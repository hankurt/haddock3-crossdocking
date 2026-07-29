#!/bin/bash
#
# HADDOCK3 ENSEMBLE CROSS-DOCKING: All Homology Models in ONE run
# Protocol: BINDING SITE RESIDUES approach
# ================================================================
#
# Difference vs. the per-model (one-by-one) script:
#   - All receptor homology models (same sequence/length, different
#     conformations, selected by DOPE) are concatenated into a single
#     multi-MODEL ensemble PDB and docked together as molecule 1.
#   - HADDOCK distributes sampling across the ensemble members.
#   - Pipeline: rigidbody -> seletop (by score) -> flexref -> emref ->
#     ilrmsdmatrix + clustrmsd -> seletopclusts.
#   - emref (energy-minimization refinement) relieves the w_vdw=0 clashes that
#     flexref (mdsteps=0) leaves, giving meaningful (negative) HADDOCK scores.
#     Alex Bonvin advises against emref for small molecules, but that advice
#     assumes a FULL flexref (MD in all 4 stages, flexible ligand). Here flexref
#     is deliberately gutted (mdsteps_rigid/cool1 = 0, and in rigid mode also
#     mdsteps_cool3 = 0 + nseg2 = 0), so nothing else removes the clashes.
#     Use --no-emref to cluster straight after flexref instead.
#   - Clustering is iL-RMSD based (ilrmsdmatrix + clustrmsd), matching HADDOCK3's
#     own protein-ligand examples. FCC is available via --fcc, but for a small
#     ligand pre-placed in ONE pocket every pose touches the same handful of BS
#     residues, so FCC saturates and lumps nearly everything into one cluster.
#   - caprieval reports i-RMSD / l-RMSD / i-l-RMSD / fnat / DockQ per pose.
#
# Everything else (restraints, cmrest/cmtight, w_vdw=0, tolerance=20,
# semi-flexible BS segments, ligand flexibility) is identical to the
# one-by-one script.
#
# Usage:
#   ./haddock_run_ensemble_crossdocking_bs.sh <models_dir> <redocking_dir> [ncores] \
#       [--rigid|--flex] [--centered] [--no-emref] [--fcc] [--min-pop N] [--clust-cutoff X]
#
# Example:
#   ./haddock_run_ensemble_crossdocking_bs.sh input_receptors redocking 32
#
# --centered (HADDOCK3 >= 2026-04-15, Issue #1491)
#   Sets `separate = false` in [rigidbody]. By default HADDOCK re-orients each
#   non-fixed molecule and pushes the molecules apart before docking, so the
#   ligand always starts outside the protein. With separate = false it keeps the
#   coordinates it is given, and only randrot (random rotation about its own
#   centre of mass) is applied.
#   No preparation is needed here: ligand_haddock.pdb already holds the crystal
#   pose in the holo frame, so it is already in the binding site. The
#   center-of-mass restraints (cmrest/cmtight), which exist only to drag a
#   scattered ligand back toward the receptor, are dropped.
#   NOTE: HADDOCK docks ONE ligand position against the whole ensemble, so the
#   models must share the holo reference frame. Checked before running.
#
# --rigid | --flex
#   Ligand flexibility in flexref/emref. Defaults to whatever the redocking
#   setup recorded in <redocking_dir>/docking_mode.txt, else rigid. Outputs go
#   to ensemble_crossdocking_bs_<mode>/ so the two never overwrite each other.
#
# --no-emref
#   Drop the [emref] stage and cluster directly on the flexref output. Only do
#   this if flexref is actually relaxing the structures (i.e. --flex, or after
#   restoring the mdsteps_* defaults); with the rigid/gutted flexref the
#   w_vdw = 0 pocket clashes survive and the HADDOCK scores stay meaningless.
#
# --fcc
#   Cluster with [clustfcc] (fraction of common contacts) instead of the default
#   [ilrmsdmatrix] + [clustrmsd]. Only useful when the ligand is large enough
#   that its contact fingerprint actually varies between poses.
#
# --min-pop N   (default 4)
#   Minimum cluster size. NOTE this is NOT only a size filter: clustfcc and
#   clustrmsd both pass min_population to rank_clusters() as the number of top
#   models averaged into the cluster score (libclust.py:rank_clusters). So
#   min_pop = 1 ranks every cluster by its single best model, which lets a
#   2-member cluster outrank a well-populated one on a single lucky pose.
#   Clusters smaller than N are dropped from the pipeline (kept only in the
#   clustering step's own structure list), which is the price of a sane ranking.
#
# --clust-cutoff X
#   Clustering cutoff. iL-RMSD (A, default 1.5) for clustrmsd; FCC fraction
#   (default 0.6) for clustfcc.
#

set -e

CENTERED=0
MODE=""
USE_EMREF=1
CLUST_METHOD="rmsd"              # rmsd (ilrmsdmatrix+clustrmsd) | fcc (clustfcc)
MIN_POP=""                       # default filled in below (4)
CLUST_CUTOFF=""                  # default filled in below (1.5 A rmsd / 0.6 fcc)
IRMSD_CUTOFF_ARG=""              # default filled in below (3.5 A)
POSARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --centered)     CENTERED=1 ;;
        --rigid)        MODE="rigid" ;;
        --flex)         MODE="flex" ;;
        --no-emref)     USE_EMREF=0 ;;
        --emref)        USE_EMREF=1 ;;
        --fcc)          CLUST_METHOD="fcc" ;;
        --rmsd-clust)   CLUST_METHOD="rmsd" ;;
        --min-pop)      MIN_POP="$2"; shift ;;
        --clust-cutoff) CLUST_CUTOFF="$2"; shift ;;
        --irmsd-cutoff) IRMSD_CUTOFF_ARG="$2"; shift ;;
        *)              POSARGS+=("$1") ;;
    esac
    shift
done

# Clustering defaults depend on the method chosen above.
if [[ "${CLUST_METHOD}" == "fcc" ]]; then
    CLUST_CUTOFF=${CLUST_CUTOFF:-0.6}    # fraction of common contacts
else
    CLUST_CUTOFF=${CLUST_CUTOFF:-1.5}    # iL-RMSD cophenetic distance (A)
fi
# HADDOCK's cluster score is the mean of the top-4 models, so anything below 4
# is flagged `under_eval=yes` by caprieval and its ranking is not comparable.
MIN_POP=${MIN_POP:-4}

MODELS_DIR=${POSARGS[0]:-"input_receptors"}
REDOCKING_DIR=${POSARGS[1]:-"redocking"}
NCORES=${POSARGS[2]:-32}
MAX_BS_DIST=10                   # max ligand-to-BS-centre distance in centered mode (A)

if [[ ! -d "${MODELS_DIR}" ]]; then
    echo "ERROR: Models directory ${MODELS_DIR} not found"
    exit 1
fi

# Check required files from redocking
for f in ligand_haddock.pdb ligand.top ligand.param ambig_active.tbl reference.pdb; do
    if [[ ! -f "${REDOCKING_DIR}/${f}" ]]; then
        echo "ERROR: ${REDOCKING_DIR}/${f} not found"
        echo "       Run redocking setup first!"
        exit 1
    fi
done

# ===========================================
# Ligand flexibility: --rigid (default) or --flex
# nseg2 = 0 pins the ligand's conformer; leaving nseg2 unset keeps HADDOCK's
# default (-1 = semi-flexible segments auto-defined from contacts), so flexref
# refines the ligand torsions. With no explicit flag, the mode recorded by the
# redocking setup is reused, so the two stages stay consistent.
# ===========================================
if [[ -z "${MODE}" ]]; then
    if [[ -s "${REDOCKING_DIR}/docking_mode.txt" ]]; then
        MODE=$(tr -d '[:space:]' < "${REDOCKING_DIR}/docking_mode.txt")
        MODE_SRC="${REDOCKING_DIR}/docking_mode.txt"
    else
        MODE="rigid"
        MODE_SRC="default"
    fi
else
    MODE_SRC="command line"
fi
if [[ "${MODE}" != "rigid" && "${MODE}" != "flex" ]]; then
    echo "ERROR: unknown ligand mode '${MODE}' (expected rigid or flex)"
    exit 1
fi

if [[ "${MODE}" == "flex" ]]; then
    FLEX_NOTE="FLEXIBLE (nseg2 auto, ligand torsions refined)"
    FLEXREF_LIG=""
    EMREF_LIG=""
else
    FLEX_NOTE="RIGID (nseg2 = 0, input conformer kept)"
    FLEXREF_LIG="mdsteps_cool3 = 0
nseg2 = 0"
    EMREF_LIG="nseg2 = 0"
fi

# caprieval interface-definition cutoff (--irmsd-cutoff, default 3.5).
#
# This is the distance (A) used to decide which residues count as interface,
# NOT a quality threshold. i-RMSD, l-RMSD, i-l-RMSD, fnat and DockQ are all
# reported regardless of its value. HADDOCK's own default is 10.0.
#
# It drives BOTH metrics: caprieval passes params["irmsd_cutoff"] to calc_irmsd
# AND to calc_ilrmsd, so it sets the i-RMSD atom set and the i-l-RMSD *fit*
# basis at the same time, and the two respond in OPPOSITE directions.
#
# The i-RMSD atom set is receptor-interface + ligand. The ligand contributes a
# fixed number of atoms while the receptor interface grows with the cutoff, so
# the ligand's share of the i-RMSD signal falls as the cutoff rises. For a small
# molecule this is a big effect: tight cutoff -> ligand is a large fraction of
# the atoms and i-RMSD behaves like a binding-mode metric; at HADDOCK's 10.0 the
# receptor backbone dominates, which across an ensemble of different receptor
# conformers measures model spread rather than pose accuracy.
#
# i-l-RMSD (fit on the receptor interface, measure the ligand) is 100 % ligand
# at any cutoff, so it never suffers that dilution -- but a larger cutoff gives
# it a broader, more stable superposition basis. Raising the cutoff therefore
# improves i-l-RMSD while degrading i-RMSD. See the README for a snippet that
# counts the atoms at each cutoff for your own reference.pdb.
#
# NOTE: this does NOT affect clustering. [ilrmsdmatrix] has its own
# contact_distance_cutoff (default 5.0 A), independent of this parameter.
IRMSD_CUTOFF=${IRMSD_CUTOFF_ARG:-3.5}

OUTDIR="ensemble_crossdocking_bs_${MODE}"

# Count models and derive scaled sampling
NMODELS=$(ls "${MODELS_DIR}"/*.pdb 2>/dev/null | wc -l | tr -d ' ')
if [[ "${NMODELS}" -eq 0 ]]; then
    echo "ERROR: No .pdb models found in ${MODELS_DIR}"
    exit 1
fi
# Rigidbody sampling = 250 per conformer, exactly like the one-by-one runs.
# The iL-RMSD matrix sits AFTER seletop, so the 10k cap applies to SELETOP_FLEX,
# not to this number.
SAMPLING=$((250 * NMODELS))

# --- Stage-1 (AFTER rigidbody): SCORE-based funnel to flexref (no iL-RMSD) ---
# Keep 50 per conformer (50*N), matching one-by-one's 250 rigid -> 50 flex.
# NOTE: this is a GLOBAL top-(50*N) by score, not a guaranteed 50-per-conformer
# (HADDOCK ensembles can't do per-conformer seletop) -- check the per-model
# spread in the caprieval tables.
SELETOP_FLEX=$((50 * NMODELS))   # = 5050 for N=101 -> flexref

# ilrmsdmatrix refuses to build a matrix above max_models (default 10000).
if [[ "${CLUST_METHOD}" == "rmsd" && "${SELETOP_FLEX}" -gt 10000 ]]; then
    echo "ERROR: seletop = ${SELETOP_FLEX} exceeds ilrmsdmatrix max_models (10000)."
    echo "       Lower the 50-per-conformer factor, or pass --fcc to cluster on contacts."
    exit 1
fi

# --- Stage-2 clustering ---
# Default: iL-RMSD (ilrmsdmatrix + clustrmsd), as in HADDOCK3's own
# docking-protein-ligand examples. This is what actually separates binding
# modes; FCC (--fcc) saturates when a small ligand sits in a single pocket,
# because every pose then contacts the same BS residues.
if [[ "${CLUST_METHOD}" == "fcc" ]]; then
    CLUST_NOTE="clustfcc (contacts), fcc_cutoff=${CLUST_CUTOFF}"
else
    CLUST_NOTE="ilrmsdmatrix + clustrmsd (distance/average), cutoff=${CLUST_CUTOFF} A"
fi

echo "=============================================="
echo "HADDOCK3 ENSEMBLE CROSS-DOCKING: BS Residues"
echo "Ligand:       ${FLEX_NOTE} [${MODE_SRC}]"
echo "=============================================="
echo "Models dir:   ${MODELS_DIR} (N=${NMODELS} conformers)"
echo "Redocking:    ${REDOCKING_DIR}"
echo "Output:       ${OUTDIR}"
echo "Cores:        ${NCORES}"
echo "----------------------------------------------"
echo "  Rigidbody:   sampling = ${SAMPLING}"
if [[ "${CENTERED}" -eq 1 ]]; then
    echo "  Placement:   CENTERED (ligand kept in the BS, separate=false)"
    echo "  w_vdw=0.0, tolerance=20, cmrest/cmtight OFF (not needed when centered)"
else
    echo "  Placement:   SCATTERED (HADDOCK default, separate=true)"
    echo "  w_vdw=0.0, tolerance=20, cmrest=true, cmtight=true"
fi
echo "  Stage-1 (post-rigidbody): seletop = ${SELETOP_FLEX} (by HADDOCK score) -> flexref"
if [[ "${USE_EMREF}" -eq 1 ]]; then
    echo "  Refinement: flexref -> emref (EM) -> meaningful (negative) scores"
    echo "  Stage-2 (post-emref):  ${CLUST_NOTE}, min_pop=${MIN_POP}"
else
    echo "  Refinement: flexref only (emref disabled via --no-emref)"
    if [[ "${MODE}" == "rigid" ]]; then
        echo "  WARNING: rigid mode gives a near no-op flexref (mdsteps 0, nseg2 = 0)."
        echo "           Without emref the w_vdw=0 clashes survive and HADDOCK scores"
        echo "           stay large and positive. Check the score range in caprieval."
    fi
    echo "  Stage-2 (post-flexref): ${CLUST_NOTE}, min_pop=${MIN_POP}"
fi
echo "  Analysis: caprieval only (irmsd_cutoff=${IRMSD_CUTOFF} = interface definition)"
echo "=============================================="
echo ""

mkdir -p "${OUTDIR}"

# ===========================================
# Extract BS residues from ambig_active.tbl
# (identical to one-by-one script)
# ===========================================
echo "==== DETERMINING BS RESIDUES ===="
# Preferred: the plain list written by haddock_run_redocking_bs.sh.
# Fallback: parse the segid-A selections out of the AIR table (BSD/GNU-safe).
if [[ -s "${REDOCKING_DIR}/bs_residues.txt" ]]; then
    BS_RESIDUES=$(tr '\n' ' ' < "${REDOCKING_DIR}/bs_residues.txt")
    echo "  Source: ${REDOCKING_DIR}/bs_residues.txt"
else
    BS_RESIDUES=$(grep -oE 'resid[[:space:]]+[0-9]+[[:space:]]+and[[:space:]]+segid[[:space:]]+A' \
                    "${REDOCKING_DIR}/ambig_active.tbl" | \
                  grep -oE '[0-9]+' | sort -un | tr '\n' ' ')
    echo "  Source: ${REDOCKING_DIR}/ambig_active.tbl"
fi
echo "  BS residues: ${BS_RESIDUES}"
if [[ -z "${BS_RESIDUES}" ]]; then
    echo "ERROR: No BS residues found in ${REDOCKING_DIR}"
    exit 1
fi

# ===========================================
# Generate semi-flexible segment definitions
# (identical to one-by-one script)
# ===========================================
echo "==== GENERATING SEMI-FLEXIBLE SEGMENTS ===="
SEMIFLEX_BLOCK=$(python3 << PYEOF
bs_residues = [${BS_RESIDUES// /, }]

ranges = []
for r in sorted(set(bs_residues)):
    ranges.append((r - 1, r + 1))

merged = [ranges[0]]
for start, end in ranges[1:]:
    prev_start, prev_end = merged[-1]
    if start <= prev_end + 1:
        merged[-1] = (prev_start, max(prev_end, end))
    else:
        merged.append((start, end))

lines = []
lines.append("# Semi-flexible segments for receptor (molecule 1)")
lines.append(f"# {len(merged)} merged segments from {len(bs_residues)} BS residues")
for i, (start, end) in enumerate(merged, start=1):
    contained = [r for r in sorted(set(bs_residues)) if start <= r <= end]
    lines.append(f"# Segment {i}: res {','.join(map(str, contained))}")
    lines.append(f"seg_sta_1_{i} = {start}")
    lines.append(f"seg_end_1_{i} = {end}")

print('\n'.join(lines))
PYEOF
)
echo "${SEMIFLEX_BLOCK}"
echo ""

# ===========================================
# BUILD THE RECEPTOR ENSEMBLE
# Each model -> ATOM-only, chain A; then all
# models concatenated into one multi-MODEL PDB.
# ===========================================
echo "==== BUILDING RECEPTOR ENSEMBLE ===="

WORKDIR="${OUTDIR}"
RECEPTOR_PARTS=()
idx=0
for MODEL_PDB in "${MODELS_DIR}"/*.pdb; do
    idx=$((idx + 1))
    PART="${WORKDIR}/receptor_$(printf '%03d' ${idx}).pdb"

    grep "^ATOM" "${MODEL_PDB}" > "${WORKDIR}/_tmp.pdb" || true
    if command -v pdb_chain >/dev/null 2>&1; then
        pdb_chain -A "${WORKDIR}/_tmp.pdb" > "${PART}" 2>/dev/null || \
            sed 's/^\(.\{21\}\)./\1A/' "${WORKDIR}/_tmp.pdb" > "${PART}"
    else
        sed 's/^\(.\{21\}\)./\1A/' "${WORKDIR}/_tmp.pdb" > "${PART}"
    fi
    rm -f "${WORKDIR}/_tmp.pdb"
    RECEPTOR_PARTS+=("${PART}")
done

ENSEMBLE="${WORKDIR}/receptor_ensemble.pdb"
if command -v pdb_mkensemble >/dev/null 2>&1; then
    pdb_mkensemble "${RECEPTOR_PARTS[@]}" > "${ENSEMBLE}"
else
    # Fallback: wrap each model in MODEL/ENDMDL, drop END/TER between
    echo "  (pdb_mkensemble not found -- using python fallback)"
    python3 - "${ENSEMBLE}" "${RECEPTOR_PARTS[@]}" << 'PYEOF'
import sys
out_path = sys.argv[1]
parts = sys.argv[2:]
with open(out_path, "w") as out:
    for i, p in enumerate(parts, start=1):
        out.write(f"MODEL     {i:>4}\n")
        with open(p) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM", "ANISOU", "TER")):
                    out.write(line)
        out.write("ENDMDL\n")
    out.write("END\n")
PYEOF
fi

NENS=$(grep -c '^MODEL' "${ENSEMBLE}" || true)
echo "  Ensemble written: ${ENSEMBLE} (${NENS} models)"
if [[ "${NENS}" -ne "${NMODELS}" ]]; then
    echo "WARNING: ensemble model count (${NENS}) != input count (${NMODELS})"
fi

# ===========================================
# Copy redocking inputs into the run dir
# ===========================================
cp "${REDOCKING_DIR}/ligand_haddock.pdb" "${WORKDIR}/"
cp "${REDOCKING_DIR}/ligand.top"         "${WORKDIR}/"
cp "${REDOCKING_DIR}/ligand.param"       "${WORKDIR}/"
cp "${REDOCKING_DIR}/ambig_active.tbl"   "${WORKDIR}/"
cp "${REDOCKING_DIR}/reference.pdb"      "${WORKDIR}/"

# ===========================================
# CENTERED MODE
# ligand_haddock.pdb already holds the crystal pose in the holo frame, so it is
# the centered ligand as-is: separate = false simply tells HADDOCK to leave it
# there. The center-of-mass restraints, which exist only to pull a scattered
# ligand back toward the receptor, are switched off.
# The only requirement is that the models share the holo reference frame; the
# check below catches the case where they do not.
# ===========================================
if [[ "${CENTERED}" -eq 1 ]]; then
    PLACEMENT_BLOCK="separate = false"

    echo ""
    echo "==== CHECKING LIGAND SITS IN THE BINDING SITE ===="
    for PART in "${RECEPTOR_PARTS[@]}"; do
        awk -v bs="${BS_RESIDUES}" -v lig="${WORKDIR}/ligand_haddock.pdb" \
            -v label="${PART}" -v maxd="${MAX_BS_DIST}" '
            BEGIN { n = split(bs, a, " "); for (i = 1; i <= n; i++) want[a[i] + 0] = 1 }
            substr($0,1,4) == "ATOM" && (substr($0,23,4) + 0) in want {
                bx += substr($0,31,8); by += substr($0,39,8); bz += substr($0,47,8); bn++
            }
            END {
                if (bn == 0) {
                    printf("ERROR: no binding-site atoms found in %s\n", label); exit 1
                }
                while ((getline l < lig) > 0)
                    if (substr(l,1,6) ~ /^(ATOM|HETATM)/) {
                        lx += substr(l,31,8); ly += substr(l,39,8); lz += substr(l,47,8); ln++
                    }
                if (ln == 0) { printf("ERROR: no atoms in %s\n", lig); exit 1 }
                d = sqrt((lx/ln - bx/bn)^2 + (ly/ln - by/bn)^2 + (lz/ln - bz/bn)^2)
                if (d > maxd) {
                    printf("  %s: ligand is %.1f A from the BS centre (limit %.1f)\n", label, d, maxd)
                    exit 1
                }
                printf("  %s: OK (%.1f A from the BS centre)\n", label, d)
            }' "${PART}" || {
            echo ""
            echo "ERROR: the ligand does not sit in this model's binding site."
            echo "       In centered mode HADDOCK docks ONE ligand position against the"
            echo "       whole ensemble, so every model must share the reference frame of"
            echo "       ${REDOCKING_DIR}/receptor.pdb (the holo structure)."
            echo "       Fix: superimpose the models onto the holo structure and re-run,"
            echo "       or drop --centered, or raise MAX_BS_DIST (currently ${MAX_BS_DIST} A)."
            exit 1
        }
    done
else
    PLACEMENT_BLOCK="cmrest = true
cmtight = true"
fi

# ===========================================
# Assemble the optional [emref] stage.
#
# There is deliberately NO [caprieval] between emref and the clustering step:
# clustering does not touch coordinates or scores, it only fills in the
# cluster_id / cluster_ranking / model-cluster_ranking columns. A caprieval on
# either side of it therefore produces byte-identical per-model metrics, and the
# post-clustering table is a strict superset. One caprieval after clustering.
# ===========================================
if [[ "${USE_EMREF}" -eq 1 ]]; then
    EMREF_BLOCK="# ===========================================
# Energy-minimization refinement (emref)
# Relieves the w_vdw=0 pocket-penetration clashes that
# flexref (mdsteps=0) leaves behind, so HADDOCK scores
# become meaningful (negative for good complexes).
# NOTE: EM cannot cross energy barriers -> severely
# clashed (non-native, high-RMSD) poses may stay
# positive; switch to [mdref] if you must relax them.
# Confirm params with 'haddock3-cfg -m emref'.
# ===========================================
[emref]
tolerance = 5
ambig_fname = \"ambig_active.tbl\"
ligand_param_fname = \"ligand.param\"
ligand_top_fname = \"ligand.top\"
${EMREF_LIG}
randremoval = false
"
else
    EMREF_BLOCK="# emref disabled (--no-emref): clustering runs on the flexref output.
"
fi

# ===========================================
# Assemble the clustering stage.
# ===========================================
if [[ "${CLUST_METHOD}" == "fcc" ]]; then
    CLUSTER_BLOCK="# ===========================================
# STAGE-2 CLUSTERING -- FCC (contact based)
# Only worth using when the ligand is big enough that
# its residue-contact fingerprint actually differs
# between poses. For a small ligand pre-placed in one
# pocket, every pose contacts the same BS residues, FCC
# saturates near 1.0 and you get one giant cluster.
# NOTE: confirm params with 'haddock3-cfg -m clustfcc'
# ===========================================
[clustfcc]
clust_cutoff = ${CLUST_CUTOFF}
min_population = ${MIN_POP}
"
else
    CLUSTER_BLOCK="# ===========================================
# STAGE-2 CLUSTERING -- iL-RMSD (HADDOCK3 default for
# protein-ligand; see examples/docking-protein-ligand).
# ilrmsdmatrix superposes on the receptor interface and
# measures the ligand, so genuinely different binding
# modes end up in different clusters.
# clust_cutoff is a cophenetic distance in Angstrom.
# NOTE: confirm params with 'haddock3-cfg -m clustrmsd'
# ===========================================
[ilrmsdmatrix]

[clustrmsd]
criterion = 'distance'
linkage = 'average'
clust_cutoff = ${CLUST_CUTOFF}
min_population = ${MIN_POP}
"
fi

# ===========================================
# WRITE THE SINGLE ENSEMBLE CONFIG
# ===========================================
cat > "${WORKDIR}/haddock3_ensemble_bs_${MODE}.cfg" << CFGEOF
# HADDOCK3 ENSEMBLE Cross-docking: ${NMODELS} homology models
# Protocol: Binding Site Residues approach - RIGID LIGAND
run_dir = "run_ensemble_crossdock_bs_${MODE}"
ncores = ${NCORES}
mode = "local"

molecules = [
    "receptor_ensemble.pdb",
    "ligand_haddock.pdb"
]

# ===========================================
# Topology generation
# (generated for every ensemble member)
# ===========================================
[topoaa]
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
delenph = false
autohis = true
set_bfactor = false

# ===========================================
# Rigid-body docking
# sampling scaled by N so each conformer is
# sampled ~250 times.
# ===========================================
[rigidbody]
sampling = ${SAMPLING}
tolerance = 20
ambig_fname = "ambig_active.tbl"
w_vdw = 0.0
inter_rigid = 0.001
mol_fix_origin_1 = true
${PLACEMENT_BLOCK}
randremoval = false
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# STAGE-1 (after rigidbody): SCORE-based funnel
# Keep the top ${SELETOP_FLEX} rigidbody models by
# HADDOCK score and pass them to flexref. This also
# keeps the later iL-RMSD matrix under its 10k cap.
# ===========================================
[seletop]
select = ${SELETOP_FLEX}

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Flexible refinement -- ligand: ${FLEX_NOTE}
# semi-flexible BS side-chains only
# ===========================================
[flexref]
tolerance = 50
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
mdsteps_rigid = 0
mdsteps_cool1 = 0
${FLEXREF_LIG}
randremoval = false
${SEMIFLEX_BLOCK}

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

${EMREF_BLOCK}
${CLUSTER_BLOCK}
# ===========================================
# Cluster ranking + per-pose metrics.
# capri_clt.tsv holds the cluster table; capri_ss.tsv
# holds every pose with its cluster_id / cluster_ranking.
# ===========================================
[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Keep clusters / models for reporting.
# This ONLY trims which cluster PDBs get copied
# out -- it does NOT change the ranking above.
# NOTE: models in clusters smaller than min_population
# are NOT carried past the clustering step, so they
# cannot reappear here. Use the pre-clustering
# caprieval tables if you need every pose.
# ===========================================
[seletopclusts]
top_clusters = 1000
top_models = 9999

# ===========================================
# Final all-atom CAPRI evaluation: also scores
# protein side chains in alignment + i/iL-RMSD.
# ===========================================
[caprieval]
allatoms = true
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}
CFGEOF

echo ""
echo "  Config written: ${WORKDIR}/haddock3_ensemble_bs_${MODE}.cfg"

# ===========================================
# RUN (single ensemble run)
# ===========================================
echo ""
echo "==== RUNNING HADDOCK3 (ENSEMBLE) ===="
STARTTIME=$(date +%s)
if (cd "${WORKDIR}" && haddock3 haddock3_ensemble_bs_${MODE}.cfg > haddock3_ensemble_${MODE}.log 2>&1); then
    ELAPSED=$(( $(date +%s) - STARTTIME ))
    echo "  DONE (${ELAPSED}s)"
else
    ELAPSED=$(( $(date +%s) - STARTTIME ))
    echo "  FAILED (${ELAPSED}s) -- check ${WORKDIR}/haddock3_ensemble_${MODE}.log"
    exit 1
fi

echo ""
echo "=============================================="
echo "ENSEMBLE RUN COMPLETE"
echo "  Results: ${WORKDIR}/run_ensemble_crossdock_bs"
echo "  Ranked poses:    .../analysis/*_caprieval/capri_ss.tsv"
echo "  Ranked clusters: .../analysis/*_caprieval/capri_clt.tsv"
echo "=============================================="
