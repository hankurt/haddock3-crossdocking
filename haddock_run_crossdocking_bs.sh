#!/bin/bash
#
# HADDOCK3 CROSS-DOCKING: Batch for All Homology Models
# Protocol: BINDING SITE RESIDUES approach
# =======================================================
#
# Strategy for closed binding sites:
#   - cmrest/cmtight: center-of-mass restraints keep ligand near receptor
#   - w_vdw=0.0: removes vdW penalty to allow pocket penetration
#   - tolerance=20: loose restraints to avoid surface trapping
#   - Semi-flexible BS residues in flexref (side-chain flexibility only)
#   - ambig_active.tbl used for both rigidbody and flexref
#   - Ligand rigid (nseg2=0) or flexible (nseg2 auto), see --rigid/--flex
#   - emref after flexref: with w_vdw=0 and mdsteps_*=0 nothing else relieves
#     the pocket clashes, so without it the HADDOCK scores stay large and
#     positive. This also keeps the scores on the same scale as
#     haddock_run_ensemble_crossdocking_bs.sh, which refines the same way.
#
# Key settings:
#   - Rigidbody: sampling=250, w_vdw=0.0, tolerance=20, cmrest/cmtight
#   - Flexref: mdsteps_rigid=0, mdsteps_cool1=0, semi-flexible BS segments
#   - emref: tolerance=5, randremoval=false (as in redocking/ensemble)
#   - Pipeline: rigidbody -> seletop(50) -> flexref -> emref -> seletop(5)
#   - Scoring/analysis: caprieval only (i-RMSD, l-RMSD, i-l-RMSD, fnat, DockQ)
#
# Usage:
#   ./haddock_run_crossdocking_bs.sh <models_dir> <redocking_dir> [ncores] [--rigid|--flex] [--centered]
#
# Example:
#   ./haddock_run_crossdocking_bs.sh input_receptors redocking 32
#
# --centered (HADDOCK3 >= 2026-04-15, Issue #1491)
#   Sets `separate = false` in [rigidbody]. By default HADDOCK re-orients each
#   non-fixed molecule and pushes the molecules apart before docking, so the
#   ligand always starts outside the protein. With separate = false it keeps the
#   coordinates it is given, and only randrot (random rotation about its own
#   centre of mass) is applied.
#   No preparation is needed here: ligand_haddock.pdb already holds the crystal
#   pose in the holo frame, so it is already in the binding site, provided the
#   models share that frame (checked before running). The center-of-mass
#   restraints (cmrest/cmtight), which exist only to drag a scattered ligand
#   back toward the receptor, are dropped.
#
# --rigid | --flex
#   Ligand flexibility in flexref. Defaults to whatever the redocking setup
#   recorded in <redocking_dir>/docking_mode.txt, else rigid. Outputs go to
#   crossdocking_runs_bs_<mode>/ so the two never overwrite each other.
#

set -e

CENTERED=0
MODE=""
IRMSD_CUTOFF_ARG=""
POSARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --centered)      CENTERED=1 ;;
        --rigid)         MODE="rigid" ;;
        --flex)          MODE="flex" ;;
        --irmsd-cutoff)  IRMSD_CUTOFF_ARG="$2"; shift ;;
        *)               POSARGS+=("$1") ;;
    esac
    shift
done

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

OUTDIR="crossdocking_runs_bs_${MODE}"

# Count models
NMODELS=$(ls "${MODELS_DIR}"/*.pdb 2>/dev/null | wc -l)
echo "=============================================="
echo "HADDOCK3 CROSS-DOCKING: BS Residues Approach"
echo "Ligand:       ${FLEX_NOTE} [${MODE_SRC}]"
echo "=============================================="
echo "Models dir:   ${MODELS_DIR} (${NMODELS} models)"
echo "Redocking:    ${REDOCKING_DIR}"
echo "Output:       ${OUTDIR}"
echo "Cores:        ${NCORES}"
echo "=============================================="
echo "  Rigidbody: sampling=250, w_vdw=0.0, tolerance=20"
if [[ "${CENTERED}" -eq 1 ]]; then
    echo "  Placement: CENTERED (ligand kept in the BS, separate=false)"
    echo "  cmrest/cmtight OFF (not needed when centered)"
else
    echo "  Placement: SCATTERED (HADDOCK default, separate=true)"
    echo "  cmrest=true, cmtight=true"
fi
echo "  Flexref: semi-flexible BS, mdsteps_rigid=0, mdsteps_cool1=0"
echo "  Analysis: caprieval (irmsd_cutoff=${IRMSD_CUTOFF} = interface definition)"
echo "=============================================="
echo ""

mkdir -p "${OUTDIR}"

# ===========================================
# Extract BS residues from ambig_active.tbl
# Parse "resi X" from segid A (receptor only)
# Filter out resi 1 (ligand, segid B)
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
# for flexref/emref from BS residues
# Each BS residue gets +/-1 padding, then
# consecutive/overlapping ranges are merged
# ===========================================
echo "==== GENERATING SEMI-FLEXIBLE SEGMENTS ===="

SEMIFLEX_BLOCK=$(python3 << PYEOF
import sys

bs_residues = [${BS_RESIDUES// /, }]

# Create padded ranges and merge overlapping
ranges = []
for r in sorted(set(bs_residues)):
    ranges.append((r - 1, r + 1))

# Merge overlapping/adjacent ranges
merged = [ranges[0]]
for start, end in ranges[1:]:
    prev_start, prev_end = merged[-1]
    if start <= prev_end + 1:
        merged[-1] = (prev_start, max(prev_end, end))
    else:
        merged.append((start, end))

# Generate segment definitions
lines = []
lines.append("# Semi-flexible segments for receptor (molecule 1)")
lines.append(f"# {len(merged)} merged segments from {len(bs_residues)} BS residues")
for i, (start, end) in enumerate(merged, start=1):
    # Collect the original BS residues in this segment
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
# SETUP: Create all working directories
# ===========================================
echo ""
echo "==== SETTING UP ALL MODELS ===="

for MODEL_PDB in "${MODELS_DIR}"/*.pdb; do
    MODEL_NAME=$(basename "${MODEL_PDB}" .pdb)
    WORKDIR="${OUTDIR}/model_${MODEL_NAME}"
    mkdir -p "${WORKDIR}"

    # Prepare receptor
    grep "^ATOM" "${MODEL_PDB}" > "${WORKDIR}/receptor_tmp.pdb" || true
    pdb_chain -A "${WORKDIR}/receptor_tmp.pdb" > "${WORKDIR}/receptor.pdb" 2>/dev/null || {
        sed 's/^\(.\{21\}\)./\1A/' "${WORKDIR}/receptor_tmp.pdb" > "${WORKDIR}/receptor.pdb"
    }
    rm -f "${WORKDIR}/receptor_tmp.pdb"
    echo "END" >> "${WORKDIR}/receptor.pdb"

    # Copy redocking files
    cp "${REDOCKING_DIR}/ligand_haddock.pdb" "${WORKDIR}/"
    cp "${REDOCKING_DIR}/ligand.top"         "${WORKDIR}/"
    cp "${REDOCKING_DIR}/ligand.param"       "${WORKDIR}/"
    cp "${REDOCKING_DIR}/ambig_active.tbl"   "${WORKDIR}/"
    cp "${REDOCKING_DIR}/reference.pdb"      "${WORKDIR}/"

    # --- Ligand placement ---
    # In centered mode the ligand is left exactly where ligand_haddock.pdb puts
    # it (the crystal pose, in the holo frame), so the center-of-mass restraints
    # that exist to pull a scattered ligand back are switched off. This only
    # holds if the model shares the holo frame, hence the check.
    if [[ "${CENTERED}" -eq 1 ]]; then
        PLACEMENT_BLOCK="separate = false"
        PLACEMENT_NOTE="ligand kept in the BS (separate = false)"

        awk -v bs="${BS_RESIDUES}" -v lig="${WORKDIR}/ligand_haddock.pdb" \
            -v label="${MODEL_NAME}" -v maxd="${MAX_BS_DIST}" '
            BEGIN { n = split(bs, a, " "); for (i = 1; i <= n; i++) want[a[i] + 0] = 1 }
            substr($0,1,4) == "ATOM" && (substr($0,23,4) + 0) in want {
                bx += substr($0,31,8); by += substr($0,39,8); bz += substr($0,47,8); bn++
            }
            END {
                if (bn == 0) {
                    printf("    ERROR: no binding-site atoms found in %s\n", label); exit 1
                }
                while ((getline l < lig) > 0)
                    if (substr(l,1,6) ~ /^(ATOM|HETATM)/) {
                        lx += substr(l,31,8); ly += substr(l,39,8); lz += substr(l,47,8); ln++
                    }
                if (ln == 0) { printf("    ERROR: no atoms in %s\n", lig); exit 1 }
                d = sqrt((lx/ln - bx/bn)^2 + (ly/ln - by/bn)^2 + (lz/ln - bz/bn)^2)
                if (d > maxd) {
                    printf("    %s: ligand is %.1f A from the BS centre (limit %.1f)\n", label, d, maxd)
                    exit 1
                }
            }' "${WORKDIR}/receptor.pdb" || {
            echo ""
            echo "ERROR: the ligand does not sit in this model's binding site."
            echo "       Centered mode reuses ligand_haddock.pdb as-is, so the models"
            echo "       must share the reference frame of the holo structure used for"
            echo "       the redocking setup."
            echo "       Fix: superimpose the models onto the holo structure and re-run,"
            echo "       or drop --centered, or raise MAX_BS_DIST (currently ${MAX_BS_DIST} A)."
            exit 1
        }
    else
        PLACEMENT_BLOCK="cmrest = true
cmtight = true"
        PLACEMENT_NOTE="cmrest/cmtight handle ligand positioning"
    fi

    # --- Config: RIGID LIGAND ---
    cat > "${WORKDIR}/haddock3_crossdock_bs_${MODE}.cfg" << CFGEOF
# HADDOCK3 Cross-docking: Model ${MODEL_NAME}
# Protocol: Binding Site Residues approach - RIGID LIGAND
# ${PLACEMENT_NOTE}
run_dir = "run_crossdock_bs_${MODE}"
ncores = ${NCORES}

molecules = [
    "receptor.pdb",
    "ligand_haddock.pdb"
]

# ===========================================
# Topology generation
# ===========================================
[topoaa]
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
delenph = false
autohis = true
set_bfactor = false

# ===========================================
# Rigid-body docking
# BS residues = ACTIVE, Ligand = ACTIVE
# Closed BS strategy: w_vdw=0, tolerance=20,
# ${PLACEMENT_NOTE}
# ===========================================
[rigidbody]
sampling = 250
tolerance = 20
ambig_fname = "ambig_active.tbl"
w_vdw = 0.0
inter_rigid = 0.001
mol_fix_origin_1 = true
${PLACEMENT_BLOCK}
randremoval = false
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"

# ===========================================
# CAPRI evaluation after rigid-body
# ===========================================
[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Selection of top 50
# ===========================================
[seletop]
select = 50

# ===========================================
# CAPRI evaluation after selection
# ===========================================
[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Flexible refinement
# BS residues = ACTIVE, Ligand = ACTIVE
# Semi-flexible BS side-chains only
# nseg2 = 0: LIGAND KEPT RIGID
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

# ===========================================
# CAPRI evaluation after flexref
# ===========================================
[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Energy-minimization refinement -- ligand: ${FLEX_NOTE}
# Relieves the w_vdw = 0 pocket-penetration clashes that
# the gutted flexref (mdsteps_rigid/cool1 = 0, and in rigid
# mode also mdsteps_cool3 = 0 + nseg2 = 0) leaves behind.
# Without it the HADDOCK scores stay large and positive and
# are NOT on the same scale as the ensemble cross-docking
# runs, which refine with emref. Matches the emref stage in
# haddock_run_redocking_bs.sh and in the ensemble script.
# ===========================================
[emref]
tolerance = 5
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
${EMREF_LIG}
randremoval = false

# ===========================================
# CAPRI evaluation after emref
# ===========================================
[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Final selection and evaluation
# ===========================================
[seletop]
select = 5

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}
CFGEOF

    echo "  Setup: model_${MODEL_NAME}"
done

echo ""
echo "  All ${NMODELS} models set up in ${OUTDIR}/"

# ===========================================
# RUN: Sequential execution
# ===========================================
echo ""
echo "==== RUNNING HADDOCK3 FOR ALL MODELS ===="
echo ""

FAILED=0
SUCCESS=0
STARTTIME=$(date +%s)

for MODEL_PDB in "${MODELS_DIR}"/*.pdb; do
    MODEL_NAME=$(basename "${MODEL_PDB}" .pdb)
    WORKDIR="${OUTDIR}/model_${MODEL_NAME}"

    if [[ -d "${WORKDIR}/run_crossdock_bs_rigid" ]]; then
        echo "[SKIP] model_${MODEL_NAME} (already exists)"
    else
        echo "[RUN]  model_${MODEL_NAME} ..."
        TSTART=$(date +%s)
        if (cd "${WORKDIR}" && haddock3 haddock3_crossdock_bs_${MODE}.cfg > haddock3_${MODE}.log 2>&1); then
            TEND=$(date +%s)
            ELAPSED=$((TEND - TSTART))
            echo "       DONE (${ELAPSED}s)"
            SUCCESS=$((SUCCESS + 1))
        else
            TEND=$(date +%s)
            ELAPSED=$((TEND - TSTART))
            echo "       FAILED (${ELAPSED}s) -- check ${WORKDIR}/haddock3_${MODE}.log"
            FAILED=$((FAILED + 1))
        fi
    fi
done

ENDTIME=$(date +%s)
TOTAL_ELAPSED=$((ENDTIME - STARTTIME))

echo ""
echo "=============================================="
echo "BATCH COMPLETE"
echo "=============================================="
echo "  Success: ${SUCCESS}"
echo "  Failed:  ${FAILED}"
echo "  Total time: ${TOTAL_ELAPSED}s"
echo "=============================================="
