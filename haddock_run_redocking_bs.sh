#!/bin/bash
#
# HADDOCK3 REDOCKING SETUP - BINDING SITE RESIDUES Approach
# =========================================================
#
# INPUT:  Holo structure (protein + ligand complex) + binding-site residue list
# OUTPUT: A self-contained <outdir>/ (default: redocking/) holding every input
#         the cross-docking runners need, plus a ready-to-run redocking config.
#
# Protocol (HADDOCK protein-ligand binding site tutorial):
#   https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-binding-sites/
#   - BS residues = ACTIVE, ligand = ACTIVE  (AIRs draw the ligand into the pocket)
#   - 2-body docking (no shape molecule)
#   - Ligand topology/parameters from the PRODRG server
#
# LIGAND FLEXIBILITY (--rigid | --flex)
#   --rigid (default)  nseg2 = 0 in flexref/emref: the ligand has no
#                      semi-flexible segments, so it keeps its input conformer
#                      throughout. Comparable to rigid-ligand AutoDock/Vina.
#   --flex             nseg2 left at its HADDOCK default (-1), so semi-flexible
#                      segments are defined automatically from contacts and the
#                      ligand torsions are refined during flexref/emref. This is
#                      the standard HADDOCK protein-ligand protocol.
#   The mode only changes the config; every shared input file is identical, so
#   one setup directory can serve both. The mode is recorded in
#   <outdir>/docking_mode.txt and picked up by the cross-docking runners.
#
# Requirements:
#   - Python 3 with: rdkit, requests (not HADDOCK3 deps -- install separately)
#     numpy, scipy and pdb-tools come with HADDOCK3
#   - pdb_chain (from pdb-tools; a sed fallback is used if absent)
#   - haddock3 (only to actually run the docking, see --run)
#
# Usage:
#   ./haddock_run_redocking_bs.sh <holo_complex.pdb> <ligand_resname> <bs_residues> [ncores] [outdir] [--rigid|--flex] [--run]
#
# Example:
#   ./haddock_run_redocking_bs.sh 6dpt_A.pdb H7M "64 67 120 152 211 316 318" 16
#
# The resulting directory is what the cross-docking runners consume:
#   ./haddock_run_ensemble_crossdocking_bs.sh <models_dir> redocking 32
#

set -e

# ===========================================
# Argument parsing
# ===========================================
RUN_HADDOCK=0
MODE="rigid"
IRMSD_CUTOFF_ARG=""
POSARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run)           RUN_HADDOCK=1 ;;
        --rigid)         MODE="rigid" ;;
        --flex)          MODE="flex" ;;
        --irmsd-cutoff)  IRMSD_CUTOFF_ARG="$2"; shift ;;
        *)               POSARGS+=("$1") ;;
    esac
    shift
done

# caprieval interface-definition cutoff (--irmsd-cutoff, default 3.5).
#
# This is the distance (A) used to decide which residues count as interface,
# NOT a quality threshold. i-RMSD, l-RMSD, i-l-RMSD, fnat and DockQ are all
# reported regardless of its value. HADDOCK's own default is 10.0.
#
# It drives BOTH metrics: caprieval passes params["irmsd_cutoff"] to calc_irmsd
# AND to calc_ilrmsd, so it sets the i-RMSD atom set and the i-l-RMSD *fit*
# basis at the same time. A tight cutoff makes i-RMSD ligand-dominated (and so
# usable as a binding-mode metric); at HADDOCK's 10.0 the receptor backbone
# dominates it. i-l-RMSD (fit on the receptor interface, measure the ligand) is
# the binding-mode metric either way, and a larger cutoff gives it a more
# stable fit basis. Keep this value identical across redocking and
# cross-docking or the tables are not comparable.
IRMSD_CUTOFF=${IRMSD_CUTOFF_ARG:-3.5}

HOLO_PDB=${POSARGS[0]:-""}
LIGAND_RESNAME=${POSARGS[1]:-""}
BS_RESIDUES=${POSARGS[2]:-""}
NCORES=${POSARGS[3]:-32}
OUTDIR=${POSARGS[4]:-"redocking"}

if [[ -z "${HOLO_PDB}" ]] || [[ -z "${LIGAND_RESNAME}" ]] || [[ -z "${BS_RESIDUES}" ]]; then
    echo "Usage: $0 <holo_complex.pdb> <ligand_resname> <bs_residues> [ncores] [outdir] [--rigid|--flex] [--run]"
    echo ""
    echo "Arguments:"
    echo "  holo_complex.pdb  - PDB file with protein (ATOM) + ligand (HETATM)"
    echo "  ligand_resname    - 3-letter residue name of the ligand (e.g. H7M)"
    echo "  bs_residues       - Space-separated binding site residue numbers (quoted)"
    echo "  ncores            - Number of CPU cores (default: 32)"
    echo "  outdir            - Output directory (default: redocking)"
    echo "  --rigid | --flex  - Ligand flexibility in flexref/emref (default: --rigid)"
    echo "  --run             - Also launch haddock3 on the generated config"
    echo ""
    echo "Example:"
    echo "  $0 6dpt_A.pdb H7M \"64 67 120 152 211 316 318\" 16"
    exit 1
fi

if [[ ! -f "${HOLO_PDB}" ]]; then
    echo "ERROR: ${HOLO_PDB} not found"
    exit 1
fi

# Absolute path: everything below runs inside ${OUTDIR}
HOLO_ABS=$(cd "$(dirname "${HOLO_PDB}")" && pwd)/$(basename "${HOLO_PDB}")

echo "=============================================="
echo "HADDOCK3 REDOCKING SETUP - Binding Site Residues"
echo "=============================================="
echo "Input:       ${HOLO_PDB}"
echo "Ligand:      ${LIGAND_RESNAME}"
echo "BS residues: ${BS_RESIDUES}"
echo "Cores:       ${NCORES}"
echo "Output dir:  ${OUTDIR}/"
if [[ "${MODE}" == "flex" ]]; then
    echo "Flexibility: FLEXIBLE (nseg2 auto, torsions refined in flexref/emref)"
else
    echo "Flexibility: RIGID (nseg2 = 0, input conformer kept)"
fi
echo "=============================================="
echo ""

mkdir -p "${OUTDIR}/prodrg"
cp "${HOLO_ABS}" "${OUTDIR}/holo_input.pdb"
cd "${OUTDIR}"

# ===========================================
# STAGE 1: EXTRACT LIGAND + PRODRG TOPOLOGY
# ===========================================
echo "==== STAGE 1: EXTRACT LIGAND ===="

echo "  [1/3] Extracting ligand ${LIGAND_RESNAME} from complex..."
# Column-based selection (cols 18-20 = resName). Keeps the first copy only and
# drops alternate locations other than blank/A.
awk -v rn="${LIGAND_RESNAME}" '
    substr($0,1,6) == "HETATM" {
        res = substr($0,18,3); gsub(/ /, "", res)
        if (res != rn) next
        alt = substr($0,17,1)
        if (alt != " " && alt != "A") next
        key = substr($0,22,1) substr($0,23,5)      # chain + resSeq + iCode
        if (first == "") first = key
        if (key != first) { skipped++; next }
        print substr($0,1,16) " " substr($0,18)    # blank the altLoc indicator
    }
    END { if (skipped) printf("        NOTE: skipped %d atoms from additional %s copies\n", skipped, rn) > "/dev/stderr" }
' holo_input.pdb > prodrg/ligand_ref_ligand_raw.pdb

NATOMS=$(wc -l < prodrg/ligand_ref_ligand_raw.pdb | tr -d ' ')
if [[ "${NATOMS}" -eq 0 ]]; then
    echo "ERROR: No HETATM records found for residue ${LIGAND_RESNAME}"
    exit 1
fi
echo "        Extracted: ${NATOMS} atoms"

echo "  [2/3] Canonicalizing ligand for PRODRG..."
python3 << 'PYEOF'
from rdkit import Chem

def format_pdb_line(atom_count, atom_name, element, coords, chain='B', resname='UNK', resnum=1):
    x, y, z = coords
    if len(element) == 1:
        atom_name_field = f" {atom_name:<3s}"
    else:
        atom_name_field = f"{atom_name:<4s}"

    line = (
        f"ATOM  {atom_count:5d} {atom_name_field} {resname:>3s} {chain}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  "
    )
    return line

mol = Chem.MolFromPDBFile('prodrg/ligand_ref_ligand_raw.pdb', removeHs=False)
if mol is None:
    print("ERROR: Could not parse ligand")
    exit(1)

# Canonicalize atom order so the PRODRG naming is reproducible
_ = Chem.MolToSmiles(mol)
mol_order = list(mol.GetPropsAsDict(True, True).get("_smilesAtomOutputOrder", range(mol.GetNumAtoms())))
mol_canonical = Chem.RenumberAtoms(mol, mol_order)

conf = mol_canonical.GetConformer()
lines = []
atom_type_counts = {}

for atom_idx in range(mol_canonical.GetNumAtoms()):
    atom = mol_canonical.GetAtomWithIdx(atom_idx)
    pos = conf.GetAtomPosition(atom_idx)

    element = atom.GetSymbol().upper()
    atom_type_counts[element] = atom_type_counts.get(element, 0) + 1
    atom_name = f"{element}{atom_type_counts[element]}"

    lines.append(format_pdb_line(
        atom_count=atom_idx + 1,
        atom_name=atom_name,
        element=element,
        coords=(pos.x, pos.y, pos.z),
    ))

with open('prodrg/ligand_for_prodrg.pdb', 'w') as f:
    f.write('\n'.join(lines) + '\nEND\n')

print(f"        Canonicalized: {mol_canonical.GetNumAtoms()} atoms")
print(f"        Written: prodrg/ligand_for_prodrg.pdb")
PYEOF

echo "  [3/3] Calling PRODRG API..."
if [[ -s prodrg/prodrg_output.json ]]; then
    echo "        Reusing existing prodrg/prodrg_output.json (delete it to force a new call)"
else
    python3 << 'PYEOF'
import requests
import sys
import time

url = 'https://wenmr.science.uu.nl/prodrg'
max_retries = 3
timeout = 120

for attempt in range(1, max_retries + 1):
    try:
        print(f"        Attempt {attempt}/{max_retries}...", flush=True)
        start = time.time()

        with open('prodrg/ligand_for_prodrg.pdb', 'rb') as f:
            files = {'input': ('ligand_for_prodrg.pdb', f, 'text/plain')}
            data = {'passkey': 'paralig4haddock'}
            r = requests.post(url, files=files, data=data, timeout=timeout)

        elapsed = time.time() - start
        print(f"        Response received in {elapsed:.1f}s", flush=True)

        if r.status_code == 200 and len(r.text) > 100:
            with open('prodrg/prodrg_output.json', 'w') as out:
                out.write(r.text)
            print(f"        PRODRG output saved ({len(r.text)} bytes)", flush=True)
            sys.exit(0)
        else:
            print(f"        Invalid response (status: {r.status_code}, size: {len(r.text)})", flush=True)
    except requests.exceptions.Timeout:
        print(f"        Timeout after {timeout}s", flush=True)
    except Exception as e:
        print(f"        Error: {e}", flush=True)

    if attempt < max_retries:
        print(f"        Retrying in 5 seconds...", flush=True)
        time.sleep(5)

print("ERROR: PRODRG API failed after all retries")
sys.exit(1)
PYEOF
fi

echo "  STAGE 1 COMPLETE"

# ===========================================
# STAGE 2: EXTRACT PRODRG FILES
# ===========================================
echo ""
echo "==== STAGE 2: EXTRACT PRODRG FILES ===="

python3 << 'PYEOF'
import json

with open('prodrg/prodrg_output.json', 'r') as f:
    data = json.load(f)

files = data.get('files', {})
print(f"  Found {len(files)} files in PRODRG response")

file_keys = list(files.keys())

param_file = [k for k in file_keys if k.endswith('.par')]
if not param_file:
    print("ERROR: No parameter file (.par) found")
    exit(1)
print(f"  Found parameters: {param_file[0]}")
with open('ligand.param', 'w') as f:
    f.write(files[param_file[0]])

top_file = [k for k in file_keys if k.endswith('.top')]
if not top_file:
    print("ERROR: No topology file (.top) found")
    exit(1)
print(f"  Found topology: {top_file[0]}")
with open('ligand.top', 'w') as f:
    f.write(files[top_file[0]])

# Heavy-atom PDB (_noh.pdb) preferred
pdb_file = [k for k in file_keys if '_noh.pdb' in k]
if not pdb_file:
    pdb_file = [k for k in file_keys if k.endswith('.pdb') and '_poh' not in k]
if not pdb_file:
    print("ERROR: No PDB file found")
    exit(1)
print(f"  Found PDB: {pdb_file[0]}")
with open('prodrg/ligand_prodrg.pdb', 'w') as f:
    f.write(files[pdb_file[0]])

print("  Written: ligand.top, ligand.param, prodrg/ligand_prodrg.pdb")
PYEOF

echo "  STAGE 2 COMPLETE"

# ===========================================
# STAGE 3: BUILD THE HADDOCK LIGAND
# PRODRG atom names + crystal coordinates
# ===========================================
echo ""
echo "==== STAGE 3: PREPARE LIGAND FOR REDOCKING ===="
echo "  Using CRYSTAL POSE (single conformer)"

python3 << 'PYEOF'
import re
import numpy as np
from collections import defaultdict
from scipy.optimize import linear_sum_assignment

def extract_topology_atom_names(top_file):
    """Extract heavy atom names from the topology file, in order."""
    with open(top_file, 'r') as f:
        content = f.read()

    residue_match = re.search(r'Residue\s+UNK.*?(?=\nResidue|\Z)', content, re.DOTALL)
    if not residue_match:
        return []

    atom_names = re.findall(r'ATOM\s+(\S+)\s+TYPE=', residue_match.group(0))
    return [name for name in atom_names if not name.startswith('H')]

def get_element_from_line(line):
    """Element from PDB columns 77-78, falling back to the atom name."""
    if len(line) >= 78:
        elem = line[76:78].strip()
        if elem and elem[0].isalpha():
            return elem.upper()
    atname = line[12:16].strip()
    name = atname.lstrip("0123456789'\"")
    if len(name) >= 2 and name[:2].upper() in ['CL', 'BR', 'FE', 'ZN', 'MG', 'CA', 'NA']:
        return name[:2].upper()
    return name[0].upper() if name else 'C'

def is_hydrogen_from_line(line):
    return get_element_from_line(line) == 'H'

def format_pdb_line(atom_count, atom_name, element, coords, chain='B', resname='UNK', resnum=1):
    x, y, z = coords
    if len(element) == 1:
        atom_name_field = f" {atom_name:<3s}"
    else:
        atom_name_field = f"{atom_name:<4s}"

    line = (
        f"ATOM  {atom_count:5d} {atom_name_field} {resname:>3s} {chain}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  "
    )
    return line

def read_heavy_atoms(path):
    atoms = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith(('HETATM', 'ATOM')) and not is_hydrogen_from_line(line):
                atoms.append({
                    'name': line[12:16].strip(),
                    'coords': np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                    'element': get_element_from_line(line),
                })
    return atoms

topology_atoms = extract_topology_atom_names('ligand.top')
print(f"  Topology heavy atoms: {len(topology_atoms)}")

prodrg_atoms = read_heavy_atoms('prodrg/ligand_prodrg.pdb')
print(f"  PRODRG heavy atoms:   {len(prodrg_atoms)}")

crystal_atoms = read_heavy_atoms('prodrg/ligand_ref_ligand_raw.pdb')
print(f"  Crystal heavy atoms:  {len(crystal_atoms)}")

# =====================================================
# Match PRODRG atoms to crystal atoms per element with
# the Hungarian algorithm (optimal global assignment).
# =====================================================
prodrg_coords = np.array([a['coords'] for a in prodrg_atoms])
crystal_coords = np.array([a['coords'] for a in crystal_atoms])

prodrg_centered = prodrg_coords - prodrg_coords.mean(axis=0)
crystal_centered = crystal_coords - crystal_coords.mean(axis=0)

prodrg_by_element = defaultdict(list)
crystal_by_element = defaultdict(list)
for i, a in enumerate(prodrg_atoms):
    prodrg_by_element[a['element']].append(i)
for j, a in enumerate(crystal_atoms):
    crystal_by_element[a['element']].append(j)

prodrg_to_crystal = {}
for element, p_indices in prodrg_by_element.items():
    c_indices = crystal_by_element.get(element, [])

    if len(c_indices) < len(p_indices):
        print(f"  ERROR: Not enough crystal {element} atoms ({len(c_indices)}) to match PRODRG ({len(p_indices)})")
        exit(1)

    cost = np.zeros((len(p_indices), len(c_indices)))
    for ii, pi in enumerate(p_indices):
        for jj, cj in enumerate(c_indices):
            cost[ii, jj] = np.linalg.norm(prodrg_centered[pi] - crystal_centered[cj])

    row_ind, col_ind = linear_sum_assignment(cost)
    for ii, jj in zip(row_ind, col_ind):
        prodrg_to_crystal[p_indices[ii]] = c_indices[jj]

    max_dist = max(cost[ii, jj] for ii, jj in zip(row_ind, col_ind))
    print(f"  Matched {len(p_indices):2d} {element} atoms (max dist: {max_dist:.3f} A)")

matched_crystal_coords = [crystal_atoms[prodrg_to_crystal[i]]['coords'] for i in range(len(prodrg_atoms))]
print(f"  Successfully matched {len(matched_crystal_coords)} atoms")

lines = [
    format_pdb_line(i + 1, atom['name'], atom['element'], coords)
    for i, (atom, coords) in enumerate(zip(prodrg_atoms, matched_crystal_coords))
]

with open('ligand_haddock.pdb', 'w') as f:
    f.write('\n'.join(lines) + '\nEND\n')

print("  Written: ligand_haddock.pdb (single crystal pose, chain B, resid 1)")
PYEOF

echo "  STAGE 3 COMPLETE"

# ===========================================
# STAGE 4: PREPARE RECEPTOR + CAPRI REFERENCE
# ===========================================
echo ""
echo "==== STAGE 4: PREPARE RECEPTOR ===="

echo "  [1/3] Extracting receptor (ATOM records only)..."
awk '
    substr($0,1,4) == "ATOM" {
        alt = substr($0,17,1)
        if (alt != " " && alt != "A") next
        print substr($0,1,16) " " substr($0,18)    # blank the altLoc indicator
    }
' holo_input.pdb > receptor_tmp.pdb

if command -v pdb_chain >/dev/null 2>&1; then
    pdb_chain -A receptor_tmp.pdb > receptor_body.pdb 2>/dev/null || \
        sed 's/^\(.\{21\}\)./\1A/' receptor_tmp.pdb > receptor_body.pdb
else
    sed 's/^\(.\{21\}\)./\1A/' receptor_tmp.pdb > receptor_body.pdb
fi
rm -f receptor_tmp.pdb

cp receptor_body.pdb receptor.pdb
printf 'TER\nEND\n' >> receptor.pdb

NRECP=$(grep -c "^ATOM" receptor.pdb)
echo "        Receptor atoms: ${NRECP} (chain A)"

echo "  [2/3] Validating BS residues against the receptor..."
MISSING=""
for R in ${BS_RESIDUES}; do
    if ! awk -v r="${R}" 'substr($0,1,4)=="ATOM" && substr($0,23,4)+0==r {found=1; exit} END{exit !found}' receptor_body.pdb; then
        MISSING="${MISSING} ${R}"
    fi
done
if [[ -n "${MISSING}" ]]; then
    echo "ERROR: BS residue(s) not present in the receptor:${MISSING}"
    echo "       Check that the numbering in the BS list matches ${HOLO_PDB}."
    exit 1
fi
echo "        All BS residues found"

echo "  [3/3] Creating CAPRI reference (receptor chain A + ligand chain B)..."
# Single END at the very end: appending after an END line hides the ligand from
# strict PDB parsers.
cat receptor_body.pdb > reference.pdb
printf 'TER\n' >> reference.pdb
grep "^ATOM" ligand_haddock.pdb >> reference.pdb
printf 'TER\nEND\n' >> reference.pdb
rm -f receptor_body.pdb
echo "        Created: reference.pdb"

echo "  STAGE 4 COMPLETE"

# ===========================================
# STAGE 5: GENERATE AIR RESTRAINTS
# BS residues = ACTIVE, ligand = ACTIVE
# ===========================================
echo ""
echo "==== STAGE 5: GENERATE AIR RESTRAINTS ===="
echo "  Protocol: HADDOCK protein-ligand binding site approach"

# Machine-readable BS list consumed by the cross-docking runners
echo "${BS_RESIDUES}" | tr ' ' '\n' | grep -v '^$' | sort -un > bs_residues.txt
echo "  Written: bs_residues.txt ($(wc -l < bs_residues.txt | tr -d ' ') residues)"

python3 << PYEOF
bs_list = sorted({int(r) for r in "${BS_RESIDUES}".split()})
lig_list = [1]

def air_block(resid, segid, targets, target_segid):
    """One AIR: resid/segid restrained to the union of targets."""
    parts = [f"(resid {t} and segid {target_segid})" for t in targets]
    sel = '\n'.join(
        ("        " if i == 0 else "     or ") + p for i, p in enumerate(parts)
    )
    return f"assign (resid {resid} and segid {segid})\n(\n{sel}\n) 2.0 2.0 0.0"

restraints = [air_block(r, 'A', lig_list, 'B') for r in bs_list]
restraints += [air_block(r, 'B', bs_list, 'A') for r in lig_list]

with open('ambig_active.tbl', 'w') as f:
    f.write('! HADDOCK AIR restraints - binding site approach\n')
    f.write('! BS residues (segid A) = ACTIVE, ligand (segid B, resid 1) = ACTIVE\n')
    f.write('! Draws the ligand into the binding pocket\n!\n')
    f.write('\n\n'.join(restraints) + '\n')

print(f"  Written: ambig_active.tbl ({len(restraints)} restraints, {len(bs_list)} BS residues)")
PYEOF

echo "  STAGE 5 COMPLETE"

# ===========================================
# STAGE 6: GENERATE HADDOCK3 CONFIG
# ===========================================
echo ""
echo "==== STAGE 6: GENERATE HADDOCK3 CONFIG ===="

# Ligand flexibility blocks: nseg2 = 0 pins the ligand's conformer, while
# leaving nseg2 unset keeps HADDOCK's default (-1 = segments auto-defined from
# contacts), which lets flexref/emref refine the ligand torsions.
if [[ "${MODE}" == "flex" ]]; then
    FLEX_NOTE="FLEXIBLE (nseg2 auto, ligand torsions refined)"
    FLEXREF_LIG=""
    EMREF_LIG=""
else
    FLEX_NOTE="RIGID (nseg2 = 0, input conformer kept)"
    FLEXREF_LIG="nseg2 = 0"
    EMREF_LIG="nseg2 = 0"
fi

echo "${MODE}" > docking_mode.txt

cat > "haddock3_redocking_bs_${MODE}.cfg" << CFGEOF
# HADDOCK3 Protein-Ligand REDOCKING
# Protocol: BINDING SITE RESIDUES approach
# Generated: $(date)
#
# Mode: REDOCKING (holo receptor, single crystal ligand conformer).
# Reference: https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-binding-sites/
#
# Restraints: ambig_active.tbl (BS = active, ligand = active) is used in every
# docking stage. Full vdW (w_vdw = 1.0) is kept because the holo pocket is open;
# the cross-docking runners drop it to 0.0 for closed/homology-model pockets.
#
# Ligand flexibility: ${FLEX_NOTE}
#
# irmsd_cutoff is the distance (A) used to define which residues count as
# interface, not a quality threshold. caprieval reports i-RMSD, l-RMSD,
# i-l-RMSD, fnat and DockQ regardless.

run_dir = "run_redocking_bs_${MODE}"
ncores = ${NCORES}
mode = "local"

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
# ===========================================
[rigidbody]
sampling = 250
tolerance = 5
ambig_fname = "ambig_active.tbl"
w_vdw = 1.0
inter_rigid = 0.001
mol_fix_origin_1 = true
randremoval = false
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Keep the top 250 -> 50 rigid-body models
# ===========================================
[seletop]
select = 50

# ===========================================
# Flexible refinement -- ligand: ${FLEX_NOTE}
# ===========================================
[flexref]
tolerance = 5
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
mdsteps_rigid = 0
mdsteps_cool1 = 0
${FLEXREF_LIG}
randremoval = false

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = ${IRMSD_CUTOFF}

# ===========================================
# Energy-minimization refinement -- ligand: ${FLEX_NOTE}
# ===========================================
[emref]
tolerance = 5
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
${EMREF_LIG}
randremoval = false

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

echo "  Written: haddock3_redocking_bs_${MODE}.cfg (ligand: ${FLEX_NOTE})"
echo "  STAGE 6 COMPLETE"

# ===========================================
# OPTIONAL RUN
# ===========================================
if [[ "${RUN_HADDOCK}" -eq 1 ]]; then
    echo ""
    echo "==== RUNNING HADDOCK3 (REDOCKING) ===="
    STARTTIME=$(date +%s)
    if haddock3 "haddock3_redocking_bs_${MODE}.cfg" > "haddock3_redocking_${MODE}.log" 2>&1; then
        echo "  DONE ($(( $(date +%s) - STARTTIME ))s)"
    else
        echo "  FAILED ($(( $(date +%s) - STARTTIME ))s) -- check ${OUTDIR}/haddock3_redocking_${MODE}.log"
        exit 1
    fi
fi

# ===========================================
# SUMMARY
# ===========================================
echo ""
echo "=============================================="
echo "SETUP COMPLETE -> ${OUTDIR}/"
echo "=============================================="
echo ""
echo "Shared inputs (consumed by the cross-docking runners):"
echo "  ligand_haddock.pdb        Ligand, crystal pose, chain B resid 1"
echo "  ligand.top                PRODRG topology"
echo "  ligand.param              PRODRG parameters"
echo "  ambig_active.tbl          AIRs (BS = active, ligand = active)"
echo "  reference.pdb             CAPRI reference complex"
echo "  bs_residues.txt           Binding site residue numbers"
echo ""
echo "Redocking-only files:"
echo "  receptor.pdb              Holo receptor, chain A"
echo "  holo_input.pdb            Copy of the input complex"
echo "  haddock3_redocking_bs_${MODE}.cfg   HADDOCK3 config (${MODE} ligand)"
echo "  docking_mode.txt          Mode picked up by the cross-docking runners"
echo "  prodrg/                   Raw PRODRG request/response artefacts"
echo ""
if [[ "${RUN_HADDOCK}" -eq 0 ]]; then
    echo "To run the redocking:"
    echo "  (cd ${OUTDIR} && haddock3 haddock3_redocking_bs_${MODE}.cfg)"
    echo ""
fi
echo "To cross-dock an ensemble of receptor models:"
echo "  ./haddock_run_ensemble_crossdocking_bs.sh <models_dir> ${OUTDIR} ${NCORES} --${MODE}"
echo "  ./haddock_run_crossdocking_bs.sh          <models_dir> ${OUTDIR} ${NCORES} --${MODE}"
echo "=============================================="
