# HADDOCK3 protein-ligand binding-site cross-docking workflow

Three scripts to dock a ligand into an ensemble of receptor conformations
(homology models, MD snapshots, apo structures) with
**[HADDOCK3](https://www.bonvinlab.org/haddock3/)**, using the
[protein-ligand binding-site protocol](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-binding-sites/),
with the ligand either **rigid** or **flexible**.

Everything starts from one reference redocking run, which builds a shared
`redocking/` directory (ligand topology/parameters, AIR restraints, CAPRI
reference). The cross-docking runners then reuse those inputs against any set of
receptor models.

This repository is HADDOCK only: docking plus HADDOCK's own `caprieval`
analysis. No external post-processing.

---

## Requirements

- **HADDOCK3** -> this also provides `numpy`, `scipy` and `pdb-tools`, which are
  hard HADDOCK3 dependencies, so they need no separate install
- **`rdkit` and `requests`** -> required by the redocking setup only, for ligand
  canonicalization and the PRODRG API call. Neither is a HADDOCK3 dependency,
  so install them into the same environment: `pip install rdkit requests`
- Network access to the PRODRG server (`wenmr.science.uu.nl`) for ligand topology

The two cross-docking runners need nothing beyond HADDOCK3 and the standard
library.

Nothing is hard-coded; every script takes its input paths as arguments.

---

## 1. Redocking setup (always first)

```bash
./haddock_run_redocking_bs.sh <holo_complex.pdb> <ligand_resname> "<bs_residues>" [ncores] [outdir] [--rigid|--flex] [--irmsd-cutoff X] [--run]
```

Example, using the included AmpC beta-lactamase structure:

```bash
./haddock_run_redocking_bs.sh 6dpt_A.pdb H7M "64 67 120 152 211 316 318" 16
```

What it does:

| Stage | Action |
|---|---|
| 1 | Extracts the ligand from the complex, canonicalises atom order (RDKit), submits it to PRODRG |
| 2 | Pulls `ligand.top`, `ligand.param` and the heavy-atom PDB out of the PRODRG response |
| 3 | Maps PRODRG atom names onto the **crystal coordinates** (per-element Hungarian assignment) → `ligand_haddock.pdb` |
| 4 | Builds `receptor.pdb` (chain A), validates the BS numbering against it, builds the CAPRI `reference.pdb` |
| 5 | Writes the AIR restraints `ambig_active.tbl` (BS = active, ligand = active) and `bs_residues.txt` |
| 6 | Writes `haddock3_redocking_bs_<mode>.cfg` and `docking_mode.txt` |

Everything lands in `redocking/` (override with the 5th argument):

```
redocking/
├── ligand_haddock.pdb          # shared with cross-docking
├── ligand.top                  # shared
├── ligand.param                # shared
├── ambig_active.tbl            # shared
├── reference.pdb               # shared
├── bs_residues.txt             # shared
├── receptor.pdb                # redocking only
├── holo_input.pdb              # copy of the input complex
├── docking_mode.txt            # rigid | flex, read by the cross-docking runners
├── haddock3_redocking_bs_<mode>.cfg
└── prodrg/                     # raw PRODRG request/response
```

The redocking itself is **not** run by default — the cross-docking runners only
need the files above. To run it as a protocol sanity check:

```bash
./haddock_run_redocking_bs.sh 6dpt_A.pdb H7M "64 67 120 152 211 316 318" 16 redocking --run
```

or later, by hand:

```bash
cd redocking && haddock3 haddock3_redocking_bs_rigid.cfg
```

Re-running the setup reuses an existing `prodrg/prodrg_output.json` instead of
hitting the PRODRG server again; delete that file to force a fresh call.

---

## 2. Cross-docking

Both runners take a directory of receptor `.pdb` models plus the `redocking/`
directory, and launch HADDOCK3 themselves.

### One run per model

```bash
./haddock_run_crossdocking_bs.sh <models_dir> redocking 32 [--rigid|--flex] [--centered] [--irmsd-cutoff X]
```

Pipeline: `rigidbody → seletop(50) → flexref → emref → seletop(5)`, one run per
model, no clustering.

Output: `crossdocking_runs_bs_<mode>/<model>/run_crossdock_bs_<mode>/`

### All models as a single ensemble (preferred)

```bash
./haddock_run_ensemble_crossdocking_bs.sh <models_dir> redocking 32 [--rigid|--flex] [--centered] [--irmsd-cutoff X] [--fcc] [--min-pop N] [--clust-cutoff X] [--no-emref]
```

Concatenates all models into one multi-MODEL receptor, scales rigid-body
sampling by the number of conformers, and runs
`rigidbody → seletop → flexref → emref → ilrmsdmatrix → clustrmsd → seletopclusts`.

Output: `ensemble_crossdocking_bs_<mode>/run_ensemble_crossdock_bs_<mode>/`
Rankings: `.../analysis/*_caprieval/capri_ss.tsv` (poses) and `capri_clt.tsv`
(clusters).

Both runners refine identically through `emref`, so their HADDOCK scores are on
the same scale and can be compared directly.

### Sampling and selection

Sampling scales with the number of conformers `N`; the funnel ratio matches
HADDOCK's own protein-ligand example (1000 → 200, i.e. 5:1):

| | one run per model | ensemble | HADDOCK example |
|---|---|---|---|
| `rigidbody sampling` | 250 | 250 × N | 1000 |
| `seletop` → `flexref` | 50 | 50 × N | 200 |
| ratio | 5:1 | 5:1 | 5:1 |

Everything selected is refined; there is no second `seletop` before `emref`.
Note the absolute per-receptor sampling is 4× below the example (250 vs 1000),
which is a deliberate cost trade when `N` is large. Raise the multipliers at the
top of the script if enrichment is poor.

---

## 3. Centered placement -> `--centered`

Both cross-docking runners accept `--centered`, which uses the HADDOCK3 feature
added on 2026-04-15 ([Issue #1491](https://github.com/haddocking/haddock3/issues/1491),
example configs `examples/docking-protein-ligand/docking-protein-ligand-centered-*.cfg`):

```bash
./haddock_run_ensemble_crossdocking_bs.sh <models_dir> redocking 32 [--rigid|--flex] --centered
./haddock_run_crossdocking_bs.sh          <models_dir> redocking 32 --centered
```

By default `[rigidbody]` runs `separate.cns`, which re-orients each non-fixed
molecule along its principal axes and pushes the molecules apart before docking
— so the ligand always starts outside the protein and the restraints have to
drag it back. `separate = false` skips that step entirely: the ligand keeps the
coordinates it is given, and only `randrot` (random rotation about its own
centre of mass) is applied, so orientation is still sampled but position is not
thrown away. A second effect matters here too: with `separate = true` the
initial rotational rigid-body minimisation runs with `flag excl vdw elec end`,
whereas with `separate = false` those terms stay on, which is what stops a
pre-placed ligand being minimised through the protein.

What the flag changes:

| | default | `--centered` |
|---|---|---|
| ligand start | scattered away from the receptor | left in the pocket |
| `separate` | `true` (HADDOCK default) | `false` |
| `cmrest` / `cmtight` | `true` | dropped — they only exist to pull a scattered ligand back |

Everything else (AIRs, the ligand file, `w_vdw = 0.0`, `tolerance = 20`,
`mol_fix_origin_1`, ligand flexibility, semi-flexible segments, the whole
downstream pipeline) is unchanged, so the two modes are directly comparable on the same
models.

### No ligand preparation is needed

`redocking/ligand_haddock.pdb` already holds the crystal pose in the holo
structure's frame — it *is* the centered ligand. As long as your receptor models
share that frame, which they must anyway for the downstream RMSD analysis, the
ligand is already sitting in every model's binding site and `separate = false`
simply tells HADDOCK to leave it there. No repositioning step, no extra tooling.

The only failure mode is models built or saved in a different frame, in which
case the ligand would start somewhere in solvent. Both runners guard against
that before launching: they compare the ligand's centre of mass with the
centroid of the model's BS residues and abort if it exceeds `MAX_BS_DIST`
(10 Å, set at the top of each script). Superimpose your models onto the holo
structure and re-run, or drop `--centered`.

Requires a HADDOCK3 with the `separate` parameter in `[rigidbody]`; check with:

```bash
haddock3-cfg -m rigidbody | grep -A3 "^separate"
```

Redocking deliberately does **not** use this: the ligand there is already at its
crystal pose, so pre-placing it would hand HADDOCK the answer and destroy the
redocking's value as a control.

---

## Clustering and refinement -> ensemble runner

The ensemble runner clusters with **iL-RMSD** (`[ilrmsdmatrix]` + `[clustrmsd]`),
matching HADDOCK3's own `examples/docking-protein-ligand/*-full.cfg`.

```bash
./haddock_run_ensemble_crossdocking_bs.sh models redocking 32 \
    [--fcc] [--min-pop N] [--clust-cutoff X] [--no-emref]
```

| flag | default | effect |
|---|---|---|
| `--min-pop N` | `4` | minimum cluster size — **see the warning below** |
| `--clust-cutoff X` | `1.5` Å (`0.6` with `--fcc`) | clustering cutoff |
| `--fcc` | off | cluster on fraction of common contacts instead of iL-RMSD |
| `--no-emref` | off | drop `[emref]`, cluster straight after `flexref` |

### Why iL-RMSD and not FCC

FCC counts residue–residue contacts. For a small ligand pre-placed in a single
pocket, every pose contacts essentially the same BS residues regardless of its
orientation, so FCC saturates near 1.0 and the run collapses into one giant
cluster holding most of the poses. iL-RMSD superposes on the receptor interface
and measures the ligand, so genuinely different binding modes separate.

`--fcc` is kept for ligands large enough that their contact fingerprint actually
varies between poses. If you use it, check the cluster populations in
`capri_clt.tsv`: one cluster holding most of the models means the cutoff has
saturated and the clustering is telling you nothing.

### `min_population` is not only a size filter

Both `clustfcc` and `clustrmsd` pass `min_population` to `rank_clusters()` as
**the number of top models averaged into the cluster score**
(`libs/libclust.py`). So `--min-pop 1` ranks every cluster by its single best
pose, which lets a 2-member cluster outrank a well-populated one on one lucky
model, and makes the `cluster_rank` column disagree with the `score` column in
`capri_clt.tsv` (caprieval computes that with `clt_threshold`, default 4).

Keep it at `4` unless you have a specific reason. The cost is that models in
clusters smaller than `N` are **not carried past the clustering step** — they
survive only in the pre-clustering `caprieval` tables.

### emref

`emref` runs by default in both cross-docking runners. Alex Bonvin advises
against it for small molecules, but that advice assumes a full `flexref` with MD
in all four stages; here `flexref` is deliberately gutted
(`mdsteps_rigid = 0`, `mdsteps_cool1 = 0`, plus `mdsteps_cool3 = 0` and
`nseg2 = 0` in rigid mode), so nothing else relieves the `w_vdw = 0` pocket
clashes. Without `emref` the residual clashes leave HADDOCK scores large and
positive, and because the clash penalty swamps the interaction terms the
score-based ranking degrades: near-native poses sink down the list.

Use `--no-emref` only together with `--flex` or restored `mdsteps_*` defaults.
To check whether it is earning its keep on your system, compare the `score`
column of the `caprieval` steps before and after it: if the pre-`emref` range
still runs to large positive values, `emref` is doing necessary work.

---

## Ligand flexibility -> `--rigid` / `--flex`

All three scripts take the same pair of flags. `--rigid` is the default, and
reproduces the previous behaviour exactly.

| | `--rigid` (default) | `--flex` |
|---|---|---|
| `nseg2` in `flexref` / `emref` | `0` — no semi-flexible segments, ligand conformer pinned | unset → HADDOCK default `-1`, segments auto-defined from contacts |
| `mdsteps_cool3` in `flexref` | `0` (cross-docking only) | unset → default `1000`, so the torsion-angle cooling phase actually runs |
| comparable to | rigid-ligand AutoDock / Vina | the standard HADDOCK protein-ligand protocol |

`mdsteps_rigid = 0` and `mdsteps_cool1 = 0` are kept in both modes — those are
the *rigid-body* hot and first-cooling phases, and dropping them matches the
upstream protein-ligand example. Ligand flexibility lives in `cool2`/`cool3`,
which is why `--flex` restores `mdsteps_cool3`.

### Keeping the stages consistent

The redocking setup records the mode in `redocking/docking_mode.txt`. If you run
a cross-docking script **without** a flag, it reads that file and reports where
the mode came from:

```
Ligand:       RIGID (nseg2 = 0, input conformer kept) [redocking/docking_mode.txt]
```

An explicit `--rigid` / `--flex` on the command line always wins.

Nothing collides between the two modes: the shared inputs (`ligand.top`,
`ligand.param`, `ambig_active.tbl`, `reference.pdb`, `ligand_haddock.pdb`,
`bs_residues.txt`) are identical, and only the configs, run directories and
output directories carry a `_rigid` / `_flex` suffix. One `redocking/` directory
therefore serves both, and you can run rigid and flex side by side:

```bash
./haddock_run_redocking_bs.sh 6dpt_A.pdb H7M "64 67 120 152 211 316 318" 16 redocking --rigid
./haddock_run_redocking_bs.sh 6dpt_A.pdb H7M "64 67 120 152 211 316 318" 16 redocking --flex
./haddock_run_ensemble_crossdocking_bs.sh models redocking 32 --rigid
./haddock_run_ensemble_crossdocking_bs.sh models redocking 32 --flex
```

---

## Analysis -> `caprieval` 

Scoring and RMSD come entirely from HADDOCK's `caprieval` module; there is no
external post-processing step. Every `caprieval` block reports **i-RMSD,
l-RMSD, i-l-RMSD, fnat and DockQ** — they are all on by default, so there is
nothing extra to enable.

Results land in `.../analysis/*_caprieval/`:

- `capri_ss.tsv` — one row per pose
- `capri_clt.tsv` — one row per cluster

### Which metric measures the binding mode

The three RMSDs differ in what they superpose on and what they then measure:

| metric | fit on | measure | |
|---|---|---|---|
| **i-RMSD** | BS interface **+ ligand** | the same set | fit set == measured set |
| **i-l-RMSD** | receptor BS interface only | **ligand only** | the binding-mode metric |
| **l-RMSD** | whole receptor | ligand only | classic docking RMSD |

**Report i-l-RMSD (and l-RMSD) as the binding-mode metrics.** i-RMSD averages
receptor backbone into the number, and across an ensemble of different homology
models that partly measures model spread rather than pose accuracy.

### `irmsd_cutoff` is not a quality threshold

`irmsd_cutoff` (`--irmsd-cutoff`, default `3.5` in all three scripts) is the
**distance in Å used to decide which residues count as interface**, based on
intermolecular contacts. HADDOCK's own default is `10.0`. It does *not* filter or
accept poses by i-RMSD, and there is no `lrmsd_cutoff` or `ilrmsd_cutoff`
parameter — the only other cutoff is `fnat_cutoff` (contact definition, default
`5.0`).

It drives **both** metrics: `caprieval` passes `irmsd_cutoff` to `calc_irmsd`
*and* to `calc_ilrmsd`, so it sets the i-RMSD atom set and the i-l-RMSD *fit
basis* at the same time.

That matters because the two metrics respond to it in opposite directions. The
i-RMSD atom set is *receptor interface + ligand*; the ligand contributes a fixed
number of atoms while the receptor interface grows with the cutoff, so the
ligand's share of the i-RMSD signal falls as the cutoff rises. For a
small molecule the effect is large: at a tight cutoff the ligand can be roughly
half the atoms, and at HADDOCK's `10.0` only a few percent, at which point
i-RMSD mostly reports how similar the receptor interface backbones are. Across
an ensemble of different receptor conformers that is model spread, not pose
accuracy.

i-l-RMSD is 100 % ligand at any cutoff, so it never suffers that dilution — but
it is *fit* on the receptor interface, and a larger cutoff gives it a broader,
more stable superposition basis. Raising the cutoff therefore improves i-l-RMSD
while degrading i-RMSD, which is the main argument for treating i-l-RMSD as the
primary metric and leaving the cutoff at HADDOCK's default.

To see where your own system sits, count the atoms at each cutoff against your
`redocking/reference.pdb` (receptor = chain A, ligand = chain B):

```bash
python3 -c '
rec, lig = [], []
for l in open("redocking/reference.pdb"):
    if l.startswith(("ATOM","HETATM")) and (l[76:78].strip() or l[12:16].strip()[0]) != "H":
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        (lig if l[21] == "B" else rec).append((l[22:26].strip(), l[12:16].strip(), xyz))
for cut in (3.5, 5.0, 10.0):
    ifres = {r for r, n, x in rec
             if any((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2 <= cut*cut
                    for _, _, y in lig)}
    bb = sum(1 for r, n, x in rec if r in ifres and n in {"C","CA","N","O"})
    print(f"cutoff {cut:5.1f} -> {len(ifres):3d} interface res, "
          f"{bb:4d} rec bb + {len(lig):3d} lig atoms, "
          f"ligand = {100*len(lig)/(bb+len(lig)):5.1f}% of i-RMSD")
'
```

Whichever you pick, keep it identical across redocking and cross-docking or the
tables are not comparable:

```bash
./haddock_run_ensemble_crossdocking_bs.sh models redocking 32 --irmsd-cutoff 10.0
```

Clustering is unaffected — `[ilrmsdmatrix]` has its own
`contact_distance_cutoff` (default `5.0` Å), independent of this parameter.

---

## Protocol notes

- **Restraints.** `ambig_active.tbl` (BS = active, ligand = active) is used in
  every docking stage of every script. Binding-site residues are read from
  `bs_residues.txt`, falling back to parsing the AIR table.
- **Ligand flexibility.** See the `--rigid` / `--flex` section above.
- **Open vs. closed pockets.** Redocking keeps full vdW (`w_vdw = 1.0`) because
  the holo pocket is open. Cross-docking uses `w_vdw = 0.0` and
  `tolerance = 20`, so the ligand can penetrate the closed pockets of
  apo/homology models; semi-flexible segments around the BS residues are
  generated automatically for `flexref`. Getting the ligand to the pocket is
  handled either by `cmrest`/`cmtight` (default) or by pre-placement
  (`--centered`, see above).
- **Clustering.** The ensemble runner clusters by **iL-RMSD**
  (`[ilrmsdmatrix]` + `[clustrmsd]`), as in HADDOCK3's protein-ligand examples.
  `--fcc` switches to fraction of common contacts, which suits larger ligands
  but saturates for a small ligand in a single pocket. See the clustering
  section above. The RMSD metrics themselves always come from `caprieval`.
- **emref.** Energy minimisation after `flexref` relieves the clashes that
  `w_vdw = 0` leaves behind, so HADDOCK scores become meaningful (negative for
  good complexes). EM cannot cross energy barriers, so badly clashed non-native
  poses may stay positive. Run in both cross-docking runners and in redocking,
  with `randremoval = false` throughout so the AIR energy is evaluated on the
  full restraint set at every stage.

---

## Example data

`6dpt_A.pdb` -> chain A of PDB **6DPT** (AmpC beta-lactamase) with the fragment
`H7M` bound. `run-for-redocking_exptBS.sh` is the one-line runner for it.

---

## License

MIT — see [LICENSE](LICENSE).
