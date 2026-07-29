#!/bin/bash
#
# Example: AmpC beta-lactamase (PDB 6DPT, chain A) with fragment H7M.
# Binding site = experimentally described BS residues.
#
# Produces redocking/ with everything the cross-docking runners need.
# Add --run at the end to also launch the redocking itself.

./haddock_run_redocking_bs.sh 6dpt_A.pdb H7M "64 67 120 152 211 316 318" 16
