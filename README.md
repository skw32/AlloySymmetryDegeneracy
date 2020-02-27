## Overview
Workflow used to determine the symmetry degeneracy of alloy supercells with substitutions amongst transition metal sites. This information is used for computing configurational entropy.

## Dependencies

- Python3
- ase
- spglib

## Getting started

Modify `INPUTS` section in `process_dataList.py` and run with `python process_dataList.py`.

The code runs in parallel. During execution, the user will be shown how many processors are available to them to use and asked how many the program should use. It is advisable not to use too many of those you have available or your machine will become sluggish!

Inputs:
- File containing list of directories where unrelaxed POSCAR's can be found (note here they are called 'POSCAR_orig')
- File containing data from previous step of workflow (to be appended by this step)
- Settings such as `threshold` (tolerance spglib will use for assigning space groups) and `scaling` which is used to determine number of random configurations to attempt when searching for equivalent structures (total attempts is total_combination_space*scaling to increase likelihood of sampling most of the possible substitutions). The final count of symmetrically degenerate structures for each input structure is divided by `scaling`.
