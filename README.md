# novel-materials-screening (A-lab)
This repository contains the code/scripts/notebooks/data used to generate and filter novel materials for trial synthesis by the Ceder Group's autonomous solid-state synthesis laboratory ("A-lab"). The provided results were utilized in the experimental study in the corresponding manuscript.

## Notebooks

The `target_screening.ipynb` notebook is a Jupyter notebook containing the code used to identify thermodynamically stable (i.e., "on the hull") candidate materials from the Materials Project. These candidate materials were further filtered and targeted by our A-lab study.

`target_screening_revisions.ipynb` is a condensed version of the above notebook with support for identifying metastable targets, as demonstrated with a <0.050 eV/atom above hull cutoff.

## Scripts
`Air-Stability_Filter/air_stable.py`: this script can be used to filter candidates by those that are air-stable (a current requirement of the A-lab). Please see the README provided in the folder for further details.

## Results

`screened_targets.xlsx`: the output of `target_screening.ipynb`. These candidates were not yet filtered for air stability or checked against the literature.

`screened_targets_50meV.xlsx`: the output of `target_screening_revisions.ipynb`. These candidates were not yet filtered for air stability or checked against the literature.

`final_targets_manually_cleaned_11_11_22.xlsx`: the output of `target_screening.ipynb`, checked manually against the literature. These candidates were not yet filtered for air stability.

## Data
All other provided files (e.g., `icsd_compounds.json.gz`) were used in the filtering steps within `target_screening.ipynb`.



