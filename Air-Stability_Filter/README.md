
A list of novel compounds should be placed in a file named ```Candidates```. An example of this file is given in the current repo.

To filter this list of compounds by their stability in air: 

```
python air_stable.py --energ_allowance=50
```

Where ```energ_allowance``` is used to specify the maximum driving force at which each compound reacts with a gaseous species in air. The current version of the script accounts for reactions with O2, CO2, and H2O. If not specified at runtime, the energy allowance defaults to 50 meV/atom.

After running this script, it will output a list of compounds that are sufficiently air-stable, in addition to their associated reaction energies with each gaseous species.
