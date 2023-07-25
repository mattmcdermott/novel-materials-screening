from pymatgen.analysis import interface_reactions
import pymatgen.analysis.phase_diagram as pd
from pymatgen.core import Composition
from mp_api.client import MPRester
from rxn_help import energetics
import time
import sys


"""
Author: Nathan J. Szymanski
Email: nathan_szymanski@berkeley.edu
Description: this script is used to filter materials
    by their stability in air. It does so by checking for
    reactivity with O2, CO2, and H2O.
"""

mpr = MPRester('You MP API key here')

energ_allowance = 50.0 # meV/atom
for arg in sys.argv:
    if '--energ_allowance' in arg:
        energ_allowance = float(arg.split('=')[1])
energ_allowance /= 1000.0 # Convert to eV/atom

all_cmpds = []
with open('Candidates') as f:
    for line in f.readlines():
        all_cmpds.append(line[:-1])

final_cmpds = {}

for i, cmpd in enumerate(all_cmpds):

    time.sleep(10)

    # Track progress
    print('%s/%s' % (i, len(all_cmpds)))

    # Standardize chemical formula and parse elems
    target_formula = Composition(cmpd).reduced_formula
    elems = [str(el) for el in Composition(cmpd).elements]

    # Add H, C, O
    if 'H' not in elems:
        elems += 'H'
    if 'C' not in elems:
        elems += 'C'
    if 'O' not in elems:
        elems += 'O'

    # Exp. conditions
    T = 800. + 273.
    atmos = 'air'

    # Compute phase diagram
    standard_pd, grand_pd = energetics.get_phase_diagrams(elems, T, atmos)

    """
    Check for stability versus O2
    """
    air_stable = False
    energs = []
    for entry in list(grand_pd.entries):

        # Get composition
        dict = entry.as_dict()
        current_formula = Composition(dict['name']).reduced_formula

        # Check for match
        if current_formula == target_formula:

            # Get energy above hull
            e_above_hull = grand_pd.get_e_above_hull(entry)

            # Allow 50 meV/atom above hull
            if e_above_hull <= energ_allowance:
                air_stable = True
                energs.append(e_above_hull)

    if air_stable:

        # Save oxidation energy
        O_energ = min(energs)

        """
        Check for stability versus CO2
        """
        rxns = interface_reactions.InterfacialReactivity(Composition(cmpd), Composition('CO2'), standard_pd)

        # Allow 50 meV/atom rxn energy (CO2 uptake)
        if rxns.minimum[1] > -energ_allowance:
            CO2_stable = True
            CO2_energ = abs(round(rxns.minimum[1], 3))
        else:
            CO2_stable = False

        if CO2_stable:

            """
            Check for stability versus H2O
            """
            rxns = interface_reactions.InterfacialReactivity(Composition(cmpd), Composition('H2O'), standard_pd)

            # Allow 50 meV/atom rxn energy (H2O uptake)
            if rxns.minimum[1] > -energ_allowance:
                H2O_stable = True
                H2O_energ = abs(round(rxns.minimum[1], 3))
            else:
                H2O_stable = False

            if H2O_stable:

                final_cmpds[target_formula] = {'O': O_energ, 'CO2': CO2_energ, 'H2O': H2O_energ}

for cmpd in final_cmpds.keys():

    print(cmpd, final_cmpds[cmpd]['O'], final_cmpds[cmpd]['CO2'], final_cmpds[cmpd]['H2O'])
