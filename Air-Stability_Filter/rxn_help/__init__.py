from pymatgen.analysis import interface_reactions
import pymatgen.analysis.phase_diagram as pd
from itertools import combinations as comb
from pymatgen.core import Composition
from rxn_help import energetics
from mp_api.client import MPRester
import pymatgen as mg
import numpy as np


mpr = MPRester('YdSP1Z8Tv8vfrLSaRCWdwp3EWvazS0vf')

def get_rxns(pair, standard_pd, grand_pd, only_max_dG=False, energ_cutoff=20.0):

    c1 = Composition(pair[0])
    c2 = Composition(pair[1])
    ir = interface_reactions.GrandPotentialInterfacialReactivity(c1, c2, grand_pd, standard_pd)

    kinks = list(ir.get_kinks())
    eqtns = [info[3] for info in kinks]
    energies = [info[2] for info in kinks]

    favorable_energs, favorable_products = [], []

    if only_max_dG:

        energ = 1000*min(energies)
        eqtn = eqtns[np.argmin(energies)]
        products = [comp.reduced_formula for comp in eqtn.products]

        favorable_energs.append(energ)
        favorable_products.append(products)

    else:

        for i, energ in enumerate(energies):

            if abs(1000*energ) > energ_cutoff:

                eqtn = eqtns[i]
                products = [comp.reduced_formula for comp in eqtn.products]

                favorable_energs.append(1000*energ)
                favorable_products.append(products)

    return favorable_energs, favorable_products

def filter_additives(reactants, target, additive_elems, T, all_anions, atmos):

    # Allow 2 additive elems to be queried at once
    additive_sublists = [additive_elems[i:i+2] for i in range(0, len(additive_elems), 2)]

    # Iterate through each group of additives
    for additive_elems in additive_sublists:

        for anion in all_anions:

            # Get elements in reactants
            all_elems = []
            for ph in reactants:
                ph_elems = [str(el) for el in Composition(ph).elements]
                all_elems += ph_elems

            # Include additive elements
            all_elems += additive_elems

            # Include additional anions
            if anion in ['OH', 'CO2']:
                all_elems += [anion[0]]
                all_elems += [anion[1]]
            else:
                all_elems += [anion]

            # Ensure all elements are unique
            unique_elems = list(set(all_elems))

            # Build phase diagram from unique elements
            all_entries = mpr.get_entries_in_chemsys(all_elems)
            standard_pd, grand_pd = energetics.get_phase_diagrams(all_elems, T, atmos)

            for elem in additive_elems:

                # Get all stable binaries
                if anion in ['OH', 'CO2']:
                    current_elems = [elem, anion[0], anion[1]]
                else:
                    current_elems = [elem, anion]
                stable_ph_info = mpr.summary.search(elements=current_elems, num_elements=[len(current_elems), len(current_elems)],
                    energy_above_hull=[0.0, 0.0], fields=['composition'])
                stable_formulae = [cmpd.composition.reduced_formula for cmpd in stable_ph_info]

                # For now, only work with two reactants
                assert len(reactants) == 2, 'Must have two reactants'

                # Iterate through rxn pathways
                all_indices = [0, 1]
                for formula in stable_formulae:
                    indices = all_indices.copy()
                    for i in all_indices:
                        pair = [formula, reactants[i]]
                        first_energies, first_psets = get_rxns(pair, standard_pd, grand_pd)
                        for energ_1, products_1 in zip(first_energies, first_psets):
                            if 'O2' in products_1:
                                products_1.remove('O2')
                            for cmpd in products_1:
                                if i == 0:
                                    pair = [cmpd, reactants[i+1]]
                                if i == 1:
                                    pair = [cmpd, reactants[i-1]]
                                second_energies, second_psets = get_rxns(pair, standard_pd, grand_pd)
                                for energ_2, products_2 in zip(second_energies, second_psets):
                                    if 'O2' in products_2:
                                        products_2.remove('O2')
                                    if target in products_2:
                                        print('Additive: %s' % formula)
                                        print('Intermediate: %s (%s meV/atom)' % (products_1, round(energ_1, 0)))
                                        print('Final: %s (%s meV/atom)' % (products_2, round(energ_2, 0)))
