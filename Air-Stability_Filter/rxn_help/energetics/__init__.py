from pymatgen.core.composition import Composition
import pymatgen.entries.computed_entries as ce
import pymatgen.analysis.phase_diagram as pd
from pymatgen.core import Element
from mp_api.client import MPRester
import numpy as np


def get_phase_diagrams(available_precursors, temperature, atmos, open_el='O'):
    """
    Get the standard phase diagram and the grand potential
    phase diagram, open w.r.t. gaseous species (O2).
    """

    mpr = MPRester('YdSP1Z8Tv8vfrLSaRCWdwp3EWvazS0vf')

    # Partial pressures of gaseous species in air
    if atmos == 'air':
        p_O2 = 21200
        p_CO2 = 4050
        p_NH3 = 16.
        p_H2O = 2300.

    # Estimation based on 1e-6 partial pressure
    elif atmos == 'inert':
        p_O2 = 0.1
        p_CO2 = 0.1
        p_NH3 = 0.1
        p_H2O = 0.1

    if open_el == 'O':
        chempot = -4.9480 # DFT-calculated energy (0 K)
        chempot += get_chempot_correction('O', temperature, p_O2)

        C_chempot = -9.2268 # DFT-calculated energy (0 K)
        C_entry = ce.ComputedEntry(Composition('C'), C_chempot)
        CO2_chempot = -8.1443583 + get_chempot_correction('CO2', temperature, p_CO2)
        CO2_entry = ce.ComputedEntry(Composition('CO2'), CO2_chempot*3)

        H2O_chempot = -5.19275 + get_chempot_correction('H2O', temperature, p_H2O)
        H2O_entry = ce.ComputedEntry(Composition('H2O'), H2O_chempot*3)

    else:
        assert False, 'Only oxygen implemented as open element'

    # Determine elements from precursors
    elems = []
    for cmpd in available_precursors:
        comp = Composition(cmpd)
        elems += [str(el) for el in comp.elements]
    elems = list(set(elems))

    # Get all entries in the chemical space
    entries = mpr.get_entries_in_chemsys(elems)
    entries += [C_entry, CO2_entry, H2O_entry]

    # Build phase diagrams
    standard_pd = pd.PhaseDiagram(entries)
    grand_pd = pd.GrandPotentialPhaseDiagram(entries, {Element('O'): chempot})

    return standard_pd, grand_pd

def get_chempot_correction(element, temp, pres):
    """
    Get the normalized correction term Δμ for chemical potential of a gas
    phase consisting of element at given temperature and pressure,
    referenced to that in the standard state (T_std = 298.15 K,
    T_std = 1 bar). The gas phase is limited to be one of O2, N2, Cl2,
    F2, H2. Calculation formula can be found in the documentation of
    Materials Project website.

    Args:
        element (string): The string representing the element.
        temp (float): The temperature of the gas phase.
        pres (float): The pressure of the gas phase.

    Returns:
        The correction of chemical potential in eV/atom of the gas
        phase at given temperature and pressure.
    """

    EV_TO_KJ_PER_MOL = 96.4853

    if element not in ['O', 'N', 'Cl', 'F', 'H', 'CO2', 'NH3', 'H2O']:
        return 0
    std_temp = 298.15
    std_pres = 1E5
    ideal_gas_const = 8.3144598
    # Cp and S at standard state in J/(K.mol). Data from
    # https://janaf.nist.gov/tables/O-029.html
    # https://janaf.nist.gov/tables/N-023.html
    # https://janaf.nist.gov/tables/Cl-073.html
    # https://janaf.nist.gov/tables/F-054.html
    # https://janaf.nist.gov/tables/H-050.html
    Cp_dict = {'O': 29.376,
               'N': 29.124,
               'Cl': 33.949,
               'F': 31.302,
               'H': 28.836,
               'CO2': 37.129,
               'NH3': 35.640,
               'H2O': 33.22}

    S_dict = {'O': 205.147,
              'N': 191.609,
              'Cl': 223.079,
              'F': 202.789,
              'H': 130.680,
              'CO2': 213.79,
              'NH3': 192.80,
              'H2O': 194.10}
    Cp_std = Cp_dict[element]
    S_std = S_dict[element]
    PV_correction = ideal_gas_const * temp * np.log(pres / std_pres)
    TS_correction = - Cp_std * (temp * np.log(temp) - std_temp * np.log(std_temp)) \
        + Cp_std * (temp - std_temp) * (1 + np.log(std_temp)) \
        - S_std * (temp - std_temp)

    dG = PV_correction + TS_correction

    # Convert to eV/molecule unit.
    dG /= 1000 * EV_TO_KJ_PER_MOL

    # Normalize by number of atoms in the gas molecule
    if element == 'H2O':
        dG /= 3
    if element == 'CO2':
        dG /= 3
    if element == 'NH3':
        dG /= 4
    if element == 'O':
        dG /= 2

    return dG
