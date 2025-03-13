import sys
from itertools import repeat
from numpy import full
from numpy import array
from numpy import dot

from resources import *

#Sprachabfrage
language = str(input('Select your language (EN/DE): '))
if language == 'EN':
    from functions_EN import *
elif language == 'DE':
    from functions_DE import *
else:
    sys.exit("I don't understand Klingon!")
Investigate = get_inputs_inverted()

#Berechnung
investigate_massmin = Investigate[0][0]
maxinten = 0
for q in range(len(Investigate)):
    if Investigate[q][1] > maxinten:
        maxinten = Investigate[q][1]
possible_molecules = invert(Investigate, investigate_massmin, maxinten, CarbonBase, HydrogenBase, OxygenBase, NitrogenBase, FluorineBase, ChlorineBase, BromineBase, IodineBase, SulfurBase, PhosphorusBase, SiliciumBase, BoronBase)

#Genauigkeit
mol_rms_min = 1e200
for x in range(len(possible_molecules)):
    MoleculeNorm = calculate_unprecise(possible_molecules[x], investigate_massmin, 1, CarbonBase, HydrogenBase, OxygenBase, NitrogenBase, FluorineBase, ChlorineBase, BromineBase, IodineBase, SulfurBase, PhosphorusBase, SiliciumBase, BoronBase, maxinten)

    mol_rms = 0
    for mz in range(smaller(len(MoleculeNorm), len(Investigate))):
        mol_rms += (Investigate[mz][1]/MoleculeNorm[mz][1] - 1)**2
    mol_rms = mol_rms**0.5
    if mol_rms <= mol_rms_min:
        mol_rms_min = mol_rms
        best_molecule = possible_molecules[x]

#Ausgabe
print_inverted(mol_rms_min, best_molecule)
