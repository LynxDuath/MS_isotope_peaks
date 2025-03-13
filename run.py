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

inputs = get_inputs()
"""
formula = inputs[0]
charge = inputs[1]
precition = inputs[2]
graphics = inputs[3]
"""
atoms = decipher_formula(inputs[0], 0, 0)[:-1]
"""
carbon = atoms[0]
hydrogen = atoms[1]
oxygen = atoms[2]
nitrogen = atoms[3]
fluorine = atoms[4]
chlorine = atoms[5]
bromine = atoms[6]
iodine = atoms[7]
sulfur = atoms[8]
phosphorus = atoms[9]
silicium = atoms[10]
boron = atoms[11]
"""

#Ungenaues Ergebnis
if inputs[2] == 'u' or inputs[2] == 'up':
    massmin = dot(Massmin_unprecise, atoms)

    #Berechnung
    MoleculeNorm = calculate_unprecise(atoms, massmin, inputs[1], CarbonBase, HydrogenBase, OxygenBase, NitrogenBase, FluorineBase, ChlorineBase, BromineBase, IodineBase, SulfurBase, PhosphorusBase, SiliciumBase, BoronBase)
        
    #Ausgabe
    for mz in MoleculeNorm:
        print(*mz, sep='\t')
    graph(inputs[3], MoleculeNorm)
    
elif inputs[2] == 'p':
    None
else:
    unprecition()

#Genaues Ergebnis
if inputs[2] == 'p' or inputs[2] == 'up':
    wait()
    massmin = dot(Massmin_precise, atoms)-0.00054857990906516*inputs[1]

    #Berechnung
    MoleculeNorm = calculate_precise(atoms, massmin, inputs[1], IsotopesMass, CarbonBase, HydrogenBase, OxygenBase, NitrogenBase, FluorineBase, ChlorineBase, BromineBase, IodineBase, SulfurBase, PhosphorusBase, SiliciumBase, BoronBase)
    
    #Ausgabe
    for mz in MoleculeNorm:
        print(*mz, sep='\t')
    graph(inputs[3], MoleculeNorm)
    
elif inputs[2] == 'u':
    None
else:
    unprecition()
