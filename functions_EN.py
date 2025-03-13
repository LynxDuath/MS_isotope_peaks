import sys
from itertools import repeat
from numpy import full
from numpy import array
from numpy import dot

#Inputs
def get_inputs() -> (str, int, str, str):
    """
    Asks for the required inputs.
    @return formula: str, formula of the molecule.
    @return charge: int, charge of the molecule.
    @return precition: str, preferred precition.
    @return graphics: str, preferred graphic output.
    """
    print("This programm only supports the following elements:\nC, H, O, N, F, Cl, Br, I, S, P, Si, B\nAdditional:\nMe, Et, Bu, Ph, Bz, Ts, Ms and parenthesis")
    formula = str(input("Please insert your formula: "))
    charge = input("Please insert your charge: ")
    try:
        charge = int(charge)
    except ValueError:
        if charge[-1] == '-':
            charge = -1*int(charge[:-1])
        else:
            charge = int(charge[:-1])
    precition = str(input("Do you want the result to be precise or unprecise (u/p/up)? "))
    graphics = str(input("Do you want graphic output (Y/N)? "))
    print()
    return formula, charge, precition, graphics

def get_inputs_inverted() -> list:
    """
    Asks for the required inputs.
    @return Investigate: list, distribution of the unknown molecule.
    """
    investigating = str(input("Which molecule shall be investigated? "))
    with open(investigating + '.py', 'r') as investigate:
        Investigate = []
        for line in investigate:
            Investigate.append([float(line.strip('\n').split()[0]), float(line.strip('\n').split()[1])])
    return Investigate

#Funktion zur Dechiffrierung der Summenformel
def decipher_formula(formula: str, level: int, element, carbon=0, hydrogen=0, oxygen=0, nitrogen=0, fluorine=0, chlorine=0, bromine=0, iodine=0, sulfur=0, phosphorus=0, silicium=0, boron=0) -> array:
    """
    Deciphers the formula.
    @param formula: str, formula of the molecule.
    @param level: int, level of the formula for parenthesis.
    @params carbon, hydrogen,...: int, initial number of atoms in the molecule, default 0.
    @returns carbon, hydrogen,...: int, total number of atoms in the molecule.
    """
    digits = ('0', '1', '2', '3', '4', '5', '6', '7', '8',' 9')
    upper_case_letters = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '(', ')')
    while element < len(formula):
        number = 0
        amount = ''
        if formula[element] == 'A':
            if formula[element+1] == 'c':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += 2*int(amount)
                    hydrogen += 3*int(amount)
                    oxygen += int(amount)
                else:
                    carbon += 2
                    hydrogen += 3
                    oxygen += 1
        elif formula[element] == 'B':
            if formula[element+1] == 'r':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    bromine += int(amount)
                else:
                    bromine += 1
            elif formula[element+1] == 'u':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += 4*int(amount)
                    hydrogen += 9*int(amount)
                else:
                    carbon += 4
                    hydrogen += 9
            elif formula[element+1] == 'z':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += 7*int(amount)
                    hydrogen += 7*int(amount)
                else:
                    carbon += 7
                    hydrogen += 7
            else:
                if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                    while formula[element+1+number] in digits:
                        amount += formula[element+1+number]
                        number += 1
                        if element+1+number == len(formula):
                            break
                    boron += int(amount)
                else:
                    boron += 1
        elif formula[element] == 'C':
            if formula[element+1] == 'l':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    chlorine += int(amount)
                else:
                    chlorine += 1
            else:
                if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                    while formula[element+1+number] in digits:
                        amount += formula[element+1+number]
                        number += 1
                        if element+1+number == len(formula):
                            break
                    carbon += int(amount)
                else:
                    carbon += 1
        elif formula[element] == 'E':
            if formula[element+1] == 't':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += 2*int(amount)
                    hydrogen += 5*int(amount)
                else:
                    carbon += 2
                    hydrogen += 5
        elif formula[element] == 'F':
            if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                while formula[element+1+number] in digits:
                    amount += formula[element+1+number]
                    number += 1
                    if element+1+number == len(formula):
                        break
                fluorine += int(amount)
            else:
                fluorine += 1
        elif formula[element] == 'H':
            if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                while formula[element+1+number] in digits:
                    amount += formula[element+1+number]
                    number += 1
                    if element+1+number == len(formula):
                        break
                hydrogen += int(amount)
            else:
                hydrogen += 1
        elif formula[element] == 'I':
            if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                while formula[element+1+number] in digits:
                    amount += formula[element+1+number]
                    number += 1
                    if element+1+number == len(formula):
                        break
                iodine += int(amount)
            else:
                iodine += 1
        elif formula[element] == 'M':
            if formula[element+1] == 'e':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += int(amount)
                    hydrogen += 3*int(amount)
                else:
                    carbon += 1
                    hydrogen += 3
            if formula[element+1] == 's':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += int(amount)
                    hydrogen += 3*int(amount)
                    oxygen += 2*int(amount)
                    sulfur += int(amount)
                else:
                    carbon += 1
                    hydrogen += 3
                    oxygen += 2
                    sulfur += 1
        elif formula[element] == 'N':
            if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                while formula[element+1+number] in digits:
                    amount += formula[element+1+number]
                    number += 1
                    if element+1+number == len(formula):
                        break
                nitrogen += int(amount)
            else:
                nitrogen += 1
        elif formula[element] == 'O':
            if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                while formula[element+1+number] in digits:
                    amount += formula[element+1+number]
                    number += 1
                    if element+1+number == len(formula):
                        break
                oxygen += int(amount)
            else:
                oxygen += 1
        elif formula[element] == 'P':
            if formula[element+1] == 'h':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += 6*int(amount)
                    hydrogen += 5*int(amount)
                else:
                    carbon += 6
                    hydrogen += 5
            else:
                if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                    while formula[element+1+number] in digits:
                        amount += formula[element+1+number]
                        number += 1
                        if element+1+number == len(formula):
                            break
                    phosphorus += int(amount)
                else:
                    phosphorus += 1
        elif formula[element] == 'S':
            if formula[element+1] == 'i':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    silicium += int(amount)
                else:
                    silicium += 1
            else:
                if element+1+number < len(formula) and formula[element+1+number] not in upper_case_letters:
                    while formula[element+1+number] in digits:
                        amount += formula[element+1+number]
                        number += 1
                        if element+1+number == len(formula):
                            break
                    sulfur += int(amount)
                else:
                    sulfur += 1
        elif formula[element] == 'T':
            if formula[element+1] == 's':
                if element+2+number < len(formula) and formula[element+2+number] not in upper_case_letters:
                    while formula[element+2+number] in digits:
                        amount += formula[element+2+number]
                        number += 1
                        if element+2+number == len(formula):
                            break
                    carbon += 7*int(amount)
                    hydrogen += 7*int(amount)
                    oxygen += 2*int(amount)
                    sulfur += int(amount)
                else:
                    carbon += 7
                    hydrogen += 7
                    oxygen += 2
                    sulfur += 1
        elif formula[element] == '(':
            exec('carbon' + str(level) + ' = carbon')
            exec('hydrogen' + str(level) + ' = hydrogen')
            exec('oxygen' + str(level) + ' = oxygen')
            exec('nitrogen' + str(level) + ' = nitrogen')
            exec('fluorine' + str(level) + ' = fluorine')
            exec('chlorine' + str(level) + ' = chlorine')
            exec('bromine' + str(level) + ' = bromine')
            exec('iodine' + str(level) + ' = iodine')
            exec('sulfur' + str(level) + ' = sulfur')
            exec('phosphorus' + str(level) + ' = phosphorus')
            exec('silicium' + str(level) + ' = silicium')
            exec('boron' + str(level) + ' = boron')
            level += 1
            carbon, hydrogen, oxygen, nitrogen, fluorine, chlorine, bromine, iodine, sulfur, phosphorus, silicium, boron, element = decipher_formula(formula, level, element + 1, carbon, hydrogen, oxygen, nitrogen, fluorine, chlorine, bromine, iodine, sulfur, phosphorus, silicium, boron)
            level -= 1
            if element+number < len(formula) and formula[element+number] not in upper_case_letters:
                while formula[element+number] in digits:
                    amount += formula[element+number]
                    number += 1
                    if element+number == len(formula):
                        break
                carbon += (carbon - eval('carbon' + str(level)))*(int(amount) - 1)
                hydrogen += (hydrogen - eval('hydrogen' + str(level)))*(int(amount) - 1)
                oxygen += (oxygen - eval('oxygen' + str(level)))*(int(amount) - 1)
                nitrogen += (nitrogen - eval('nitrogen' + str(level)))*(int(amount) - 1)
                fluorine += (fluorine - eval('fluorine' + str(level)))*(int(amount) - 1)
                chlorine += (chlorine - eval('chlorine' + str(level)))*(int(amount) - 1)
                bromine += (bromine - eval('bromine' + str(level)))*(int(amount) - 1)
                iodine += (iodine - eval('iodine' + str(level)))*(int(amount) - 1)
                sulfur += (sulfur - eval('sulfur' + str(level)))*(int(amount) - 1)
                phosphorus += (phosphorus - eval('phosphorus' + str(level)))*(int(amount) - 1)
                silicium += (silicium - eval('silicium' + str(level)))*(int(amount) - 1)
                boron += (boron - eval('boron' + str(level)))*(int(amount) - 1)
        elif formula[element] == ')':
            element += 1
            break
        element += 1
    return array([carbon, hydrogen, oxygen, nitrogen, fluorine, chlorine, bromine, iodine, sulfur, phosphorus, silicium, boron, element])

#Funktion zur Berechnung der Verteilung der Atomsorten
def distribute(Atomtype: list, atomtype: int, AtomtypeBase: list) -> list:
    """
    Calculates the distribution of additional masses for the selected atom type.
    @param Atomtype: list, distribution of the selected atom type before calculations, default [1].
    @param atomtype: int, number of atoms of selected type in the molecule.
    @param AtomtypeBase: list, base distribution for one atom of the selected type.
    @return Atomtype: list, distribution of the selected atom type after calculations.
    """
    for a in range(1, atomtype+1):
        AtomtypeMem = Atomtype
        Atomtype = list(repeat(0,len(AtomtypeBase)-1))
        for r in range(len(AtomtypeMem)):
            Atomtype.append(0)
            for d in range(len(AtomtypeBase)):
                Atomtype[r+d] += AtomtypeMem[r]*AtomtypeBase[d]
    return Atomtype

#Normieren der Intensitäten
def norm(Molecule: list, maxprob: float, rounding: int, threshold: float, maxinten=100) -> list:
    """
    Scaling of the intensity distribution.
    @param Molecule: list, initial distribution of the molecule.
    @param maxprob: float, maximum propability for scaling.
    @param rounding: int, decimal place to round the intensities.
    @param threshold: float, threshold of the intensity to show the peaks.
    @param maxinten: float, maximum intensity to be achieved, default 100.
    @return MoleculeNorm: list, scaled distribution of the molecule.
    """
    MoleculeNorm=[]
    for l in range(len(Molecule)):
        massrel = round(Molecule[l][1]/maxprob*maxinten,rounding)
        if massrel >= threshold:
            MoleculeNorm.append([Molecule[l][0],massrel])
    return MoleculeNorm

#Graphische Ausgabe
def graph(graphics: str, MoleculeNorm: list) -> None:
    """
    Graphic outputof the results.
    @param graphics: str, preferred graphic output.
    @param MoleculeNorm: list, normed distribution of the molecule.
    """
    if graphics == 'Y':
        print("\n---Graphic Output---\n")
        Graphics = full((11,len(MoleculeNorm)),'       ')
        for g1 in range(10):
            for g2 in range(len(MoleculeNorm)):
                if round(MoleculeNorm[g2][1]/10) >= 10-g1:
                    Graphics[g1][g2] = 'iiiiiii'
        for g3 in range(len(MoleculeNorm)):
            Graphics[10][g3] = str(MoleculeNorm[g3][0])
        for intensity in Graphics:
            print(*intensity, sep='\t')
        print()

#Output inverted
def print_inverted(mol_rms_min: float, best_molecule: list) -> None:
    """
    Prints out the given molecule as a formula.
    @param mol_rms_min: float, minimal discrepancy from the input.
    @param best_molecule: list, stom distribution of the best approximation.
    """
    output = 'With a match of ' + str((1-mol_rms_min)*100) + '% your molecule is '
    if best_molecule[0] != 0:
        output += 'C'
        if best_molecule[0] > 1:
            output += str(best_molecule[0])
    if best_molecule[1] != 0:
        output += 'H'
        if best_molecule[1] > 1:
            output += str(best_molecule[1])
    if best_molecule[2] != 0:
        output += 'O'
        if best_molecule[2] > 1:
            output += str(best_molecule[2])
    if best_molecule[3] != 0:
        output += 'N'
        if best_molecule[3] > 1:
            output += str(best_molecule[3])
    if best_molecule[4] != 0:
        output += 'F'
        if best_molecule[4] > 1:
            output += str(best_molecule[4])
    if best_molecule[5] != 0:
        output += 'Cl'
        if best_molecule[5] > 1:
            output += str(best_molecule[5])
    if best_molecule[6] != 0:
        output += 'Br'
        if best_molecule[6] > 1:
            output += str(best_molecule[6])
    if best_molecule[7] != 0:
        output += 'I'
        if best_molecule[7] > 1:
            output += str(best_molecule[7])
    if best_molecule[8] != 0:
        output += 'S'
        if best_molecule[8] > 1:
            output += str(best_molecule[8])
    if best_molecule[9] != 0:
        output += 'P'
        if best_molecule[9] > 1:
            output += str(best_molecule[9])
    if best_molecule[10] != 0:
        output += 'Si'
        if best_molecule[10] > 1:
            output += str(best_molecule[10])
    if best_molecule[11] != 0:
        output += 'B'
        if best_molecule[11] > 1:
            output += str(best_molecule[11])
    output += '.'
    print(output)

#Additional
def wait():
    print("---Please hold the line!---\n")

def unprecition():
    sys.exit("Your answer is quite too unprecise!")

def smaller(MoleculeNorm, Investigate):
    if MoleculeNorm >= Investigate:
        return Investigate
    else:
        return MoleculeNorm

#Makrofunktionen
def calculate_unprecise(atoms: list, massmin: int, charge: int, CarbonBase, HydrogenBase, OxygenBase, NitrogenBase, FluorineBase, ChlorineBase, BromineBase, IodineBase, SulfurBase, PhosphorusBase, SiliciumBase, BoronBase, maxinten=100) -> list:
    """
    Calculates the precise distribution of the Molecule.
    @param atoms: list, number of the atoms in the molecule.
    @param massmin: float, mass of the lightest isotope of the molecule.
    @param charge: int, charge of the molecule.
    @param maxinten: float, maximum intensity to be achieved, default 100.
    @return MoleculeNorm: list, scaled distribution of the molecule. 
    """
    #Addition
    Molecule = [1]
    Molecule = distribute(Molecule, atoms[0], CarbonBase)
    Molecule = distribute(Molecule, atoms[1], HydrogenBase)
    Molecule = distribute(Molecule, atoms[2], OxygenBase)
    Molecule = distribute(Molecule, atoms[3], NitrogenBase)
    Molecule = distribute(Molecule, atoms[4], FluorineBase)
    Molecule = distribute(Molecule, atoms[5], ChlorineBase)
    Molecule = distribute(Molecule, atoms[6], BromineBase)
    Molecule = distribute(Molecule, atoms[7], IodineBase)
    Molecule = distribute(Molecule, atoms[8], SulfurBase)
    Molecule = distribute(Molecule, atoms[9], PhosphorusBase)
    Molecule = distribute(Molecule, atoms[10], SiliciumBase)
    Molecule = distribute(Molecule, atoms[11], BoronBase)
    
    for m in range(len(Molecule)):
        Molecule[m] = [round((m + massmin)/abs(charge), 2), Molecule[m]]
    
    #Normieren der Wahrscheinlichkeiten auf relative Intensitäten
    maxprob = 0
    for q in range(len(Molecule)):
        if Molecule[q][1] > maxprob:
            maxprob = Molecule[q][1]
    
    MoleculeNorm = norm(Molecule, maxprob, 1, 1.0, maxinten)
    return MoleculeNorm

def calculate_precise(atoms: list, massmin: float, charge: int, IsotopesMass: list, CarbonBase, HydrogenBase, OxygenBase, NitrogenBase, FluorineBase, ChlorineBase, BromineBase, IodineBase, SulfurBase, PhosphorusBase, SiliciumBase, BoronBase, maxinten=100) -> list:
    """
    Calculates the precise distribution of the Molecule.
    @param atoms: list, number of the atoms in the molecule.
    @param massmin: float, mass of the lightest isotope of the molecule.
    @param charge: int, charge of the molecule.
    @param IsotopesMass: list, base information of the isotope masses.
    @param maxinten: float, maximum intensity to be achieved, default 100.
    @return MoleculeNorm: list, scaled distribution of the molecule. 
    """
    #Distribution and Addition
    Carbon = distribute([1], atoms[0], CarbonBase)
    Hydrogen = distribute([1], atoms[1], HydrogenBase)
    Oxygen = distribute([1], atoms[2], OxygenBase)
    Nitrogen = distribute([1], atoms[3], NitrogenBase)
    Fluorine = distribute([1], atoms[4], FluorineBase)
    Chlorine = distribute([1], atoms[5], ChlorineBase)
    Bromine = distribute([1], atoms[6], BromineBase)
    Iodine = distribute([1], atoms[7], IodineBase)
    Sulfur = distribute([1], atoms[8], SulfurBase)
    Phosphorus = distribute([1], atoms[9], PhosphorusBase)
    Silicium = distribute([1], atoms[10], SiliciumBase)
    Boron = distribute([1], atoms[11], BoronBase)
    
    Molecule = []    
    for c1 in range(atoms[0]+1):
        for h1 in range(atoms[1]+1):
            for o1 in range(atoms[2]+1):
                for o2 in range(atoms[2]-o1+1):
                    for n1 in range(atoms[3]+1):
                        for cl2 in range(atoms[5]+1):
                            for br2 in range(atoms[6]+1):
                                for s1 in range(atoms[8]+1):
                                    for s2 in range(atoms[8]-s1+1):
                                        for s4 in range(atoms[8]-s1-s2+1):
                                            for si1 in range(atoms[10]+1):
                                                for si2 in range(atoms[10]-si1+1):
                                                    for b1 in range(atoms[11]+1):
                                                        Isotopes = array([c1, h1, o1, o2, n1, cl2, br2, s1, s2, s4, si1, si2, b1])
                                                        mass = massmin+dot(Isotopes, IsotopesMass)
                                                        Molecule.append([mass/abs(charge),Carbon[c1]*Hydrogen[h1]*Oxygen[o1+2*o2]*Nitrogen[n1]*Chlorine[2*cl2]*Bromine[2*br2]*Sulfur[s1+2*s2+4*s4]*Silicium[si1+2*si2]*Boron[b1]])
    
    #Scaling the propabilities to relative intensities
    maxprob = 0
    Molecule = sorted(Molecule, key=lambda mass:mass[0])
    Molecule[0][0] = round(Molecule[0][0],5)
    for j1 in range(1,len(Molecule)):
        j2 = 0
        if round(Molecule[j1][0],5) == Molecule[j1-1-j2][0]:
            Molecule[j1-1-j2][1] += round(Molecule[j1][1])
            j2 += 1
        else:
            Molecule[j1][0] = round(Molecule[j1][0],5)
            j2 = 0
        if Molecule[j1-1-j2][1] > maxprob:
            maxprob = Molecule[j1-1-j2][1]
    
    MoleculeNorm = norm(Molecule, maxprob, 3, 0.001)
    return MoleculeNorm

def invert(Investigate: list, investigate_massmin: int, maxinten, CarbonBase, HydrogenBase, OxygenBase, NitrogenBase, FluorineBase, ChlorineBase, BromineBase, IodineBase, SulfurBase, PhosphorusBase, SiliciumBase, BoronBase) -> list:
    """
    Calculates possible compositions of an unknown molecule with a certain mass distribution.
    @param Investigate: list, distribution of the unknown molecule.
    @param investigate_massmin: int: mass of the lowest isotope of the unknown molecule.
    @param maxinten: float, maximum intensity to be achieved.
    @return possible_molecules: list, possible compositions of the unknown molecule.
    """
    if len(Investigate) >1:
        if Investigate[1][1]/Investigate[0][1] >= 1:
            sys.exit("Your molecule is too big to calculate,a cluster/fullerene, an overlap of two molecules (e.g. acidity), or contains boron.\nEither way, this program won't work.")

    #Anzahl Stickstoff
    if investigate_massmin%2 == 0:
        nitrogen_proto = 'even'
    else:
        nitrogen_proto = 'odd'

    #Auswertung des [M+2] Peaks
    chlorine = 0
    bromine = 0
    ssi_proto = 'impossible'
    if len(Investigate) > 2:
        if Investigate[2][1]/Investigate[0][1] >= 0.03:
            ssi_proto = 'possible'
            if Investigate[2][1]/Investigate[0][1] >= 0.24:
                hal = (len(Investigate)-1)//2
                hal_rms_min = 1e200
                for cl in range(hal+1):
                    br = hal - cl
                    Hal = [1]
                    Hal = distribute(Hal, cl, ChlorineBase)
                    Hal = distribute(Hal, br, BromineBase)
                    for m in range(len(Hal)):
                        Hal[m] = [m, Hal[m]]
                    maxprob = 0
                    for q in range(len(Hal)):
                        if Hal[q][1] > maxprob:
                            maxprob = Hal[q][1]
                    HalNorm = norm(Hal, maxprob, 1, 0, maxinten)
                    hal_rms = 0
                    for hal1 in range(len(HalNorm)):
                        if hal1%2 == 0:
                            hal_rms += (Investigate[hal1][1]/HalNorm[hal1][1] - 1)**2
                    hal_rms = hal_rms**0.5
                    if hal_rms <= hal_rms_min:
                        hal_rms_min = hal_rms
                        chlorine = cl
                        bromine = br
                hal -= 1
                for cl in range(hal+1):
                    br = hal - cl
                    Hal = [1]
                    Hal = distribute(Hal, cl, ChlorineBase)
                    Hal = distribute(Hal, br, BromineBase)
                    for m in range(len(Hal)):
                        Hal[m] = [m, Hal[m]]
                    maxprob = 0
                    for q in range(len(Hal)):
                        if Hal[q][1] > maxprob:
                            maxprob = Hal[q][1]
                    HalNorm = norm(Hal, maxprob, 1, 0, maxinten)
                    hal_rms = 0
                    for hal1 in range(len(HalNorm)):
                        if hal1%2 == 0:
                            hal_rms += (Investigate[hal1][1]/HalNorm[hal1][1] - 1)**2
                    hal_rms = hal_rms**0.5
                    if hal_rms <= hal_rms_min:
                        hal_rms_min = hal_rms
                        chlorine = cl
                        bromine = br
    
    #Auswertung des [M+1] Peaks
    if len(Investigate) > 1:
        carbon_proto = Investigate[1][1]/Investigate[0][1]/0.011
    
    #Mögliche sonstige Atome
    restmass = investigate_massmin - 35*chlorine - 79*bromine
    possible_molecules = []
    for c in range(int(restmass//12)+1):
        if c <= carbon_proto+1:
            for n in range(int((restmass-12*c)//14)+1):
                if (nitrogen_proto == 'even' and n%2 == 0) or (nitrogen_proto == 'odd' and n%2 == 1):
                    if c+n/3 - carbon_proto >= -1 and c+n/3 - carbon_proto <= 1:
                        for i in range(int((restmass-12*c-14*n)//127)+1):
                            for p in range(int((restmass-12*c-14*n-127*i)//31)+1):
                                for f in range(int((restmass-12*c-14*n-127*i-31*p)//17)+1):
                                    for o in range(int((restmass-12*c-14*n-127*i-31*p-17*f)//16)+1):
                                        if ssi_proto == 'possible':
                                            for s in range(int((restmass-12*c-14*n-127*i-31*p-17*f-16*o)//32)+1):
                                                for si in range(int((restmass-12*c-14*n-127*i-31*p-17*f-16*o-32*s)//28)+1):
                                                    h = int(restmass-12*c-14*n-127*i-31*p-17*f-16*o-32*s-28*si)
                                                    possible_molecules.append([c, h, o, n, f, chlorine, bromine, i, s, p, si, 0])
                                        elif ssi_proto == 'impossible':
                                            h = int(restmass-12*c-14*n-127*i-31*p-17*f-16*o)
                                            possible_molecules.append([c, h, o, n, f, chlorine, bromine, i, 0, p, 0, 0])

    return possible_molecules
