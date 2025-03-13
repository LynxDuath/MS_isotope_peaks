import sys
from itertools import repeat
from numpy import full
from numpy import array
from numpy import dot

#Relevante Basisdaten
CarbonBase = [0.9890, 0.0110]
HydrogenBase = [0.99985, 0.00015]
OxygenBase = [0.99762, 0.00038, 0.00200]
NitrogenBase = [0.99634, 0.00366]
FluorineBase = [1]
ChlorineBase = [0.7577, 0, 0.2423]
BromineBase = [0.5069, 0, 0.4931]
IodineBase = [1]
SulfurBase = [0.9502, 0.0075, 0.0421, 0, 0.0002]
PhosphorusBase = [1]
SiliciumBase = [0.9223, 0.0467, 0.0310]
BoronBase = [0.199, 0.801]

Massmin_unprecise = array([12, 1, 16, 14, 19, 35, 79, 127, 32, 31, 28, 10])
Massmin_precise = array([12.0, 1.007825032239, 15.9949146195717, 14.0030740044320, 18.99984031627392, 34.96885268237, 78.918337614, 126.90447, 31.972071174414, 30.9737619984270, 27.9769265346544, 10.0129369541])
IsotopesMass = array([1.0033548350723, 1.006276746, 1.004217137, 2.004244993, 0.9970348945, 1.99704992, 1.9979521, 0.9993877354, 1.99579583, 3.995009538, 0.9995681303, 1.996843602, 0.9963684104])
