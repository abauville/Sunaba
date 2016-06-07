# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json
from InputDef import *


print("\n"*5)

grid = Grid()
numerics = Numerics()
particles = Particles()
physics = Physics()
visu = Visu()
char = Char()
bcStokes = BCStokes()
bcThermal = BCThermal()



Phase0 = Material()
Phase1 = Material()

Phase0.name = "Matrix"


Phase1.name = "Inclusion"
Phase1.eta0 = 100*Phase0.eta0
Phase1.rho0 = 5*Phase0.rho0

numerics.nTimeSteps = -1
bcStokes.backStrainRate = -0.0E-1
numerics.CFL_fac = 0.5

MatProps = {'0': Phase0.__dict__, '1': Phase1.__dict__}

Description = "This is a test input file. Which defines to materials: a matrix and an inclusion 100 times stronger in a square box in pure shear"

#myJsonFile = dict(Description = Description, Grid = grid.__dict__, Numerics = numerics.__dict__, Particles = particles.__dict__, Physics = physics.__dict__, Visu = visu.__dict__, MatProps = MatProps, Char = char.__dict__, BCStokes = bcStokes.__dict__, BCThermal = bcThermal.__dict__);
myJsonFile = dict(Grid = grid.__dict__, Numerics = numerics.__dict__, Particles = particles.__dict__, Physics = physics.__dict__, Visu = visu.__dict__, MatProps = MatProps, Char = char.__dict__, BCStokes = bcStokes.__dict__, BCThermal = bcThermal.__dict__);

outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, separators=(',', ': '), ensure_ascii=False)

#print(myJsonFile["Description"])
