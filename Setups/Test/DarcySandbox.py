# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json
from InputDef import *
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *

print("\n"*5)

##             Description
## =====================================
Description = "This is a test input file. Which defines to materials: a matrix and an inclusion 100 times stronger in a square box in pure shear"

##      Declare singleton objects
## =====================================
Grid = Grid()
Numerics = Numerics()
Particles = Particles()
Physics = Physics(True)
Visu = Visu()
Char = Char()
BCStokes = BCStokes()
BCThermal = BCThermal()
Geometry = {}


##          Material properties
## =====================================
Phase0 = Material("StickyAir")
Phase1   = Material("Sediments")
PhaseRef = Phase1

#Phase0.eta0 = 1e19
#Phase0.G    = 1e10
Phase0.rho0 = 1000.0
Phase0.eta0 = 1e18

Phase1.eta0 = 1e23
Phase1.G    = 1e11
Phase0.G    = 1e11

Backphi = 0.001
RefPerm = 5e-20
Phase0.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))
RefPerm = 5e-20
Phase1.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))

Phase1.isRef = True

MatProps = {'0': Phase0.__dict__, '1': Phase1.__dict__}



##              Some info
## =====================================
Backphi = 0.01
RefPerm = Phase0.perm0*(Backphi * Backphi *Backphi  *  (1.0-Backphi)*(1.0-Backphi))
CompactionLength = sqrt(4/3*RefPerm/Physics.eta_f * (PhaseRef.eta0/Backphi))
DeltaRho = (Physics.rho_f-PhaseRef.rho0)
CompactionVelocity = abs(RefPerm/(Physics.eta_f*Backphi) * (1-Backphi) * DeltaRho*Physics.gy)
CompactionTime = CompactionLength/CompactionVelocity

print("CompactionLength: " + str(CompactionLength) + " m")
print("CompactionTime: " + str(CompactionTime/(3600*24*365*1e6)) + " Myr")

##                 BC
## =====================================
#BCStokes.SetupType = "PureShear"
BCStokes.SetupType = "SandBox"
BCThermal.SetupType = "SandBox"









##              Particles
## =====================================
Particles.noiseFactor = 0.95
Particles.nPCX = 3
Particles.nPCY = 3





##              Grid
## =====================================
#Grid.xmin = -32*CompactionLength
#Grid.xmax =  32*CompactionLength
#Grid.ymin =  -1.0*(Grid.xmax-Grid.xmin)
#Grid.ymax =  0.0*(Grid.xmax-Grid.xmin)

#RefinementFac = 2.0
#Grid.nyC = round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
#Grid.nxC = round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.xmin = -10e3
Grid.xmax =  0.0
Grid.ymin =  0.0
Grid.ymax = 4e3;
Grid.nxC = 128#round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 96#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.fixedBox = False

##              Non Dim
## =====================================
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)
Char.set_based_on_lithostatic_pressure(PhaseRef,BCThermal,Physics,Grid)






##              Geometry
## =====================================
H = Grid.ymax-Grid.ymin
L = Grid.xmax-Grid.xmin
Hsed = 1e3
Physics.y_oceanSurface = Hsed+7e3
#DepthWater = H/2.0
#TopWater = Hsed+DepthWater

air = 0
#water = 1
sediments = 1

i = 0
#Geometry["%05d_rect" % i] = Geom_Rect(sediments,0.0,Hsed/2.0,L,Hsed/2.0)
#i+=1
Geometry["%05d_line" % i] = Geom_Line(sediments,0.0,Hsed,"y","<",Grid.xmin,Grid.xmax)
#Geometry["%05d_sine" % i] = Geom_Sine(sediments,Hsed,Hsed/20.0,0.,L/5.0,"y","<",Grid.xmin,Grid.xmax)


#plt.axis([Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax])

#for key in Geometry:
#    Geometry[key].plot()

#plt.show()

#make dict of geometry
for key in Geometry:
   Geometry[key] = vars(Geometry[key])




##            Visualization
## =====================================
Visu.showParticles = False
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.height = 1 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "StrainRate"
Visu.writeImages = False
Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
#Visu.outputFolder = "/Users/abauville/GoogleDrive/Output/"
Visu.transparency = True



##              Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -1.0e-15
Numerics.CFL_fac = 0.2
Numerics.nLineSearch = 10
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 1
Numerics.maxNonLinearIter = 3

Numerics.absoluteTolerance = 1e-5

Numerics.etaMin = 1e-5



##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);
outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)


