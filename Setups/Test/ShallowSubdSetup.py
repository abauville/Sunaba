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
Phase2   = Material("Sediments") # oc_crust
Phase3   = Material("Sediments") # lith mantle
Phase4   = Material("Sediments") # asth mantle
Phase5   = Material("Sediments") # cont crust
PhaseRef = Phase1

Phase2.name = "Basement"

#Phase0.eta0 = 1e19
#Phase0.G    = 1e10
Phase0.rho0 = 1000.0


Phase2.rho0 = 3000.0
Phase3.rho0 = 3000.0
Phase4.rho0 = 2950.0
Phase5.rho0 = 2600.0


#Phase1.frictionAngle = 5*pi/180
Phase0.eta0 = 1e19
Phase1.eta0 = 1e21
Phase2.eta0 = 1e22
Phase3.eta0 = 1e24
Phase4.eta0 = 1e20
Phase5.eta0 = 1e19

Phase0.G    = 1e10
Phase1.G    = 1e10
Phase2.G    = 1e10
Phase3.G    = 1e10
Phase4.G    = 1e10
Phase5.G    = 1e10



Backphi = 0.0001
RefPerm = 1e-19
Phase0.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))
RefPerm = 1e-19
Phase1.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))

RefPerm = 1e-21
Phase2.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))


Phase1.isRef = True

MatProps = {'0': Phase0.__dict__, '1': Phase1.__dict__, '2': Phase2.__dict__, '3': Phase3.__dict__, '4': Phase4.__dict__, '5': Phase5.__dict__}



##              Some info
## =====================================
Backphi = 0.001
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
#BCStokes.SetupType = "SandBox"
#BCThermal.SetupType = "SandBox"









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

Grid.xmin = -1055.0e3
Grid.xmax =  1055.0e3
Grid.ymin =  -200.0e3
Grid.ymax =  35.0e3;
Grid.nxC = 1024#round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 512#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.fixedBox = True

##              Non Dim
## =====================================
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)
Char.set_based_on_lithostatic_pressure(PhaseRef,BCThermal,Physics,Grid)






##              Geometry
## =====================================
H = Grid.ymax-Grid.ymin
L = Grid.xmax-Grid.xmin
Hsed = 1e3
Hcrust = 5e3
Hcontcrust = 30e3
inc_cy = -25e3

Physics.y_oceanSurface = Hsed+7e3
#DepthWater = H/2.0
#TopWater = Hsed+DepthWater

A = 1.0*2.0*Hsed/10.0

Leff = L#-L/15.0

air = 0
#water = 1
sediments = 1
oc_crust = 2
oc_mantle = 3
asth_mantle = 4
concrust = 5

i = 0
#Geometry["%05d_rect" % i] = Geom_Rect(sediments,0.0,Hsed/2.0,L,Hsed/2.0)

Geometry["%05d_line" % i] = Geom_Line(sediments,0.0,Hsed,"y","<",Grid.xmin,Grid.xmax)
 
i+=1
Geometry["%05d_line" % i] = Geom_Line(oc_crust,0.0,0.0,"y","<",0.0,Grid.xmax)

i+=1
Geometry["%05d_line" % i] = Geom_Line(concrust,0.0,0.0,"y","<",Grid.xmin,0.0)

i+=1
Geometry["%05d_line" % i] = Geom_Line(oc_mantle,0.0,-Hcrust,"y","<",0.0,Grid.xmax)

i+=1
Geometry["%05d_line" % i] = Geom_Line(oc_mantle,0.0,-Hcontcrust,"y","<",Grid.xmin,0.0)

i+=1
Geometry["%05d_line" % i] = Geom_Line(asth_mantle,0.0,-50.0e3,"y","<",0.0,Grid.xmax)

i+=1
Geometry["%05d_line" % i] = Geom_Line(asth_mantle,0.0,-100.0e3,"y","<",Grid.xmin,0.0)

#i+=1
#Geometry["%05d_circle" % i] = Geom_Circle(inclusion,0.0, inc_cy, Hcrust/2.0)


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

Visu.height = 1/3 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "SIIOvYield"
Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output/"
Visu.transparency = True



##              Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -1.0e-15
Numerics.CFL_fac = 0.1
Numerics.nLineSearch = 10
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 1
Numerics.maxNonLinearIter = 10

Numerics.absoluteTolerance = 5e-4

Numerics.etaMin = 1e-5



##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);
outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)


