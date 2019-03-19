#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD
import sys
import os
sys.path.insert(0, '../../src/UserInput')
sys.path.insert(0, '../../PostProcessing/Paper_Decollement/CriticalTaper')
#import json
#from InputDef import *
import InputDef as Input
import MaterialsDef as material
from CritTaper_utils import Taper
import numpy as np
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *
from math import pi, sqrt, tan, sin, cos, log10, log2, log
print("\n"*5)

##             Units
## =====================================
m       = 1.0
s       = 1.0
K       = 1.0
kg      = 1.0

cm      = 0.01      * m
km      = 1000.0    * m

mn      = 60        * s
hour    = 60        * mn
day     = 24        * hour
yr      = 365       * day
Myr     = 1e6       * yr

Pa      = kg/m/s/s
MPa     = 1e6       * Pa


##      Declare singleton objects
## =====================================
Setup = Input.Setup(isDimensional=True)
Grid = Setup.Grid
Numerics = Setup.Numerics
Particles = Setup.Particles
Physics = Setup.Physics
Visu = Setup.Visu
Char = Setup.Char
BCStokes = Setup.BC.Stokes
BCThermal = Setup.BC.Thermal
ICThermal = Setup.IC.Thermal
ICDarcy = Setup.IC.Darcy
MatProps = Setup.MatProps
Geometry = Setup.Geometry
Output = Setup.Output





##             Main Parameters
## =====================================
ProductionMode = False


Bottom_type = "inactive"
#Bottom_type = "fixed"
#Bottom_type = "weakenable"

Hc_nd = 1.0/2.0
#Hc_nd = 1.0/1.0
#Hc_nd = 1.0/8.0


Lambda = 0.9
#weakFac = 0.4
PfWeakFac = 0.00
frictionWeakFac = 0.0
cohesionWeakFac = 0.99
Lambda_b_Fac = 0.0


maxElasticStrain = 0.05



timeFac = .75
beta        =  0.0 * pi/180.0 # place holder




#alpha = 13.0*pi/180.0
#alpha = 25.0*pi/180.0

## ============= RefTaper =================    
rho_w = 1000.0
rho = 2500.0*(1.0-Lambda)
phiRef   = 30.0*pi/180.0
LambdaRef=Lambda

thisTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                 Lambda=LambdaRef, Lambda_b=LambdaRef,
                 rho_w=rho_w, rho=rho)
thisTaper.computeAlphaVsBeta(n=2010)

betaMinRef = np.min(thisTaper.beta_all)
betaMaxRef = np.max(thisTaper.beta_all)

alpha  = thisTaper.findAlpha(beta,"upper")
## ========================================

L = 9.0
Lwedge = L

Hwedge = 1.0#Lwedge * tan(alpha)

Htotal = Hwedge + 1.25
shFac = Hwedge*Lwedge/2.0  

print("Lambda = %.2f, alpha = %.2f deg, shFac = %.2f" % (Lambda, alpha*180.0/pi, shFac))


if ProductionMode:
    nGrid_H = 64
else:
    nGrid_H = 32
    
nGrid_H = 48

Setup.Description = "Hc = %.5e, Lambda = %.5e, weakFac = %.5e, Beta = %.5e, alpha = %.5e, shFac = %.5e, nGrid_H = %i" % (Hc_nd, Lambda, PfWeakFac, beta, alpha, shFac, nGrid_H)


localMachineIndex = 0 # 0: Mac, 1: Desktop Linux, 2: DA System
runMachineIndex = 2 # 0: Mac, 1: Desktop Linux, 2: DA System
if localMachineIndex==0:
    localPreBaseFolder = "/Users/abauville/Output/"
elif localMachineIndex==1:
    localPreBaseFolder = "/home/abauvill/Output/"
elif localMachineIndex==2:
    localPreBaseFolder = "/work/G10501/abauville/"
else:
    raise ValueError("Unknown localMachineIndex");
    
if runMachineIndex==0:
    runPreBaseFolder = "/Users/abauville/Output/"
elif runMachineIndex==1:
    runPreBaseFolder = "/home/abauvill/Output/"
elif runMachineIndex==2:
    runPreBaseFolder = "/work/G10501/abauville/"
else:
    raise ValueError("Unknown runMachineIndex");
    
##          Material properties
## =====================================
StickyAir   = Input.Material("StickyAir")
Sediment    = Input.Material("Sediments")
Basement    = Input.Material("Sediments")
Backstop   = Input.Material("Sediments")
WeakChannel   = Input.Material("Sediments")


Setup.MatProps = {"0":StickyAir,"1":Sediment,"2":Basement, "3":Backstop, "4":WeakChannel}

Numerics.stickyAirSwitchPhaseTo = 1
Numerics.stickyAirSwitchPassiveTo = 1


PhaseRef = Sediment
#PhaseRef = StickyAir
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Sediment.name = "Sediment"
Basement.name = "Basement"
Backstop.name = "Backstop"

Sediment.vDiff = material.DiffusionCreep       ("Off")
Basement.vDiff = material.DiffusionCreep       ("Off")
Backstop.vDiff = material.DiffusionCreep       ("Off")

StickyAir.rho0 = rho_w



Sediment.use_dtMaxwellLimit = False

Numerics.invariantComputationType = 1



# Static Fluid Pressure Factor
Sediment.staticPfFac = Lambda

# StrainWeakening
Sediment.staticPfFacWeakFac = PfWeakFac
Sediment.frictionAngleWeakFac = frictionWeakFac
Sediment.cohesionWeakFac = cohesionWeakFac#weakFac # 0.0 is none weak, 1.0 is fully weakened Cfinal = Cini*(1-CweakFac)

Sediment.strainWeakStart = 0.025
Sediment.strainWeakEnd = .2

    



Backstop.G = 5e8*100.0
WeakChannel.G  = 5e8
Basement.G  = Backstop.G*10.0



##              Numerics
## =====================================
Numerics.etaMin = 1e-8
Numerics.etaMax = 1e8
Numerics.nTimeSteps = 15000
Numerics.CFL_fac_Stokes = .25
#if weakFac>=0.6:
#    Numerics.CFL_fac_Stokes = .05
Numerics.CFL_fac_Darcy = 1000.0
Numerics.CFL_fac_Thermal = 10000.0
Numerics.nLineSearch = 1
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 20
Numerics.maxNonLinearIter = 20
#if Bottom_type!="inactive":
#    Numerics.maxNonLinearIter = 4
Numerics.dtAlphaCorr = .3
Numerics.absoluteTolerance = 1e-4
Numerics.relativeTolerance  = 1e-3


Numerics.dtMaxwellFac_EP_ov_E  = .5   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .5   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress


Numerics.stressSubGridDiffFac = 1.0

Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.5


Numerics.yieldComputationType = 1


## Main parameters for this setup
## =====================================

Sediment.frictionAngle  = 30.0/180.0*pi
Backstop.frictionAngle = 30.0/180.0*pi
Basement.frictionAngle  = Sediment.frictionAngle
WeakChannel.frictionAngle  = 30.0/180.0*pi

Backstop.cohesion = 50000.0e6
WeakChannel.cohesion = 5.0e6
Basement.cohesion = 50*1e6


Sediment.rho0 = rho


HFac        = 2.0
Hsed        = HFac*1.0e3


Hc = Hc_nd*Hsed
g = 9.81
C = Hc*((1.0-Lambda)*Sediment.rho0*g*tan(Sediment.frictionAngle))

#Lambda_hydro = 0.4
#Lambda_ov = (Lambda-Lambda_hydro)/(1.0-Lambda_hydro)
#rhoWater = 1000.0
#rho_e = Sediment.rho0-rhoWater
#C = Hc*((1.0-Lambda_ov)*rho_e*g*tan(Sediment.frictionAngle))

Sediment.cohesion = C
StickyAir.cohesion = 1.0*Sediment.cohesion


BCStokes.Bottom_frictionAngle = Sediment.frictionAngle
BCStokes.Bottom_cohesion = Sediment.cohesion
BCStokes.Bottom_staticPfFac = Lambda + Lambda_b_Fac * (1.0-Lambda)
BCStokes.Bottom_type = Bottom_type # values can be: "inactive", "fixed", "weakenable"


print("cohesion = %.4f MPa" % (Sediment.cohesion/1e6))

LWRatio = L/Htotal



Grid.xmin = -Htotal*Hsed*LWRatio
Grid.ymax =  Htotal*Hsed

Grid.xmax = 0.0e3
Grid.ymin = 0.0e3

#if ProductionMode:
Grid.nxC = round(Htotal*nGrid_H*LWRatio) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = round(Htotal*nGrid_H)#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)
    
Grid.fixedBox = True

print("Grid.nxC = %i, Grid.nyC = %i" % (Grid.nxC, Grid.nyC))

VatBound = - 10.0 * cm/yr
dx = (Grid.xmax-Grid.xmin)/Grid.nxC
dy = (Grid.ymax-Grid.ymin)/Grid.nyC
BCStokes.backStrainRate = VatBound / (Grid.xmax-Grid.xmin)

##Numerics.maxTime = shFac*Hsed/abs(VatBound)
#alpha = 5.0*np.pi/180.0
#Lwedge = 17.0
#Numerics.maxTime = (Lwedge*Hsed)**2*np.tan(alpha)/(2.0*np.abs(VatBound)*(Hwedge*Hsed)) # time necessary to create a wedge of length Lwedge and of angle alpha


Plitho = Sediment.rho0 * abs(Physics.gy) * 1.0*Hsed
Sigma_y = Sediment.cohesion*cos(Sediment.frictionAngle) + sin(Sediment.frictionAngle)*(1.0-Lambda)*Plitho
print("RefViscBrittle = %.2e Pa.s" % (Sigma_y/abs(BCStokes.backStrainRate)))
print("backStrainRate = %.2e, Sigma_y = %.2e MPa" % (BCStokes.backStrainRate, Sigma_y/1e6))


RefVisc =  (Sigma_y/abs(BCStokes.backStrainRate))
RefViscSurf =  (Sediment.cohesion/abs(BCStokes.backStrainRate))

RefVisc *= 1
StickyAir.vDiff = material.DiffusionCreep(eta0=RefViscSurf/10000.0)
Sediment.vDisl = material.DislocationCreep     (eta0=RefVisc*100.0, n=1)
Backstop.vDisl = material.DislocationCreep    (eta0=RefVisc*1, n=1)
Basement.vDisl = material.DislocationCreep     (eta0=RefVisc*100, n=1)



BoxTilt = -beta
Physics.gx = -9.81*sin(BoxTilt);
Physics.gy = -9.81*cos(BoxTilt);








##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin

Hbase = 0


Wseamount = .15e3*HFac
xseamount = Grid.xmin + 1e3

i = 0
AirPhase = 0
SedPhase = 1
BasementPhase = 2
BackPhase = 3
ChannelPhase = 4

Lweak = Grid.xmax-Grid.xmin
Hweak = .24e3*HFac
ThickWeak = .05e3*HFac



HHorst = Hsed*0.25
WHorst = 2.0*HHorst

Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,0.0,Hsed,"y","<",Grid.xmin,Grid.xmax)

    
    
i+=1
slope = 0.0*pi/180.0
Lwedge = 4.0
Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hsed - slope*(L-Lwedge)*Hsed,"y","<",Grid.xmin,Grid.xmax)

HSFac = 1
BCStokes.Sandbox_TopSeg00 = Hbase + 0*Hbase + 0*dy + 0*HSFac*dy
BCStokes.Sandbox_TopSeg01 = BCStokes.Sandbox_TopSeg00+HSFac*dy#0.405e3*HFac





##                 BC
## =====================================
BCStokes.SetupType = "Sandbox_InternalBC"
BCStokes.SetupType = "Sandbox"

BCStokes.Sandbox_NoSlipWall = True


##                 IC
## =====================================
ICThermal.age = 100*Myr

ICDarcy.background = 0.0#Numerics.phiMin
ICDarcy.Amp = 0.0
ICDarcy.xc = Grid.xmin+(Grid.xmax-Grid.xmin)/2
ICDarcy.yc = Grid.ymin+(Grid.ymax-Grid.ymin)/2
ICDarcy.wx = (Grid.xmax-Grid.xmin)/24.0
ICDarcy.wy = (Grid.xmax-Grid.xmin)/24.0

##              Non Dim
## =====================================
SedVisc = Sediment.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
BaseVisc = Basement.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))


##              Non Dim
## =====================================
Char.length =  Hsed#(Grid.xmax-Grid.xmin)/2
Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0


EII = abs(VatBound) / dx
eta = Sediment.getRefVisc(0.0,Char.temperature,EII)

S3 = Setup.Physics.Pback
C = Sediment.cohesion
phi = Sediment.frictionAngle


Sy_back = ( C*cos(phi) + (1.0-Lambda)*Plitho*sin(phi) ) / (1.0-sin(phi))

#myRefTime = 4*yr

Sediment.G = Sy_back/2.0 / maxElasticStrain # Choose G such that the strain needing to reach the yield is a given value (e.g. 0.5%)
G = Sediment.G
StickyAir.G = Sediment.G*1.0
RefTime  =  eta/G * log(2.0*eta*EII / (2.0*eta*EII - Sy_back )); # time at which stress has built up to the 
print("Gref = %.2e Pa, C = %.2e MPa, RefTime = %.2f yr\n" % (Sediment.G , Sediment.cohesion/MPa, RefTime/yr))


Visu.colorMap.Strain.max=Sediment.strainWeakEnd/2.0


Char.time = RefTime

   
  
#CharVisc = 1.0/(1.0/eta+1.0/(G*Char.time))
#CharStress = CharVisc/Char.time
CharStress = Plitho

Char.mass   = CharStress*Char.time*Char.time*Char.length

   
Numerics.dtIni = timeFac*RefTime
Numerics.dtMin = timeFac*RefTime
Numerics.dtMax = timeFac*RefTime

#
#if Bottom_type!="inactive":
#    Numerics.dtIni = 10.0*RefTime
#    Numerics.dtMin = 10.0*RefTime
#    Numerics.dtMax = 10.0*RefTime



###              Output
### =====================================

postBaseFolder = "ListricDecollement/Test_local/Lambda%03d_Hc%03d_Weak%03d_GFac%03d/" % (Lambda*100, Hc_nd*100, cohesionWeakFac*100, maxElasticStrain*100)

baseFolder = localPreBaseFolder + postBaseFolder


Output.folder = (runPreBaseFolder + postBaseFolder + "Output/" )
if ProductionMode: 
    ResFac = 0
    Output.folder = (runPreBaseFolder + postBaseFolder + "Output/" )

    Output.strain     = True
    Output.phase = True

    Output.particles_pos = True            
    Output.particles_strain   = True

    
    Output.particles_posIni = True
    Output.timeFrequency = Numerics.dtMax*800.0
    
#    if Bottom_type!="inactive":
#        Output.timeFrequency = RefTime*1600.0
    
    if Lambda==0.6 and beta==0.0 and Bottom_type=="inactive":
        Output.timeFrequency = Numerics.dtMax*200.0
        Output.P = True
        Output.sigma_xx = True
        Output.sigma_xy = True

        Output.Vx = True
        Output.Vy = True
        
    #
    Output.breakpointRealTimeFrequency = 24.0*hour
    Output.restartAfterBreakpoint = True


##              Info
## =====================================
print("\n"*5)
CharExtra = Input.CharExtra(Char)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

CharVisc2 = 1.0/(1.0/StickyAirVisc+1.0/(StickyAir.G*Char.time))


print("RefVisc = %.2e" % RefVisc)
print("Sediment Visc = %.2e" % SedVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)
print("BaseVisc = %.2e" %  BaseVisc)

print("Lc = " + str(  (Sediment.cohesion*cos(Sediment.frictionAngle)) / (Sediment.rho0 * abs(Physics.gy) * sin(Sediment.frictionAngle))  ))




##            Visualization
## =====================================

Particles.passiveGeom = "Grid_w_Layers"

Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/10 / (Grid.ymax/Hsed)
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.shaderFolder = "../Shaders/Sandbox_w_Layers_Backstop" # Relative path from the running folder (of StokesFD)
#Visu.shaderFolder = "../Shaders/Default" # Relative path from the running folder (of StokesFD)

#Visu.type = "StrainRate"
Visu.type = "StrainRate"
#Visu.renderFrequency = 10#round(128*yr/Numerics.dtMin)
#        Visu.renderTimeFrequency = 32*yr
#Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/StokesFD_Output/Test_NewRotation"
#Visu.outputFolder = ("/Users/abauville/Output/Sandbox_NumericalConvergenceTest_NewRHS/dt_%.0fyr/ResFac_%.1f" % (Numerics.dtMin/yr, ResFac) )
Visu.outputFolder = (baseFolder + "Visu/")
Visu.transparency = True

Visu.glyphMeshType = "TensorCross"
Visu.glyphType = "DeviatoricStressTensor"
#Visu.showGlyphs = True
#Visu.glyphScale = 8.0/(abs(VatBound)/(Char.length/Char.time))
Visu.glyphScale = 15.0
Visu.glyphSamplingRateX = nGrid_H/4.0
Visu.glyphSamplingRateY = nGrid_H/4.0
#glyphSpacing = (Grid.ymax-Grid.ymin)/64 #50 * km
#Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
#Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = .75 * Visu.height
Visu.width = 1.25 * Visu.width


Visu.filter = "Nearest"

Visu.shiftFacY = -0.5
Visu.shiftFacZ = -0.1

print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0
Visu.colorMap.Stress.scale  = Sediment.cohesion/CharExtra.stress
Visu.colorMap.Stress.center = 0.0
Visu.colorMap.Stress.max    = 1.00
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5
Visu.colorMap.Porosity.max = 1.0

Visu.colorMap.Pressure.scale  = .5*Plitho/CharExtra.stress
Visu.colorMap.Pressure.center = 2.0
Visu.colorMap.Pressure.max    = 4.00


Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 10.0

Visu.colorMap.Velocity.log10on = True
Visu.colorMap.Velocity.scale = (10.0*cm/yr) / (Char.length/Char.time)#abs(VatBound) / (Char.length/Char.time)
Visu.colorMap.Velocity.center = 0
Visu.colorMap.Velocity.max = 2.0

Visu.colorMap.Vorticity.max = 0.0005/yr /  (1.0/Char.time) # in rad/yr



Visu.colorMap.VelocityDiv.scale = 1e-4



Visu.colorMap.POvPlitho.log10on = True
Visu.colorMap.POvPlitho.center = 0.0
Visu.colorMap.POvPlitho.max = log10(2.0)

Visu.closeAtTheEndOfSimulation = False

Visu.colorMap.VxRes.scale = 1e-6
Visu.colorMap.VyRes.scale = 1e-6
Visu.colorMap.PRes.scale = 1e-6

###          Write the Input file
### =====================================
Input.writeInputFile(Setup)

#
if Visu.writeImages or Output.write():
    os.system("mkdir " + baseFolder )
    os.system("mkdir " + baseFolder + "Output")
    os.system("mkdir " + baseFolder + "Visu")
    os.system("mkdir " + baseFolder + "Input")    
    os.system("cp ../input.json " + baseFolder + "Input/input.json")

#if Output.breakpointFrequency > 0:
os.system("mkdir " + baseFolder + "Breakpoint")

if Bottom_type=='inactive':
    outJobFile = 'B%02d_W%02d_L%02d' % (int(beta*180.0/pi*10.0), int(Sediment.staticPfFacWeakFac*100),int(Lambda*100))
elif Bottom_type=='fixed':
    outJobFile = 'BotFixed_B%02d_W%02d_L%02d_Lb%02d' % (int(beta*180.0/pi*10.0), int(Sediment.staticPfFacWeakFac*100),int(Lambda*100),int(Lambda_b_Fac*100))
elif Bottom_type=='weakenable':
    outJobFile = 'BotWeak_B%02d_W%02d_L%02d_Lb%02d' % (int(beta*180.0/pi*10.0), int(Sediment.staticPfFacWeakFac*100),int(Lambda*100),int(Lambda_b_Fac*100))
#    os.system("mkdir " + Visu.outputFolder)
#    os.system("mkdir " + Output.folder)
if Lambda>0.6:
    memsize = 6
else:
    memsize = 8    
restartNumber = -1
# Write a job submission file
JobFileContent = """#!/bin/csh
#PBS -q l                                       # batch queue 
#PBS -b 1                                       # Number of jobs per request 
#PBS -r n                                       # rerunning disable
#PBS -l elapstim_req=26:00:00                   # Elapsed time per request
#PBS -l cpunum_job=4                            # Number of CPU cores per job
#PBS -l memsz_job=%igb                          # Memory size per job
#PBS -v OMP_NUM_THREADS=4                       # Number of threads per process
#PBS -o /home/G10501/abauville/Jobs/%s.o.%%s.%%j                              # standard output to outJobFileName.<reqID>.<jobNo>
#PBS -e /home/G10501/abauville/Jobs/%s.e.%%s.%%j                              # standard error to outJobFileName.<reqID>.<jobNo>
/work/G10501/abauville/Software/StokesFD/ReleaseDA/StokesFD /work/G10501/abauville/%s/input.json %05d""" % (memsize, outJobFile, outJobFile, postBaseFolder + "Input", restartNumber)
   

file = open(baseFolder + "Input/job.sh","w") 
file.write(JobFileContent)
file.close()
file = open(baseFolder + "Input/jobRestart.sh","w") 
file.write(JobFileContent)
file.close()
#os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ./input.json")
