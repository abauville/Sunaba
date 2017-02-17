#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
#import json
#from InputDef import *
import InputDef as Input
import MaterialsDef as material
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *
from math import pi, sqrt, tan, sin, cos

import os


from subprocess import call


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

Pa      = 1.0
MPa     = 1e6       * Pa

degree = pi/180




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








##                 BC
## =====================================
BCStokes.SetupType = "Sandbox"

BCStokes.Sandbox_NoSlipWall = False







##              Numerics
## =====================================


Numerics.phiMin = 1e-5
Numerics.phiMax = 0.9

Numerics.etaMin = 1e-6

Numerics.nTimeSteps = -10
BCStokes.backStrainRate = -1.0e-14
Numerics.CFL_fac_Stokes = 0.45
Numerics.CFL_fac_Darcy = 0.5
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 1
Numerics.maxNonLinearIter = 2
Numerics.absoluteTolerance = 1e-6





Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.9





##          Material properties
## =====================================
StickyAir   = Input.Material("StickyAir")
Sediment    = Input.Material("Sediments")
Basement    = Input.Material("Sediments")

Setup.MatProps = {"0":StickyAir,"1":Sediment,"2":Basement}

PhaseRef = StickyAir
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Sediment.name = "Sediment"
Basement.name = "Basement"

Sediment.vDiff = material.DiffusionCreep       ("Off")
Basement.vDiff = material.DiffusionCreep       ("Off")
#Basement.vDiff = material.DiffusionCreep       (eta0 = 1e23)

Sediment.vDisl = material.DislocationCreep     (eta0=1E90, n=10)
Basement.vDisl = material.DislocationCreep     (eta0=1E150, n=10)

Sediment.vDisl = material.DislocationCreep     (eta0=5E26, n=1)
Basement.vDisl = material.DislocationCreep     (eta0=5E29, n=1)


Basement.cohesion = 100 * MPa
Basement.frictionAngle = 50 * degree

#StickyAir.rho0 = 1000.0
StickyAir.rho0 = 0000.00


Sediment.G = 1e10
Basement.G = 1e12
StickyAir.G = 1e10

StickyAir.vDiff = material.DiffusionCreep(eta0=1E17)




             


##                   Systematics loop - creating the input files
## ============================================================================


#Syst_HFac               = [0.001, 0.01, 0.1, 1.0, 10.0]
#Syst_surfaceAngleAngle  = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30] # in degrees
#Syst_frictionAngle      = [4, 8, 16, 32]    # in degrees
#Syst_cohesion           = [0.001*MPa, 0.01*MPa, 0.1*MPa, 1.0*MPa, 10.0*MPa, 100.0*MPa]

Syst_HFac               = [1.0]
Syst_surfaceAngleAngle  = [0] # in degrees
Syst_frictionAngle      = [16]    # in degrees
Syst_cohesion           = [1.0*MPa]


no_simulation = len(Syst_HFac) * len(Syst_surfaceAngleAngle) * len(Syst_frictionAngle) * len(Syst_cohesion)
              
Visu.height = round (0.5* Visu.height)
Visu.width = 1* Visu.width
                
                
print("number of simulations: " + str(no_simulation))

for thisHFac in Syst_HFac:
    for thisSurfAngle in Syst_surfaceAngleAngle:
        for thisCohesion in Syst_cohesion:
            for thisFrictionAngle in Syst_frictionAngle:

                Sediment.frictionAngle = thisFrictionAngle   * degree
                Sediment.cohesion = thisCohesion
                
                surface_angle = thisSurfAngle * degree
                slope = tan(surface_angle)
                
                
                
                
                HFac = thisHFac
                
                
                
                
                
                
                StickyAir.cohesion = 1.0*Sediment.cohesion/5.0
                
                
                
                
                ##              Grid
                ## =====================================
                
                LWRatio = 7
                
                Grid.xmin = HFac* -(1.28+.32) * km *LWRatio
                Grid.xmax = HFac*  0.0e3
                Grid.ymin = HFac* 0.0e3
                Grid.ymax = HFac* (1.28+.96) * km * 1.0
                Grid.nxC = round(1/1*(128*LWRatio)) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
                Grid.nyC = round(1/1*(128))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)
                
                Grid.fixedBox = True
                
                #print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))
                ##              Non Dim
                ## =====================================
                
                L = (Grid.ymax-Grid.ymin)/2.0
                
                Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)
                
                
                
                
                
                
                
                
                
                ##              Geometry
                ## =====================================
                
                W = Grid.xmax-Grid.xmin
                H = Grid.ymax-Grid.ymin
                Hsed = HFac*1.0e3
                Hbase = HFac*0.1e3
                
                
                
                i = 0
                SedPhase = 1
                BasementPhase = 2
                Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hsed - slope*W,"y","<",Grid.xmin,Grid.xmax)
                #i+=1
                #Geometry["%05d_line" % i] = Input.Geom_Line(BasementPhase,slope,Hbase - slope*W,"y","<",Grid.xmin,Grid.xmax)
                
                
                BCStokes.Sandbox_TopSeg00 = 0.45e3*HFac
                BCStokes.Sandbox_TopSeg01 = 0.50e3*HFac
                
                #i+=1
                #Geometry["%05d_line" % i] = Input.Geom_Line(BasementPhase,0.0,Grid.xmax-Hbase,"x",">",BCStokes.Sandbox_TopSeg01,Hsed)
                
                
                ##              Output
                ## =====================================
                Output.folder = "/Users/abauville/StokesFD_Output/FaultAngle_Systematics_Batch02/" + "HSed" + str(round(thisHFac*1000)) + "_SA" + str(round(thisSurfAngle)) + "_C" + str(round(thisCohesion)) + "_FA" + str(round(thisFrictionAngle))
#                Output.khi          = True
#                Output.strainRate   = True
#                Output.sigma_xx     = True
#                Output.sigma_xy     = True
#                Output.P            = True
#                #Output.phase        = True
#                Output.frequency    = Numerics.nTimeSteps-1
#                Output.saveFirstStep = False
                
                print(Output.folder)
                
                ## Description
                ## =====================================
                Setup.Description = str({'surface_angle':surface_angle, 'Hsed':Hsed, 'friction_angle':Sediment.frictionAngle, 'cohesion':Sediment.cohesion})
                
                
                
                
                
                
                
                
                
                
                
                
                
                             
                
                ##            Visualization
                ## =====================================
                Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/16
                Particles.passiveDx = Particles.passiveDy
                
                Visu.showParticles = False
                Visu.filter = "Nearest"
                Visu.particleMeshRes = 6
                Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC
                
                
                Visu.type = "StrainRate"
                Visu.writeImages = False
                #Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
                Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_SandboxNew/"
                Visu.transparency = True
                
                Visu.showGlyphs = False
                Visu.glyphMeshType = "Triangle"
                Visu.glyphScale = 0.1/(BCStokes.refValue/(Char.length/Char.time))
                glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
                Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
                Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))
                
                
                #Visu.filter = "Linear"
                Visu.filter = "Nearest"
                
                Visu.shiftFacY = -0.0
                
                
                print("\n"*5)
                CharExtra = Input.CharExtra(Char)
                RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
                SedVisc = Sediment.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
                BaseVisc = Basement.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
                
                #StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/1000.0)
                
                StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
                
                print("RefVisc = %.2e" % RefVisc)
                print("Sediment Visc = %.2e" % SedVisc)
                print("StickyAirVisc = %.2e" % StickyAirVisc)
                print("BaseVisc = %.2e" %  BaseVisc)
                
                
                
                
                RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0
                
                Visu.colorMap.Stress.scale  = 20.0e6/CharExtra.stress
                Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
                Visu.colorMap.Stress.max    = 1.0
                Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
                Visu.colorMap.Viscosity.max = 4.0
                Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
                Visu.colorMap.StrainRate.max = 1.5
                
                
                
                Visu.colorMap.Pressure.scale  = 250e6/CharExtra.stress
                Visu.colorMap.Pressure.center = 0.0
                Visu.colorMap.Pressure.max    = 1.00
                Visu.colorMap.CompactionPressure.scale  = 2e6/CharExtra.stress
                Visu.colorMap.CompactionPressure.center = 0.0
                Visu.colorMap.CompactionPressure.max    = 1.0
                Visu.colorMap.FluidPressure.scale  = 10e6/CharExtra.stress
                Visu.colorMap.FluidPressure.center = 0.0
                Visu.colorMap.FluidPressure.max    = 1.00
                
                
                Visu.colorMap.VelocityDiv.scale = 1e-1
                
                Visu.colorMap.Khi.max = 5.0
                Visu.colorMap.Khib.max = 5.0

             
                
                
                
                
                
                
                
                
                
                
                #print(Setup.Description)
                
                ###          Write the Input file
                ### =====================================
                Input.writeInputFile(Setup)
                #os.system("mkdir " + Output.folder)
                #os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ../input.json")
                
                
                
                
