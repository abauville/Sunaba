#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD
import sys
import os
sys.path.insert(0, './../../src/UserInput')
#import json
#from InputDef import *
import InputDef as Input
import MaterialsDef as material
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *
from math import pi, sqrt, tan, sin, cos, exp, log
#print("\n"*5)





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
MPa     = 1e6
GPa     = 1e9

deg     = pi/180




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

## Description
## =====================================
Setup.Description = "Paper on Dynamic stress and Strain softening. Fig. dtDependence"

ProductionMode = False


Numerics.etaMin = 1e-5
Numerics.etaMax = 1e+5

##          Material properties
## =====================================
StickyAir   = Input.Material("StickyAir")
Matrix      = Input.Material("Sediments")
Inclusion   = Input.Material("Sediments")


Setup.MatProps = {"0":StickyAir,"1":Matrix,"2":Inclusion}

PhaseRef = Matrix#Inclusion
PhaseRef.isRef = True

Matrix.name     = "Sediment"
Inclusion.name  = "Inclusion"

Matrix.vDiff    = material.DiffusionCreep       ("Off")
Inclusion.vDiff = material.DiffusionCreep       ("Off")


Matrix.use_dtMaxwellLimit = True

Matrix.vDisl    = material.DislocationCreep     (eta0=5E23, n=1)
Inclusion.vDisl = material.DislocationCreep     (eta0=0.1*5E23, n=1)
StickyAir.vDiff = material.DiffusionCreep       (eta0=5E23/10000.0)


Matrix.rho0     = 0.0*2700  * kg/(m**3)
Inclusion.rho0  = 0.0*2700  * kg/(m**3)

Matrix.cohesion     = 50    * MPa
Inclusion.cohesion  = 50    * MPa

Matrix.frictionAngle    = 30 * deg
Inclusion.frictionAngle = 30 * deg

Matrix.G                = 1.0 * GPa
Inclusion.G             = 0.1*Matrix.G 
StickyAir.G             = Matrix.G/10000.0
StickyAir.cohesion      = Matrix.cohesion

#StickyAir.cohesion = Matrix.cohesion



##              Grid
## =====================================

if ProductionMode:
    RFac = 2
else:
    RFac = 1; # Resolution Factor
HFac = 1.0


H = HFac * 1 * km
HStickyAir = H*.28
Grid.ymin =  0.0
Grid.ymax =  H + HStickyAir
Grid.nyC = round(128*RFac)
dy = (Grid.ymax-Grid.ymin)/(Grid.nyC+1)



r = 6*dy# H/8.0         # inclusion radius
d = 2.0*r
theta = 33/180*pi # effective shear zone angle
#W = r*cos(45/180*pi) + (H-r*sin(45/180*pi))/tan(theta) # takes into account that the shear zone starts at 45 degree on the inclusion perimeter
#W = W*1.5
W = 192/100 * H



Grid.xmin = 0.0
Grid.xmax = W
#Grid.xmin = -W
#Grid.xmax = 0.0


Grid.nxC = round(128*W/(H+HStickyAir)*RFac)

dx = (Grid.xmax-Grid.xmin)/(Grid.nxC+1)




Grid.fixedBox = True

Physics.Pback = 100 * MPa


##              Numerics
## =====================================
Numerics.nTimeSteps = 1000000
BCStokes.backStrainRate = -1.0e-15
Numerics.dtAlphaCorr = 1.0
Numerics.CFL_fac_Stokes = 0.25
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 1
if ProductionMode:
    Numerics.maxNonLinearIter = 150
else: 
    Numerics.maxNonLinearIter = 1000

Numerics.absoluteTolerance = 1e-6
Numerics.relativeTolerance = 1e-3 # time current residual


Numerics.dtMaxwellFac_EP_ov_E  = 0.5;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = 0.0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = 0.5;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = False



Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.0



##              Geometry
## =====================================



inclusion_w = d/2
inclusion_h = d


slope = tan(0*pi/180)

i = 0
MatrixPhase = 1
Geometry["%05d_line" % i] = Input.Geom_Line(MatrixPhase,slope,H,"y","<",Grid.xmin,Grid.xmax)
InclusionPhase = 2
i+=1
Geometry["%05d_line" % i] = Input.Geom_Line(InclusionPhase,0.0,inclusion_w,"y","<",Grid.xmin,Grid.xmin+inclusion_w)



dt_stressFacList = [9.0e-2] # used only for the scaling
#[1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
for dt_stressFac in dt_stressFacList:    
    ##              Non Dim
    ## =====================================
    Char.length =  (Grid.xmax-Grid.xmin)/2
    Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
#    CharStress =    Matrix.cohesion*cos(Matrix.frictionAngle) + Physics.Pback *sin((Matrix.frictionAngle))
    
    #if ProductionMode:
    #    a_f = 100.0
    #else:    
    #    a_f = 100.0
    #    
    
        
    #DeltaSigma = CharStress*dt_stressFac ;
    G = Matrix.G
    EII = abs(BCStokes.backStrainRate)
    eta = Matrix.getRefVisc(0.0,Char.temperature,EII)
    
    S3 = Setup.Physics.Pback
    C = Matrix.cohesion
    phi = Matrix.frictionAngle
    S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
    Sy_back = (S1-S3)/2.0
    P_Lim = (S1+S3)/2.0
#    P = Setup.Physics.Pback
#    Sy_back = C*cos(phi) + P*sin(phi)
    RefTime  = eta/G * log(2*eta*EII / (2*eta*EII - Sy_back )); # time at which stress has built up to the 
    Char.time = RefTime*dt_stressFac
    
    CharVisc = 1.0/(1.0/eta+1.0/(G*Char.time))
    CharStress = CharVisc/Char.time
    
    Char.mass   = CharStress*Char.time*Char.time*Char.length
    
    Numerics.dtMin = 60 * s #Char.time #* 1e-12
    Numerics.dtMax = Char.time #* 1.0
    
    ####### !!!!!!!!!
    Numerics.dt_stressFac = dt_stressFac#0.05 # Used for the computation
    ####### !!!!!!!!!
    
    
    min_nTimeSteps = 100
    Numerics.maxTime        = max(RefTime*1.5 , min_nTimeSteps*Numerics.dtMin)
    dt = Numerics.dtMin
    nSteps = round(Numerics.maxTime/Numerics.dtMin)
    print("Fac = %.1e, maxTime= %.2f Myrs, nSteps= %i, CharStress*Char.time = %.2e"  % (dt_stressFac ,Numerics.maxTime/Myr, nSteps, CharStress*Char.time))
    print("eta = %.2e, eta_ve = %.2e, ratio_eta_Gdt = %.2e" % (eta, CharVisc, eta/G/dt))
    
    etaStickyAir = StickyAir.getRefVisc(0.0,Char.temperature,EII)
    print("etaStickyAir = %.2e, etaStickyAir_ve = %.2e, ratio_eta_Gdt = %.2e" % (etaStickyAir, 1.0/(1.0/etaStickyAir+1.0/(StickyAir.G*Char.time)), etaStickyAir/G/dt))
    
    etaInclusion = Inclusion.getRefVisc(0.0,Char.temperature,EII)
    print("etaInclusion = %.2e, etaInclusion_ve = %.2e, ratio_eta_Gdt = %.2e" % (etaInclusion, 1.0/(1.0/etaInclusion+1.0/(Inclusion.G*Char.time)), etaInclusion/G/dt))






    ###              Output
    ### =====================================
    #Output.folder = "/Users/abauville/Output_Paper_DynDecollement/DynStress_PureShear/nx_%i_ny_%i_G_%.2e_C_%.2e_fric_%.2e_Pref_%.2e" % (Grid.nxC, Grid.nyC, Matrix.G, Matrix.cohesion, Matrix.frictionAngle*180/pi, Physics.Pback)

    if sys.platform == "linux" or sys.platform == "linux2":
        # linux
        if ProductionMode:
            Output.folder = "/home/abauvill/StokesFD_Output/Paper_DynStress/Output/dtDependence/Production/dt_stressFac_%.1e" % Numerics.dt_stressFac      
        else:
            Output.folder = "/home/abauvill/StokesFD_Output/Paper_DynStress/Output/dtDependence/Test_NoAdv_Interp_adaptative/dt_stressFac_%.1e" % Numerics.dt_stressFac

    elif sys.platform == "darwin":
        # OS X
        if ProductionMode:
            Output.folder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Production/dt_stressFac_%.1e" % Numerics.dt_stressFac      
        else:
            Output.folder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_WeakInclusion_NoAdv_NoInterp_adaptative_NEW/dt_stressFac_%.1e" % Numerics.dt_stressFac



#elif sys.platform == "win32":
    # Windows...

   
    Output.strainRate = True
    Output.sigma_II = True
    Output.sigma_xx = True
    Output.sigma_xy_node = True
    Output.khi = True
    Output.P = True
#    Output.Z = True
    Output.strainRate = True
    
    Output.frequency = 1#timeFac
    
    
    
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
    Visu.outputFolder = "/Users/abauville/GoogleDrive/Output/"
    Visu.transparency = False
    
    Visu.showGlyphs =  False
    Visu.glyphMeshType = "Triangle"
    Visu.glyphScale = 250.0
    glyphSpacing = (Grid.ymax-Grid.ymin)/24 #50 * km
    Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
    Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))
    
    Visu.height = 1.0 * Visu.height
    Visu.width = 1* Visu.width
    
    #Visu.filter = "Linear"
    Visu.filter = "Nearest"
    
    Visu.shiftFacY = 0.0
    
    
    CharExtra = Input.CharExtra(Char)
    eta = Matrix.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
    #RefVisc = CharExtra.viscosity
    
    print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))
    
    RefP = PhaseRef.rho0*abs(Physics.gy)*H/2.0
    
    Visu.colorMap.Stress.scale  = 100.0e6/CharExtra.stress
    Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
    Visu.colorMap.Stress.max    = 1.0
    Visu.colorMap.Viscosity.scale = 1.0#RefVisc/CharExtra.visc
    Visu.colorMap.Viscosity.max = 4.0
    Visu.colorMap.EffectiveViscosity.max = 4.0
    Visu.colorMap.StrainRate.scale = 1.0#abs(BCStokes.backStrainRate/(1.0/Char.time))
    Visu.colorMap.StrainRate.max = 3.0
    Visu.colorMap.Velocity.scale = (abs(BCStokes.backStrainRate) * (Grid.xmax-Grid.xmin) ) / (Char.length/Char.time)
    Visu.colorMap.Velocity.center = 0.0
    Visu.colorMap.Velocity.max = 2.0
    Visu.colorMap.Velocity.log10on = True
    
    
    Visu.colorMap.Pressure.scale  = Physics.Pback/CharExtra.stress
    Visu.colorMap.Pressure.center = 1.0
    Visu.colorMap.Pressure.max    = 2.0
    
    
    Visu.colorMap.VelocityDiv.scale = 1e-1
    
    Visu.colorMap.Khi.max = 5.0
    Visu.colorMap.Khib.max = 5.0
    
    
    
    
    ###          Write the Input file
    ### =====================================
    os.system("mkdir " + Output.folder)
    os.system("cd " + Output.folder + "\n rm -r *")
    if (Visu.writeImages):
        os.system("mkdir " + Visu.outputFolder)
    Input.writeInputFile(Setup)
    print("platform = " + sys.platform)
    if sys.platform == "linux" or sys.platform == "linux2":
        # linux
        os.system("/home/abauvill/mySoftwares/StokesFD/ReleaseLinux/StokesFD ./../input.json")
        print("should have launched")
    elif sys.platform == "darwin":
        # OS X
        os.system("/Users/abauville/JAMSTEC/StokesFD/Release/StokesFD ./../input.json")
    #elif sys.platform == "win32":
        # Windows...
