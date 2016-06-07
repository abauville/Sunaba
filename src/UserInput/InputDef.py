# ============================================
#
# Classes of User Input
#
# Arthur Bauville, June 3, 2016
#
# ============================================
from math import *

class Frozen(object): # A metaclass that prevents the creation of new attributes
    __List = []
    def __setattr__(self, key, value):
        setIsOK = False
        for item in self.__List:
            if key == item:
                setIsOK = True
                
        if setIsOK == True:
            object.__setattr__(self, key, value)
        else:
            raise TypeError( "%r has no attributes %r" % (self, key) )



class Grid(Frozen):
    _Frozen__List = ["xmin","xmax","ymin","ymax","nxC","nyC"]
    def __init__(self):
        self.xmin   = -1.0
        self.xmax   =  1.0
        self.ymin   = -1.0
        self.ymax   =  1.0

        self.nxC    = 64
        self.nyC    = 64


     
class Numerics(Frozen):
    _Frozen__List = ["nTimeSteps", "nLineSearch", "maxNonLinearIter", "minNonLinearIter", "relativeTolerance", "absoluteTolerance","maxCorrection","CFL_fac","etaMin","etaMax"]
    def __init__(self):
        self.nTimeSteps  = 1 #  negative value for infinite
        self.nLineSearch = 1
        self.maxNonLinearIter = 1 # should always be greater than the number of line searches
        self.minNonLinearIter = 1 # should always be greater than the number of line searches
        self.relativeTolerance = 3E-5 # relative tolerance to the one of this time step
        self.absoluteTolerance = 3E-5 # relative tolerance to the first one of the simulation
        self.maxCorrection = 1.0

        self.CFL_fac = 0.5

        self.etaMin = 1E-4
        self.etaMax = 1E4

        


class Material(Frozen):
    _Frozen__List = ["name","material","n","cohesion","frictionAngle","rho0","eta0","alpha","beta","k","G"]
    def __init__(self):
            self.name = ""
            self.material = "Default"
            self.n = 1.0
            self.cohesion = 1E20
            self.frictionAngle = 30.0/180*pi
            self.rho0 = 1.0
            self.eta0 = 1.0
            
            self.alpha = 1.0
            self.beta = 1.0

            self.k = 1.0
            self.G = 1E20
            


class Particles(Frozen):
    _Frozen__List = ["nPCX","nPCY","noiseFactor"]
    def __init__(self):
        self.nPCX = 4
        self.nPCY = 4
        self.noiseFactor = 0.0



class Physics(Frozen):
    _Frozen__List = ["Cp","gx","gy"]
    def __init__(self):
        self.Cp = 1.0
        self.gx = 0.0
        self.gy = -1.0



class Visu(Frozen):
    _Frozen__List = ["type","typeParticles","showParticles","shiftFacX","shiftFacY","shiftFacZ","writeImages","transparency","alphaOnValue","showGlyphs","glyphType","glyphMeshType","glyphScale","glyphSamplingRateX","glyphSamplingRateY","width","height","outputFolder","retinaScale","particleMeshRes"]
    def __init__(self):
        self.type 	    = "StrainRate"; # Default
        self.typeParticles  = "Phase"; # Default
        self.showParticles  = True;
        self.shiftFacX      = 0.00;
        self.shiftFacY 	    = 0.00;
        self.shiftFacZ 	    = +.05;
        self.writeImages 	= False;
        self.transparency 	= False;
        self.alphaOnValue 	= False;
        self.showGlyphs 	= False;
        self.glyphType		= "StokesVelocity";
        self.glyphMeshType	= "ThinArrow";
        self.glyphScale		= 0.05;
        self.glyphSamplingRateX  = 3;
        self.glyphSamplingRateY  = 6;

        self.width              = 1024;
        self.height             = 1024;

        self.outputFolder = "../../OutputStokesFD"

        self.retinaScale = 2
        
        Visu.particleMeshRes = 6; # minimum 3, higher number make the particle voronoi diagram look smoother


class Char(Frozen):
    _Frozen__List = ["length","mass","time","temperature"]
    def __init__(self):
         self.length = 1.0
         self.mass = 1.0
         self.time = 1.0
         self.temperature = 1.0


class BCStokes(Frozen):
    _Frozen__List = ["backStrainRate","SetupType"]
    def __init__(self):
        self.backStrainRate = -1.0
        self.SetupType  = "PureShear"

class BCThermal(Frozen):
    _Frozen__List = ["TT","TB","SetupType"]
    def __init__(self):
        self.TT = 0
        self.TB = 0
        self.SetupType  = "PureShear"
        



