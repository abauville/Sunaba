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

        #self.ySeg   = self.xmin
        #self.ySeg   = self.xmin

     
class Numerics(Frozen):
    _Frozen__List = ["nTimeSteps", "nLineSearch", "maxNonLinearIter", "minNonLinearIter", "relativeTolerance", "absoluteTolerance","maxCorrection","CFL_fac","etaMin","etaMax","dtMin","dtMax"]
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

        self.dtMin = 0.
        self.dtMax = 1e100

        


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
            
            self.alpha = 0.0
            self.beta = 0.0

            self.k = 1.0
            self.G = 1E20
            


class Particles(Frozen):
    _Frozen__List = ["nPCX","nPCY","noiseFactor","minPartPerCellFactor","maxPartPerCellFactor"]
    def __init__(self):
        self.nPCX = 4
        self.nPCY = 4
        self.noiseFactor = 0.0

        self.minPartPerCellFactor = 0.2
        self.maxPartPerCellFactor = 3.0



class Physics(Frozen):
    _Frozen__List = ["Cp","gx","gy"]
    def __init__(self,Dimensional):
        if Dimensional == True:
            self.Cp = 0.5
            self.gx = 0
            self.gy = -9.81
        else:
            self.Cp = 1.0
            self.gx = 0.0
            self.gy = -1.0
     



class Visu(Frozen):
    _Frozen__List = ["type","typeParticles","showParticles","shiftFacX","shiftFacY","shiftFacZ","writeImages","transparency","alphaOnValue","showGlyphs","glyphType","glyphMeshType","glyphScale","glyphSamplingRateX","glyphSamplingRateY","width","height","outputFolder","retinaScale","particleMeshRes","particleMeshSize","filter"]
    def __init__(self):
        self.type 	    = "StrainRate" # Default
        self.typeParticles  = "Phase" # Default
        self.showParticles  = True
        self.shiftFacX      = 0.00
        self.shiftFacY 	    = 0.00
        self.shiftFacZ 	    = -0.05
        self.writeImages 	= False
        self.transparency 	= False
        self.alphaOnValue 	= False
        self.showGlyphs 	= False
        self.glyphType		= "StokesVelocity"
        self.glyphMeshType	= "ThinArrow"
        self.glyphScale		= 0.05
        self.glyphSamplingRateX  = 3
        self.glyphSamplingRateY  = 6

        self.width              = 1024
        self.height             = 1024

        self.outputFolder = "../../OutputStokesFD"

        self.retinaScale = 2
        
        self.particleMeshRes = 6 # minimum 3, higher number make the particle voronoi diagram look smoother

        self.particleMeshSize = 0.05

        self.filter = "Linear"

class Char(Frozen):
    _Frozen__List = ["length","mass","time","temperature"]
    def __init__(self):
         self.length = 1.0
         self.mass = 1.0
         self.time = 1.0
         self.temperature = 1.0

    def set_based_on_strainrate(self,PhaseRef,BCStokes,BCThermal,Grid):
        self.time   = abs(1.0/BCStokes.backStrainRate)
        self.length = (Grid.xmax-Grid.xmin)/2

        CharStress = 2*PhaseRef.eta0
        self.mass   = CharStress*self.time*self.time*self.length
        self.temperature = (BCThermal.TB + BCThermal.TT)/2.0
        
        
    def set_based_on_lithostatic_pressure(self,PhaseRef,BCThermal,Physics,Grid):
        self.length = (Grid.ymax-Grid.ymin)/2.0
        CharStress  = PhaseRef.rho0*abs(Physics.gy)*self.length
        CharVisc    = PhaseRef.eta0

        self.time   = CharVisc/CharStress
        self.mass   = CharStress*self.time*self.time*self.length
        self.temperature = (BCThermal.TB + BCThermal.TT)/2.0
        

        
    



class BCStokes(Frozen):
    _Frozen__List = ["backStrainRate","SetupType"]
    def __init__(self):
        self.backStrainRate = -1.0
        self.SetupType  = "PureShear"

class BCThermal(Frozen):
    _Frozen__List = ["TT","TB","SetupType"]
    def __init__(self):
        self.TT = 1.0
        self.TB = 1.0
        self.SetupType  = "PureShear"
        


## Geometry ##
class Geom_Circle(object):
    def __init__(self,phase,cx,cy,radius):
        self.phase  = phase
        self.cx     = cx
        self.cy     = cy
        self.radius = radius

class Geom_Rect(object):
    def __init__(self,phase,llx,lly,width,height):
        self.phase = phase
        self.llx = llx
        self.lly = lly
        self.width = width
        self.height = height

class Geom_Line(object):
    def __init__(self,phase,a,b,definedFor,condition,Min,Max):
        self.phase = phase
        self.a = a
        self.b = b
        self.definedFor = definedFor # "y" or "x"
        self.condition = condition # ">", "<"
        self.min = Min
        self.max = Max


class Geom_Sine(object):
    # if definedFor "y":
    # y = base + amplitude*sin(wavelength*x*2*pi+ wavephase)
    def __init__(self,phase,base,amplitude,wavephase,wavelength,definedFor,condition,Min,Max):
        self.phase = phase
        self.amplitude = amplitude
        self.wavephase = base
        self.base  = base
        self.wavelength = wavelength
        self.definedFor = definedFor
        self.condition = condition
        self.min = Min
        self.max = Max
         
class Geom_Polygon(object):
    # x and y should be arrays of vertices
    # the test for polygon is not ready yet
    def __init__(self,phase,x,y):
        self.phase = phase
        self.x = x
        self.y = y
