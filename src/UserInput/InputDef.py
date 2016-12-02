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
    _Frozen__List = ["xmin","xmax","ymin","ymax","nxC","nyC","fixedBox"]
    def __init__(self):
        self.xmin   = -1.0
        self.xmax   =  1.0
        self.ymin   = -1.0
        self.ymax   =  1.0

        self.nxC    = 64
        self.nyC    = 64

        #self.ySeg   = self.xmin
        #self.ySeg   = self.xmin
    
        self.fixedBox = False
class Numerics(Frozen):
    _Frozen__List = ["nTimeSteps", "nLineSearch", "maxNonLinearIter", "minNonLinearIter", "relativeTolerance", "absoluteTolerance","maxCorrection","CFL_fac_Stokes","CFL_fac_Thermal","CFL_fac_Darcy","etaMin","etaMax","dtMin","dtMax"]
    def __init__(self):
        self.nTimeSteps  = 1 #  negative value for infinite
        self.nLineSearch = 1
        self.maxNonLinearIter = 1 # should always be greater than the number of line searches
        self.minNonLinearIter = 1 # should always be greater than the number of line searches
        self.relativeTolerance = 3E-5 # relative tolerance to the one of this time step
        self.absoluteTolerance = 3E-5 # relative tolerance to the first one of the simulation
        self.maxCorrection = 1.0

        self.CFL_fac_Stokes  = 0.75
        self.CFL_fac_Thermal = 10.0
        self.CFL_fac_Darcy   = 0.5

        self.etaMin = 1E-4
        self.etaMax = 1E4

        self.dtMin = 0.
        self.dtMax = 1e100




class Material(Frozen):
    _Frozen__List = ["name","material","n","cohesion","frictionAngle","rho0","eta0","alpha","beta","k","G","perm0","eta_b","B","isAir","isWater", "isRef"]
    def __init__(self,material="Default",name=""):
        self.isRef    = False
        if material == "Default":
            self.name = name
            self.material = "Default"
            self.n = 1.0
            self.cohesion = 1E100
            self.frictionAngle = 30.0/180*pi
            self.rho0 = 1.0
            self.eta0 = 1.0

            self.alpha = 0.01
            self.beta = 0.01

            self.k = 1.0
            self.G = 1E100

            self.perm0  = 0.0001
            self.eta_b  = 1.0
            self.B      = 1E20

            self.isAir = False
            self.isWater = False
            

        elif material == "StickyAir":
            self.name = name
            self.material = "StickyAir"
            self.n = 1.0
            self.cohesion = 1E100
            self.frictionAngle = 30.0/180*pi
            self.rho0 = 1.0
            self.eta0 = 1E17

            self.alpha = 1E-5
            self.beta  = 1E-11

            self.k = 3.0
            self.G = 1E100

            self.perm0  = 1E-5
            self.eta_b  = 1E25
            self.B      = 1E23

            self.isAir = True
            self.isWater = False


        elif material == "StickyWater":
            self.name = name
            self.material = "StickyWater"
            self.n = 1.0
            self.cohesion = 1E100
            self.frictionAngle = 30.0/180*pi
            self.rho0 = 1000.0
            self.eta0 = 1E17

            self.alpha = 1E-5
            self.beta  = 1E-11

            self.k = 3.0
            self.G = 1E100

            self.perm0  = 1E-5
            self.eta_b  = 1E25
            self.B      = 1E23

            self.isAir = False
            self.isWater = True

        elif material == "Sediments":
            self.name = name
            self.material = "Sediments"
            self.n = 1.0
            self.cohesion = 10E6
            self.frictionAngle = 30.0/180*pi
            self.rho0 = 2500.0
            self.eta0 = 1E21

            self.alpha = 1E-5
            self.beta  = 1E-11

            self.k = 3.0
            self.G = 1E11

            self.perm0  = 5E-9
            self.eta_b  = 1E23
            self.B      = 1E20

            self.isAir = False
            self.isWater = False
















class Particles(Frozen):
    _Frozen__List = ["nPCX","nPCY","noiseFactor","minPartPerCellFactor","maxPartPerCellFactor","passiveGeom","passiveRes"]
    def __init__(self):
        self.nPCX = 4
        self.nPCY = 4
        self.noiseFactor = 0.0

        self.minPartPerCellFactor = 0.2
        self.maxPartPerCellFactor = 3.0

        self.passiveGeom = "Grid"
        self.passiveRes = 1.0/8.0

class Physics(Frozen):
    _Frozen__List = ["Cp","gx","gy","eta_f","rho_f","y_oceanSurface"]
    def __init__(self,Dimensional):
        if Dimensional == True:
            self.Cp = 1000.0
            self.gx = 0.0
            self.gy = -9.81
            self.eta_f = 100.0
            self.rho_f = 1000.0
            self.y_oceanSurface = 0.0
        else:
            self.Cp = 1000.0
            self.gx = 0.0
            self.gy = -1.0
            self.eta_f = 0.0001
            self.rho_f = 0.3
            self.y_oceanSurface = 0.0


            
            
            
class SingleColorMap(Frozen):
    _Frozen__List = ["type","scale","center","max","log10on"]
    def __init__(self,colormapType="Manual",colormap="Default",scale=1.0,center=0.0,maxValue=1.0,log10on=False):
        self.type       = colormapType # "automatic would go from min to max values"
        self.colormap   = colormap
        self.scale      = scale
        self.center     = center # centered value (scaled)
        self.max        = maxValue # maximum value (scaled) from the center
        self.log10on    = log10on        

        
class ColorMapList(Frozen):
    _Frozen__List = ["type","scale","center","max","log10on"]
    def __init__(self):
        self.Viscosity      = SingleColorMap(log10on=True)
        self.Khi            = SingleColorMap(log10on=True)
        self.Khib           = SingleColorMap(log10on=True)
        self.StrainRate     = SingleColorMap()
        self.Stress         = SingleColorMap()
        self.Velocity       = SingleColorMap()
        self.VelocityDiv    = SingleColorMap()
        self.SIIOvYield     = SingleColorMap()
        self.PeOvYield      = SingleColorMap()
        self.Pressure       = SingleColorMap()
        self.Density        = SingleColorMap()
        self.Temperature    = SingleColorMap()
        self.FluidPressure  = SingleColorMap()
        self.CompactionPressure   = SingleColorMap()
        self.Permeability   = SingleColorMap()
        self.Porosity       = SingleColorMap()
        self.Phase          = SingleColorMap()
        self.VxRes          = SingleColorMap()
        self.VyRes          = SingleColorMap()
        self.PRes           = SingleColorMap()
        self.PfRes          = SingleColorMap()
        self.PcRes          = SingleColorMap()
        self.TRes           = SingleColorMap()

class Visu(Frozen):
    _Frozen__List = ["type","typeParticles","showParticles","shiftFacX","shiftFacY","shiftFacZ","writeImages","transparency","alphaOnValue","showGlyphs","glyphType","glyphMeshType","glyphScale","glyphSamplingRateX","glyphSamplingRateY","width","height","outputFolder","retinaScale","particleMeshRes","particleMeshSize","filter"]
    def __init__(self):
        self.type 	    = "StrainRate" # Default
        self.typeParticles  = "PartPhase" # Default
        self.showParticles  = True
        self.shiftFacX      = 0.00
        self.shiftFacY 	= 0.00
        self.shiftFacZ 	 = -0.05
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
        
        self.colormaps
    
        
        
        
        
        
        
        
        
        
        
        
        

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


    def set_based_on_lithostatic_pressure(self,PhaseRef,BCThermal,Physics,Grid,Length=0):
        if (Length == 0):
          self.length = (Grid.ymax-Grid.ymin)/2.0
        else:
          self.length = Length

        CharStress  = PhaseRef.rho0*abs(Physics.gy)*self.length
        CharVisc    = PhaseRef.eta0

        self.time   = CharVisc/CharStress
        self.mass   = CharStress*self.time*self.time*self.length
        self.temperature = (BCThermal.TB + BCThermal.TT)/2.0
        
    def set_based_on_corner_flow(self,PhaseRef,BCStokes,BCThermal,Physics,Grid,Length=0):
        if (Length == 0):
          self.length = (Grid.ymax-Grid.ymin)/2.0
        else:
          self.length = Length
          
        #CharStress  = PhaseRef.rho0*abs(Physics.gy)*self.length
        CharVisc    = PhaseRef.eta0
        CharStress  = CharVisc*BCStokes.refValue/self.length;
        self.time   = CharVisc/CharStress
        self.mass   = CharStress*self.time*self.time*self.length
        self.temperature = (BCThermal.TB + BCThermal.TT)/2.0
        






class BCStokes(Frozen):
    _Frozen__List = ["backStrainRate","SetupType","refValue"]
    def __init__(self):
        self.backStrainRate = -1.0
        self.SetupType  = "PureShear"
        self.refValue = 1.0;

class BCThermal(Frozen):
    _Frozen__List = ["TT","TB","SetupType","refValue"]
    def __init__(self):
        self.TT = 1.0
        self.TB = 1.0
        self.SetupType  = "TT_TB_LRNoFlux"
        self.refValue = 1.0;


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
        self.wavephase = wavephase
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
