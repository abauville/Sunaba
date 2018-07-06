
# ============================================
#
# Classes of User Input
#
# Arthur Bauville, June 3, 2016
#
# ============================================
from math import pi, cos, sin, tan, pow

import gc
import json
import copy
import os

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

    def setFromDict(self, dictionary):
        """Constructor"""
        for key in dictionary:
            setattr(self, key, dictionary[key])
            
            
            
from MaterialsDef import Material


          
class Setup(Frozen):
    _Frozen__List = ["Description","Physics","Grid","Numerics","Particles","Char","Visu","BC","IC","Geometry","MatProps","Output"]
    def __init__(self,isDimensional=False):
        self.Description    = ""
        self.Physics        = Physics(isDimensional)
        self.Grid           = Grid()
        self.Numerics       = Numerics()
        self.Particles      = Particles()
        self.Char           = Char()
        self.Visu           = Visu()
        self.BC             = BC() 
        self.IC             = IC() 
        self.Geometry       = {}
        self.MatProps       = {}
        self.Output         = Output()

class Grid(Frozen):
    _Frozen__List = ["xmin","xmax","ymin","ymax","nxC","nyC","fixedBox"]
    def __init__(self):
        self.xmin   = -.5
        self.xmax   =  .5
        self.ymin   = -.5
        self.ymax   =  .5

        self.nxC    = 64
        self.nyC    = 64

        #self.ySeg   = self.xmin
        #self.ySeg   = self.xmin
    
        self.fixedBox = False
class Numerics(Frozen):
    _Frozen__List = ["nTimeSteps", "maxTime", "nLineSearch", "maxNonLinearIter", "minNonLinearIter", "relativeTolerance", "absoluteTolerance","maxCorrection","CFL_fac_Stokes","CFL_fac_Thermal","CFL_fac_Darcy","etaMin","etaMax","phiMin","phiMax","phiCrit","dtMin","dtMax","use_dtMaxwellLimit","stickyAirSwitchingDepth","stickyAirSwitchPhaseTo","stickyAirSwitchPassiveTo","stickyAirTimeSwitchPassive","dtAlphaCorr","dtVep","dtMaxwellFac_EP_ov_E","dtMaxwellFac_VP_ov_E","dtMaxwellFac_VP_ov_EP","dt_stressFac","deltaSigmaMin","dt_plasticFac","yieldComputationType", "invariantComputationType","stressSubGridDiffFac","dtIni"]
    def __init__(self):
        self.nTimeSteps  = 1 #  negative value for infinite
        self.maxTime     = 14*1e9*(3600*24*365) #  in s, by default 14Gyrs
        self.nLineSearch = 1
        self.maxNonLinearIter = 1 # should always be greater than the number of line searches
        self.minNonLinearIter = 1 # should always be greater than the number of line searches
        self.relativeTolerance = 1E-18 # relative tolerance to the one of this time step
        self.absoluteTolerance = 1E-6 # relative tolerance to the first one of the simulation
        self.maxCorrection = 1.0

        self.CFL_fac_Stokes  = 0.75
        self.CFL_fac_Thermal = 10.0
        self.CFL_fac_Darcy   = 0.5
        
        self.use_dtMaxwellLimit = True

        self.etaMin = 1E-6
        self.etaMax = 1E6
        
        self.phiMin     = 1e-5
        self.phiMax     = 0.8
        self.phiCrit    = 1e-4

        self.dtMin = 0.
        self.dtMax = 1e100

        self.stickyAirSwitchingDepth = -1e100; # effectively switched off by default
        self.stickyAirSwitchPhaseTo  = 0
        self.stickyAirSwitchPassiveTo  = 0
        self.stickyAirTimeSwitchPassive  = 1e100

        self.dtAlphaCorr = 1.0 # correction factor for the time step size (0<alpha<=1)

        self.dtVep   = 0.0 # if value is 0 then use dtAdv for the computation, otherwise use the value given

        
        self.dtMaxwellFac_EP_ov_E  = 0.0;   # lowest,       ElastoPlasticVisc   /   G
        self.dtMaxwellFac_VP_ov_E  = 0.5;   # intermediate, ViscoPlasticVisc    /   G
        self.dtMaxwellFac_VP_ov_EP = 0.5;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress

        self.dt_stressFac          = 1.0;
        self.dt_plasticFac         = 0.5    # between 0 and 1; 0 = EP/E limit; 1 = VP/EP limit
        
        self.deltaSigmaMin         = 0.0 * 1e6; # 5 MPa by default
        
        self.yieldComputationType       = 0 # 0: Cell-and-Node, 1:Cell-interp2Node, 2: Markers
        self.invariantComputationType   = 0 # 0: (interp)^2, 1: interp(^2)



        self.stressSubGridDiffFac = 1.0
        self.dtIni = -1.0 # Default Dummy, in the fucntion write_inputfile it is assigned to the value of Char.time if the dummy value has not been overwritten
        

class Particles(Frozen):
    _Frozen__List = ["nPCX","nPCY","noiseFactor","minPartPerCellFactor","maxPartPerCellFactor","passiveGeom","passiveDx","passiveDy"]
    def __init__(self):
        self.nPCX = 4
        self.nPCY = 4
        self.noiseFactor = 0.0

        self.minPartPerCellFactor = 0.2
        self.maxPartPerCellFactor = 3.0

        self.passiveGeom = "Grid"
        self.passiveDx = 1.0
        self.passiveDy = 1.0

class Physics(Frozen):
    _Frozen__List = ["Cp","gx","gy","eta_f","rho_f","y_oceanSurface","Pback"]
    def __init__(self,Dimensional):
        if Dimensional == True:
            self.Cp = 1000.0
            self.gx = 0.0
            self.gy = -9.81
            self.eta_f = 100.0
            self.rho_f = 1000.0
            self.y_oceanSurface = 0.0
            self.Pback = 0.0
        else:
            self.Cp = 1000.0
            self.gx = 0.0
            self.gy = -1.0
            self.eta_f = 0.0001
            self.rho_f = 0.3
            self.y_oceanSurface = 0.0
            self.Pback = 0.0


            
            
            
class SingleColorMap(Frozen):
    _Frozen__List = ["type","colorMap","scale","center","max","log10on","A0number","colorMapRes","alphaAbsThreshold"]
    def __init__(self,colormapType="Manual",colormap="Default",scale=1.0,center=0.0,maxValue=1.0,log10on=False,number=0,alphaAbsThreshold=-1.0):
        self.A0number     = number
        self.type       = colormapType # "automatic would go from min to max values"
        self.colorMapRes= 0
        self.colorMap   = colormap
        self.scale      = scale
        self.center     = center   # centered value (scaled)
        self.max        = maxValue # maximum value (scaled)
        self.log10on    = log10on    
        self.alphaAbsThreshold = alphaAbsThreshold # absolute value of the threshold for transparecny (not affected by log10on), negative values effectively put it off
        if (self.center>self.max):
            raise ValueError( "%r has a max value lower than its center value (%.2e < %.2e)" % (self, self.max,self.center) )

        
class ColorMapList(Frozen):
    _Frozen__List = ["Viscosity","Khi","Khib","StrainRate","Stress","Velocity","VelocityDiv","SIIOvYield","PeOvYield","Pressure","Density","Temperature",
    "FluidPressure","CompactionPressure","Permeability","Porosity","Phase","VxRes","VyRes","PRes","PfRes","PcRes","TRes","Strain","Vorticity","POvPlitho", 
    "EffectiveViscosity", "ShearModulus","ExtraField"]
    def __init__(self):
        self.Viscosity          = SingleColorMap(log10on=True,  number= 1)
        self.StrainRate         = SingleColorMap(log10on=True,  number= 2)
        self.Velocity           = SingleColorMap(               number= 3)
        self.Pressure           = SingleColorMap(               number= 4)
        self.Density            = SingleColorMap(               number= 5)
        self.Temperature        = SingleColorMap(               number= 6)
        self.Stress             = SingleColorMap(               number= 7)
        self.FluidPressure      = SingleColorMap(               number= 8)
        self.Permeability       = SingleColorMap(log10on=True,  number= 9)
        self.Porosity           = SingleColorMap(               number=10)
        self.CompactionPressure = SingleColorMap(               number=11)
        self.Phase              = SingleColorMap(               number=12)
        self.VxRes              = SingleColorMap(log10on=True,  number=13, scale=1e-6, maxValue=2.0)
        self.VyRes              = SingleColorMap(log10on=True,  number=14, scale=1e-6, maxValue=2.0)
        self.PRes               = SingleColorMap(log10on=True,  number=15, scale=1e-6, maxValue=2.0)
        self.PfRes              = SingleColorMap(log10on=True,  number=16, scale=1e-6, maxValue=2.0)
        self.PcRes              = SingleColorMap(log10on=True,  number=17, scale=1e-6, maxValue=2.0)
        self.TRes               = SingleColorMap(log10on=True,  number=18, scale=1e-6, maxValue=2.0)
        self.VelocityDiv        = SingleColorMap(               number=19, scale=1e-6, maxValue=2.0)
        self.SIIOvYield         = SingleColorMap(               number=20,             maxValue=1.5, center=1.0)
        self.PeOvYield          = SingleColorMap(               number=21,             maxValue=1.5, center=1.0)
        self.Khi                = SingleColorMap(log10on=True,  number=22)
        self.Khib               = SingleColorMap(log10on=True,  number=23)
        self.Strain             = SingleColorMap(               number=24)
        self.Vorticity          = SingleColorMap(               number=25)
        self.POvPlitho          = SingleColorMap(               number=26, scale=1.0 , maxValue=2.0, center=1.0)
        self.EffectiveViscosity = SingleColorMap(log10on=True,  number=27)
        self.ShearModulus       = SingleColorMap(log10on=True,  number=28)
        self.ExtraField         = SingleColorMap(               number=29)
    
class Visu(Frozen):
    _Frozen__List = ["type","typeParticles","showParticles","shiftFacX","shiftFacY","shiftFacZ","writeImages","transparency","alphaOnValue","showGlyphs","glyphType","glyphMeshType","glyphScale","glyphSamplingRateX","glyphSamplingRateY","width","height","outputFolder","retinaScale","particleMeshRes","particleMeshSize","filter","colorMap","typeNumber","shaderFolder","renderFrequency", "renderTimeFrequency","closeAtTheEndOfSimulation"]
    def __init__(self):
        self.type           = "StrainRate" # Default
        self.typeNumber         = 0
        self.typeParticles  = "PartPhase" # Default
        self.showParticles  = True
        self.shiftFacX      = 0.00
        self.shiftFacY      = 0.00
        self.shiftFacZ      = -0.05
        self.writeImages    = False
        self.transparency   = False
        self.alphaOnValue   = False
        self.showGlyphs     = False
        self.glyphType      = "StokesVelocity"
        self.glyphMeshType  = "ThinArrow"
        self.glyphScale     = 0.05
        self.glyphSamplingRateX  = 3
        self.glyphSamplingRateY  = 6

        self.width              = 1024
        self.height             = 1024

        self.outputFolder = "../../OutputStokesFD"
        
        self.shaderFolder = "../Shaders/Default" # Relative path from the running folder (of StokesFD)

        self.renderFrequency = 1 # Render every this number of steps
        self.renderTimeFrequency = 0 # Render every this amount of time, if 0: then renderFrequency is used instead

        self.retinaScale = 2

        self.particleMeshRes = 6 # minimum 3, higher number make the particle voronoi diagram look smoother

        self.particleMeshSize = 0.05

        self.filter = "Linear"
        
        self.colorMap = ColorMapList()
    
        self.closeAtTheEndOfSimulation = True
    
    def dictionarize(self):
        self.colorMap = vars(self.colorMap)
        for key in self.colorMap:
            self.colorMap[key] = vars(self.colorMap[key])
         
    def finalize(self):
        self.dictionarize()
        ListOfTypes = ("Blank", "Viscosity", "StrainRate", "Velocity", "Pressure", "Density", "Temperature", "Stress", "FluidPressure", "Permeability", "Porosity", "CompactionPressure", "Phase",
                       "VxRes", "VyRes", "PRes", "PfRes", "PcRes", "TRes", "VelocityDiv","SIIOvYield", "PeOvYield", "Khi", "Khib","Strain","Vorticity","POvPlitho", "EffectiveViscosity", "ShearModulus", "ExtraField")
        self.typeNumber = ListOfTypes.index(self.type)
        #Here goes the automatic computation of colormapRes
    
        
        
        
        
        
        

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
        self.temperature = (BCThermal.TB + BCThermal.TT)/2.0
        CharVisc = PhaseRef.getRefVisc(0,self.temperature,abs(BCStokes.backStrainRate))

        CharStress = 2*CharVisc*abs(BCStokes.backStrainRate)
        self.mass   = CharStress*self.time*self.time*self.length
        
        if (PhaseRef.isRef == False):
            raise ValueError("PhaseRef.isRef == False")

    def set_based_on_lithostatic_pressure(self,PhaseRef,BCStokes,BCThermal,Physics,Grid,Length=0):
        if (Length == 0):
            self.length = (Grid.ymax-Grid.ymin)/2.0
        else:
            self.length = Length

        self.temperature = (BCThermal.TB + BCThermal.TT)/2.0
        CharVisc = PhaseRef.getRefVisc(0,self.temperature,abs(BCStokes.backStrainRate))
          
          
        CharStress  = PhaseRef.rho0*abs(Physics.gy)*self.length


        self.time   = CharVisc/CharStress
        self.mass   = CharStress*self.time*self.time*self.length
        if (PhaseRef.isRef == False):
            raise ValueError("PhaseRef.isRef == False")
        
    def set_based_on_corner_flow(self,PhaseRef,BCStokes,BCThermal,Physics,Grid,Length=0):
        if (Length == 0):
            self.length = (Grid.ymax-Grid.ymin)/2.0
        else:
            self.length = Length
            
        self.temperature = (BCThermal.TB + BCThermal.TT)/2.0
        CharVisc = PhaseRef.getRefVisc(0,self.temperature,abs(BCStokes.backStrainRate))

        
        CharStress  = CharVisc*BCStokes.refValue/self.length;
        self.time   = CharVisc/CharStress
        self.mass   = CharStress*self.time*self.time*self.length
        if (PhaseRef.isRef == False):
            raise ValueError("PhaseRef.isRef == False")
        

            
            
            
class CharExtra(Frozen):
    _Frozen__List = ["visc","stress","strainrate"]
    def __init__(self,Char):
        kg = Char.mass
        m  = Char.length
        s  = Char.time
        #K  = Char.temperature
        self.visc = kg/m/s
        self.stress = kg/m/s/s
        self.strainrate = 1/s     
            

            

class BC(Frozen):
    _Frozen__List = ["Stokes","Thermal"]
    def __init__(self):
        self.Stokes = BCStokes()
        self.Thermal = BCThermal()



class BCStokes(Frozen):
    _Frozen__List = ["backStrainRate","SetupType","refValue","DeltaL","Sandbox_TopSeg00","Sandbox_TopSeg01","Sandbox_NoSlipWall","Corner_SubductionAngle"]
    def __init__(self):
        self.backStrainRate = -1.0
        self.SetupType  = "PureShear"
        self.refValue = 1.0
        self.DeltaL = 1.0
        
        self.Sandbox_TopSeg00 = 0.0
        self.Sandbox_TopSeg01 = 0.0
        self.Sandbox_NoSlipWall = False
        
        self.Corner_SubductionAngle = 30/180*pi
        

class BCThermal(Frozen):
    _Frozen__List = ["TT","TB","SetupType","refValue","DeltaL"]
    def __init__(self):
        self.TT = 1.0
        self.TB = 1.0
        self.SetupType  = "TT_TB_LRNoFlux"
        self.refValue = 1.0
        self.DeltaL = 1.0

class BCDarcy(Frozen):
    _Frozen__List = ["backStrainRate","SetupType","refValue","DeltaL"]
    def __init__(self):
        self.PfT_type = "Dirichlet"
        self.PfT_val  = 0.0
        self.PfB_type = "Dirichlet"
        self.PfB_val  = 0.0
        self.PfL_type = "Neumann"
        self.PfL_val  = 0.0
        self.PfR_type = "Neuman"
        self.PfR_val  = 0.0
       
        
        
        
class IC(Frozen):
    _Frozen__List = ["Darcy","Thermal"]
    def __init__(self):
        self.Thermal = IC_HSC()
        self.Darcy   = IC_Gaussian(background=1e-3)
        
class IC_HSC(Frozen):
    _Frozen__List = ["A0type","Tm","age","noise"]
    def __init__(self,noise = 0.0, Tm = 1300.0+273.0, age = 0.0):
        self.A0type = "HSC"
        self.noise   = noise # in K
        self.Tm      = Tm
        self.age     = age
        
    
class IC_Gaussian(Frozen):
    _Frozen__List = ["A0type","noise","background","Amp","xc","yc","wx","wy"]
    def __init__(self,noise = 0.0, background = 0.0, Amp = 0.0, xc =0.0, yc = 0.0, wx = 0.5, wy = 0.5):
        self.A0type  = "Gaussian";
        self.noise   = noise # in the unit of the grandeur (i.e. K for T, dimensionless for porosity)
        self.background = background
        self.Amp     = Amp
        self.xc      = xc
        self.yc      = yc
        self.wx      = wx
        self.wy      = wy
        
        

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
    def __init__(self,phase,base,amplitude,wavephase,wavelength,definedFor,condition,Min,Max,slope = 0.0):
        self.phase = phase
        self.amplitude = amplitude
        self.wavephase = wavephase
        self.base  = base
        self.wavelength = wavelength
        self.definedFor = definedFor
        self.condition = condition
        self.min = Min
        self.max = Max
        self.slope = slope
        if (slope>pi/2.0 or slope<-pi/2.0):
            raise ValueError("the slope is not in the range -pi/2 -> pi/2. Hint: the slope should be given in radians")


class Geom_Polygon(object):
    # x and y should be arrays of vertices
    # the test for polygon is not ready yet
    def __init__(self,phase,x,y):
        self.phase = phase
        self.x = x
        self.y = y

        
        


class Output(Frozen):
    _Frozen__List = ["folder","Vx","Vy","P","Pf","Pc","eta","porosity","Z","G","khi","sigma_xx","sigma_xy","sigma_xx0","sigma_xy0","sigma_II","strainRate","strain","temperature", "phase", "frequency", "timeFrequency", "saveFirstStep", "particles_pos","particles_posIni","particles_phase","particles_passive","particles_T","particles_stress","particles_phi","sigma_xy_node","particles_strain","particles_timeLastPlastic"]
    def __init__(self, folder= "./Output/", Vx=False, Vy=False, P=False, Pf=False, Pc=False, eta=False, porosity=False, Z=False, G=False, khi=False, sigma_xx=False, sigma_xy=False, sigma_xy_node=False, sigma_xx0=False, sigma_xy0=False, sigma_II=False, strainRate=False, strain=False, temperature=False, phase=False, saveFirstStep=True, frequency=1, timeFrequency=0.0):
        self.folder         = folder
        self.Vx             = Vx
        self.Vy             = Vy
        self.P              = P
        self.Pf             = Pf
        self.Pc             = Pc
        self.eta            = eta
        self.porosity       = porosity
        self.Z              = Z
        self.G              = G
        self.khi            = khi
        self.sigma_xx       = sigma_xx
        self.sigma_xy       = sigma_xy
        self.sigma_xy_node  = sigma_xy_node
        self.sigma_xx0      = sigma_xx0
        self.sigma_xy0      = sigma_xy0
        self.sigma_II       = sigma_II
        self.strainRate     = strainRate
        self.strain         = strain
        self.temperature    = temperature
        self.phase          = phase
        self.particles_pos      = False
        self.particles_posIni   = False
        self.particles_phase    = False
        self.particles_passive  = False
        self.particles_T        = False
        self.particles_stress   = False
        self.particles_phi      = False
        self.particles_strain   = False
        self.particles_timeLastPlastic   = False
        self.frequency      = frequency
        self.saveFirstStep  = saveFirstStep
        self.timeFrequency  = timeFrequency # in seconds. note: if timeFrequency is larger than 0 it will use time frequency instead of frequency
        
        
    def write(self):
        # returns True if at least on of the Output is set to True
        write = False
        for attr, value in self.__dict__.items():
            if type(value) == bool:
                write = write or value
        return write
        
        
        
        
def writeInputFile(Setup,Filename='default'):
    if (Filename=='default'):
        print("A")
        # Look for the root directory
        # The root directory contains the git ignore file
        # So the idea is to look recursively upward for this file
        folder = './'
        while (os.path.isfile(folder + '.gitignore')==False):
            folder = folder + '../'
        Filename = folder + "Setups/input.json"
        

    
    CSetup = copy.deepcopy(Setup) # to avoid transforming the elements of Setup in dicts, which is annoying for parameter checking
    #make dicts
    CSetup.Visu.finalize()
    
    # Some error check, should be moved
    if CSetup.Output.folder[-1]!="/" :
        CSetup.Output.folder = CSetup.Output.folder + "/"
    
    if CSetup.Visu.outputFolder[-1]!="/" :
        CSetup.Visu.outputFolder = CSetup.Visu.outputFolder + "/"
    
    
    print(CSetup.Output.folder)
    
    if CSetup.Output.frequency == 0:
        raise ValueError("CSetup.Output.frequency == 0, should be at least 1")
    
    
    if CSetup.Visu.showGlyphs:
        if (CSetup.Visu.glyphSamplingRateX==0):
            raise ValueError("CSetup.Visu.glyphSamplingRateX == 0")
        if (CSetup.Visu.glyphSamplingRateY==0):
            raise ValueError("CSetup.Visu.glyphSamplingRateY == 0")    
    
    if (CSetup.BC.Stokes.Sandbox_TopSeg00 > CSetup.BC.Stokes.Sandbox_TopSeg01):
        raise ValueError("CSetup.BC.Stokes.Sandbox_TopSeg00 > CSetup.BC.Stokes.Sandbox_TopSeg00Sandbox_TopSeg01 , should be <=")
    
    
    
    
    
    
    for key in CSetup.Geometry:
        try:
            CSetup.MatProps[str(CSetup.Geometry[key].phase)]
        except KeyError:
            raise ValueError("The geometry object %s was CSetup with phase # %i, however this phase is not defined in CSetup.MatProps." % (key,CSetup.Geometry[key].phase))
        CSetup.Geometry[key] = vars(CSetup.Geometry[key])
       
    for key in CSetup.MatProps:
        if CSetup.MatProps[key].vDisl.B == 0:
            raise ValueError("vDisl.B=0 in Phase %s, use the function correctUnitsAndComputeB() to compute it (also, be careful with units)" % CSetup.MatProps[key]);
        if CSetup.MatProps[key].vDiff.B == 0:
            raise ValueError("vDiff.B=0 in Phase %s, use the function correctUnitsAndComputeB() to compute it (also, be careful with units)" % CSetup.MatProps[key]);
        if CSetup.MatProps[key].vPei.B == 0:
            raise ValueError("vPei.B=0 in Phase %s, use the function correctUnitsAndComputeB() to compute it (also, be careful with units)" % CSetup.MatProps[key]);
                
        CSetup.MatProps[key].vDisl = vars(CSetup.MatProps[key].vDisl)
        CSetup.MatProps[key].vDiff = vars(CSetup.MatProps[key].vDiff)
        CSetup.MatProps[key].vPei  = vars(CSetup.MatProps[key].vPei)
        CSetup.MatProps[key] = vars(CSetup.MatProps[key])

            
    CSetup.BC.Stokes   = vars(CSetup.BC.Stokes)
    CSetup.BC.Thermal  = vars(CSetup.BC.Thermal)
    
    CSetup.IC.Thermal  = vars(CSetup.IC.Thermal)
    CSetup.IC.Darcy    = vars(CSetup.IC.Darcy)
    
    
    
    
    
    if CSetup.Numerics.dtIni == -1.0:
        CSetup.Numerics.dtIni = CSetup.Char.time
    
    
    
    
    myJsonFile = dict(Description = CSetup.Description, Grid = vars(CSetup.Grid), Numerics = vars(CSetup.Numerics), Particles = vars(CSetup.Particles), Physics = vars(CSetup.Physics), Visu = vars(CSetup.Visu), MatProps = CSetup.MatProps, Char = vars(CSetup.Char), BC = vars(CSetup.BC), IC = vars(CSetup.IC), Geometry = CSetup.Geometry, Output = vars(CSetup.Output));

    json.dump(myJsonFile, open(Filename, 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)


