# ============================================
#
# Classes of User Input
#
# Arthur Bauville, June 3, 2016
#
# ============================================



class Grid(Model):
    __slots__("xmin","xmax","ymin","ymax","nxC","nyC")
    def __init__(self):
        self.xmin   = -1
        self.xmax   =  1
        self.ymin   = -1
        self.ymax   =  1

        self.nxC    = 64
        self.nyC    = 64
        

class Numerics(object):
    def __init__(self):
        self.nTimeSteps  = 1, #  negative value for infinite
        self.nLineSearch = 1,
        self.maxNonLinearIter = 1, # should always be greater than the number of line searches
        self.minNonLinearIter = 1, # should always be greater than the number of line searches
        self.relativeTolerance = 3E-5, # relative tolerance to the one of this time step
        self.absoluteTolerance = 3E-5, # relative tolerance to the first one of the simulation
        self.maxCorrection = 1.0

        self.CLF_fac = 0.5
        


class Material(object):
    def __init__(self):
            name = ""
            Material = "Default"
            n = 1.0
            cohesion = 1.0
            frictionAngle = 30.0*180/pi
            rho0 = 1
            eta0 = 1
            
            alpha = 1
            beta = 1
            


class Particles(object):
    def __init__(self):
        self.nPCX = 4
        self.nPCY = 4
        self.noiseFactor = 0.0



class Physics(object):
    def __init__(self):
        self.Cp = 1.0
        self.gx = 0.0
        self.gy = 1.0



class Visu(object):
    def __init__(self):
        self.type 			= StrainRate; # Default
        self.typeParticles	= Phase; # Default
        self.showParticles  = true;
        self.shiftFac[0]    = 0.00;
        self.shiftFac[1] 	= 0.00;
        self.shiftFac[2] 	= +.05;
        self.writeImages 	= true;
        self.transparency 	= true;
        self.alphaOnValue 	= true;
        self.showGlyphs 	= true;
        self.glyphType		= StokesVelocity;
        self.glyphMeshType	= ThinArrow;
        self.glyphScale		= 0.05;
        self.glyphSamplingRateX  = 3;
        self.glyphSamplingRateY  = 6;

        self.outputFolder = "../../OutputStokesFD"

        self.retinaScale = 1
        
        Visu.particleMeshRes = 6; # minimum 3, higher number make the particle voronoi diagram look smoother


class Char(object):
    def __init__(self):
         self.length = 1
         self.mass = 1
         self.time = 1
         self.temp = 1

    def compute(self,Grid,Physics,MatProps):
        pass





