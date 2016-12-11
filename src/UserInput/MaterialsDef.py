#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 16:53:49 2016

@author: abauville
"""

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

class Material(Frozen):
    _Frozen__List = ["name","material","n","cohesion","frictionAngle","rho0","eta0","alpha","beta","k","G","perm0","eta_b",
    "B","isAir","isWater", "isRef","vDisl","vDiff","vPei"]
    def __init__(self,material="Default",name=""):
        self.isRef    = False
        if material == "Default":
            self.name = name
            self.material = "Default"
            self.isAir = False
            self.isWater = False
            
            
            # General
            self.rho0 = 1.0
            
            
            # Temperature
            self.alpha = 0.01
            self.beta = 0.01
            self.k = 1.0
            
            
            # Rheology
            # Plasticity
            self.cohesion = 1E100
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E100
            
            
            # Viscosity
            self.vDisl = DislocationCreep("Default")
            self.vDiff = PeierlsCreep("Default")
            self.vPei  = PeierlsCreep("Default")

            self.eta0 = 1.0
            self.n = 1.0
            
            
            
            # Darcy
            self.perm0  = 0.0001
            self.eta_b  = 1.0
            self.B      = 1E20

            
            
        

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

            
            
            
            
# The definition of the flow laws and the material compilation has been borrowed from LaMEM (Kaus, Popov et al.)    
            
# We assume that the creep law has the form:
# Diffusion:   eII = Ad*Tau   * C_OH^r * d^-p *exp( - (Ed + P*Vd)/(R*T))
# Dislocation: eII = An*Tau^n * C_OH^r        *exp( - (En + P*Vn)/(R*T))
# In addition, we take into account that the creep-laws are typically measured under uniaxial or simple shear,
# whereas we need them in tensorial format (self.tensorCorrection and F2) as defined in T. Gerya book.
#
# The resulting expressions for effective viscosity:
# Diffusion:   inv_eta_diff = 2 * [Bd * exp(-(Ed + P*Vd)/(R*T))]
# Dislocation: inv_eta_disl = 2 * [Bn * exp(-(En + P*Vn)/(R*T))]^(1/n) * eII^(1-1/n)
#
# In LaMEM we include the effect of grain size, H2O and tensor correction in the pre-factor (Bd,Bn) such that:
# Diffusion:   Bd  = (2*F2)^(-1) * Ad [Pa] * d^-p * C_OH^r
# Dislocation: Bn  = (2*F2)^(-n) * An [Pa]        * C_OH^r
#
#   eII     -   strain rate             [1/s]
#   Tau     -   stress                  [Pa]
#   P       -   pressure                [Pa]
#   R       -   gas constant
#   Ad, An  -   prefactor (Bn before taking into account grain size and water fugacity) [Pa^(-n)s^(-1)]
#   Bd, Bn  -   prefactor               [Pa^(-n)s^(-1)]
#   n       -   power-law exponent (n=1 for diffusion creep)
#   Ed, En  -   activation Energy       [J/self.MPa/mol]
#   Vd, Vn  -   activation volume       [m^3/mol]
#   d       -   grain size              [in micro-meters (1e-6 meter)]
#   p       -   exponent of grain size
#   C_OH    -   water fugacity in H/10^6 Si  (see Hirth & Kohlstedt 2003 for a description)
#   r       -   power-law exponent of C_OH term
#   MPa     -   transform units: 0 - units in Pa 1 - units in self.MPa
    
            
            
            
            
class DislocationCreep(Frozen):
    _Frozen__List = ["Bn","n","En","Vn","tensorCorrection","MPa","C_OH_0","r"]
    def __init__(self,material="Default",name=""):
        if flowLaw == "Dry_Olivine-Ranalli_1995":	
            # after Ranalli 1995
            self.Bn                 =   2.5e4
            self.n                  =   3.5
            self.En                 =   532e3
            self.Vn                 =   17e-6
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Wet_Olivine-Ranalli_1995":	
            # after Ranalli 1995
            self.Bn                 =   2.0e3
            self.n                  =   4.0
            self.En                 =   471e3
            self.Vn                 =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  ß=   0
	

        elif flowLaw == "Quartz_Diorite-Hansen_Carter_1982":
            # taken from Carter and Tsenn (1986). Flow properties of continental lithosphere - page 18.
            self.Bn            =   pow(10,-1.5)
            self.n             =   2.4
            self.En            =   212e3
            self.Vn            =   0
            self.tensorCorrection =   "SimpleShear"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Diabase-Caristan_1982":
            # Taken from J. de Bremond d'Ars et al./Tectonophysics (1999). Hydrothermalism and Diapirism in the Archaean: gravitational instability constrains. - page 5
            self.Bn            =   6e-2
            self.n             =   3.05
            self.En            =   276e3
            self.Vn            =   1
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   False
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Tumut_Pond_Serpentinite-Raleigh_Paterson_1965":
            # Taken from J. de Bremond d'Ars et al./Tectonophysics (1999). Hydrothermalism and Diapirism in the Archaean: gravitational instability constrains. - page 5
            self.Bn            =   6.3e-7
            self.n             =   2.8
            self.En            =   66e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Wet_Quarzite-Ranalli_1995":
            # used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
            # Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
            # Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
            # in Ranalli 1995 (page 334 Table 10.3)
            self.Bn            =   3.2e-4
            self.n             =   2.3
            self.En            =   154e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
            

        elif flowLaw == "Quarzite-Ranalli_1995":
            # used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
            # Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
            # Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
            # in Ranalli 1995 (page 334 Table 10.3)
            self.Bn            =   6.7e-6
            self.n             =   2.4
            self.En            =   156e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Mafic_Granulite-Ranalli_1995":
            # used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
            # Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
            # Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
            # in Ranalli 1995 (page 334 Table 10.3)
            self.Bn            =   1.4e4
            self.n             =   4.2
            self.En            =   445e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Maryland_strong_diabase-Mackwell_et_al_1998":
            # Mackwell, Zimmerman & Kohlstedt (1998). High-temperature deformation
            # of dry diabase with application to tectonics on Venus. JGR 103. B1. 975-984. page 980
            self.Bn            =   8
            self.n             =   4.7
            self.En            =   485e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Wet_Quarzite-Ueda_et_al_2008":
            # Parameters used in Ueda et al (PEPI 2008)
            self.Bn            =   pow(10,-3.5)
            self.n             =   2.3
            self.En            =   154e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Diabase-Huismans_et_al_2001":
            # parameters used in Huismans et al 2001
            self.Bn            =   3.2e-20
            self.n             =   3.05
            self.En            =   276e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   False
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Granite-Huismans_et_al_2001":
            # parameters used in Huismans et al 2001
            self.Bn            =   3.16e-26
            self.n             =   3.3
            self.En            =   186.5e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   False
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Dry_Upper_Crust-Schmalholz_Kaus_Burg_2009":
            # granite - Burg And Podladchikov (1999)
            self.Bn            =   3.16e-26
            self.n             =   3.3
            self.En            =   190e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   False
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Weak_Lower_Crust-Schmalholz_Kaus_Burg_2009":
            # diabase - Burg And Podladchikov (1999)
            self.Bn            =   3.2e-20
            self.n             =   3.0
            self.En            =   276e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   False
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Plagioclase_An75-Ranalli_1995":
            self.Bn            =   3.3e-4
            self.n             =   3.2
            self.En            =   238e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	
        elif flowLaw == "Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
            self.Bn            =   1600
            self.n             =   3.5
            self.En            =   520e3
            self.Vn            =   22e-6
            self.tensorCorrection =   "SimpleShear"
            self.MPa              =   True
            self.C_OH_0           =   1000
            self.r                =   1.2
	
        elif flowLaw == "Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003_constant_C_OH":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "wet dislocation (constant C_OH)" parameters
            self.Bn            =   90
            self.n             =   3.5
            self.En            =   480e3
            self.Vn            =   11e-6
            self.tensorCorrection =   "SimpleShear"
            self.MPa              =   True
            self.C_OH_0           =   1000
            self.r                =   1.2
	

        elif flowLaw == "Dry_Olivine_disl_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
            self.Bn            =   1.1e5
            self.n             =   3.5
            self.En            =   530e3
            self.Vn            =   15e-6
            self.tensorCorrection =   "SimpleShear"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	

        elif flowLaw == "Olivine-Burg_Podladchikov_1999":
            # after Burg and Podladchikov 1999
            self.Bn            =   7.1e-14
            self.n             =   3.0
            self.En            =   510e3
            self.Vn            =   0
            self.tensorCorrection =   "SimpleShear"
            self.MPa              =   False
            self.C_OH_0           =   1
            self.r                =   0
	
        elif flowLaw == "Wet_Upper_Mantle-Burg_Schmalholz_2008":
            # used in  SchmalholzKausBurg(2009), Geology (wet olivine)
            self.Bn            =   2e-21
            self.n             =   4.0
            self.En            =   471e3
            self.Vn            =   0
            self.tensorCorrection =   "SimpleShear"
            self.MPa              =   False
            self.C_OH_0           =   1
            self.r                =   0

        elif flowLaw == "Granite-Tirel_et_al_2008":
            # used in  SchmalholzKausBurg(2009), Geology
            self.Bn            =   1.25e-9
            self.n             =   3.2
            self.En            =   123e3
            self.Vn            =   0
            self.tensorCorrection =   "SimpleShear"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
	
        elif flowLaw == "Ara_rocksalt-Urai_et_al.(2008)":
            # Ara rocksalt as published in Urai et al.(2008)
            self.Bn            =   1.82e-9
            self.n             =   5
            self.En            =   32.4e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self. r                =   0

        elif flowLaw == "Polycrystalline_Anhydrite-Mueller_and_Briegel(1978)":
            self.Bn            =   3.16228e1
            self.n             =   2
            self.En            =   152.3e3
            self.Vn            =   0
            self.tensorCorrection =   "UniAxial"
            self.MPa              =   True
            self.C_OH_0           =   1
            self.r                =   0
    
    
        else
		raise ValueError("No such dislocation creep profile: %s! " % flowLaw)
	
 
 
 

            
class DiffusionCreep(Frozen):
    _Frozen__List = ["self.tensorCorrection","self.MPa","self.d0","p","self.C_OH_0","r"]
    def __init__(self,flowLaw="Default"):
        self.flowLaw = flowLaw
        if flowLaw == "Default"
            self.Bd                 =   1.0
            self.Ed                 =   0.0
            self.Vd                 =   0.0
            self.tensorCorrection   =   "None"
            self.MPa                =   false
            self.d0                 =   1.0
            self.p                  =   1.0
            self.C_OH_0             =   1.0
            self.r                  =   0.0

        elif flowLaw == "Dry_Olivine_diff_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            self.Bd                 =   1.5e9
            self.Ed                 =   375e3
            self.Vd                 =   5e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.self.d0            =   10e3
            self.p                  =   3
            self.self.C_OH_0        =   1
            self.r                  =   0
	

        elif flowLaw == "Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003_constant_C_OH":	
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            self.Bd                 =   1.0e6
            self.Ed                 =   335e3
            self.Vd                 =   4e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.d0                 =   10e3
            self.p                  =   3
            self.C_OH_0             =   1000
            self.r                  =   1
	

        elif flowLaw == "Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            self.Bd                 =   2.5e7
            self.Ed                 =   375e3
            self.Vd                 =   10e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.d0                 =   10e3
            self.p                  =   3
            self.C_OH_0             =   1000
            self.r                  =   0.8
	

        else:
		raise ValueError( "No such diffusion creep profile: %s! " % flowLaw)













	

//---------------------------------------------------------------------------
// set Peierls creep profiles from literature
#undef __FUNCT__
#define __FUNCT__ "SetPeirProfile"
PetscErrorCode SetPeirProfile(Material_t *m, char name[])

# We assume that the creep law has the form:
# Peierls:   eII = Bp * exp( - (EP + P*VP)/(R*T)*(1-gamma)^q) * (Tau/gamma/taup)^s
#            s   = (Ep+p*Vp)/(R*T)*(1-gamma)^(q-1)*q*gamma
#
# where:
# Bp         - pre-exponential constant for the Peierls mechanism [1/s]
# Ep         - activation energy [J/mol K]
# Vp         - activation volume [m3/mol ]
# taup       - Peierl stress [Pa]
# gamma      - adjustable constant [-]
# q          - stress dependence for Peierls creep [-]
# s          - Peierls creep exponent (typical values between 7-11) [-]

	PetscFunctionBegin

	flowLaw == "Olivine_Peierls-Kameyama_1999":
	
	# used in Kameyama et al 1999 (EPSL), vol 168., pp. 159-172
	# original source: Guyot and Dorn (1967) and Poirier (1985)
		self.Bp            = 5.7e11
		self.Ep            = 5.4e5
		self.Vp            = 0.0
		self.taup          = 8.5e9
		self.gamma         = 0.1
		self.q             = 2
	

	else
	
		SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "No such Peierls creep profile: %s! ",name)
	

	PetscFunctionReturn(0)

//---------------------------------------------------------------------------
// set tensor and units correction for rheological profiles
#undef __FUNCT__
#define __FUNCT__ "SetProfileCorrection"
PetscErrorCode SetProfileCorrection(PetscScalar *B, PetscScalar n, self.tensorCorrection self.tensorCorrection, PetscInt self.self.MPa)

	PetscScalar F2, Bi
# Lab. experiments are typically done under simple shear or uni-axial
# compression, which requires a correction in order to use them in tensorial format.
# An explanation is given in the textbook of Taras Gerya, chapter 6, p. 71-78.

	PetscFunctionBegin

	Bi = *B

# Tensor correction
# In LaMEM this is added to the pre-factor and not to the effective viscosity as in T. Gerya
	if      (self.tensorCorrection == "UniAxial")    F2 = pow(0.5,(n-1)/n) / pow(3,(n+1)/(2*n)) //  F2 = 1/2^((n-1)/n)/3^((n+1)/2/n)
	else if (self.tensorCorrection == "SimpleShear") F2 = pow(0.5,(2*n-1)/n)                    //  F2 = 1/2^((2*n-1)/n)
	else if (self.tensorCorrection == "None")        F2 = 0.5
	else
	
		 SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Unknown tensor correction in creep mechanism profile!")
	

# Units correction from [self.MPa^(-n)s^(-1)] to [Pa^(-n)s^(-1)] if required
	if (self.MPa) Bi = pow(2*F2,-n) * pow(1e6*pow(Bi,-1/n),-n)
	else     Bi = pow(2*F2,-n) * Bi

	(*B) = Bi

	PetscFunctionReturn(0)

class PeierlsCreep(Frozen):
    _Frozen__List = [""]
    def __init__(self,material="Default",name=""):
        print("")            
            