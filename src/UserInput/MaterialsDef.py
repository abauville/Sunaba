#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 16:53:49 2016

@author: abauville
"""

from math import pi, pow, exp

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
    _Frozen__List = ["name","material","cohesion","frictionAngle","rho0","alpha","beta","k","G","perm0",
    "isAir","isWater", "isRef","vDisl","vDiff","vPei"]
    def __init__(self,material="Default",name=""):
        self.isRef    = False
        self.name = name
        self.material = material
        self.isAir = False
        self.isWater = False
        if material == "Default":
            # Density
            self.rho0 = 1.0
            
            #  Heat
            self.alpha = 0.0
            self.beta = 0.0
            self.k = 1.0
            
            # Rheology
            # Plasticity
            self.cohesion = 1E100
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E100
            
            # Viscosity
            self.vDisl = DislocationCreep   ()
            self.vDiff = DiffusionCreep     ("Off")
            self.vPei  = PeierlsCreep       ("Off")
            
            # Darcy
            self.perm0  = 0.0001

            
            
        

        elif material == "StickyAir":
            self.isAir = True;
            # Density
            self.rho0 = 0.1
            self.alpha = 0.0
            self.beta  = 0.0
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 50E6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E20
            
            # Viscosity
            #self.vDisl = DislocationCreep   (eta0=1E17)
            #self.vDiff = DiffusionCreep     ("Off")
            self.vDisl = DislocationCreep   ("Off")
            self.vDiff = DiffusionCreep     (eta0=1e17)
            self.vPei  = PeierlsCreep       ("Off")

            
            # Darcy
            self.perm0  = 1E-5
            

        elif material == "StickyWater":
            self.isWater = True;
            # Density
            self.rho0 = 1.0
            self.alpha = 0.0
            self.beta  = 0.0
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 50E6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E20
            
            # Viscosity
            self.vDisl = DislocationCreep   ("Off")
            self.vDiff = DiffusionCreep     (eta0=1E17)
            self.vPei  = PeierlsCreep       ("Off")

            # Darcy
            self.perm0  = 1E-5

            
        elif material == "Sediments":
            # Density
            self.rho0 = 2500
            self.alpha = 1E-5
            self.beta  = 1E-11
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 10e6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E11
            
            # Viscosity
            self.vDisl = DislocationCreep   ("Off")
            self.vDiff = DiffusionCreep     (eta0=1E21)
            self.vPei  = PeierlsCreep       ("Off")

            # Darcy
            self.perm0  = 5E-9
            
        
        elif material == "Wet_Quartzite":
            # Density
            self.rho0 = 2800
            self.alpha = 1E-5
            self.beta  = 1E-11
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 10e6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E11
            
            # Viscosity
            self.vDisl = DislocationCreep   ("Wet_Quarzite-Ueda_et_al_2008")
            self.vDiff = DiffusionCreep     ("Off")
            self.vPei  = PeierlsCreep       ("Off")

            # Darcy
            self.perm0  = 5E-9
            
        elif material == "Quartzite":
            # Density
            self.rho0 = 2800
            self.alpha = 1E-5
            self.beta  = 1E-11
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 10e6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E11
            
            # Viscosity
            self.vDisl = DislocationCreep   ("Quarzite-Ranalli_1995")
            self.vDiff = DiffusionCreep     ("Off")
            self.vPei  = PeierlsCreep       ("Off")

            # Darcy
            self.perm0  = 5E-9
            
            
            

        elif material == "Dry_Olivine":
            # Density
            self.rho0 = 3300
            self.alpha = 1E-5
            self.beta  = 1E-11
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 10e6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E11
            
            # Viscosity
            self.vDisl = DislocationCreep   ("Dry_Olivine-Ranalli_1995")
            self.vDiff = DiffusionCreep     ("Dry_Olivine_diff_creep-Hirth_Kohlstedt_2003")
            self.vPei  = PeierlsCreep       ("Olivine_Peierls-Kameyama_1999")

            # Darcy
            self.perm0  = 5E-9


        elif material == "Wet_Olivine":
            # Density
            self.rho0 = 3300
            self.alpha = 1E-5
            self.beta  = 1E-11
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 10e6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E11
            
            # Viscosity
            self.vDisl = DislocationCreep   ("Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003_constant_C_OH")
            self.vDiff = DiffusionCreep     ("Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003_constant_C_OH")
            self.vPei  = PeierlsCreep       ("Olivine_Peierls-Kameyama_1999")

            # Darcy
            self.perm0  = 5E-9
            
            
        elif material == "Diabase":
            # Density
            self.rho0 = 2900
            self.alpha = 1E-5
            self.beta  = 1E-11
            
            # Heat
            self.k = 3.0
            
            # Rheology
            # Plasticity
            self.cohesion = 10e6
            self.frictionAngle = 30.0/180*pi
            
            # Elasticity
            self.G = 1E11
            
            # Viscosity
            self.vDisl = DislocationCreep   ("Maryland_strong_diabase-Mackwell_et_al_1998")
            self.vDiff = DiffusionCreep     ("Off")
            self.vPei  = PeierlsCreep       ("Off")

            # Darcy
            self.perm0  = 5E-9
        
        else :
            raise ValueError("No such dislocation material: %s " % (material) )
            
            
            
        
            
    def dictionarize(self):
        self.vDisl = vars(self.vDisl)
        self.vDiff = vars(self.vDiff)
        self.vPei  = vars(self.vPei)
        
        
    def getRefVisc(self,P,T,Eii):
        R = 8.3144598
        invEtaDiff = 0.0
        invEtaDisl = 0.0
        # note: P*V is omitted in the dislocation creep, and etaPei is also omitted
        if self.vDisl.isActive:
            invEtaDisl = (2.0*pow(self.vDisl.B,1.0/self.vDisl.n)*pow(abs(Eii),(-1.0/self.vDisl.n+1.0))) * exp( - (self.vDisl.E+P*self.vDisl.V) / (self.vDisl.n*R*T))
        if self.vDiff.isActive:
            invEtaDiff = (2.0*self.vDiff.B * exp( - (self.vDiff.E+P*self.vDiff.V)  / (R*T)))
        
        print(  (self.vDisl.n*R*T))
        print((self.vDisl.E+P*self.vDisl.V) )
        print (exp( - (self.vDisl.E+P*self.vDisl.V) / (self.vDisl.n*R*T) ))
        print("Diff:" + str(invEtaDiff))
        print("Disl:" + str(invEtaDisl))
        return 1.0/(invEtaDisl+invEtaDiff)
            
# The definition of the flow laws and the material compilation has been borrowed from LaMEM (Kaus, Popov et al.)    
            
# We assume that the creep law has the form:
# Diffusion:   eII = Ad*Tau   * C_OH^r * d^-p *exp( - (Ed + P*Vd)/(R*T))
# Dislocation: eII = An*Tau^n * C_OH^r        *exp( - (En + P*Vn)/(R*T))
# Peierls:     eII = Bp * exp( - (EP + P*VP)/(R*T)*(1-gamma)^q) * (Tau/gamma/taup)^s
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
#   Ed, En  -   activation Energy       [J/mol]
#   Vd, Vn  -   activation volume       [m^3/mol]
#   d       -   grain size              [in micro-meters (1e-6 meter)]
#   p       -   exponent of grain size
#   C_OH    -   water fugacity in H/10^6 Si  (see Hirth & Kohlstedt 2003 for a description)
#   r       -   power-law exponent of C_OH term
#   MPa     -   transform units: 0 - units in Pa 1 - units in self.MPa
    
# For Peierls:
#  s   = (Ep+p*Vp)/(R*T)*(1-gamma)^(q-1)*q*gamma
# Bp         - pre-exponential constant for the Peierls mechanism [1/s]
# Ep         - activation energy [J/mol]
# Vp         - activation volume [m^3/mol]
# taup       - Peierl stress [Pa]
# gamma      - adjustable constant [-]
# q          - stress dependence for Peierls creep [-]
# s          - Peierls creep exponent (typical values between 7-11) [-]            
            
            
            
class DislocationCreep(Frozen):
    _Frozen__List = ["flowLaw","B","A","n","E","V","tensorCorrection","MPa","C_OH_0","r","isActive"]
    def __init__(self,flowLaw="Default",eta0=1.0,n=1.0):
        self.flowLaw = flowLaw
        self.isActive = True
        self.B = 0;
        if flowLaw == "Off":
            self.isActive = False
            flowLaw = "Default"
            
            
        if flowLaw == "Default":	
            self.A                  =   0.5/eta0
            self.n                  =   n
            self.E                  =   0.0
            self.V                  =   0.0
            self.tensorCorrection   =   "None"
            self.MPa                =   False
            self.C_OH_0             =   1.0
            self.r                  =   0.0
        
            
        elif flowLaw == "Dry_Olivine-Ranalli_1995":	
            # after Ranalli 1995
            self.A                  =   2.5e4
            self.n                  =   3.5
            self.E                  =   532e3
            self.V                  =   17e-6
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Wet_Olivine-Ranalli_1995":	
            # after Ranalli 1995
            self.A                  =   2.0e3
            self.n                  =   4.0
            self.E                  =   471e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Quartz_Diorite-Hansen_Carter_1982":
            # taken from Carter and Tsenn (1986). Flow properties of continental lithosphere - page 18.
            self.A                  =   pow(10,-1.5)
            self.n                  =   2.4
            self.E                  =   212e3
            self.V                  =   0
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Diabase-Caristan_1982":
            # Taken from J. de Bremond d'Ars et al./Tectonophysics (1999). Hydrothermalism and Diapirism in the Archaean: gravitational instability constrains. - page 5
            self.A                  =   6e-2
            self.n                  =   3.05
            self.E                  =   276e3
            self.V                  =   1
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   False
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Tumut_Pond_Serpentinite-Raleigh_Paterson_1965":
            # Taken from J. de Bremond d'Ars et al./Tectonophysics (1999). Hydrothermalism and Diapirism in the Archaean: gravitational instability constrains. - page 5
            self.A                  =   6.3e-7
            self.n                  =   2.8
            self.E                  =   66e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Wet_Quarzite-Ranalli_1995":
            # used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
            # Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
            # Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
            # in Ranalli 1995 (page 334 Table 10.3)
            self.A                  =   3.2e-4
            self.n                  =   2.3
            self.E                  =   154e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
            

        elif flowLaw == "Quarzite-Ranalli_1995":
            # used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
            # Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
            # Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
            # in Ranalli 1995 (page 334 Table 10.3)
            self.A                  =   6.7e-6
            self.n                  =   2.4
            self.E                  =   156e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Mafic_Granulite-Ranalli_1995":
            # used in LI, Z. H., GERYA, T. V. and BURG, J.-P. (2010),
            # Influence of tectonic overpressure on PT paths of HPUHP rocks in continental collision zones: thermomechanical modelling.
            # Journal of Metamorphic Geology, 28: 227247. doi: 10.1111/j.1525-1314.2009.00864.x Table 2
            # in Ranalli 1995 (page 334 Table 10.3)
            self.A                  =   1.4e4
            self.n                  =   4.2
            self.E                  =   445e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Maryland_strong_diabase-Mackwell_et_al_1998":
            # Mackwell, Zimmerman & Kohlstedt (1998). High-temperature deformation
            # of dry diabase with application to tectonics on Venus. JGR 103. B1. 975-984. page 980
            self.A                  =   8
            self.n                  =   4.7
            self.E                  =   485e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Wet_Quarzite-Ueda_et_al_2008":
            # Parameters used in Ueda et al (PEPI 2008)
            self.A                  =   pow(10,-3.5)
            self.n                  =   2.3
            self.E                  =   154e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Diabase-Huismans_et_al_2001":
            # parameters used in Huismans et al 2001
            self.A                  =   3.2e-20
            self.n                  =   3.05
            self.E                  =   276e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   False
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Granite-Huismans_et_al_2001":
            # parameters used in Huismans et al 2001
            self.A                  =   3.16e-26
            self.n                  =   3.3
            self.E                  =   186.5e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   False
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Dry_Upper_Crust-Schmalholz_Kaus_Burg_2009":
            # granite - Burg And Podladchikov (1999)
            self.A                  =   3.16e-26
            self.n                  =   3.3
            self.E                  =   190e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   False
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Weak_Lower_Crust-Schmalholz_Kaus_Burg_2009":
            # diabase - Burg And Podladchikov (1999)
            self.A                  =   3.2e-20
            self.n                  =   3.0
            self.E                  =   276e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   False
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Plagioclase_An75-Ranalli_1995":
            self.A                  =   3.3e-4
            self.n                  =   3.2
            self.E                  =   238e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	
        elif flowLaw == "Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
            self.A                  =   1600
            self.n                  =   3.5
            self.E                  =   520e3
            self.V                  =   22e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.C_OH_0             =   1000
            self.r                  =   1.2
	
        elif flowLaw == "Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003_constant_C_OH":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "wet dislocation (constant C_OH)" parameters
            self.A                  =   90
            self.n                  =   3.5
            self.E                  =   480e3
            self.V                  =   11e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.C_OH_0             =   1000
            self.r                  =   1.2
	

        elif flowLaw == "Dry_Olivine_disl_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
            self.A                  =   1.1e5
            self.n                  =   3.5
            self.E                  =   530e3
            self.V                  =   15e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	

        elif flowLaw == "Olivine-Burg_Podladchikov_1999":
            # after Burg and Podladchikov 1999
            self.A                  =   7.1e-14
            self.n                  =   3.0
            self.E                  =   510e3
            self.V                  =   0
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   False
            self.C_OH_0             =   1
            self.r                  =   0
	
        elif flowLaw == "Wet_Upper_Mantle-Burg_Schmalholz_2008":
            # used in  SchmalholzKausBurg(2009), Geology (wet olivine)
            self.A                  =   2e-21
            self.n                  =   4.0
            self.E                  =   471e3
            self.V                  =   0
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   False
            self.C_OH_0             =   1
            self.r                  =   0

        elif flowLaw == "Granite-Tirel_et_al_2008":
            # used in  SchmalholzKausBurg(2009), Geology
            self.A                  =   1.25e-9
            self.n                  =   3.2
            self.E                  =   123e3
            self.V                  =   0
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
	
        elif flowLaw == "Ara_rocksalt-Urai_et_al.(2008)":
            # Ara rocksalt as published in Urai et al.(2008)
            self.A                  =   1.82e-9
            self.n                  =   5
            self.E                  =   32.4e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self. r                 =   0

        elif flowLaw == "Polycrystalline_Anhydrite-Mueller_and_Briegel(1978)":
            self.A                  =   3.16228e1
            self.n                  =   2
            self.E                  =   152.3e3
            self.V                  =   0
            self.tensorCorrection   =   "UniAxial"
            self.MPa                =   True
            self.C_OH_0             =   1
            self.r                  =   0
    
    
        
        
            
        else:
            raise ValueError("No such dislocation creep profile: %s! " % flowLaw)
 
        
            
        self.correctUnitsAndComputeB()
        
        
            
    def correctUnitsAndComputeB(self):
        n = self.n
        if self.MPa == True:
            self.A *= 1e6**(-n)
            
        
        self.B = self.A * self.C_OH_0**(self.r);
        
        if      (self.tensorCorrection == "UniAxial"):
            F2 = pow(0.5,(n-1)/n) / pow(3,(n+1)/(2*n)) #  F2 = 1/2^((n-1)/n)/3^((n+1)/2/n)
        elif (self.tensorCorrection == "SimpleShear"):
            F2 = pow(0.5,(2*n-1)/n)                    #  F2 = 1/2^((2*n-1)/n)
        elif (self.tensorCorrection == "None"):
            F2 = 0.5
        else:
            raise ValueError("Unknown tensor correction in vDisl!")
            
        self.B *= (2.0*F2)**(-n)
            

            
class DiffusionCreep(Frozen):
    _Frozen__List = ["flowLaw","A","B","E","V","tensorCorrection","MPa","d0","p","C_OH_0","r","isActive"]
    def __init__(self,flowLaw="Default",eta0=1.0):
        self.isActive = True
        self.flowLaw = flowLaw
        self.B = 0
        if flowLaw == "Off":
            self.isActive = False
            flowLaw = "Default"
            
            
            
        if flowLaw == "Default":
            self.A                 =   0.5/eta0
            self.E                 =   0.0
            self.V                 =   0.0
            self.tensorCorrection   =  "None"
            self.MPa                =   False
            self.d0                 =   1.0
            self.p                  =   0.0
            self.C_OH_0             =   1.0
            self.r                  =   0.0

        elif flowLaw == "Dry_Olivine_diff_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            self.A                 =   1.5e9
            self.E                 =   375e3
            self.V                 =   5e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.d0            =   10e3
            self.p                  =   3
            self.C_OH_0        =   1
            self.r                  =   0
	

        elif flowLaw == "Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003_constant_C_OH":	
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            self.A                 =   1.0e6
            self.E                 =   335e3
            self.V                 =   4e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.d0                 =   10e3
            self.p                  =   3
            self.C_OH_0             =   1000
            self.r                  =   1
	

        elif flowLaw == "Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003":
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            self.A                 =   2.5e7
            self.E                 =   375e3
            self.V                 =   10e-6
            self.tensorCorrection   =   "SimpleShear"
            self.MPa                =   True
            self.d0                 =   10e3
            self.p                  =   3
            self.C_OH_0             =   1000
            self.r                  =   0.8


            
        else:
            raise ValueError( "No such diffusion creep profile: %s!" % flowLaw)

        
        self.correctUnitsAndComputeB()
        
        
            
    def correctUnitsAndComputeB(self):
        n = 1.0
        if self.MPa == True:
            self.A *= 1e6**(-n) # [MPa^-n s^-1] to [Pa^-n s^-1]
            
        #self.d0 /= 1e6; # [1e-6m] to [m]
        

            
        self.B = self.A * self.C_OH_0**(self.r) * self.d0**(-self.p)
        
        if      (self.tensorCorrection == "UniAxial"):
            F2 = pow(0.5,(n-1)/n) / pow(3,(n+1)/(2*n)) #  F2 = 1/2^((n-1)/n)/3^((n+1)/2/n)
        elif (self.tensorCorrection == "SimpleShear"):
            F2 = pow(0.5,(2*n-1)/n)                    #  F2 = 1/2^((2*n-1)/n)
        elif (self.tensorCorrection == "None"):
            F2 = 0.5
        else:
            raise ValueError("Unknown tensor correction in vDisl!")
            
        self.B *= (2.0*F2)**(-n)




class PeierlsCreep(Frozen):
    _Frozen__List = ["flowLaw","B","E","V","tau","gamma","q","isActive"]
    def __init__(self,flowLaw="Default"):
        self.flowLaw = flowLaw     
        self.isActive = True
        if flowLaw == "Off":
            self.isActive = False
            flowLaw = "Default"    
        
        if flowLaw == "Default":
            self.B            = 1.0
            self.E            = 0.0
            self.V            = 0.0
            self.tau          = 1
            self.gamma         = 1
            self.q             = 1
	
        elif flowLaw == "Olivine_Peierls-Kameyama_1999":
            # used in Kameyama et al 1999 (EPSL), vol 168., pp. 159-172
            # original source: Guyot and Dorn (1967) and Poirier (1985)
            self.B           = 5.7e11
            self.E            = 5.4e5
            self.V            = 0.0
            self.tau          = 8.5e9
            self.gamma         = 0.1
            self.q             = 2

        else :
            raise ValueError(  "No such Peierls creep profile: %s!" % flowLaw);



