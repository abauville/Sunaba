/*
 * MatProps.c
 *
 *  Created on: Dec 12, 2016
 *      Author: abauville
 */



/*
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
*/
