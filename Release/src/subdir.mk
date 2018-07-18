################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/BC.c \
../src/Breakpoint.c \
../src/Char.c \
../src/Darcy.c \
../src/EqSystem.c \
../src/FreeSurface.c \
../src/Grid.c \
../src/IC.c \
../src/Input.c \
../src/Interp.c \
../src/LocalStencil.c \
../src/Main_Init.c \
../src/MatProps.c \
../src/Numbering.c \
../src/Numerics.c \
../src/Output.c \
../src/Particles.c \
../src/Physics.c \
../src/Physics_CellVal.c \
../src/Physics_Eta.c \
../src/Utils.c \
../src/Visu.c \
../src/main.c 

OBJS += \
./src/BC.o \
./src/Breakpoint.o \
./src/Char.o \
./src/Darcy.o \
./src/EqSystem.o \
./src/FreeSurface.o \
./src/Grid.o \
./src/IC.o \
./src/Input.o \
./src/Interp.o \
./src/LocalStencil.o \
./src/Main_Init.o \
./src/MatProps.o \
./src/Numbering.o \
./src/Numerics.o \
./src/Output.o \
./src/Particles.o \
./src/Physics.o \
./src/Physics_CellVal.o \
./src/Physics_Eta.o \
./src/Utils.o \
./src/Visu.o \
./src/main.o 

C_DEPS += \
./src/BC.d \
./src/Breakpoint.d \
./src/Char.d \
./src/Darcy.d \
./src/EqSystem.d \
./src/FreeSurface.d \
./src/Grid.d \
./src/IC.d \
./src/Input.d \
./src/Interp.d \
./src/LocalStencil.d \
./src/Main_Init.d \
./src/MatProps.d \
./src/Numbering.d \
./src/Numerics.d \
./src/Output.d \
./src/Particles.d \
./src/Physics.d \
./src/Physics_CellVal.d \
./src/Physics_Eta.d \
./src/Utils.d \
./src/Visu.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	/usr/local/bin/gcc-7 -I/usr/local/include/ -O3 -Ofast -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


