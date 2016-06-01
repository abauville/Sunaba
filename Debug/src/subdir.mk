################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/BC.c \
../src/Char.c \
../src/Darcy.c \
../src/EqSystem.c \
../src/Grid.c \
../src/Numbering.c \
../src/Numerics.c \
../src/Particles.c \
../src/Physics.c \
../src/Utils.c \
../src/Visualization.c \
../src/main.c 

OBJS += \
./src/BC.o \
./src/Char.o \
./src/Darcy.o \
./src/EqSystem.o \
./src/Grid.o \
./src/Numbering.o \
./src/Numerics.o \
./src/Particles.o \
./src/Physics.o \
./src/Utils.o \
./src/Visualization.o \
./src/main.o 

C_DEPS += \
./src/BC.d \
./src/Char.d \
./src/Darcy.d \
./src/EqSystem.d \
./src/Grid.d \
./src/Numbering.d \
./src/Numerics.d \
./src/Particles.d \
./src/Physics.d \
./src/Utils.d \
./src/Visualization.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	/usr/local/bin/gcc-5 -I/usr/local/include/ -O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


