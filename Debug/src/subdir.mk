################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/BC.c \
../src/char.c \
../src/darcy.c \
../src/eqSystem.c \
../src/grid.c \
../src/main.c \
../src/memory.c \
../src/numbering.c \
../src/particles.c \
../src/physics.c \
../src/utils.c \
../src/visualization.c 

OBJS += \
./src/BC.o \
./src/char.o \
./src/darcy.o \
./src/eqSystem.o \
./src/grid.o \
./src/main.o \
./src/memory.o \
./src/numbering.o \
./src/particles.o \
./src/physics.o \
./src/utils.o \
./src/visualization.o 

C_DEPS += \
./src/BC.d \
./src/char.d \
./src/darcy.d \
./src/eqSystem.d \
./src/grid.d \
./src/main.d \
./src/memory.d \
./src/numbering.d \
./src/particles.d \
./src/physics.d \
./src/utils.d \
./src/visualization.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	/usr/local/bin/gcc-5 -O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


