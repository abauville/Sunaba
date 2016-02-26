################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/BC.c \
../src/Physics.c \
../src/eqSystem.c \
../src/main.c \
../src/memory.c \
../src/numbering.c \
../src/particles.c \
../src/utils.c \
../src/visualization.c 

OBJS += \
./src/BC.o \
./src/Physics.o \
./src/eqSystem.o \
./src/main.o \
./src/memory.o \
./src/numbering.o \
./src/particles.o \
./src/utils.o \
./src/visualization.o 

C_DEPS += \
./src/BC.d \
./src/Physics.d \
./src/eqSystem.d \
./src/main.d \
./src/memory.d \
./src/numbering.d \
./src/particles.d \
./src/utils.d \
./src/visualization.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	/usr/local/bin/gcc-5 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


