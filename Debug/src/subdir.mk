################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/main.c \
../src/memory.c \
../src/particles.c \
../src/visualization.c 

OBJS += \
./src/main.o \
./src/memory.o \
./src/particles.o \
./src/visualization.o 

C_DEPS += \
./src/main.d \
./src/memory.d \
./src/particles.d \
./src/visualization.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	/usr/local/bin/gcc-5 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


