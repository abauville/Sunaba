################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/MiscLibraries/jsmn.c 

OBJS += \
./src/MiscLibraries/jsmn.o 

C_DEPS += \
./src/MiscLibraries/jsmn.d 


# Each subdirectory must supply rules for building sources it contributes
src/MiscLibraries/%.o: ../src/MiscLibraries/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	icc -I/usr/local/include -O0 -g -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


