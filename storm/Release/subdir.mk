################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Alignment.cpp \
../Config.cpp \
../DoubleKeyHashTable.cpp \
../GeneralQueryDatasetFilter.cpp \
../HashTable.cpp \
../HashTableMethod.cpp \
../OmegaGraphConstructor.cpp \
../OmegaHashTable.cpp \
../OverlapGraphConstructor.cpp \
../PairEndReadMerger.cpp \
../PairedEndRead.cpp \
../PairedEndReadsMerger.cpp \
../QueryDataset.cpp \
../QueryDatasetFilter.cpp \
../QueryRead.cpp \
../SingleKeyHashTable.cpp \
../SubjectDataset.cpp \
../SubjectRead.cpp \
../main.cpp 

OBJS += \
./Alignment.o \
./Config.o \
./DoubleKeyHashTable.o \
./GeneralQueryDatasetFilter.o \
./HashTable.o \
./HashTableMethod.o \
./OmegaGraphConstructor.o \
./OmegaHashTable.o \
./OverlapGraphConstructor.o \
./PairEndReadMerger.o \
./PairedEndRead.o \
./PairedEndReadsMerger.o \
./QueryDataset.o \
./QueryDatasetFilter.o \
./QueryRead.o \
./SingleKeyHashTable.o \
./SubjectDataset.o \
./SubjectRead.o \
./main.o 

CPP_DEPS += \
./Alignment.d \
./Config.d \
./DoubleKeyHashTable.d \
./GeneralQueryDatasetFilter.d \
./HashTable.d \
./HashTableMethod.d \
./OmegaGraphConstructor.d \
./OmegaHashTable.d \
./OverlapGraphConstructor.d \
./PairEndReadMerger.d \
./PairedEndRead.d \
./PairedEndReadsMerger.d \
./QueryDataset.d \
./QueryDatasetFilter.d \
./QueryRead.d \
./SingleKeyHashTable.d \
./SubjectDataset.d \
./SubjectRead.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -fopenmp -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


