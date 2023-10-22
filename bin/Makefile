# Modified by P. Loskot in January 2022
#
# 1. all *.mk files collected into this single makefile
# 2. warnings removed by modifying the code or by -Wno-xx flags
# 3. incremental rebuilds now enabled

# COMMANDS
GG = g++ -O3 -Wall -std=c++11 -fmessage-length=0 -MMD -MP
CC = gcc -O3 -Wall -std=c11 -fmessage-length=0 -MMD -MP
WARNS = -Wno-deprecated-declarations
RM = rm -rf

# VARIABLES
OBJS := 
C_SRCS := 
CC_SRCS := 
C_DEPS := 
CC_DEPS := 

###############################################
# src/nauty24
###############################################

C_SRCS += \
../src/nauty24/nausparse.c \
../src/nauty24/nautil.c \
../src/nauty24/nauty.c 

OBJS += \
./src/nauty24/nausparse.o \
./src/nauty24/nautil.o \
./src/nauty24/nauty.o 

C_DEPS += \
./src/nauty24/nausparse.d \
./src/nauty24/nautil.d \
./src/nauty24/nauty.d 

src/nauty24/%.o: ../src/nauty24/%.c
	@echo 'Building file: $<'
	$(CC) $(WARNS) -Wno-unused-variable -Wno-unused-function \
		-c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFutil
###############################################

CC_SRCS += \
../src/NFutil/MTrand/mtrand.cpp \
../src/NFutil/MTrand/mttest.cpp 

OBJS += \
./src/NFutil/MTrand/mtrand.o \
./src/NFutil/MTrand/mttest.o 

CC_DEPS += \
./src/NFutil/MTrand/mtrand.d \
./src/NFutil/MTrand/mttest.d 

src/NFutil/MTrand/%.o: ../src/NFutil/MTrand/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFutil/conversion.cpp \
../src/NFutil/random.cpp \
../src/NFutil/stringOperations.cpp 

OBJS += \
./src/NFutil/conversion.o \
./src/NFutil/random.o \
./src/NFutil/stringOperations.o 

CC_DEPS += \
./src/NFutil/conversion.d \
./src/NFutil/random.d \
./src/NFutil/stringOperations.d 

src/NFutil/%.o: ../src/NFutil/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFtest
###############################################

CC_SRCS += \
../src/NFtest/transcription/transcription.cpp 

OBJS += \
./src/NFtest/transcription/transcription.o 

CC_DEPS += \
./src/NFtest/transcription/transcription.d 

src/NFtest/transcription/%.o: ../src/NFtest/transcription/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFtest/tlbr/tlbr.cpp 

OBJS += \
./src/NFtest/tlbr/tlbr.o 

CC_DEPS += \
./src/NFtest/tlbr/tlbr.d 

src/NFtest/tlbr/%.o: ../src/NFtest/tlbr/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFtest/simple_system/simple_system.cpp 

OBJS += \
./src/NFtest/simple_system/simple_system.o 

CC_DEPS += \
./src/NFtest/simple_system/simple_system.d 

src/NFtest/simple_system/%.o: ../src/NFtest/simple_system/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFtest/agentcell/cell/cell.cpp \
../src/NFtest/agentcell/cell/environment.cpp \
../src/NFtest/agentcell/cell/util.cpp 

OBJS += \
./src/NFtest/agentcell/cell/cell.o \
./src/NFtest/agentcell/cell/environment.o \
./src/NFtest/agentcell/cell/util.o 

CC_DEPS += \
./src/NFtest/agentcell/cell/cell.d \
./src/NFtest/agentcell/cell/environment.d \
./src/NFtest/agentcell/cell/util.d 

src/NFtest/agentcell/cell/%.o: ../src/NFtest/agentcell/cell/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFtest/agentcell/agentcell.cpp 

OBJS += \
./src/NFtest/agentcell/agentcell.o 

CC_DEPS += \
./src/NFtest/agentcell/agentcell.d 

src/NFtest/agentcell/%.o: ../src/NFtest/agentcell/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFscheduler 
###############################################

CC_SRCS += \
../src/NFscheduler/NFstream.cpp \
../src/NFscheduler/Scheduler.cpp 

OBJS += \
./src/NFscheduler/NFstream.o \
./src/NFscheduler/Scheduler.o 

CC_DEPS += \
./src/NFscheduler/NFstream.d \
./src/NFscheduler/Scheduler.d 

src/NFscheduler/%.o: ../src/NFscheduler/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFreactions
###############################################

CC_SRCS += \
../src/NFreactions/transformations/moleculeCreator.cpp \
../src/NFreactions/transformations/speciesCreator.cpp \
../src/NFreactions/transformations/transformation.cpp \
../src/NFreactions/transformations/transformationSet.cpp 

OBJS += \
./src/NFreactions/transformations/moleculeCreator.o \
./src/NFreactions/transformations/speciesCreator.o \
./src/NFreactions/transformations/transformation.o \
./src/NFreactions/transformations/transformationSet.o 

CC_DEPS += \
./src/NFreactions/transformations/moleculeCreator.d \
./src/NFreactions/transformations/speciesCreator.d \
./src/NFreactions/transformations/transformation.d \
./src/NFreactions/transformations/transformationSet.d 

src/NFreactions/transformations/%.o: ../src/NFreactions/transformations/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFreactions/reactions/DORreaction.cpp \
../src/NFreactions/reactions/reaction.cpp 

OBJS += \
./src/NFreactions/reactions/DORreaction.o \
./src/NFreactions/reactions/reaction.o 

CC_DEPS += \
./src/NFreactions/reactions/DORreaction.d \
./src/NFreactions/reactions/reaction.d 

src/NFreactions/reactions/%.o: ../src/NFreactions/reactions/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFreactions/reactantLists/reactantList.cpp \
../src/NFreactions/reactantLists/reactantTree.cpp 

OBJS += \
./src/NFreactions/reactantLists/reactantList.o \
./src/NFreactions/reactantLists/reactantTree.o 

CC_DEPS += \
./src/NFreactions/reactantLists/reactantList.d \
./src/NFreactions/reactantLists/reactantTree.d 

src/NFreactions/reactantLists/%.o: ../src/NFreactions/reactantLists/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFreactions/mappings/mapping.cpp \
../src/NFreactions/mappings/mappingGenerator.cpp \
../src/NFreactions/mappings/mappingSet.cpp 

OBJS += \
./src/NFreactions/mappings/mapping.o \
./src/NFreactions/mappings/mappingGenerator.o \
./src/NFreactions/mappings/mappingSet.o 

CC_DEPS += \
./src/NFreactions/mappings/mapping.d \
./src/NFreactions/mappings/mappingGenerator.d \
./src/NFreactions/mappings/mappingSet.d 

src/NFreactions/mappings/%.o: ../src/NFreactions/mappings/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFreactions/NFreactions.cpp 

OBJS += \
./src/NFreactions/NFreactions.o 

CC_DEPS += \
./src/NFreactions/NFreactions.d 

src/NFreactions/%.o: ../src/NFreactions/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFoutput
###############################################

CC_SRCS += \
../src/NFoutput/NFoutput.cpp 

OBJS += \
./src/NFoutput/NFoutput.o 

CC_DEPS += \
./src/NFoutput/NFoutput.d 

src/NFoutput/%.o: ../src/NFoutput/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFinput
###############################################

CC_SRCS += \
../src/NFinput/TinyXML/tinystr.cpp \
../src/NFinput/TinyXML/tinyxml.cpp \
../src/NFinput/TinyXML/tinyxmlerror.cpp \
../src/NFinput/TinyXML/tinyxmlparser.cpp 

OBJS += \
./src/NFinput/TinyXML/tinystr.o \
./src/NFinput/TinyXML/tinyxml.o \
./src/NFinput/TinyXML/tinyxmlerror.o \
./src/NFinput/TinyXML/tinyxmlparser.o 

CC_DEPS += \
./src/NFinput/TinyXML/tinystr.d \
./src/NFinput/TinyXML/tinyxml.d \
./src/NFinput/TinyXML/tinyxmlerror.d \
./src/NFinput/TinyXML/tinyxmlparser.d 

src/NFinput/TinyXML/%.o: ../src/NFinput/TinyXML/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFinput/NFinput.cpp \
../src/NFinput/commandLineParser.cpp \
../src/NFinput/parseFuncXML.cpp \
../src/NFinput/parseSymRxns.cpp \
../src/NFinput/rnfRunner.cpp \
../src/NFinput/walk.cpp 

OBJS += \
./src/NFinput/NFinput.o \
./src/NFinput/commandLineParser.o \
./src/NFinput/parseFuncXML.o \
./src/NFinput/parseSymRxns.o \
./src/NFinput/rnfRunner.o \
./src/NFinput/walk.o 

CC_DEPS += \
./src/NFinput/NFinput.d \
./src/NFinput/commandLineParser.d \
./src/NFinput/parseFuncXML.d \
./src/NFinput/parseSymRxns.d \
./src/NFinput/rnfRunner.d \
./src/NFinput/walk.d 

src/NFinput/%.o: ../src/NFinput/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFfunction
###############################################

CC_SRCS += \
../src/NFfunction/muParser/muParser.cpp \
../src/NFfunction/muParser/muParserBase.cpp \
../src/NFfunction/muParser/muParserBytecode.cpp \
../src/NFfunction/muParser/muParserCallback.cpp \
../src/NFfunction/muParser/muParserComplex.cpp \
../src/NFfunction/muParser/muParserError.cpp \
../src/NFfunction/muParser/muParserInt.cpp \
../src/NFfunction/muParser/muParserTokenReader.cpp 

OBJS += \
./src/NFfunction/muParser/muParser.o \
./src/NFfunction/muParser/muParserBase.o \
./src/NFfunction/muParser/muParserBytecode.o \
./src/NFfunction/muParser/muParserCallback.o \
./src/NFfunction/muParser/muParserComplex.o \
./src/NFfunction/muParser/muParserError.o \
./src/NFfunction/muParser/muParserInt.o \
./src/NFfunction/muParser/muParserTokenReader.o 

CC_DEPS += \
./src/NFfunction/muParser/muParser.d \
./src/NFfunction/muParser/muParserBase.d \
./src/NFfunction/muParser/muParserBytecode.d \
./src/NFfunction/muParser/muParserCallback.d \
./src/NFfunction/muParser/muParserComplex.d \
./src/NFfunction/muParser/muParserError.d \
./src/NFfunction/muParser/muParserInt.d \
./src/NFfunction/muParser/muParserTokenReader.d 

src/NFfunction/muParser/%.o: ../src/NFfunction/muParser/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -Wno-strict-aliasing -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFfunction/compositeFunction.cpp \
../src/NFfunction/funcParser.cpp \
../src/NFfunction/function.cpp \
../src/NFfunction/localFunction.cpp 

OBJS += \
./src/NFfunction/compositeFunction.o \
./src/NFfunction/funcParser.o \
./src/NFfunction/function.o \
./src/NFfunction/localFunction.o 

CC_DEPS += \
./src/NFfunction/compositeFunction.d \
./src/NFfunction/funcParser.d \
./src/NFfunction/function.d \
./src/NFfunction/localFunction.d 

src/NFfunction/%.o: ../src/NFfunction/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

###############################################
# src/NFcore
###############################################

CC_SRCS += \
../src/NFcore/reactionSelector/directSelector.cpp \
../src/NFcore/reactionSelector/logClassSelector.cpp 

OBJS += \
./src/NFcore/reactionSelector/directSelector.o \
./src/NFcore/reactionSelector/logClassSelector.o 

CC_DEPS += \
./src/NFcore/reactionSelector/directSelector.d \
./src/NFcore/reactionSelector/logClassSelector.d 

src/NFcore/reactionSelector/%.o: ../src/NFcore/reactionSelector/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFcore/moleculeLists/moleculeList.cpp 

OBJS += \
./src/NFcore/moleculeLists/moleculeList.o 

CC_DEPS += \
./src/NFcore/moleculeLists/moleculeList.d 

src/NFcore/moleculeLists/%.o: ../src/NFcore/moleculeLists/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

CC_SRCS += \
../src/NFcore/complex.cpp \
../src/NFcore/complexList.cpp \
../src/NFcore/molecule.cpp \
../src/NFcore/moleculeType.cpp \
../src/NFcore/observable.cpp \
../src/NFcore/reactionClass.cpp \
../src/NFcore/system.cpp \
../src/NFcore/templateMolecule.cpp 

OBJS += \
./src/NFcore/complex.o \
./src/NFcore/complexList.o \
./src/NFcore/molecule.o \
./src/NFcore/moleculeType.o \
./src/NFcore/observable.o \
./src/NFcore/reactionClass.o \
./src/NFcore/system.o \
./src/NFcore/templateMolecule.o 

CC_DEPS += \
./src/NFcore/complex.d \
./src/NFcore/complexList.d \
./src/NFcore/molecule.d \
./src/NFcore/moleculeType.d \
./src/NFcore/observable.d \
./src/NFcore/reactionClass.d \
./src/NFcore/system.d \
./src/NFcore/templateMolecule.d 

src/NFcore/%.o: ../src/NFcore/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -Wno-unused-variable -c "$<" -o "$@"
	@echo ' '

###############################################
# src/ ...
###############################################

CC_SRCS += \
../src/NFsim.cpp 

OBJS += \
./src/NFsim.o 

CC_DEPS += \
./src/NFsim.d 

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	$(GG) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" $(WARNS) -c "$<" -o "$@"
	@echo ' '

all: NFsim

# linking
NFsim: $(OBJS)
	@echo 'Building target: $@'
	g++  -o "NFsim" $(OBJS)
	@echo 'Finished building target: $@'

clean:
	-$(RM) $(OBJS)$(C_DEPS)$(CC_DEPS) NFsim
	-@echo ' '

.PHONY: all clean
