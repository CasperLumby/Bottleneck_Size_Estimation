#NOTE: Some classes needs functionality from gcc version 4.5 or higher. Solve by running 'module load gcc/4.9.2' on Cambridge HPC cluster.
#NOTE: Some classes needs functionality from gcc version 4.5 or higher. Solve by running 'module load gcc/5.3.0' on Cambridge HPC cluster.
CC		= g++
CC_FLAGS	= -fopenmp -std=c++0x -g -O3 -Wall -I/usr/local/Cluster-Apps/gsl/1.9/include -lgsl -lgslcblas -I../cubature-master
LD_FLAGS	= -L/usr/local/Cluster-Apps/gsl/1.9/lib
ALL_FILES	= anaParam.o analyser.o analyserPH.o analyserPHWH.o analyserPHMG.o analyserPHMGWH.o analyserPHMR.o data.o dataPH.o dataPHgen.o dataPHMG.o dataPHMGgen.o dataPHMR.o diploidSequence.o ../cubature-master/hcubature.c inputParser.o misc.o path.o pathPH.o pathPHMG.o pathPHMR.o ../cubature-master/pcubature.c model.o modelMR.o sequence.o simParam.o simParamPH.o
SimCmdInput	= ${ALL_FILES} simDataCmdInput.cpp
SimRealObsCount	= ${ALL_FILES} simRealObsCount.cpp
SimRepCmdInput	= ${ALL_FILES} simRepDataCmdInput.cpp
OpenMPCheck	= openMPCheck.cpp
IntegralEval	= ${ALL_FILES} integralEvaluator.cpp
REALDATA	= ${ALL_FILES} realData.cpp
REALDATAWH	= ${ALL_FILES} realDataWH.cpp
MMATRIXCHECK	= ${ALL_FILES} mMatrixChecker.cpp

all: siminput simrealobs simrepinput openmpcheck inteval realdata
siminput: $(SimCmdInput)
	$(CC) $(LD_FLAGS) $(SimCmdInput) $(CC_FLAGS) -o SimDataCmdInput.run
simrealobs: $(SimRealObsCount)
	$(CC) $(LD_FLAGS) $(SimRealObsCount) $(CC_FLAGS) -o SimRealObsCount.run
simrepinput: $(SimRepCmdInput)
	$(CC) $(LD_FLAGS) $(SimRepCmdInput) $(CC_FLAGS) -o SimRepDataCmdInput.run
openmpcheck: ${OpenMPCheck}
	${CC} ${LD_FLAGS} ${OpenMPCheck} ${CC_FLAGS} -o openMPCheck.run
inteval: ${IntegralEval}
	${CC} ${LD_FLAGS} ${IntegralEval} ${CC_FLAGS} -o integralEvaluator.run
realdata: ${REALDATA}
	$(CC) $(LD_FLAGS) $(REALDATA) $(CC_FLAGS) -o realData.run
realdatawh: ${REALDATAWH}
	$(CC) $(LD_FLAGS) $(REALDATAWH) $(CC_FLAGS) -o realDataWH.run
mmatrix: $(MMATRIXCHECK)
	$(CC) $(LD_FLAGS) $(MMATRIXCHECK) $(CC_FLAGS) -o mMatrixChecker.run


anaParam.o: anaParam.cpp
	$(CC) $(CC_FLAGS) -c anaParam.cpp
analyser.o: analyser.cpp
	$(CC) $(CC_FLAGS) -c analyser.cpp
analyserPH.o: analyserPH.cpp
	$(CC) $(CC_FLAGS) -c analyserPH.cpp
analyserPHWH.o: analyserPHWH.cpp
	$(CC) $(CC_FLAGS) -c analyserPHWH.cpp
analyserPHMG.o: analyserPHMG.cpp
	$(CC) $(CC_FLAGS) -c analyserPHMG.cpp
analyserPHMGWH.o: analyserPHMGWH.cpp
	$(CC) $(CC_FLAGS) -c analyserPHMGWH.cpp
analyserPHMR.o: analyserPHMR.cpp
	$(CC) $(CC_FLAGS) -c analyserPHMR.cpp
data.o: data.cpp
	$(CC) $(CC_FLAGS) -c data.cpp
dataPH.o: dataPH.cpp
	$(CC) $(CC_FLAGS) -c dataPH.cpp
dataPHgen.o: dataPHgen.cpp
	$(CC) $(CC_FLAGS) -c dataPHgen.cpp
dataPHMG.o: dataPHMG.cpp
	$(CC) $(CC_FLAGS) -c dataPHMG.cpp
dataPHMGgen.o: dataPHMGgen.cpp
	$(CC) $(CC_FLAGS) -c dataPHMGgen.cpp
dataPHMR.o: dataPHMR.cpp
	$(CC) $(CC_FLAGS) -c dataPHMR.cpp
diploidSequence.o: diploidSequence.cpp
	$(CC) $(CC_FLAGS) -c diploidSequence.cpp
inputParser.o: inputParser.cpp
	$(CC) $(CC_FLAGS) -c inputParser.cpp
misc.o: misc.cpp
	${CC} ${CC_FLAGS} -c misc.cpp
path.o: path.cpp
	$(CC) $(CC_FLAGS) -c path.cpp
pathPH.o: pathPH.cpp
	$(CC) $(CC_FLAGS) -c pathPH.cpp
pathPHMG.o: pathPHMG.cpp
	$(CC) $(CC_FLAGS) -c pathPHMG.cpp
pathPHMR.o: pathPHMR.cpp
	$(CC) $(CC_FLAGS) -c pathPHMR.cpp
model.o: model.cpp
	$(CC) $(CC_FLAGS) -c model.cpp
modelMR.o: modelMR.cpp
	$(CC) $(CC_FLAGS) -c modelMR.cpp
sequence.o: sequence.cpp
	$(CC) $(CC_FLAGS) -c sequence.cpp
simParam.o: simParam.cpp
	$(CC) $(CC_FLAGS) -c simParam.cpp
simParamPH.o: simParamPH.cpp
	$(CC) $(CC_FLAGS) -c simParamPH.cpp



clean:
	\rm *.o *~ *.run

