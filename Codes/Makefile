CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I/usr/local/include/ -std=c++11
LD_FLAGS	= -L/usr/local/lib -lm  -lgsl
QSTARSTAR   = qstarstar.cpp
BOTTLENECK   = bottleneck.cpp
AVG_BOTTLENECK   = avg_bottleneck.cpp

#Run files
all: qstarstar bottleneck avg_bottleneck mlhaprec
qstarstar:	$(QSTARSTAR)
	$(CC) $(CC_FLAGS) $(QSTARSTAR) -o run_qstarstar $(LD_FLAGS)
bottleneck:	$(BOTTLENECK)
	$(CC) $(CC_FLAGS) $(BOTTLENECK) -o run_bottleneck $(LD_FLAGS)
avg_bottleneck:	$(AVG_BOTTLENECK)
	$(CC) $(CC_FLAGS) $(AVG_BOTTLENECK) -o run_avg_bottleneck $(LD_FLAGS)
mlhaprec:
	cd MLHapRec && make

#Other files


#Clean up all compiled files
clean: cleanbottleneck cleanmlhaprec
cleanbottleneck:
	\rm -f *.o *~ run_*
cleanmlhaprec:
	cd MLHapRec && make clean
