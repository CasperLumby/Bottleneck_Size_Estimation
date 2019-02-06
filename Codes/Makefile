CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I/usr/local/include/ -std=c++11
LD_FLAGS	= -L/usr/local/lib -lm  -lgsl
QSTARSTAR   = qstarstar.o
BOTTLENECK   = bottleneck.o
AVG_BOTTLENECK   = avg_bottleneck.o
	
qstarstar:	$(QSTARSTAR)
	$(CC) $(CC_FLAGS) $(QSTARSTAR) -o run_qstarstar $(LD_FLAGS)
qstarstar.o: qstarstar.cpp
	$(CC) $(CC_FLAGS) -c qstarstar.cpp
bottleneck:	$(BOTTLENECK)
	$(CC) $(CC_FLAGS) $(BOTTLENECK) -o run_bottleneck $(LD_FLAGS)
bottleneck.o: bottleneck.cpp
	$(CC) $(CC_FLAGS) -c bottleneck.cpp
avg_bottleneck:	$(AVG_BOTTLENECK)
	$(CC) $(CC_FLAGS) $(AVG_BOTTLENECK) -o run_avg_bottleneck $(LD_FLAGS)
avg_bottleneck.o: avg_bottleneck.cpp
	$(CC) $(CC_FLAGS) -c avg_bottleneck.cpp


