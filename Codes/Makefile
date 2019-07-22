CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I/usr/local/include/ -std=c++11
LD_FLAGS	= -L/usr/local/lib -lm  -lgsl
QSTARSTAR   = qstarstar.cpp
BOTTLENECK   = bottleneck.cpp
AVG_BOTTLENECK   = avg_bottleneck.cpp

#Run files
all: qstarstar bottleneck avg_bottleneck	
qstarstar:	$(QSTARSTAR)
	$(CC) $(CC_FLAGS) $(QSTARSTAR) -o qstarstar.run $(LD_FLAGS)
bottleneck:	$(BOTTLENECK)
	$(CC) $(CC_FLAGS) $(BOTTLENECK) -o bottleneck.run $(LD_FLAGS)
avg_bottleneck:	$(AVG_BOTTLENECK)
	$(CC) $(CC_FLAGS) $(AVG_BOTTLENECK) -o avg_bottleneck.run $(LD_FLAGS)

#Other files


#Clean up all compiled files
clean:
	\rm -f *.o *~ *.run
