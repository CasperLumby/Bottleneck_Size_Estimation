CC	      = g++
CC_FLAGS	= -g3 -O3 -Wall -I/usr/local/include/ -std=c++11
LD_FLAGS	= -L/usr/local/lib -lm  -lgsl
MLHAPREC   = MLHapRec.cpp
	
#Run files
MLHapRec:	$(MLHAPREC)
	$(CC) $(CC_FLAGS) $(MLHAPREC) -o run_MLHapRec $(LD_FLAGS)

#Other files


#Clean up all compiled files
clean:
	\rm -f *.o *~ run_*

