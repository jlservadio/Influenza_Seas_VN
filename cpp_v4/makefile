# driver: driver.o rkf.o
# 	lcc  driver.o rkf.o -o driver -lm
# makefile
# Tabs *must* be used for the indentations below;
# spaces cause make syntax errors.

CC=g++
#CFLAGS=-fast -xO4 -xdepend -xarch=v8plusa -xprefetch -xvector -xunroll=8 -fsimple=2 -xsafe=mem
LIBS=-lm
GSLLIBS=-lgsl -lgslcblas 

ode:		
		$(CC) -O3 -c -o predict.o predict.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o func.o func.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o readdata.o readdata.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o prms.o prms.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o main_for_odesim.o main_for_odesim.cpp $(LIBS) $(GSLLIBS)
		$(CC) $(CFLAGS) -o odesim main_for_odesim.o predict.o func.o readdata.o prms.o  $(LIBS) $(GSLLIBS) 


clean:
		rm -f *.o core odesim*  *~




