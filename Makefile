CC = g++
CFLAGS = -g -O3 -Wfatal-errors -Wno-deprecated -D__USE_MINGW_ANSI_STDIO
INCLUDEDIRS = -I./SOCO_SI -I./fastFractal-CEC2008
OBJS = lmcma.o 
F7OBJS = RanQD1.o  RanTable.o  FastFractal.o  DoubleDip.o  FractalFunction1D.o UnitFunction1D.o
FUNCTOBJS = funsoft.o

lshade.exe : $(OBJS)  $(FUNCTOBJS) $(F7OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(FUNCTOBJS) $(F7OBJS) -o lmcma.exe

$(OBJS) : lmcma.cpp 
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c lmcma.cpp 
	
$(FUNCTOBJS) : SOCO_SI/funsoft.cpp
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c SOCO_SI/funsoft.cpp
  
$(F7OBJS) : ./fastFractal-CEC2008/RanQD1.cpp  ./fastFractal-CEC2008/RanTable.cpp  ./fastFractal-CEC2008/FastFractal.cpp  ./fastFractal-CEC2008/DoubleDip.cpp  ./fastFractal-CEC2008/FractalFunction1D.cpp ./fastFractal-CEC2008/UnitFunction1D.cpp
  $(CC) $(CFLAGS) $(INCLUDEDIRS) -c ./fastFractal-CEC2008/RanQD1.cpp  ./fastFractal-CEC2008/RanTable.cpp  ./fastFractal-CEC2008/FastFractal.cpp  ./fastFractal-CEC2008/DoubleDip.cpp  ./fastFractal-CEC2008/FractalFunction1D.cpp ./fastFractal-CEC2008/UnitFunction1D.cpp

	
.PHONY: clean

clean:
	rm -rf *.o