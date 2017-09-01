SHELL = /bin/sh
PROJECT = swingmp;

CC = icpc
TARGET = swingmp
CFLAGS = -fopenmp -qopt-report=5 -xCORE-AVX2 -axMIC-AVX512 -O3  
LFLAGS = -qopt-report-phase=loop,vec
OFLAGS = -qopt-report-file=stderr

TFILE = $(PROJECT).tgz
ZFILE = $(PROJECT).zip
PACKFILES = Makefile swingmp.cpp


all:    $(TARGET)

run:    all
	- ./$(TARGET)

swingmp.o:        swingmp.cpp Makefile



$(TARGET):      
	$(CC) $(CFLAGS)  $(LFLAGS) $(OFLAGS) -o $(TARGET) $(TARGET).cpp 
clean:
	- /bin/rm -f *.o $(TARGET) *.optrpt 
#$(TFILE) $(ZFILE)



pack:
	tar zcvf $(TFILE) $(PACKFILES)
	zip -r $(ZFILE) $(PACKFILES)

###############################################################################
