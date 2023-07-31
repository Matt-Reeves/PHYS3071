SOURCEFILES = schrodinger.cpp #list of c++ files to compile
SUBDIRS = ./data #Where data will be stored
CC = /opt/homebrew/bin/g++-13 #location of our c++ compiler
INCL = -I/opt/homebrew/include #include directory
LIBS = -L/opt/homebrew/lib  -lfftw3 -lm #library directory and compiler library flags 
OPTFLAGS =  -O3 #optional compiler flags 
BINARY = schrodinger.out #executable file name

build : $(SOURCEFILES)
	$(CC) $(SOURCEFILES) $(OPTFLAGS) -o $(BINARY) $(INCL) $(LIBS)
	$(shell mkdir -p $(SUBDIRS)) 

.PHONY : clean
clean :
	$(RM) -r schrodinger.out *.dat *.txt ./data

rebuild: clean build


