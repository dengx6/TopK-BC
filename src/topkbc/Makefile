
#help for make
#http://www.cnblogs.com/wang_yb/p/3990952.html
#https://segmentfault.com/a/1190000000349917
#http://blog.csdn.net/cuiyifang/article/details/7910268

#to use gprof to analyse efficience of the program:
#http://blog.chinaunix.net/uid-25194149-id-3215487.html
#to use doxygen+graphviz+htmlhelp to generate document from source code:
#http://www.doxygen.nl/
#(also include good comments norm)
#http://blog.csdn.net/u010740725/article/details/51387810

#TODO:the dependences are not complete!

#compile parameters

CC = g++
MPICC = mpicxx


#the optimazition level of gcc/g++
#http://blog.csdn.net/hit_090420216/article/details/44900215
#NOTICE: -O2 is recommended, while -O3 is dangerous
#when developing, not use -O because it will disturb the normal 
#routine. use it for test and release.
CFLAGS = -c -Wall -g -pg #-O2

CFLAGS += -DNO_THREAD
#CFLAGS += -DDEBUG_TRACK
#CFLAGS += -DMAIN_LOG
#CFLAGS += -DSTD_TRACK
#CFLAGS += -DSTD_MAIN
CFLAGS += -DQUERYSTR
CFLAGS += -DPTR_MATCH
CFLAGS += -DSPAN
#CFLAGS += -DI_AND_D
#CFLAGS += -DTUPLE
#CFLAGS += -DUNCOMP
#CFLAGS += -DCOLORING
#CFLAGS += -DREDCC
#CFLAGS += -DGLOBAL_COMMENT
#CFLAGS += -DSPACE_LOG
#CFLAGS += -DWHOLE_STR
#CFLAGS += -DALL_I
#CFLAGS += -DINVALID_READ
#CFLAGS += -DCYCLE_TRACK
#CFLAGS += -DSELECTIVE



EXEFLAG = -g -O2 -pg #-O2
LDFLAG = -pthread

# paths

homedir = XXXXXX

topkbcdir = $(homedir)TopkBC/

edgedir = $(homedir)edge/

gstreamdir = $(homedir)gstream/

objdir = $(homedir)objs/

exedir = $(homedir)bin/

# objects

topkbcobj = $(objdir)topkbc.o $(objdir)edge.o $(objdir)gstream.o

objfile = $(topkbcobj)

def64IO = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

all: topkbc

#topkbc
topkbc: $(objdir)topkbc_main.o $(objfile)
	$(CC) $(EXEFLAG) ${LDFLAG} -o $(exedir)topkbc $(objdir)topkbc_main.o $(objfile)


#<obj:topkbc>

$(objdir)topkbc_main.o: $(topkbcdir)topkbc_main.cpp 
	$(CC) $(CFLAGS) $(topkbcdir)topkbc_main.cpp -o $(objdir)topkbc_main.o 

$(objdir)topkbc.o: $(topkbcdir)topkbc.cpp $(topkbcdir)topkbc.h 
	$(CC) $(CFLAGS) $(topkbcdir)topkbc.cpp -o $(objdir)topkbc.o
	
$(objdir)edge.o: $(edgedir)edge.cpp $(edgedir)edge.h 
	$(CC) $(CFLAGS) $(edgedir)edge.cpp -o $(objdir)edge.o
	
$(objdir)gstream.o: $(gstreamdir)gstream.cpp $(gstreamdir)gstream.h 
	$(CC) $(CFLAGS) $(gstreamdir)gstream.cpp -o $(objdir)gstream.o

#</obj:topkbc>

#<obj:other>

#$(objdir)util.o: $(utildir)util.h $(utildir)util.cpp $(utildir)list_T.h
#	$(CC) $(CFLAGS) $(utildir)util.cpp -o $(objdir)util.o

#$(objdir)list.o:  $(utildir)list_T.h 
#	$(CC) $(CFLAGS) $(utildir)list_T.h -o $(objdir)list.o

#$(objdir)runtime.o:  $(utildir)runtime.h $(utildir)runtime.cpp
#	$(CC) $(CFLAGS) $(utildir)runtime.cpp -o $(objdir)runtime.o

#</obj:other>

clean:
	rm -rf $(exedir)* $(objdir)*.o
	#rm -rf .project .cproject .settings   just for eclipse
	#rm -rf cscope* just for vim


