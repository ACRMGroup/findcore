CC = gcc -O3 -ansi -Wall -pedantic
COPT = -I$(HOME)/include
LOPT = -L$(HOME)/lib
LIBS = -lbiop -lgen -lm -lxml2

TARGETS = profitcore findcore findcora

all : $(TARGETS)

profitcore : profitcore.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

findcore : findcore.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

findcora : findcora.o
	$(CC) $(LOPT) -o $@ $< $(LIBS)

.c.o :
	$(CC) $(COPT) -o $@ -c $<

clean :
	\rm -f *.o

distclean : clean
	\rm -f $(TARGETS)


