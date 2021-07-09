CC = gcc

RM = rm -f
LAPACK_DIR = /opt/libflame/gcc-7.5.0/2.2/

LDFLAGS = -mcmodel=medium -L${LAPACK_DIR}/lib -llapack -lm
CFLAGS = -c -g -O2 -Wall -mcmodel=medium -I${LAPACK_DIR}/include
BUILD = afmpd
OBJS = afmpd.o

.PHONY: all clean

default: afmpd
all: $(BUILD)
clean:
	$(RM) *.o

afmpd: $(OBJS)
	$(CC) $^	$(LDLIBS) $(LDFLAGS) -o $@

afmpd.o:	$($@:.o=.c)	Makefile

.c.o:
	$(CC) $(CFLAGS) -o $@ $<
