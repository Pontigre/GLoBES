
prefix = /usr/local
exec_prefix = /usr/local
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O4

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
preopen_modules:=$(shell $(globesconf) --dlpreopen)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)
ifdef preopen_modules
predefs = -DGLB_STATIC
endif

BIN = main th13delta th23dm32 deltasig events test_all test1 test2 test3 test4 1Dproj
OBJ = $(BIN).o

all: $(BIN)

main: main.o
	gcc main.o -o main $(LDFLAGS) $(local_LDFLAGS)

th13delta: th13delta.o
	gcc th13delta.o -o th13delta $(LDFLAGS) $(local_LDFLAGS)

th23dm32: th23dm32.o
	gcc th23dm32.o -o th23dm32 $(LDFLAGS) $(local_LDFLAGS)

deltasig: deltasig.o
	gcc deltasig.o -o deltasig $(LDFLAGS) $(local_LDFLAGS)

events: events.o
	gcc events.o -o events $(LDFLAGS) $(local_LDFLAGS)

1Dproj: 1Dproj.o
	gcc 1Dproj.o -o 1Dproj $(LDFLAGS) $(local_LDFLAGS)


test_all: test_all.o
	gcc test_all.o -o test_all $(LDFLAGS) $(local_LDFLAGS)

test1: test1.o
	gcc test1.o -o test1 $(LDFLAGS) $(local_LDFLAGS)

test2: test2.o
	gcc test2.o -o test2 $(LDFLAGS) $(local_LDFLAGS)

test3: test3.o
	gcc test3.o -o test3 $(LDFLAGS) $(local_LDFLAGS)

test4: test4.o
	gcc test4.o -o test4 $(LDFLAGS) $(local_LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm -f $(BIN) *.o *~ deltasig.dat th13delta.dat th13delta_05xerror.dat events_* init_flux* *.eps rates* test1.dat test2.dat test3.dat test4a.dat test4b.dat test4c.dat test4d.dat th13.dat th23.dat th23dm32.dat th23dm32_05xerror.dat

