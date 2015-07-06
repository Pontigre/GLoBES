
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

BIN = th13delta deltasig
OBJ = $(BIN).o

all: $(BIN)

th13delta: th13delta.o
	gcc th13delta.o -o th13delta $(LDFLAGS) $(local_LDFLAGS)

deltasig: deltasig.o
	gcc deltasig.o -o deltasig $(LDFLAGS) $(local_LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm -f $(BIN) *.o *~ deltasig.dat th13delta.dat events.dat *.eps
