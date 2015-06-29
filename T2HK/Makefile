
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

BIN = th13delta
OBJ = $(BIN).o

$(BIN): $(BIN).o
	gcc $(BIN).o -o $(BIN) $(LDFLAGS) $(local_LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ) *~ deltasig.dat th13delta.dat
