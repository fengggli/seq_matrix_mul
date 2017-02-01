## Compiler, tools and options
CC      = cc
CCFLAGS =-O0 -Wall -Wextra


## PAPI Setup
PAPI_ROOT  = /opt/cray/papi/5.3.0
PAPI_LIB   = $(PAPI_ROOT)/lib/libpapi.a

INCPATH    = -I. -I$(PAPI_ROOT)/include
LIBS       = -lm $(PAPI_LIB)

## Files
OBJECTS = seq_matrix_mul.o papi_timer.o
TARGET  = run_mul

## Implicit rules
.SUFFIXES: .c

.c.o:
	$(CC) -c $(CCFLAGS) $(INCPATH) $<

## Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core

