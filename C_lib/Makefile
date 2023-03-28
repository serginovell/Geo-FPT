CC    = gcc #run without -mp-6 on splinter
FLAGS        = #-std=gnu99 -Iinclude -I/usr/local/include
CFLAGS       = -fPIC -g -pedantic -Wall -Wextra -ggdb3 -O3 -ffast-math -fopenmp -I/usr/include/gsl -I/usr/lib/x86_64-linux-gnu/   


LFLAGS       =  -lgomp -lm -fopenmp -L/usr/lib/x86_64-linux-gnu/ -L/usr/include/gsl -lgsl -lgslcblas  
LDFLAGS      = -shared

DEBUGFLAGS   = -O0 -D _DEBUG
RELEASEFLAGS = -O2 -D NDEBUG -combine -fwhole-program

TARGET  = clib.so
SOURCES = $(wildcard *.c)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.c=.o)


all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LFLAGS)
