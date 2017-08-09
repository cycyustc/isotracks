
OS="$(shell uname)"
ARCH="$(shell uname -m)"
UBUNTU="$(shell lsb_release -rs)"
ISTPAGB="$(shell ls tpagb || echo 0)"

ifeq ($(OS),"Linux")
#Linux:
  LIBRARIES=/usr/lib/libm.a #/usr/lib/libcfitsio.a
  ifeq ($(ARCH),"x86_64")
    LIBRARIES=/usr/lib/x86_64-linux-gnu/libm.a #/usr/lib/x86_64-linux-gnu/libcfitsio.a
  endif
  ifeq ($(UBUNTU),"12.04")
    LIBRARIES += -L/usr/lib -lm #-lpthread -lz
  endif
  ifeq ($(UBUNTU),"14.04")
    LIBRARIES += -L/usr/lib -lm #-lpthread -lz
  endif
  ifeq ($(UBUNTU),"16.10")
    LIBRARIES=-L/usr/lib/x86_64-linux-gnu/ -lm
  endif
# LIBRARIES=
#Mac:
else 
 LIBRARIES=-L/usr/lib -lm
endif

#Linux without graphics: 
LIBS=  

# COMPILER FLAGS 
#CFLAGS= -g
# for optimization:
#CFLAGS= -O3 
# for optimization + cfitsio:
CFLAGS= -g -O3 -w 

# LINKER FLAGS 
#LFLAGS= -g -O 
# for optimization:
#LFLAGS= -O3 --static
# for optimization + cfitsio:
LFLAGS= -g -O3 

# DIR CONTAINING FILES COMMON TO CMD AND TRILEGAL, 
# does not work properly, needs to split header files 
# instead of having symbolic links there
#VPATH = .:./common
#ifneq ($(ISTPAGB),"0") # adds commons for tpagb
#	VPATH= .:./common:./tpagb
#endif
#$(error VPATH is $(VPATH))

TARGET= main	
HEAD= parametri.h
COMMONHEAD= #abort.h #colors.h dataFileReader.h 
HEADER= $(HEAD) $(COMMONHEAD)
SRCS= main.c
COMMONSRCS=
OLDSRCS=
SOURCES= $(SRCS) $(COMMONSRCS)

OBJS= $(SOURCES:.c=.o)

#################################################
all: .depend $(TARGET)

.depend:
	$(CC) $(CFLAGS) -M $(SRCS) >> .depend

ifeq ($wildcard .depend,.depend)
 include .depend	
endif	

$(TARGET): $(OBJS) Makefile
	$(CC) -o $(TARGET) $(LFLAGS) $(OBJS) $(LIBS) $(LIBRARIES)

$(OBJS): $(HEADER) Makefile

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f *.o $(TARGET) .depend
