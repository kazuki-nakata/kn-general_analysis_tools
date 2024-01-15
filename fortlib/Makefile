SYSDIR=sys
include ARCH
include $(SYSDIR)/Sysdep_$(SYS)

SRCDIR=.$(DELIMITER)src$(DELIMITER)
BINDIR=.$(DELIMITER)bin$(DELIMITER)
LIBDIR=.$(DELIMITER)lib$(DELIMITER)
PYTHONDIR=..$(DELIMITER)knool$(DELIMITER)
INCLUDEDIR=.$(DELIMITER)include$(DELIMITER)

OBJS = $(SRCDIR)*.o
MODS = $(SRCDIR)*.mod
TARGET = fortlib
LIB  = lib$(TARGET).a

MAIN = filter.F90 similarity.F90 coregistration.F90 sensor_geometry.F90 grid_data.F90 pmw_processor.F90 pmw_processor_l1b.F90 deblur.F90
SUB = sort.o #pmw_processor.o #filter.o

MAINOBJS = $(addprefix $(SRCDIR), $(MAIN))
MODOBJS = $(addprefix $(SRCDIR), $(SUB))


.SUFFIXES:
.SUFFIXES: .F90 .o
.SUFFIXES: .F90 .mod
.SUFFIXES: .F90 .f90

.F90.o:
	$(FC) $(FFLAGS) $(OPT) -o $@ -J $(SRCDIR) -c $<

.F90.mod:
	$(FC) $(FFLAGS) $(OPT) -o $@ -J $(SRCDIR) -c $< 

.F90.f90:
	$(CPP) $(OPT)  $< > $@

all:$(TARGET)

$(MAINOBJS):$(LIB)

$(TARGET):$(MODOBJS) $(MAINOBJS)
	$(PC) -c $(PFFLAGS) $(LDFLAGS) $(OPT) -m $@ $(MAINOBJS) $(MODOBJS)

$(LIB):$(MODOBJS)
	$(AR) crvf $(SRCDIR)$(LIB) $^
	$(RANLIB) $(SRCDIR)$(LIB)

$(BINDIR):
	$(MD) $(BINDIR)

$(LIBDIR):
	$(MD) $(LIBDIR)

$(INCLUDEDIR):
	$(MD) $(INCLUDEDIR)

install:
	$(MV) $(SRCDIR)$(LIB) $(LIBDIR) > nul
	$(MV) $(TARGET)* $(PYTHONDIR) > nul
	$(MV) $(MODS) $(INCLUDEDIR) > nul
	$(MV) $(OBJS) $(INCLUDEDIR) > nul

clean:
	$(RM) $(OBJS) $(MODS) $(SRCDIR)$(LIB) $(TARGET).pyd

uninstall:
	$(RM) $(LIBDIR)$(LIB) $(PYTHONDIR)$(TARGET).pyd $(INCLUDEDIR)*.o $(INCLUDEDIR)*.mod
#	rmdir $(BINDIR)
#	rmdir $(LIBDIR)
#	rmdir $(INCLUDEDIR)
