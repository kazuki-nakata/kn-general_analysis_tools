SYSDIR=sys
include ARCH
include $(SYSDIR)/Sysdep_$(SYS)

SRCDIR=.$(DELIMITER)src$(DELIMITER)
BINDIR=.$(DELIMITER)bin$(DELIMITER)
LIBDIR=.$(DELIMITER)lib$(DELIMITER)
INCLUDEDIR=.$(DELIMITER)include$(DELIMITER)

OBJS = $(SRCDIR)*.o
MODS = $(SRCDIR)*.mod
TARGET = respack
LIB  = lib$(TARGET).a

MAIN = timeseries_analysis.F90
SUB = io_geotiff.o filter.o arc_array_arrangement.o sort.o sort2.o object_module.o object_module3.o data_extraction.o util_shadow_module.o util_texture_module.o

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
	$(AR) $(SRCDIR)$(LIB) $^
	$(RANLIB) $(SRCDIR)$(LIB)

$(BINDIR):
	$(MD) $(BINDIR)

$(LIBDIR):
	$(MD) $(LIBDIR)

$(INCLUDEDIR):
	$(MD) $(INCLUDEDIR)

install:
	$(MV) $(SRCDIR)$(LIB) $(LIBDIR) > nul
	$(MV) $(TARGET).pyd $(LIBDIR) > nul
	$(MV) $(MODS) $(INCLUDEDIR) > nul
	$(MV) $(OBJS) $(INCLUDEDIR) > nul

clean:
	$(RM) $(OBJS) $(MODS) $(SRCDIR)$(LIB) $(TARGET).pyd

uninstall:
	$(RM) $(LIBDIR)$(LIB) $(LIBDIR)$(TARGET).pyd $(INCLUDEDIR)*.o $(INCLUDEDIR)*.mod
#	rmdir $(BINDIR)
#	rmdir $(LIBDIR)
#	rmdir $(INCLUDEDIR)
