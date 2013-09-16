SHELL=/bin/sh

.SUFFIXES:

AR   = ar cr
AT   = @
ECHO = $(AT)echo
ifdef IFORT
  F77 = ifort
else
  F77 = gfortran
endif

SOURCES = $(wildcard ?????.f)
# some stuff missing in libblas.a
ifndef MKL
  SOURCES += $(wildcard ??????.f)
endif

OBJDIR  = ./obj
DOBJDIR = ./dobj
OUTDIR  = ../../lib
OBJS    = $(patsubst %.f, $(OBJDIR)/%.o, $(SOURCES))
DOBJS   = $(patsubst %.f, $(DOBJDIR)/%.o, $(SOURCES))
FCOMMON_FLAGS = -fPIC
FFLAGS  = -O3 $(FCOMMON_FLAGS)
FFLAGSD = -g $(FCOMMON_FLAGS)

ifdef EM64T
  ifdef MKL
    LIBNAME  = libftn_mkl_em64t.a
    DLIBNAME = libftn_mkl_em64t_debug.a
  else
    LIBNAME  = libftn_em64t.a
    DLIBNAME = libftn_em64t_debug.a
  endif
else  
  ifdef MKL
    LIBNAME  = libftn_mkl.a
    DLIBNAME = libftn_mkl_debug.a
  else
    LIBNAME  = libftn.a
    DLIBNAME = libftn_debug.a
  endif
endif

all:     debug release
release: $(OUTDIR)/$(LIBNAME)
debug:   $(OUTDIR)/$(DLIBNAME)

$(OUTDIR)/$(LIBNAME): $(OBJS)
	$(ECHO) Archiving: $@
	$(AT)$(AR) $@ $(OBJS)

$(OUTDIR)/$(DLIBNAME): $(DOBJS)
	$(ECHO) Archiving: $@
	$(AT)$(AR) $@ $(DOBJS)

$(OBJDIR)/%.o: %.f
	$(AT)if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(ECHO) Compiling release: $<
	$(AT)$(F77) $(FFLAGS) -c -o $@ $<

# ifort optimizer is messing up with complex numbers returned from DOT functions, thus turning it off
ifdef IFORT
$(OBJDIR)/cmplx.o: cmplx.f
	$(AT)if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(ECHO) Compiling debug: $<
	$(AT)$(F77) $(FFLAGSD) -c -o $@ $<
endif

$(DOBJDIR)/%.o: %.f
	$(AT)if [ ! -d $(DOBJDIR) ]; then mkdir $(DOBJDIR); fi
	$(ECHO) Compiling debug: $<
	$(AT)$(F77) $(FFLAGSD) -c -o $@ $<

clean:
	$(ECHO) Cleaning...
	$(AT)rm -rf $(OBJDIR) $(DOBJDIR) $(OUTDIR)/libftn*
