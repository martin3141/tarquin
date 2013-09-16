SHELL=/bin/sh

.SUFFIXES:

LIBS=-lpthread -lpthread_nonshared
DLIBS=-lpthread -lpthread_nonshared

AR   = ar cr
AT   = @
ECHO = $(AT)echo
ifdef ICCT
  CC=icpc
else
  CC=g++
endif

COMMON_FLAGS = -Wall
ifdef CVM_FLOAT
  COMMON_FLAGS += -DCVM_FLOAT
endif

ifdef CVM_STLPORT
  COMMON_FLAGS += -I/usr/include/stlport
  LIBS  += -lstlport
  DLIBS += -lstlportg
endif

ifndef ICCT
  LIBS  += -lgfortran
  DLIBS += -lgfortran
endif

DOBJDIR=./dobj
OBJDIR=./obj
OUTDIR=../../lib

CVM_LIB_PATH=-L$(OUTDIR)
ifndef MKL_PATH
  MKL_PATH=/opt/intel/mkl/10.0.011/lib/
endif

ifdef EM64T
  LIBS  += -lcvm_em64t
  DLIBS += -lcvm_em64t_debug
  OUTNAME_RELEASE=$(OUTDIR)/regtest_cvm_em64t
  OUTNAME_DEBUG=$(OUTDIR)/regtest_cvm_em64t_debug
  ifdef STATIC_ONLY
    ifdef MKL
      LIBS  += -lftn_mkl_em64t
      DLIBS += -lftn_mkl_em64t_debug
    else
      LIBS  += -lftn_em64t
      DLIBS += -lftn_em64t_debug
    endif
  endif
else
  LIBS  += -lcvm
  DLIBS += -lcvm_debug
  OUTNAME_RELEASE=$(OUTDIR)/regtest_cvm_ia32
  OUTNAME_DEBUG=$(OUTDIR)/regtest_cvm_ia32_debug
  ifdef STATIC_ONLY
    ifdef MKL
      LIBS  += -lftn_mkl
      DLIBS += -lftn_mkl_debug
    else
      LIBS  += -lftn
      DLIBS += -lftn_debug
    endif
  endif
endif

ifdef STATIC_ONLY
  ifdef MKL
    # Intel MKL 8.1 library:
    ifdef EM64T
      INTEL_MKL_PATH=-L$(MKL_PATH)/em64t
      LIBS  += -lmkl_lapack -lmkl_em64t
      DLIBS += -lmkl_lapack -lmkl_em64t
    else
      INTEL_MKL_PATH=-L$(MKL_PATH)/32
      LIBS  += -lmkl_lapack -lmkl_ia32
      DLIBS += -lmkl_lapack -lmkl_ia32
    endif
  else
    LIBS  += -llapack -lblas
    DLIBS += -llapack -lblas
#    DLIBS += -llapack_debug -lblas_debug
  endif
endif

CFLAGS=$(COMMON_FLAGS)
CFLAGSD=$(COMMON_FLAGS) -g -DDEBUG

ifdef PROF
  CFLAGSD += -pg
  LFLAGSD += -pg
endif

#$(call make-depend,source-file,object-file,depend-file)
define make-depend
  $(AT)$(CC) -MM -MF $3 -MP -MT $2 $(COMMON_FLAGS) $1
endef

SOURCES=$(wildcard *.cpp)
OBJS=$(patsubst %.cpp, $(OBJDIR)/%.o, $(SOURCES))
DOBJS=$(patsubst %.cpp, $(DOBJDIR)/%.o, $(SOURCES))


all: debug release
release : $(OUTNAME_RELEASE)
debug : $(OUTNAME_DEBUG)


$(OUTNAME_RELEASE) : $(OBJS)
	$(ECHO) Linking: $@
	$(AT)$(CC) -o $@ $(OBJS)  $(LFLAGS) $(CVM_LIB_PATH) $(INTEL_MKL_PATH) -Wl,--start-group $(LIBS) -Wl,--end-group

$(OUTNAME_DEBUG) : $(DOBJS)
	$(ECHO) Linking: $@
	$(AT)$(CC) $(LFLAGSD) -o $@ $(DOBJS) $(LFLAGSD) $(CVM_LIB_PATH) $(INTEL_MKL_PATH) -Wl,--start-group $(DLIBS) -Wl,--end-group

$(OBJDIR)/%.o: %.cpp
	$(AT)if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(ECHO) Compiling: $<
	$(AT)$(CC) $(CFLAGS) -c -o $@ $<

$(DOBJDIR)/%.o: %.cpp
	$(AT)if [ ! -d $(DOBJDIR) ]; then mkdir $(DOBJDIR); fi
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(ECHO) Compiling: $<
	$(AT)$(CC) $(CFLAGSD) -c -o $@ $<

clean:
	$(ECHO) Cleaning...
	$(AT)-rm -rf $(OBJDIR) $(DOBJDIR) $(OUTDIR)/regtest_cvm* $(OUTDIR)/*testout*.txt

ifneq "$(MAKECMDGOALS)" "clean"
  -include $(OBJDIR)/*.d $(DOBJDIR)/*.d
endif
