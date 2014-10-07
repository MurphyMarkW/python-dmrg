# System-specific parameters
# --------------- ----------
  # Automatic platform identification.
  PLATFORM = $(shell uname)

  # C and C++ compilers, respectively.
  GCC = gcc
  CXX = g++

  # Extra compiler flags to use.
  GCCEXTRA = -pthread -D WITH_LAPACK -D WITH_COMPLEX -funroll-loops -O3 -fstrict-aliasing -fnogcse -D_REENTRANT
  CXXEXTRA = -pthread -ftemplate-depth-150 -D WITH_LAPACK -D WITH_COMPLEX -funroll-loops -O3 -fstrict-aliasing -fno-gcse -D_REENTRANT

  # Library directors and libraries to use.
  LIBDIRS = ./ /excport/home/mark/lib/ /export/home/mark/lib64/
  LIBS = m c lapack cblas f77blas gfortran atlas

  # Source directors to use.
  INCDIRS = ./ /export/home/mark/include/ /export/home/mark/src/dmrg /export/home/mark/src/Tensor

  # Clean-up command.
  RM = rm -Rf


# Targets
# -------
  # C targets:
  CHEADERS = # C header files. *.h
  CSOURCES = # C source files. *.c
  COBJECTS = $(patsubst %.c,%.o,$(CSOURCES)) $(patsubst %.h,%.o,$(CHEADERS))

  # C++ targets:
  CPPHEADERS = # C++ header files. *.hpp
  CPPSOURCES = main.cpp # C++ source files. *.cpp
  CPPOBJECTS = $(patsubst %.cpp,%.o,$(CPPSOURCES)) $(patsubst %.hpp,%.o,$(CPPHEADERS))

  OBJECTS = $(COBJECTS) $(CPPOBJECTS)

  ONAME = main # Output file name.

  SHARED = False # Shared object flag. True/False


# Compiler / linker / loader flags
# -------- - ------ - ------ -----
LDFLAGS = $(patsubst %,-L%,$(LIBDIRS)) $(patsubst %,-l%,$(LIBS))
INCLUDE = $(patsubst %,-I%,$(INCDIRS))


ifeq ($(PLATFORM),Linux) # Linux system specific options.
  GCCFLAGS = $(INCLUDE) -DWITH_CBLAS
  CXXFLAGS = $(INCLUDE) -DWITH_CBLAS

  ifeq ($(SHARED),True) # Shared object flags.
    # Flag for compiling with position-independent-code.
    GCCFLAGS = $(GCCFLAGS) -fPIC
    CXXFLAGS = $(CXXFLAGS) -fPIC
    # Flag for compiling as a shared object.
    SYSFLAGS = $(SYSFLAGS) -shared
  endif
endif

ifeq ($(PLATFORM),Darwin) # OS X system specific options.
  # Turn all warnings on.
  GCCFLAGS = $(INCLUDE) -Wall
  CXXFLAGS = $(INCLUDE) -Wall

  ifeq ($(SHARED),True) # Shared object flags.
    # Flag for compiling with position-independent-code.
    GCCFLAGS = $(GCCFLAGS) -fPIC
    CXXFLAGS = $(CXXFLAGS) -fPIC
    # Flag for compiling as a dynamic library.
    SYSFLAGS = $(SYSFLAGS) -dynamiclib -Wl,-undefined,dynamic_lookup
  endif
endif

  

# Make Commands
# ---- --------
all: $(ONAME)

$(ONAME): $(OBJECTS)
	$(CXX) $(OBJECTS) $(SYSFLAGS) $(LDFLAGS) -o $(ONAME)

%.o: %.c
	$(GCC) $(GCCFLAGS) $(GCCEXTRA) -g -c $< -o $@

%.o: %.h
	$(GCC) $(GCCFLAGS) $(GCCEXTRA) -g -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CXXEXTRA) -g -c $< -o $@

%.o: %.hpp
	$(CXX) $(CXXFLAGS) $(CXXEXTRA) -g -c $< -o $@

clean:
	$(RM) $(OBJECTS) $(ONAME)

