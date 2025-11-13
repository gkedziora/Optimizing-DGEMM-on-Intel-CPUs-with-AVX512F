BINARY_NAME = dgemm_x86
CC          = gcc
ICX         = icx
CRAYCC      = cc

# Base flags
BASE_FLAGS = -lpthread -fopenmp
BASE_FLAGS_ICX = -lpthread -qopenmp

# MKL configuration
#MKLPATH     = /opt/intel/mkl
LDFLAGS     = -L/nasa/intel/oneapi-2025.2.1/mkl/2025.2/lib -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -DMKL_ILP64 -m64
INCFLAGS    = -I/nasa/intel/oneapi-2025.2.1/mkl/2025.2/include

# GCC optimization flags
GCC_OPT_O0  = -O0 -march=native $(BASE_FLAGS)
GCC_OPT_O2  = -O2 -march=native -ftree-vectorize $(BASE_FLAGS)
GCC_OPT_O3  = -O3 -march=native -ffast-math -funroll-loops -ftree-vectorize $(BASE_FLAGS)

# ICX optimization flags
ICX_OPT_O0  = -O0 -march=znver5 $(BASE_FLAGS_ICX)
ICX_OPT_O2  = -O2 -march=common-avx512 $(BASE_FLAGS_ICX)
ICX_OPT_O3  = -O3 -march=common-avx512 $(BASE_FLAGS_ICX)
 
# CRAYCC optimization flags
CRAYCC_OPT_O0 = -O0 -march=znver5 $(BASE_FLAGS)
CRAYCC_OPT_O2 = -O2 -march=cray-znver5 -fvectorize $(BASE_FLAGS)
#CRAYCC_OPT_O3 = -O3 -march=znver5 -ffast-math -funroll-loops -fvectorize $(BASE_FLAGS)
CRAYCC_OPT_O3 = -O3 -march=cray-znver5  -ffast-math -funroll-loops -fvectorize $(BASE_FLAGS)

SRC         = $(wildcard *.c)

# Default build (no optimization)
build: gcc-O0

# GCC builds with different optimization levels
gcc-O0: $(SRC)
	$(CC) $(GCC_OPT_O0) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_gcc_O0

gcc-O2: $(SRC)
	$(CC) $(GCC_OPT_O2) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_gcc_O2

gcc-O3: $(SRC)
	$(CC) $(GCC_OPT_O3) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_gcc_O3

# ICX builds with different optimization levels
ICX-O0: $(SRC)
	$(ICX) $(ICX_OPT_O0) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_ICX_O0

ICX-O2: $(SRC)
	$(ICX) $(ICX_OPT_O2) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_ICX_O2

ICX-O3: $(SRC)
	$(ICX) $(ICX_OPT_O3) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_ICX_O3

# CRAYCC builds with different optimization levels
CRAYCC-O0: $(SRC)
	$(CRAYCC) $(CRAYCC_OPT_O0) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_CRAYCC_O0

CRAYCC-O2: $(SRC)
	$(CRAYCC) $(CRAYCC_OPT_O2) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_CRAYCC_O2

CRAYCC-O3: $(SRC)
	$(CRAYCC) $(CRAYCC_OPT_O3) $(LDFLAGS) $(INCFLAGS) $(SRC) -o $(BINARY_NAME)_CRAYCC_O3

# Build all variants
all: gcc-O0 gcc-O2 gcc-O3 ICX-O0 ICX-O2 ICX-O3 CRAYCC-O0 CRAYCC-O2 CRAYCC-O3

# Clean all builds
clean:
	rm -f $(BINARY_NAME)_* *.optrpt

# Print compiler versions
compiler-info:
	@echo "GCC version:"
	@$(CC) --version
	@echo "\nICX version:"
	@$(ICX) --version 2>/dev/null || echo "ICX not available"
	@echo "\nCRAYCC version:"
	@$(CRAYCC) --version 2>/dev/null || echo "CRAYCC not available"
