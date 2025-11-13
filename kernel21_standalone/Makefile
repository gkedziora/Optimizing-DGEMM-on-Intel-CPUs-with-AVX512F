# Compiler options
CC          = gcc
ICC         = icc
ICX         = icx
CRAYCC      = cc

# Base flags
BASE_FLAGS  = -lpthread -fopenmp
BASE_FLAGS_ICX = -lpthread -qopenmp

# MKL (optional) - override by sourcing environment or setting MKL_INC/MKL_LDFLAGS
MKL_INC ?= -I/nasa/intel/oneapi-2025.2.1/mkl/2025.2/include
MKL_LDFLAGS ?= -L/nasa/intel/oneapi-2025.2.1/mkl/2025.2/lib -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -DMKL_ILP64 -m64

# GCC optimization flags
GCC_FLAGS_O0    = -O0 -march=native $(BASE_FLAGS)
GCC_FLAGS_O2    = -O2 -march=native -ftree-vectorize $(BASE_FLAGS)
GCC_FLAGS_O3    = -O3 -march=native -ffast-math -funroll-loops -ftree-vectorize $(BASE_FLAGS)
GCC_FLAGS_ZNVER5 = -O3 -march=znver5 -ffast-math -funroll-loops -ftree-vectorize $(BASE_FLAGS)

# ICC optimization flags
ICC_FLAGS_O0    = -O0 -xHost $(BASE_FLAGS)
ICC_FLAGS_O2    = -O2 -xHost -qopt-report=5 $(BASE_FLAGS)
ICC_FLAGS_O3    = -O3 -xHost -qopt-report=5 -ffast-math -unroll-aggressive $(BASE_FLAGS)

# ICX optimization flags
#ICX_OPT_O0  = -O0 -march=znver5 $(BASE_FLAGS_ICX)
ICX_FLAGS_O0  = -O0 -march=common-avx512 $(BASE_FLAGS_ICX)
ICX_FLAGS_O2  = -O2 -march=common-avx512 -qopt-report $(BASE_FLAGS_ICX)
ICX_FLAGS_O3  = -O3 -march=common-avx512 $(BASE_FLAGS_ICX)
ICX_FLAGS_Ofast  = -Ofast -march=common-avx512 $(BASE_FLAGS_ICX)

# CRAYCC optimization flags
CRAYCC_FLAGS_O3  = -O3 -march=cray-znver5 $(BASE_FLAGS)
CRAYCC_FLAGS_fast  = -march=cray-znver5 -fast $(BASE_FLAGS)

SRC = dgemm_test.c utils.c

# Default build
build: gcc-O3

# GCC builds
gcc-O0: $(SRC)
	$(CC) -DUSE_MKL $(GCC_FLAGS_O0) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_gcc_O0

gcc-O2: $(SRC)
	$(CC) -DUSE_MKL $(GCC_FLAGS_O2) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_gcc_O2

gcc-O3: $(SRC)
	$(CC) -DUSE_MKL $(GCC_FLAGS_O3) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_gcc_O3

gcc-znver4: $(SRC)
	$(CC) -DUSE_MKL $(GCC_FLAGS_ZNVER4) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_gcc_znver4

# ICC builds
icc-O0: $(SRC)
	$(ICC) -DUSE_MKL $(ICC_FLAGS_O0) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_icc_O0

icc-O2: $(SRC)
	$(ICC) -DUSE_MKL $(ICC_FLAGS_O2) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_icc_O2

icc-O3: $(SRC)
	$(ICC) -DUSE_MKL $(ICC_FLAGS_O3) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_icc_O3

# ICX builds
icx-O2: $(SRC)
	$(ICX) -DUSE_MKL $(ICX_FLAGS_O2) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_icx_O2

icx-O3: $(SRC)
	$(ICX) -DUSE_MKL $(ICX_FLAGS_O3) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_icx_O3

icx-Ofast: $(SRC)
	$(ICX) -DUSE_MKL $(ICX_FLAGS_O3) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_icx_Ofast


# CRAYCC builds
craycc-O3: $(SRC)
	$(CRAYCC) $(CRAYCC_FLAGS_O3) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_CRAYCC_O3

craycc-fast: $(SRC)
	$(CRAYCC) $(CRAYCC_FLAGS_fast) $(MKL_INC) $(MKL_LDFLAGS) $(SRC) -o dgemm_CRAYCC_fast

# Build all variants
all: gcc-O0 gcc-O2 gcc-O3 gcc-znver4 icc-O0 icc-O2 icc-O3 CRAYCC-O0 CRAYCC-O2 CRAYCC-O3

# Clean
clean:
	rm -f dgemm_* *.optrpt

# Print compiler versions
info:
	@echo "GCC version:"
	@$(CC) --version
	@echo "\nICC version:"
	@$(ICC) --version 2>/dev/null || echo "ICC not available"
	@echo "\nCRAYCC version:"
	@$(CRAYCC) --version 2>/dev/null || echo "CRAYCC not available"
