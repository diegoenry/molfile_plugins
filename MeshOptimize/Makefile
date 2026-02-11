# Makefile for mesh_optimize — STL mesh smoother, decimator, and HDF5 converter
#
# Usage:
#   make                    # Build with HDF5 + OpenMP (auto-detected)
#   make HDF5=0             # Build without HDF5 (STL output only)
#   make OPENMP=0           # Build without OpenMP
#   make config             # Show detected configuration
#   make clean              # Remove built files
#
# Configuration:
#   Override variables via environment or command line:
#   make CXX=g++-15 HDF5_PREFIX=/opt/homebrew

# ──────────────── Platform detection ─────────────────────

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# ──────────────── Compiler detection ─────────────────────
# Apple Clang lacks OpenMP; prefer Homebrew GCC on macOS.

# Only auto-detect if CXX was not explicitly set by the user
ifeq ($(origin CXX),default)
    ifeq ($(UNAME_S),Darwin)
        # Apple Clang lacks OpenMP; prefer Homebrew GCC
        ifneq ($(shell which g++-15 2>/dev/null),)
            CXX := g++-15
        else ifneq ($(shell which g++-14 2>/dev/null),)
            CXX := g++-14
        else ifneq ($(shell which g++-13 2>/dev/null),)
            CXX := g++-13
        endif
    endif
endif

CXXFLAGS = -O2 -std=c++17 -Wall

# ──────────────── OpenMP detection ───────────────────────

OPENMP ?= 1
ifeq ($(OPENMP),1)
    # Test if the compiler supports -fopenmp
    OPENMP_OK := $(shell echo 'int main(){}' | $(CXX) -fopenmp -x c++ - -o /dev/null 2>/dev/null && echo 1)
    ifeq ($(OPENMP_OK),1)
        CXXFLAGS += -fopenmp
        LDFLAGS  += -fopenmp
    else
        $(info Note: Compiler does not support OpenMP, building single-threaded.)
    endif
endif

# ──────────────── HDF5 detection ─────────────────────────

HDF5 ?= 1
ifeq ($(HDF5),1)
    ifdef HDF5_PREFIX
        HDF5_INCLUDE ?= $(HDF5_PREFIX)/include
        HDF5_LIB     ?= $(HDF5_PREFIX)/lib
        HDF5_FOUND   := 1
    else
        # Try pkg-config
        HDF5_CFLAGS_PC  := $(shell pkg-config --cflags hdf5 2>/dev/null)
        HDF5_LDFLAGS_PC := $(shell pkg-config --libs hdf5 2>/dev/null)

        ifneq ($(HDF5_CFLAGS_PC),)
            HDF5_FOUND := 1
        else
            # Try common locations
            ifeq ($(UNAME_S),Darwin)
                ifeq ($(UNAME_M),arm64)
                    HDF5_INCLUDE ?= /opt/homebrew/include
                    HDF5_LIB     ?= /opt/homebrew/lib
                else
                    HDF5_INCLUDE ?= /usr/local/include
                    HDF5_LIB     ?= /usr/local/lib
                endif
            else
                HDF5_INCLUDE ?= /usr/include
                HDF5_LIB     ?= /usr/lib
            endif
            # Verify hdf5.h exists
            HDF5_FOUND := $(shell test -f "$(HDF5_INCLUDE)/hdf5.h" && echo 1)
        endif
    endif

    ifeq ($(HDF5_FOUND),1)
        ifneq ($(HDF5_CFLAGS_PC),)
            CXXFLAGS += $(HDF5_CFLAGS_PC)
            LDFLAGS  += $(HDF5_LDFLAGS_PC)
        else
            CXXFLAGS += -I"$(HDF5_INCLUDE)"
            LDFLAGS  += -L"$(HDF5_LIB)" -lhdf5
        endif
    else
        $(info Note: HDF5 not found, building without .h5mesh output support.)
        $(info       Set HDF5_PREFIX or install HDF5 to enable.)
    endif
endif

# ──────────────── Target ─────────────────────────────────

TARGET = mesh_optimize

all: $(TARGET)

$(TARGET): mesh_optimize.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) mesh_optimize.cpp $(LDFLAGS)

clean:
	rm -f $(TARGET)

# ──────────────── Install ────────────────────────────────

PREFIX ?= /usr/local
install: $(TARGET)
	install -d $(PREFIX)/bin
	install -m 755 $(TARGET) $(PREFIX)/bin/

# ──────────────── Configuration ──────────────────────────

config:
	@echo "=== Build Configuration ==="
	@echo "Platform:     $(UNAME_S) $(UNAME_M)"
	@echo "Compiler:     $(CXX)"
	@echo "CXXFLAGS:     $(CXXFLAGS)"
	@echo "LDFLAGS:      $(LDFLAGS)"
ifeq ($(OPENMP_OK),1)
	@echo "OpenMP:       yes"
else
	@echo "OpenMP:       no"
endif
ifeq ($(HDF5_FOUND),1)
	@echo "HDF5:         yes"
  ifneq ($(HDF5_CFLAGS_PC),)
	@echo "  (via pkg-config)"
  else
	@echo "  include:    $(HDF5_INCLUDE)"
	@echo "  lib:        $(HDF5_LIB)"
  endif
else
	@echo "HDF5:         no (STL output only)"
endif
	@echo ""
	@echo "Override: make CXX=g++-15 HDF5_PREFIX=/opt/homebrew"

.PHONY: all clean install config
