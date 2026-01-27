# Makefile for OpenMiChroM CNDB VMD plugin
#
# Usage:
#   make                    # Build the plugin
#   make install            # Install to VMD plugin directory
#   make clean              # Remove built files
#
# Configuration:
#   Override variables via environment or command line:
#   make VMD_DIR=/path/to/vmd HDF5_PREFIX=/usr/local
#
# Required environment for non-standard installations:
#   VMD_DIR      - Path to VMD installation directory
#   HDF5_PREFIX  - Prefix where HDF5 is installed (contains include/ and lib/)
#   Or set HDF5_INCLUDE and HDF5_LIB separately

# Detect operating system and architecture
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Platform detection
ifeq ($(UNAME_S),Darwin)
    PLATFORM := MACOSX
    ifeq ($(UNAME_M),arm64)
        ARCH := ARM64
    else
        ARCH := X86_64
    endif
else ifeq ($(UNAME_S),Linux)
    PLATFORM := LINUX
    ifeq ($(UNAME_M),x86_64)
        ARCH := AMD64
    else ifeq ($(UNAME_M),aarch64)
        ARCH := ARM64
    else
        ARCH := $(UNAME_M)
    endif
else
    $(error Unsupported platform: $(UNAME_S))
endif

VMD_ARCH := $(PLATFORM)$(ARCH)

# VMD plugin directories - try to auto-detect common locations
ifndef VMD_DIR
    ifeq ($(UNAME_S),Darwin)
        # Try common macOS locations
        VMD_DIR := $(wildcard /Applications/VMD*.app/Contents/vmd)
        ifndef VMD_DIR
            VMD_DIR := $(wildcard $(HOME)/Applications/VMD*.app/Contents/vmd)
        endif
        ifndef VMD_DIR
            VMD_DIR := $(wildcard $(HOME)/software/VMD*.app/Contents/vmd*)
        endif
    else
        # Try common Linux locations
        VMD_DIR := $(wildcard /usr/local/lib/vmd)
        ifndef VMD_DIR
            VMD_DIR := $(wildcard $(HOME)/vmd)
        endif
    endif
    ifndef VMD_DIR
        $(warning VMD_DIR not set and could not be auto-detected. Set VMD_DIR environment variable.)
        VMD_DIR := /path/to/vmd
    endif
endif

# VMD directory structure varies - try both common layouts
ifeq ($(wildcard $(VMD_DIR)/lib/plugins/include),)
    # Direct plugins directory (VMD 1.9.4+)
    VMD_PLUGIN_INCLUDE = $(VMD_DIR)/plugins/include
    VMD_PLUGIN_DIR = $(VMD_DIR)/plugins/$(VMD_ARCH)/molfile
else
    # lib/plugins structure (older VMD)
    VMD_PLUGIN_INCLUDE = $(VMD_DIR)/lib/plugins/include
    VMD_PLUGIN_DIR = $(VMD_DIR)/lib/plugins/$(VMD_ARCH)/molfile
endif

# HDF5 detection - try pkg-config first, then common locations
ifdef HDF5_PREFIX
    HDF5_INCLUDE ?= $(HDF5_PREFIX)/include
    HDF5_LIB ?= $(HDF5_PREFIX)/lib
else
    # Try pkg-config
    HDF5_CFLAGS := $(shell pkg-config --cflags hdf5 2>/dev/null)
    HDF5_LDFLAGS := $(shell pkg-config --libs hdf5 2>/dev/null)

    ifdef HDF5_CFLAGS
        # pkg-config found HDF5
        HDF5_INCLUDE :=
        HDF5_LIB :=
    else
        # Try common locations based on platform
        ifndef HDF5_INCLUDE
            ifeq ($(UNAME_S),Darwin)
                # Homebrew paths
                ifeq ($(UNAME_M),arm64)
                    HDF5_INCLUDE ?= /opt/homebrew/include
                    HDF5_LIB ?= /opt/homebrew/lib
                else
                    HDF5_INCLUDE ?= /usr/local/include
                    HDF5_LIB ?= /usr/local/lib
                endif
            else
                # Linux standard paths
                HDF5_INCLUDE ?= /usr/include
                HDF5_LIB ?= /usr/lib
            endif
        endif
    endif
endif

# Compiler and flags
CC ?= cc
CFLAGS = -O2 -fPIC -Wall

# Add VMD include (quote to handle spaces in paths)
CFLAGS += -I"$(VMD_PLUGIN_INCLUDE)"

# Add HDF5 includes
ifdef HDF5_CFLAGS
    CFLAGS += $(HDF5_CFLAGS)
else
    CFLAGS += -I"$(HDF5_INCLUDE)"
endif

# Platform-specific linker flags
ifeq ($(UNAME_S),Darwin)
    LDFLAGS = -bundle -flat_namespace -undefined suppress
else
    LDFLAGS = -shared
endif

# Add HDF5 library (quote to handle spaces in paths)
ifdef HDF5_LDFLAGS
    LDFLAGS += $(HDF5_LDFLAGS)
else
    LDFLAGS += -L"$(HDF5_LIB)" -lhdf5
endif

# Target
PLUGIN = cndbplugin.so

# Default target
all: $(PLUGIN)

# Show detected configuration
config:
	@echo "=== Build Configuration ==="
	@echo "Platform:          $(PLATFORM)"
	@echo "Architecture:      $(ARCH)"
	@echo "VMD Architecture:  $(VMD_ARCH)"
	@echo "VMD_DIR:           $(VMD_DIR)"
	@echo "VMD_PLUGIN_INCLUDE: $(VMD_PLUGIN_INCLUDE)"
	@echo "VMD_PLUGIN_DIR:    $(VMD_PLUGIN_DIR)"
	@echo "HDF5_INCLUDE:      $(HDF5_INCLUDE)"
	@echo "HDF5_LIB:          $(HDF5_LIB)"
	@echo "CC:                $(CC)"
	@echo "CFLAGS:            $(CFLAGS)"
	@echo "LDFLAGS:           $(LDFLAGS)"
	@echo ""
	@echo "To override, use: make VMD_DIR=/path/to/vmd HDF5_PREFIX=/path/to/hdf5"

$(PLUGIN): cndbplugin.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(PLUGIN) cndbplugin.c

clean:
	rm -f $(PLUGIN) *.o

install: $(PLUGIN)
	@echo "Installing plugin to $(VMD_PLUGIN_DIR)"
	@mkdir -p $(VMD_PLUGIN_DIR)
	cp $(PLUGIN) $(VMD_PLUGIN_DIR)/

test: $(PLUGIN)
	@echo "Testing plugin compilation..."
	@echo "Plugin compiled successfully: $(PLUGIN)"
	@file $(PLUGIN)

.PHONY: all clean install test config
