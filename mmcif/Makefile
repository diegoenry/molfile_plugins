# Makefile for mmCIF VMD plugin
#
# Usage:
#   make                    # Build the plugin
#   make install            # Install to VMD plugin directory
#   make clean              # Remove built files
#   make config             # Show detected configuration
#   make test               # Verify compiled plugin binary
#
# Configuration:
#   Override variables via environment or command line:
#   make VMD_DIR=/path/to/vmd

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
        VMD_DIR := $(wildcard /Applications/VMD*.app/Contents/vmd)
        ifndef VMD_DIR
            VMD_DIR := $(wildcard $(HOME)/Applications/VMD*.app/Contents/vmd)
        endif
        ifndef VMD_DIR
            VMD_DIR := $(wildcard $(HOME)/software/VMD*.app/Contents/vmd*)
        endif
    else
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
# Use shell test to handle spaces in paths correctly
ifndef VMD_PLUGIN_INCLUDE
    VMD_HAS_LIB := $(shell test -d "$(VMD_DIR)/lib/plugins/include" && echo yes)
    ifeq ($(VMD_HAS_LIB),yes)
        VMD_PLUGIN_INCLUDE = $(VMD_DIR)/lib/plugins/include
        VMD_PLUGIN_DIR = $(VMD_DIR)/lib/plugins/$(VMD_ARCH)/molfile
    else
        VMD_PLUGIN_INCLUDE = $(VMD_DIR)/plugins/include
        VMD_PLUGIN_DIR = $(VMD_DIR)/plugins/$(VMD_ARCH)/molfile
    endif
else
    VMD_PLUGIN_DIR ?= $(VMD_DIR)/plugins/$(VMD_ARCH)/molfile
endif

# Compiler and flags
CXX ?= c++
CXXFLAGS = -O2 -fPIC -Wall -std=c++11

# Platform-specific linker flags
ifeq ($(UNAME_S),Darwin)
    LDFLAGS = -bundle -flat_namespace -undefined dynamic_lookup
else
    LDFLAGS = -shared
endif

# Target
PLUGIN = mmcifplugin.so

# Source files
SOURCES = mmcifplugin.cpp cifparse.cpp

# Default target
all: $(PLUGIN)

$(PLUGIN): $(SOURCES) cifparse.h periodic_table.h
	$(CXX) $(CXXFLAGS) -I'$(VMD_PLUGIN_INCLUDE)' -o $(PLUGIN) $(SOURCES) $(LDFLAGS)

clean:
	rm -f $(PLUGIN) *.o

install: $(PLUGIN)
	@echo "Installing plugin to $(VMD_PLUGIN_DIR)"
	@mkdir -p $(VMD_PLUGIN_DIR)
	cp $(PLUGIN) $(VMD_PLUGIN_DIR)/

# Show detected configuration
config:
	@echo "=== Build Configuration ==="
	@echo "Platform:           $(PLATFORM)"
	@echo "Architecture:       $(ARCH)"
	@echo "VMD Architecture:   $(VMD_ARCH)"
	@echo "VMD_DIR:            $(VMD_DIR)"
	@echo "VMD_PLUGIN_INCLUDE: $(VMD_PLUGIN_INCLUDE)"
	@echo "VMD_PLUGIN_DIR:     $(VMD_PLUGIN_DIR)"
	@echo "CXX:                $(CXX)"
	@echo "CXXFLAGS:           $(CXXFLAGS)"
	@echo "LDFLAGS:            $(LDFLAGS)"
	@echo ""
	@echo "To override, use: make VMD_DIR=/path/to/vmd"

# Test the compiled plugin binary
test: $(PLUGIN)
	@echo "Testing plugin compilation..."
	@echo "Plugin compiled successfully: $(PLUGIN)"
	@file $(PLUGIN)
	@echo ""
	@echo "Plugin size:"
	@ls -la $(PLUGIN)

.PHONY: all clean install test config
