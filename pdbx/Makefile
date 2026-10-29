# Makefile for the PDBx/mmCIF VMD molfile plugin (pdbx 1.0)
#
#   make                 build pdbxplugin.so
#   make install         copy it into VMD's molfile plugin directory
#   make test            build the test driver and run it on tests/data/*.cif
#   make asan            same driver built with AddressSanitizer + UBSan
#   make config          show the detected paths
#   make clean
#
# Override detection with e.g.  make VMD_DIR=/path/to/vmd
# or point straight at the headers with  make VMD_PLUGIN_INCLUDE=/path/include

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),arm64)
        VMD_ARCH := MACOSXARM64
    else
        VMD_ARCH := MACOSXX86_64
    endif
else ifeq ($(UNAME_S),Linux)
    ifeq ($(UNAME_M),x86_64)
        VMD_ARCH := LINUXAMD64
    else ifeq ($(UNAME_M),aarch64)
        VMD_ARCH := LINUXARM64
    else
        VMD_ARCH := LINUX
    endif
else
    $(error Unsupported platform: $(UNAME_S))
endif

# VMD installation. App bundle paths contain spaces, so detection uses the
# shell rather than $(wildcard).
ifndef VMD_DIR
    ifeq ($(UNAME_S),Darwin)
        VMD_DIR := $(shell for d in /Applications/VMD*.app/Contents/vmd* \
                     $(HOME)/Applications/VMD*.app/Contents/vmd* \
                     $(HOME)/software/VMD*.app/Contents/vmd*; do \
                     if [ -d "$$d" ]; then echo "$$d"; break; fi; done)
    else
        VMD_DIR := $(shell for d in /usr/local/lib/vmd $(HOME)/vmd; do \
                     if [ -d "$$d" ]; then echo "$$d"; break; fi; done)
    endif
endif

ifndef VMD_PLUGIN_INCLUDE
    VMD_HAS_LIB := $(shell test -d "$(VMD_DIR)/lib/plugins/include" && echo yes)
    ifeq ($(VMD_HAS_LIB),yes)
        VMD_PLUGIN_INCLUDE = $(VMD_DIR)/lib/plugins/include
        VMD_PLUGIN_DIR ?= $(VMD_DIR)/lib/plugins/$(VMD_ARCH)/molfile
    else
        VMD_PLUGIN_INCLUDE = $(VMD_DIR)/plugins/include
        VMD_PLUGIN_DIR ?= $(VMD_DIR)/plugins/$(VMD_ARCH)/molfile
    endif
else
    VMD_PLUGIN_DIR ?= $(VMD_DIR)/plugins/$(VMD_ARCH)/molfile
endif

CXX      ?= c++
CXXFLAGS ?= -O2 -Wall -Wextra
CXXFLAGS += -std=c++11 -fPIC

ifeq ($(UNAME_S),Darwin)
    LDFLAGS_PLUGIN = -bundle -undefined dynamic_lookup
else
    LDFLAGS_PLUGIN = -shared
endif

PLUGIN  = pdbxplugin.so
DRIVER  = tests/pdbx_dump
ASAN    = tests/pdbx_dump_asan
FUZZER  = tests/fuzz_pdbx
TESTCIF = $(wildcard tests/data/*.cif)

# Coverage-guided fuzzing needs a clang that ships libFuzzer (Apple's does
# not; Homebrew's llvm does). Corpus and crash files go to FUZZ_DIR.
FUZZ_CXX  ?= $(shell for c in /opt/homebrew/opt/llvm/bin/clang++ /usr/local/opt/llvm/bin/clang++; do \
               if [ -x "$$c" ]; then echo "$$c"; exit 0; fi; done; echo clang++)
FUZZ_DIR  ?= $(or $(TMPDIR),/tmp)/pdbx_fuzz
FUZZ_TIME ?= 300
FUZZ_JOBS ?= 4

all: $(PLUGIN)

$(PLUGIN): pdbxplugin.C
	$(CXX) $(CXXFLAGS) -I'$(VMD_PLUGIN_INCLUDE)' -o $@ pdbxplugin.C $(LDFLAGS_PLUGIN)

$(DRIVER): tests/pdbx_dump.cpp pdbxplugin.C
	$(CXX) $(CXXFLAGS) -I'$(VMD_PLUGIN_INCLUDE)' -o $@ tests/pdbx_dump.cpp

$(ASAN): tests/pdbx_dump.cpp pdbxplugin.C
	$(CXX) -std=c++11 -g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined \
	    -I'$(VMD_PLUGIN_INCLUDE)' -o $@ tests/pdbx_dump.cpp

# Each fixture's summary line must match tests/expected.txt; then every
# fixture that loads must survive a write/read round trip unchanged; then the
# error paths (compressed input, read errors, pipes, console output).
test: $(DRIVER)
	@for f in $(TESTCIF); do \
	    printf '%s ' "$$f"; ./$(DRIVER) -q "$$f" 2>/dev/null | grep '^RESULT' || echo "RESULT failed"; \
	done > tests/actual.txt
	@if diff -u tests/expected.txt tests/actual.txt; then \
	    echo "summary lines: all $$(wc -l < tests/expected.txt | tr -d ' ') fixtures match"; \
	else echo "summary lines: MISMATCH (see diff above)"; exit 1; fi
	@for b in 1 2 3 7; do for f in $(TESTCIF); do \
	    printf '%s ' "$$f"; ./$(DRIVER) -q -b $$b "$$f" 2>/dev/null | grep '^RESULT' || echo "RESULT failed"; \
	done > tests/actual.txt; \
	if ! diff -u tests/expected.txt tests/actual.txt; then \
	    echo "summary lines with a $$b-byte read buffer: MISMATCH"; exit 1; fi; done; \
	echo "summary lines: same with 1, 2, 3 and 7-byte read buffers"
	@sh tests/roundtrip.sh ./$(DRIVER) $(TESTCIF)
	@sh tests/messages.sh ./$(DRIVER)
	@rm -f tests/actual.txt

asan: $(ASAN)
	@for f in $(TESTCIF); do \
	    printf '%-40s ' "$$f"; ./$(ASAN) -q "$$f" 2>&1 | grep -E '^RESULT|ERROR' || echo FAILED; \
	done

$(FUZZER): tests/fuzz_pdbx.cpp pdbxplugin.C
	$(FUZZ_CXX) -std=c++11 -g -O1 -fno-omit-frame-pointer \
	    -fsanitize=fuzzer,address,undefined -fno-sanitize-recover=undefined \
	    -I'$(VMD_PLUGIN_INCLUDE)' -o $@ tests/fuzz_pdbx.cpp

# Runs FUZZ_JOBS fuzzers for FUZZ_TIME seconds, seeded with tests/data.
# Inputs that crash, hang (>10 s) or use more than 2 GB are saved in FUZZ_DIR.
fuzz: $(FUZZER)
	@mkdir -p '$(FUZZ_DIR)/corpus' '$(FUZZ_DIR)/crashes'
	./$(FUZZER) -max_total_time=$(FUZZ_TIME) -jobs=$(FUZZ_JOBS) -workers=$(FUZZ_JOBS) \
	    -max_len=16384 -timeout=10 -rss_limit_mb=2048 -dict=tests/fuzz.dict \
	    -artifact_prefix='$(FUZZ_DIR)/crashes/' '$(FUZZ_DIR)/corpus' tests/data

install: $(PLUGIN)
	@echo "Installing $(PLUGIN) into $(VMD_PLUGIN_DIR)"
	mkdir -p '$(VMD_PLUGIN_DIR)'
	cp $(PLUGIN) '$(VMD_PLUGIN_DIR)/'

config:
	@echo "VMD_ARCH:           $(VMD_ARCH)"
	@echo "VMD_DIR:            $(VMD_DIR)"
	@echo "VMD_PLUGIN_INCLUDE: $(VMD_PLUGIN_INCLUDE)"
	@echo "VMD_PLUGIN_DIR:     $(VMD_PLUGIN_DIR)"
	@echo "CXX / CXXFLAGS:     $(CXX) $(CXXFLAGS)"

clean:
	rm -f $(PLUGIN) $(DRIVER) $(ASAN) $(FUZZER) tests/actual.txt fuzz-*.log
	rm -rf tests/*.dSYM

.PHONY: all test asan fuzz install config clean
