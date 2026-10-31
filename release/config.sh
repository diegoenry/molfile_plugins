# Settings shared by the release scripts (sourced, not run).

# Libraries built from source and linked into the plugins. A download is used
# only if its SHA-256 matches.
HDF5_VER=2.2.0
HDF5_URLS="https://github.com/HDFGroup/hdf5/releases/download/$HDF5_VER/hdf5-$HDF5_VER.tar.gz"
HDF5_SHA256=1a1ab8209b35586fbc1aa279ba76d102130b95badcb20ca329587219112d8c16
ZLIB_VER=1.3.2
ZLIB_URLS="https://zlib.net/zlib-$ZLIB_VER.tar.gz https://zlib.net/fossils/zlib-$ZLIB_VER.tar.gz"
ZLIB_SHA256=bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16
# MinGW's thread library ends up in the Windows plugins; its license must ship with them.
WINPTHREADS_LICENSE_URLS="https://raw.githubusercontent.com/mingw-w64/mingw-w64/v14.0.0/mingw-w64-libraries/winpthreads/COPYING"
WINPTHREADS_LICENSE_SHA256=63263614cdd29f2f93cba85e992f041b31f9fc7b4033692f31269489a8a1b177

# Oldest systems the binaries run on
MACOS_ARM64_MIN=11.0
MACOS_X86_64_MIN=10.13
GLIBC_MIN=2.17

ALL_ARCHES="MACOSXARM64 MACOSXX86_64 LINUXAMD64 LINUXARM64 WIN64"

# Work folder: sources, build trees, static libraries, test results. The Linux
# tests mount it into containers, so with colima or Docker Desktop it must be
# under your home folder.
W=${W:-$HOME/.cache/molfile-release}
export W
