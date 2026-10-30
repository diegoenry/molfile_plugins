/*
 * pdbxplugin.C -- PDBx/mmCIF reader and writer plugin for VMD
 *
 * A rewrite of VMD's "pdbx" molfile plugin, and the single mmCIF plugin of
 * this repository: it also replaces the separate "mmcif" reader. It reads
 * text mmCIF and BinaryCIF, plain or gzip-compressed. It registers as
 * "pdbx" (extensions .cif, .mmcif, .bcif) and, for scripts that name the old
 * reader, as "mmcif", both with a higher version number than the plugins they
 * replace, so VMD picks this one when several are loaded. No code is shared
 * with either original implementation.
 *
 * What changed compared with the original plugin (0.14):
 *   - Atoms whose element symbol is not in the periodic table (coarse-grained
 *     beads, dummy atoms, "X", ...) no longer abort the load. They get atomic
 *     number 0 (VMD element "X") and VMD guesses their mass and radius from
 *     the atom name, as it does for any file without element data. The
 *     original symbol stays available in the atom "type" field.
 *   - A real CIF tokenizer instead of fixed line/column assumptions: rows may
 *     span lines, quoted values may contain spaces, semicolon text fields,
 *     comments anywhere, any column order, tag/value (non-loop) categories,
 *     case-insensitive tags. The file is streamed, never held in memory.
 *   - A row with missing or extra values no longer shifts every later row.
 *     Damaged rows are read line by line and matched to their columns by
 *     value type; a row that cannot be placed unambiguously is skipped with
 *     a warning rather than loaded with wrong coordinates.
 *   - Residue numbers and chains come from the author fields (auth_seq_id,
 *     auth_asym_id), matching the PDB-format file; label_* are fallbacks.
 *   - Multi-model files (NMR ensembles) load as trajectory frames.
 *   - Alternate locations (label_alt_id) and insertion codes are kept.
 *   - Explicit bonds from _struct_conn (disulfides, covalent links, metal
 *     coordination: what a PDB file has as CONECT records) in addition to the
 *     original plugin's _pdbx_validate_rmsd_bond/_angle ones.
 *   - Chemical component files (_chem_comp_atom/_chem_comp_bond) load with
 *     their full bond list and bond orders.
 *   - Unit cell from _cell (the 1 x 1 x 1 NMR/EM placeholder is ignored).
 *   - The writer emits valid, quoted, locale-independent mmCIF, one model per
 *     frame.
 *   - gzip-compressed files (.cif.gz, as the PDB distributes them) are read
 *     and written directly when built with zlib (_USE_ZLIB, as VMD's other
 *     plugins do).
 *   - Only some chains can be loaded: "file.cif?chain=A,B" or MMCIF_CHAIN
 *     (from the former mmcif plugin).
 *   - BinaryCIF (.bcif, MessagePack) is read with all seven of its column
 *     encodings and loads exactly like the text file it was made from.
 *
 * Robustness: memory stays proportional to the useful content of the file.
 * A single value is limited to 256 MiB and the number of models to one
 * million; models that do not match the first one keep no coordinates. Read
 * errors and damaged or truncated compressed data fail the load instead of
 * truncating it, other compressed formats are recognized and named, every
 * size and offset in a BinaryCIF file is checked, and text from the file is
 * filtered before it is printed.
 *
 * Reader field mapping (atom_site):
 *   name      <- label_atom_id        (fallback auth_atom_id, then type_symbol)
 *   type      <- type_symbol          (fallback: the atom name)
 *   element   <- type_symbol          (0 / "X" when not a known element)
 *   resname   <- label_comp_id        (fallback auth_comp_id)
 *   resid     <- auth_seq_id          (fallback label_seq_id)
 *   chain     <- first character of auth_asym_id (fallback label_asym_id);
 *                the molfile ABI stores a single character
 *   segname   <- the full chain ID (up to 7 characters)
 *   altloc    <- label_alt_id
 *   insertion <- pdbx_PDB_ins_code
 *   occupancy <- occupancy, beta <- B_iso_or_equiv,
 *   charge    <- pdbx_formal_charge
 *
 * Chain selection: "file.cif?chain=A,B" (the part after the last '?', unless
 * a file with the whole name exists) or the MMCIF_CHAIN environment variable
 * loads only atoms whose chain ID (the one stored in segname: the full ID, or
 * the first 7 characters of a longer one) is listed. IDs are separated by
 * commas or blanks and are case-sensitive. The result is exactly those atoms
 * of the unselected load.
 *
 * Environment variables:
 *   VMDPDBXUSELABELIDS   when set, chain/segname/resid come from label_asym_id
 *                        and label_seq_id instead of the author fields.
 *   MMCIF_CHAIN          chains to load when the file name selects none.
 *
 * Author: Diego Enry Barreto Gomes
 * License: University of Illinois Open Source License
 */

/* vmdplugin.h includes windows.h on Windows: keep its min/max macros out */
#if defined(_WIN32) && !defined(NOMINMAX)
#define NOMINMAX
#endif
#include "molfile_plugin.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <map>
#include <new>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>

#if defined(_WIN32) || defined(_WIN64)
typedef __int64 pdbx_off_t;
#define PDBX_FSEEK _fseeki64
#else
typedef off_t pdbx_off_t;
#define PDBX_FSEEK fseeko
#endif

#ifdef _USE_ZLIB
#include <zlib.h>
#endif

#define PDBX_MAJOR_VERSION 2
#define PDBX_MINOR_VERSION 1

namespace {

/* ========================================================================
 * Small string and number helpers (ASCII only, locale independent)
 * ======================================================================== */

inline bool is_ws(char c) {
  /* CIF whitespace is space, tab and line ends; \v, \f and NUL are accepted
   * as separators too so that damaged files cannot glue tokens together. */
  return c == ' ' || c == '\t' || c == '\n' || c == '\r' ||
         c == '\v' || c == '\f' || c == '\0';
}

inline bool is_eol(char c) { return c == '\n' || c == '\r'; }
inline bool is_digit(char c) { return c >= '0' && c <= '9'; }
inline bool is_alpha(char c) {
  return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}
inline char to_lower(char c) { return (c >= 'A' && c <= 'Z') ? (char)(c + 32) : c; }
inline char to_upper(char c) { return (c >= 'a' && c <= 'z') ? (char)(c - 32) : c; }

std::string lowercase(const char *p, size_t n) {
  std::string s(p, n);
  for (size_t i = 0; i < s.size(); i++) s[i] = to_lower(s[i]);
  return s;
}

/* case-insensitive comparison of a counted string with a C string */
bool iequals(const char *p, size_t n, const char *lit) {
  size_t i = 0;
  for (; i < n && lit[i]; i++)
    if (to_lower(p[i]) != to_lower(lit[i])) return false;
  return i == n && lit[i] == '\0';
}

bool iequals(const std::string &s, const char *lit) {
  return iequals(s.data(), s.size(), lit);
}

bool iprefix(const char *p, size_t n, const char *lit) {
  size_t i = 0;
  for (; lit[i]; i++)
    if (i >= n || to_lower(p[i]) != to_lower(lit[i])) return false;
  return true;
}

const double kPow10[] = {
  1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,  1e10, 1e11,
  1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22
};

/* Parse a CIF number: [+-]digits[.digits][e[+-]digits] with an optional
 * standard uncertainty in parentheses, e.g. "12.345(6)". A comma is accepted
 * as the decimal separator because some writers honor the process locale.
 * Does not depend on LC_NUMERIC. Returns false for anything else. */
bool parse_number(const char *s, size_t n, double *out) {
  size_t i = 0;
  bool neg = false;
  if (i < n && (s[i] == '+' || s[i] == '-')) { neg = (s[i] == '-'); i++; }

  unsigned long long mant = 0;
  int ndig = 0, exp10 = 0;
  bool any = false;
  for (; i < n && is_digit(s[i]); i++) {
    int d = s[i] - '0';
    any = true;
    if (ndig < 19) {
      if (mant || d) { mant = mant * 10 + (unsigned)d; ndig++; }
    } else {
      exp10++;                     /* drop digits beyond double precision */
    }
  }
  if (i < n && (s[i] == '.' || s[i] == ',')) {
    i++;
    for (; i < n && is_digit(s[i]); i++) {
      int d = s[i] - '0';
      any = true;
      if (ndig < 19) {
        if (mant || d) { mant = mant * 10 + (unsigned)d; ndig++; }
        exp10--;
      }
    }
  }
  if (!any) return false;

  if (i < n && (s[i] == 'e' || s[i] == 'E')) {
    size_t j = i + 1;
    bool eneg = false;
    if (j < n && (s[j] == '+' || s[j] == '-')) { eneg = (s[j] == '-'); j++; }
    if (j < n && is_digit(s[j])) {
      int e = 0;
      for (; j < n && is_digit(s[j]); j++)
        if (e < 100000) e = e * 10 + (s[j] - '0');
      exp10 += eneg ? -e : e;
      i = j;
    }
  }

  /* only a standard uncertainty "(...)" may follow */
  if (i < n && s[i] != '(') return false;

  double v = (double)mant;
  if (v != 0.0) {
    if (exp10 > 0)
      v = (exp10 <= 22) ? v * kPow10[exp10] : v * std::pow(10.0, (double)exp10);
    else if (exp10 < 0)
      v = (-exp10 <= 22) ? v / kPow10[-exp10] : v * std::pow(10.0, (double)exp10);
  }
  if (!std::isfinite(v)) return false;
  *out = neg ? -v : v;
  return true;
}

bool parse_float(const std::string &s, float *out) {
  double d;
  if (!parse_number(s.data(), s.size(), &d)) return false;
  if (d > FLT_MAX || d < -FLT_MAX) return false;
  *out = (float)d;
  return true;
}

/* Strict integer parse with saturation: [+-]digits and nothing else. */
bool parse_int(const std::string &s, int *out) {
  size_t i = 0, n = s.size();
  bool neg = false;
  if (i < n && (s[i] == '+' || s[i] == '-')) { neg = (s[i] == '-'); i++; }
  if (i >= n) return false;
  long long v = 0;
  for (; i < n; i++) {
    if (!is_digit(s[i])) return false;
    if (v < 10000000000LL) v = v * 10 + (s[i] - '0');
  }
  if (neg) v = -v;
  /* INT_MIN is reserved as the "no residue number" marker (kSeqMissing) */
  if (v > INT_MAX) v = INT_MAX;
  if (v < INT_MIN + 1) v = INT_MIN + 1;
  *out = (int)v;
  return true;
}

/* Formal charge: "1", "-2", "+1", and also the "2+" / "1-" / "+" spellings. */
bool parse_charge(const std::string &s, float *out) {
  if (parse_float(s, out)) return true;
  size_t n = s.size();
  if (n == 0) return false;
  char sign = s[n - 1];
  if (sign != '+' && sign != '-') return false;
  int mag = 1;
  if (n > 1 && !parse_int(s.substr(0, n - 1), &mag)) return false;
  *out = (float)((sign == '-') ? -mag : mag);
  return true;
}

/* Residue number: an integer, optionally followed by one letter that some
 * writers append as the insertion code ("100A"). */
bool parse_seq(const std::string &s, int *out, char *ins) {
  if (ins) *ins = 0;
  if (parse_int(s, out)) return true;
  size_t n = s.size();
  if (n >= 2 && is_alpha(s[n - 1]) && is_digit(s[n - 2]) && parse_int(s.substr(0, n - 1), out)) {
    if (ins) *ins = s[n - 1];
    return true;
  }
  return false;
}

/* Copy into a fixed molfile string field; returns true if truncated. */
bool copy_field(char *dst, size_t dstsize, const std::string &src) {
  size_t n = src.size();
  bool trunc = false;
  if (n > dstsize - 1) { n = dstsize - 1; trunc = true; }
  memcpy(dst, src.data(), n);
  dst[n] = '\0';
  return trunc;
}

/* File text made safe for the console: printable ASCII only (no terminal
 * escape sequences), at most maxlen characters. */
std::string printable(const std::string &s, size_t maxlen = 40) {
  std::string out;
  for (size_t i = 0; i < s.size(); i++) {
    if (out.size() >= maxlen) { out += "..."; break; }
    unsigned char c = (unsigned char)s[i];
    out += (c >= 32 && c < 127) ? (char)c : '?';
  }
  return out;
}

unsigned fnv1a(unsigned h, const char *p, size_t n) {
  for (size_t i = 0; i < n; i++) { h ^= (unsigned char)p[i]; h *= 16777619u; }
  h ^= 0xffu; h *= 16777619u;               /* field separator */
  return h;
}

/* ========================================================================
 * Element symbols
 *
 * Index = atomic number. The list stops at Rg (111) because that is where
 * VMD's own periodic table ends; heavier symbols are treated as unknown so
 * VMD never sees an atomic number it cannot look up.
 * ======================================================================== */

const char *const kElements[] = {
  "X",
  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
  "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
  "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
  "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
  "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
  "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
  "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",
  "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
  "Mt", "Ds", "Rg"
};
const int kMaxElement = (int)(sizeof(kElements) / sizeof(kElements[0])) - 1;

struct ElementTable {
  unsigned char z[26][27];   /* [first letter][0 = none, 1..26 = second] */
  ElementTable() {
    memset(z, 0, sizeof(z));
    for (int i = 1; i <= kMaxElement; i++) {
      const char *s = kElements[i];
      int c0 = to_upper(s[0]) - 'A';
      int c1 = s[1] ? (to_upper(s[1]) - 'A' + 1) : 0;
      z[c0][c1] = (unsigned char)i;
    }
  }
};

const ElementTable &element_table() {
  static const ElementTable table;   /* thread-safe init (C++11) */
  return table;
}

/* Atomic number for a type_symbol value: > 0 known element, 0 unknown symbol,
 * -1 empty value. Matching is case-insensitive ("FE" and "Fe"), and an ionic
 * charge suffix is accepted ("Fe2+", "O1-", "Na+"). "D" and "T" are hydrogen
 * isotopes. Anything else ("Cg", "Dm", "SC1", "X", "DUM") is unknown. */
int element_from_symbol(const std::string &symbol) {
  /* blanks around the symbol (a quoted 'Fe ') do not matter */
  size_t b = 0, e = symbol.size();
  while (b < e && (symbol[b] == ' ' || symbol[b] == '\t')) b++;
  while (e > b && (symbol[e - 1] == ' ' || symbol[e - 1] == '\t')) e--;
  const char *sym = symbol.data() + b;
  size_t n = e - b;
  if (n == 0) return -1;
  size_t letters = 0;
  while (letters < n && is_alpha(sym[letters])) letters++;
  if (letters == 0 || letters > 2) return 0;
  if (letters < n) {
    /* suffix must be digits followed by exactly one sign */
    if (sym[n - 1] != '+' && sym[n - 1] != '-') return 0;
    for (size_t k = letters; k < n - 1; k++)
      if (!is_digit(sym[k])) return 0;
  }
  int c0 = to_upper(sym[0]) - 'A';
  if (letters == 1 && (c0 == 'D' - 'A' || c0 == 'T' - 'A')) return 1;
  int c1 = (letters == 2) ? (to_upper(sym[1]) - 'A' + 1) : 0;
  return element_table().z[c0][c1];
}

/* ========================================================================
 * Streaming CIF tokenizer
 *
 * Reads the file in chunks; only the token being scanned is kept in memory.
 * Token text points into the internal buffer and stays valid until the next
 * call to next() or restore().
 * ======================================================================== */

enum TokKind { TK_EOF, TK_DATA, TK_LOOP, TK_SAVE, TK_GLOBAL, TK_STOP, TK_TAG, TK_VALUE };

/* Longest single token kept in memory. Real values are at most a few MB
 * (long sequences); only an unterminated text field or binary data gets here. */
const size_t kMaxTokenBytes = (size_t)1 << 28;

/* Initial read buffer; it grows only for a longer value. The test programs
 * define PDBX_TEST_HOOKS to shrink it, so that tokens straddle refills. */
#ifdef PDBX_TEST_HOOKS
size_t g_read_buffer_bytes = (size_t)1 << 16;
#else
const size_t g_read_buffer_bytes = (size_t)1 << 16;
#endif

/* Compressed data can only be read forward: going back to re-read a damaged
 * loop decompresses again from the start of the file. The total of that work
 * is capped at a few times the data read, plus a fixed allowance, so that a
 * file with many damaged loops cannot take quadratic time; past the cap the
 * strict reading of a loop is kept, as for a pipe. */
#ifdef PDBX_TEST_HOOKS
long long g_gz_rewind_factor = 4;
long long g_gz_rewind_extra = (long long)64 << 20;
#elif defined(_USE_ZLIB)
const long long g_gz_rewind_factor = 4;
const long long g_gz_rewind_extra = (long long)64 << 20;
#endif

struct Token {
  TokKind kind;
  const char *p;     /* payload: block name, tag without '_', or value */
  size_t n;
  bool quoted;       /* value was quoted or a text field */
  bool bol;          /* first token on its line */
  long line;
};

inline bool token_is_null(const Token &t) {
  return !t.quoted && t.n == 1 && (t.p[0] == '.' || t.p[0] == '?');
}

struct ReaderState {
  pdbx_off_t offset;
  bool prev_eol;
  bool at_start;
  long line;
  long unterminated_quotes, unterminated_text;
};

class CifReader {
 public:
  CifReader()
    : unterminated_quotes(0), unterminated_text(0), fp_(NULL), pos_(0), end_(0),
      mark_(0), base_(0), eof_(false), prev_eol_(true), at_start_(true), line_(1),
      gz_(false), seekable_(false) {
#ifdef _USE_ZLIB
    memset(&zs_, 0, sizeof(zs_));
    zs_init_ = zin_eof_ = zend_ = false;
    zpos_ = zhigh_ = zrewound_ = 0;
#endif
  }
  ~CifReader() {
#ifdef _USE_ZLIB
    if (zs_init_) inflateEnd(&zs_);
#endif
    if (fp_) fclose(fp_);
  }

  bool open(const char *path) {
    fp_ = fopen(path, "rb");
    if (!fp_) return false;
    seekable_ = (PDBX_FSEEK(fp_, 0, SEEK_CUR) == 0);      /* false for a pipe */
    buf_.resize(g_read_buffer_bytes > 0 ? g_read_buffer_bytes : 1);
    while (end_ < 8 && !eof_) refill();
#ifdef _USE_ZLIB
    if (end_ >= 2 && (unsigned char)buf_[0] == 0x1f && (unsigned char)buf_[1] == 0x8b) {
      start_gzip();
      while (end_ < 8 && !eof_) refill();
    }
#endif
    /* skip a UTF-8 byte order mark */
    if (end_ >= 3 && (unsigned char)buf_[0] == 0xEF &&
        (unsigned char)buf_[1] == 0xBB && (unsigned char)buf_[2] == 0xBF)
      pos_ = mark_ = 3;
    return true;
  }

  /* true if the file is gzip-compressed and being decompressed */
  bool compressed() const { return gz_; }

  /* BinaryCIF: a MessagePack map with string keys ("version", "encoder",
   * "dataBlocks"); text mmCIF never starts with these bytes. */
  bool binarycif() const {
    const unsigned char *b = (const unsigned char *)(buf_.empty() ? "" : &buf_[0]);
    return end_ >= 2 && ((b[0] >= 0x80 && b[0] <= 0x8f) || b[0] == 0xde || b[0] == 0xdf) &&
           ((b[1] >= 0xa0 && b[1] <= 0xbf) || b[1] == 0xd9);
  }

  /* All of the data, for a format that is read whole (BinaryCIF). Call
   * right after open(). Fails beyond cap bytes. */
  void read_all(std::vector<char> &out, size_t cap) {
    out.assign(buf_.begin(), buf_.begin() + (std::ptrdiff_t)end_);
    struct stat st;
    if (!gz_ && fp_ && fstat(fileno(fp_), &st) == 0 && st.st_size > 0 &&
        (unsigned long long)st.st_size <= cap)
      out.reserve((size_t)st.st_size);         /* no regrowth for a plain file */
    const size_t chunk = std::max(g_read_buffer_bytes, (size_t)1);
    while (!eof_) {
      const size_t old = out.size();
      out.resize(old + chunk);
      const size_t got = raw_read(&out[old], chunk);
      out.resize(old + got);
      if (got == 0) eof_ = true;
      if (out.size() > cap)
        throw std::runtime_error("BinaryCIF files larger than 1 GiB are not supported");
    }
  }

  /* Name of a compressed or binary format the data starts with (after
   * decompression), or NULL. Call right after open(). */
  const char *binary_format() const {
    const unsigned char *b = (const unsigned char *)(buf_.empty() ? "" : &buf_[0]);
    size_t n = end_;
    if (n >= 2 && b[0] == 0x1f && b[1] == 0x8b) return "gzip-compressed";
    if (n >= 4 && b[0] == 'B' && b[1] == 'Z' && b[2] == 'h' && b[3] >= '1' && b[3] <= '9')
      return "bzip2-compressed";
    if (n >= 6 && memcmp(b, "\xfd" "7zXZ\0", 6) == 0) return "xz-compressed";
    if (n >= 4 && memcmp(b, "\x28\xb5\x2f\xfd", 4) == 0) return "zstd-compressed";
    if (n >= 4 && memcmp(b, "PK\x03\x04", 4) == 0) return "a zip archive";
    /* BinaryCIF is a MessagePack map with string keys */
    if (n >= 2 && ((b[0] >= 0x80 && b[0] <= 0x8f) || b[0] == 0xde || b[0] == 0xdf) &&
        ((b[1] >= 0xa0 && b[1] <= 0xbf) || b[1] == 0xd9))
      return "BinaryCIF (MessagePack), not text mmCIF";
    return NULL;
  }

  ReaderState state() const {
    ReaderState s;
    s.offset = base_ + (pdbx_off_t)pos_;
    s.prev_eol = prev_eol_;
    s.at_start = at_start_;
    s.line = line_;
    s.unterminated_quotes = unterminated_quotes;
    s.unterminated_text = unterminated_text;
    return s;
  }

  /* Whether restore(s) is sure to succeed `times` times, whatever is read in
   * between. A loop is re-read only if the reader can also come back to
   * its normal reading; otherwise that reading would be lost. */
  bool can_restore(const ReaderState &s, int times) const {
    if (!fp_) return false;
    const bool buffered = s.offset >= base_ && s.offset <= base_ + (pdbx_off_t)end_;
    if (buffered && eof_) return true;      /* nothing is discarded any more */
    if (!seekable_) return false;
#ifdef _USE_ZLIB
    if (gz_) return (long long)zrewound_ + times * (long long)s.offset <= gz_budget();
#else
    (void)times;
#endif
    return true;
  }

  /* Rewind to a saved state (used to re-read a malformed loop). Fails on
   * non-seekable input such as a pipe, unless the state is still buffered. */
  bool restore(const ReaderState &s) {
    if (!fp_) return false;
    if (s.offset >= base_ && s.offset <= base_ + (pdbx_off_t)end_) {
      pos_ = mark_ = (size_t)(s.offset - base_);     /* still in the buffer */
    } else {
#ifdef _USE_ZLIB
      if (gz_) {
        if (!gz_seek(s.offset)) return false;
      } else
#endif
      {
        if (PDBX_FSEEK(fp_, s.offset, SEEK_SET) != 0) return false;
        clearerr(fp_);
      }
      base_ = s.offset;
      pos_ = end_ = mark_ = 0;
      eof_ = false;
    }
    prev_eol_ = s.prev_eol;
    at_start_ = s.at_start;
    line_ = s.line;
    unterminated_quotes = s.unterminated_quotes;   /* don't count twice */
    unterminated_text = s.unterminated_text;
    return true;
  }

  Token next();

  long unterminated_quotes;
  long unterminated_text;

 private:
  FILE *fp_;
  std::vector<char> buf_;
  size_t pos_, end_, mark_;   /* mark_ = start of data that must be kept */
  pdbx_off_t base_;           /* file offset of buf_[0] */
  bool eof_;
  bool prev_eol_;             /* char before pos_ is a line end (or BOF) */
  bool at_start_;             /* no token returned yet */
  long line_;
  bool gz_;                   /* decompressing gzip data */
  bool seekable_;             /* the file can be repositioned (not a pipe) */

#ifdef _USE_ZLIB
  z_stream zs_;
  bool zs_init_;
  std::vector<unsigned char> zin_;   /* compressed input */
  bool zin_eof_;              /* no more compressed input */
  bool zend_;                 /* last gzip member finished */
  pdbx_off_t zpos_;           /* offset of the next decompressed byte */
  pdbx_off_t zhigh_;          /* farthest offset decompressed so far */
  pdbx_off_t zrewound_;       /* bytes decompressed again by gz_seek */

  /* Switch to decompression; the bytes read so far are compressed input. */
  void start_gzip() {
    zin_.assign(std::max(std::max(g_read_buffer_bytes, (size_t)16), end_), 0);
    memcpy(&zin_[0], &buf_[0], end_);
    if (inflateInit2(&zs_, 15 + 16) != Z_OK)        /* gzip wrapper only */
      throw std::runtime_error("cannot initialize zlib");
    zs_init_ = true;
    zs_.next_in = &zin_[0];
    zs_.avail_in = (uInt)end_;
    zin_eof_ = eof_;
    zend_ = false;
    zpos_ = zhigh_ = zrewound_ = 0;
    gz_ = true;
    pos_ = end_ = mark_ = 0;
    base_ = 0;
    eof_ = false;
  }

  /* More compressed input after what is left; false at the end of file. */
  bool zfill() {
    if (zin_eof_) return false;
    size_t keep = zs_.avail_in;
    if (keep >= zin_.size()) return true;
    if (keep > 0 && zs_.next_in != &zin_[0]) memmove(&zin_[0], zs_.next_in, keep);
    size_t want = zin_.size() - keep;
    size_t got = fread(&zin_[keep], 1, want, fp_);
    if (got < want && ferror(fp_))
      throw std::runtime_error(std::string("read error: ") + strerror(errno));
    if (got == 0) zin_eof_ = true;
    zs_.next_in = &zin_[0];
    zs_.avail_in = (uInt)(keep + got);
    return got > 0;
  }

  /* Decompress up to want bytes; 0 only at the end of the data. */
  size_t gz_read(char *dst, size_t want) {
    if (zend_ || want == 0) return 0;
    if (want > ((size_t)1 << 30)) want = (size_t)1 << 30;
    zs_.next_out = (Bytef *)dst;
    zs_.avail_out = (uInt)want;
    while (zs_.avail_out == (uInt)want) {          /* until some output */
      if (zs_.avail_in == 0 && !zfill())
        throw std::runtime_error("the compressed data ends early (incomplete download "
                                 "or truncated .gz file)");
      int r = inflate(&zs_, Z_NO_FLUSH);
      if (r == Z_STREAM_END) {
        /* another gzip member may follow (concatenated .gz files) */
        if (zs_.avail_in < 2) zfill();
        if (zs_.avail_in >= 2 && zs_.next_in[0] == 0x1f && zs_.next_in[1] == 0x8b) {
          if (inflateReset(&zs_) != Z_OK) throw std::runtime_error("zlib error");
          continue;
        }
        /* After the last member only zero padding may follow. Anything else
         * is most likely a damaged member header, and ignoring it would load
         * part of the file as if it were all of it. */
        for (;;) {
          for (uInt i = 0; i < zs_.avail_in; i++)
            if (zs_.next_in[i] != 0)
              throw std::runtime_error("the compressed data is damaged (unexpected bytes "
                                       "after a gzip member)");
          zs_.avail_in = 0;
          if (!zfill()) break;
        }
        zend_ = true;
        break;
      }
      if (r == Z_BUF_ERROR && zs_.avail_in == 0) continue;   /* needs more input */
      if (r != Z_OK)
        throw std::runtime_error(std::string("the compressed data is damaged (") +
                                 (zs_.msg ? zs_.msg : "zlib error") + ")");
    }
    size_t got = want - zs_.avail_out;
    zpos_ += (pdbx_off_t)got;
    if (zpos_ > zhigh_) zhigh_ = zpos_;
    return got;
  }

  /* how many bytes may be decompressed again, in total */
  long long gz_budget() const {
    return g_gz_rewind_factor * (long long)zhigh_ + g_gz_rewind_extra;
  }

  /* Continue decompressing at offset off of the decompressed data. */
  bool gz_seek(pdbx_off_t off) {
    if (off < zpos_) {
      if ((long long)zrewound_ + (long long)off > gz_budget()) return false;
      if (PDBX_FSEEK(fp_, 0, SEEK_SET) != 0) return false;   /* e.g. a pipe */
      clearerr(fp_);
      if (inflateReset(&zs_) != Z_OK) return false;
      zs_.next_in = &zin_[0];
      zs_.avail_in = 0;
      zin_eof_ = zend_ = false;
      zpos_ = 0;
      zrewound_ += off;
    }
    std::vector<char> skip((size_t)1 << 16);
    while (zpos_ < off) {
      size_t n = (size_t)std::min((pdbx_off_t)skip.size(), off - zpos_);
      if (gz_read(&skip[0], n) == 0) return false;
    }
    return true;
  }
#endif

  size_t raw_read(char *dst, size_t want) {
#ifdef _USE_ZLIB
    if (gz_) return gz_read(dst, want);
#endif
    size_t got = fread(dst, 1, want, fp_);
    if (got < want && ferror(fp_))   /* never mistake an I/O error for the end */
      throw std::runtime_error(std::string("read error: ") + strerror(errno));
    return got;
  }

  /* Keep [mark_, end_), read more data. True if pos_ < end_ afterwards. */
  bool refill() {
    if (eof_ || !fp_) return pos_ < end_;
    if (mark_ > 0) {
      size_t keep = end_ - mark_;
      if (keep) memmove(&buf_[0], &buf_[mark_], keep);
      pos_ -= mark_;
      end_ = keep;
      base_ += (pdbx_off_t)mark_;
      mark_ = 0;
    }
    if (end_ == buf_.size()) {        /* one token larger than the buffer */
      if (buf_.size() >= kMaxTokenBytes)
        throw std::runtime_error("a single value is larger than 256 MiB (an unterminated "
                                 "text field?); the file is damaged");
      buf_.resize(buf_.size() * 2);
    }
    size_t want = buf_.size() - end_;
    size_t got = raw_read(&buf_[end_], want);
    if (got == 0) {
      eof_ = true;
      return pos_ < end_;
    }
    end_ += got;
    return pos_ < end_;
  }

  bool avail() { return pos_ < end_ || refill(); }

  /* character after pos_, or -1 at end of input */
  int peek1() {
    if (pos_ + 1 >= end_) refill();
    return (pos_ + 1 < end_) ? (unsigned char)buf_[pos_ + 1] : -1;
  }

  void adv() {
    char c = buf_[pos_++];
    if (c == '\n') line_++;
    prev_eol_ = is_eol(c);
  }

  Token make(TokKind k, size_t from, size_t to, bool quoted, bool bol, long line) {
    Token t;
    t.kind = k;
    t.p = buf_.empty() ? "" : &buf_[0] + from;
    t.n = to - from;
    t.quoted = quoted;
    t.bol = bol;
    t.line = line;
    return t;
  }
};

Token CifReader::next() {
  bool nl = at_start_;

  /* whitespace and comments */
  for (;;) {
    mark_ = pos_;
    if (!avail()) {
      Token t = make(TK_EOF, pos_, pos_, false, nl, line_);
      return t;
    }
    char c = buf_[pos_];
    if (c == '\n') { pos_++; line_++; prev_eol_ = true; nl = true; continue; }
    if (c == '\r') { pos_++; prev_eol_ = true; nl = true; continue; }
    if (is_ws(c)) { pos_++; prev_eol_ = false; continue; }
    if (c == '#') {
      for (;;) {
        while (pos_ < end_ && !is_eol(buf_[pos_])) pos_++;
        if (pos_ < end_) break;
        mark_ = pos_;
        if (!refill()) break;
      }
      prev_eol_ = false;
      continue;
    }
    break;
  }

  at_start_ = false;
  long line = line_;
  char c = buf_[pos_];

  /* semicolon text field: ';' in the first column */
  if (c == ';' && prev_eol_) {
    adv();
    for (;;) {
      if (!avail()) {
        unterminated_text++;
        return make(TK_VALUE, mark_ + 1, pos_, true, nl, line);
      }
      if (buf_[pos_] == ';' && prev_eol_) {
        size_t vend = pos_;
        if (vend > mark_ + 1 && buf_[vend - 1] == '\n') vend--;
        if (vend > mark_ + 1 && buf_[vend - 1] == '\r') vend--;
        adv();
        return make(TK_VALUE, mark_ + 1, vend, true, nl, line);
      }
      adv();
    }
  }

  /* quoted string: ends at a matching quote followed by whitespace; a quote
   * cannot span lines, so an unterminated one ends at the line end */
  if (c == '\'' || c == '"') {
    char q = c;
    adv();
    for (;;) {
      if (!avail()) {
        unterminated_quotes++;
        return make(TK_VALUE, mark_ + 1, pos_, true, nl, line);
      }
      char d = buf_[pos_];
      if (is_eol(d)) {
        unterminated_quotes++;
        return make(TK_VALUE, mark_ + 1, pos_, true, nl, line);
      }
      if (d == q) {
        int nx = peek1();
        if (nx < 0 || is_ws((char)nx)) {
          size_t vend = pos_;
          adv();
          return make(TK_VALUE, mark_ + 1, vend, true, nl, line);
        }
      }
      adv();
    }
  }

  /* bare word */
  for (;;) {
    while (pos_ < end_ && !is_ws(buf_[pos_])) pos_++;
    if (pos_ < end_ || !refill()) break;
  }
  prev_eol_ = false;
  const char *w = &buf_[mark_];
  size_t n = pos_ - mark_;

  if (w[0] == '_')
    return make(TK_TAG, mark_ + 1, pos_, false, nl, line);
  if (n >= 5 && iprefix(w, n, "data_"))
    return make(TK_DATA, mark_ + 5, pos_, false, nl, line);
  if (iequals(w, n, "loop_"))
    return make(TK_LOOP, pos_, pos_, false, nl, line);
  if (n >= 5 && iprefix(w, n, "save_"))
    return make(TK_SAVE, mark_ + 5, pos_, false, nl, line);
  if (iequals(w, n, "global_"))
    return make(TK_GLOBAL, pos_, pos_, false, nl, line);
  if (iequals(w, n, "stop_"))
    return make(TK_STOP, pos_, pos_, false, nl, line);
  return make(TK_VALUE, mark_, pos_, false, nl, line);
}

/* ========================================================================
 * Loop consumers ("sinks")
 * ======================================================================== */

class Sink {
 public:
  virtual ~Sink() {}
  /* item names are lowercase, without the category prefix */
  virtual void begin(const std::vector<std::string> &items) = 0;
  /* col can reach past the last item for a row with extra values */
  virtual void value(size_t col, const char *p, size_t n, bool null) = 0;
  /* nvals differs from the item count when the row was damaged */
  virtual void end_row(size_t nvals) = 0;
  virtual void reset() = 0;
  /* whether a tag/value item is used at all (others are not stored) */
  virtual bool wants(const std::string &item) const { (void)item; return true; }
  /* on while a damaged loop is re-read line by line */
  virtual void set_recovery(bool on) { (void)on; }
  /* rows that look shifted or could not be used */
  virtual long bad_rows() const { return 0; }
  virtual size_t nrows() const = 0;
};

/* The items of a category that a consumer needs. Values are stored back to
 * back in one string (nulls as ""), so memory stays close to the file size. */
struct Table {
  std::vector<std::string> items;
  std::string text;
  std::vector<size_t> ends;     /* end of each value in text, row by row */
  size_t nrows;
  Table() : nrows(0) {}
  void clear_rows() { text.clear(); ends.clear(); nrows = 0; }
  int col(const char *name) const {
    for (size_t i = 0; i < items.size(); i++)
      if (items[i] == name) return (int)i;
    return -1;
  }
  std::string get(size_t row, int c) const {
    if (c < 0 || row >= nrows || (size_t)c >= items.size()) return std::string();
    size_t i = row * items.size() + (size_t)c;
    size_t b = i ? ends[i - 1] : 0;
    return text.substr(b, ends[i] - b);
  }
};

class TableSink : public Sink {
 public:
  Table table;
  /* wanted: NULL-terminated item names; max_rows: rows kept (e.g. 1 for _cell) */
  explicit TableSink(const char *const *wanted, size_t max_rows = (size_t)-1)
    : wanted_(wanted), max_rows_(max_rows) {}
  bool wants(const std::string &item) const {
    for (const char *const *w = wanted_; *w; w++)
      if (item == *w) return true;
    return false;
  }
  void begin(const std::vector<std::string> &items) {
    std::vector<std::string> kept;
    sel_.assign(items.size(), -1);
    for (size_t i = 0; i < items.size(); i++)
      if (wants(items[i])) { sel_[i] = (int)kept.size(); kept.push_back(items[i]); }
    if (kept != table.items) { table.items = kept; table.clear_rows(); }
    row_.assign(kept.size(), std::string());
  }
  void value(size_t col, const char *p, size_t n, bool null) {
    if (col < sel_.size() && sel_[col] >= 0) {
      std::string &v = row_[(size_t)sel_[col]];
      if (null) v.clear(); else v.assign(p, n);
    }
  }
  void end_row(size_t) {
    if (!row_.empty() && table.nrows < max_rows_) {
      for (size_t k = 0; k < row_.size(); k++) {
        table.text += row_[k];
        table.ends.push_back(table.text.size());
      }
      table.nrows++;
    }
    for (size_t k = 0; k < row_.size(); k++) row_[k].clear();
  }
  void reset() { table.clear_rows(); }
  size_t nrows() const { return table.nrows; }
 private:
  const char *const *wanted_;
  size_t max_rows_;
  std::vector<int> sel_;
  std::vector<std::string> row_;
};

/* An array that grows in fixed-size chunks. Unlike a vector it never
 * reallocates and copies its contents as it grows, so millions of atoms need
 * no transient second copy, and no large freed blocks are left behind for
 * the allocator to keep. */
template <class T> class ChunkedArray {
 public:
  ChunkedArray() : n_(0) {}
  void push_back(const T &v) {
    if (chunks_.empty() || chunks_.back().size() == kChunk) {
      chunks_.push_back(std::vector<T>());
      if (chunks_.size() > 1) chunks_.back().reserve(kChunk);   /* small files stay small */
    }
    chunks_.back().push_back(v);
    n_++;
  }
  T &operator[](size_t i) { return chunks_[i / kChunk][i % kChunk]; }
  const T &operator[](size_t i) const { return chunks_[i / kChunk][i % kChunk]; }
  size_t size() const { return n_; }
  bool empty() const { return n_ == 0; }
  void clear() { std::vector<std::vector<T> >().swap(chunks_); n_ = 0; }   /* frees */
  void swap(ChunkedArray &o) { chunks_.swap(o.chunks_); std::swap(n_, o.n_); }
  void copy_to(T *dst) const {
    for (size_t c = 0; c < chunks_.size(); c++) {
      if (!chunks_[c].empty()) memcpy(dst, &chunks_[c][0], chunks_[c].size() * sizeof(T));
      dst += chunks_[c].size();
    }
  }
 private:
  static const size_t kChunk = 65536;
  std::vector<std::vector<T> > chunks_;
  size_t n_;
};

/* String interning for chain IDs (few distinct values, many atoms). */
class StringPool {
 public:
  StringPool() : last_(-1) {}
  int intern(const std::string &s) {
    if (last_ >= 0 && strs_[(size_t)last_] == s) return last_;
    std::unordered_map<std::string, int>::const_iterator it = map_.find(s);
    if (it != map_.end()) return last_ = it->second;
    int id = (int)strs_.size();
    strs_.push_back(s);
    map_.insert(std::make_pair(s, id));
    return last_ = id;
  }
  int find(const std::string &s) const {
    std::unordered_map<std::string, int>::const_iterator it = map_.find(s);
    return (it == map_.end()) ? -1 : it->second;
  }
  const std::string &str(int id) const { return strs_[(size_t)id]; }
  size_t size() const { return strs_.size(); }
  void clear() { map_.clear(); strs_.clear(); last_ = -1; }
 private:
  std::unordered_map<std::string, int> map_;
  std::vector<std::string> strs_;
  int last_;
};

enum AtomCol {
  AC_GROUP, AC_TYPE, AC_LATOM, AC_AATOM, AC_ALT, AC_LCOMP, AC_ACOMP,
  AC_LASYM, AC_AASYM, AC_LSEQ, AC_ASEQ, AC_INS, AC_X, AC_Y, AC_Z,
  AC_OCC, AC_B, AC_CHARGE, AC_MODEL, AC_MODEL2,
  AC_ID, AC_ENTITY,          /* only used to check damaged rows */
  AC_NUM
};

/* columns holding names or codes (missing values there are usually empty
 * strings a writer left out) */
inline bool is_text_col(int k) {
  return k == AC_TYPE || k == AC_LATOM || k == AC_AATOM || k == AC_ALT ||
         k == AC_LCOMP || k == AC_ACOMP || k == AC_LASYM || k == AC_AASYM || k == AC_INS;
}

struct ColAlias { const char *item; int col; };

const ColAlias kAtomSiteCols[] = {
  {"group_pdb", AC_GROUP},          {"type_symbol", AC_TYPE},
  {"label_atom_id", AC_LATOM},      {"auth_atom_id", AC_AATOM},
  {"label_alt_id", AC_ALT},         {"label_comp_id", AC_LCOMP},
  {"auth_comp_id", AC_ACOMP},       {"label_asym_id", AC_LASYM},
  {"auth_asym_id", AC_AASYM},       {"label_seq_id", AC_LSEQ},
  {"auth_seq_id", AC_ASEQ},         {"pdbx_pdb_ins_code", AC_INS},
  {"cartn_x", AC_X},                {"cartn_y", AC_Y},
  {"cartn_z", AC_Z},                {"occupancy", AC_OCC},
  {"b_iso_or_equiv", AC_B},         {"pdbx_formal_charge", AC_CHARGE},
  {"pdbx_pdb_model_num", AC_MODEL}, {"ihm_model_id", AC_MODEL2},
  {"id", AC_ID},                    {"label_entity_id", AC_ENTITY},
  {NULL, 0}
};

/* IHM starting-model coordinates use shorter item names. Different
 * starting_model_id values are templates for different parts of the system,
 * not conformations of the same atoms, so they are all loaded together (as
 * the original plugin did) rather than split into frames. */
const ColAlias kIhmStartCols[] = {
  {"group_pdb", AC_GROUP},   {"type_symbol", AC_TYPE},
  {"atom_id", AC_LATOM},     {"comp_id", AC_LCOMP},
  {"asym_id", AC_LASYM},     {"seq_id", AC_LSEQ},
  {"cartn_x", AC_X},         {"cartn_y", AC_Y},
  {"cartn_z", AC_Z},         {"b_iso_or_equiv", AC_B},
  {"id", AC_ID},             {"entity_id", AC_ENTITY},
  {NULL, 0}
};

const int kSeqMissing = INT_MIN;

/* More distinct model ids than any real ensemble or trajectory has; beyond
 * this, rows of new models are ignored so bookkeeping stays bounded. */
const size_t kMaxModels = 1000000;

/* Chains to load ("file.cif?chain=A,B" or MMCIF_CHAIN). Matched against the
 * chain ID stored in segname (auth_asym_id, or label_asym_id with
 * VMDPDBXUSELABELIDS): the full ID, or for a longer ID its first 7
 * characters as segname shows them. Case-sensitive, as chain IDs are. */
struct ChainFilter {
  std::vector<std::string> ids;     /* sorted, unique */
  std::vector<char> hit;            /* ids[i] occurs in the file */
  std::string source;               /* where the selection came from */

  ChainFilter() : have_last_(false), last_k_(-1) {}
  bool active() const { return !ids.empty(); }

  /* add IDs separated by commas or blanks */
  void add(const std::string &list) {
    size_t i = 0;
    while (i < list.size()) {
      while (i < list.size() && (list[i] == ',' || is_ws(list[i]))) i++;
      size_t b = i;
      while (i < list.size() && list[i] != ',' && !is_ws(list[i])) i++;
      if (i > b) ids.push_back(list.substr(b, i - b));
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    hit.assign(ids.size(), 0);
    have_last_ = false;
  }

  /* index of chain in ids, or -1. A chain ID longer than segname holds is
   * also selected by its first characters, which is what VMD shows. */
  int find(const std::string &chain) {
    if (have_last_ && chain == last_) return last_k_;       /* rows come in runs */
    const size_t seglen = sizeof(molfile_atom_t::segid) - 1;
    last_k_ = index_of(chain);
    if (last_k_ < 0 && chain.size() > seglen) last_k_ = index_of(chain.substr(0, seglen));
    last_ = chain;
    have_last_ = true;
    return last_k_;
  }

  /* selected? (and remember that it occurs) */
  bool accepts(const std::string &chain) {
    int k = find(chain);
    if (k >= 0) hit[(size_t)k] = 1;
    return k >= 0;
  }

  /* "A, B" -- all IDs, or only those (not) found in the file */
  std::string list(int which = 0) const {
    std::string s;
    size_t shown = 0;
    for (size_t i = 0; i < ids.size(); i++) {
      if ((which > 0 && !hit[i]) || (which < 0 && hit[i])) continue;
      if (shown++ == 20) { s += ", ..."; break; }
      if (!s.empty()) s += ", ";
      s += printable(ids[i], 20);
    }
    return s;
  }

 private:
  bool have_last_;
  int last_k_;
  std::string last_;

  int index_of(const std::string &id) const {
    std::vector<std::string>::const_iterator it =
        std::lower_bound(ids.begin(), ids.end(), id);
    return (it != ids.end() && *it == id) ? (int)(it - ids.begin()) : -1;
  }
};

/* A later model. Its rows are compared with the first model as they arrive;
 * once one differs, its coordinates are dropped, since it cannot be a frame. */
struct ExtraModel {
  std::string key;
  ChunkedArray<float> xyz;
  size_t rows;
  bool matches;
  ExtraModel() : rows(0), matches(true) {}
};

/* Streams _atom_site rows into compact per-atom arrays. The first model
 * defines the structure; each further model only stores coordinates. */
class AtomSink : public Sink {
 public:
  explicit AtomSink(const ColAlias *aliases, bool use_label)
    : aliases_(aliases), use_label_(use_label), filter_(NULL), ncols_(0) { reset(); }

  /* load only rows of these chains (NULL: all) */
  void set_filter(ChainFilter *f) {
    filter_ = (f && f->active()) ? f : NULL;
    found.assign(filter_ ? filter_->ids.size() : 0, 0);
  }

  /* model 1 (structure) */
  ChunkedArray<molfile_atom_t> atoms;
  ChunkedArray<float> xyz;
  ChunkedArray<unsigned> ident;
  /* identifiers kept for bond lookup */
  StringPool asyms;
  ChunkedArray<int> auth_asym, label_asym;  /* pool ids, -1 if absent */
  ChunkedArray<int> auth_seq, label_seq;    /* kSeqMissing if absent */
  ChunkedArray<char> insertion;             /* author insertion code, 0 = none */
  std::string model0_key;
  std::vector<ExtraModel> extra;

  /* which columns were seen */
  bool has_col[AC_NUM];

  /* statistics for diagnostics */
  long bad_coord_rows, bad_group_rows;
  long bad_field_rows;             /* a value that cannot be right for its column */
  long fixed_rows, dropped_rows;   /* damaged rows: placed / skipped */
  long n_known_elem, n_unknown_elem, n_no_elem;
  long multichar_chain, long_chain, long_altloc, long_name, long_resname;
  long missing_resid;
  long dropped_model_rows;      /* rows of models beyond kMaxModels */
  long filtered_rows;           /* rows of chains not selected */
  std::vector<char> found;      /* with a filter: ids[i] has atoms in the structure */
  std::vector<std::string> chains_seen;   /* with a filter: chains of the first model */
  bool more_chains;             /* chains_seen is incomplete */
  bool any_altloc, any_insertion, any_charge;
  std::map<std::string, long> unknown_symbols;

  void reset() {
    atoms.clear(); xyz.clear(); ident.clear(); asyms.clear();
    auth_asym.clear(); label_asym.clear(); auth_seq.clear(); label_seq.clear();
    insertion.clear();
    model0_key.clear(); extra.clear(); model_index_.clear();
    cur_model_ = -1;
    for (int i = 0; i < AC_NUM; i++) { has_col[i] = false; have_[i] = false; }
    bad_coord_rows = bad_group_rows = bad_field_rows = fixed_rows = dropped_rows = 0;
    n_known_elem = n_unknown_elem = n_no_elem = 0;
    multichar_chain = long_chain = long_altloc = long_name = long_resname = 0;
    missing_resid = 0;
    dropped_model_rows = 0;
    filtered_rows = 0;
    found.assign(found.size(), 0);
    chains_seen.clear();
    more_chains = false;
    any_altloc = any_insertion = any_charge = false;
    unknown_symbols.clear();
    recovering_ = false;
    prof_.clear();
    prof_rows_ = 0;
  }

  bool wants(const std::string &item) const {
    for (const ColAlias *a = aliases_; a->item; a++)
      if (item == a->item) return true;
    return false;
  }

  void begin(const std::vector<std::string> &items) {
    ncols_ = items.size();
    colmap_.assign(items.size(), -1);
    for (int k = 0; k < AC_NUM; k++) srccol_[k] = -1;
    for (size_t i = 0; i < items.size(); i++) {
      for (const ColAlias *a = aliases_; a->item; a++) {
        if (items[i] == a->item) {
          colmap_[i] = a->col;
          has_col[a->col] = true;
          srccol_[a->col] = (int)i;          /* a repeated item: the last wins */
          break;
        }
      }
    }
    if (raw_.size() < ncols_) { raw_.resize(ncols_); rawnull_.resize(ncols_); }
    for (int k = 0; k < AC_NUM; k++) have_[k] = false;
  }

  void value(size_t col, const char *p, size_t n, bool null) {
    if (col >= 2 * ncols_) return;          /* far too many values: ignored */
    if (col >= raw_.size()) { raw_.resize(col + 1); rawnull_.resize(col + 1); }
    rawnull_[col] = null ? 1 : 0;
    if (!null) raw_[col].assign(p, n);
  }

  void set_recovery(bool on) {
    recovering_ = on;
    prof_.assign(on ? ncols_ : 0, ColProfile());
    prof_rows_ = 0;
  }

  void end_row(size_t nvals) {
    if (nvals == ncols_) {
      if (recovering_) learn();
      for (int k = 0; k < AC_NUM; k++) {
        int c = srccol_[k];
        src_[k] = (c >= 0 && !rawnull_[(size_t)c]) ? c : -1;
      }
    } else if (align(nvals)) {
      fixed_rows++;
    } else {
      dropped_rows++;
      return;
    }
    load_fields();
    process_row();
  }

  /* used to choose between the strict and the line-by-line reading of a
   * damaged loop */
  long bad_rows() const {
    return bad_coord_rows + bad_group_rows + bad_field_rows + dropped_rows;
  }
  /* rows with coordinates, including those of chains not selected: the
   * first data block with coordinates is the one read, whatever the filter */
  size_t nrows() const { return atoms.size() + (size_t)filtered_rows; }

 private:
  struct ColProfile {
    long nnull, nnum, ndot;
    ColProfile() : nnull(0), nnum(0), ndot(0) {}
  };

  /* How the values of a damaged row map onto the columns: a gap of `gap`
   * missing values at position `at` plus `tail` missing at the end, or
   * `extra` surplus values dropped at position `at`. */
  struct Placement {
    size_t at, gap, tail, extra, ncols;
    long raw_index(size_t c) const {          /* -1: a missing value */
      if (extra) return (long)(c < at ? c : c + extra);
      if (c < at) return (long)c;
      if (c < at + gap) return -1;
      if (c < ncols - tail) return (long)(c - gap);
      return -1;
    }
    bool end_padding() const { return !extra && gap == 0; }
  };

  const ColAlias *aliases_;
  bool use_label_;
  ChainFilter *filter_;
  size_t ncols_;
  std::vector<int> colmap_;                /* column -> AtomCol, or -1 */
  int srccol_[AC_NUM];                     /* AtomCol -> column, or -1 */
  std::vector<std::string> raw_;           /* values of the current row */
  std::vector<char> rawnull_;
  int src_[AC_NUM];                        /* AtomCol -> raw value, -1 = null */
  std::string val_[AC_NUM];
  bool have_[AC_NUM];
  int cur_model_;                                   /* -1 = none yet */
  std::unordered_map<std::string, int> model_index_;
  bool recovering_;
  std::vector<ColProfile> prof_;           /* complete rows of this loop */
  long prof_rows_;

  void load_fields() {
    for (int k = 0; k < AC_NUM; k++) {
      if (src_[k] < 0) { have_[k] = false; continue; }
      val_[k].swap(raw_[(size_t)src_[k]]);   /* raw_ is rewritten next row */
      have_[k] = true;
    }
    if (have_[AC_TYPE]) {                     /* a quoted 'Fe ' is the symbol Fe */
      std::string &t = val_[AC_TYPE];
      size_t b = 0, e = t.size();
      while (b < e && (t[b] == ' ' || t[b] == '\t')) b++;
      while (e > b && (t[e - 1] == ' ' || t[e - 1] == '\t')) e--;
      if (b > 0 || e < t.size()) t = t.substr(b, e - b);
      if (t.empty()) have_[AC_TYPE] = false;
    }
  }

  static bool is_number(const std::string &v) {
    double d;
    return parse_number(v.data(), v.size(), &d);
  }

  void learn() {
    prof_rows_++;
    for (size_t c = 0; c < ncols_; c++) {
      ColProfile &p = prof_[c];
      if (rawnull_[c]) { p.nnull++; continue; }
      if (is_number(raw_[c])) {
        p.nnum++;
        if (raw_[c].find_first_of(".,") != std::string::npos) p.ndot++;
      }
    }
  }

  const std::string *current_model_key() const {
    if (cur_model_ < 0) return NULL;
    return (cur_model_ == 0) ? &model0_key : &extra[(size_t)cur_model_ - 1].key;
  }

  /* Score of a placement (higher is more plausible), or -1 if some value
   * does not fit its column. Fills src with the raw value of each AtomCol. */
  int evaluate(const Placement &pl, int *src) const {
    static const int coords[3] = { AC_X, AC_Y, AC_Z };
    for (int i = 0; i < 3; i++) {             /* cheap test first */
      int c = srccol_[coords[i]];
      if (c < 0) return -1;
      long r = pl.raw_index((size_t)c);
      if (r < 0 || rawnull_[(size_t)r] || !is_number(raw_[(size_t)r])) return -1;
      if (prof_rows_ > 0 && prof_[(size_t)c].ndot == prof_rows_ &&
          raw_[(size_t)r].find_first_of(".,") == std::string::npos)
        return -1;                            /* "2" where "2.000" belongs */
    }
    const long kMinRowsForNulls = 2;
    int score = 0;
    for (int k = 0; k < AC_NUM; k++) src[k] = -1;
    for (size_t c = 0; c < ncols_; c++) {
      long r = pl.raw_index(c);
      int kc = colmap_[c];
      if (r < 0) {
        if (kc >= 0) { src[kc] = -1; if (is_text_col(kc)) score++; }
        continue;
      }
      bool null = rawnull_[(size_t)r] != 0;
      if (prof_rows_ > 0) {                   /* compare with complete rows */
        const ColProfile &p = prof_[c];
        long nonnull = prof_rows_ - p.nnull;
        if (null) {
          if (prof_rows_ >= kMinRowsForNulls && p.nnull == 0) return -1;
        } else {
          if (prof_rows_ >= kMinRowsForNulls && nonnull == 0) return -1;
          bool num = is_number(raw_[(size_t)r]);
          if (nonnull > 0 && p.nnum == nonnull && !num) return -1;
          if (nonnull > 0 && p.nnum == 0 && num) return -1;
        }
      }
      if (kc >= 0) src[kc] = null ? -1 : (int)r;
    }

    /* values that cannot be right for their column */
    if (src[AC_GROUP] >= 0 && !iequals(raw_[(size_t)src[AC_GROUP]], "ATOM") &&
        !iequals(raw_[(size_t)src[AC_GROUP]], "HETATM"))
      return -1;
    static const int numeric[] = { AC_ID, AC_ENTITY, AC_OCC, AC_B, AC_MODEL, AC_MODEL2 };
    for (size_t i = 0; i < sizeof(numeric) / sizeof(numeric[0]); i++)
      if (src[numeric[i]] >= 0 && !is_number(raw_[(size_t)src[numeric[i]]])) return -1;
    int seq;
    float q;
    if (src[AC_LSEQ] >= 0 && !parse_seq(raw_[(size_t)src[AC_LSEQ]], &seq, NULL)) return -1;
    if (src[AC_ASEQ] >= 0 && !parse_seq(raw_[(size_t)src[AC_ASEQ]], &seq, NULL)) return -1;
    if (src[AC_CHARGE] >= 0 && !parse_charge(raw_[(size_t)src[AC_CHARGE]], &q)) return -1;
    if (src[AC_INS] >= 0 && raw_[(size_t)src[AC_INS]].size() != 1) return -1;
    if (src[AC_TYPE] >= 0) {
      const std::string &t = raw_[(size_t)src[AC_TYPE]];
      size_t b = t.find_first_not_of(" \t");
      if (b == std::string::npos || !is_alpha(t[b])) return -1;
    }

    /* evidence that the values sit in the right columns */
    if (src[AC_LCOMP] >= 0 && src[AC_ACOMP] >= 0 &&
        raw_[(size_t)src[AC_LCOMP]] == raw_[(size_t)src[AC_ACOMP]])
      score += 2;
    if (src[AC_LATOM] >= 0 && src[AC_AATOM] >= 0 &&
        raw_[(size_t)src[AC_LATOM]] == raw_[(size_t)src[AC_AATOM]])
      score += 2;
    if (src[AC_LASYM] >= 0 && src[AC_AASYM] >= 0 &&
        raw_[(size_t)src[AC_LASYM]] == raw_[(size_t)src[AC_AASYM]])
      score++;                                /* the usual case for polymers */
    int nm = (src[AC_LATOM] >= 0) ? src[AC_LATOM] : src[AC_AATOM];
    if (src[AC_TYPE] >= 0 && nm >= 0) {       /* "C" and "CA", "Fe" and "FE" */
      const std::string &t = raw_[(size_t)src[AC_TYPE]], &a = raw_[(size_t)nm];
      size_t b = t.find_first_not_of(" \t"), i = 0;
      while (b + i < t.size() && i < a.size() && is_alpha(t[b + i]) &&
             to_upper(t[b + i]) == to_upper(a[i]))
        i++;
      if (i > 0 && (b + i == t.size() || !is_alpha(t[b + i]))) score++;
    }
    int mk = (src[AC_MODEL] >= 0) ? src[AC_MODEL] : src[AC_MODEL2];
    const std::string *cur = current_model_key();
    if (mk < 0 || !cur || cur->empty() || *cur == raw_[(size_t)mk]) score++;
    if (pl.end_padding()) score++;          /* a cut-off row is the usual case */
    return score;
  }

  bool same_fields(const int *a, const int *b) const {
    for (int k = 0; k < AC_NUM; k++) {
      if (k == AC_ID || k == AC_ENTITY) continue;
      if ((a[k] < 0) != (b[k] < 0)) return false;
      if (a[k] >= 0 && raw_[(size_t)a[k]] != raw_[(size_t)b[k]]) return false;
    }
    return true;
  }

  /* A row with missing or extra values. Missing values are tried at the end
   * (a cut-off row), as one gap anywhere, or (up to 3 of them) as a gap plus
   * the end; up to 3 extra values as one run anywhere. Each placement is
   * checked against the column types and the complete rows read so far. The
   * row is used only if one reading is the most plausible: wrong coordinates
   * are worse than a missing atom. */
  bool align(size_t nvals) {
    const size_t kMaxShift = 3, kMaxCols = 64;
    if (ncols_ == 0 || nvals == 0 || nvals > raw_.size()) return false;
    int best = -1;
    bool tied = false;
    int best_src[AC_NUM], src[AC_NUM];
    Placement pl;
    pl.ncols = ncols_;
    std::vector<Placement> cands;
    if (nvals < ncols_) {
      size_t m = ncols_ - nvals;
      pl.extra = 0;
      pl.at = nvals; pl.gap = 0; pl.tail = m;
      cands.push_back(pl);                    /* missing at the end */
      if (ncols_ <= kMaxCols) {
        for (size_t g = 1; g <= m; g++) {
          if (g < m && m > kMaxShift) break;
          for (size_t at = 0; at < nvals; at++) {
            pl.at = at; pl.gap = g; pl.tail = m - g;
            cands.push_back(pl);
          }
        }
      }
    } else {
      size_t e = nvals - ncols_;
      if (e > kMaxShift || ncols_ > kMaxCols) return false;
      pl.gap = pl.tail = 0;
      pl.extra = e;
      for (size_t at = 0; at <= ncols_; at++) { pl.at = at; cands.push_back(pl); }
    }
    for (size_t i = 0; i < cands.size(); i++) {
      int sc = evaluate(cands[i], src);
      if (sc < 0) continue;
      if (sc > best) {
        best = sc;
        tied = false;
        memcpy(best_src, src, sizeof(src));
      } else if (sc == best && !same_fields(src, best_src)) {
        tied = true;
      }
    }
    if (best < 0 || tied) return false;
    memcpy(src_, best_src, sizeof(src_));
    return true;
  }

  const std::string &field(int k) const {
    static const std::string empty;
    return have_[k] ? val_[k] : empty;
  }

  int select_model() {
    const std::string *key = NULL;
    if (have_[AC_MODEL]) key = &val_[AC_MODEL];
    else if (have_[AC_MODEL2]) key = &val_[AC_MODEL2];

    if (cur_model_ < 0) {                     /* first row */
      cur_model_ = 0;
      if (key) { model0_key = *key; model_index_[*key] = 0; }
      return 0;
    }
    if (!key) return cur_model_;              /* missing model id: same model */

    std::string &cur_key = (cur_model_ == 0) ? model0_key
                                             : extra[(size_t)cur_model_ - 1].key;
    if (cur_key == *key) return cur_model_;
    if (cur_key.empty()) {                    /* adopt the first real id */
      cur_key = *key;
      model_index_[*key] = cur_model_;
      return cur_model_;
    }
    std::unordered_map<std::string, int>::const_iterator it = model_index_.find(*key);
    if (it != model_index_.end()) return cur_model_ = it->second;
    if (extra.size() >= kMaxModels) return -1;
    extra.push_back(ExtraModel());
    extra.back().key = *key;
    cur_model_ = (int)extra.size();
    model_index_[*key] = cur_model_;
    return cur_model_;
  }

  void process_row() {
    float x, y, z;
    if (!(have_[AC_X] && have_[AC_Y] && have_[AC_Z] &&
          parse_float(val_[AC_X], &x) && parse_float(val_[AC_Y], &y) &&
          parse_float(val_[AC_Z], &z))) {
      bad_coord_rows++;
      return;
    }
    if (have_[AC_GROUP] && !iequals(val_[AC_GROUP], "ATOM") &&
        !iequals(val_[AC_GROUP], "HETATM"))
      bad_group_rows++;
    bool field_ok = !(have_[AC_INS] && val_[AC_INS].size() > 1) &&
                    !(have_[AC_TYPE] && !is_alpha(val_[AC_TYPE][0]));

    /* identity fields shared by the structure and later models */
    const std::string &name = have_[AC_LATOM] ? val_[AC_LATOM]
                            : have_[AC_AATOM] ? val_[AC_AATOM]
                            : field(AC_TYPE);
    const std::string &resname = have_[AC_LCOMP] ? val_[AC_LCOMP] : field(AC_ACOMP);
    int chain_primary = use_label_ ? AC_LASYM : AC_AASYM;
    int chain_second  = use_label_ ? AC_AASYM : AC_LASYM;
    const std::string &chain = have_[chain_primary] ? val_[chain_primary]
                                                    : field(chain_second);
    int seq_primary = use_label_ ? AC_LSEQ : AC_ASEQ;
    int seq_second  = use_label_ ? AC_ASEQ : AC_LSEQ;
    int resid = 0;
    char seq_ins = 0;                /* letter appended to the number, "100A" */
    bool have_resid =
        (have_[seq_primary] && parse_seq(val_[seq_primary], &resid, &seq_ins)) ||
        (have_[seq_second] && parse_seq(val_[seq_second], &resid, &seq_ins));
    if (!have_resid) resid = 0;
    if (have_[seq_primary] && !have_resid) field_ok = false;
    const std::string &alt = field(AC_ALT);

    unsigned id = 2166136261u;
    id = fnv1a(id, name.data(), name.size());
    id = fnv1a(id, resname.data(), resname.size());
    id = fnv1a(id, chain.data(), chain.size());
    id = fnv1a(id, (const char *)&resid, sizeof(resid));
    id = fnv1a(id, alt.data(), alt.size());

    float occ = 1.0f, bfac = 0.0f, charge = 0.0f;
    bool has_charge = false;
    if (have_[AC_OCC] && !parse_float(val_[AC_OCC], &occ)) field_ok = false;
    if (have_[AC_B] && !parse_float(val_[AC_B], &bfac)) field_ok = false;
    if (have_[AC_CHARGE]) {
      if (parse_charge(val_[AC_CHARGE], &charge)) has_charge = true;
      else field_ok = false;
    }
    if (!field_ok) bad_field_rows++;

    int model = select_model();
    if (model < 0) { dropped_model_rows++; return; }

    /* Chain selection: after the counts above, which judge how the loop was
     * read and must not depend on it, and after the model of the row is
     * known, so that the first model still defines the structure: the
     * selection gives exactly the selected atoms of the unselected load.
     * Later models are compared on the selected atoms only. */
    if (filter_) {
      int k = filter_->find(chain);
      if (model == 0) {
        if (k >= 0) found[(size_t)k] = 1;
        if (chains_seen.empty() || chains_seen.back() != chain) {
          if (std::find(chains_seen.begin(), chains_seen.end(), chain) == chains_seen.end()) {
            if (chains_seen.size() < 50) chains_seen.push_back(chain);
            else more_chains = true;
          }
        }
      }
      if (k < 0) { filtered_rows++; return; }
    }
    if (model > 0) {
      ExtraModel &m = extra[(size_t)model - 1];
      size_t k = m.rows++;
      if (!m.matches) return;
      if (k >= ident.size() || ident[k] != id) {
        m.matches = false;
        m.xyz.clear();
        return;
      }
      m.xyz.push_back(x); m.xyz.push_back(y); m.xyz.push_back(z);
      return;
    }

    if (atoms.size() >= (size_t)(INT_MAX / 4))
      throw std::runtime_error("too many atoms");

    molfile_atom_t a;
    memset(&a, 0, sizeof(a));
    if (copy_field(a.name, sizeof(a.name), name)) long_name++;
    copy_field(a.type, sizeof(a.type), have_[AC_TYPE] ? val_[AC_TYPE] : name);
    if (copy_field(a.resname, sizeof(a.resname), resname)) long_resname++;
    a.resid = resid;
    if (!have_resid && (has_col[AC_ASEQ] || has_col[AC_LSEQ])) missing_resid++;

    /* the ABI's chain field holds one character; the full ID goes to segid */
    if (!chain.empty()) {
      a.chain[0] = chain[0];
      if (chain.size() > 1) multichar_chain++;
    }
    if (copy_field(a.segid, sizeof(a.segid), chain)) long_chain++;

    if (!alt.empty()) {
      a.altloc[0] = alt[0];
      any_altloc = true;
      if (alt.size() > 1) long_altloc++;
    }
    /* insertion codes belong to the author numbering; label_seq_id is unique
     * without them */
    char ins = ' ';
    if (have_[AC_INS] && !val_[AC_INS].empty() && val_[AC_INS][0] != ' ')
      ins = val_[AC_INS][0];
    else if (seq_ins && !use_label_)
      ins = seq_ins;
    a.insertion[0] = use_label_ ? ' ' : ins;
    if (a.insertion[0] != ' ') any_insertion = true;

    a.occupancy = occ;
    a.bfactor = bfac;
    a.charge = charge;
    if (has_charge) any_charge = true;

    int zn = have_[AC_TYPE] ? element_from_symbol(val_[AC_TYPE]) : -1;
    if (zn > 0) {
      n_known_elem++;
      a.atomicnumber = zn;
    } else if (zn == 0) {
      n_unknown_elem++;
      a.atomicnumber = 0;
      if (unknown_symbols.size() < 32 || unknown_symbols.count(val_[AC_TYPE]))
        unknown_symbols[val_[AC_TYPE]]++;
    } else {
      n_no_elem++;
      a.atomicnumber = 0;
    }
    /* mass and radius are left for VMD to derive from the element (or, for
     * atomic number 0, from the atom name) */

    atoms.push_back(a);
    xyz.push_back(x); xyz.push_back(y); xyz.push_back(z);
    ident.push_back(id);

    insertion.push_back(ins == ' ' ? '\0' : ins);
    auth_asym.push_back(have_[AC_AASYM] ? asyms.intern(val_[AC_AASYM]) : -1);
    label_asym.push_back(have_[AC_LASYM] ? asyms.intern(val_[AC_LASYM]) : -1);
    int s;
    auth_seq.push_back((have_[AC_ASEQ] && parse_seq(val_[AC_ASEQ], &s, NULL)) ? s : kSeqMissing);
    label_seq.push_back((have_[AC_LSEQ] && parse_int(val_[AC_LSEQ], &s)) ? s : kSeqMissing);
  }
};

/* ========================================================================
 * BinaryCIF input
 *
 * BinaryCIF (as served by RCSB, PDBe and Mol*) is the same data model as
 * mmCIF, stored column by column in MessagePack and compressed with a chain
 * of simple encodings. The file is read whole, parsed into a flat list of
 * MessagePack values that point into the file buffer, and only the columns
 * the sinks use are decoded. Each value is then written back as the text a
 * CIF file holds, so everything after this (elements, chains, models, bonds,
 * the chain selection) works exactly as for text. Every length, offset and
 * count in the file is checked, and the decoded size is capped relative to
 * the file size, so a small file cannot claim billions of values.
 * ======================================================================== */

const size_t kMaxBcifBytes = (size_t)1 << 30;   /* the whole file, decompressed */
const size_t kMaxMpValues = (size_t)1 << 20;    /* MessagePack values in the file */
const int kMaxMpDepth = 32;                     /* BinaryCIF nests about 8 deep */
const long long kMaxBcifRows = (long long)1 << 28;

void bcif_fail(const std::string &what) {
  throw std::runtime_error("BinaryCIF: " + what);
}

/* One MessagePack value. Strings and binary data point into the file
 * buffer; an array's elements (a map's keys and values, alternating) are
 * `count` consecutive entries of the document from `first` on. */
struct MpVal {
  enum Type { NIL, BOOL, INT, FLOAT, STR, BIN, ARRAY, MAP };
  unsigned char type;
  bool b;
  long long i;
  double f;
  const char *p;
  size_t n;
  size_t first, count;
  MpVal() : type(NIL), b(false), i(0), f(0), p(NULL), n(0), first(0), count(0) {}
};

class MpDoc {
 public:
  std::vector<MpVal> v;                       /* v[0] is the top-level value */

  void parse(const char *data, size_t size) {
    p_ = (const unsigned char *)data;
    n_ = size;
    pos_ = 0;
    v.assign(1, MpVal());
    fill(0, 0);
  }

  const MpVal &root() const { return v[0]; }
  const MpVal &elem(const MpVal &a, size_t k) const { return v[a.first + k]; }

  /* the value of a map entry, or NULL */
  const MpVal *get(const MpVal &m, const char *key) const {
    if (m.type != MpVal::MAP) return NULL;
    const size_t kl = strlen(key);
    for (size_t k = 0; k < m.count; k++) {
      const MpVal &key_val = v[m.first + 2 * k];
      if (key_val.type == MpVal::STR && key_val.n == kl && memcmp(key_val.p, key, kl) == 0)
        return &v[m.first + 2 * k + 1];
    }
    return NULL;
  }

 private:
  const unsigned char *p_;
  size_t n_, pos_;

  void need(unsigned long long k) {
    if (k > n_ - pos_) bcif_fail("the data ends in the middle of a value (truncated file?)");
  }
  unsigned long long be(int bytes) {          /* big-endian unsigned */
    need((unsigned long long)bytes);
    unsigned long long u = 0;
    for (int k = 0; k < bytes; k++) u = (u << 8) | p_[pos_++];
    return u;
  }
  static long long sign_extend(unsigned long long u, int bytes) {
    if (bytes < 8 && (u >> (8 * bytes - 1)) & 1) u |= ~0ULL << (8 * bytes);
    long long s;
    memcpy(&s, &u, sizeof(s));
    return s;
  }

  void fill(size_t idx, int depth) {
    if (depth > kMaxMpDepth) bcif_fail("values are nested too deeply");
    need(1);
    const unsigned c = p_[pos_++];
    MpVal &x = v[idx];
    if (c <= 0x7f) { x.type = MpVal::INT; x.i = (long long)c; return; }
    if (c >= 0xe0) { x.type = MpVal::INT; x.i = (long long)c - 256; return; }
    if (c >= 0xa0 && c <= 0xbf) { bytes(idx, MpVal::STR, c & 0x1f); return; }
    if (c >= 0x90 && c <= 0x9f) { container(idx, MpVal::ARRAY, c & 0x0f, depth); return; }
    if (c >= 0x80 && c <= 0x8f) { container(idx, MpVal::MAP, c & 0x0f, depth); return; }
    switch (c) {
      case 0xc0: x.type = MpVal::NIL; return;
      case 0xc2: case 0xc3: x.type = MpVal::BOOL; x.b = (c == 0xc3); return;
      case 0xcc: x.type = MpVal::INT; x.i = (long long)be(1); return;
      case 0xcd: x.type = MpVal::INT; x.i = (long long)be(2); return;
      case 0xce: x.type = MpVal::INT; x.i = (long long)be(4); return;
      case 0xcf: {
        unsigned long long u = be(8);
        if (u > (unsigned long long)LLONG_MAX) bcif_fail("an integer is too large");
        v[idx].type = MpVal::INT; v[idx].i = (long long)u;
        return;
      }
      case 0xd0: x.type = MpVal::INT; x.i = sign_extend(be(1), 1); return;
      case 0xd1: x.type = MpVal::INT; x.i = sign_extend(be(2), 2); return;
      case 0xd2: x.type = MpVal::INT; x.i = sign_extend(be(4), 4); return;
      case 0xd3: x.type = MpVal::INT; x.i = sign_extend(be(8), 8); return;
      case 0xca: {
        uint32_t u = (uint32_t)be(4);
        float f;
        memcpy(&f, &u, sizeof(f));
        v[idx].type = MpVal::FLOAT; v[idx].f = f;
        return;
      }
      case 0xcb: {
        unsigned long long u = be(8);
        double d;
        memcpy(&d, &u, sizeof(d));
        v[idx].type = MpVal::FLOAT; v[idx].f = d;
        return;
      }
      case 0xd9: bytes(idx, MpVal::STR, be(1)); return;
      case 0xda: bytes(idx, MpVal::STR, be(2)); return;
      case 0xdb: bytes(idx, MpVal::STR, be(4)); return;
      case 0xc4: bytes(idx, MpVal::BIN, be(1)); return;
      case 0xc5: bytes(idx, MpVal::BIN, be(2)); return;
      case 0xc6: bytes(idx, MpVal::BIN, be(4)); return;
      case 0xdc: container(idx, MpVal::ARRAY, be(2), depth); return;
      case 0xdd: container(idx, MpVal::ARRAY, be(4), depth); return;
      case 0xde: container(idx, MpVal::MAP, be(2), depth); return;
      case 0xdf: container(idx, MpVal::MAP, be(4), depth); return;
      default:
        bcif_fail("unsupported MessagePack type (BinaryCIF uses no extension types)");
    }
  }

  void bytes(size_t idx, unsigned char type, unsigned long long len) {
    need(len);
    v[idx].type = type;
    v[idx].p = (const char *)p_ + pos_;
    v[idx].n = (size_t)len;
    pos_ += (size_t)len;
  }

  void container(size_t idx, unsigned char type, unsigned long long count, int depth) {
    const unsigned long long slots = (type == MpVal::MAP) ? 2 * count : count;
    need(slots);                              /* every value takes a byte or more */
    if (v.size() + slots > kMaxMpValues)
      bcif_fail("too many values (more than any BinaryCIF file has)");
    const size_t first = v.size();
    v.resize(first + (size_t)slots);         /* x would dangle: index from here */
    v[idx].type = type;
    v[idx].first = first;
    v[idx].count = (size_t)count;
    for (size_t k = 0; k < (size_t)slots; k++) fill(first + k, depth + 1);
  }
};

/* A decoded array: integers (optionally fixed point), reals, or strings. */
struct BcifArray {
  enum Kind { INTS, REALS, STRINGS };
  Kind kind;
  std::vector<int32_t> ints;       /* INTS; STRINGS: index into the table, -1 = "" */
  int decimals;                    /* INTS: value = ints / 10^decimals; -1: integers */
  std::vector<double> reals;
  bool single;                     /* REALS hold float32 values */
  std::vector<const char *> sp;    /* STRINGS: the table, pointing into the file */
  std::vector<size_t> sn;
  BcifArray() : kind(INTS), decimals(-1), single(false) {}
  size_t size() const { return kind == REALS ? reals.size() : ints.size(); }
};

class BcifDecoder {
 public:
  BcifDecoder(const MpDoc &d, size_t file_bytes)
    : d_(d), spent_(0),
      budget_((long long)64 * (long long)file_bytes + ((long long)16 << 20)),
      emitted_(0),
      emit_budget_((long long)64 * (long long)file_bytes + ((long long)64 << 20)) {}

  /* Bytes of text given to the sinks, capped like the decoded values: an
   * index column can repeat a long string in every row, and a category can
   * claim many rows, so the text is not bounded by the file otherwise. */
  void emit(unsigned long long bytes) {
    if (bytes > (unsigned long long)emit_budget_ ||
        (emitted_ += (long long)bytes) > emit_budget_)
      bcif_fail("the decoded text is implausibly large for the file size");
  }

  /* Decode {data, encoding} into exactly `len` values. */
  BcifArray column(const MpVal &enc_data, size_t len, const std::string &what) {
    what_ = what;
    if (enc_data.type != MpVal::MAP) fail("the column data is not a map");
    const MpVal *data = d_.get(enc_data, "data");
    const MpVal *enc = d_.get(enc_data, "encoding");
    if (!data || (data->type != MpVal::BIN && data->type != MpVal::STR))
      fail("no binary data");
    if (!enc || enc->type != MpVal::ARRAY) fail("no list of encodings");
    BcifArray a = chain(*data, *enc, len);
    if (a.size() != len) fail("the number of values does not match the row count");
    return a;
  }

 private:
  const MpDoc &d_;
  long long spent_;
  const long long budget_;
  long long emitted_;
  const long long emit_budget_;
  std::string what_;

  void fail(const std::string &why) const { bcif_fail(what_ + ": " + why); }

  /* Account for n new values; also caps one step at a few times the rows. */
  void spend(unsigned long long n, size_t len) {
    if (n > 8ULL * len + 1024) fail("an encoding step produces far more values than there are rows");
    spent_ += (long long)n;
    if (spent_ > budget_) fail("the decoded data is implausibly large for the file size");
  }

  long long num_int(const MpVal &e, const char *key) const {
    const MpVal *x = d_.get(e, key);
    if (x && x->type == MpVal::INT) return x->i;
    if (x && x->type == MpVal::FLOAT && x->f == std::floor(x->f) && std::fabs(x->f) < 9.0e15)
      return (long long)x->f;
    fail(std::string("an encoding has no integer '") + key + "'");
    return 0;
  }
  double num_real(const MpVal &e, const char *key) const {
    const MpVal *x = d_.get(e, key);
    if (x && x->type == MpVal::INT) return (double)x->i;
    if (x && x->type == MpVal::FLOAT && std::isfinite(x->f)) return x->f;
    fail(std::string("an encoding has no number '") + key + "'");
    return 0;
  }
  static bool kind_is(const MpVal *k, const char *name) {
    return k && k->type == MpVal::STR && k->n == strlen(name) && memcmp(k->p, name, k->n) == 0;
  }
  void need_ints(const BcifArray &a, const char *step) const {
    if (a.kind != BcifArray::INTS || a.decimals >= 0)
      fail(std::string(step) + " applied to values that are not integers");
  }
  static void push_checked(std::vector<int32_t> &out, long long x, const BcifDecoder *self) {
    if (x < INT32_MIN || x > INT32_MAX) self->fail("an integer does not fit in 32 bits");
    out.push_back((int32_t)x);
  }

  /* Undo the encodings, last first. The first step takes the raw bytes. */
  BcifArray chain(const MpVal &raw, const MpVal &enc, size_t len) {
    if (enc.count == 0) fail("empty list of encodings");
    BcifArray a;
    bool have_bytes = true;
    for (size_t k = enc.count; k-- > 0;) {
      const MpVal &e = d_.elem(enc, k);
      if (e.type != MpVal::MAP) fail("an encoding is not a map");
      const MpVal *kind = d_.get(e, "kind");
      static const char *const known[] = { "ByteArray", "StringArray", "FixedPoint",
        "IntervalQuantization", "RunLength", "Delta", "IntegerPacking", NULL };
      bool is_known = false;
      for (const char *const *kn = known; *kn && !is_known; kn++) is_known = kind_is(kind, *kn);
      if (!is_known) {
        if (kind && kind->type == MpVal::STR)
          fail("unsupported encoding '" + printable(std::string(kind->p, kind->n)) + "'");
        fail("an encoding has no kind");
      }
      if (kind_is(kind, "ByteArray")) {
        if (!have_bytes) fail("ByteArray is not the first step");
        byte_array(raw, num_int(e, "type"), len, a);
        have_bytes = false;
      } else if (kind_is(kind, "StringArray")) {
        if (!have_bytes || k != 0) fail("StringArray is not the only step");
        string_array(raw, e, len, a);
        have_bytes = false;
      } else {
        if (have_bytes) fail("the first step is not ByteArray or StringArray");
        if (kind_is(kind, "FixedPoint")) fixed_point(e, a);
        else if (kind_is(kind, "IntervalQuantization")) interval(e, a);
        else if (kind_is(kind, "RunLength")) run_length(e, len, a);
        else if (kind_is(kind, "Delta")) delta(e, a);
        else integer_packing(e, len, a);
      }
    }
    if (have_bytes) fail("no decoding step");
    return a;
  }

  void byte_array(const MpVal &raw, long long type, size_t len, BcifArray &a) {
    size_t size = 0;
    switch (type) {
      case 1: case 4: size = 1; break;
      case 2: case 5: size = 2; break;
      case 3: case 6: case 32: size = 4; break;
      case 33: size = 8; break;
      default: fail("unknown ByteArray type");
    }
    if (raw.n % size != 0) fail("the byte count is not a whole number of values");
    const size_t n = raw.n / size;
    spend(n, len);
    const unsigned char *b = (const unsigned char *)raw.p;
    a = BcifArray();
    if (type == 32 || type == 33) {
      a.kind = BcifArray::REALS;
      a.single = (type == 32);
      a.reals.resize(n);
      for (size_t i = 0; i < n; i++) {
        unsigned long long u = 0;
        for (size_t k = size; k-- > 0;) u = (u << 8) | b[i * size + k];   /* little-endian */
        if (type == 32) {
          uint32_t u32 = (uint32_t)u;
          float f;
          memcpy(&f, &u32, sizeof(f));
          a.reals[i] = f;
        } else {
          double dv;
          memcpy(&dv, &u, sizeof(dv));
          a.reals[i] = dv;
        }
      }
      return;
    }
    a.ints.reserve(n);
    for (size_t i = 0; i < n; i++) {
      unsigned long long u = 0;
      for (size_t k = size; k-- > 0;) u = (u << 8) | b[i * size + k];
      long long x;
      if (type <= 3) {                          /* signed */
        const unsigned long long sign = 1ULL << (8 * size - 1);
        x = (u & sign) ? (long long)u - (long long)(sign << 1) : (long long)u;
      } else {
        x = (long long)u;
      }
      push_checked(a.ints, x, this);
    }
  }

  void fixed_point(const MpVal &e, BcifArray &a) {
    need_ints(a, "FixedPoint");
    const double factor = num_real(e, "factor");
    if (!(factor > 0)) fail("FixedPoint factor is not positive");
    const bool single = d_.get(e, "srcType") && num_int(e, "srcType") == 32;
    for (int k = 0; k <= 9; k++) {
      if (factor == kPow10[k]) {                /* exact decimals: keep the integers */
        a.decimals = k;
        a.single = single;
        return;
      }
    }
    BcifArray r;
    r.kind = BcifArray::REALS;
    r.single = single;
    r.reals.resize(a.ints.size());
    for (size_t i = 0; i < a.ints.size(); i++) r.reals[i] = a.ints[i] / factor;
    a.ints.clear();
    a = r;
  }

  void interval(const MpVal &e, BcifArray &a) {
    need_ints(a, "IntervalQuantization");
    const double lo = num_real(e, "min"), hi = num_real(e, "max");
    const long long steps = num_int(e, "numSteps");
    if (steps < 2) fail("IntervalQuantization has fewer than 2 steps");
    const double delta = (hi - lo) / (double)(steps - 1);
    BcifArray r;
    r.kind = BcifArray::REALS;
    r.single = d_.get(e, "srcType") && num_int(e, "srcType") == 32;
    r.reals.resize(a.ints.size());
    for (size_t i = 0; i < a.ints.size(); i++) r.reals[i] = lo + delta * a.ints[i];
    a = r;
  }

  void run_length(const MpVal &e, size_t len, BcifArray &a) {
    need_ints(a, "RunLength");
    const long long src = num_int(e, "srcSize");
    if (src < 0) fail("RunLength srcSize is negative");
    if (a.ints.size() % 2 != 0) fail("RunLength data is not (value, count) pairs");
    spend((unsigned long long)src, len);
    std::vector<int32_t> out;
    out.reserve((size_t)src);
    for (size_t i = 0; i < a.ints.size(); i += 2) {
      const long long count = a.ints[i + 1];
      if (count < 0 || (long long)out.size() + count > src)
        fail("RunLength counts do not add up to srcSize");
      out.insert(out.end(), (size_t)count, a.ints[i]);
    }
    if ((long long)out.size() != src) fail("RunLength counts do not add up to srcSize");
    a.ints.swap(out);
  }

  void delta(const MpVal &e, BcifArray &a) {
    need_ints(a, "Delta");
    long long acc = num_int(e, "origin");
    for (size_t i = 0; i < a.ints.size(); i++) {
      acc += a.ints[i];
      if (acc < INT32_MIN || acc > INT32_MAX) fail("Delta values do not fit in 32 bits");
      a.ints[i] = (int32_t)acc;
    }
  }

  void integer_packing(const MpVal &e, size_t len, BcifArray &a) {
    need_ints(a, "IntegerPacking");
    const long long bytes = num_int(e, "byteCount");
    const MpVal *u = d_.get(e, "isUnsigned");
    const bool is_unsigned = u && u->type == MpVal::BOOL && u->b;
    const long long src = num_int(e, "srcSize");
    if (bytes != 1 && bytes != 2) fail("IntegerPacking byteCount is not 1 or 2");
    if (src < 0) fail("IntegerPacking srcSize is negative");
    spend((unsigned long long)src, len);
    const long long upper = is_unsigned ? (bytes == 1 ? 0xFF : 0xFFFF) : (bytes == 1 ? 0x7F : 0x7FFF);
    const long long lower = -upper - 1;
    std::vector<int32_t> out;
    out.reserve((size_t)src);
    const size_t n = a.ints.size();
    for (size_t i = 0; i < n; i++) {
      long long value = 0;
      long long t = a.ints[i];
      /* a limit value continues the number into the next one */
      while (t == upper || (!is_unsigned && t == lower)) {
        value += t;
        if (++i >= n) fail("IntegerPacking data ends inside a value");
        t = a.ints[i];
      }
      value += t;
      if ((long long)out.size() >= src) fail("IntegerPacking gives more values than srcSize");
      push_checked(out, value, this);
    }
    if ((long long)out.size() != src) fail("IntegerPacking gives fewer values than srcSize");
    a.ints.swap(out);
  }

  /* String offsets count characters; convert them to byte offsets when the
   * strings are not plain ASCII (UTF-8: count the bytes that start one). */
  void string_array(const MpVal &raw, const MpVal &e, size_t len, BcifArray &a) {
    const MpVal *sd = d_.get(e, "stringData");
    const MpVal *offs = d_.get(e, "offsets");
    const MpVal *oenc = d_.get(e, "offsetEncoding");
    const MpVal *denc = d_.get(e, "dataEncoding");
    if (!sd || sd->type != MpVal::STR) fail("StringArray has no stringData");
    if (!offs || (offs->type != MpVal::BIN && offs->type != MpVal::STR))
      fail("StringArray has no offsets");
    if (!oenc || oenc->type != MpVal::ARRAY || !denc || denc->type != MpVal::ARRAY)
      fail("StringArray has no offset or data encoding");

    const std::string saved = what_;
    what_ = saved + " (string offsets)";
    BcifArray off = chain(*offs, *oenc, sd->n + 1);
    what_ = saved + " (string indices)";
    BcifArray idx = chain(raw, *denc, len);
    what_ = saved;
    need_ints(off, "StringArray offsets");
    need_ints(idx, "StringArray indices");
    if (off.ints.empty()) fail("StringArray has no offsets");

    std::vector<size_t> byte_at;               /* character index -> byte offset */
    bool ascii = true;
    for (size_t i = 0; i < sd->n && ascii; i++)
      if ((unsigned char)sd->p[i] >= 0x80) ascii = false;
    if (!ascii) {
      for (size_t i = 0; i < sd->n; i++)
        if (((unsigned char)sd->p[i] & 0xC0) != 0x80) byte_at.push_back(i);
      byte_at.push_back(sd->n);
    }
    const size_t nchars = ascii ? sd->n : byte_at.size() - 1;

    a = BcifArray();
    a.kind = BcifArray::STRINGS;
    const size_t nstr = off.ints.size() - 1;
    a.sp.resize(nstr);
    a.sn.resize(nstr);
    for (size_t k = 0; k < nstr; k++) {
      const long long b0 = off.ints[k], b1 = off.ints[k + 1];
      if (b0 < 0 || b1 < b0 || (size_t)b1 > nchars) fail("a string offset is out of range");
      const size_t s0 = ascii ? (size_t)b0 : byte_at[(size_t)b0];
      const size_t s1 = ascii ? (size_t)b1 : byte_at[(size_t)b1];
      if (s1 - s0 > kMaxTokenBytes) fail("a string is larger than 256 MiB");
      a.sp[k] = sd->p + s0;
      a.sn[k] = s1 - s0;
    }
    for (size_t i = 0; i < idx.ints.size(); i++)
      if (idx.ints[i] < -1 || idx.ints[i] >= (long long)nstr) fail("a string index is out of range");
    a.ints.swap(idx.ints);
  }
};

/* A wanted column of a BinaryCIF category. */
struct BcifColumn {
  std::string item;                /* lowercase */
  BcifArray values;
  std::vector<unsigned char> mask; /* empty, or 0 = value, 1 = '.', 2 = '?' */
};

/* Text of an integer: digits, or fixed point with `decimals` places. */
size_t format_fixed(long long v, int decimals, char *buf) {
  char tmp[32];
  size_t n = 0;
  const bool neg = v < 0;
  unsigned long long u = neg ? 0ULL - (unsigned long long)v : (unsigned long long)v;
  int digits = 0;
  do {                              /* digits from the right; the loop runs past */
    tmp[n++] = (char)('0' + (int)(u % 10));   /* the point once, so 5 with 3 */
    u /= 10;                                  /* decimals gives "0.005" */
    if (++digits == decimals) tmp[n++] = '.';
  } while (u > 0 || digits <= decimals);
  size_t o = 0;
  if (neg) buf[o++] = '-';
  while (n > 0) buf[o++] = tmp[--n];
  return o;
}

/* Shortest text that reads back as the same number: a value that came from
 * a short decimal ("17.047") is written that way again, so BinaryCIF reads
 * exactly like the text file it was made from. */
size_t format_real(double v, bool single, char *buf) {
  for (int d = 0; d <= 6; d++) {
    const double s = v * kPow10[d];
    if (!(std::fabs(s) < 9.0e15)) break;
    const double r = (double)std::llround(s);
    const double back = r / kPow10[d];
    if (single ? ((float)back == (float)v) : (back == v)) {
      if (r == 0 && std::signbit(v)) {              /* keep the sign of -0.0 */
        buf[0] = '-';
        return 1 + format_fixed(0, d, buf + 1);
      }
      return format_fixed((long long)r, d, buf);
    }
  }
  int n = snprintf(buf, 40, single ? "%.9g" : "%.17g", v);
  if (n < 0 || n >= 40) return 0;
  for (int i = 0; i < n; i++)
    if (buf[i] == ',') buf[i] = '.';                 /* a comma-decimal locale */
  return (size_t)n;
}

/* ========================================================================
 * Parser: drives the tokenizer, routes categories to sinks
 * ======================================================================== */

class Parser {
 public:
  explicit Parser(CifReader &r)
    : stray_values(0), extra_blocks(false), reader_(r), has_pending_(false),
      seen_data_(false) {}

  void add_sink(const char *category, Sink *s) { sinks_[category] = s; }

  /* Parse data blocks until one has content according to has_content();
   * blocks without content are discarded. */
  template <class Pred> void run(Pred has_content) {
    for (;;) {
      begin_block();
      bool more = parse_block();
      end_block();
      if (has_content()) { extra_blocks = more; return; }
      if (!more) return;
    }
  }

  /* The same for BinaryCIF: each category is a loop whose rows are given to
   * the sink as text, column by column, as a CIF file would hold them. */
  template <class Pred> void run_binary(const MpDoc &doc, size_t file_bytes,
                                        Pred has_content) {
    BcifDecoder dec(doc, file_bytes);
    const MpVal *blocks = doc.get(doc.root(), "dataBlocks");
    if (!blocks || blocks->type != MpVal::ARRAY) bcif_fail("no list of data blocks");
    for (size_t b = 0; b < blocks->count; b++) {
      begin_block();
      const MpVal &blk = doc.elem(*blocks, b);
      const MpVal *cats = doc.get(blk, "categories");
      if (!cats || cats->type != MpVal::ARRAY) bcif_fail("a data block has no categories");
      for (size_t c = 0; c < cats->count; c++)
        binary_category(doc, dec, doc.elem(*cats, c));
      end_block();
      if (has_content()) { extra_blocks = (b + 1 < blocks->count); return; }
    }
  }

  /* loops re-read line by line: (category, line of the first value) */
  std::vector<std::pair<std::string, long> > recovered;
  long stray_values;
  bool extra_blocks;

 private:
  struct KVTable {
    std::vector<std::string> items;
    std::vector<std::string> vals;
    std::vector<char> nulls;
    std::map<std::string, size_t> index;    /* a repeated tag: the last wins */
  };

  CifReader &reader_;
  std::map<std::string, Sink *> sinks_;
  std::map<std::string, KVTable> kv_;
  std::set<std::string> looped_;
  Token pending_;
  ReaderState pending_state_;
  ReaderState tok_state_;      /* reader state before the token last returned */
  bool has_pending_;
  bool seen_data_;

  Token get() {
    if (has_pending_) {
      has_pending_ = false;
      tok_state_ = pending_state_;
      return pending_;
    }
    tok_state_ = reader_.state();
    return reader_.next();
  }
  void unget(const Token &t) {
    pending_ = t;
    pending_state_ = tok_state_;
    has_pending_ = true;
  }

  Sink *sink_for(const std::string &cat) {
    std::map<std::string, Sink *>::iterator it = sinks_.find(cat);
    return (it == sinks_.end()) ? NULL : it->second;
  }

  /* one BinaryCIF category: decode the columns its sink uses, feed the rows */
  void binary_category(const MpDoc &doc, BcifDecoder &dec, const MpVal &cat) {
    const MpVal *name = doc.get(cat, "name");
    if (!name || name->type != MpVal::STR) bcif_fail("a category has no name");
    std::string cname = lowercase(name->p, name->n);
    if (!cname.empty() && cname[0] == '_') cname.erase(0, 1);   /* "_atom_site" */
    Sink *s = sink_for(cname);
    if (!s) return;
    const MpVal *rc = doc.get(cat, "rowCount");
    const MpVal *cols = doc.get(cat, "columns");
    if (!rc || rc->type != MpVal::INT || rc->i < 0 || rc->i > kMaxBcifRows)
      bcif_fail("_" + printable(cname) + " has no valid row count");
    if (!cols || cols->type != MpVal::ARRAY) bcif_fail("_" + printable(cname) + " has no columns");
    const size_t nrows = (size_t)rc->i;

    std::vector<BcifColumn> use;
    for (size_t k = 0; k < cols->count; k++) {
      const MpVal &col = doc.elem(*cols, k);
      const MpVal *cn = doc.get(col, "name");
      if (!cn || cn->type != MpVal::STR) bcif_fail("_" + printable(cname) + " has a column without a name");
      const std::string item = lowercase(cn->p, cn->n);
      if (!s->wants(item)) continue;
      const std::string what = "_" + printable(cname) + "." + printable(item);
      const MpVal *data = doc.get(col, "data");
      if (!data) bcif_fail(what + ": no data");
      use.push_back(BcifColumn());
      use.back().item = item;
      use.back().values = dec.column(*data, nrows, what);
      const MpVal *mask = doc.get(col, "mask");
      if (mask && mask->type != MpVal::NIL) {
        BcifArray m = dec.column(*mask, nrows, what + " (mask)");
        if (m.kind != BcifArray::INTS || m.decimals >= 0) bcif_fail(what + ": the mask is not integers");
        std::vector<unsigned char> &mk = use.back().mask;
        mk.resize(m.ints.size());
        for (size_t r = 0; r < m.ints.size(); r++)
          mk[r] = (unsigned char)(m.ints[r] == 0 ? 0 : m.ints[r] == 1 ? 1 : 2);
      }
    }

    /* every row and value costs at least a byte, charged before any row is
     * fed, so a huge row count fails at once (also with no columns used) */
    dec.emit((unsigned long long)nrows * (use.size() + 1));
    looped_.insert(cname);
    std::vector<std::string> items(use.size());
    for (size_t k = 0; k < use.size(); k++) items[k] = use[k].item;
    s->begin(items);
    char buf[48];
    for (size_t r = 0; r < nrows; r++) {
      for (size_t k = 0; k < use.size(); k++) {
        const BcifColumn &col = use[k];
        if (!col.mask.empty() && col.mask[r] != 0) {
          s->value(k, col.mask[r] == 1 ? "." : "?", 1, true);
          continue;
        }
        const BcifArray &a = col.values;
        if (a.kind == BcifArray::STRINGS) {
          const int32_t ix = a.ints[r];
          if (ix < 0) { s->value(k, "", 0, false); continue; }   /* no string: "" */
          dec.emit(a.sn[(size_t)ix]);
          s->value(k, a.sp[(size_t)ix], a.sn[(size_t)ix], false);
        } else if (a.kind == BcifArray::REALS) {
          if (!std::isfinite(a.reals[r])) { s->value(k, "?", 1, true); continue; }
          size_t n = format_real(a.reals[r], a.single, buf);
          dec.emit(n);
          s->value(k, buf, n, n == 0);
        } else {
          size_t n = format_fixed(a.ints[r], a.decimals > 0 ? a.decimals : 0, buf);
          dec.emit(n);
          s->value(k, buf, n, false);
        }
      }
      s->end_row(use.size());
    }
  }

  static void split_tag(const char *p, size_t n, std::string &cat, std::string &item) {
    const char *dot = (const char *)memchr(p, '.', n);
    if (dot) {
      cat = lowercase(p, (size_t)(dot - p));
      item = lowercase(dot + 1, n - (size_t)(dot - p) - 1);
    } else {
      cat = lowercase(p, n);       /* DDL1-style tag without a category dot */
      item.clear();
    }
  }

  void begin_block() {
    for (std::map<std::string, Sink *>::iterator it = sinks_.begin(); it != sinks_.end(); ++it)
      it->second->reset();
    kv_.clear();
    looped_.clear();
    recovered.clear();
  }

  /* feed tag/value categories to their sinks as one-row loops */
  void end_block() {
    for (std::map<std::string, KVTable>::iterator it = kv_.begin(); it != kv_.end(); ++it) {
      if (looped_.count(it->first)) continue;
      Sink *s = sink_for(it->first);
      if (!s) continue;
      const KVTable &kv = it->second;
      s->begin(kv.items);
      for (size_t i = 0; i < kv.vals.size(); i++)
        s->value(i, kv.vals[i].data(), kv.vals[i].size(), kv.nulls[i] != 0);
      s->end_row(kv.vals.size());
    }
  }

  /* returns true if another data block follows */
  bool parse_block() {
    for (;;) {
      Token t = get();
      switch (t.kind) {
        case TK_EOF:
          return false;
        case TK_DATA:
          if (seen_data_) return true;
          seen_data_ = true;
          break;
        case TK_LOOP:
          parse_loop();
          break;
        case TK_TAG:
          parse_item(t);
          break;
        case TK_VALUE:
          stray_values++;
          break;
        default:      /* save_, global_, stop_: frame structure is ignored */
          break;
      }
    }
  }

  void parse_item(const Token &tag) {
    std::string cat, item;
    split_tag(tag.p, tag.n, cat, item);
    Token v = get();
    if (v.kind != TK_VALUE) { unget(v); return; }   /* tag without value */
    Sink *s = sink_for(cat);
    if (!s || !s->wants(item)) return;
    KVTable &kv = kv_[cat];
    std::map<std::string, size_t>::iterator it = kv.index.find(item);
    if (it == kv.index.end()) {
      kv.index[item] = kv.items.size();
      kv.items.push_back(item);
      kv.vals.push_back(std::string(v.p, v.n));
      kv.nulls.push_back(token_is_null(v) ? 1 : 0);
    } else {
      kv.vals[it->second].assign(v.p, v.n);
      kv.nulls[it->second] = token_is_null(v) ? 1 : 0;
    }
  }

  void parse_loop() {
    std::vector<std::string> items;
    std::string cat, c, item;
    Token t = get();
    while (t.kind == TK_TAG) {
      split_tag(t.p, t.n, c, item);
      if (items.empty()) cat = c;
      items.push_back(item);
      t = get();
    }
    if (items.empty() || t.kind != TK_VALUE) { unget(t); return; }

    /* a category looped twice in one block (invalid, but seen in concatenated
     * files) appends to the first loop; re-reading would discard that one */
    const bool repeated = !looped_.insert(cat).second;
    Sink *s = sink_for(cat);
    if (!s) {                                   /* not needed: skip fast */
      while (t.kind == TK_VALUE) t = get();
      unget(t);
      return;
    }

    const ReaderState start = tok_state_;       /* before the first value */
    const long first_line = t.line;
    const size_t ncols = items.size();

    long misaligned = 0;
    size_t col = strict_pass(s, items, t, &misaligned);
    const bool partial = (col != 0);          /* values do not fill whole rows */
    if (partial) s->end_row(col);
    const long strict_bad = s->bad_rows();

    /* A loop whose values do not divide into whole rows, or whose rows look
     * shifted, may have rows with missing values. Re-read it so that each
     * row starts on a new line, as virtually every mmCIF writer emits them,
     * and keep that reading only if it is cleaner: a valid loop that packs
     * several rows per line must keep its normal reading. This needs two
     * rewinds (there and possibly back), so it is done only when both are
     * sure to work: a pipe, or compressed input past its rewind budget,
     * keeps the normal reading. */
    if (!repeated && misaligned > 0 && (partial || strict_bad > 0) &&
        reader_.can_restore(start, 2) && reader_.restore(start)) {
      has_pending_ = false;
      s->reset();
      s->begin(items);
      s->set_recovery(true);
      reread_by_line(s, ncols);
      s->set_recovery(false);
      const long reread_bad = s->bad_rows();
      if (reread_bad < strict_bad || (partial && reread_bad == strict_bad) ||
          !reader_.restore(start)) {
        recovered.push_back(std::make_pair(cat, first_line));
        return;
      }
      has_pending_ = false;
      s->reset();
      t = get();
      col = strict_pass(s, items, t, &misaligned);
      if (col != 0) s->end_row(col);
    }
    unget(t);
  }

  /* Rows by value count, as the CIF grammar defines them. Leaves t at the
   * first token after the loop and returns the number of values in an
   * unfinished last row. */
  size_t strict_pass(Sink *s, const std::vector<std::string> &items, Token &t,
                     long *misaligned) {
    const size_t ncols = items.size();
    s->begin(items);
    size_t col = 0;
    long rows = 0;
    *misaligned = 0;
    while (t.kind == TK_VALUE) {
      if (col == 0 && rows > 0 && !t.bol) (*misaligned)++;
      s->value(col, t.p, t.n, token_is_null(t));
      if (++col == ncols) { s->end_row(ncols); col = 0; rows++; }
      t = get();
    }
    return col;
  }

  /* Line-aware row assembly: tokens are grouped by physical line; a line
   * that would overflow the current row starts a new row instead, and the
   * unfinished row is passed on as short. A line with more values than one
   * row but fewer than two is one row with extra values (an unquoted value
   * with a blank, say). Lines holding several rows and rows wrapped over
   * several lines are still handled. The sink decides where missing or
   * extra values were. */
  void reread_by_line(Sink *s, size_t ncols) {
    std::vector<std::string> group;
    std::vector<char> gnull;
    size_t rowvals = 0;

    Token t = get();
    for (;;) {
      bool done = (t.kind != TK_VALUE);
      if (done || t.bol) {
        if (!group.empty()) {
          if (rowvals > 0 && rowvals + group.size() > ncols) {
            s->end_row(rowvals);
            rowvals = 0;
          }
          if (rowvals == 0 && group.size() > ncols && group.size() < 2 * ncols) {
            for (size_t k = 0; k < group.size(); k++)
              s->value(k, group[k].data(), group[k].size(), gnull[k] != 0);
            s->end_row(group.size());
          } else {
            for (size_t k = 0; k < group.size(); k++) {
              s->value(rowvals, group[k].data(), group[k].size(), gnull[k] != 0);
              if (++rowvals == ncols) { s->end_row(ncols); rowvals = 0; }
            }
          }
          group.clear();
          gnull.clear();
        }
      }
      if (done) break;
      group.push_back(std::string(t.p, t.n));
      gnull.push_back(token_is_null(t) ? 1 : 0);
      t = get();
    }
    if (rowvals > 0) s->end_row(rowvals);
    unget(t);
  }
};

/* ========================================================================
 * Plugin handle and post-processing
 * ======================================================================== */

struct PdbxHandle {
  int natoms;
  int optflags;
  ChunkedArray<molfile_atom_t> atoms;
  std::vector<ChunkedArray<float> > frames;
  size_t next_frame;
  bool has_cell;
  float cell[6];
  std::vector<int> bond_from, bond_to;
  std::vector<float> bond_order;
  bool bonds_complete;        /* full bond list (chem_comp), not just special */
  bool have_orders;
  std::vector<molfile_graphics_t> graphics;

  PdbxHandle() : natoms(0), optflags(0), next_frame(0), has_cell(false),
                 bonds_complete(false), have_orders(false) {
    cell[0] = cell[1] = cell[2] = 0.0f;
    cell[3] = cell[4] = cell[5] = 90.0f;
  }
};

/* the items read from each table; nothing else is stored */
const char *const kCellItems[] = {
  "length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma", NULL
};
const char *const kStructConnItems[] = {
  "conn_type_id", "pdbx_value_order", "ptnr1_symmetry", "ptnr2_symmetry",
  "ptnr1_auth_asym_id", "ptnr1_auth_seq_id", "pdbx_ptnr1_pdb_ins_code",
  "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id",
  "pdbx_ptnr1_label_alt_id", "ptnr1_label_comp_id",
  "ptnr2_auth_asym_id", "ptnr2_auth_seq_id", "pdbx_ptnr2_pdb_ins_code",
  "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id",
  "pdbx_ptnr2_label_alt_id", "ptnr2_label_comp_id", NULL
};
const char *const kValidateItems[] = {
  "pdb_model_num",
  "auth_asym_id_1", "auth_seq_id_1", "pdb_ins_code_1", "auth_atom_id_1",
  "label_alt_id_1", "auth_comp_id_1",
  "auth_asym_id_2", "auth_seq_id_2", "pdb_ins_code_2", "auth_atom_id_2",
  "label_alt_id_2", "auth_comp_id_2",
  "auth_asym_id_3", "auth_seq_id_3", "pdb_ins_code_3", "auth_atom_id_3",
  "label_alt_id_3", "auth_comp_id_3", NULL
};
const char *const kChemCompAtomItems[] = {
  "comp_id", "atom_id", "type_symbol", "charge", "model_cartn_x", "model_cartn_y",
  "model_cartn_z", "pdbx_model_cartn_x_ideal", "pdbx_model_cartn_y_ideal",
  "pdbx_model_cartn_z_ideal", NULL
};
const char *const kChemCompBondItems[] = {
  "comp_id", "atom_id_1", "atom_id_2", "value_order", NULL
};
const char *const kSphereItems[] = {
  "cartn_x", "cartn_y", "cartn_z", "object_radius", "model_id", "asym_id", NULL
};

struct ParseContext {
  AtomSink atom_site;
  AtomSink ihm_start;
  TableSink cell, struct_conn, val_bond, val_angle, cc_atom, cc_bond, spheres;
  explicit ParseContext(bool use_label)
    : atom_site(kAtomSiteCols, use_label), ihm_start(kIhmStartCols, use_label),
      cell(kCellItems, 1), struct_conn(kStructConnItems), val_bond(kValidateItems),
      val_angle(kValidateItems), cc_atom(kChemCompAtomItems),
      cc_bond(kChemCompBondItems), spheres(kSphereItems) {}
};

void read_cell(const Table &t, PdbxHandle *h) {
  if (t.nrows == 0) return;
  float a = 0, b = 0, c = 0, al = 90, be = 90, ga = 90;
  parse_float(t.get(0, t.col("length_a")), &a);
  parse_float(t.get(0, t.col("length_b")), &b);
  parse_float(t.get(0, t.col("length_c")), &c);
  parse_float(t.get(0, t.col("angle_alpha")), &al);
  parse_float(t.get(0, t.col("angle_beta")), &be);
  parse_float(t.get(0, t.col("angle_gamma")), &ga);
  /* NMR and EM entries carry a 1 x 1 x 1 placeholder meaning "no cell" */
  if (a == 1.0f && b == 1.0f && c == 1.0f && al == 90.0f && be == 90.0f && ga == 90.0f)
    return;
  if (a > 0 && b > 0 && c > 0 && al > 0 && al < 180 && be > 0 && be < 180 &&
      ga > 0 && ga < 180) {
    h->has_cell = true;
    h->cell[0] = a; h->cell[1] = b; h->cell[2] = c;
    h->cell[3] = al; h->cell[4] = be; h->cell[5] = ga;
  }
}

/* ---- explicit bonds from struct_conn / validation tables ---- */

struct Partner {
  bool auth;            /* true: auth_asym/auth_seq/ins; false: label ids */
  std::string asym;
  int seq;              /* kSeqMissing = any residue of that asym */
  char ins;             /* 0 = none */
  std::string atom;
  char alt;             /* 0 = unspecified */
  std::string comp;     /* residue name, used to break ties; may be empty */
};

struct BondRequest {
  Partner a, b;
  float order;
};

bool make_partner(const Table &t, size_t r, int c_aasym, int c_aseq, int c_ins,
                  int c_lasym, int c_lseq, int c_atom, int c_alt, Partner *p) {
  const std::string &atom = t.get(r, c_atom);
  if (atom.empty()) return false;
  /* atom names are stored truncated to the molfile field size */
  p->atom = atom.substr(0, sizeof(((molfile_atom_t *)0)->name) - 1);
  p->comp.clear();
  const std::string &alt = t.get(r, c_alt);
  p->alt = alt.empty() ? 0 : alt[0];
  const std::string &ins = t.get(r, c_ins);
  p->ins = (ins.empty() || ins[0] == ' ') ? 0 : ins[0];
  int seq, lseq;
  char seq_ins = 0;
  const std::string &aasym = t.get(r, c_aasym);
  const std::string &lasym = t.get(r, c_lasym);
  const bool label_exact = !lasym.empty() && parse_int(t.get(r, c_lseq), &lseq);
  /* without an insertion-code column the author number is ambiguous (52 or
   * 52A?), while the label number is not */
  if (!aasym.empty() && parse_seq(t.get(r, c_aseq), &seq, &seq_ins) &&
      (c_ins >= 0 || seq_ins || !label_exact)) {
    p->auth = true;
    p->asym = aasym;
    p->seq = seq;
    if (!p->ins) p->ins = seq_ins;
    return true;
  }
  if (!lasym.empty()) {
    p->auth = false;
    p->asym = lasym;
    p->seq = parse_int(t.get(r, c_lseq), &seq) ? seq : kSeqMissing;
    p->ins = 0;
    return true;
  }
  return false;
}

void collect_struct_conn(const Table &t, std::vector<BondRequest> &out, long *skipped_symm) {
  if (t.nrows == 0) return;
  int c_type = t.col("conn_type_id");
  int c_order = t.col("pdbx_value_order");
  int c_sym1 = t.col("ptnr1_symmetry"), c_sym2 = t.col("ptnr2_symmetry");
  for (size_t r = 0; r < t.nrows; r++) {
    /* the connection types the PDB format writes as CONECT records; hydrogen
     * bonds, salt bridges, mismatches and modres entries are not bonds */
    const std::string &ty = t.get(r, c_type);
    if (ty.empty() || (!iequals(ty, "covale") && !iequals(ty, "covale_base") &&
        !iequals(ty, "covale_phosph") && !iequals(ty, "covale_sugar") &&
        !iequals(ty, "disulf") && !iequals(ty, "metalc")))
      continue;
    const std::string &s1 = t.get(r, c_sym1), &s2 = t.get(r, c_sym2);
    if ((!s1.empty() && s1 != "1_555") || (!s2.empty() && s2 != "1_555")) {
      (*skipped_symm)++;     /* partner is a symmetry mate, not in the file */
      continue;
    }
    BondRequest b;
    if (!make_partner(t, r, t.col("ptnr1_auth_asym_id"), t.col("ptnr1_auth_seq_id"),
                      t.col("pdbx_ptnr1_pdb_ins_code"), t.col("ptnr1_label_asym_id"),
                      t.col("ptnr1_label_seq_id"), t.col("ptnr1_label_atom_id"),
                      t.col("pdbx_ptnr1_label_alt_id"), &b.a))
      continue;
    b.a.comp = t.get(r, t.col("ptnr1_label_comp_id"));
    if (!make_partner(t, r, t.col("ptnr2_auth_asym_id"), t.col("ptnr2_auth_seq_id"),
                      t.col("pdbx_ptnr2_pdb_ins_code"), t.col("ptnr2_label_asym_id"),
                      t.col("ptnr2_label_seq_id"), t.col("ptnr2_label_atom_id"),
                      t.col("pdbx_ptnr2_label_alt_id"), &b.b))
      continue;
    b.b.comp = t.get(r, t.col("ptnr2_label_comp_id"));
    b.order = 1.0f;
    const std::string &o = t.get(r, c_order);
    if (iequals(o, "doub")) b.order = 2.0f;
    else if (iequals(o, "trip")) b.order = 3.0f;
    else if (iequals(o, "quad")) b.order = 4.0f;
    out.push_back(b);
  }
}

/* _pdbx_validate_rmsd_bond (atoms 1-2) and _pdbx_validate_rmsd_angle (1-2 and
 * 2-3) list covalent bonds with unusual geometry, which VMD's distance-based
 * bond search can miss. */
void collect_validate(const Table &t, int natoms_in_angle,
                      std::vector<BondRequest> &out, const std::string &model0) {
  if (t.nrows == 0) return;
  int c_model = t.col("pdb_model_num");
  int ca[3], cs[3], ci[3], cn[3], cl[3], cc[3];
  for (int k = 0; k < natoms_in_angle; k++) {
    char buf[64];
    snprintf(buf, sizeof(buf), "auth_asym_id_%d", k + 1); ca[k] = t.col(buf);
    snprintf(buf, sizeof(buf), "auth_seq_id_%d", k + 1); cs[k] = t.col(buf);
    snprintf(buf, sizeof(buf), "pdb_ins_code_%d", k + 1); ci[k] = t.col(buf);
    snprintf(buf, sizeof(buf), "auth_atom_id_%d", k + 1); cn[k] = t.col(buf);
    snprintf(buf, sizeof(buf), "label_alt_id_%d", k + 1); cl[k] = t.col(buf);
    snprintf(buf, sizeof(buf), "auth_comp_id_%d", k + 1); cc[k] = t.col(buf);
  }
  for (size_t r = 0; r < t.nrows; r++) {
    const std::string &m = t.get(r, c_model);
    if (!m.empty() && !model0.empty() && m != model0) continue;
    Partner p[3];
    bool ok = true;
    for (int k = 0; k < natoms_in_angle && ok; k++) {
      ok = make_partner(t, r, ca[k], cs[k], ci[k], -1, -1, cn[k], cl[k], &p[k]);
      if (ok) p[k].comp = t.get(r, cc[k]);
    }
    if (!ok) continue;
    for (int k = 0; k + 1 < natoms_in_angle; k++) {
      BondRequest b;
      b.a = p[k];
      b.b = p[k + 1];
      b.order = 1.0f;
      out.push_back(b);
    }
  }
}

std::string partner_key(bool auth, const std::string &asym, int seq, char ins,
                        const char *atom) {
  std::string k;
  k.reserve(asym.size() + strlen(atom) + 24);
  k += auth ? 'a' : 'l';
  k += '\x1f'; k += asym;
  k += '\x1f';
  if (seq == kSeqMissing) k += '*';
  else { char b[16]; snprintf(b, sizeof(b), "%d", seq); k += b; }
  k += '\x1f'; if (ins) k += ins;
  k += '\x1f'; k += atom;
  return k;
}

/* Two residue types can share one position (microheterogeneity, told apart
 * by altloc). When the bond record names the residue, prefer atoms of that
 * residue; if none match (label/auth naming differs), keep them all. */
std::vector<int> by_residue_name(const AtomSink &as, const std::vector<int> &c,
                                 const std::string &comp) {
  if (comp.empty() || c.size() < 2) return c;
  std::vector<int> out;
  for (size_t k = 0; k < c.size(); k++)
    if (strncmp(as.atoms[(size_t)c[k]].resname, comp.c_str(),
                sizeof(as.atoms[0].resname) - 1) == 0)
      out.push_back(c[k]);
  return out.empty() ? c : out;
}

/* With a chain selection, a record naming a chain of which no atom was
 * loaded is counted as outside the selection; every other record that does
 * not match is unmatched, as without a selection. */
void resolve_bonds(const AtomSink &as, const std::vector<BondRequest> &req,
                   PdbxHandle *h, long *unresolved, bool filtered, long *outside) {
  if (req.empty()) return;
  long *no_chain = filtered ? outside : unresolved;  /* a partner's chain not loaded */
  const size_t natoms = as.atoms.size();

  /* lookup keys needed by the requests, and which chains are involved */
  std::unordered_map<std::string, int> slots;
  std::vector<std::vector<int> > cand;
  std::vector<char> wild;          /* slot names a chain but no residue number */
  std::vector<char> rel_auth(as.asyms.size(), 0), rel_label(as.asyms.size(), 0);
  bool any_label_wild = false;
  std::vector<std::pair<int, int> > req_slots(req.size(), std::make_pair(-1, -1));
  for (size_t i = 0; i < req.size(); i++) {
    const Partner *ps[2] = { &req[i].a, &req[i].b };
    int sl[2] = { -1, -1 };
    for (int k = 0; k < 2; k++) {
      const Partner &p = *ps[k];
      int aid = as.asyms.find(p.asym);
      if (aid < 0) continue;
      if (p.auth) rel_auth[(size_t)aid] = 1;
      else { rel_label[(size_t)aid] = 1; if (p.seq == kSeqMissing) any_label_wild = true; }
      std::string key = partner_key(p.auth, p.asym, p.seq, p.ins, p.atom.c_str());
      std::unordered_map<std::string, int>::iterator it = slots.find(key);
      if (it == slots.end()) {
        it = slots.insert(std::make_pair(key, (int)cand.size())).first;
        cand.push_back(std::vector<int>());
        wild.push_back((!p.auth && p.seq == kSeqMissing) ? 1 : 0);
      }
      sl[k] = it->second;
    }
    req_slots[i] = std::make_pair(sl[0], sl[1]);
  }
  if (cand.empty()) { *no_chain += (long)req.size(); return; }

  /* one pass over the atoms of the involved chains. Atoms sharing a key are
   * alternate locations (rarely more than 4); the cap bounds the pairing work
   * a crafted file can cause. */
  const size_t kMaxCand = 8;
  for (size_t i = 0; i < natoms; i++) {
    const molfile_atom_t &a = as.atoms[i];
    char ins = as.insertion[i];
    int aa = as.auth_asym[i];
    if (aa >= 0 && rel_auth[(size_t)aa] && as.auth_seq[i] != kSeqMissing) {
      std::unordered_map<std::string, int>::const_iterator it =
          slots.find(partner_key(true, as.asyms.str(aa), as.auth_seq[i], ins, a.name));
      if (it != slots.end() && cand[(size_t)it->second].size() < kMaxCand)
        cand[(size_t)it->second].push_back((int)i);
    }
    int la = as.label_asym[i];
    if (la >= 0 && rel_label[(size_t)la]) {
      if (as.label_seq[i] != kSeqMissing) {
        std::unordered_map<std::string, int>::const_iterator it =
            slots.find(partner_key(false, as.asyms.str(la), as.label_seq[i], 0, a.name));
        if (it != slots.end() && cand[(size_t)it->second].size() < kMaxCand)
          cand[(size_t)it->second].push_back((int)i);
      }
      if (any_label_wild) {
        std::unordered_map<std::string, int>::const_iterator it =
            slots.find(partner_key(false, as.asyms.str(la), kSeqMissing, 0, a.name));
        if (it != slots.end() && cand[(size_t)it->second].size() < kMaxCand)
          cand[(size_t)it->second].push_back((int)i);
      }
    }
  }

  /* a partner given only as "this atom name somewhere in this chain" (e.g. a
   * water O, label_seq_id '.') is usable only if one residue has it */
  for (size_t sl = 0; sl < cand.size(); sl++) {
    if (!wild[sl] || cand[sl].size() < 2) continue;
    const size_t i0 = (size_t)cand[sl][0];
    for (size_t k = 1; k < cand[sl].size(); k++) {
      const size_t i = (size_t)cand[sl][k];
      if (as.atoms[i].resid != as.atoms[i0].resid || as.insertion[i] != as.insertion[i0] ||
          as.auth_asym[i] != as.auth_asym[i0]) {
        cand[sl].clear();
        break;
      }
    }
  }

  /* pair candidates; alternate locations bond within the same conformer */
  std::set<std::pair<int, int> > seen;
  for (size_t i = 0; i < req.size(); i++) {
    int s1 = req_slots[i].first, s2 = req_slots[i].second;
    if (s1 < 0 || s2 < 0) { (*no_chain)++; continue; }
    bool made = false;
    const std::vector<int> c1 = by_residue_name(as, cand[(size_t)s1], req[i].a.comp);
    const std::vector<int> c2 = by_residue_name(as, cand[(size_t)s2], req[i].b.comp);
    for (size_t x = 0; x < c1.size(); x++) {
      char alt1 = as.atoms[(size_t)c1[x]].altloc[0];
      if (req[i].a.alt && alt1 && alt1 != req[i].a.alt) continue;
      for (size_t y = 0; y < c2.size(); y++) {
        char alt2 = as.atoms[(size_t)c2[y]].altloc[0];
        if (req[i].b.alt && alt2 && alt2 != req[i].b.alt) continue;
        if (alt1 && alt2 && alt1 != alt2) continue;
        int u = c1[x], v = c2[y];
        if (u == v) continue;
        std::pair<int, int> key = (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
        if (!seen.insert(key).second) { made = true; continue; }
        h->bond_from.push_back(u + 1);         /* molfile bonds are 1-based */
        h->bond_to.push_back(v + 1);
        h->bond_order.push_back(req[i].order);
        if (req[i].order != 1.0f) h->have_orders = true;
        made = true;
      }
    }
    if (!made) (*unresolved)++;
  }
}

/* ---- chemical component dictionary files ---- */

bool build_from_chem_comp(const Table &atoms, const Table &bonds, PdbxHandle *h,
                          long *n_unknown, std::map<std::string, long> *unknown_syms) {
  if (atoms.nrows == 0) return false;
  int c_comp = atoms.col("comp_id"), c_atom = atoms.col("atom_id");
  int c_type = atoms.col("type_symbol"), c_charge = atoms.col("charge");
  const char *ideal[3] = { "pdbx_model_cartn_x_ideal", "pdbx_model_cartn_y_ideal",
                           "pdbx_model_cartn_z_ideal" };
  const char *model[3] = { "model_cartn_x", "model_cartn_y", "model_cartn_z" };
  int ci[3], cm[3];
  for (int k = 0; k < 3; k++) { ci[k] = atoms.col(ideal[k]); cm[k] = atoms.col(model[k]); }

  size_t n = atoms.nrows;
  if (n >= (size_t)(INT_MAX / 4)) throw std::runtime_error("too many atoms");
  std::vector<float> xi(3 * n, 0.0f), xm(3 * n, 0.0f);
  size_t ok_i = 0, ok_m = 0;
  bool any_charge = false, any_type = false;
  std::map<std::pair<std::string, std::string>, int> index;

  for (size_t r = 0; r < n; r++) {
    molfile_atom_t a;
    memset(&a, 0, sizeof(a));
    const std::string &name = atoms.get(r, c_atom);
    std::string sym = atoms.get(r, c_type);
    while (!sym.empty() && (sym[0] == ' ' || sym[0] == '\t')) sym.erase(0, 1);
    while (!sym.empty() && (sym[sym.size() - 1] == ' ' || sym[sym.size() - 1] == '\t'))
      sym.erase(sym.size() - 1);
    copy_field(a.name, sizeof(a.name), name);
    copy_field(a.type, sizeof(a.type), sym.empty() ? name : sym);
    copy_field(a.resname, sizeof(a.resname), atoms.get(r, c_comp));
    a.resid = 1;
    a.insertion[0] = ' ';
    a.occupancy = 1.0f;
    if (parse_charge(atoms.get(r, c_charge), &a.charge)) any_charge = true;
    int zn = element_from_symbol(sym);
    if (zn >= 0) any_type = true;
    if (zn > 0) a.atomicnumber = zn;
    else if (zn == 0) {
      (*n_unknown)++;
      if (unknown_syms->size() < 32 || unknown_syms->count(sym)) (*unknown_syms)[sym]++;
    }
    bool gi = true, gm = true;
    for (int k = 0; k < 3; k++) {
      gi = parse_float(atoms.get(r, ci[k]), &xi[3 * r + k]) && gi;
      gm = parse_float(atoms.get(r, cm[k]), &xm[3 * r + k]) && gm;
    }
    if (gi) ok_i++;
    if (gm) ok_m++;
    index[std::make_pair(atoms.get(r, c_comp), name)] = (int)r;
    h->atoms.push_back(a);
  }

  /* prefer the complete ideal coordinates, then the model coordinates */
  const std::vector<float> *use = &xi;
  if (ok_i == n) {
    use = &xi;
  } else if (ok_m == n) {
    use = &xm;
  } else {
    printf("pdbxplugin) Warning: chemical component has incomplete coordinates; "
           "missing values set to 0\n");
    use = (ok_i >= ok_m) ? &xi : &xm;
  }
  h->frames.push_back(ChunkedArray<float>());
  for (size_t k = 0; k < use->size(); k++) h->frames.back().push_back((*use)[k]);

  int b_comp = bonds.col("comp_id"), b_a1 = bonds.col("atom_id_1");
  int b_a2 = bonds.col("atom_id_2"), b_ord = bonds.col("value_order");
  for (size_t r = 0; r < bonds.nrows; r++) {
    const std::string &comp = bonds.get(r, b_comp);
    std::map<std::pair<std::string, std::string>, int>::const_iterator i1 =
        index.find(std::make_pair(comp, bonds.get(r, b_a1)));
    std::map<std::pair<std::string, std::string>, int>::const_iterator i2 =
        index.find(std::make_pair(comp, bonds.get(r, b_a2)));
    if (i1 == index.end() || i2 == index.end() || i1->second == i2->second) continue;
    const std::string &o = bonds.get(r, b_ord);
    float order = 1.0f;
    if (iequals(o, "doub")) order = 2.0f;
    else if (iequals(o, "trip")) order = 3.0f;
    else if (iequals(o, "quad")) order = 4.0f;
    else if (iequals(o, "arom") || iequals(o, "delo")) order = 1.5f;
    h->bond_from.push_back(i1->second + 1);
    h->bond_to.push_back(i2->second + 1);
    h->bond_order.push_back(order);
  }
  h->bonds_complete = true;
  h->have_orders = !h->bond_from.empty();

  h->optflags = MOLFILE_OCCUPANCY;
  if (any_type) h->optflags |= MOLFILE_ATOMICNUMBER;
  if (any_charge) h->optflags |= MOLFILE_CHARGE;
  return true;
}

/* ---- IHM spheres as raw graphics (same behavior as the original plugin) ---- */

/* with a chain selection, spheres are matched by their asym_id */
void build_spheres(const Table &t, PdbxHandle *h, ChainFilter *filter, long *filtered) {
  if (t.nrows == 0) return;
  int cx = t.col("cartn_x"), cy = t.col("cartn_y"), cz = t.col("cartn_z");
  int cr = t.col("object_radius"), cmod = t.col("model_id"), casym = t.col("asym_id");
  const std::string first_model = t.get(0, cmod);
  for (size_t r = 0; r < t.nrows; r++) {
    if (t.get(r, cmod) != first_model) continue;
    if (filter && filter->active() && !filter->accepts(t.get(r, casym))) {
      (*filtered)++;
      continue;
    }
    molfile_graphics_t g;
    memset(&g, 0, sizeof(g));
    float rad = 1.0f;
    if (!parse_float(t.get(r, cx), &g.data[0]) || !parse_float(t.get(r, cy), &g.data[1]) ||
        !parse_float(t.get(r, cz), &g.data[2]))
      continue;
    if (!parse_float(t.get(r, cr), &rad) || !(rad > 0.0f)) rad = 1.0f;
    g.type = MOLFILE_SPHERE;
    g.style = 12;       /* sphere tessellation resolution */
    g.size = rad;
    h->graphics.push_back(g);
  }
}

/* ---- diagnostics ---- */

void report_recovered_loops(const Parser &ps, const AtomSink *as, const char *atom_cat) {
  long atom_line = -1;
  for (size_t i = 0; i < ps.recovered.size(); i++) {
    const std::string &cat = ps.recovered[i].first;
    if (as && cat == atom_cat) { atom_line = ps.recovered[i].second; continue; }
    printf("pdbxplugin) Warning: _%s has rows with missing or extra values (loop starts "
           "at line %ld); each row was read from its own line\n",
           printable(cat).c_str(), ps.recovered[i].second);
  }
  if (!as || (as->fixed_rows == 0 && as->dropped_rows == 0)) return;
  printf("pdbxplugin) Warning: %ld _%s row(s) have missing or extra values",
         as->fixed_rows + as->dropped_rows, atom_cat);
  if (atom_line > 0) printf(" (loop starts at line %ld)", atom_line);
  printf(".\n");
  if (as->fixed_rows > 0)
    printf("pdbxplugin)   %ld of them were matched to their columns by value type; "
           "missing values are unknown.\n", as->fixed_rows);
  if (as->dropped_rows > 0)
    printf("pdbxplugin)   %ld could not be matched with confidence and were skipped.\n",
           as->dropped_rows);
}

void report_atom_warnings(const AtomSink &as, const char *cat) {
  if (as.bad_coord_rows > 0)
    printf("pdbxplugin) Warning: skipped %ld _%s row(s) without numeric coordinates\n",
           as.bad_coord_rows, cat);
  if (as.n_unknown_elem > 0) {
    std::string list;
    int shown = 0;
    for (std::map<std::string, long>::const_iterator it = as.unknown_symbols.begin();
         it != as.unknown_symbols.end(); ++it) {
      if (shown++ == 10) { list += ", ..."; break; }
      char b[96];
      snprintf(b, sizeof(b), "%s%s (%ld)", list.empty() ? "" : ", ",
               printable(it->first, 20).c_str(), it->second);
      list += b;
    }
    printf("pdbxplugin) Note: %ld atom(s) have element symbols that are not in VMD's "
           "periodic table: %s.\n"
           "pdbxplugin)   They are loaded as element X; the symbol is kept in the atom "
           "type and mass/radius are guessed from the atom name.\n",
           as.n_unknown_elem, list.c_str());
  }
  if (as.multichar_chain > 0)
    printf("pdbxplugin) Note: chain IDs longer than one character are present; "
           "'chain' holds the first character and 'segname' the full ID\n");
  if (as.long_chain > 0)
    printf("pdbxplugin) Warning: %ld atom(s) have chain IDs longer than 7 characters; "
           "segname was truncated\n", as.long_chain);
  if (as.long_name > 0 || as.long_resname > 0)
    printf("pdbxplugin) Warning: truncated %ld atom name(s) and %ld residue name(s) "
           "that exceed the plugin field sizes\n", as.long_name, as.long_resname);
  if (as.missing_resid > 0)
    printf("pdbxplugin) Warning: %ld atom(s) have no numeric residue number; set to 0\n",
           as.missing_resid);
}

/* ========================================================================
 * Reader callbacks
 * ======================================================================== */

/* "file.cif?chain=A,B": options after the last '?' of a name that is not
 * itself an existing file. Without a chain option, MMCIF_CHAIN applies.
 * Returns false with a message for an option it does not know. */
bool split_options(const char *filepath, std::string *path, ChainFilter *filter,
                   std::string *err) {
  *path = filepath;
  const char *q = strrchr(filepath, '?');
  struct stat st;
  /* Options are "key=value" pairs. A known key starts them, whatever its
   * value holds (a chain ID may contain '/'); any other "word=" is an
   * unknown option only if no path separator follows, since otherwise the
   * '?' is part of a directory name. */
  if (q) {
    size_t k = 0;
    while (is_alpha(q[1 + k]) || is_digit(q[1 + k]) || q[1 + k] == '_') k++;
    const std::string key = lowercase(q + 1, k);
    const bool known = q[1 + k] == '=' && (key == "chain" || key == "chains");
    const bool option_like = k > 0 && q[1 + k] == '=' && !strchr(q, '/') && !strchr(q, '\\');
    if (!known && !option_like) q = NULL;
  }
  if (q && stat(filepath, &st) != 0) {
    path->assign(filepath, (size_t)(q - filepath));
    const std::string opts(q + 1);
    size_t i = 0;
    while (i < opts.size()) {
      size_t amp = opts.find('&', i);
      if (amp == std::string::npos) amp = opts.size();
      const std::string kv = opts.substr(i, amp - i);
      i = amp + 1;
      if (kv.empty()) continue;
      size_t eq = kv.find('=');
      const std::string key = lowercase(kv.data(), eq == std::string::npos ? kv.size() : eq);
      if ((key == "chain" || key == "chains") && eq != std::string::npos) {
        filter->add(kv.substr(eq + 1));
        filter->source = "the file name";
      } else {
        *err = "unknown option '" + printable(kv) + "' after '?' in the file name "
               "(supported: ?chain=A,B)";
        return false;
      }
    }
    if (!filter->active()) {
      *err = "no chain IDs after '?chain=' in the file name";
      return false;
    }
  }
  if (!filter->active()) {
    const char *env = getenv("MMCIF_CHAIN");
    if (env && *env) {
      filter->add(env);
      if (filter->active()) filter->source = "MMCIF_CHAIN";
    }
  }
  return true;
}

void *open_pdbx_read(const char *filepath_arg, const char *filetype, int *natoms) {
  (void)filetype;
  if (!filepath_arg || !natoms) return NULL;
  PdbxHandle *h = NULL;
  std::string path_buf;
  const char *filepath = filepath_arg;           /* for messages until split */
  try {
    ChainFilter filter;
    std::string err;
    if (!split_options(filepath_arg, &path_buf, &filter, &err)) {
      printf("pdbxplugin) %s: %s\n", printable(filepath_arg, 200).c_str(), err.c_str());
      return NULL;
    }
    filepath = path_buf.c_str();
    CifReader reader;
    if (!reader.open(filepath)) {
      printf("pdbxplugin) Unable to open file %s: %s\n", filepath, strerror(errno));
      return NULL;
    }
    const bool binary = reader.binarycif();
    if (const char *fmt = binary ? NULL : reader.binary_format()) {
      const char *hint = "";
      if (strstr(fmt, "compressed")) hint = " (decompress it first)";
#ifndef _USE_ZLIB
      if (strstr(fmt, "gzip")) hint = " (decompress it first: this build has no zlib)";
#endif
      if (reader.compressed())
        printf("pdbxplugin) %s: the decompressed data is %s; this plugin reads text "
               "mmCIF files\n", filepath, fmt);
      else
        printf("pdbxplugin) %s is %s; this plugin reads text mmCIF files%s\n", filepath,
               fmt, hint);
      return NULL;
    }
    h = new PdbxHandle();
    const bool use_label = (getenv("VMDPDBXUSELABELIDS") != NULL);
    ParseContext ctx(use_label);
    ctx.atom_site.set_filter(&filter);
    ctx.ihm_start.set_filter(&filter);

    Parser parser(reader);
    parser.add_sink("atom_site", &ctx.atom_site);
    parser.add_sink("ihm_starting_model_coord", &ctx.ihm_start);
    parser.add_sink("cell", &ctx.cell);
    parser.add_sink("struct_conn", &ctx.struct_conn);
    parser.add_sink("pdbx_validate_rmsd_bond", &ctx.val_bond);
    parser.add_sink("pdbx_validate_rmsd_angle", &ctx.val_angle);
    parser.add_sink("chem_comp_atom", &ctx.cc_atom);
    parser.add_sink("chem_comp_bond", &ctx.cc_bond);
    parser.add_sink("ihm_sphere_obj_site", &ctx.spheres);

    const auto has_content = [&ctx]() {
      return ctx.atom_site.nrows() > 0 || ctx.ihm_start.nrows() > 0 ||
             ctx.cc_atom.nrows() > 0 || ctx.spheres.nrows() > 0;
    };
    if (binary) {
      std::vector<char> data;                    /* freed once the sinks hold the values */
      reader.read_all(data, kMaxBcifBytes);
      MpDoc doc;
      doc.parse(data.empty() ? "" : &data[0], data.size());
      parser.run_binary(doc, data.size(), has_content);
    } else {
      parser.run(has_content);
    }

    if (reader.unterminated_quotes > 0 || reader.unterminated_text > 0)
      printf("pdbxplugin) Warning: %ld unterminated quoted value(s) and %ld unterminated "
             "text field(s)\n", reader.unterminated_quotes, reader.unterminated_text);
    if (parser.stray_values > 0)
      printf("pdbxplugin) Warning: ignored %ld value(s) that belong to no tag or loop\n",
             parser.stray_values);
    if (parser.extra_blocks)
      printf("pdbxplugin) Note: file has more than one data block; only the first "
             "one with coordinates was read\n");

    read_cell(ctx.cell.table, h);
    long filtered_spheres = 0;
    build_spheres(ctx.spheres.table, h, &filter, &filtered_spheres);

    AtomSink *as = NULL;
    const char *cat = "atom_site";
    if (ctx.atom_site.nrows() > 0) as = &ctx.atom_site;
    else if (ctx.ihm_start.nrows() > 0) { as = &ctx.ihm_start; cat = "ihm_starting_model_coord"; }

    report_recovered_loops(parser, as, cat);

    long n_bonds_req = 0, unresolved = 0, skipped_symm = 0, outside = 0;
    if (as) {
      report_atom_warnings(*as, cat);

      /* frames: model 1, then every later model that matches it atom for atom */
      size_t n = as->atoms.size();
      h->frames.push_back(ChunkedArray<float>());
      h->frames.back().swap(as->xyz);
      size_t skipped = 0;
      for (size_t m = 0; m < as->extra.size(); m++) {
        ExtraModel &em = as->extra[m];
        if (em.matches && em.rows == n && em.xyz.size() == 3 * n) {
          h->frames.push_back(ChunkedArray<float>());
          h->frames.back().swap(em.xyz);
        } else if (em.rows == 0) {
          /* a model without rows of the selected chains: not a frame */
        } else if (++skipped <= 3) {
          printf("pdbxplugin) Warning: model %s (%lu atoms) does not match the %lu "
                 "atoms of model %s; skipped\n", printable(em.key).c_str(),
                 (unsigned long)em.rows, (unsigned long)n,
                 as->model0_key.empty() ? "1" : printable(as->model0_key).c_str());
        }
        em.xyz.clear();
      }
      if (skipped > 3)
        printf("pdbxplugin) Warning: %lu more model(s) skipped for the same reason; "
               "only models with the same atoms as the first load as frames\n",
               (unsigned long)(skipped - 3));
      if (as->dropped_model_rows > 0)
        printf("pdbxplugin) Warning: more than %lu models; %ld row(s) of later models "
               "were ignored\n", (unsigned long)kMaxModels, as->dropped_model_rows);
      as->ident.clear();

      std::vector<BondRequest> req;
      collect_struct_conn(ctx.struct_conn.table, req, &skipped_symm);
      collect_validate(ctx.val_bond.table, 2, req, as->model0_key);
      collect_validate(ctx.val_angle.table, 3, req, as->model0_key);
      n_bonds_req = (long)req.size();
      resolve_bonds(*as, req, h, &unresolved, filter.active(), &outside);

      int fl = MOLFILE_INSERTION;
      if (as->n_known_elem + as->n_unknown_elem > 0) fl |= MOLFILE_ATOMICNUMBER;
      if (as->has_col[AC_OCC]) fl |= MOLFILE_OCCUPANCY;
      if (as->has_col[AC_B]) fl |= MOLFILE_BFACTOR;
      if (as->any_charge) fl |= MOLFILE_CHARGE;
      if (as->any_altloc) fl |= MOLFILE_ALTLOC;
      if (!as->any_insertion) fl &= ~MOLFILE_INSERTION;
      h->optflags = fl;
      h->atoms.swap(as->atoms);
    } else {
      long n_unknown = 0;
      std::map<std::string, long> unknown_syms;
      if (build_from_chem_comp(ctx.cc_atom.table, ctx.cc_bond.table, h,
                               &n_unknown, &unknown_syms)) {
        if (n_unknown > 0)
          printf("pdbxplugin) Note: %ld atom(s) have element symbols that are not in "
                 "VMD's periodic table; loaded as element X\n", n_unknown);
        if (filter.active())
          printf("pdbxplugin) Note: the chain selection (%s, from %s) does not apply to "
                 "a chemical component file; all atoms were loaded\n",
                 filter.list().c_str(), filter.source.c_str());
      }
    }

    /* the molfile API counts in int */
    if (h->bond_from.size() > (size_t)INT_MAX || h->graphics.size() > (size_t)INT_MAX)
      throw std::runtime_error("too many bonds or graphics objects");
    h->natoms = (int)h->atoms.size();
    const bool filtered = filter.active() && !h->bonds_complete;
    const long n_filtered = (as ? as->filtered_rows : 0) + filtered_spheres;
    /* selected chains with atoms in the structure, or spheres (filter.hit
     * holds the sphere matches) */
    if (filtered && as)
      for (size_t i = 0; i < filter.hit.size() && i < as->found.size(); i++)
        if (as->found[i]) filter.hit[i] = 1;
    const bool multi_model = as && !as->extra.empty();
    if (h->natoms == 0 && h->graphics.empty()) {
      long bad = ctx.atom_site.bad_coord_rows + ctx.ihm_start.bad_coord_rows +
                 ctx.atom_site.dropped_rows + ctx.ihm_start.dropped_rows;
      if (filtered && n_filtered > 0) {
        std::string have;
        if (as) {
          for (size_t i = 0; i < as->chains_seen.size(); i++) {
            if (i == 20) { have += ", ..."; break; }
            if (!have.empty()) have += ", ";
            have += as->chains_seen[i].empty() ? std::string("(none)")
                                                : printable(as->chains_seen[i], 20);
          }
          if (as->more_chains && as->chains_seen.size() <= 20) have += ", ...";
        }
        printf("pdbxplugin) No atoms loaded from %s: chain(s) %s (selected by %s) have "
               "no atoms with coordinates in the %s%s%s\n", filepath, filter.list().c_str(),
               filter.source.c_str(), multi_model ? "first model" : "file",
               have.empty() ? "" : ", whose chains are ", have.c_str());
      } else if (bad > 0)
        printf("pdbxplugin) No atoms loaded from %s: none of the %ld atom row(s) has "
               "usable coordinates\n", filepath, bad);
      else
        printf("pdbxplugin) No atom coordinates found in %s "
               "(looked for _atom_site, _ihm_starting_model_coord, _chem_comp_atom)\n",
               filepath);
      delete h;
      return NULL;
    }
    if (h->natoms == 0) h->frames.clear();

    printf("pdbxplugin) Read %d atom(s), %lu model(s)", h->natoms,
           (unsigned long)h->frames.size());
    if (!h->bond_from.empty())
      printf(", %lu %s bond(s)", (unsigned long)h->bond_from.size(),
             h->bonds_complete ? "component" : "explicit");
    if (!h->graphics.empty())
      printf(", %lu IHM sphere(s) as graphics", (unsigned long)h->graphics.size());
    if (use_label) printf(" [label_asym_id/label_seq_id numbering]");
    if (binary) printf(" [BinaryCIF]");
    if (reader.compressed()) printf(" [gzip]");
    printf("\n");
    if (filtered) {
      printf("pdbxplugin) Note: loaded only chain(s) %s (selected by %s)",
             filter.list(1).c_str(), filter.source.c_str());
      if (n_filtered > 0) printf("; skipped %ld row(s) of other chains", n_filtered);
      printf("\n");
      std::string missing = filter.list(-1);
      if (!missing.empty())
        printf("pdbxplugin) Warning: chain(s) %s have no atoms with coordinates in the %s\n",
               missing.c_str(), multi_model ? "first model" : "file");
    }
    if (unresolved > 0)
      printf("pdbxplugin) Note: %ld of %ld explicit bond record(s) did not match %s\n",
             unresolved, n_bonds_req, filtered ? "the loaded atoms" : "atoms in the file");
    if (outside > 0)
      printf("pdbxplugin) Note: skipped %ld explicit bond record(s) to chains that were "
             "not loaded\n", outside);
    if (skipped_symm > 0)
      printf("pdbxplugin) Note: skipped %ld bond(s) to symmetry-related copies\n",
             skipped_symm);

    *natoms = h->natoms;
    return h;
  } catch (const std::bad_alloc &) {
    printf("pdbxplugin) Out of memory reading %s\n", filepath);
  } catch (const std::exception &e) {
    printf("pdbxplugin) Error reading %s: %s\n", filepath, e.what());
  } catch (...) {
    printf("pdbxplugin) Unknown error reading %s\n", filepath);
  }
  delete h;
  return NULL;
}

int read_pdbx_structure(void *v, int *optflags, molfile_atom_t *atoms) {
  PdbxHandle *h = (PdbxHandle *)v;
  *optflags = MOLFILE_NOOPTIONS;
  if (h->natoms == 0) return MOLFILE_NOSTRUCTUREDATA;
  if (h->atoms.size() != (size_t)h->natoms) {
    printf("pdbxplugin) read_structure called more than once\n");
    return MOLFILE_ERROR;
  }
  h->atoms.copy_to(atoms);
  /* VMD keeps its own copy; release ours (~80 bytes/atom) right away */
  h->atoms.clear();
  int fl = h->optflags;
  if (!h->bond_from.empty() && !h->bonds_complete) fl |= MOLFILE_BONDSSPECIAL;
  *optflags = fl;
  return MOLFILE_SUCCESS;
}

int read_pdbx_bonds(void *v, int *nbonds, int **fromptr, int **toptr,
                    float **bondorderptr, int **bondtype, int *nbondtypes,
                    char ***bondtypename) {
  PdbxHandle *h = (PdbxHandle *)v;
  *nbonds = (int)h->bond_from.size();
  *fromptr = h->bond_from.empty() ? NULL : &h->bond_from[0];
  *toptr = h->bond_to.empty() ? NULL : &h->bond_to[0];
  *bondorderptr = (h->have_orders && !h->bond_order.empty()) ? &h->bond_order[0] : NULL;
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;
  return MOLFILE_SUCCESS;
}

int read_pdbx_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  PdbxHandle *h = (PdbxHandle *)v;
  if (h->next_frame >= h->frames.size()) return MOLFILE_EOF;
  if (natoms != h->natoms) return MOLFILE_ERROR;
  if (ts) {
    const ChunkedArray<float> &f = h->frames[h->next_frame];
    if (f.size() != 3 * (size_t)natoms) return MOLFILE_ERROR;
    f.copy_to(ts->coords);
    ts->A = h->cell[0]; ts->B = h->cell[1]; ts->C = h->cell[2];
    ts->alpha = h->cell[3]; ts->beta = h->cell[4]; ts->gamma = h->cell[5];
  }
  /* release frames as they are consumed */
  h->frames[h->next_frame].clear();
  h->next_frame++;
  return MOLFILE_SUCCESS;
}

int read_pdbx_rawgraphics(void *v, int *nelem, const molfile_graphics_t **data) {
  PdbxHandle *h = (PdbxHandle *)v;
  *nelem = (int)h->graphics.size();
  *data = h->graphics.empty() ? NULL : &h->graphics[0];
  return MOLFILE_SUCCESS;
}

void close_pdbx_read(void *v) {
  delete (PdbxHandle *)v;
}

/* ========================================================================
 * Writer
 * ======================================================================== */

/* Writer output: buffered text, gzip-compressed for a name ending in .gz. */
class Out {
 public:
  Out() : fp_(NULL), gz_(false), error_(false) {
#ifdef _USE_ZLIB
    memset(&zs_, 0, sizeof(zs_));
    zinit_ = false;
#endif
  }
  ~Out() {
#ifdef _USE_ZLIB
    if (zinit_) deflateEnd(&zs_);
#endif
    if (fp_) fclose(fp_);
  }

  /* false with a message in *err */
  bool open(const char *path, bool gz, std::string *err) {
#ifndef _USE_ZLIB
    if (gz) { *err = "this build has no zlib, so it cannot write .gz files"; return false; }
#endif
    fp_ = fopen(path, "wb");
    if (!fp_) { *err = strerror(errno); return false; }
#ifdef _USE_ZLIB
    if (gz) {
      if (deflateInit2(&zs_, 6, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        *err = "cannot initialize zlib";
        return false;
      }
      zinit_ = true;
      zout_.resize((size_t)1 << 16);
    }
#endif
    gz_ = gz;
    return true;
  }

  void put(const char *s) { put(s, strlen(s)); }
  void put(const char *s, size_t n) {
    buf_.append(s, n);
    if (buf_.size() >= kFlushBytes) flush(false);
  }
  void put(char c) {
    buf_ += c;
    if (buf_.size() >= kFlushBytes) flush(false);
  }
  void format(const char *fmt, ...) {
    char b[256];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(b, sizeof(b), fmt, ap);
    va_end(ap);
    if (n < 0) { error_ = true; return; }
    if ((size_t)n < sizeof(b)) { put(b, (size_t)n); return; }
    std::vector<char> big((size_t)n + 1);
    va_start(ap, fmt);
    vsnprintf(&big[0], big.size(), fmt, ap);
    va_end(ap);
    put(&big[0], (size_t)n);
  }

  bool error() const { return error_; }

  /* flush everything and close; false if anything failed */
  bool close() {
    if (!fp_) return false;
    flush(true);
#ifdef _USE_ZLIB
    if (zinit_) { deflateEnd(&zs_); zinit_ = false; }
#endif
    if (fclose(fp_) != 0) error_ = true;
    fp_ = NULL;
    return !error_;
  }

 private:
  static const size_t kFlushBytes = (size_t)1 << 16;
  FILE *fp_;
  bool gz_, error_;
  std::string buf_;
#ifdef _USE_ZLIB
  z_stream zs_;
  bool zinit_;
  std::vector<unsigned char> zout_;
#endif

  void flush(bool final) {
    if (!fp_ || error_) { buf_.clear(); return; }
#ifdef _USE_ZLIB
    if (gz_) {
      zs_.next_in = (Bytef *)(buf_.empty() ? NULL : &buf_[0]);
      zs_.avail_in = (uInt)buf_.size();
      for (;;) {
        zs_.next_out = &zout_[0];
        zs_.avail_out = (uInt)zout_.size();
        int r = deflate(&zs_, final ? Z_FINISH : Z_NO_FLUSH);
        if (r == Z_STREAM_ERROR) { error_ = true; break; }
        size_t n = zout_.size() - zs_.avail_out;
        if (n && fwrite(&zout_[0], 1, n, fp_) != n) { error_ = true; break; }
        if (final ? (r == Z_STREAM_END) : (zs_.avail_in == 0 && zs_.avail_out > 0)) break;
      }
      buf_.clear();
      return;
    }
#else
    (void)final;
#endif
    if (!buf_.empty() && fwrite(buf_.data(), 1, buf_.size(), fp_) != buf_.size())
      error_ = true;
    buf_.clear();
  }
};

struct PdbxWriter {
  Out out;
  std::string block;
  int natoms;
  int optflags;
  std::vector<molfile_atom_t> atoms;
  int nframes;
  bool io_error;
  PdbxWriter() : natoms(0), optflags(0), nframes(0), io_error(false) {}
};

/* The text the writer stores for a field: without the padding molfile fields
 * may carry, and with characters CIF 1.1 does not allow (control characters,
 * line breaks, non-ASCII bytes) replaced by '_'. Empty means null. */
std::string cif_text(const char *s) {
  size_t n = strlen(s);
  while (n > 0 && s[n - 1] == ' ') n--;
  size_t b = 0;
  while (b < n && s[b] == ' ') b++;
  std::string v(s + b, n - b);
  for (size_t i = 0; i < v.size(); i++) {
    unsigned char c = (unsigned char)v[i];
    if (c < 32 || c >= 127) v[i] = '_';
  }
  return v;
}

/* whether quote character q can delimit v: inside a value, q may not be
 * followed by a blank, which would end the value */
bool quote_fits(const std::string &v, char q) {
  for (size_t i = 0; i + 1 < v.size(); i++)
    if (v[i] == q && v[i + 1] == ' ') return false;
  return true;
}

/* Write a value, quoting only when CIF requires it. */
void put_value(Out &o, const char *s, const char *null_tok) {
  std::string v = cif_text(s);
  if (v.empty()) { o.put(null_tok); return; }

  bool quote = false;
  char c0 = v[0];
  if (c0 == '_' || c0 == '#' || c0 == '$' || c0 == '\'' || c0 == '"' ||
      c0 == '[' || c0 == ']' || c0 == '{' || c0 == '}' || c0 == ';')
    quote = true;
  if (v == "." || v == "?") quote = true;
  if (iprefix(v.data(), v.size(), "data_") || iprefix(v.data(), v.size(), "save_") ||
      iequals(v, "loop_") || iequals(v, "stop_") || iequals(v, "global_"))
    quote = true;
  for (size_t i = 0; i < v.size() && !quote; i++)
    if (v[i] == ' ') quote = true;

  if (!quote) { o.put(v.c_str()); return; }
  if (quote_fits(v, '\'') || quote_fits(v, '"')) {
    char q = quote_fits(v, '\'') ? '\'' : '"';
    o.put(q);
    o.put(v.c_str());
    o.put(q);
    return;
  }
  /* both quote kinds occur before a blank: only a text field holds it */
  o.put("\n;");
  o.put(v.c_str());
  o.put("\n;");
}

/* Fixed-point formatting independent of LC_NUMERIC. */
void put_fixed(Out &o, double v, int decimals) {
  if (!std::isfinite(v)) { o.put("?"); return; }
  double scaled = v * kPow10[decimals];
  if (std::fabs(scaled) >= 9.0e18) {
    /* too large for fixed point: mantissa and exponent, still locale free */
    int e = (int)std::floor(std::log10(std::fabs(v)));
    double m = v / std::pow(10.0, (double)e);
    put_fixed(o, m, 6);
    o.format("e%d", e);
    return;
  }
  long long iv = std::llround(scaled);
  bool neg = iv < 0;
  unsigned long long u = (unsigned long long)(neg ? -iv : iv);
  unsigned long long p = (unsigned long long)kPow10[decimals];
  o.format("%s%llu.%0*llu", neg ? "-" : "", u / p, decimals, u % p);
}

bool is_polymer_residue(const char *resname) {
  static const char *const names[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "ASX", "GLX", "UNK",
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN",
    "A", "C", "G", "U", "I", "N", "DA", "DC", "DG", "DT", "DI", "DU", "DN",
    NULL
  };
  for (int i = 0; names[i]; i++)
    if (strcmp(resname, names[i]) == 0) return true;
  return false;
}

std::string writer_chain_id(const molfile_atom_t &a) {
  std::string chain(a.chain, strnlen(a.chain, sizeof(a.chain)));
  std::string seg(a.segid, strnlen(a.segid, sizeof(a.segid)));
  while (!seg.empty() && seg[seg.size() - 1] == ' ') seg.erase(seg.size() - 1);
  if (chain == " ") chain.clear();
  /* the reader stores the full chain ID in segid; use it when consistent */
  if (!seg.empty() && seg.find(' ') == std::string::npos &&
      (chain.empty() || seg[0] == chain[0]))
    return seg;
  return chain;                       /* may be empty: written as '?' */
}

std::string writer_element(const molfile_atom_t &a, int optflags) {
  if (a.atomicnumber > 0 && a.atomicnumber <= kMaxElement) {
    std::string s(kElements[a.atomicnumber]);
    for (size_t i = 0; i < s.size(); i++) s[i] = to_upper(s[i]);
    return s;
  }
  /* element given as unknown (0): the type may hold the file's symbol, e.g.
   * a coarse-grained bead type, so write that -- but only if it reads back as
   * unknown too. A force-field type or name such as "CA" would come back as
   * calcium. */
  if ((optflags & MOLFILE_ATOMICNUMBER) && a.type[0]) {
    std::string t(a.type, strnlen(a.type, sizeof(a.type)));
    if (element_from_symbol(cif_text(t.c_str())) == 0) return t;   /* as written */
  }
  return std::string();               /* no element information: '?' */
}

void *open_pdbx_write(const char *path, const char *filetype, int natoms) {
  (void)filetype;
  if (natoms < 1) {
    printf("pdbxplugin) Cannot write a file with %d atoms\n", natoms);
    return NULL;
  }
  PdbxWriter *w = NULL;
  try {
    w = new PdbxWriter();
    const size_t plen = strlen(path);
    const bool gz = plen > 3 && iequals(path + plen - 3, 3, ".gz");
    const size_t stem = gz ? plen - 3 : plen;
    if (stem >= 5 && iequals(path + stem - 5, 5, ".bcif")) {
      printf("pdbxplugin) Cannot write %s: this plugin writes text mmCIF only "
             "(use a name ending in .cif or .cif.gz)\n", path);
      delete w;
      return NULL;
    }
    std::string err;
    if (!w->out.open(path, gz, &err)) {
      printf("pdbxplugin) Unable to open %s for writing: %s\n", path, err.c_str());
      delete w;
      return NULL;
    }
    w->natoms = natoms;

    /* data block name from the file name, without directory and extension */
    const char *base = path;
    for (const char *p = path; *p; p++)
      if (*p == '/' || *p == '\\') base = p + 1;
    std::string name(base);
    if (gz) name.erase(name.size() - std::min(name.size(), (size_t)3));
    size_t dot = name.rfind('.');
    if (dot != std::string::npos && dot > 0) name.erase(dot);
    for (size_t i = 0; i < name.size(); i++)
      if (is_ws(name[i]) || (unsigned char)name[i] < 32) name[i] = '_';
    w->block = name.empty() ? std::string("vmd") : name;
    return w;
  } catch (...) {
    delete w;                 /* Out closes the file */
    return NULL;
  }
}

int write_pdbx_structure(void *v, int optflags, const molfile_atom_t *atoms) {
  PdbxWriter *w = (PdbxWriter *)v;
  try {
    w->atoms.assign(atoms, atoms + w->natoms);
  } catch (...) {
    return MOLFILE_ERROR;
  }
  w->optflags = optflags;
  return MOLFILE_SUCCESS;
}

void write_header(PdbxWriter *w, const molfile_timestep_t *ts) {
  Out &o = w->out;
  o.put("data_");
  o.put(w->block.c_str());
  o.format("\n# written by the VMD pdbx plugin %d.%d\n#\n",
          PDBX_MAJOR_VERSION, PDBX_MINOR_VERSION);
  if (ts && ts->A > 0 && ts->B > 0 && ts->C > 0) {
    o.put("_cell.length_a    "); put_fixed(o, ts->A, 3); o.put('\n');
    o.put("_cell.length_b    "); put_fixed(o, ts->B, 3); o.put('\n');
    o.put("_cell.length_c    "); put_fixed(o, ts->C, 3); o.put('\n');
    o.put("_cell.angle_alpha "); put_fixed(o, ts->alpha, 2); o.put('\n');
    o.put("_cell.angle_beta  "); put_fixed(o, ts->beta, 2); o.put('\n');
    o.put("_cell.angle_gamma "); put_fixed(o, ts->gamma, 2); o.put('\n');
    o.put("#\n_symmetry.space_group_name_H-M 'P 1'\n#\n");
  }
  static const char *const items[] = {
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
    "B_iso_or_equiv", "pdbx_formal_charge", "auth_seq_id", "auth_comp_id",
    "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num", NULL
  };
  o.put("loop_\n");
  for (int i = 0; items[i]; i++) o.format("_atom_site.%s\n", items[i]);
}

void write_atoms(PdbxWriter *w, const float *coords, int model) {
  Out &o = w->out;
  const int fl = w->optflags;
  for (int i = 0; i < w->natoms; i++) {
    const molfile_atom_t &a = w->atoms[(size_t)i];
    std::string resname(a.resname, strnlen(a.resname, sizeof(a.resname)));
    std::string name(a.name, strnlen(a.name, sizeof(a.name)));
    std::string chain = writer_chain_id(a);
    std::string elem = writer_element(a, fl);
    bool polymer = is_polymer_residue(resname.c_str());
    char alt[2] = { (fl & MOLFILE_ALTLOC) ? a.altloc[0] : '\0', '\0' };
    char ins[2] = { (fl & MOLFILE_INSERTION) ? a.insertion[0] : '\0', '\0' };
    char num[32];

    o.put(polymer ? "ATOM   " : "HETATM ");
    o.format("%lld ", (long long)(model - 1) * w->natoms + i + 1);
    put_value(o, elem.c_str(), "?");               o.put(' ');
    put_value(o, name.c_str(), "?");               o.put(' ');
    put_value(o, alt, ".");                        o.put(' ');
    put_value(o, resname.c_str(), "?");            o.put(' ');
    put_value(o, chain.c_str(), "?");              o.put(" ? ");
    if (polymer) o.format("%d ", a.resid); else o.put(". ");
    put_value(o, ins, "?");                        o.put(' ');
    if (coords) {
      put_fixed(o, coords[3 * i], 3);     o.put(' ');
      put_fixed(o, coords[3 * i + 1], 3); o.put(' ');
      put_fixed(o, coords[3 * i + 2], 3); o.put(' ');
    } else {
      o.put("0.000 0.000 0.000 ");
    }
    put_fixed(o, (fl & MOLFILE_OCCUPANCY) ? a.occupancy : 1.0, 2); o.put(' ');
    put_fixed(o, (fl & MOLFILE_BFACTOR) ? a.bfactor : 0.0, 2);     o.put(' ');
    /* only integral charges are formal charges */
    if ((fl & MOLFILE_CHARGE) && std::isfinite(a.charge) &&
        std::fabs(a.charge - std::floor(a.charge + 0.5f)) < 1e-3f &&
        std::fabs(a.charge) < 100.0f) {
      snprintf(num, sizeof(num), "%d ", (int)std::floor(a.charge + 0.5f));
      o.put(num);
    } else {
      o.put("? ");
    }
    o.format("%d ", a.resid);
    put_value(o, resname.c_str(), "?");            o.put(' ');
    put_value(o, chain.c_str(), "?");              o.put(' ');
    put_value(o, name.c_str(), "?");               o.put(' ');
    o.format("%d\n", model);
  }
  if (o.error()) w->io_error = true;
}

int write_pdbx_timestep(void *v, const molfile_timestep_t *ts) {
  PdbxWriter *w = (PdbxWriter *)v;
  if (w->atoms.size() != (size_t)w->natoms || !ts || !ts->coords)
    return MOLFILE_ERROR;
  try {                       /* no C++ exception may reach VMD */
    if (w->nframes == 0) write_header(w, ts);
    w->nframes++;
    write_atoms(w, ts->coords, w->nframes);
  } catch (...) {
    w->io_error = true;
  }
  return w->io_error ? MOLFILE_ERROR : MOLFILE_SUCCESS;
}

void close_pdbx_write(void *v) {
  PdbxWriter *w = (PdbxWriter *)v;
  if (w) {
    try {
      if (w->nframes == 0 && w->atoms.size() == (size_t)w->natoms) {
        /* structure without coordinates: keep the atom records anyway */
        write_header(w, NULL);
        write_atoms(w, NULL, 1);
      }
      if (w->nframes > 0 || !w->atoms.empty()) w->out.put("#\n");
    } catch (...) {
      w->io_error = true;
    }
    if (!w->out.close() || w->io_error)
      printf("pdbxplugin) Error while writing the mmCIF file\n");
  }
  delete w;
}

molfile_plugin_t plugin;
molfile_plugin_t mmcif_alias;    /* same reader and writer, for "type mmcif" */

} /* namespace */

/* ========================================================================
 * Registration
 * ======================================================================== */

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "pdbx";
  plugin.prettyname = "mmCIF/PDBx";
  plugin.author = "Diego Enry Barreto Gomes";
  plugin.majorv = PDBX_MAJOR_VERSION;
  plugin.minorv = PDBX_MINOR_VERSION;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "cif,mmcif,bcif";
  plugin.open_file_read = open_pdbx_read;
  plugin.read_structure = read_pdbx_structure;
  plugin.read_bonds = read_pdbx_bonds;
  plugin.read_next_timestep = read_pdbx_timestep;
  plugin.read_rawgraphics = read_pdbx_rawgraphics;
  plugin.close_file_read = close_pdbx_read;
  plugin.open_file_write = open_pdbx_write;
  plugin.write_structure = write_pdbx_structure;
  plugin.write_timestep = write_pdbx_timestep;
  plugin.close_file_write = close_pdbx_write;

  /* The repository's former "mmcif" reader is merged into this plugin. The
   * alias keeps "mol new ... type mmcif" working and, with its higher
   * version, shadows an old mmcifplugin.so that is still installed. */
  mmcif_alias = plugin;
  mmcif_alias.name = "mmcif";
  mmcif_alias.prettyname = "mmCIF (same as PDBx)";
  mmcif_alias.filename_extension = "mmcif";
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  (*cb)(v, (vmdplugin_t *)&mmcif_alias);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}
