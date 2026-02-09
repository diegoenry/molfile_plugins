/*
 * cifparse.h - Minimal CIF parser for mmCIF VMD plugin
 *
 * Reads a CIF file into in-memory tables that can be queried by category name.
 * Handles loop_ constructs, tag-value pairs, quoted strings, semicolon text
 * fields, and comments. Case-insensitive category/column matching.
 *
 * Author: Diego Enry Barreto Gomes
 */

#ifndef CIFPARSE_H
#define CIFPARSE_H

#include <string>
#include <vector>
#include <stdexcept>

namespace cif {

struct Table {
    std::string category;                         /* e.g. "atom_site" */
    std::vector<std::string> columns;             /* column names without prefix */
    std::vector<std::vector<std::string>> rows;

    /* Case-insensitive column lookup; returns -1 if not found */
    int find_column(const std::string &name) const;
};

struct DataBlock {
    std::string name;
    std::vector<Table> tables;

    /* Case-insensitive category lookup; returns nullptr if not found */
    const Table *find_table(const std::string &category) const;
};

/* Parse a CIF file and return the first data block.
 * Throws std::runtime_error on I/O or syntax errors. */
DataBlock parse_file(const char *filename);

} /* namespace cif */

#endif /* CIFPARSE_H */
