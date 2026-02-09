/*
 * cifparse.cpp - Minimal CIF parser implementation
 *
 * Tokenizes a CIF file character-by-character and builds in-memory tables.
 * Supports CIF 1.1 features needed for mmCIF: data_ blocks, loop_ constructs,
 * tag-value pairs, single/double quoted strings, semicolon text fields, and
 * # comments. Keywords are case-insensitive.
 *
 * Author: Diego Enry Barreto Gomes
 */

#include "cifparse.h"

#include <cstdio>
#include <cstring>
#include <cctype>

namespace cif {

/* ------------------------------------------------------------------
 * Table / DataBlock helpers
 * ------------------------------------------------------------------ */

static bool iequals(const std::string &a, const std::string &b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (tolower((unsigned char)a[i]) != tolower((unsigned char)b[i]))
            return false;
    return true;
}

int Table::find_column(const std::string &name) const {
    for (size_t i = 0; i < columns.size(); i++)
        if (iequals(columns[i], name))
            return (int)i;
    return -1;
}

const Table *DataBlock::find_table(const std::string &cat) const {
    for (size_t i = 0; i < tables.size(); i++)
        if (iequals(tables[i].category, cat))
            return &tables[i];
    return nullptr;
}

/* ------------------------------------------------------------------
 * Tokenizer
 * ------------------------------------------------------------------ */

enum TokenType { T_DATA, T_LOOP, T_TAG, T_VALUE, T_EOF };

struct Token {
    TokenType type;
    std::string text;   /* payload: block name, tag, or value */
};

/* CIF whitespace: HT(9), LF(10), CR(13), SPACE(32) */
static bool is_ws(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

class Tokenizer {
public:
    Tokenizer(const char *data, size_t len)
        : buf(data), end(data + len), pos(data), has_peeked(false) {}

    Token next() {
        if (has_peeked) {
            has_peeked = false;
            return peeked;
        }
        return read_token();
    }

    TokenType peek_type() {
        if (!has_peeked) {
            peeked = read_token();
            has_peeked = true;
        }
        return peeked.type;
    }

private:
    const char *buf;
    const char *end;
    const char *pos;
    bool has_peeked;
    Token peeked;

    void skip_ws_and_comments() {
        while (pos < end) {
            if (is_ws(*pos)) {
                pos++;
            } else if (*pos == '#') {
                /* skip to end of line */
                while (pos < end && *pos != '\n') pos++;
            } else {
                break;
            }
        }
    }

    Token read_token() {
        skip_ws_and_comments();
        if (pos >= end) return Token{T_EOF, ""};

        /* Semicolon text field: must be at start of line */
        if (*pos == ';' && (pos == buf || *(pos - 1) == '\n' || *(pos - 1) == '\r')) {
            pos++;  /* skip opening ';' */
            const char *start = pos;
            /* Find closing ';' at start of a line */
            while (pos < end) {
                if (*pos == ';' && (*(pos - 1) == '\n' || *(pos - 1) == '\r')) {
                    std::string val(start, pos - start);
                    /* Trim trailing newline(s) */
                    while (!val.empty() && (val.back() == '\n' || val.back() == '\r'))
                        val.pop_back();
                    pos++;  /* skip closing ';' */
                    return Token{T_VALUE, val};
                }
                pos++;
            }
            /* Unterminated - return what we have */
            return Token{T_VALUE, std::string(start, end - start)};
        }

        /* Quoted string */
        if (*pos == '\'' || *pos == '"') {
            char quote = *pos;
            pos++;
            const char *start = pos;
            /* Quote terminates at matching quote followed by whitespace or EOF */
            while (pos < end) {
                if (*pos == quote && (pos + 1 >= end || is_ws(*(pos + 1)))) {
                    std::string val(start, pos - start);
                    pos++;  /* skip closing quote */
                    return Token{T_VALUE, val};
                }
                pos++;
            }
            return Token{T_VALUE, std::string(start, end - start)};
        }

        /* Bare word: tag, keyword, or unquoted value */
        const char *start = pos;
        while (pos < end && !is_ws(*pos)) pos++;
        std::string word(start, pos - start);

        /* Check for keywords (case-insensitive) */
        if (word.size() >= 5) {
            std::string lower = word;
            for (auto &c : lower) c = tolower((unsigned char)c);

            if (lower.compare(0, 5, "data_") == 0)
                return Token{T_DATA, word.substr(5)};
            if (lower == "loop_")
                return Token{T_LOOP, ""};
            /* Skip save_ and stop_ blocks */
            if (lower.compare(0, 5, "save_") == 0 || lower == "stop_")
                return read_token();  /* skip and get next */
        }

        /* Tag: starts with '_' */
        if (!word.empty() && word[0] == '_')
            return Token{T_TAG, word.substr(1)};  /* strip leading '_' */

        return Token{T_VALUE, word};
    }
};

/* ------------------------------------------------------------------
 * Parser
 * ------------------------------------------------------------------ */

/* Split a tag "category.column" into (category, column).
 * If no dot, category is empty. */
static void split_tag(const std::string &tag, std::string &cat, std::string &col) {
    size_t dot = tag.find('.');
    if (dot == std::string::npos) {
        cat.clear();
        col = tag;
    } else {
        cat = tag.substr(0, dot);
        col = tag.substr(dot + 1);
    }
    /* Lowercase the category for consistent matching */
    for (auto &c : cat) c = tolower((unsigned char)c);
}

static DataBlock parse_data(Tokenizer &tok) {
    DataBlock block;

    while (tok.peek_type() != T_EOF) {
        TokenType tt = tok.peek_type();

        if (tt == T_DATA) {
            /* If we already have a block, stop (we return first block only) */
            if (!block.name.empty()) break;
            Token t = tok.next();
            block.name = t.text;
            continue;
        }

        if (tt == T_LOOP) {
            tok.next();  /* consume loop_ */

            /* Collect all tags */
            Table table;
            while (tok.peek_type() == T_TAG) {
                Token tag = tok.next();
                std::string cat, col;
                split_tag(tag.text, cat, col);
                if (table.category.empty())
                    table.category = cat;
                table.columns.push_back(col);
            }

            /* Collect values into rows */
            size_t ncols = table.columns.size();
            if (ncols == 0) continue;

            std::vector<std::string> row;
            row.reserve(ncols);

            while (tok.peek_type() == T_VALUE) {
                Token v = tok.next();
                row.push_back(v.text);
                if (row.size() == ncols) {
                    table.rows.push_back(std::move(row));
                    row.clear();
                    row.reserve(ncols);
                }
            }

            block.tables.push_back(std::move(table));
            continue;
        }

        if (tt == T_TAG) {
            /* Tag-value pair(s) outside loop: group consecutive same-category tags */
            Token tag = tok.next();
            std::string cat, col;
            split_tag(tag.text, cat, col);

            /* Find or create table for this category */
            Table *tbl = nullptr;
            for (auto &t : block.tables) {
                if (iequals(t.category, cat)) {
                    tbl = &t;
                    break;
                }
            }
            if (!tbl) {
                block.tables.push_back(Table());
                tbl = &block.tables.back();
                tbl->category = cat;
            }

            /* Add column and value */
            tbl->columns.push_back(col);
            Token val = tok.next();
            if (tbl->rows.empty())
                tbl->rows.push_back(std::vector<std::string>());
            tbl->rows[0].push_back(val.text);
            continue;
        }

        /* Skip unexpected tokens */
        tok.next();
    }

    return block;
}

/* ------------------------------------------------------------------
 * Public API
 * ------------------------------------------------------------------ */

DataBlock parse_file(const char *filename) {
    FILE *f = fopen(filename, "rb");
    if (!f)
        throw std::runtime_error(std::string("Cannot open file: ") + filename);

    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);

    std::string contents(sz, '\0');
    size_t nread = fread(&contents[0], 1, sz, f);
    fclose(f);
    contents.resize(nread);

    Tokenizer tok(contents.data(), contents.size());
    return parse_data(tok);
}

} /* namespace cif */
