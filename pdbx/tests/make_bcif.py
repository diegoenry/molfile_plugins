"""make_bcif.py -- generate the BinaryCIF test files in tests/bcif/

    uv run --with mmcif --with msgpack python tests/make_bcif.py

The generated files are committed, so `make test` needs no Python. Run this
again only to change them. It writes:

  <name>.bcif     text fixtures from tests/data converted by RCSB's py-mmcif,
                  the encoder behind https://models.rcsb.org/<id>.bcif. Each
                  must load exactly like tests/data/<name>.cif.
  encodings.bcif  a file written here that uses every BinaryCIF encoding,
  encodings.cif   including those py-mmcif never emits, and its text twin.
                  py-mmcif's own decoder checks that the two hold the same
                  values, so this encoder is not the only reading of the spec.
  bad_*.bcif      damaged files; bad.txt lists the message each must give.
"""
import copy
import os
import struct
import sys

import msgpack
from mmcif.io.BinaryCifReader import BinaryCifReader
from mmcif.io.BinaryCifWriter import BinaryCifWriter
from mmcif.io.IoAdapterPy import IoAdapterPy

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "bcif")

# text fixtures that py-mmcif reads as this plugin does (the damaged ones
# are read differently by design, so they are left out)
CONVERT = [
    "altloc_bonds", "chains", "chem_comp", "charges_isotopes", "ihm_spheres",
    "ihm_starting", "label_vs_auth", "microhet_bond", "models", "models_mismatch",
    "multiblock", "no_type_symbol", "partial_type_symbol", "quoting",
    "single_atom_kv",
]  # (py-mmcif 1.2.1 cannot encode quoting.cif; it is skipped with a note)

# ---- encoders (the inverse of the decoding steps in pdbxplugin.C) ----------

TYPES = {1: "b", 2: "h", 3: "i", 4: "B", 5: "H", 6: "I", 32: "f", 33: "d"}


def byte_array(values, typ):
    return struct.pack("<%d%s" % (len(values), TYPES[typ]), *values), {"kind": "ByteArray", "type": typ}


def integer_packing(values, byte_count, unsigned):
    upper = (0xFF if byte_count == 1 else 0xFFFF) if unsigned else (0x7F if byte_count == 1 else 0x7FFF)
    lower = -upper - 1
    out = []
    for v in values:
        if v >= 0:
            while v >= upper:
                out.append(upper)
                v -= upper
        else:
            assert not unsigned
            while v <= lower:
                out.append(lower)
                v -= lower
        out.append(v)
    typ = {(1, True): 4, (1, False): 1, (2, True): 5, (2, False): 2}[(byte_count, unsigned)]
    return out, typ, {"kind": "IntegerPacking", "byteCount": byte_count,
                      "isUnsigned": unsigned, "srcSize": len(values)}


def delta(values):
    if not values:
        return [], {"kind": "Delta", "origin": 0, "srcType": 3}
    out = [0] + [values[i] - values[i - 1] for i in range(1, len(values))]
    return out, {"kind": "Delta", "origin": values[0], "srcType": 3}


def run_length(values):
    out = []
    i = 0
    while i < len(values):
        j = i
        while j < len(values) and values[j] == values[i]:
            j += 1
        out += [values[i], j - i]
        i = j
    return out, {"kind": "RunLength", "srcType": 3, "srcSize": len(values)}


def int_chain(values, steps):
    """steps: names applied in order, ending with a packing or ByteArray."""
    encs = []
    for s in steps:
        if s == "delta":
            values, e = delta(values)
        elif s == "rle":
            values, e = run_length(values)
        elif s.startswith("pack"):        # pack1u, pack1s, pack2u, pack2s
            values, typ, e = integer_packing(values, int(s[4]), s[5] == "u")
            encs.append(e)
            data, e = byte_array(values, typ)
            encs.append(e)
            return data, encs
        elif s.startswith("bytes"):       # bytes:<type>
            data, e = byte_array(values, int(s.split(":")[1]))
            encs.append(e)
            return data, encs
        encs.append(e)
    raise ValueError("chain must end with packing or bytes")


def string_column(strings, data_steps, offset_steps):
    """strings: list of str or None (masked)."""
    table = []
    for s in strings:
        if s is not None and s not in table:
            table.append(s)
    offsets = [0]
    for s in table:
        offsets.append(offsets[-1] + len(s))     # characters, as py-mmcif counts
    indices = [-1 if s is None else table.index(s) for s in strings]
    idata, ienc = int_chain(indices, data_steps)
    odata, oenc = int_chain(offsets, offset_steps)
    return idata, [{"kind": "StringArray", "dataEncoding": ienc, "stringData": "".join(table),
                    "offsetEncoding": oenc, "offsets": odata}]


def column(name, data, encoding, mask=None):
    m = None
    if mask is not None and any(mask):
        mdata, menc = int_chain(mask, ["rle", "pack1u"])
        m = {"data": mdata, "encoding": menc}
    return {"name": name, "data": {"data": data, "encoding": encoding}, "mask": m}


def mask_of(values):
    return [0 if v not in (".", "?") else (1 if v == "." else 2) for v in values]


# ---- the hand-made file ----------------------------------------------------

ATOMS = [  # group, type, name, alt, comp, asym, entity, seq, ins, occ, charge, aseq, aasym
    ("ATOM",   "N",  "N",  ".", "CYS", "A", 1, 1,   "?", 100, "?", 1,    "A"),
    ("ATOM",   "S",  "SG", ".", "CYS", "A", 1, 1,   "?", 100, "?", 1,    "A"),
    ("ATOM",   "S",  "SG", "A", "CYS", "A", 1, 2,   "?", 50,  "?", 300,  "A"),
    ("ATOM",   "S",  "SG", "B", "CYS", "A", 1, 2,   "?", 50,  "?", 300,  "A"),
    ("HETATM", "C",  "Ä1", ".", "LIG", "B", 2, ".", "?", 100, 1,   -5,   "B"),
    ("HETATM", "O",  "O",  ".", "HOH", "C", 3, ".", "B", 100, "?", 1000, "W"),
]
# x: fixed point (1000) with a jump of more than 32.767 A between atoms
X = [[12345, 12845, -40000, -40010, 1500, 0], [12346, 12846, -39999, -40009, 1501, 1]]
Y = [[1.5, -2.25, 0.1, 0.5, 3.75, -0.125], [1.5, -2.25, 0.1, 0.5, 3.75, 0.125]]   # float32
Z = [[-10.0, -9.5, 0.0, 0.5, 9.5, 10.0], [-9.0, -8.5, 1.0, 1.5, 8.5, 9.0]]          # 0.5 steps
B = [12.34, 0.0, 99.99, 100.0, 5.5, 20.0]                                           # float64


def atom_site_rows():
    rows = []
    for m in range(2):
        for i, a in enumerate(ATOMS):
            rows.append((m, i, a))
    return rows


def encodings_text():
    lines = ["data_first", "_entry.id FIRST", "#", "data_second", "#",
             "_cell.length_a 50.5", "_cell.length_b 60.25", "_cell.length_c 70.125",
             "_cell.angle_alpha 90.0", "_cell.angle_beta 100.5", "_cell.angle_gamma 90.0", "#",
             "loop_"]
    items = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
             "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x",
             "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
             "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]
    lines += ["_atom_site." + i for i in items]
    for m, i, a in atom_site_rows():
        g, ty, nm, alt, comp, asym, ent, seq, ins, occ, chg, aseq, aasym = a
        x = X[m][i]
        lines.append(" ".join(str(v) for v in [
            g, m * 6 + i + 1, ty, nm, alt, comp, asym, ent, seq, ins,
            "%s%d.%03d" % ("-" if x < 0 else "", abs(x) // 1000, abs(x) % 1000),
            repr(Y[m][i]), repr(Z[m][i]), "%d.%02d" % (occ // 100, occ % 100), repr(B[i]),
            chg, aseq, comp, aasym, nm, m + 1]))
    lines += ["#", "loop_", "_struct_conn.id", "_struct_conn.conn_type_id",
              "_struct_conn.ptnr1_auth_asym_id", "_struct_conn.ptnr1_auth_seq_id",
              "_struct_conn.ptnr1_label_atom_id", "_struct_conn.ptnr2_auth_asym_id",
              "_struct_conn.ptnr2_auth_seq_id", "_struct_conn.ptnr2_label_atom_id",
              "disulf1 disulf A 1 SG A 300 SG", "#",
              "data_third", "loop_", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
              "_atom_site.Cartn_z", "9 9 9", "#", ""]
    return "\n".join(lines)


def encodings_doc():
    rows = atom_site_rows()
    col = lambda k: [r[2][k] for r in rows]
    n = len(rows)
    cols = []
    d, e = string_column(col(0), ["rle", "pack1u"], ["delta", "pack1u"])
    cols.append(column("group_PDB", d, e))
    d, e = int_chain(list(range(1, n + 1)), ["delta", "rle", "pack1u"])
    cols.append(column("id", d, e))
    d, e = string_column(col(1), ["bytes:3"], ["bytes:6"])                 # Int32 / Uint32
    cols.append(column("type_symbol", d, e))
    d, e = string_column(col(2), ["delta", "pack2s"], ["delta", "rle", "pack1u"])
    cols.append(column("label_atom_id", d, e))
    alts = col(3)
    d, e = string_column([None if v == "." else v for v in alts], ["rle", "pack1s"], ["pack1u"])
    cols.append(column("label_alt_id", d, e, mask_of(alts)))
    d, e = string_column(col(4), ["delta", "rle", "pack2s"], ["delta", "pack2u"])
    cols.append(column("label_comp_id", d, e))
    d, e = string_column(col(5), ["bytes:1"], ["bytes:5"])                 # Int8 / Uint16
    cols.append(column("label_asym_id", d, e))
    d, e = int_chain(col(6), ["bytes:5"])                                    # Uint16
    cols.append(column("label_entity_id", d, e))
    seqs = col(7)
    d, e = int_chain([0 if v == "." else v for v in seqs], ["bytes:3"])     # Int32
    cols.append(column("label_seq_id", d, e, mask_of(seqs)))
    ins = col(8)
    d, e = string_column([None if v == "?" else v for v in ins], ["pack1s"], ["pack1u"])
    cols.append(column("pdbx_PDB_ins_code", d, e, mask_of(ins)))
    xs = [X[m][i] for m, i, _ in rows]
    d, e = int_chain(xs, ["delta", "pack2s"])
    cols.append(column("Cartn_x", d, [{"kind": "FixedPoint", "factor": 1000, "srcType": 33}] + e))
    d, e = byte_array([Y[m][i] for m, i, _ in rows], 32)                    # Float32
    cols.append(column("Cartn_y", d, [e]))
    q = [int(round((Z[m][i] + 10.0) / 0.5)) for m, i, _ in rows]
    d, e = int_chain(q, ["bytes:5"])
    cols.append(column("Cartn_z", d, [{"kind": "IntervalQuantization", "min": -10.0, "max": 10.0,
                                       "numSteps": 41, "srcType": 32}] + e))
    d, e = int_chain(col(9), ["rle", "pack1u"])
    cols.append(column("occupancy", d, [{"kind": "FixedPoint", "factor": 100, "srcType": 32}] + e))
    d, e = byte_array([B[i] for _, i, _ in rows], 33)                       # Float64
    cols.append(column("B_iso_or_equiv", d, [e]))
    chg = col(10)
    d, e = int_chain([0 if v == "?" else v for v in chg], ["bytes:1"])     # Int8
    cols.append(column("pdbx_formal_charge", d, e, mask_of(chg)))
    d, e = int_chain(col(11), ["delta", "pack1s"])                          # multi-unit values
    cols.append(column("auth_seq_id", d, e))
    d, e = string_column(col(4), ["rle", "pack1u"], ["delta", "pack1u"])
    cols.append(column("auth_comp_id", d, e))
    d, e = string_column(col(12), ["rle", "pack1u"], ["pack1u"])
    cols.append(column("auth_asym_id", d, e))
    d, e = string_column(col(2), ["bytes:3"], ["bytes:3"])
    cols.append(column("auth_atom_id", d, e))
    d, e = int_chain([m + 1 for m, _, _ in rows], ["bytes:6"])              # Uint32
    cols.append(column("pdbx_PDB_model_num", d, e))

    def f64col(name, v):
        d, e = byte_array([v], 33)
        return column(name, d, [e])

    def strcol(name, vals):
        d, e = string_column(vals, ["pack1u"], ["pack1u"])
        return column(name, d, e)

    def intcol(name, vals):
        d, e = int_chain(vals, ["pack2u"])
        return column(name, d, e)

    cell = {"name": "_cell", "rowCount": 1, "columns": [
        f64col("length_a", 50.5), f64col("length_b", 60.25), f64col("length_c", 70.125),
        f64col("angle_alpha", 90.0), f64col("angle_beta", 100.5), f64col("angle_gamma", 90.0)]}
    # no leading '_' and upper-case column names: both are normalized
    conn = {"name": "struct_conn", "rowCount": 1, "columns": [
        strcol("ID", ["disulf1"]), strcol("CONN_TYPE_ID", ["disulf"]),
        strcol("PTNR1_AUTH_ASYM_ID", ["A"]), intcol("PTNR1_AUTH_SEQ_ID", [1]),
        strcol("PTNR1_LABEL_ATOM_ID", ["SG"]), strcol("PTNR2_AUTH_ASYM_ID", ["A"]),
        intcol("PTNR2_AUTH_SEQ_ID", [300]), strcol("PTNR2_LABEL_ATOM_ID", ["SG"])]}
    third = {"name": "_atom_site", "rowCount": 1, "columns": [
        f64col("Cartn_x", 9.0), f64col("Cartn_y", 9.0), f64col("Cartn_z", 9.0)]}
    return {"version": "0.3.0", "encoder": "pdbx tests/make_bcif.py", "dataBlocks": [
        {"header": "first", "categories": [{"name": "_entry", "rowCount": 1,
                                             "columns": [strcol("id", ["FIRST"])]}]},
        {"header": "second", "categories": [cell, {"name": "_atom_site", "rowCount": n,
                                                   "columns": cols}, conn]},
        {"header": "third", "categories": [third]}]}


def check_with_pymmcif(bcif_path, cif_path):
    """py-mmcif's decoder must read the same values as its text reader."""
    b = {c.getName(): c for c in BinaryCifReader().deserialize(bcif_path)}
    t = {c.getName(): c for c in IoAdapterPy().readFile(cif_path)}
    assert sorted(b) == sorted(t), (sorted(b), sorted(t))
    for block in t:
        for cat in t[block].getObjNameList():
            # py-mmcif drops the first character of a category name, assuming
            # the leading '_': struct_conn (written without it) is "truct_conn"
            tc = t[block].getObj(cat)
            bc = b[block].getObj(cat) or b[block].getObj(cat[1:])
            assert bc is not None, (block, cat)
            for item in tc.getAttributeList():
                for r in range(tc.getRowCount()):
                    tv, bv = tc.getValue(item, r), bc.getValue(item.upper() if cat == "struct_conn" else item, r)
                    if tv in (".", "?"):
                        ok = bv in (None, ".", "?", tv)
                    else:
                        try:
                            ok = abs(float(tv) - float(bv)) <= 1e-6 * max(1.0, abs(float(tv)))
                        except (TypeError, ValueError):
                            ok = str(bv) == tv
                    assert ok, (block, cat, item, r, tv, bv)


# ---- damaged files ---------------------------------------------------------

def atom_site_of(doc):
    return [c for c in doc["dataBlocks"][1]["categories"] if c["name"] == "_atom_site"][0]


def col_of(doc, name):
    return [c for c in atom_site_of(doc)["columns"] if c["name"] == name][0]


def bad_files(doc, packed):
    bad = []

    def add(name, data, message):
        bad.append((name, data, message))

    add("bad_truncated", packed[: len(packed) // 2], "ends in the middle of a value")
    d = copy.deepcopy(doc); atom_site_of(d)["rowCount"] += 1
    add("bad_rowcount", msgpack.packb(d), "does not match the row count")
    d = copy.deepcopy(doc); col_of(d, "id")["data"]["encoding"][1]["srcSize"] = 10 ** 9
    add("bad_runlength", msgpack.packb(d), "far more values than there are rows")
    d = copy.deepcopy(doc); col_of(d, "type_symbol")["data"]["encoding"][0]["stringData"] = "N"
    add("bad_offsets", msgpack.packb(d), "string offset is out of range")
    d = copy.deepcopy(doc); c = col_of(d, "type_symbol")
    c["data"]["data"] = struct.pack("<12i", *([99] * 12))
    add("bad_index", msgpack.packb(d), "string index is out of range")
    d = copy.deepcopy(doc); col_of(d, "B_iso_or_equiv")["data"]["encoding"][0]["kind"] = "Zip"
    add("bad_kind", msgpack.packb(d), "unsupported encoding 'Zip'")
    d = copy.deepcopy(doc); col_of(d, "B_iso_or_equiv")["data"]["data"] += b"\0"
    add("bad_bytes", msgpack.packb(d), "not a whole number of values")
    d = copy.deepcopy(doc); c = col_of(d, "auth_seq_id")
    c["data"]["data"] = c["data"]["data"][:-1] + struct.pack("<b", 127)
    add("bad_packing", msgpack.packb(d), "IntegerPacking data ends inside a value")
    d = copy.deepcopy(doc); c = col_of(d, "id"); c["data"]["encoding"][0]["origin"] = 2 ** 31 - 3
    add("bad_delta", msgpack.packb(d), "Delta values do not fit in 32 bits")
    add("bad_noblocks", msgpack.packb({"version": "0.3.0", "encoder": "x"}), "no list of data blocks")
    nested = {"version": "0.3.0", "dataBlocks": []}
    for _ in range(40):
        nested = {"version": [nested]}
    add("bad_nesting", msgpack.packb(nested), "nested too deeply")
    add("bad_ext", msgpack.packb({"version": "0.3.0", "dataBlocks": [msgpack.ExtType(5, b"x")]}),
        "unsupported MessagePack type")
    d = copy.deepcopy(doc); atom_site_of(d)["rowCount"] = 2 ** 27
    c = col_of(d, "pdbx_PDB_model_num")       # decoded first: 2^27 values from 8 bytes
    c["data"]["encoding"] = [{"kind": "RunLength", "srcType": 3, "srcSize": 2 ** 27},
                             {"kind": "ByteArray", "type": 3}]
    c["data"]["data"] = struct.pack("<2i", 1, 2 ** 27)
    cols = atom_site_of(d)["columns"]
    cols.insert(0, cols.pop(cols.index(c)))
    add("bad_bomb", msgpack.packb(d), "implausibly large for the file size")
    return bad


def main():
    os.makedirs(OUT, exist_ok=True)
    converted = 0
    for name in CONVERT:
        src = os.path.join(DATA, name + ".cif")
        dst = os.path.join(OUT, name + ".bcif")
        if os.path.exists(dst):
            os.remove(dst)
        ok = BinaryCifWriter(dictionaryApi=None, useFloat64=True).serialize(
            dst, IoAdapterPy().readFile(src))
        if ok is False or not os.path.exists(dst):
            print("note: py-mmcif could not encode %s.cif; skipped" % name)
            if os.path.exists(dst):
                os.remove(dst)
            continue
        converted += 1

    doc = encodings_doc()
    packed = msgpack.packb(doc, use_bin_type=True)
    with open(os.path.join(OUT, "encodings.bcif"), "wb") as f:
        f.write(packed)
    with open(os.path.join(OUT, "encodings.cif"), "w", encoding="utf-8") as f:
        f.write(encodings_text())
    check_with_pymmcif(os.path.join(OUT, "encodings.bcif"), os.path.join(OUT, "encodings.cif"))

    lines = []
    for name, data, message in bad_files(doc, packed):
        with open(os.path.join(OUT, name + ".bcif"), "wb") as f:
            f.write(data)
        lines.append("%s.bcif|%s" % (name, message))
    with open(os.path.join(OUT, "bad.txt"), "w") as f:
        f.write("# file|message it must give\n" + "\n".join(lines) + "\n")
    print("wrote %d files to %s" % (converted + 2 + len(lines) + 1, OUT))


if __name__ == "__main__":
    sys.exit(main())
