#!/usr/bin/env python3
"""
filter_ermes_res.py

Reads an ERMES/GiD post-process results file (*.res, frequency domain) and
writes a new one that keeps the mesh (coordinates, elements, Gauss points)
but only the results you select.

The list of results is read from the file itself, so any result ERMES wrote
can be selected: rE, iE, mod(E), Poynting_vector, Electron_density, Bext,
Plasma_permittivity_*, and any others present in your case.

Ways to select results (case-insensitive, comma or space separated):
    3            -> result number 3 in the printed list
    2-5          -> results 2, 3, 4 and 5
    rE           -> one result by its exact name
    E            -> group: real, imaginary and modulus of E (rE, iE, mod(E));
                    works for any field written as rX / iX / mod(X)
    S  or  P     -> shortcut for Poynting_vector
    all          -> every result (useful only to re-compress / copy)

Works with binary GiDPostEx files (gzip, ERMES default) and ASCII .res files.
Only the Python standard library is used.

HOW TO RUN
----------
The script never stops to ask anything: it takes the parameters from the
command line if given, otherwise from the USER SETTINGS block below.
If something is missing or wrong it prints an error (and the list of
available results) and stops.

From Spyder (F5):
    Edit INPUT_FILE, SELECTION (and optionally OUTPUT_FILE) and run.

From a console (arguments override the USER SETTINGS):
    python filter_ermes_res.py input.res --list            (only list results)
    python filter_ermes_res.py input.res E
    python filter_ermes_res.py input.res E S Bext
    python filter_ermes_res.py input.res "1-3, 6" -o out.res
    python filter_ermes_res.py                             (USER SETTINGS / asks)

From the Spyder console, after running the script once:
    filter_res(r"D:\\ON-WORK\\ERMES\\case.res", "E, Bext")
"""

#==============================================================================
#  USER SETTINGS  (used when running from Spyder or without arguments)
#==============================================================================

# Input flavia.res file
INPUT_FILE = r"./MASTU-2D-t20_p20.gid/MASTU-2D-t20_p20.flavia.res"

# e.g. "E", "S", "E, Bext", "1-3, 6", "all"
SELECTION = "S"   

# "" -> <input>_<selection>.res in the input folder
OUTPUT_FILE = "only_S.flavia.res"       

#==============================================================================

import gzip
import os
import re
import struct
import sys

SHORTCUTS = {"s": "Poynting_vector", "p": "Poynting_vector",
             "poynting": "Poynting_vector"}

# Binary GiDPostEx markers
END_VALUES_BIN = b"\xff\xff\xff\xff" + struct.pack("<i", 11) + b"End Values\x00"
RESULT_HDR_RE = re.compile(rb'Result "([^"]*)"[^\x00]*\x00')
BAR = "=" * 70

# --------------------------------------------------------------------------
# Reading / splitting
# --------------------------------------------------------------------------
def read_raw(path):
    """Return (bytes, was_gzipped)."""
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw), True
    return raw, False


def is_binary_gid(data):
    return b"GiDPostEx" in data[:64]


def split_binary(data):
    """
    Split a binary GiDPostEx stream into
        head   : everything before the first result (header + mesh + GPs)
        blocks : list of (result_name, header_text, bytes)
    Strings are stored as <int32 length><chars>\\0, and each 'Values' block
    ends with an id of -1 followed by the string 'End Values'.
    """
    first = None
    for m in RESULT_HDR_RE.finditer(data):
        s = m.start()
        if s >= 4 and struct.unpack("<i", data[s - 4:s])[0] == m.end() - s:
            first = s - 4
            break
    if first is None:
        raise ValueError("No results found in the file.")

    head, blocks, pos = data[:first], [], first
    while pos < len(data):
        (n,) = struct.unpack("<i", data[pos:pos + 4])
        hdr = data[pos + 4:pos + 4 + n]
        m = RESULT_HDR_RE.match(hdr)
        if not m:
            raise ValueError(f"Unexpected record at byte {pos}: {hdr[:60]!r}")
        end = data.find(END_VALUES_BIN, pos)
        if end < 0:
            raise ValueError(f"'End Values' not found for {m.group(1)!r}")
        end += len(END_VALUES_BIN)
        blocks.append((m.group(1).decode("latin-1"),
                       hdr[:-1].decode("latin-1"), data[pos:end]))
        pos = end
    return head, blocks


def split_ascii(data):
    lines = data.decode("latin-1").splitlines(keepends=True)
    head, blocks, i = [], [], 0
    while i < len(lines):
        m = re.match(r'\s*Result\s+"([^"]*)"', lines[i], re.IGNORECASE)
        if not m:
            if not blocks:
                head.append(lines[i])
            i += 1
            continue
        j = i
        while j < len(lines) and not re.match(r"\s*End\s+Values", lines[j], re.I):
            j += 1
        blocks.append((m.group(1), lines[i].strip(),
                       "".join(lines[i:j + 1]).encode("latin-1")))
        i = j + 1
    return "".join(head).encode("latin-1"), blocks


def load(in_path):
    data, gz = read_raw(in_path)
    binary = is_binary_gid(data)
    head, blocks = split_binary(data) if binary else split_ascii(data)
    fmt = ("binary GiDPostEx" if binary else "ASCII") + (" (gzip)" if gz else "")
    return head, blocks, gz, fmt


# --------------------------------------------------------------------------
# Result list and selection
# --------------------------------------------------------------------------
def find_groups(names):
    """Fields written as rX / iX / mod(X) -> {X: [names...]}."""
    groups = {}
    for n in names:
        m = re.fullmatch(r"(?:r|i)(.+)|mod\((.+)\)", n)
        if m:
            groups.setdefault(m.group(1) or m.group(2), []).append(n)
    return {k: v for k, v in groups.items() if len(v) >= 2}


def print_results(blocks):
    names = [b[0] for b in blocks]
    print("Available results:")
    for k, (name, hdr, _) in enumerate(blocks, 1):
        rtype = next((t for t in ("Scalar", "Vector", "Matrix", "ComplexScalar",
                                  "ComplexVector") if f" {t} " in f" {hdr} "), "")
        print(f"  {k:3d}) {name:28s} {rtype}")
    groups = find_groups(names)
    if groups:
        print("Groups:")
        for g, members in groups.items():
            print(f"   {g:8s} -> {', '.join(members)}")
    if "Poynting_vector" in names:
        print("   S or P   -> Poynting_vector")


def parse_selection(text, names):
    """Return the list of selected result names (in file order) or raise."""
    lower = {n.lower(): n for n in names}
    groups = {g.lower(): m for g, m in find_groups(names).items()}
    chosen = set()
    tokens = [t for t in re.split(r"[,\s;]+", str(text).strip()) if t]
    if not tokens:
        raise ValueError("Empty selection.")
    for t in tokens:
        tl = t.lower()
        if tl == "all":
            chosen.update(names)
        elif re.fullmatch(r"\d+", t):
            k = int(t)
            if not 1 <= k <= len(names):
                raise ValueError(f"Number {k} out of range 1-{len(names)}.")
            chosen.add(names[k - 1])
        elif re.fullmatch(r"\d+-\d+", t):
            a, b = map(int, t.split("-"))
            if not (1 <= a <= b <= len(names)):
                raise ValueError(f"Range {t} out of range 1-{len(names)}.")
            chosen.update(names[a - 1:b])
        elif tl in lower:                      # exact name (wins over groups)
            chosen.add(lower[tl])
        elif tl in groups:
            chosen.update(groups[tl])
        elif tl in SHORTCUTS and SHORTCUTS[tl] in names:
            chosen.add(SHORTCUTS[tl])
        else:
            raise ValueError(f"'{t}' is not a result, group or number "
                             f"available in this file.")
    return [n for n in names if n in chosen]


def auto_suffix(selection_text, selected):
    tokens = [t for t in re.split(r"[,\s;]+", str(selection_text).strip()) if t]
    if all(not re.fullmatch(r"[\d-]+", t) for t in tokens) and len(tokens) <= 4:
        parts = ["Poynting" if t.lower() in SHORTCUTS else t for t in tokens]
    else:
        parts = selected[:3] + (["etc"] if len(selected) > 3 else [])
    return "_".join(re.sub(r"[^\w]+", "", p) for p in parts) or "filtered"


# --------------------------------------------------------------------------
# Main work
# --------------------------------------------------------------------------
def filter_res(in_path, selection, out_path=""):
    """Filter one file (no interactive questions)."""
    in_path = os.path.abspath(in_path)
    if not os.path.isfile(in_path):
        raise FileNotFoundError(f"File not found: {in_path}")

    print(BAR)
    print(" ERMES .res filter (frequency domain)")
    print(BAR)
    print(f" Input file : {in_path}")
    print("Reading ...")
    head, blocks, gz, fmt = load(in_path)
    names = [b[0] for b in blocks]
    print(f"  format: {fmt},  {len(names)} results\n")
    print_results(blocks)

    if not str(selection).strip():
        raise ValueError("No results selected. Set SELECTION in the script "
                         "or give it after the file name, e.g. 'E, Bext'.")
    selected = parse_selection(selection, names)

    root, ext = os.path.splitext(in_path)
    out_path = os.path.abspath(
        out_path or f"{root}_{auto_suffix(selection, selected)}{ext or '.res'}")
    if out_path == in_path:
        raise ValueError("Output file would overwrite the input file.")

    print("\n" + BAR)
    print(f" Input file : {in_path}")
    print(f" Selected   : {', '.join(selected)}")
    print(f" Output file: {out_path}")
    print(BAR)

    out = head + b"".join(b[2] for b in blocks if b[0] in selected)
    if gz:
        with gzip.open(out_path, "wb", compresslevel=6) as f:
            f.write(out)
    else:
        with open(out_path, "wb") as f:
            f.write(out)

    print(f"Kept {len(selected)} of {len(names)} results")
    print(BAR)
    print(" OUTPUT FILE WRITTEN:")
    print(f"   {out_path}")
    print(f"   folder : {os.path.dirname(out_path)}")
    print(f"   size   : {os.path.getsize(out_path) / 1e6:.1f} MB "
          f"(input was {os.path.getsize(in_path) / 1e6:.1f} MB)")
    print(BAR)
    return out_path


def parse_cli(argv):
    """Return (input, selection, output, list_only) from the command line."""
    inp, sel, out, list_only = "", [], "", False
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-o", "--output") and i + 1 < len(argv):
            out = argv[i + 1]
            i += 1
        elif a in ("-l", "--list"):
            list_only = True
        elif a.startswith("-") and not re.fullmatch(r"-?\d+(-\d+)?", a):
            pass                                  # ignore unknown flags (IPython)
        elif not inp:
            inp = a
        else:
            sel.append(a)
        i += 1
    return inp, " ".join(sel), out, list_only


def main():
    inp, sel, out, list_only = parse_cli(sys.argv[1:])

    in_path = inp or INPUT_FILE

    try:
        if not in_path:
            raise ValueError("No input file. Set INPUT_FILE in the script "
                             "or give it as the first argument.")
        if not os.path.isfile(in_path):
            raise FileNotFoundError(f"File not found: {os.path.abspath(in_path)}")
        if list_only:
            print(f"File: {os.path.abspath(in_path)}")
            _, blocks, _, fmt = load(in_path)
            print(f"Format: {fmt}\n")
            print_results(blocks)
            return
        filter_res(in_path, sel or SELECTION, out or OUTPUT_FILE)
    except Exception as err:          # no sys.exit: keeps Spyder's console alive
        print(f"\nERROR: {err}")


if __name__ == "__main__":
    main()
