"""Phase 2 shim generator: replaces each of the six original RPI toolchain files
with a thin shim that re-exports, from rpi_core.py, EXACTLY the top-level names
originally defined in that specific file (public and private, tests included) --
never a wildcard `from rpi_core import *`, since rpi_core's namespace holds all
six modules' names and an unscoped wildcard would make e.g. `python -m unittest
reassortment_detector` discover and run every other module's tests too.

The per-file name list is computed by reading each of the six files AS THEY
CURRENTLY STAND in the working tree -- Phase 1 added rpi_core.py alongside them
but did not modify them, so at the point this script runs they are still each
module's own real, pre-shim implementation.

Run once, after rpi_core.py (Phase 1) is committed and its own test suite is
verified green. Overwrites the six original files -- recoverable via git only,
this being Phase 2's own commit boundary. NOT idempotent past the first run
(a second run would read already-shimmed files and produce empty shims) --
this is intentional, matching Phase 2 being a one-time transformation.
"""
import ast
import textwrap

# name in this file -> name it now has in rpi_core.py (only entries that differ)
RENAMES = {
    "reassortment_detector.py": {"main": "_detector_main"},
    "rpi_calibration.py": {"main": "_calibration_main", "build_panel": "_calibration_build_panel"},
    "rpi_discrimination.py": {"main": "_discrimination_main"},
    "rpi_discrimination_stats.py": {},
    "validate_reassortment_detector.py": {"main": "_validate_main", "build_panel": "_validate_build_panel"},
    "pipeline.py": {"main": "_pipeline_main"},
}

GUARD_BODY = {
    "reassortment_detector.py": "sys.exit(main())",
    "rpi_calibration.py": "sys.exit(main())",
    "rpi_discrimination.py": "sys.exit(main())",
    "rpi_discrimination_stats.py": "unittest.main()",
    "validate_reassortment_detector.py": "sys.exit(main())",
    "pipeline.py": "sys.exit(main())",
}

GUARD_IMPORT = {
    "reassortment_detector.py": "sys",
    "rpi_calibration.py": "sys",
    "rpi_discrimination.py": "sys",
    "rpi_discrimination_stats.py": "unittest",
    "validate_reassortment_detector.py": "sys",
    "pipeline.py": "sys",
}

FILES = list(RENAMES)


def top_level_names_in_order(fname):
    with open(fname, encoding="utf-8") as f:
        src = f.read()
    tree = ast.parse(src, filename=fname)
    names = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.append(node.name)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    names.append(t.id)
        elif isinstance(node, ast.AnnAssign):
            if isinstance(node.target, ast.Name):
                names.append(node.target.id)
    return names


def format_import_block(plain_names, aliased_pairs):
    """`from rpi_core import (a, b, c, ...)` for plain names, one alias-import
    line per renamed name (kept separate and explicit rather than mixed into
    the parenthesized list, so a renamed symbol is visually obvious)."""
    lines = []
    if plain_names:
        wrapped = textwrap.wrap(
            ", ".join(plain_names), width=76,
            break_long_words=False, break_on_hyphens=False)
        lines.append("from rpi_core import (")
        for w in wrapped:
            lines.append(f"    {w}")
        lines.append(")")
    for orig, new in aliased_pairs:
        lines.append(f"from rpi_core import {new} as {orig}")
    return "\n".join(lines)


for fname in FILES:
    names = top_level_names_in_order(fname)
    renames = RENAMES[fname]
    plain = [n for n in names if n not in renames]
    aliased = [(n, renames[n]) for n in names if n in renames]

    import_block = format_import_block(plain, aliased)
    guard_import = GUARD_IMPORT[fname]
    guard_body = GUARD_BODY[fname]

    module_title = fname
    content = (
        f'"""{module_title} -- thin shim.\n\n'
        f"Re-exports this module's original public surface from rpi_core.py, the\n"
        f"single consolidated module the RPI toolchain now lives in (branch\n"
        f"05-module-consolidation). See notes_consolidation.md for exactly what\n"
        f"moved and why; see rpi_core.py's own banner comments for this module's\n"
        f"section within it.\n"
        f'"""\n'
        f"import {guard_import}\n\n"
        f"{import_block}\n\n\n"
        f'if __name__ == "__main__":\n'
        f"    {guard_body}\n"
    )

    with open(fname, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"wrote shim {fname}: {len(plain)} plain + {len(aliased)} aliased re-exports")
