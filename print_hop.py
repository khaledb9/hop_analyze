#!/usr/bin/env python3
import re
import sys

HEADER_RE = re.compile(
    r'^Hopping\s+<a\|H\|b>\s+between\s+(?P<pair>.+?)\s+in\s+sphere\s+#\s*(?P<sphere>\d+)'
    r'\s+with\s+radius\s+(?P<radius>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    r'(?:\s+--\s*(?P<idx>\d+):)?\s*$'
)
RADIUS_VEC_RE = re.compile(r'^\s*Radius\s+vector\s+is:\s*(.*)$')
NUM_RE = re.compile(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?')

def clean_pair_label(raw_pair: str) -> str:
    s = re.sub(r'\s*\([^)]*\)', '', raw_pair)     # drop all “( … )”
    s = re.sub(r'\s*<-->\s*', '<-->', s)          # normalize arrow spacing
    return s.strip()

def parse_blocks(lines):
    """
    Yields blocks with keys: pair, radius, matrix_lines.
    Robust to: blank separators, lone '--', and back-to-back headers.
    """
    current = None
    after_radius_vec = False

    for raw in lines:
        line = raw.rstrip("\n")

        # If we see a new header, flush previous (if any)
        m = HEADER_RE.match(line)
        if m:
            if current is not None and current['matrix_lines']:
                yield current
            current = {
                'pair': clean_pair_label(m.group('pair').strip()),
                'radius': float(m.group('radius')),
                'matrix_lines': [],
            }
            after_radius_vec = False
            continue

        if current is None:
            continue

        if RADIUS_VEC_RE.match(line):
            after_radius_vec = True
            continue

        # End block on blank / naked separator
        if not line.strip() or line.strip() == '--':
            if current['matrix_lines']:
                yield current
                current = None
            continue

        # Collect numeric matrix lines after the radius-vector marker
        if after_radius_vec and NUM_RE.search(line):
            current['matrix_lines'].append(line)

    # Tail flush
    if current is not None and current['matrix_lines']:
        yield current

def max_abs_from_lines(matrix_lines):
    max_val = 0.0
    for line in matrix_lines:
        for s in NUM_RE.findall(line):
            v = abs(float(s))
            if v > max_val:
                max_val = v
    return max_val

def group_by_distance(rows, gap=2.0, eps=1e-9):
    """
    rows: list[(pair:str, hop:float, dist:float)]
    New group starts only if dist - last_dist > gap + eps
    """
    rows_sorted = sorted(rows, key=lambda r: r[2])
    groups, cur, last_d = [], [], None
    for r in rows_sorted:
        d = r[2]
        if last_d is None or (d - last_d) <= (gap + eps):
            cur.append(r)
        else:
            if cur: groups.append(cur)
            cur = [r]
        last_d = d
    if cur:
        groups.append(cur)
    return groups

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input_file> [--gap FLOAT]")
        sys.exit(1)

    path = sys.argv[1]
    gap = 2.0
    if '--gap' in sys.argv:
        try:
            gap = float(sys.argv[sys.argv.index('--gap') + 1])
        except Exception:
            print("Invalid --gap value; using default 2.0", file=sys.stderr)

    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        blocks = list(parse_blocks(f))

    rows = []
    for b in blocks:
        rows.append((b['pair'], max_abs_from_lines(b['matrix_lines']), float(b['radius'])))

    groups = group_by_distance(rows, gap=gap, eps=1e-9)

    # Pretty print (3 columns, by your requested order)
    for gi, grp in enumerate(groups, 1):
        dmin, dmax = min(r[2] for r in grp), max(r[2] for r in grp)
        print(f"=== Distance group {gi}  [{dmin:.6f} .. {dmax:.6f}] ===")
        print(f"{'PAIR':40s}\t{'HOPPING':>12s}\t{'DISTANCE':>12s}")
        # sort each group by distance, then by |hopping| desc
        for pair, hop, dist in sorted(grp, key=lambda r: (r[2], -r[1])):
            print(f"{pair:40s}\t{hop:12.6f}\t{dist:12.6f}")
        print()

if __name__ == '__main__':
    main()
