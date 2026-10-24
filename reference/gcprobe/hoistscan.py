# ABOUTME: Scans a log of Julia's JuliaLICM optimization remarks for allocations hoisted out of
# ABOUTME: loops whose fields hold object references but are written by a non-reference ("bits") store.
import re, sys, collections
def scan(path):
    lines = open(path, errors='replace').read().split('\n')
    blocks, i = [], 0
    while i < len(lines):
        l = lines[i]
        if l.startswith('remark:') and 'escape analysis for' in l:
            loc = l.split(': escape analysis')[0][len('remark: '):]
            j = i + 1; body = []
            while j < len(lines) and not lines[j].startswith('remark:'):
                body.append(lines[j]); j += 1
            hoisted = j < len(lines) and 'hoisting gc allocation' in lines[j]
            blocks.append((loc, hoisted, body)); i = j; continue
        i += 1
    hoisted = [b for b in blocks if b[1]]
    danger = collections.Counter()
    for loc, _, body in hoisted:
        # walk MemOps fields: a field is dangerous when hasobjref: 1 and one access is a bits store
        field_obj, prev_store, bad = False, False, False
        for k, l in enumerate(body):
            s = l.strip()
            if s.startswith('hasobjref:'): field_obj = s.endswith('1')
            if s.startswith('store ') or ' = store' in s: prev_store = True; continue
            if s in ('bits', 'objref') and prev_store:
                if field_obj and s == 'bits': bad = True
                prev_store = False
            elif s.startswith('%') or s.startswith('call'): prev_store = False
        if bad: danger[loc] += 1
    print(f"{path}\n  escape analyses: {len(blocks)}  hoisted allocations: {len(hoisted)}  "
          f"hoisted with a bits store into a reference field: {sum(danger.values())}")
    for loc, n in danger.most_common(): print(f"    DANGER {n}x at {loc}")
    for loc, n in collections.Counter(b[0] for b in hoisted).most_common(15): print(f"    hoist {n}x at {loc}")
for p in sys.argv[1:]: scan(p)
