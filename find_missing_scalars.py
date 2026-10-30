import sys
sys.path.insert(0, '/home/wukenton/dgswem/src')
from gen_context import parse_scalars
import re

def check_file(fname):
    with open(fname) as f:
        content = f.read()
    tp = set()
    for line in content.splitlines():
        m = re.search(r'THREADPRIVATE\s*\((.*?)\)', line, re.IGNORECASE)
        if m:
            for v in m.group(1).split(','):
                item = v.strip().strip('/').lower()
                if item: tp.add(item)
    ints, reals = parse_scalars(fname)
    missing_ints = [s[0] for s in ints if s[0].lower() not in tp]
    missing_reals = [s[0] for s in reals if s[0].lower() not in tp]
    return missing_ints, missing_reals

dg_i, dg_r = check_file('/home/wukenton/dgswem/src/dg.F')
gl_i, gl_r = check_file('/home/wukenton/dgswem/src/global.F')

with open('/home/wukenton/dgswem/missing_scalars.txt', 'w') as out:
    out.write(f'dg.F missing reals: {dg_r}\n')
    out.write(f'global.F missing ints: {gl_i}\n')
    out.write(f'global.F missing reals: {gl_r}\n')
print('Done!')
