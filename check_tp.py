import re

for mod_file in ['/home/wukenton/dgswem/src/dg.F', '/home/wukenton/dgswem/src/global.F']:
    with open(mod_file) as f:
        text = f.read()
    
    spec = text.split('CONTAINS')[0]
    
    tp_vars = set()
    for m in re.finditer(r'!\$OMP\s+THREADPRIVATE\s*\((.*?)\)', spec, re.IGNORECASE):
        for v in m.group(1).split(','):
            tp_vars.add(v.strip().lower().strip('/'))
            
    declared = set()
    for line in spec.splitlines():
        line = line.strip()
        if line.startswith('!') or line.startswith('C') or line.startswith('c') or line.startswith('#'):
            continue
        m = re.match(r'^(?:INTEGER|REAL\s*\(.*?\)|LOGICAL|CHARACTER\s*\(.*?\))\s*(?:::)?\s*(.*)', line, re.IGNORECASE)
        if m:
            rest = m.group(1)
            if 'PARAMETER' in line.upper():
                continue
            for item in rest.split(','):
                vname = item.strip().split('(')[0].split('=')[0].strip()
                if vname and not vname.startswith('!'):
                    declared.add(vname.lower())
                    
    not_tp = [v for v in declared if v not in tp_vars]
    print(f'=== {mod_file} ===')
    print(f'Declared: {len(declared)}, Not TP: {len(not_tp)}')
    print('Sample not TP:', not_tp[:30])
