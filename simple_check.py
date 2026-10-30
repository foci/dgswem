for fname in ['/home/wukenton/dgswem/src/dg.F', '/home/wukenton/dgswem/src/global.F']:
    with open(fname) as f:
        lines = f.readlines()
    
    tp_count = 0
    tp_vars = set()
    for l in lines:
        if 'THREADPRIVATE' in l.upper():
            tp_count += 1
            idx1 = l.find('(')
            idx2 = l.rfind(')')
            if idx1 != -1 and idx2 != -1:
                for v in l[idx1+1:idx2].split(','):
                    tp_vars.add(v.strip().strip('/').lower())
                    
    print(f"{fname}: {tp_count} THREADPRIVATE directives, {len(tp_vars)} items")
