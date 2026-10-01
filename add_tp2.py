import re
import sys
import os

def join_lines(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    joined = []
    curr = ''
    for line in lines:
        raw = line.rstrip('\r\n')
        s = raw.strip()
        if s.startswith('#'):
            if curr:
                joined.append((curr, False))
                curr = ''
            joined.append((s, True))
            continue
        if s.startswith(('c', 'C', '*')) and (len(raw) > 0 and raw[0] in ('c', 'C', '*')):
            continue
        clean = raw.split('!')[0]
        if not clean.strip():
            continue
        is_cont = False
        if len(clean) > 5 and clean[5] not in (' ', '0') and clean[:5].strip() == '':
            is_cont = True
            part = clean[6:]
        elif clean.strip().startswith(('&', '$')):
            is_cont = True
            part = clean.strip()[1:]
        elif curr.endswith('&'):
            is_cont = True
            curr = curr[:-1]
            part = clean.strip()
        
        if is_cont and curr:
            curr += ' ' + part.strip()
        else:
            if curr:
                joined.append((curr, False))
            curr = clean.strip()
    if curr:
        joined.append((curr, False))
    return joined

def get_variables(filename):
    vars_list = []
    joined = join_lines(filename)
    pattern = re.compile(r'^\s*(INTEGER|REAL\s*\(\s*SZ\s*\))', re.IGNORECASE)
    for line, is_dir in joined:
        if is_dir: continue
        match = pattern.search(line)
        if match:
            if 'PARAMETER' in line.upper(): continue
            # Find the part after :: if it exists, else after type
            if '::' in line:
                vars_part = line.split('::', 1)[1]
            else:
                # Fortran 77 style: INTEGER A, B, C
                vars_part = re.sub(r'^\s*(INTEGER|REAL\s*\(\s*SZ\s*\))\s*', '', line, flags=re.IGNORECASE)
            
            var_defs = re.split(r',\s*(?![^()]*\))', vars_part)
            for vd in var_defs:
                vd = vd.strip()
                if not vd: continue
                # remove array bounds
                vname = vd.split('(')[0].strip()
                # remove initialization
                vname = vname.split('=')[0].strip()
                vmatch = re.match(r'([A-Za-z0-9_]+)', vname)
                if vmatch:
                    vars_list.append(vmatch.group(1))
    return vars_list

def add_threadprivate(filename):
    vars_list = get_variables(filename)
    # Filter out variables already threadprivate
    with open(filename, 'r') as f:
        content = f.read()
    
    already_tp = set()
    for match in re.finditer(r'!\$OMP\s+THREADPRIVATE\s*\(([^)]+)\)', content, re.IGNORECASE):
        for v in match.group(1).split(','):
            already_tp.add(v.strip().lower())
            
    new_tp = [v for v in vars_list if v.lower() not in already_tp]
    
    if not new_tp:
        return
        
    print(f"Adding {len(new_tp)} THREADPRIVATE vars to {filename}")
    
    lines = content.split('\n')
    out_lines = []
    for line in lines:
        out_lines.append(line)
        if 'IMPLICIT NONE' in line.upper() or ('SAVE' in line.upper() and '!' not in line.split('S')[0]):
            pass # wait, I will just append at the end of the variable declarations
    
    # Actually, the safest place to append is right before CONTAINS or END MODULE
    insert_idx = -1
    for i, line in enumerate(lines):
        if re.match(r'^\s*CONTAINS', line, re.IGNORECASE):
            insert_idx = i
            break
    if insert_idx == -1:
        for i, line in reversed(list(enumerate(lines))):
            if re.match(r'^\s*END\s+MODULE', line, re.IGNORECASE):
                insert_idx = i
                break
                
    # group by 5 to avoid long lines
    tp_lines = []
    for i in range(0, len(new_tp), 5):
        chunk = new_tp[i:i+5]
        tp_lines.append(f"      !$OMP THREADPRIVATE({','.join(chunk)})")
        
    lines = lines[:insert_idx] + tp_lines + lines[insert_idx:]
    with open(filename, 'w') as f:
        f.write('\n'.join(lines))

add_threadprivate('src/global.F')
add_threadprivate('src/dg.F')
