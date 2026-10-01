import re

def get_scalars(filename):
    scalars = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.split('!')[0].split('C')[0].split('c')[0]  # strip comments roughly
        if re.match(r'^\s*(?:INTEGER|REAL\s*\(\s*SZ\s*\))\s*(?:,\s*[A-Za-z0-9_]+)*\s*(?:::\s*)?(.*)', line, re.IGNORECASE):
            if 'POINTER' in line.upper() or 'ALLOCATABLE' in line.upper() or 'PARAMETER' in line.upper(): continue
            match = re.match(r'^\s*(?:INTEGER|REAL\s*\(\s*SZ\s*\))\s*(?:,\s*[A-Za-z0-9_]+)*\s*(?:::\s*)?(.*)', line, re.IGNORECASE)
            if match:
                vars_part = match.group(1).strip()
                if not vars_part: continue
                # Split by commas not inside parens
                var_defs = re.split(r',\s*(?![^()]*\))', vars_part)
                for vd in var_defs:
                    vd = vd.strip()
                    if not vd: continue
                    if '(' in vd: continue  # it's an array
                    if '=' in vd:
                        vd = vd.split('=')[0].strip() # remove initialization
                    if vd:
                        scalars.append(vd)
    return scalars

print("GLOBAL:")
print(get_scalars('src/global.F'))
print("DG:")
print(get_scalars('src/dg.F'))
