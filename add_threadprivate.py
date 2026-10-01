import re
import os

def process_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    out_lines = []
    pointers = []
    
    in_module = False
    module_name = ""
    contains_idx = -1
    
    for i, line in enumerate(lines):
        out_lines.append(line)
        if re.match(r'^\s*MODULE\s+([A-Za-z0-9_]+)', line, re.IGNORECASE):
            in_module = True
            module_name = re.match(r'^\s*MODULE\s+([A-Za-z0-9_]+)', line, re.IGNORECASE).group(1)
        if re.match(r'^\s*CONTAINS', line, re.IGNORECASE):
            contains_idx = i
            
        match = re.search(r'([A-Za-z0-9_\(\)]+)\s*,\s*POINTER\s*::(.*)', line, re.IGNORECASE)
        if match and not re.match(r'^\s*[c!*]', line, re.IGNORECASE):
            vars_part = match.group(2).strip()
            var_defs = re.split(r',\s*(?![^()]*\))', vars_part)
            for vd in var_defs:
                vd = vd.strip()
                if not vd: continue
                vmatch = re.match(r'([A-Za-z0-9_]+)(\(.*?\))?', vd.replace(' ', ''))
                if vmatch:
                    name = vmatch.group(1)
                    pointers.append(name)
                    
    if in_module and len(pointers) > 0:
        # insert before CONTAINS or at the end if no CONTAINS
        tp_lines = []
        for p in pointers:
            tp_lines.append(f"      !$OMP THREADPRIVATE({p})\n")
            
        if contains_idx != -1:
            out_lines = out_lines[:contains_idx] + tp_lines + out_lines[contains_idx:]
        else:
            # find END MODULE
            end_idx = -1
            for i in range(len(out_lines)-1, -1, -1):
                if re.match(r'^\s*END\s+MODULE', out_lines[i], re.IGNORECASE):
                    end_idx = i
                    break
            if end_idx != -1:
                out_lines = out_lines[:end_idx] + tp_lines + out_lines[end_idx:]
            else:
                out_lines.extend(tp_lines)
                
    with open(filename, 'w') as f:
        f.writelines(out_lines)

process_file('src/dg.F')
process_file('src/global.F')
