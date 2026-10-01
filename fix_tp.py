import re

def remove_allocatable_tp(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # Find allocatables
    allocatables = []
    lines = content.split('\n')
    for line in lines:
        if 'ALLOCATABLE' in line.upper():
            match = re.search(r'::\s*(.*)', line)
            if match:
                vars_part = match.group(1).strip()
                vdefs = re.split(r',\s*(?![^()]*\))', vars_part)
                for vd in vdefs:
                    vd = vd.strip()
                    vname = vd.split('(')[0].strip()
                    if vname:
                        allocatables.append(vname.lower())

    # Now parse all THREADPRIVATE lines and remove the allocatables
    new_lines = []
    for line in lines:
        match = re.search(r'!\$OMP\s+THREADPRIVATE\s*\((.*)\)', line, re.IGNORECASE)
        if match:
            vars_part = match.group(1)
            vdefs = [v.strip() for v in vars_part.split(',')]
            filtered = [v for v in vdefs if v.lower() not in allocatables]
            if not filtered:
                continue # remove line
            else:
                new_lines.append(line[:match.start()] + f"!$OMP THREADPRIVATE({','.join(filtered)})" + line[match.end():])
        else:
            new_lines.append(line)

    with open(filename, 'w') as f:
        f.write('\n'.join(new_lines))

remove_allocatable_tp('src/global.F')
remove_allocatable_tp('src/dg.F')
