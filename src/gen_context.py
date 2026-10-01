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

def parse_pointers(filename):
    pointers = []
    joined_lines = join_lines(filename)
    pattern = re.compile(r'([A-Za-z0-9_\(\)]+)\s*,\s*POINTER\s*::(.*)', re.IGNORECASE)
    for line, is_dir in joined_lines:
        if is_dir:
            continue
        if re.match(r'^\s*CONTAINS\b', line, re.IGNORECASE):
            break
        match = pattern.search(line)
        if match:
            type_decl = match.group(1).strip()
            vars_part = match.group(2).strip()
            
            var_defs = re.split(r',\s*(?![^()]*\))', vars_part)
            for vd in var_defs:
                vd = vd.strip()
                if not vd: continue
                vmatch = re.match(r'([A-Za-z0-9_]+)(\(.*?\))?', vd.replace(' ', ''))
                if vmatch and vmatch.group(1).lower() not in ['dx', 'dy', 'y1', 'y2', 'initial_value', 'alternate_value', 'function']:
                    name = vmatch.group(1)
                    shape = vmatch.group(2) if vmatch.group(2) else ""
                    rank = 0
                    if shape:
                        rank = shape.count(':')
                    pointers.append((name, type_decl, rank))
    return pointers

def parse_allocatables(filename):
    allocs = []
    joined_lines = join_lines(filename)
    pattern = re.compile(r'\bALLOCATABLE\b.*?::\s*(.*)', re.IGNORECASE)

    for line, is_dir in joined_lines:
        if is_dir:
            if re.match(r'#\s*if', line, re.IGNORECASE):
                allocs.append((line, '__DIRECTIVE__'))
            elif re.match(r'#\s*endif', line, re.IGNORECASE):
                allocs.append((line, '__DIRECTIVE__'))
            elif re.match(r'#\s*else', line, re.IGNORECASE):
                allocs.append((line, '__DIRECTIVE__'))
            elif re.match(r'#\s*elif', line, re.IGNORECASE):
                allocs.append((line, '__DIRECTIVE__'))
            continue

        if re.match(r'^\s*CONTAINS\b', line, re.IGNORECASE):
            break

        match = pattern.search(line)
        if match:
            vars_part = match.group(1).strip()
            var_defs = re.split(r',\s*(?![^()]*\))', vars_part)
            for vd in var_defs:
                vd = vd.strip()
                if not vd:
                    continue
                vmatch = re.match(r'([A-Za-z0-9_]+)', vd)
                if vmatch and vmatch.group(1).lower() not in ['dx', 'dy', 'y1', 'y2', 'initial_value', 'alternate_value', 'function']:
                    name = vmatch.group(1)
                    allocs.append((name, 'VAR'))
    return allocs

def parse_scalars(filename):
    int_scalars = []
    real_scalars = []
    joined_lines = join_lines(filename)
    pattern = re.compile(r'^\s*(INTEGER|REAL\s*\(\s*SZ\s*\))', re.IGNORECASE)
    for line, is_dir in joined_lines:
        if is_dir: continue
        if re.match(r'^\s*CONTAINS\b', line, re.IGNORECASE):
            break
        match = pattern.search(line)
        if match:
            if 'PARAMETER' in line.upper() or 'POINTER' in line.upper() or 'ALLOCATABLE' in line.upper(): continue
            is_int = 'INTEGER' in match.group(1).upper()
            if '::' in line:
                vars_part = line.split('::', 1)[1]
            else:
                vars_part = re.sub(r'^\s*(INTEGER|REAL\s*\(\s*SZ\s*\))\s*', '', line, flags=re.IGNORECASE)
            
            var_defs = re.split(r',\s*(?![^()]*\))', vars_part)
            for vd in var_defs:
                vd = vd.strip()
                if not vd: continue
                vname = vd.split('=')[0].strip()
                vmatch = re.match(r'([A-Za-z0-9_]+)(\(.*?\))?', vname.replace(' ', ''))
                if vmatch and vmatch.group(1).lower() not in ['dx', 'dy', 'initial_value', 'alternate_value', 'function', 'float_type', 'np_g', 'int_initial_value', 'num_fd_records', 'num_items_per_record', 'num_records_this', 'specifier', 'dummy_buf']:
                    var_name = vmatch.group(1)
                    var_shape = vmatch.group(2) if vmatch.group(2) else ""
                    if is_int:
                        int_scalars.append((var_name, var_shape))
                    else:
                        real_scalars.append((var_name, var_shape))
    return int_scalars, real_scalars

def main():
    if len(sys.argv) >= 3:
        input_files = sys.argv[1:-1]
        output_file = sys.argv[-1]
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_files = [
            os.path.join(script_dir, 'sizes.F'),
            os.path.join(script_dir, 'dg.F'),
            os.path.join(script_dir, 'global.F'),
            os.path.join(script_dir, 'nodalattr.F')
        ]
        output_file = os.path.join(script_dir, 'dagswem_state.F90')

    all_ptrs = []
    all_allocs = []
    int_scalars = []
    real_scalars = []
    
    state_ptr_files = [f for f in input_files if os.path.basename(f) in ['sizes.F', 'dg.F', 'global.F']]
    for infile in state_ptr_files:
        all_ptrs.extend(parse_pointers(infile))
        ints, reals = parse_scalars(infile)
        int_scalars.extend(ints)
        real_scalars.extend(reals)
        
    for infile in input_files:
        all_allocs.extend(parse_allocatables(infile))
    
    seen_ptrs = set()
    unique_ptrs = []
    for p in all_ptrs:
        if p[0].lower() not in seen_ptrs:
            seen_ptrs.add(p[0].lower())
            unique_ptrs.append(p)
            
    seen_ints = set()
    unique_ints = []
    for s in int_scalars:
        if s[0].lower() not in seen_ints:
            seen_ints.add(s[0].lower())
            unique_ints.append(s)
            
    seen_reals = set()
    unique_reals = []
    for s in real_scalars:
        if s[0].lower() not in seen_reals:
            seen_reals.add(s[0].lower())
            unique_reals.append(s)

    with open(output_file, 'w') as f:
        f.write("MODULE DAGSWEM_STATE\n")
        f.write("  USE SIZES\n")
        f.write("  USE GLOBAL\n")
        f.write("  USE DG\n")
        f.write("  USE NodalAttributes, ONLY: SwanWaveRefrac, STARTDRY, FRIC, TAU0VAR, TAU0BASE, &\n")
        f.write("                             z0land, vcanopy, BridgePilings, Chezy, ManningsN, &\n")
        f.write("                             GeoidOffset, EVM, EVC\n")
        f.write("  USE FSTARPU_MOD\n")
        f.write("  USE ISO_C_BINDING\n")
        f.write("  IMPLICIT NONE\n")
        f.write("\n")
        
        # Calculate total element count for int scalars/fixed arrays
        total_int_entries = 0
        for name, shape in unique_ints:
            if shape:
                dims = [int(d) for d in re.findall(r'\d+', shape)]
                prod = 1
                for d in dims: prod *= d
                total_int_entries += prod
            else:
                total_int_entries += 1

        total_real_entries = 0
        for name, shape in unique_reals:
            if shape:
                dims = [int(d) for d in re.findall(r'\d+', shape)]
                prod = 1
                for d in dims: prod *= d
                total_real_entries += prod
            else:
                total_real_entries += 1

        # Total handles = unique_ptrs + 2 (one for int scalars, one for real scalars)
        num_state_handles = len(unique_ptrs) + 2
        num_int_scalars = len(unique_ptrs) + total_int_entries + 1
        
        f.write(f"  INTEGER, PARAMETER :: NUM_STATE_HANDLES = {num_state_handles}\n")
        f.write(f"  INTEGER, PARAMETER :: NUM_INT_SCALARS = {num_int_scalars}\n")
        f.write(f"  INTEGER, PARAMETER :: NUM_REAL_SCALARS = {total_real_entries}\n")
        f.write("  REAL(SZ), TARGET, SAVE :: DUMMY_BUF(1) = 0.0_SZ\n")
        f.write("\n")
        f.write("  CONTAINS\n")
        f.write("\n")
        
        f.write("  SUBROUTINE DGSWEM_STATE_REGISTER(handles)\n")
        f.write("    TYPE(C_PTR), INTENT(OUT) :: handles(NUM_STATE_HANDLES)\n")
        f.write("    INTEGER, POINTER :: SCALAR_INT_BUF(:)\n")
        f.write("    REAL(SZ), POINTER :: SCALAR_REAL_BUF(:)\n")
        f.write("    ALLOCATE(SCALAR_INT_BUF(NUM_INT_SCALARS))\n")
        f.write("    ALLOCATE(SCALAR_REAL_BUF(NUM_REAL_SCALARS))\n")
        
        for idx, (name, _, rank) in enumerate(unique_ptrs):
            lb = ",".join(["LBOUND(" + name + "," + str(i+1) + ")" for i in range(rank)])
            if rank == 1:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 1\n")
                f.write(f"      CALL fstarpu_vector_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 0\n")
                f.write(f"      CALL fstarpu_vector_data_register(handles({idx+1}), 0, C_LOC(DUMMY_BUF(1)), 1, C_SIZEOF(DUMMY_BUF(1)))\n")
                f.write(f"    END IF\n")
            elif rank == 2:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 1\n")
                f.write(f"      CALL fstarpu_matrix_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), SIZE({name},1), SIZE({name},2), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 0\n")
                f.write(f"      CALL fstarpu_matrix_data_register(handles({idx+1}), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))\n")
                f.write(f"    END IF\n")
            elif rank == 3:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 1\n")
                f.write(f"      CALL fstarpu_block_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), SIZE({name},1)*SIZE({name},2), SIZE({name},1), SIZE({name},2), SIZE({name},3), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 0\n")
                f.write(f"      CALL fstarpu_block_data_register(handles({idx+1}), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))\n")
                f.write(f"    END IF\n")
            elif rank == 4:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 1\n")
                f.write(f"      CALL fstarpu_tensor_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), SIZE({name},1)*SIZE({name},2), SIZE({name},1)*SIZE({name},2)*SIZE({name},3), SIZE({name},1), SIZE({name},2), SIZE({name},3), SIZE({name},4), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      SCALAR_INT_BUF({idx+1}) = 0\n")
                f.write(f"      CALL fstarpu_tensor_data_register(handles({idx+1}), 0, C_LOC(DUMMY_BUF(1)), 1, 1, 1, 1, 1, 1, 1, C_SIZEOF(DUMMY_BUF(1)))\n")
                f.write(f"    END IF\n")
            f.write(f"    NULLIFY({name})\n")
        
        # Now register the scalars and fixed arrays
        num_ptrs = len(unique_ptrs)
        offset = num_ptrs
        for name, shape in unique_ints:
            if shape:
                dims = [int(d) for d in re.findall(r'\d+', shape)]
                prod = 1
                for d in dims: prod *= d
                if len(dims) == 1:
                    for idx_dim in range(1, dims[0] + 1):
                        offset += 1
                        f.write(f"    SCALAR_INT_BUF({offset}) = {name}({idx_dim})\n")
                elif len(dims) == 2:
                    for d2 in range(1, dims[1] + 1):
                        for d1 in range(1, dims[0] + 1):
                            offset += 1
                            f.write(f"    SCALAR_INT_BUF({offset}) = {name}({d1},{d2})\n")
            else:
                offset += 1
                f.write(f"    SCALAR_INT_BUF({offset}) = {name}\n")
        f.write("    IF (vertexslope) THEN\n")
        f.write(f"      SCALAR_INT_BUF({num_int_scalars}) = 1\n")
        f.write("    ELSE\n")
        f.write(f"      SCALAR_INT_BUF({num_int_scalars}) = 0\n")
        f.write("    END IF\n")
        f.write(f"    CALL fstarpu_vector_data_register(handles({num_ptrs + 1}), 0, C_LOC(SCALAR_INT_BUF(1)), NUM_INT_SCALARS, C_SIZEOF(SCALAR_INT_BUF(1)))\n")
        f.write("    NULLIFY(SCALAR_INT_BUF)\n")
        
        offset = 0
        for name, shape in unique_reals:
            if shape:
                dims = [int(d) for d in re.findall(r'\d+', shape)]
                prod = 1
                for d in dims: prod *= d
                if len(dims) == 1:
                    for idx_dim in range(1, dims[0] + 1):
                        offset += 1
                        f.write(f"    SCALAR_REAL_BUF({offset}) = {name}({idx_dim})\n")
                elif len(dims) == 2:
                    for d2 in range(1, dims[1] + 1):
                        for d1 in range(1, dims[0] + 1):
                            offset += 1
                            f.write(f"    SCALAR_REAL_BUF({offset}) = {name}({d1},{d2})\n")
            else:
                offset += 1
                f.write(f"    SCALAR_REAL_BUF({offset}) = {name}\n")
        f.write(f"    CALL fstarpu_vector_data_register(handles({num_ptrs + 2}), 0, C_LOC(SCALAR_REAL_BUF(1)), NUM_REAL_SCALARS, C_SIZEOF(SCALAR_REAL_BUF(1)))\n")
        f.write("    NULLIFY(SCALAR_REAL_BUF)\n")
        
        f.write("  END SUBROUTINE DGSWEM_STATE_REGISTER\n")
        f.write("\n")
        
        f.write("  SUBROUTINE DGSWEM_STATE_ACTIVATE(buffers)\n")
        f.write("    TYPE(C_PTR), VALUE, INTENT(IN) :: buffers\n")
        f.write("    TYPE(C_PTR) :: curr_ptr\n")
        f.write("    INTEGER, POINTER :: SCALAR_INT_BUF(:)\n")
        f.write("    REAL(SZ), POINTER :: SCALAR_REAL_BUF(:)\n")
        num_ptrs = len(unique_ptrs)
        
        # First unpack the scalars
        f.write(f"    curr_ptr = fstarpu_vector_get_ptr(buffers, {num_ptrs})\n")
        f.write(f"    CALL c_f_pointer(curr_ptr, SCALAR_INT_BUF, shape=[NUM_INT_SCALARS])\n")
        offset = num_ptrs
        for name, shape in unique_ints:
            if shape:
                dims = [int(d) for d in re.findall(r'\d+', shape)]
                prod = 1
                for d in dims: prod *= d
                if len(dims) == 1:
                    for idx_dim in range(1, dims[0] + 1):
                        offset += 1
                        f.write(f"    {name}({idx_dim}) = SCALAR_INT_BUF({offset})\n")
                elif len(dims) == 2:
                    for d2 in range(1, dims[1] + 1):
                        for d1 in range(1, dims[0] + 1):
                            offset += 1
                            f.write(f"    {name}({d1},{d2}) = SCALAR_INT_BUF({offset})\n")
            else:
                offset += 1
                f.write(f"    {name} = SCALAR_INT_BUF({offset})\n")
        f.write(f"    vertexslope = (SCALAR_INT_BUF({num_int_scalars}) == 1)\n")
            
        f.write(f"    curr_ptr = fstarpu_vector_get_ptr(buffers, {num_ptrs + 1})\n")
        f.write(f"    CALL c_f_pointer(curr_ptr, SCALAR_REAL_BUF, shape=[NUM_REAL_SCALARS])\n")
        offset = 0
        for name, shape in unique_reals:
            if shape:
                dims = [int(d) for d in re.findall(r'\d+', shape)]
                prod = 1
                for d in dims: prod *= d
                if len(dims) == 1:
                    for idx_dim in range(1, dims[0] + 1):
                        offset += 1
                        f.write(f"    {name}({idx_dim}) = SCALAR_REAL_BUF({offset})\n")
                elif len(dims) == 2:
                    for d2 in range(1, dims[1] + 1):
                        for d1 in range(1, dims[0] + 1):
                            offset += 1
                            f.write(f"    {name}({d1},{d2}) = SCALAR_REAL_BUF({offset})\n")
            else:
                offset += 1
                f.write(f"    {name} = SCALAR_REAL_BUF({offset})\n")
        
        # Now unpack the pointers based on SCALAR_INT_BUF mask
        for idx, (name, _, rank) in enumerate(unique_ptrs):
            if rank == 1:
                f.write(f"    IF (SCALAR_INT_BUF({idx+1}) == 1) THEN\n")
                f.write(f"      curr_ptr = fstarpu_vector_get_ptr(buffers, {idx})\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_vector_get_nx(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
            elif rank == 2:
                f.write(f"    IF (SCALAR_INT_BUF({idx+1}) == 1) THEN\n")
                f.write(f"      curr_ptr = fstarpu_matrix_get_ptr(buffers, {idx})\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_matrix_get_nx(buffers, {idx}), fstarpu_matrix_get_ny(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
            elif rank == 3:
                f.write(f"    IF (SCALAR_INT_BUF({idx+1}) == 1) THEN\n")
                f.write(f"      curr_ptr = fstarpu_block_get_ptr(buffers, {idx})\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_block_get_nx(buffers, {idx}), fstarpu_block_get_ny(buffers, {idx}), fstarpu_block_get_nz(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
            elif rank == 4:
                f.write(f"    IF (SCALAR_INT_BUF({idx+1}) == 1) THEN\n")
                f.write(f"      curr_ptr = fstarpu_tensor_get_ptr(buffers, {idx})\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_tensor_get_nx(buffers, {idx}), fstarpu_tensor_get_ny(buffers, {idx}), fstarpu_tensor_get_nz(buffers, {idx}), fstarpu_tensor_get_nt(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
        
        f.write("  END SUBROUTINE DGSWEM_STATE_ACTIVATE\n")
        f.write("\n")
        
        f.write("  SUBROUTINE DGSWEM_DEALLOC_ALLOCATABLES()\n")
        seen_allocs = set()
        for item, kind in all_allocs:
            if kind == '__DIRECTIVE__':
                f.write(f"{item}\n")
            elif kind == 'VAR':
                if item.lower() not in seen_allocs:
                    seen_allocs.add(item.lower())
                    f.write(f"    IF (ALLOCATED({item})) DEALLOCATE({item})\n")
        f.write("  END SUBROUTINE DGSWEM_DEALLOC_ALLOCATABLES\n")
        f.write("\n")
        
        f.write("END MODULE DAGSWEM_STATE\n")

if __name__ == '__main__':
    main()
