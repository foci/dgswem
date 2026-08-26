import re
import sys
import os

def parse_pointers(filename):
    pointers = []
    with open(filename, 'r') as f:
        content = f.read()
    
    # Simple regex to find pointers
    pattern = re.compile(r'([A-Za-z0-9_\(\)]+)\s*,\s*POINTER\s*::(.*)', re.IGNORECASE)
    for line in content.split('\n'):
        line = line.split('!')[0]
        if line.startswith('c') or line.startswith('C'):
            continue
        
        match = pattern.search(line)
        if match:
            type_decl = match.group(1).strip()
            vars_part = match.group(2).strip()
            
            var_defs = re.split(r',\s*(?![^()]*\))', vars_part)
            for vd in var_defs:
                vd = vd.strip()
                if not vd: continue
                vmatch = re.match(r'([A-Za-z0-9_]+)(\(.*?\))?', vd.replace(' ', ''))
                if vmatch:
                    name = vmatch.group(1)
                    shape = vmatch.group(2) if vmatch.group(2) else ""
                    
                    # determine rank
                    rank = 0
                    if shape:
                        rank = shape.count(':')
                    
                    pointers.append((name, type_decl, rank))
    return pointers

def main():
    if len(sys.argv) >= 3:
        input_files = sys.argv[1:-1]
        output_file = sys.argv[-1]
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_files = [os.path.join(script_dir, 'dg.F'), os.path.join(script_dir, 'global.F')]
        output_file = os.path.join(script_dir, 'dagswem_state.f90')

    all_ptrs = []
    for infile in input_files:
        all_ptrs.extend(parse_pointers(infile))
    
    # Remove duplicates
    seen = set()
    unique_ptrs = []
    for p in all_ptrs:
        if p[0].lower() not in seen:
            seen.add(p[0].lower())
            unique_ptrs.append(p)
            
    with open(output_file, 'w') as f:
        f.write("MODULE DAGSWEM_STATE\n")
        f.write("  USE SIZES\n")
        f.write("  USE GLOBAL\n")
        f.write("  USE DG\n")
        f.write("  USE FSTARPU_MOD\n")
        f.write("  USE ISO_C_BINDING\n")
        f.write("  IMPLICIT NONE\n")
        f.write("\n")
        f.write(f"  INTEGER, PARAMETER :: NUM_STATE_HANDLES = {len(unique_ptrs)}\n")
        f.write("\n")
        f.write("  CONTAINS\n")
        f.write("\n")
        
        # Generator for DGSWEM_STATE_REGISTER
        f.write("  SUBROUTINE DGSWEM_STATE_REGISTER(handles)\n")
        f.write("    TYPE(C_PTR), INTENT(OUT) :: handles(NUM_STATE_HANDLES)\n")
        
        for idx, (name, _, rank) in enumerate(unique_ptrs):
            # calculate bounds string
            lb = ",".join(["LBOUND(" + name + "," + str(i+1) + ")" for i in range(rank)])
            if rank == 0:
                pass # impossible in this codebase
            elif rank == 1:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      CALL fstarpu_vector_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      handles({idx+1}) = C_NULL_PTR\n")
                f.write(f"    END IF\n")
            elif rank == 2:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      CALL fstarpu_matrix_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), SIZE({name},1), SIZE({name},2), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      handles({idx+1}) = C_NULL_PTR\n")
                f.write(f"    END IF\n")
            elif rank == 3:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      CALL fstarpu_block_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), SIZE({name},1)*SIZE({name},2), SIZE({name},1), SIZE({name},2), SIZE({name},3), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      handles({idx+1}) = C_NULL_PTR\n")
                f.write(f"    END IF\n")
            elif rank == 4:
                f.write(f"    IF (ASSOCIATED({name})) THEN\n")
                f.write(f"      CALL fstarpu_tensor_data_register(handles({idx+1}), 0, C_LOC({name}({lb})), SIZE({name},1), SIZE({name},1)*SIZE({name},2), SIZE({name},1)*SIZE({name},2)*SIZE({name},3), SIZE({name},1), SIZE({name},2), SIZE({name},3), SIZE({name},4), C_SIZEOF({name}({lb})))\n")
                f.write(f"    ELSE\n")
                f.write(f"      handles({idx+1}) = C_NULL_PTR\n")
                f.write(f"    END IF\n")
            f.write(f"    NULLIFY({name})\n")
        f.write("  END SUBROUTINE DGSWEM_STATE_REGISTER\n")
        f.write("\n")
        
        # Generator for DGSWEM_STATE_ACTIVATE
        f.write("  SUBROUTINE DGSWEM_STATE_ACTIVATE(buffers)\n")
        f.write("    TYPE(C_PTR), VALUE, INTENT(IN) :: buffers\n")
        f.write("    TYPE(C_PTR) :: curr_ptr\n")
        for idx, (name, _, rank) in enumerate(unique_ptrs):
            if rank == 1:
                f.write(f"    curr_ptr = fstarpu_vector_get_ptr(buffers, {idx})\n")
                f.write(f"    IF (C_ASSOCIATED(curr_ptr)) THEN\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_vector_get_nx(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
            elif rank == 2:
                f.write(f"    curr_ptr = fstarpu_matrix_get_ptr(buffers, {idx})\n")
                f.write(f"    IF (C_ASSOCIATED(curr_ptr)) THEN\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_matrix_get_nx(buffers, {idx}), fstarpu_matrix_get_ny(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
            elif rank == 3:
                f.write(f"    curr_ptr = fstarpu_block_get_ptr(buffers, {idx})\n")
                f.write(f"    IF (C_ASSOCIATED(curr_ptr)) THEN\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_block_get_nx(buffers, {idx}), fstarpu_block_get_ny(buffers, {idx}), fstarpu_block_get_nz(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
            elif rank == 4:
                f.write(f"    curr_ptr = fstarpu_tensor_get_ptr(buffers, {idx})\n")
                f.write(f"    IF (C_ASSOCIATED(curr_ptr)) THEN\n")
                f.write(f"      CALL c_f_pointer(curr_ptr, {name}, shape=[fstarpu_tensor_get_nx(buffers, {idx}), fstarpu_tensor_get_ny(buffers, {idx}), fstarpu_tensor_get_nz(buffers, {idx}), fstarpu_tensor_get_nt(buffers, {idx})])\n")
                f.write(f"    ELSE\n")
                f.write(f"      NULLIFY({name})\n")
                f.write(f"    END IF\n")
        f.write("  END SUBROUTINE DGSWEM_STATE_ACTIVATE\n")
        
        f.write("END MODULE DAGSWEM_STATE\n")

if __name__ == '__main__':
    main()
