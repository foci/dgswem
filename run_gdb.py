import subprocess

gdb_cmds = """
set pagination off
set confirm off
b DG_HYDRO_TIMESTEP_STAGE_INTEGRATE
b DG_HYDRO_TIMESTEP_STAGE_POSTPROCESS
run
print "At stage integrate 1"
print ZE(1, 1:MNE, 1)
c
print "At stage postprocess 1"
print ZE(1, 1:MNE, 1)
quit
"""

with open("gdb_script.txt", "w") as f:
    f.write(gdb_cmds)

# run gdb-oneapi on dagswem
p = subprocess.run(["gdb-oneapi", "-x", "gdb_script.txt", "../../build-debug/dagswem"],
                   cwd="/home/wukenton/dgswem/tests/sbdag", capture_output=True, text=True)
print("STDOUT:\n", p.stdout[-2000:])
print("STDERR:\n", p.stderr[-1000:])
