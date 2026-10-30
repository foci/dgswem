import subprocess, hashlib, sys

for i in range(1, 21):
    r1 = subprocess.run(["/home/wukenton/dgswem/build-debug/dagswem"], cwd="/home/wukenton/dgswem/tests/slopingbeach", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p = subprocess.Popen(["/home/wukenton/dgswem/build-debug/adcpost"], cwd="/home/wukenton/dgswem/tests/slopingbeach", stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)
    p.communicate(input="post\nN\nquit\n")
    with open("/home/wukenton/dgswem/tests/slopingbeach/fort.63", "rb") as f:
        h = hashlib.md5(f.read()).hexdigest()
    print(f"Run {i}: {h}", flush=True)
