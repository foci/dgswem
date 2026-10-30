
import subprocess, hashlib, os

ref = '92ab900e5a5e6fd9b8113ff9ba4f6a07'
for run in range(1, 20):
    cmd = 'bash -c "source /etc/profile.d/modules.sh 2>/dev/null || true; spacktivate dagswem 2>/dev/null || true; cd tests/slopingbeach && ../../build-debug/src/dagswem"'
    res = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open('tests/slopingbeach/fort.63', 'rb') as fp:
        h = hashlib.md5(fp.read()).hexdigest()
    if h != ref:
        print(f'Run {run} deviated: {h}')
        subprocess.run('cp tests/slopingbeach/fort.63 tests/slopingbeach/fort.63.bad', shell=True)
        break
    else:
        print(f'Run {run} matched')
