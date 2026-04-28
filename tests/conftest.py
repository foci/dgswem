import pytest
import os
import subprocess
import glob

def dir_path(string):
    if os.path.isdir(string):
        return os.path.abspath(string)
    raise ValueError

def pytest_addoption(parser):
    parser.addoption("--binpath", help="Path to build directory",
                     type=dir_path,
                     required=True)

@pytest.fixture
def binpath(request):
    p = request.config.getoption("--binpath")
    return dir_path(p)

@pytest.fixture
def mpi_aps(monkeypatch):
    ''' Provides a fixture to quickly setup profiling any MPI execution inside a subprocess with Intel APS.'''
    try:
        subprocess.run(["vtune", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (FileNotFoundError, subprocess.CalledProcessError):
        pytest.skip("Intel VTune is not available, skipping MPI APS profiling fixture.")
        
    print('Intercepting subprocess calls with aps wrapper...')
    original_run = subprocess.run
    commands_executed = []

    def wrapped_run(*args, **kwargs):
        cmd = args[0] if args else kwargs.get('args')
        if cmd and isinstance(cmd, str) and 'mpirun' in cmd:
            print(f"Intercepted command: {cmd}, with arguments: {args[1:] if len(args) > 1 else ''}, kwargs: {kwargs}")
            print("Injecting Intel APS wrapper...")

            result = original_run('aps ' + cmd, *args[1:], **kwargs)
            if result.returncode != 0:
                print(f"Error executing command: {result.stderr.decode()}")
            else:
                working_dir = kwargs.get('cwd')
                aps_reports = [d for d in glob.glob(os.path.join(working_dir, 'aps_result*')) if os.path.isdir(d)]
                newest_report = max(aps_reports, key=os.path.getctime)
                print(f"Generating report for APS data in {newest_report}...")
                result = subprocess.run(['aps', '--report', '.'], check=True, cwd=newest_report, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif isinstance(cmd, list) and 'mpirun' in cmd[0]:
            raise ValueError("Injection into non-shell subprocesses is a work in progress.")
        else: 
            result = original_run(*args, **kwargs)
        return result

    monkeypatch.setattr(subprocess, 'run', wrapped_run)

    yield

