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
    ''' Provides a fixture to quickly setup profiling any MPI execution with Intel APS.'''
    try:
        subprocess.run(["vtune", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (FileNotFoundError, subprocess.CalledProcessError):
        raise ValueError("Intel APS is part of Intel Vtune. Please install or find Vtune to use this feature.")
        
    print('Intercepting subprocess calls with aps wrapper...')
    original_run = subprocess.run
    commands_executed = []

    def wrapped_run(*args, **kwargs):
        cmd = args[0] if args else kwargs.get('args')
        if cmd and isinstance(cmd, str) and 'mpirun' in cmd:
            print("Execution command intercepted.")
            print(f"Intercepted command: {cmd}, with arguments: {args[1:] if len(args) > 1 else ''}, kwargs: {kwargs}")
            print("Injecting Intel APS wrapper...")

            result = original_run('aps ' + cmd, *args[1:], **kwargs)
            if result.returncode != 0:
                print(f"Error executing command: {result.stderr.decode()}")
            else:
                working_dir = kwargs.get('cwd')
                aps_reports = glob.glob(os.path.join(working_dir, 'aps_report*'))
                newest_report = max(aps_reports, key=os.path.getctime)
                result = subprocess.run('aps --report .', check=True, cwd=newest_report, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            result = original_run(*args, **kwargs)
        return result

    monkeypatch.setattr(subprocess, 'run', wrapped_run)

    yield

