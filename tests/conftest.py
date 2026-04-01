import pytest
import os

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
