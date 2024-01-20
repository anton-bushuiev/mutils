import sys
import inspect
import pathlib


# Paths
# Root
MUTILS_ROOT_DIR = pathlib.Path(__file__).parent.absolute()
# Data
MUTILS_DATA_DIR = MUTILS_ROOT_DIR.parent / 'data'
MUTILS_SKEMPI2_DIR = MUTILS_DATA_DIR / 'SKEMPI2'
MUTILS_C380_DIR = MUTILS_DATA_DIR / 'C380'


def export():
    for name, obj in inspect.getmembers(sys.modules[__name__]):
        # Export paths
        if isinstance(obj, pathlib.PosixPath):
            val = str(obj)
            print(f'export {name}={val}')
        # Export string variables
        if isinstance(obj, str):
            print(f'export {name}={obj}')
