import sys, os
import pytest
import subprocess

# Directories required, modify if package structure changes
TEST_DIR       = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR     = f'{TEST_DIR}/../fetchmgs'         # Points to the package directory
TEST_DATA_DIR  = f'{TEST_DIR}/test_data/'

# Load fetchMGs functions
sys.path.insert(0, PARENT_DIR) 
from fetchMGs import cli


def test_cli_no_mode():
    """ Test for the command line interface output when no mode is defined """

    proc = subprocess.Popen(["python", f"{PARENT_DIR}/fetchMGs.py", "-h"], stdout=subprocess.PIPE)
    test_out = proc.communicate()[0].decode()
    print(test_out)

    original_out = open(f"{TEST_DATA_DIR}/cli_no_mode.txt", "r")
    print(original_out.read())

    print(test_out == original_out)
test_cli_no_mode()