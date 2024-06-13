import sys, os
import pytest
import subprocess
from fetchmgs.test_config import *

def test_cli_nomode():
    """
    Test for the command line interface output when no mode is defined 
    
    If cli arguments are change, regenerate the cli txt files in test data by 
    python3 fetchmgs.py -h > ../test/test_input/cli_nomode.38.txt (Python < 3.10)
    python3 fetchmgs.py -h > ../test/test_input/cli_nomode.310.txt (Python >= 3.10)
    """

    proc = subprocess.Popen(["python", f"{PACKAGE_DIR}/fetchmgs.py", "-h"], stdout=subprocess.PIPE)
    test_out = proc.communicate()[0].decode().split()
    if sys.version_info[1] < 10:
        original_out = open(f"{TEST_INPUT_DIR}/cli_nomode.38.txt", "r").read().split()
    else:
        original_out = open(f"{TEST_INPUT_DIR}/cli_nomode.310.txt", "r").read().split()
    assert test_out[1:] == original_out[1:]


def test_cli_calibration():
    """
    Test for the command line interface output when calibration mode is selected.
    
    If cli arguments are change, regenerate the cli txt files in test data by 
    python3 fetchmgs.py -m calibration -h > ../test/test_input/cli_calibration.38.txt (Python < 3.10)
    python3 fetchmgs.py -m calibration -h > ../test/test_input/cli_calibration.310.txt (Python >= 3.10)
    """

    proc = subprocess.Popen(["python", f"{PACKAGE_DIR}/fetchmgs.py", "-m", "calibration", "-h"], stdout=subprocess.PIPE)
    test_out = proc.communicate()[0].decode().split()
    if sys.version_info[1] < 10:
        original_out = open(f"{TEST_INPUT_DIR}/cli_calibration.38.txt", "r").read().split()
    else:
        original_out = open(f"{TEST_INPUT_DIR}/cli_calibration.310.txt", "r").read().split()
    assert test_out[1:] == original_out[1:]


def test_cli_extraction():
    """
    Test for the command line interface output when extraction mode is selected.
    
    If cli arguments are change, regenerate the cli txt files in test data by 
    python3 fetchmgs.py -m extraction -h > ../test/test_input/cli_extraction.38.txt (Python < 3.10)
    python3 fetchmgs.py -m extraction -h > ../test/test_input/cli_extraction.310.txt (Python >= 3.10)
    """

    proc = subprocess.Popen(["python", f"{PACKAGE_DIR}/fetchmgs.py", "-m", "extraction", "-h"], stdout=subprocess.PIPE)
    test_out = proc.communicate()[0].decode().split()
    if sys.version_info[1] < 10:
        original_out = open(f"{TEST_INPUT_DIR}/cli_extraction.38.txt", "r").read().split()
    else:
        original_out = open(f"{TEST_INPUT_DIR}/cli_extraction.310.txt", "r").read().split()
    assert test_out[1:] == original_out[1:]
