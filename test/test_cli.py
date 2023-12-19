import sys, os
import pytest
import subprocess
from fetchmgs.test_config import *

# Load fetchMGs functions
sys.path.insert(0, PACKAGE_DIR) 
from fetchMGs import cli

def test_cli_no_mode():
    """
    Test for the command line interface output when no mode is defined 
    
    If cli arguments are change, regenerate the cli txt files in test data by 
    python3 fetchMGs.py -h > ../test/test_data/cli_no_mode.txt
    """

    proc = subprocess.Popen(["python", f"{PACKAGE_DIR}/fetchMGs.py", "-h"], stdout=subprocess.PIPE)
    test_out = proc.communicate()[0].decode().split()
    original_out = open(f"{TEST_INPUT_DIR}/cli_no_mode.txt", "r").read().split()
    assert test_out == original_out


def test_cli_calibration():
    """
    Test for the command line interface output when calibration mode is selected.
    
    If cli arguments are change, regenerate the cli txt files in test data by 
    python3 fetchMGs.py -m calibration -h > ../test/test_data/cli_calibration.txt
    """

    proc = subprocess.Popen(["python", f"{PACKAGE_DIR}/fetchMGs.py", "-m", "calibration", "-h"], stdout=subprocess.PIPE)
    test_out = proc.communicate()[0].decode().split()
    original_out = open(f"{TEST_INPUT_DIR}/cli_calibration.txt", "r").read().split()
    assert test_out == original_out


def test_cli_extraction():
    """
    Test for the command line interface output when extraction mode is selected.
    
    If cli arguments are change, regenerate the cli txt files in test data by 
    python3 fetchMGs.py -m extraction -h > ../test/test_data/cli_extraction.txt
    """

    proc = subprocess.Popen(["python", f"{PACKAGE_DIR}/fetchMGs.py", "-m", "extraction", "-h"], stdout=subprocess.PIPE)
    test_out = proc.communicate()[0].decode().split()
    original_out = open(f"{TEST_INPUT_DIR}/cli_extraction.txt", "r").read().split()
    assert test_out == original_out
