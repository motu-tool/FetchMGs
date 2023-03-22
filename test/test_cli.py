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

    proc = subprocess.Popen(["python", "-c", "import writer; writer.write()"], stdout=subprocess.PIPE)
    out = proc.communicate()[0]
    print out.upper()