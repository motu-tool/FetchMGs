#!/usr/bin/env python

#############################################################
#
# test_mycode.py
#
# Usage:
# - call 'pytest' from the directory were this file is located
#
#############################################################

import mycode as code     # import functions from mycode.py

def test_multiply():
    """ Test our multiplication function against the standard addition operator """
    assert code.multiply(3,5) == 15     # looking for a specific output
    
    # Test with extended coverage of the function 
    for x in range(10):
        assert code.multiply(3,x) == x + x + x