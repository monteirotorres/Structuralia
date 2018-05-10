# License
###############################################################################
'''
Structuralia: A suite of python scripts to easily manipulate PDBs
Copyright (C) 2018  Pedro H. M. Torres

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact info:
Department Of Biochemistry
University of Cambridge
80 Tennis Court Road
Cambridge CB2 1GA
E-mail address: monteirotorres@gmail.com
'''

# Description
###############################################################################
'''
This is a simple function that can be imported and used as a wrapper to time
the execution of any wrapped function.

Developed by: Pedro Torres
'''
# Imports
###############################################################################
from functools import wraps
from time import time


# Functions
###############################################################################

def timed(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        start = time()
        result = f(*args, **kwds)
        elapsed = time() - start
        # timelog = open('timelog.txt', 'a+')
        # timelog.write("{0} \t {1:5.3f} seconds \n".format(f.__name__, elapsed))
        print("{0} \t {1:5.5f} seconds \n".format(f.__name__, elapsed))
        return result
    return wrapper
