# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 11:18:25 2016

@author: meekBSD
"""

import argparse
import sys

parser = argparse.ArgumentParser(description="The USAGE")     #simplifys the wording of using argparse as stated in the python tutorial
parser.add_argument("-i", "--input", action = 'store', dest='testName', required = True , help="input a string passing to testName") # allows input of the file for reading
parser.add_argument("-m", "--may_day", type = int)
parser.add_argument("-cn", "--cook_num", type = int)

try:
    options = parser.parse_args()
    if options.testName == None:
        parser.print_help()
        sys.exit(0)
except:
    sys.exit(0)

print (options.testName)
print ("Your second argument is %d, it multiplys 6 equals to %d"%(options.may_day, options.may_day *6))
print ("Your third argument is %d, it multiply 10 equals to %d"%(options.cook_num,options.cook_num *10))
print ("\n\nDone.")
