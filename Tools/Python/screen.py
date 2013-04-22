"""               screen.py   
 
 Useful utilities for flushing output to screen
 A. Mignone (mignone@to.astro.it)

""" 

import string
import os
import sys

###########################################################
def dump(string):

# 
#  Dump string to screen
#

  sys.stdout.write (string)
  sys.stdout.flush()

###########################################################
def dumpclear(string):

# 
#  Clear screen and dump string to screen
#

  os.system("clear")
  sys.stdout.write (string)
  sys.stdout.flush()

