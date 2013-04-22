import os
import shutil
import sys
import string
import re
import time
import file
import menu


print sys.argv
for x in sys.argv:
  if (x == "--with-chombo:"): 
    i = sys.argv.index(x) + 1
    a = repr(2)
    for y in sys.argv[i:]: a = a+' '+y
    print "do "+a





