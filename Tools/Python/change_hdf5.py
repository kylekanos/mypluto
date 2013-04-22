import os
import shutil
import sys
import string
import re
import time
import file
import ut
import menu

# Script to automatically change Checkpoint_interval in pluto.ini

if (len(sys.argv) == 1): fname = "pluto.ini"
else:                    fname = sys.argv[1]

for fname in sys.argv[1:]:
  print "doing "+fname
  try:
    scrh = file.string_find (fname, 'Checkpoint_interval')
    file.replace_line (fname, "Checkpoint_interval  -1.0  0\n", scrh[0])

    scrh = file.string_find (fname, 'Plot_interval')
    file.replace_line (fname, "Plot_interval         1.0  0\n", scrh[0])
  except:
    continue
  
  
  
  



