import os
import shutil
import sys
import string
import re
import time
import file
import ut
import menu

#
# Input to conversion are:
# 
#  dbl/flt,     : select via argv
#  little/big,  : pass to vtk function
#  single/multi : 
#  num. of scalars:
#  num. of vectors:
#  rec_beg:
#  rec_end:
#    
# Input to VTK:
#
#
# size,
# number of scalars,
# input file names and offsets
# geometry
# endianity 
# float/double
#     

entries = []
options = []
default = []
scrh    = []



outname = "flt.out"
line = file.read_lines (outname,0,0)

rec_beg = 0
rec_end = file.count_lines (outname)

scrh = []
entries.append('REC_BEG')
for n in range(rec_end): scrh.append(repr(n))
options.append(scrh)
default.append('0')
 
scrh = []
entries.append('REC_END')
for n in range(rec_end): scrh.append(repr(n))
options.append(scrh)
default.append(str(rec_end-1))

scrh = []
entries.append('REC_SKIP')
for n in range(rec_end): scrh.append(repr(n+1))
options.append(scrh)
default.append('1')

scrh = []
scrh = string.split(line[0])
for x in scrh[6:]: 

  vvect = 0
  bvect = 0
  y = []
  if (x == "v1"):
    i = scrh.index(x)
    vvect = 1
    if (scrh[i+1] == "v2"):
      i = i + 1
      vvect = 2 
      if (scrh[i+1] == "v3"):
        vvect = 3

  if (x == "b1"):
    i = scrh.index(x)
    bvect = 1
    if (scrh[i+1] == "b2"):
      i = i + 1
      bvect = 2 
      if (scrh[i+1] == "b3"):
        bvect = 3

  if (x == "v2" or x == "b2"): continue
  if (x == "v3" or x == "b3"): continue

  if (vvect == 0 and bvect == 0): 
    entries.append(x)
    options.append(['YES','NO'])
    default.append('YES')
  elif (vvect == 2):
    entries.append('Velocity')
    options.append(['YES','NO'])
    default.append('YES')
  elif (bvect == 2):
    entries.append('Magnetic field')
    options.append(['YES','NO'])
    default.append('YES')




print entries
selection = '' 
selection = menu.select(entries,options,default,'Input')
