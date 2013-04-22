import os
import shutil
import sys
import string
import re
import time
import file
import ut
import menu

# Script to automatically change pluto.ini from v. 2 to v. 3.


if (len(sys.argv) == 1): fname = "pluto.ini"
else:                    fname = sys.argv[1]

# change headers first

scrh = file.string_find (fname, 'DOMAIN')
file.replace_line (fname, "[Grid]\n\n", scrh[0])

scrh = file.string_find (fname, 'TIME')
file.replace_line (fname, "[Time]\n", scrh[0])

scrh = file.string_find (fname, 'SOLVER')
file.replace_line (fname, "[Solver]\n", scrh[0])

scrh = file.string_find (fname, 'BOUNDARIES')
file.replace_line (fname, "[Boundary]\n", scrh[0])

scrh = file.string_find (fname, 'OUTPUT')
file.replace_line (fname, "[Uniform Grid Output]\n", scrh[0])

scrh = file.string_find (fname, 'USER-DEF')
file.replace_line (fname, "[Parameters]\n", scrh[0])

# correct old labels - 1

old_label = ["X-domain","Y-domain","Z-domain", "CFL_Max_Var", 
             "X-Left","X-Right", "Y-Left","Y-Right","Z-Left","Z-Right"]
new_label = ["X1-grid" ,"X2-grid", "X3-grid", "CFL_max_var", 
             "X1-beg", "X1-end", "X2-beg", "X2-end","X3-beg", "X3-end"]

for label in old_label:
  indx = old_label.index(label)
  scrh = file.string_find (fname, label)
  if len(scrh) == 0: continue
  ipos = scrh[0]
  line = file.read_lines (fname, ipos, ipos+1)
  line_list = string.split(line[0])
  corrected_line = new_label[indx]

  for scrh in line_list[1:]:
    corrected_line = corrected_line + '    '+scrh

  corrected_line += '\n'
  file.replace_line (fname,corrected_line, ipos)


# correct old labels - 2
# some old pluto.ini still have "X1-low" in the boundary 
# section

old_label = ["X1-low","X1-up","X2-low", "X2-up", "X3-low","X3-up"]
new_label = ["X1-beg", "X1-end", "X2-beg", "X2-end","X3-beg", "X3-end"]

for label in old_label:
  indx = old_label.index(label)
  scrh = file.string_find (fname, label)
  if len(scrh) == 0: continue
  ipos = scrh[0]
  line = file.read_lines (fname, ipos, ipos+1)
  line_list = string.split(line[0])
  corrected_line = new_label[indx]

  for scrh in line_list[1:]:
    corrected_line = corrected_line + '    '+scrh

  corrected_line += '\n'
  file.replace_line (fname,corrected_line, ipos)

# correct output section

scrh = file.string_find (fname, 'tsave')
if (len(scrh) == 0): scrh = file.string_find(fname, 'time')

ipos = scrh[0]
line = file.read_lines (fname, ipos, ipos)
line_list = string.split(line[0])
tsave = line_list[1]

file.delete_lines (fname, ipos-1,ipos+2)

lines = ["uservar    0\n",
         "dbl       "+tsave+"  -1   single_file\n",
         "flt       -1.0  -1   single_file\n",
         "vtk       -1.0  -1   single_file\n",
         "tab       -1.0  -1   \n",
         "ppm       -1.0  -1   \n",
         "png       -1.0  -1   \n",
         "log        1 \n", 
         "analysis  -1.0  -1 \n", " \n"]

ipos = ipos - 1
for x in lines:
  file.insert (fname, x, ipos)
  ipos = ipos + 1


# add Chombo refinement section

scrh = file.string_find (fname, '[Time]')
ipos = scrh[0]

lines = ["[Chombo Refinement]\n", " \n", 
         "Levels           4\n", 
         "Ref_ratio        2 2 2 2 2\n",
         "Regrid_interval  2 2 2 2\n",
         "Refine_thresh    0.3\n",
         "Tag_buffer_size  3\n", 
         "Block_factor     4\n",
         "Max_grid_size    32\n", 
         "Fill_ratio       0.75\n", " \n"]
	  
for x in lines:
  file.insert(fname, x, ipos)
  ipos = ipos + 1


# add Chombo HDF5 output section

#scrh = file.string_find (fname, '[Time]')
#ipos = scrh[0]

scrh = file.string_find (fname, '[Parameters]')
ipos = scrh[0]

lines = ["[Chombo HDF5 output]\n", " \n",
         "Checkpoint_interval  -1.0  0\n",
         "Plot_interval         1.0  0\n", " \n" ]

for x in lines:
  file.insert(fname, x, ipos)
  ipos = ipos + 1
	 




print "Successfully converted"


