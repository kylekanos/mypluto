import os
import shutil
import sys
import string
import re
import time
import file

######################################################################
def nchange (fini, s, n, repl_str):

# 
# Search for the occurence of string "s" in file "fini".
# Replace the "n"-th word (n=0,1,2...) of the string with 
# a different string "repl_str". This function is invoked 
# frequently by other functions in this module.
#

  nl   = (file.string_find(fini,s))[0]     # find the line number
  line = (file.read_lines(fini,nl,nl))[0]  # read the whole line
  scrh = string.split(line)                # split the words into a list
  scrh[n] = repl_str                       # replace the n-th word
  scrh = string.join(scrh,'   ')           # join words, add space
  scrh = scrh+"\n"                         # add newline
  file.replace_line(fini,scrh,nl)          # overwrite

def xbeg (fini, x):
  nchange(fini, "X1-grid", 2, str(x))

def xend (fini, x):
  nchange(fini, "X1-grid", 5, str(x))

def ybeg (fini, y):
  nchange(fini, "X2-grid", 2, str(y))

def yend (fini, y):
  nchange(fini, "X2-grid", 5, str(y))

def zbeg (fini, z):
  nchange(fini, "X3-grid", 2, str(z))

def zend (fini, z):
  nchange(fini, "X3-grid", 5, str(z))
  
def xres (fini, n):
  nchange(fini, "X1-grid", 3, str(n))

def yres (fini, n):
  nchange(fini, "X2-grid", 3, str(n))

def zres (fini, n):
  nchange(fini, "X3-grid", 3, str(n))

def CFL (fini, cfl):
  scrh = file.string_find (fini, 'CFL')
  file.replace_line (fini, "CFL     "+str(cfl)+"\n", scrh[0])

def CFL_max_var(fini, q):
  scrh = file.string_find (fini, 'CFL_max_var')
  file.replace_line (fini, "CFL_max_var     "+str(q)+"\n", scrh[0])

def tstop (fini, t):
  scrh = file.string_find (fini, 'tstop')
  file.replace_line (fini, "tstop     "+str(t)+"\n", scrh[0])
  
def first_dt(fini, dt):
  scrh = file.string_find (fini, 'first_dt')
  file.replace_line (fini, "first_dt      "+str(dt)+"\n", scrh[0])

def analysis (fini, xf,nf):
  scrh = file.string_find (fini, 'analysis')
  file.replace_line (fini, "analysis    "+str(xf)+"    "+str(nf)+"\n", scrh[0])

def dbl (fini, xf, nf, mode):
  scrh = file.string_find (fini, 'dbl')
  file.replace_line (fini, "dbl    "+str(xf)+"    "+str(nf)+"    "+mode+"\n", scrh[0])

def flt (fini, xf, nf, mode):
  scrh = file.string_find (fini, 'flt')
  file.replace_line (fini, "flt    "+str(xf)+"    "+str(nf)+"    "+mode+"\n", scrh[0])

def interpolation(fdef,interp):
  scrh = file.string_find(fdef, 'INTERPOLATION')
  file.replace_line(fdef, '#define    INTERPOLATION    '+interp+'\n',scrh[0]);
  
  
  
  
  
  



