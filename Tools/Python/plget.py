import os
import shutil
import sys
import string
import re
import time
import file


def xres (fini):
  dimstr = file.string_list(fini,"X1-grid")
  dimstr = string.split(dimstr[0])
  return string.atoi(dimstr[3])


def xend (fini):
  dimstr = file.string_list(fini,"X1-grid")
  dimstr = string.split(dimstr[0])
  return string.atof(dimstr[5])


def DIMENSIONS (fdef):
  dimstr = file.string_list(fdef,"DIMENSIONS")
  dimstr = string.split(dimstr[0])
  return string.atoi(dimstr[2])
 


  
  
  
  



