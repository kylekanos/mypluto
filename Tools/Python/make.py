import os
import shutil
import sys
import string
import re
import time
import file
import ut
import menu

###############################################################
#
#              B U I L D       P R O B L E M
#
###############################################################
def problem(work_dir, pluto_path, pluto_dir, additional_files, 
            additional_flags, makefile_template, AUTO_UPDATE):

#  before proceeding we scan pre-existing definitions.h 
#  for obsolete statements in order to and make  proper replacements.

  replace_obsoletes(work_dir)
 
#  sys.exit(1)

  HD     = 1
  RHD    = 2
  MHD    = 3
  SRMHD  = 4
  RMHD   = 5
  
  WITH_CHOMBO = 0  # enable/disable Chombo AMR
  FULL        = 0
  WITH_FD     = 0  # enable/disable finite difference schemes 
  WITH_SB     = 0  # enable/disable shearing box module
  WITH_FARGO  = 0  # enable/disable FARGO module for orbital advection

# *****************************************
#      check command line arguments   
# *****************************************

  for x in sys.argv:
    if (x == "--with-chombo" or x == "--with-chombo:"):
      WITH_CHOMBO = 1
    if (x == "--full"):
      FULL = 1
    if (x == "--with-fd"):
      WITH_FD = 1
    if (x == "--with-sb"):
      WITH_SB = 1
    if (x == "--with-fargo"):
      WITH_FARGO = 1
        
# -- default makefile template -- 

  if (WITH_CHOMBO):
    makefile_template.append('/Src/Templates/makefile.chombo')
  else:
    makefile_template.append('/Src/Templates/makefile')
  
# Create user-editable menu

  entries = []
  options = []
  default = []
  scrh    = []

  entries.append('PHYSICS')
  if (WITH_FD or WITH_FARGO):
    options.append(['HD','MHD'])
    default.append('HD')
  elif (WITH_SB):
    options.append(['MHD'])
    default.append('MHD')
  else:
    options.append(['HD','RHD','MHD','RMHD'])
    default.append('HD')

  entries.append('DIMENSIONS')
  if(WITH_SB or WITH_FARGO):
    options.append(['2','3'])
    default.append('2')
  else:
    options.append(['1','2','3'])
    default.append('1')


  entries.append('COMPONENTS')
  if (WITH_SB):
    options.append(['2','3'])
    default.append('2')
  else:
    options.append(['1','2','3'])
    default.append('1')

  entries.append('GEOMETRY')
  if (WITH_CHOMBO):
#    options.append(['CARTESIAN','CYLINDRICAL','AXISYM'])
    options.append(['CARTESIAN','CYLINDRICAL'])
    default.append('CARTESIAN')
  elif (WITH_FD or WITH_SB):
    options.append(['CARTESIAN'])
    default.append('CARTESIAN')
  else:
#    options.append(['CARTESIAN','CYLINDRICAL','AXISYM','POLAR','SPHERICAL'])
    options.append(['CARTESIAN','CYLINDRICAL','POLAR','SPHERICAL'])
    default.append('CARTESIAN')

#  entries.append('UNIFORM_GRID')
#  if (WITH_CHOMBO):
#    options.append('YES')
#    default.append('YES')
#  else:
#    options.append(['YES','NO'])
#    default.append('YES')

    
  entries.append('BODY_FORCE')
  if (WITH_SB): 
    options.append(['VECTOR', 'POTENTIAL', '(VECTOR+POTENTIAL)'])
    default.append('VECTOR')
  else:
    options.append(['NO','VECTOR', 'POTENTIAL', '(VECTOR+POTENTIAL)'])
    default.append('NO')

  entries.append('COOLING')
  options.append(['NO','POWER_LAW','TABULATED','SNEq','MINEq']) # ,'H2_COOL'])
  default.append('NO')

#  entries.append('INCLUDE_PARTICLES')
#  options.append(['NO','YES'])
#  default.append('NO')

# set interpolation options

  entries.append('INTERPOLATION')
  if (WITH_CHOMBO):
    options.append(['FLAT','LINEAR','WENO3','PARABOLIC'])
    default.append('LINEAR')
  elif (FULL):
    options.append(['FLAT','LINEAR','LimO3',
                    'WENO3','PARABOLIC', 'MP5'])
    default.append('LINEAR')
  elif (WITH_FD):
    options.append(['WENO3_FD', 'WENOZ_FD', 
                    'MP5_FD','LIMO3_FD'])
    default.append('WENOZ_FD')
  else: 
    options.append(['FLAT','LINEAR','LimO3','WENO3','PARABOLIC'])  # Default
    default.append('LINEAR')


  entries.append('TIME_STEPPING')
  if (WITH_CHOMBO):
#    options.append(['HANCOCK','CHARACTERISTIC_TRACING','RK_MIDPOINT','EULER','RK2'])
    options.append(['EULER','HANCOCK','CHARACTERISTIC_TRACING','RK2'])
    default.append('HANCOCK')
  elif (WITH_FD):
    options.append(['RK3'])
    default.append('RK3')
  else:    
#    options.append(['EULER','RK_MIDPOINT','RK2','RK3','HANCOCK','CHARACTERISTIC_TRACING'])
    options.append(['EULER','RK2','RK3','HANCOCK','CHARACTERISTIC_TRACING'])
    default.append('RK2')

  entries.append('DIMENSIONAL_SPLITTING')
  if (WITH_CHOMBO or WITH_FARGO):
    options.append(['NO'])
    default.append('NO')
  else:  
    options.append(['YES','NO'])
    default.append('YES')

  entries.append('NTRACER')
  for n in range(5): scrh.append(repr(n))
  options.append(scrh)
  default.append('0')

  scrh =[]
  entries.append('USER_DEF_PARAMETERS')
  for n in range(32): scrh.append(repr(n+1))
  options.append(scrh)
  default.append('1')

# **************************************************
#   If  definitions.h already exists try to read 
#   default values from it. 
# **************************************************

  if (os.path.exists(work_dir+'/definitions.h')):
    for x in entries:
      try:
        scrh = file.string_list(work_dir+'/definitions.h',x)
        tmp  = string.split(scrh[0])
        i  = entries.index(tmp[1])
        y  = options[i]
        i2 = y.index(tmp[2]) 
        default[entries.index(x)] = tmp[2]
  
      except ValueError:
        continue

      except IndexError:
        continue

#  ****  ok, call menu  ****

  if AUTO_UPDATE == 0:
    selection = '' 
    menu.SetTitle("Setup problem")
    selection = menu.Browse(entries,default,options)
#    selection = menu.select(entries,options,default,'Setup problem')

# **************************************************
#       define Physics module sub-menus
# **************************************************

  i              = entries.index('DIMENSIONS')
  dimensions     = int(default[i])
  i              = entries.index('PHYSICS')
  physics_module = default[i]
  i              = entries.index('GEOMETRY')
  geometry       = default[i]
  i              = entries.index('INTERPOLATION')
  interpolation  = default[i]
  i              = entries.index('TIME_STEPPING')
  time_stepping  = default[i]
  i              = entries.index('DIMENSIONAL_SPLITTING')
  dim_splitting  = default[i]

# HD

  if physics_module == 'HD':
    entries_HD = []
    options_HD = []
    default_HD = []

    entries_HD.append('EOS')
    options_HD.append(['IDEAL','ISOTHERMAL'])
    default_HD.append('IDEAL')

    entries_HD.append('ENTROPY_SWITCH')
    options_HD.append(['NO','YES'])
    default_HD.append('NO')

    entries_HD.append('THERMAL_CONDUCTION')
    if (WITH_CHOMBO == 1):
      options_HD.append(['NO','EXPLICIT'])
    else:
#      options_HD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_CHEBYSHEV'])
      options_HD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING'])
      
    default_HD.append('NO')

    entries_HD.append('VISCOSITY')
    if (WITH_CHOMBO == 1):
      options_HD.append(['NO','EXPLICIT'])
    else:
#      options_HD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_CHEBYSHEV'])
      options_HD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING'])
      
    default_HD.append('NO')

    entries_HD.append('ROTATING_FRAME')
    options_HD.append(['NO','YES'])
    default_HD.append('NO')


#    if (physics_module == "HD" and dim_splitting == "YES"):
#      if (geometry == "POLAR" or geometry == "SPHERICAL"):
#        entries_HD.append('FARGO_SCHEME')
#        options_HD.append(['NO','YES'])
#        default_HD.append('NO')

#  Read HD pre-existing HD submenu defaults, if they exists

    if (os.path.exists(work_dir+'/definitions.h')):
      for x in entries_HD:
        try:
          scrh = file.string_list(work_dir+'/definitions.h',x)
          tmp  = string.split(scrh[0])
          default_HD[entries_HD.index(x)] = tmp[2]
        except IndexError:
          continue


    if AUTO_UPDATE == 0:
      selection = ''
      menu.SetTitle("HD Menu")
      selection = menu.Browse(entries_HD, default_HD, options_HD)
#      selection = menu.select(entries_HD,options_HD,default_HD,'HD MENU')

# RHD

  if physics_module == 'RHD':
    entries_RHD = []
    options_RHD = []
    default_RHD = []

    entries_RHD.append('EOS')
    options_RHD.append(['IDEAL','TAUB'])
    default_RHD.append('IDEAL')

    entries_RHD.append('USE_FOUR_VELOCITY')
    options_RHD.append(['NO','YES'])
    default_RHD.append('NO')

    entries_RHD.append('ENTROPY_SWITCH')
    options_RHD.append(['NO','YES'])
    default_RHD.append('NO')

#  Read RHD pre-existing RHD submenu defaults, if they exists

    if (os.path.exists(work_dir+'/definitions.h')):
      for x in entries_RHD:
        try:
          scrh = file.string_list(work_dir+'/definitions.h',x)
          tmp  = string.split(scrh[0])
          default_RHD[entries_RHD.index(x)] = tmp[2]
        except IndexError:
          continue


    if AUTO_UPDATE == 0:
      selection = ''
      menu.SetTitle("RHD Menu")
      selection = menu.Browse(entries_RHD, default_RHD, options_RHD)
#      selection = menu.select(entries_RHD,options_RHD,default_RHD,'RHD MENU')

# MHD

  if (physics_module == 'MHD'):
    entries_MHD = []
    options_MHD = []
    default_MHD = []

    entries_MHD.append('EOS')
    options_MHD.append(['IDEAL','BAROTROPIC','ISOTHERMAL'])
    default_MHD.append('IDEAL')

    entries_MHD.append('ENTROPY_SWITCH')
    options_MHD.append(['NO','YES'])
    default_MHD.append('NO')

    entries_MHD.append('MHD_FORMULATION')
    if (WITH_CHOMBO or WITH_FD):
      options_MHD.append(['NONE','EIGHT_WAVES','DIV_CLEANING'])
      default_MHD.append('EIGHT_WAVES')
    elif (WITH_SB or WITH_FARGO):
      options_MHD.append(['CONSTRAINED_TRANSPORT'])
      default_MHD.append('CONSTRAINED_TRANSPORT')
    else:
      options_MHD.append(['NONE','EIGHT_WAVES','DIV_CLEANING', 'CONSTRAINED_TRANSPORT'])
      default_MHD.append('EIGHT_WAVES')


    entries_MHD.append('BACKGROUND_FIELD')
    options_MHD.append(['NO','YES'])
    default_MHD.append('NO')

    entries_MHD.append('RESISTIVE_MHD')
    if (WITH_CHOMBO == 1):
      options_MHD.append(['NO','EXPLICIT'])
    else:
#      options_MHD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_CHEBYSHEV'])
      options_MHD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING'])

    default_MHD.append('NO')

    entries_MHD.append('THERMAL_CONDUCTION')
    if (WITH_CHOMBO == 1):
      options_MHD.append(['NO','EXPLICIT'])
    else:
#      options_MHD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_CHEBYSHEV'])
      options_MHD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING'])

    default_MHD.append('NO')

    entries_MHD.append('VISCOSITY')
    if (WITH_CHOMBO == 1):
      options_MHD.append(['NO','EXPLICIT']) 
    else:   
#      options_MHD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_CHEBYSHEV'])
      options_MHD.append(['NO','EXPLICIT','SUPER_TIME_STEPPING'])
    
    default_MHD.append('NO')

    entries_MHD.append('ROTATING_FRAME')
    options_MHD.append(['NO','YES'])
    default_MHD.append('NO')

#  Read MHD pre-existing MHD submenu defaults, if they exists
 
    if (os.path.exists(work_dir+'/definitions.h')):
      for x in entries_MHD:
        ix = entries_MHD.index(x)
        try:
          scrh = file.string_list(work_dir+'/definitions.h',x)
          tmp  = string.split(scrh[0])
          try:   # if the choice is not in the permitted list of options,
                 # ValueError is raised and the default choice is kept.
            i = options_MHD[ix].index(tmp[2])
            default_MHD[ix] = tmp[2]
          except ValueError:
            continue
        except IndexError:
          continue

    if AUTO_UPDATE == 0:
      selection = ''
      if (physics_module == 'MHD'):
        menu.SetTitle("MHD Menu")
        selection = menu.Browse(entries_MHD, default_MHD, options_MHD)
#        selection = menu.select(entries_MHD,options_MHD,default_MHD,'MHD MENU')


# RMHD

  if physics_module == 'RMHD':
    entries_RMHD = []
    options_RMHD = []
    default_RMHD = []

    entries_RMHD.append('EOS')
    options_RMHD.append(['IDEAL','TAUB'])
    default_RMHD.append('IDEAL')

    entries_RMHD.append('ENTROPY_SWITCH')
    options_RMHD.append(['NO','YES'])
    default_RMHD.append('NO')

    entries_RMHD.append('MHD_FORMULATION')
    if (WITH_CHOMBO == 1):
      options_RMHD.append(['NONE','EIGHT_WAVES','DIV_CLEANING'])
    else:
      options_RMHD.append(['NONE','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT'])
 
    default_RMHD.append('NONE')

#  Read RMHD pre-existing RMHD submenu defaults, if they exists

    if (os.path.exists(work_dir+'/definitions.h')):
      for x in entries_RMHD:
        ix = entries_RMHD.index(x)
        try:
          scrh = file.string_list(work_dir+'/definitions.h',x)
          tmp  = string.split(scrh[0])
          try:   # if the choice is not in the permitted list of options,
                 # ValueError is raised and the default choice is kept.
            i = options_RMHD[ix].index(tmp[2])
            default_RMHD[ix] = tmp[2]
          except ValueError:
            continue
        except IndexError:
          continue

    if AUTO_UPDATE == 0:
      selection = ''
      if (physics_module == 'RMHD'):
        menu.SetTitle("RMHD Menu")
        selection = menu.Browse(entries_RMHD, default_RMHD, options_RMHD)
#        selection = menu.select(entries_RMHD,options_RMHD,default_RMHD,'RMHD MENU')

# ********************************************************************
#    USER_DEF PARAMETERS:
#
#      Read them from definitions.h if the file exists; 
#      assign  'SCRH' otherwise.
# ********************************************************************

  i    = entries.index('USER_DEF_PARAMETERS')
  npar = int(default[i])
  u_def_par = []
  if (os.path.exists(work_dir+'/definitions.h')):
    scrh = file.word_find(work_dir+'/definitions.h','parameters')
    k0   = scrh[0] + 2
    scrh = file.read_lines(work_dir+'/definitions.h', k0, k0 + npar)
    for n in range(npar):
      try:
        x = string.split(scrh[n])
        if (x[0] == "#define"):
          u_def_par.append(x[1])
        else:          
          u_def_par.append('SCRH')
      except IndexError:
        u_def_par.append('SCRH')
#        scrh[n + 1] = '#define SCRH  x'


  else:
    for n in range(npar):
      u_def_par.append('SCRH')


  if AUTO_UPDATE == 0:
#    u_def_par = ['SCRH']
    menu.SetTitle ("User-defined parameters")
    menu.Insert(int(default[i]),u_def_par)
#    menu.put(int(default[i]),u_def_par,'USER DEF SECTION')

# -----------------------------------------------------
#   Create the list 'tmp' that will be used to write
#   the new definitions.h header file
# -----------------------------------------------------

  tmp = []
  for x in entries:
    i = entries.index(x)
    y = default[i]
    tmp.append('#define  '+x.ljust(21)+'   '+y+'\n')

  tmp.append('\n/* -- physics dependent declarations -- */\n\n');

  if (physics_module == 'MHD'):
    for x in entries_MHD:
      i = entries_MHD.index(x)
      tmp.append('#define  '+x.ljust(21)+'   '+default_MHD[i]+'\n')

  if physics_module == 'HD':
    for x in entries_HD:
      i = entries_HD.index(x)
      tmp.append('#define    '+x.ljust(21)+'   '+default_HD[i]+'\n')

  if physics_module == 'RHD':
    for x in entries_RHD:
      i = entries_RHD.index(x)
      tmp.append('#define    '+x.ljust(21)+'   '+default_RHD[i]+'\n')

  if physics_module == 'RMHD':
    for x in entries_RMHD:
      i = entries_RMHD.index(x)
      tmp.append('#define    '+x.ljust(21)+'   '+default_RMHD[i]+'\n')

# *****************************************
#     add Chombo-AMR specific flags 
# *****************************************

  if (WITH_CHOMBO and os.path.exists(work_dir+'/definitions.h')):
    scrh = file.string_list(work_dir+'/definitions.h','AMR_EN_SWITCH')
    tmp.append('\n/* -- Chombo-AMR flags -- */\n\n');
    if (len(scrh) != 0): 
      strsplit = string.split(scrh[0])
      value    = strsplit[2]
      tmp.append('#define  AMR_EN_SWITCH   '+value+'\n')
    else:
      tmp.append('#define  AMR_EN_SWITCH   NO\n')

# *******************************************************
#     add pointer names to user defined parameters
# *******************************************************

  tmp.append('\n/* -- pointers to user-def parameters -- */\n\n');

  for x in u_def_par:
    i = u_def_par.index(x)
    tmp.append('#define  '+x.ljust(16)+'   '+repr(i)+'\n')

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #   In this section, also, udpate the pluto.ini 
  #   initialization file by appending the correct 
  #   sequence of user defined parameters
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#  menu.RestoreScreen()
  plutoini = work_dir+"/pluto.ini"
  scrh = file.string_find (plutoini, 'Parameters')
  ipos = scrh[0] + 2  # start reading/writing pluto.ini at this point  
                      # to build a list of the userdef parameters
  scrh = file.read_lines(plutoini,ipos,999) # read until end of file
    
  pname  = []
  pvalue = []

  for x in scrh:
    y = string.split(x)  # pname and pvalue contain the list of parameter
    if (len(y) == 0): continue  # skip blank lines 
    pname.append (y[0])  # names and their respective values as present
    pvalue.append(y[1])  # in pluto.ini

  for x in u_def_par:    # loop on the actual list of parameters;
    par_exist = 0        # is the parameter already present in pluto.ini ?
    for y in pname:
      if (x == y):       # if yes, use the original value instead of '0'
        par_exist = 1
        i = pname.index(y)
        file.insert (plutoini, pname[i] + "    " + pvalue[i] + "\n", ipos)
        break

    if (par_exist == 0):  # if not, set default value to 0
      file.insert (plutoini, x+'    0\n',ipos)

    ipos = ipos + 1

  # delete all other lines in pluto.ini

  file.delete_lines(plutoini, ipos, ipos + 256) 

# ***************************************************************
#    check dependencies for conditional inclusion
#    of additional object files; this is useful 
#    for makefile creation
# ***************************************************************

# additional files will be compiled 
# depending on the user's choice of
# time integration scheme
 
# the next line removes all elements from additional_files
# and pluto_path so the list can be built from scratch

#  for x in additional_files: additional_files.remove(x)
#  for x in pluto_path      : pluto_path.remove(x)

  additional_files[:] = []
  additional_flags[:] = []
  pluto_path[:] = []

  divb = ' '

  # **********************************************
  #   physics module & parabolic term treatment
  # **********************************************

  sts_flag = 0 # if changed to 1, include the super-time-stepping driver
  rkc_flag = 0 # if changed to 1, include the runge-kutta-Chebyshev driver
  exp_flag = 0 # if changed to 1, include explicit parabolic time-stepping

  if (physics_module == 'HD'):
    pluto_path.append('HD/')
    entries_MOD = entries_HD
    default_MOD = default_HD

    n = entries_MOD.index('THERMAL_CONDUCTION')
    if (default_MOD[n] != 'NO'):
      pluto_path.append('Thermal_Conduction/')
      if   (default_MOD[n] == 'SUPER_TIME_STEPPING'): sts_flag = 1
      elif (default_MOD[n] == 'RK_CHEBYSHEV'):        rkc_flag = 1
      elif (default_MOD[n] == 'EXPLICIT'):            exp_flag = 1

    n = entries_MOD.index('VISCOSITY')
    if (default_MOD[n] != 'NO'):
      pluto_path.append('Viscosity/')
      if   (default_MOD[n] == 'SUPER_TIME_STEPPING'): sts_flag = 1
      elif (default_MOD[n] == 'RK_CHEBYSHEV'):        rkc_flag = 1
      elif (default_MOD[n] == 'EXPLICIT'):            exp_flag = 1
      
  elif (physics_module == 'RHD'):
    pluto_path.append('RHD/')
    entries_MOD = entries_RHD
    default_MOD = default_RHD

  elif (physics_module == 'MHD'):
    pluto_path.append('MHD/')
    entries_MOD = entries_MHD
    default_MOD = default_MHD
    n    = entries_MOD.index('MHD_FORMULATION')
    divb = default_MOD[n]
    if (divb == 'CONSTRAINED_TRANSPORT'):
      pluto_path.append('MHD/CT/')
    elif (divb == 'DIV_CLEANING'):
      pluto_path.append('MHD/GLM/')

    n = entries_MOD.index('RESISTIVE_MHD')
    if (default_MOD[n] != 'NO'):
      pluto_path.append('MHD/Resistive/')
      if   (default_MOD[n] == 'SUPER_TIME_STEPPING'): sts_flag = 1
      elif (default_MOD[n] == 'RK_CHEBYSHEV'):        rkc_flag = 1
      elif (default_MOD[n] == 'EXPLICIT'):            exp_flag = 1

    n = entries_MOD.index('THERMAL_CONDUCTION')
    if (default_MOD[n] != 'NO'):
      pluto_path.append('Thermal_Conduction/')
      if   (default_MOD[n] == 'SUPER_TIME_STEPPING'): sts_flag = 1
      elif (default_MOD[n] == 'RK_CHEBYSHEV'):        rkc_flag = 1
      elif (default_MOD[n] == 'EXPLICIT'):            exp_flag = 1
      
    n = entries_MOD.index('VISCOSITY')
    if (default_MOD[n] != 'NO'):
      pluto_path.append('Viscosity/')
      if   (default_MOD[n] == 'SUPER_TIME_STEPPING'): sts_flag = 1
      elif (default_MOD[n] == 'RK_CHEBYSHEV'):        rkc_flag = 1
      elif (default_MOD[n] == 'EXPLICIT'):            exp_flag = 1
  
    # ***************************
    #  need shearingbox module ?
    # ***************************

    if (WITH_SB == 1):
      pluto_path.append('MHD/ShearingBox/')
      additional_flags.append(' -DSHEARINGBOX')
       
  
  elif (physics_module == 'RMHD'):
    pluto_path.append('RMHD/')
    entries_MOD = entries_RMHD
    default_MOD = default_RMHD
    n    = entries_MOD.index('MHD_FORMULATION')
    divb = default_MOD[n]
    if (divb == 'CONSTRAINED_TRANSPORT'):
      pluto_path.append('MHD/CT/')
    elif (divb == 'DIV_CLEANING'):
      pluto_path.append('MHD/GLM/')

  # *****************************
  #       Interpolation   
  # *****************************

  if (interpolation == 'FLAT'):
    additional_files.append('states_flat.o')
  elif   (interpolation == 'LINEAR'):
    additional_files.append('states_plm.o')
  elif (interpolation == 'PARABOLIC'):
    additional_files.append('states_ppm.o')
  elif (interpolation == 'LimO3'):
    additional_files.append('states_limo3.o')
  elif (interpolation == 'WENO3'):
    additional_files.append('states_weno3.o')

  elif (WITH_FD):    # Finite Difference

    additional_files.append('states_fd.o')
    additional_files.append('fd_reconstruct.o')
    additional_files.append('fd_flux.o')
    additional_flags.append(' -DFINITE_DIFFERENCE')
    

  # *****************************
  #       Time Stepping
  # *****************************

  if (WITH_CHOMBO):     # ** AMR file inclusion ** 

    if (dimensions == 1): additional_files.append('PatchCTU.1D.o')
    else:                 
      if (time_stepping == 'RK_MIDPOINT'): 
        additional_files.append('PatchRK.3D.o')
      elif (time_stepping == 'EULER'):
        additional_files.append('PatchEuler.3D.o')
      elif (time_stepping == 'RK2'):
        additional_files.append('PatchEuler.3D.o')
      else:
        additional_files.append('PatchCTU.3D.o')

  else:
    if (time_stepping == 'CHARACTERISTIC_TRACING'):
       if (dim_splitting == 'YES'): additional_files.append('sweep.o')
       else:                        additional_files.append('unsplit_ctu.o')
    elif (time_stepping == 'HANCOCK'):
       if (dim_splitting == 'YES'): additional_files.append('sweep.o')
       else:                        additional_files.append('unsplit_ctu.o')
    elif (time_stepping == 'RK_MIDPOINT'):
       additional_files.append('rk_midpoint.o')
    else:
       if (dim_splitting == 'YES'): additional_files.append('sweep.o')
       else:                        additional_files.append('unsplit.o')


  if (time_stepping == 'CHARACTERISTIC_TRACING'):
    additional_files.append('char_tracing.o')
  elif (time_stepping == 'HANCOCK'):
    additional_files.append('hancock.o')


#  if (interpolation != 'LINEAR'):
#    if (time_stepping == 'CHARACTERISTIC_TRACING'):
#      additional_files.append('states_ppm_char.o')
#    else:
#      additional_files.append('rk_states.o')


  # *****************************
  #      FARGO Scheme 
  # *****************************
 
  if (WITH_FARGO):  
    pluto_path.append('Fargo/')
    additional_flags.append(' -DFARGO')

  # *****************************
  #      Geometry 
  # *****************************

 #  additional_files.append('rhs.o')
 #   additional_files.append('rhs.o')

#  if (geometry == "CARTESIAN"):
#    additional_files.append('rhs_cart.o')
#  if (geometry == "AXISYM"):
#    additional_files.append('cylsource.o')

  # *****************************
  #         Cooling 
  # *****************************

  n = entries.index('COOLING')
    
  if (default[n] == 'POWER_LAW'):
    pluto_path.append('Cooling/Power_Law/')
  elif (default[n] == 'TABULATED'):
    pluto_path.append('Cooling/Tab/')
    additional_files.append('cooling_source.o')
    additional_files.append('cooling_ode_solver.o')
  elif (default[n] == 'SNEq'):
    pluto_path.append('Cooling/SNEq/')
    additional_files.append('cooling_source.o')
    additional_files.append('cooling_ode_solver.o')
  elif (default[n] == 'MINEq'):
    pluto_path.append('Cooling/MINEq/')
    additional_files.append('cooling_source.o')
    additional_files.append('cooling_ode_solver.o')
  elif (default[n] == 'H2_COOL'):
    pluto_path.append('Cooling/H2_COOL/')
    additional_files.append('cooling_source.o')
    additional_files.append('cooling_ode_solver.o')
    

  # **********************************
  #   parabolic flux: file inclusion
  # **********************************

  if (sts_flag == 1): 
    additional_files.append('sts.o')
    additional_files.append('parabolic_rhs.o')

  if (rkc_flag == 1): 
    additional_files.append('rkc.o')
    additional_files.append('parabolic_rhs.o')

  if (exp_flag == 1): 
    additional_files.append('parabolic_flux.o')

  # **********************************************
  #          Vector potential
  # **********************************************

  if (physics_module == 'MHD' or physics_module == 'RMHD'):
    additional_files.append('vec_pot_diff.o')
    if (WITH_CHOMBO == 0): additional_files.append('vec_pot_update.o')
  
  # *****************************
  #          Particles
  # *****************************

#  n = entries.index('INCLUDE_PARTICLES')
#  if (default[n] == 'YES'):
#    pluto_path.append('Particles/')

# ****************************************************************
#       DEFINE ALL NON-USER FRIENDLY CONSTANTS
# ****************************************************************

  tmp.append('\n/* -- supplementary constants (user editable) -- */ \n\n')

  no_us_fr = []

  no_us_fr.append('#define  INITIAL_SMOOTHING     NO\n')
  no_us_fr.append('#define  WARNING_MESSAGES      NO\n')
  no_us_fr.append('#define  PRINT_TO_FILE         NO\n')
  no_us_fr.append('#define  INTERNAL_BOUNDARY     NO\n')
  no_us_fr.append('#define  SHOCK_FLATTENING      NO\n')
  if (not WITH_FD): 
    no_us_fr.append('#define  ARTIFICIAL_VISCOSITY  NO\n')
    no_us_fr.append('#define  CHAR_LIMITING         NO\n')
    no_us_fr.append('#define  LIMITER               DEFAULT\n')

# add geometry-dependent switches

  if (geometry == "AXISYM"):
    no_us_fr.append('#define  ADD_CYLSOURCE            1\n')

#  if (geometry == "POLAR" and physics_module == "HD"):
#    no_us_fr.append('#define  FARGO_SCHEME          NO\n')

# add flux ct switches for MHD or RMHD
  
  if (divb == "CONSTRAINED_TRANSPORT"):
    no_us_fr.append('#define  CT_EMF_AVERAGE           UCT_HLL\n')
    no_us_fr.append('#define  CT_EN_CORRECTION         NO\n')
    no_us_fr.append('#define  ASSIGN_VECTOR_POTENTIAL  YES\n')  
  elif (physics_module == "MHD" or physics_module == "RMHD"):
    no_us_fr.append('#define  ASSIGN_VECTOR_POTENTIAL  NO\n')  
    
     
  if (physics_module == "MHD" or physics_module == "RMHD"):
    if (not WITH_CHOMBO): 
      no_us_fr.append('#define  UPDATE_VECTOR_POTENTIAL  NO\n')  

  if (time_stepping == 'HANCOCK'):
    if (physics_module == 'RMHD'):
      no_us_fr.append('#define  PRIMITIVE_HANCOCK     NO\n')
    else:
      no_us_fr.append('#define  PRIMITIVE_HANCOCK     YES\n')

# add definition of STS_NU

  if (sts_flag == 1):
    no_us_fr.append('#define  STS_nu                0.01\n')
    
# ***************************************************************
#    Read pre-existing non-user-editable defaults;
#    Use the default value if they exist
# ***************************************************************

  if (os.path.exists(work_dir+'/definitions.h')):
    for x in no_us_fr:
      try:
        xl   = string.split(x)
        scrh = file.string_list(work_dir+'/definitions.h',xl[1])
        no_us_fr[no_us_fr.index(x)] = scrh[0]
      except IndexError:
        continue


# add elements of no_us_fr to tmp 

  for x in no_us_fr:
    tmp.append(x)

# -- create definitions.h --
    
  file.create_file(work_dir+'/definitions.h',tmp)

  return
###############################################################
#
#                  BUILD MAKEFILE
#
#
# - if no makefile is present, build one from scratch using
#   the templates in Src/Templates
#
# - if AUTO_UPDATE = 0 build a new one anyway, overwriting the 
#   old one (if present)
#
# - if AUTO_UPDATE = 1 attempt to automatically update the
#   existing makefile
#
###############################################################
def makefile(work_dir, pluto_path, pluto_dir, additional_files, 
             additional_flags, makefile_template, AUTO_UPDATE):

 mkfl_name  = work_dir+'/makefile'
 mkfl_exist = os.path.exists(mkfl_name)

# is Chombo required ?

 WITH_CHOMBO = 0
 for x in sys.argv: 
   if (x == "--with-chombo" or x == "--with-chombo:"):      
     WITH_CHOMBO = 1
   
# try to get "ARCH" from makefile

 if (mkfl_exist): 
   scrh = file.string_list(mkfl_name,'ARCH')
     
# Create a new makefile under the following conditions:

 if (AUTO_UPDATE == 0) or (not mkfl_exist) or (len(scrh) == 0):

 # define architecture

   entries = os.listdir(pluto_dir + '/Config')

   entries.sort()
   menu.SetTitle("Change makefile")
   arch        = menu.Browse(entries)  # call menu to select config. file
   arch_string = 'ARCH         = '+ arch + '\n'
 else:                           # try to update makefile automatically
   arch_string = scrh[0]
   arch        = (string.split(arch_string))[2]

# Starting print here

 row = 1
 menu.Print ("> Generating makefile... ["+arch+"]", row=row,sleep=1)
 
# if CHOMBO is required, change dir to Lib/Chombo/lib, make vars to generate 
# a list of all the variables needed by CHOMBO makefile and copy the list to 
# make.vars in local working directory.
# This file will be used later by PLUTO makefile.

 if (WITH_CHOMBO):      

 # get the number of DIMENSIONS from definitions.h
   
   scrh = file.string_list (work_dir+"/definitions.h", "DIMENSIONS")
   scrh = string.split(scrh[0])

 # build chombo configuration string

   chombo_config_string = 'DIM='+scrh[2]
   for x in sys.argv: 
     if (x == '--with-chombo:'):
       i = sys.argv.index(x) + 1
       for y in sys.argv[i:]: chombo_config_string += ' '+y
     
   row = row + 1
   menu.Print("  - Chombo config string: "+chombo_config_string,row=row) 
   row = row+1
   menu.Print("  - creating make.vars...",row=row,sleep=0) 
   os.chdir(pluto_dir+"/Lib/Chombo-3.1/lib")
   os.system("make "+chombo_config_string+" vars > make.vars\n")
   os.system("cp make.vars "+work_dir+"\n")
   os.chdir(work_dir)

# copy template
 
 shutil.copy(pluto_dir + makefile_template[0], mkfl_name)

# write main PLUTO dir

 scrh = file.word_find(mkfl_name,'PLUTO_DIR')
 ipos = scrh[0]
 file.replace_line (mkfl_name, 'PLUTO_DIR    = '+pluto_dir+'\n', ipos)

#  write architecture 

 scrh = file.word_find(mkfl_name,'ARCH')
 ipos = scrh[0]
 file.replace_line (mkfl_name, arch_string, ipos)

#  Write additional objects to makefile

 scrh = file.word_find(mkfl_name,'Additional_object_files_here')
 ipos = scrh[0] + 3

 for x in additional_files:
   file.insert(mkfl_name, 'OBJ += '+x + '\n', ipos)
   ipos = ipos + 1

# add included makefile; useful for header dependences

 for x in pluto_path:
   file.insert(mkfl_name, 'include $(SRC)/' + x + 'makefile' + '\n',ipos)
   ipos = ipos + 1

#  Write additional flags for C compiler

 for x in additional_flags:
   file.insert(mkfl_name, 'CFLAGS += '+x+'\n', ipos)
   ipos = ipos + 1

# Check for libpng

# if (WITH_CHOMBO == 0):
#   if (os.path.exists("sysconf.out")):
#     scrh = file.string_find ("sysconf.out", "LIBPNG")
#     line = file.read_lines("sysconf.out", scrh[0], scrh[0]+1)
#     line_list = string.split(line[0])
#     if (line_list[2] == 'YES'):
#       file.insert(this_makefile, 'LDFLAGS += -lpng\n', ipos)
#       ipos = ipos + 1
#       file.insert(this_makefile, 'CFLAGS += -DHAVE_LIBPNG\n', ipos)
  
 return

#####################################################################
#
#  Automatically replace old statements with new ones
#  in pre-existing definitions.h. 
#
#####################################################################
def replace_obsoletes(work_dir):

  # replace potential occurences of "INCLUDE_ROTATION" --> "ROTATING_FRAME"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','INCLUDE_ROTATION')
    tmp  = string.split(scrh[0])
    newline = "#define    ROTATING_FRAME        "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','INCLUDE_ROTATION')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  # replace potential occurences of "INCLUDE_BODY_FORCE" --> "BODY_FORCE"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','INCLUDE_BODY_FORCE')
    tmp  = string.split(scrh[0])
    newline = "#define    BODY_FORCE        "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','INCLUDE_BODY_FORCE')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  # replace potential occurences of "INCLUDE_COOLING" --> "COOLING"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','INCLUDE_COOLING')
    tmp  = string.split(scrh[0])
    newline = "#define    COOLING        "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','INCLUDE_COOLING')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  # replace potential occurences of "INCLUDE_BACKGROUND_FIELD" --> "BACKGROUND_FIELD"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','INCLUDE_BACKGROUND_FIELD')
    tmp  = string.split(scrh[0])
    newline = "#define    BACKGROUND_FIELD        "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','INCLUDE_BACKGROUND_FIELD')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  # replace potential occurences of "TIME_EVOLUTION" --> "TIME_STEPPING"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','TIME_EVOLUTION')
    tmp  = string.split(scrh[0])
    newline = "#define    TIME_STEPPING    "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','TIME_EVOLUTION')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  # replace potential occurences of "FLUX_CT" --> "CONSTRAINED_TRANSPORT"
  try:
    scrh = file.string_find(work_dir+'/definitions.h','FLUX_CT')
    file.replace_line(work_dir+'/definitions.h', 
                      "#define    MHD_FORMULATION    CONSTRAINED_TRANSPORT\n",scrh[0])
  except:
    pass


  # replace potential occurences of "RAYMOND" --> "SNEq"
  try:
    scrh = file.string_find(work_dir+'/definitions.h','RAYMOND')
    file.replace_line(work_dir+'/definitions.h', 
                      "#define    COOLING   SNEq\n",scrh[0])
  except:
    pass

  # replace potential occurences of "NEQ" --> "MINEq"
  try:
    scrh = file.string_find(work_dir+'/definitions.h','NEQ')
    file.replace_line(work_dir+'/definitions.h', 
                      "#define    COOLING   MINEq\n",scrh[0])
  except:
    pass


  # replace potential occurences of "CT_VEC_POT_INIT" --> "USE_VECTOR_POTENTIAL"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','CT_VEC_POT_INIT')
    tmp  = string.split(scrh[0])
    newline = "#define  USE_VECTOR_POTENTIAL  "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','CT_VEC_POT_INIT')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  # replace potential occurences of "USE_VECTOR_POTENTIAL" -->
  #                                 "ASSIGN_VECTOR_POTENTIAL"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','USE_VECTOR_POTENTIAL')
    tmp  = string.split(scrh[0])

    newline = "#define  ASSIGN_VECTOR_POTENTIAL  "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','USE_VECTOR_POTENTIAL')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  # replace potential occurences of "SAVE_VEC_POT" -->
  #                                 "UPDATE_VECTOR_POTENTIAL"
  try:
    scrh = file.string_list(work_dir+'/definitions.h','SAVE_VEC_POT')
    tmp  = string.split(scrh[0])
    newline = "#define  UPDATE_VECTOR_POTENTIAL  "+tmp[2]+"\n"
    n    = file.string_find(work_dir+'/definitions.h','SAVE_VEC_POT')
    n    = n[0]
    file.replace_line(work_dir+'/definitions.h', newline, n)
  except:
    pass

  return
