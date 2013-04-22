import sys
import string
import re
import menu

###############################################################
#
#              B U I L D       P R O B L E M
#
###############################################################

# Create user-editable menu

entries = []
options = []
default = []
scrh    = []

entries.append('PHYSICS')
options.append(['HD','RHD','MHD','RMHD'])
default.append('HD')

entries.append('DIMENSIONS')
options.append(['1','2','3'])
default.append('1')

entries.append('COMPONENTS')
options.append(['1','2','3'])
default.append('1')

entries.append('GEOMETRY')
options.append(['CARTESIAN','CYLINDRICAL','POLAR','SPHERICAL'])
default.append('CARTESIAN')

entries.append('INCLUDE_BODY_FORCE')
options.append(['NO','EXPLICIT'])
default.append('NO')

entries.append('INCLUDE_COOLING')
options.append(['NO','POWER_LAW','TABULATED','SNEq','MINEq']) # ,'H2_COOL'])
default.append('NO')

entries.append('INCLUDE_PARTICLES')
options.append(['NO','YES'])
default.append('NO')

entries.append('INTERPOLATION')
options.append(['FLAT','LINEAR','LimO3','WENO3','PARABOLIC'])  # Default
default.append('LINEAR')

entries.append('TIME_STEPPING')
options.append(['EULER','RK2','RK3','HANCOCK','CHARACTERISTIC_TRACING'])
default.append('RK2')

entries.append('DIMENSIONAL_SPLITTING')
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

menu.SetTitle('Test Menu')
menu.browse(entries, default, options)

