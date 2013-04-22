import os
import shutil
import string
import sys
import time
import make
import menu

#######################################################
def main_menu(pluto_dir):
#
#
#
#
#
#
#######################################################

# check if setup has been invoked with --with-chombo
    
  work_dir = os.getcwd()
  options = ['Setup problem', 'Change makefile', 
             'Auto-update','Save Setup','Quit']

  additional_files  = ['']
  additional_flags  = ['']
  pluto_path        = ['']
  makefile_template = []
  what              = ''
  
  while (what != options[4]):
  
    menu.SetTitle("Python setup (09/2012)", 
                  "Working dir: "+work_dir+"\nPLUTO dir  : "+pluto_dir)
    what = menu.Browse(options)
     
    if   (what == options[0]):   # build problem 

      if (not os.path.exists(work_dir+'/init.c')):
        shutil.copy(pluto_dir+'/Src/Templates/init.c',work_dir+'/init.c')
        
      if (not os.path.exists(work_dir+'/pluto.ini')):
        shutil.copy(pluto_dir+'/Src/Templates/pluto.ini',work_dir+'/pluto.ini')

      make.problem(work_dir, pluto_path, pluto_dir, additional_files, 
                   additional_flags, makefile_template, 0)
      make.makefile(work_dir, pluto_path, pluto_dir, additional_files, 
                    additional_flags, makefile_template, 1)

    elif (what == options[1]):   # makefile
    
      # --------------------------------------------------------
      #   Read problem definitions again;
      #   Useful when one wants to change the makefile 
      #   without necessarily going to build/change problem.
      # --------------------------------------------------------
      
      make.problem(work_dir, pluto_path, pluto_dir, additional_files,  
                   additional_flags, makefile_template, 1)
      make.makefile(work_dir, pluto_path, pluto_dir, additional_files, 
                    additional_flags, makefile_template, 0)

    elif (what == options[2]):    # automatic update

      menu.Prompt('Press Enter to Update '+work_dir)

      make.problem (work_dir, pluto_path, pluto_dir, additional_files, 
                    additional_flags, makefile_template, 1)

      make.makefile(work_dir, pluto_path, pluto_dir, additional_files, 
                    additional_flags, makefile_template, 1)
      menu.Print ("Configuration up to date",sleep=0.75)
      sys.exit(1)             
    elif (what == options[3]):   # save setup
  
      # ------------------------------------
      #  get the path basename and make 
      #  a tar archive 
      # ------------------------------------

      menu.RestoreScreen()
      os.chdir(work_dir)
      default_name = os.path.basename(work_dir)
      os.system("clear")
      setup_name = raw_input(" > Archive name ("+default_name+"): ")
      complete_backup = 0
      scrh       = raw_input(" > Backup source files ? (n): ")
      if (setup_name == ''): 
        setup_name = default_name
 
      if (scrh == "y"): 
        complete_backup = 1

      setup_name = setup_name+".tar"
      ifail = os.system ("tar cvf "+setup_name+" *.c *.ini *.h")

      # -------------------------------
      #  back up source tree
      # -------------------------------

      if (complete_backup):
        ifail = ifail or os.system("tar rvf "+setup_name+" --directory="+pluto_dir+"/   Src/ ")
        ifail = ifail or os.system("tar rvf "+setup_name+" --directory="+pluto_dir+"/   Config/")
        ifail = ifail or os.system("tar rvf "+setup_name+" makefile")
          
      ifail = ifail or os.system ("gzip -f "+setup_name)
      if (ifail == 0):
        print " > "+setup_name+".gz successfully created\n"
      else:
        print " ! Error creating "+setup_name+".gz\n"
 
      sys.exit()
      

  menu.RestoreScreen()
  print "\n\n Bye."
  sys.exit()
    


