#
# Build ARC Module List .RST page
#
# Author: Andy Gittings
# 
#
import sys
from spider import Spider

# Read the module list 
modlist = Spider()

# Open the output file
with open('modules.rst','w') as op:

# Make file the default standard output
    so = sys.stdout 
    sys.stdout = op
    
# Iterate though all the module names
    for modName in modlist.get_names():
       print (modName)
       print ("-"*(len(modName)),"\n"*3)
       print ("**Description**","\n"*2)
       print (modlist.get_key(modName,"Description")[0].replace("\n",""),"\n"*2)
       print ("**More Information**","\n"*2)
       print (modlist.get_key(modName,"URL")[0],"\n"*2)
       print ("**Available Versions**::","\n"*2)
       for version in modlist.get_all_versions(modName):
           print("   ",version)
       print ("\n"*2)

# Restore standard output
    sys.stdout = so
