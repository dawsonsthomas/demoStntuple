#!/usr/bin/env python

import os, re, string, subprocess

Import('env')
from stntuple_helper import *
#------------------------------------------------------------------------------
# last two components of the path. Ex: /not/this/but/THIS/AND_THIS
#                                      "AND_THIS" is usually "src"
#------------------------------------------------------------------------------
x = subprocess.call('scripts/build_config',shell=True)
#------------------------------------------------------------------------------
# done
#------------------------------------------------------------------------------
env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen})
env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

Export('stntuple_helper')
# Import('env')
# # Import('mu2e_helper')
# #------------------------------------------------------------------------------
# # last two components of the path. Ex: /not/this/but/THIS/AND_THIS
# #                                      "AND_THIS" is usually "src"
# #------------------------------------------------------------------------------
# # x = subprocess.call('scripts/build_config',shell=True)
# if (os.environ.get("MU2E_SATELLITE_RELEASE")) :
#     env['CPPPATH' ].append('-I'+os.environ['MU2E_SATELLITE_RELEASE']+'/include');
#     env['CXXFLAGS'].append('-I'+os.environ['MU2E_SATELLITE_RELEASE']+'/include');
# else :
#     env['CPPPATH' ].append('-I'+os.environ['MU2E_BASE_RELEASE']+'/include');
#     env['CXXFLAGS'].append('-I'+os.environ['MU2E_BASE_RELEASE']+'/include');
    
# #------------------------------------------------------------------------------
# # done
# #------------------------------------------------------------------------------
# from gianipez_helper import *
# #from site_init import *

# env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen })
# env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

# Export('gianipez_helper')
