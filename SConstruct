import os
import sys
import commands

vars = Variables()
vars.AddVariables(
  EnumVariable( 'order',
                'convergence order of the ADER-DG method',
                'none',
                allowed_values=('none', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12')
              ),
  PathVariable( 'buildDir', 'where to build the code', 'build', PathVariable.PathIsDirCreate ),
  EnumVariable( 'compileMode', 'mode of the compilation', 'release',
                allowed_values=('debug', 'release')
              ),
  EnumVariable( 'compiler',
                'Select the compiler (default: intel)',
                'intel',
                allowed_values=('intel', 'gcc')),
)

# set environment
env = Environment(variables=vars)
env['ENV'] = os.environ

# generate help text
Help(vars.GenerateHelpText(env))

# handle unknown, maybe misspelled variables
unknownVariables = vars.UnknownVariables()

# exit in the case of unknown variables
if unknownVariables:
  ConfigurationError("*** The following build variables are unknown: " + str(unknownVariables.keys()))

if env['order'] == 'none':
  ConfigurationError("*** Convergence order not set.")

#
# preprocessor, compiler and linker
#

# Basic compiler setting
if env['compiler'] == 'intel':
    env['CC'] = 'icc'
    env['CXX'] = 'icpc'
    env['F90'] = 'ifort'
elif env['compiler'] == 'gcc':
    env['CC'] = 'gcc'
    env['CXX'] = 'g++'
    env['F90'] = 'gfortran'
else:
    assert(false)

#
# Common settings
#

# enforce restrictive C/C++-Code
env.Append(CFLAGS   = ['-Wall', '-Werror', '-ansi'],
           CXXFLAGS = ['-Wall', '-Werror', '-ansi'])

#
# Compile mode settings
#

# set (pre-)compiler flags for the compile modes
if env['compileMode'] == 'debug':
  env.Append(CFLAGS  = ['-O0', '-g'],
             CXXFLAGS = ['-O0', '-g'])
elif env['compileMode'] == 'release':
  env.Append(CPPDEFINES = ['NDEBUG'])
  env.Append(CFLAGS   = ['-O3'],
             CXXFLAGS = ['-O3'])

#
# Basic preprocessor defines
#

env.Append(CPPDEFINES=['CONVERGENCE_ORDER='+env['order']])
env.Append(CPPPATH=['#/src'])

#
# setup the program name and the build directory
#
env['programName'] = 'lina'
env['programFile'] = '%s/%s' %(env['buildDir'], env['programName'])

# build directory
env['buildDir'] = '%s/build_%s' %(env['buildDir'], env['programName'])

# get the source files
env.sourceFiles = []

Export('env')
SConscript('src/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
Import('env')

# build standard version
env.Program('#/'+env['programFile'], env.sourceFiles)
