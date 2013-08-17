#****************************************************************************
# Copyright (C) 2013 SUNCAT
# This file is distributed under the terms of the
# GNU General Public License. See the file `COPYING'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#****************************************************************************

from sys import argv,executable,exit,stderr

if len(argv)<2 or len(argv)>3 or argv[1]!='install' or len(argv)==3 and argv[2][:9]!='--prefix=':
    print >>stderr, 'usage: '+executable+' '+argv[0]+' install [--prefix=installation-directory]'
    exit(1)

from os import path
from distutils.sysconfig import get_python_lib
pylib = get_python_lib()

if len(argv)==3:
    paths = [path.realpath(x) for x in [pylib,executable]]
    prefix = path.commonprefix(paths)
    instdir = argv[2][9:]
    pylib = path.join(instdir,paths[0].replace(prefix,''))
else:
    instdir = path.normpath(path.join(pylib,'../../..'))

from os import system
exit(system('make doinstall instdir='+instdir+' py='+executable+' pylib='+pylib))
