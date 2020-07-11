ADPcomp is written for Python3 and requires NumPy and SymPy. These can be
installed easily if you use an environment manager such as Anaconda, but can
also be installed by the following:

$ pip3 install numpy 
$ pip3 install sympy

To install ADPcomp, copy ADPcomp.tar.gz to a directory such as 
/applications/xray/ADPcomp and extract the files with:

$ tar -xvzf ADPcomp.tar.gz

Make the main script executable and place a link in the system path:

$ chmod +x adpcomp.py
$ ln -s $PWD/adpcomp.py /usr/local/bin

You may have to change the Python3 path in the first line of ADPcomp.py. The
correct path can be found with:

$ which python3


Further help and information are provided by calling: adpcomp -h







********************************************************************************
License and copyright information:

Copyright 2016, Daniel W. Paley
Contact: dwp2111@columbia.edu

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
