# ADPcomp

Measure atomic mean-square displacement amplitudes along arbitrary directions.

## Installation


```
$ git clone https://github.com/dwpaley/ADPcomp
$ chmod +x ADPcomp/ADPcomp.py
$ ln -s $PWD/ADPcomp/ADPcomp.py /usr/local/bin/adpcomp
```


## Usage

Help and information are provided by calling: adpcomp -h

Example input file:

```
file perovtest.cif  # load a cif file

eqiv $1 3/2-x, 1/2-y, -1/2+z  # Define a symmetry code

ip br02 pb01 pb01_$1  # measure atom br02 on the vector bisecting br02-pb01 and br02-pb01_$1 
oop br02 pb01 pb01_$1  # measure br02 along the normal to the plane of the three atoms

bond br02 pb01  # measure br02 along the axis br02-pb01
axis br02 pb01 pb01_$1  # measure br02 along the axis pb01-pb01_$1

uvw br02 0 0 1  # measure br02 along direct axis 001
hkl br02 2 0 1  # measure br02 along reciprocal axis 001
```

This file would be saved as `run.aip` and called with `$ ADPcomp.py run.aip`. The output
is `run.aop` with the requested measurements appended to each line.

An interactive mode is also available with `$ adpcomp -i`.






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
