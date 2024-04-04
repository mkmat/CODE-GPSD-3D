_______________________________________________________________________________________________________________________________

This is GPSD-3D version 1.0
written Nov 2023 by Martin Kroger, ETH Zurich, www.complexfluids.ethz.ch, Email: mk@mat.ethz.ch

Related publication (GPSD-3D): Comput. Phys. Commun. (2024) in press (see Github link below).
Related publication (GPSD-2D): Phys. Rev. E 107 (2023) 015307. Link: http://doi.org/DOI:10.1103/PhysRevE.107.015307

GPSD-3D Code available from: https://github.com/mkmat/CODE-GPSD-3D
GPSD-2D Code available from: https://github.com/mkmat/CODE-GPSD-2D
_______________________________________________________________________________________________________________________________

_______________________________________________________________________________________________________________________________

GPSD-3D.tar contains the following files
_______________________________________________________________________________________________________________________________:

GPSD-3D.f90
GPSD-3D-grid.f90
GPSD-3D-nlin.f90
convert_samarth_to_config.pl
replicate.pl
voro-parser.cpp
voro-parser-0.cpp
get-max-threads.f90
MESSAGE.txt
LICENSE.txt
TEMPLATE.pl
INSTALL.pl
GPSD-3D-benchmark-tests.pl
hard-code-software-location.pl
windows.pl
software-locations-example.txt
benchmark/benchmark-10-config
benchmark/benchmark-11-config
benchmark/benchmark-12-config
benchmark/benchmark-1-config
benchmark/benchmark-2-config
benchmark/benchmark-3-config
benchmark/benchmark-4-config
benchmark/benchmark-5-config
benchmark/benchmark-6-config
benchmark/benchmark-7-config
benchmark/benchmark-8-config
benchmark/benchmark-9-config
benchmark/benchmark-10-box
benchmark/benchmark-11-box
benchmark/benchmark-12-box
benchmark/benchmark-1-box
benchmark/benchmark-2-box
benchmark/benchmark-3-box
benchmark/benchmark-4-box
benchmark/benchmark-5-box
benchmark/benchmark-6-box
benchmark/benchmark-7-box
benchmark/benchmark-8-box
benchmark/benchmark-9-box
README.txt

_______________________________________________________________________________________________________________________________

Documentation
_______________________________________________________________________________________________________________________________

Documentation is available from the related publication and the above-mentioned github repository:
https://github.com/mkmat/CODE-GPSD-3D

After installation, help is available via: perl ./GPSD-3D

_______________________________________________________________________________________________________________________________

Installation
_______________________________________________________________________________________________________________________________

Unpack the tar-ball. Then type: perl ./INSTALL.pl

GPSD-3D requires a fortran90 and c++ compiler, as well as voro++.

In case you need to install any of the three, the INSTALL.pl script will ask for. Here are the links:
g++ is available for free download at https://gcc.gnu.org
gfortran is available for free download at https://gcc.gnu.org
voro++ is available for free download at https://math.lbl.gov/voro++

The install-script should create a file software-locations.txt. 
An example for such file is available as software-locations-example.txt.

If any of the links does not work, consult the GPSD-3D github repository for updates.

_______________________________________________________________________________________________________________________________

Questions/Ideas
_______________________________________________________________________________________________________________________________

For questions or ideas please contact mk@mat.ethz.ch
