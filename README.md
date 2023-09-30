# CODE-GPSD-3D
Generalized pore size distribution for periodic systems of monodisperse spheres

<img src="images/schematic-GPSD3D-github.png" width="40%">

## Installation

Download GPSD-3D.tar and unpack it in a new directory named GPSD-3D: 

    gunzip GPSD-3D.tar.gz; tar -x GPSD-3D.tar; 

This will create the following files: GPSD-3D.f90, and GPSD-3D-grid.f90, INSTALL.pl, TEMPLATE.pl, MESSAGE.txt. Then switch to the new directory and call the installation script via

    cd GPSD-3D; perl ./INSTALL.pl

This will produce some message on your screen like 

      
    This is GPSD-3D version 1.0 written by Martin Kroger 2023, ETH Zurich, https://www.complexfluids.ethz.ch, Email: mk@mat.ethz.ch

    Related publication (GPSD-3D): Comput. Phys. Commun. (2023) submitted
    Related publication (GPSD-2D): Phys. Rev. E 107 (2023) 015307. Link: http://doi.org/DOI:10.1103/PhysRevE.107.015307

    GPSD-3D Code available from: https://github.com/mkmat/CODE-GPSD-3D
    GPSD-2D Code available from: https://github.com/mkmat/CODE-GPSD-2D

and a new file GPSD-3D (a perl script). Ideally, it finishes with saying 

    You have successfully installed GPSD-3D
    Copy GPSD-3D (just this single file) to a directory where you'll need it or where it can be found.
    Then call (i) GPSD-3D or (ii) ./GPSD-3D or (iii) perl GPSD-3D to start GPSD-3D (includes documentation)

But, if you have voro++ not yet installed, or no fortran compiler, you will end up with an error message that provides you with a download link. In this case you need to do the installation, and call the installation script INSTALL.pl again. 

## Configuration file formats

The input required by GPSD-3D are (i) coordinates: the center positions of *N* monodisperse or polydisperse spheres, that are located inside a rectangular, periodic box, whose corners are specified by (ii) box geometry: 6 values: xlo,xhi,ylo,yhi,zlo,zhi. We offer various file formats, where each of the *N* lines contains the coordinates, and eventually also the radius of a single sphere

1. x y z (monodisperse system)
2. id x y z (monodisperse system)        
3. id x y z radius  (polydisperse system, monodisperse if all radii are equal)
4. samarth-type configuration file (do not specify -boxfile in this case)

The six values for the box can be either saved in a txt-file (single line, six values xlo xhi ylo yhi zlo zhi separated by blank or commata), or passed over on the command line.  


## How to run the code. Parameters. Command-line options

      perl GPSD-3D 
          -in=<filename>
          [-box=<filename>] OR [-xlo=.. -xhi=.. ylo=.. -yhi=.. -zlo=.. -zhi=..]
          [-rp=<value>] 
          [-ro=<value>] 
          [-rc=<value>] 
          [-q=<integer>] 
          [-quiet]

If called without any argument, GPSD-3D shows the documentation.

**-in=filename**:    name of the file containing the configuration (as described above)

**-box=< boxfile >**:      name of the file containing the box geometry (as described above, alternatively, the box size can be passed over on the command line)

**-ro= *r<sub>o</sub>***:   particle radius *r<sub>o</sub>* (required, if the particle radii are not contained in the input file. If ro is specified, existing radii in the input file are ignored)

**-rp= *r<sub>p</sub>***:   (optionally) probe particle radius *r<sub>p</sub>*  (if not specified, rp=0 is used)

**-rc= *r<sub>c</sub>***:   (optionally) shell thickness *r<sub>c</sub>* (if not specified, rc=0 is used)

**-q= *q***:    positive quality integer *q* (if not specified, quality=1 is used. The number of random shots is 10000 times *q*)



## Output

GPSD-3D reads the configuration and 

## How to run GPSD-3D



