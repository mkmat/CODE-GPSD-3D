# CODE-GPSD-3D
Generalized geometric pore size distribution (G-PSD) for periodic systems composed of spheres. This quantity was defined by Gelb and Gubbins, there are other pore size distributions such as T-PSD, as discussed in detail in the accompanying [article](#about)). The spheres constitute the 'material', which is surrounded by 'pore space'. This codes allows to calculate the distribution *P(r;r<sub>p</sub>|r<sub>c</sub>)*  of pore radii *r* for such system as function of the radius *r<sub>p</sub>* of a thought probe particle and a thought coating thickness *r<sub>c</sub>* of the material spheres. 

<img src="images/schematic-GPSD3D-github.png" width="100%">  

For monodisperse systems GPSD-3D uses the advanced (grid-free) voronoi-based algorithm. For polydisperse systems, it uses a basic grid-based algorithm, whose resolution is limited by the amount of available memory. The settings of the grid-based algorithm are hard-coded in the file GPSD-3D, and can be modified there (see [below](#hardcoded) for details).  

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

1. x y z (monodisperse system, requires specifying -ro on the command line)
2. id x y z (monodisperse system, requires specifying -ro on the command line)        
3. id x y z radius  (polydisperse system, monodisperse if all radii are equal)
4. samarth-type configuration file (do not specify box dimensions in that case, requires specifying -ro on the command line)

The six values for the box can be either saved in a txt-file (single line, six values xlo xhi ylo yhi zlo zhi separated by blank or commata), or passed over on the command line.  


## How to run the code. Parameters. Command-line options

      perl GPSD-3D
          -in=<filename>
          -box=<filename> OR -xlo=.. -xhi=.. ylo=.. -yhi=.. -zlo=.. -zhi=..
          [-ro=<value>] 
          [-rp=<value>] 
          [-rc=<value>] 
          [-q=<integer>] 
          [-o=<filename>]
          [-np=<integer>]
          [-more]
          [-info]
          [-quiet]
          [-clean]

**-in=configfilename**:    name of the file containing the configuration (as described above)

**-box=< boxfile >**:      name of the file containing the box geometry (as described above, alternatively, the box size can be passed over on the command line)

**-ro= *r<sub>o</sub>***:   particle radius *r<sub>o</sub>* (required, if the particle radii are not contained in the input file. If *r<sub>o</sub>* is specified, existing radii in the input file are ignored)

**-rp= *r<sub>p</sub>***:   (optionally) probe particle radius *r<sub>p</sub>*  (if not specified, *r<sub>p</sub>=0* is used)

**-rc= *r<sub>c</sub>***:   (optionally) shell thickness *r<sub>c</sub>* (if not specified, *r<sub>c</sub>=0* is used)

**-q= *q***:    (optionally) positive quality value *q* (if not specified, *q=1* is used. The number of random shots is 10000 times *q*)

**-o=outputfilename**:     (optionally) name of the resulting file containing a list of pore radii (if not specified, the list is saved in configfilename.gpsd)

**-np= *n<sub>p</sub>***:  (optionally) number of threads *n<sub>p</sub>* to be used (default: maximum number of threads)

**-more**:                 (optionally) add, besides pore radius (column 1), the corresponding pore center (columns 2,3,4) to the outputfile (for visualization purposes)

**-info**:                 (optionally) runtime information (cpu times etc) is collected in a file whose name ends with -info (see <a href="#info">below</a> for detailed information)

**-quiet**:                (optionally) do not produce any stdout.

**-clean**:                GPSD-3D -clean removes all temporary directories that may have been generated during previous crashes. 

### Comments: 

1. The leading 'perl' is only needed, if your system does not automatically recognize GPSD-3D to be a perl file. If GPSD-3D is not found, call it via: perl ./GPSD-3D.
2. If called without any argument, GPSD-3D displays the documentation.
3. GPSD-3D can be called in parallel. 
4. Each GPSD-3D call runs in a unique temporary directory, upon successful completion the temporary directory is removed.
5. The maximum number of threads used by OpenMP is reported during the installation and also if GPSD-3D is called without arguments. For very large systems with, say, more than 1000000 spheres, running at the maximum number of threats must not be an advantage and you should check the speed also using a single processor, using -np=1


## Output

GPSD-3D returns a list of pore radii *r* in a file, either in configfilename.gpsd or outputfilename, if the latter had explicitly been defined on the command line. 

        0.658312
        0.274754
        1.070546
        0.685289
        ...

This list of *r* values (for the chosen values *r<sub>p</sub>* and *r<sub>c</sub>*) gives rise to a distribution of pore radii, so called generalized geometric pore radius distribution *P(r;r<sub>p</sub>|r<sub>c</sub>)*. The bare radius of the particles *r<sub>o</sub>* is usually not mentioned here, as it belongs to the system. For monodisperse systems only the sum or *r<sub>o</sub>+r<sub>c</sub>* matters. For polydisperse systems each spherical particle has its own radius according to the configuration file, and *r<sub>c</sub>* can be used to effectively modify the stored particle radii, without changing the configuration file. If you call GPSD-3D with the -more option, the same file will contain four columns (no header)

        r        x    y    z
        0.658312 1.31 2.13 1.99
        0.274754 2.30 1.01 4.02
        1.070546 ...
        0.685289 ...
        ...

where *x*,*y*,*z* are the center coordinates of the pore with radius *r*.

## Test configurations and test runs

A number of configurations and corresponding box-files are available from the current respository. They are named .benchmark-x-config and .benchmark-x-box. A test call, using 10 of the available threads, 50000 Monte Carlo trials (*q=5*), for a probe sphere with zero radius, and materials spheres of radius 1.0 is 

        perl ./GPSD-3D -in=.benchmark-7-config -box=.benchmark-7-box -rp=0.0 -ro=1.0 -q=5 -np=10

As we did not suppress stdout via -quiet, it should produce the following:

       _______________________________________________________________________________________________________________________________

        This is GPSD-3D version 1.0 written by Martin Kroger 2023, ETH Zurich, https://www.complexfluids.ethz.ch, Email: mk@mat.ethz.ch

        Related publication (GPSD-3D): Comput. Phys. Commun. (2023) submitted
        Related publication (GPSD-2D): Phys. Rev. E 107 (2023) 015307. Link: http://doi.org/DOI:10.1103/PhysRevE.107.015307

        GPSD-3D Code available from: https://github.com/mkmat/CODE-GPSD-3D
        GPSD-2D Code available from: https://github.com/mkmat/CODE-GPSD-2D
        _______________________________________________________________________________________________________________________________


        [INFO] using configuration file .benchmark-7-config
        [INFO] using box file .benchmark-7-box
        [PREPARING] scanning .benchmark-7-box
        [PREPARING] recognized format (B)
        [INFO] monodisperse: 1
        [INFO] .benchmark-7-config contains 2000 particle coordinates (4 columns)
        [INFO] created files in .tmp-GPSD-3D-70888 including .parameters.
        [INFO] monodisperse system. The particle radius is taken as 1, shell thickness 0, test particle radius 0.
        [GPSD-3D] Using 49996 shots on 10 threads
        [GPSD-3D] Please stand by ..
        [GPSD-3D]                               reading box ..
        [GPSD-3D]                                          box      24.000       24.000       24.000
        [VORO++]                                  max_vertices          53
        [VORO++]                                voro_max_faces          31
        [VORO++]                    voro_max_vertices_for_face          12
        [GPSD-3D]                                    triangles      119244
        [GPSD-3D]                      parallel processes (np)          10
        [GPSD-3D]                                            N        2000
        [GPSD-3D]                                    triangles      119244
        [GPSD-3D]                              UpperPoreRadius       2.847
        [GPSD-3D]                                           ro       1.000
        [GPSD-3D]                                           rc       0.000
        [GPSD-3D]                                           rp       0.000
        [GPSD-3D]                                 reff = rc+rp       0.000
        [GPSD-3D]                    creating neighbor list ..
        [GPSD-3D]                               neighborlist_M           4            4            4
        [GPSD-3D]                            neighborlist_size       6.000        6.000        6.000
        [GPSD-3D]                                  start MC ..
        [GPSD-3D]                    cpu per 1000 shots [secs]       1.416
        [GPSD-3D]                    volume fraction phi(reff)       0.474
        [GPSD-3D]                                      V(0|rc)    7266.752
        [GPSD-3D]                              min pore radius       0.075
        [GPSD-3D]                             mean pore radius       1.793 +/-        0.002
        [GPSD-3D]                              max pore radius       3.418
        [GPSD-3D]             created a list {r} of pore radii
        [GPSD-3D]                    shots (use -q to enlarge)       49996
        [GPSD-3D]    cpu time spent in read_voro_output [secs]       0.159
        [GPSD-3D]     cpu time spent in setup_triangles [secs]       0.002
        [GPSD-3D]          cpu time spent in MonteCarlo [secs]      70.788
        [GPSD-3D]         real time spent in MonteCarlo [secs]       7.320
        [GPSD-3D] completed after 8.35455751419067 secs



and the following file (a list of roughly 14000 *r* values) should have been generated (if you do not see it, type: ls -lat): 

        .benchmark-7-config-ro=1-rp=0-rc=0.gpsd

With such list of values at hand, creating the normalized histogram (the pore radius distribution) is straightforward using any software that can bin the values, and visualize a graph. Some quantities derived from the list, such as minimum and maximum pore radius, as well as the mean pore radius including its standard error are mentioned in the above stdout. If you are not satisfied with the name of the resulting file, use the -o option. 

## DEBUGGING - REFERENCE CODE

    [GPSD-3D]                          reading box ..
    [GPSD-3D]                                     box    24.000     24.000     24.000
    [GPSD-3D]                  reading _triangles.txt
    [GPSD-3D]     reading and processing triangles ..
    [GPSD-3D]                                       N      2000
    [GPSD-3D]                               triangles    119244
    [GPSD-3D]                         UpperPoreRadius     2.847
    [GPSD-3D]                                   ro+rc     1.000
    [GPSD-3D]                                      rp     0.000
    [GPSD-3D]                            reff = rc+rp     0.000
    [GPSD-3D]              max_triangle_max_extension     3.461
    [GPSD-3D]               creating neighbor list ..
    [GPSD-3D]                          neighborlist_M         3          3          3
    [GPSD-3D]                       neighborlist_size     8.000      8.000      8.000
    [GPSD-3D]                             start MC .. 
    [GPSD-3D]               cpu per 1000 shots [secs]     2.516
    [GPSD-3D]               volume fraction phi(reff)     0.310
    [GPSD-3D]                                 V(0|rc)  9541.325
    [GPSD-3D]                         min pore radius     0.019
    [GPSD-3D]                        mean pore radius     1.542 +/- 0.003
    [GPSD-3D]                         max pore radius     2.847
    [GPSD-3D]        created a list {r} of pore radii
    [GPSD-3D]               shots (use -q to enlarge)     20000
    [GPSD-3D]   time spent in read_voro_output [secs]     0.241
    [GPSD-3D]    time spent in setup_triangles [secs]     0.002
    [GPSD-3D]         time spent in MonteCarlo [secs]    50.320
    [GPSD-3D] completed
    [GPSD-3D] created: .benchmark-7-config-ro=1-rp=0-rc=0.gpsd

This is the voro++ output

1 10.5 17.5 0.5 11 (-1,-0,0) (-0.447214,0,-0.894427) (-0.5547,-0.83205,2.51468e-16) (0,-0.948683,0.316228) (0.447214,-0.894427,0) (0.707107,-0.707107,0) (5.12411e-17,2.13504e-17,1) (9.5162e-17,1,1.23711e-16) (1,0,0) (0,0,-1) (-0.801784,0.534522,0.267261) 16 (-0.785714,-1.64286,0.0714286) (-2.5,-0.5,-1.38778e-16) (-0.5,-1.5,0.5) (0.5,0.5,-0.5) (-0.5,-1.5,-0.5) (0.5,0.5,0.5) (-1.5,0.5,-0.5) (-1.5,-1.16667,-0.5) (0.5,-0.5,0.5) (0.5,-0.5,-0.5) (-1.83333,0.5,0.5) (-2.07143,0.5,-0.214286) (-0.785714,-1.64286,-0.5) (-2.5,-0.25,-1.59595e-16) (-1,-1.5,0.5) (-2.5,-0.5,0.5) (1,13,15) (1,7,6,11,13) (1,15,14,0,12,7) (2,0,14) (2,4,12,0) (2,8,9,4) (2,14,15,10,5,8) (3,5,10,11,6) (3,9,8,5) (3,6,7,12,4,9) (10,15,13,11)
2 12.5 23.5 10.5 15 (-0.57735,-0.57735,-0.57735) (-1.11022e-16,-6.16298e-33,-1) (0.57735,-0.57735,-0.57735) (0.707107,1.10397e-16,-0.707107) (0.894427,-0.447214,6.4714e-16) (-1,1.11022e-16,-3.33067e-16) (-0.408248,-0.816497,-0.408248) (-0.707107,9.15887e-17,0.707107) (0.816497,0.408248,0.408248) (0.57735,0.57735,-0.57735) (-7.06542e-16,0.707107,-0.707107) (-9.93014e-17,0.894427,-0.447214) (-0.534522,-0.801784,0.267261) (-0.57735,-0.57735,0.57735) (0.707107,-1.96262e-17,0.707107) 21 (0.5,-0.5,-0.5) (5.55112e-17,-1,-0.5) (1,-0.5,1.22125e-15) (-0.5,-1.5,0.5) (-0.5,-0.5,-0.5) (0.5,1.5,0.5) (0.25,-2,0.75) (0.2,-2.1,1.1) (-0.5,0.5,-0.5) (0.0714286,1.78571,1.07143) (0.5,0.5,-0.5) (-0.5,1.5,0.5) (3.70074e-17,-1.5,2.77556e-16) (1.3,0.1,0.3) (1.16667,-0.166667,0.833333) (1.16667,0.5,0.166667) (0.5,0.5,1.5) (0.5,-1.5,1.5) (0.192308,-2.11538,1.03846) (0.25,-2,1.25) (-2.06779e-15,-2,1) (1,4,3,12) (1,0,10,8,4) (1,12,6,2,0) (2,13,15,10,0) (2,6,18,7,19,17,14,13) (3,4,8,11) (3,20,18,6,12) (3,11,9,16,17,19,20) (5,15,13,14,16,9) (5,10,15) (5,11,8,10) (5,9,11) (7,18,20) (7,20,19) (14,17,16)


## Polydisperse systems: Grid-based <a name="hardcoded">

For the case of polydisperse systems, the GPSD-3D script contains two lines that may be edited by users to increase or reduce the resolution further. The default setting is:  

        $min_delta_grid     = 0.005;      # USER-defined minimum grid spacing (in units of the effective particle radius ro+rc)
        $maxvoxels_grid     = 1000000;    # USER-defined upper limit for number of voxels 

## -info file <a name="info">

    N (number of material spheres)
    ro (material sphere radius)
    rp (probe sphere radius)
    rc (shell thickness)
    box volume
    total number of triangles
    total number of shots
    total number of pore radius values
    maximum pore radius
    minimum pore radius
    mean pore radius
    maximum pore radius detected
    standard error of the mean pore radius
    volume fraction phi(reff)
    cpu time used to process the voro++ output
    cpu time used to setup triangles
    cpu time spent during Monte Carlo (affected by the number of processors -np)
    real time spent during Monte Carlo
    number of neighbor cells
    np number of threads used in parallel
    
    

## About <a name="about">

Related publication (GPSD-3D): Comput. Phys. Commun. (2023) submitted

Related publication (GPSD-2D): Phys. Rev. E 107 (2023) 015307. Link: http://doi.org/DOI:10.1103/PhysRevE.107.015307

GPSD-3D Code available from: https://github.com/mkmat/CODE-GPSD-3D

GPSD-2D Code available from: https://github.com/mkmat/CODE-GPSD-2D

