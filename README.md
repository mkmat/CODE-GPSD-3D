# CODE-GPSD-3D
Generalized geometric pore size distribution (G-PSD) for periodic systems composed of spheres. This quantity was defined by Gelb and Gubbins, there are other pore size distributions such as T-PSD, as discussed in detail in the accompanying [article](#about). The spheres constitute the 'material', which is surrounded by 'pore space'. This codes allows to calculate the distribution *P(r;r<sub>p</sub>|r<sub>c</sub>)*  of pore radii *r* for such system as function of the radius *r<sub>p</sub>* of a thought probe particle and a thought coating thickness *r<sub>c</sub>* of the material spheres. 

<img src="images/schematic-GPSD3D-github.png" width="100%">  

For monodisperse systems GPSD-3D uses the advanced (grid-free) voronoi-based algorithm. For polydisperse systems, it uses a basic grid-based algorithm, whose resolution is limited by the amount of available memory. The settings of the grid-based algorithm are hard-coded in the file GPSD-3D, and can be modified there (see [below](#hardcoded) for details).  

## Installation

Create a new directory named GPSD-3D. Download GPSD-3D.tar.gz to the new directory GPSD-3D and unpack it via 

    gunzip GPSD-3D.tar.gz; tar -x GPSD-3D.tar; 

This will create the following files: GPSD-3D.f90, and GPSD-3D-grid.f90, voro-parser.cpp, get-max-threads.f90, INSTALL.pl, TEMPLATE.pl, MESSAGE.txt. Then switch to the new directory and call the installation script via

    perl ./INSTALL.pl

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

## Configuration and box file formats <a name=format>

The input required by GPSD-3D are (i) coordinates: the center positions of *N* monodisperse or polydisperse spheres, that are located inside a rectangular, periodic box, whose corners are specified by (ii) box geometry: 6 values: xlo,xhi,ylo,yhi,zlo,zhi. We offer various file formats, where each of the *N* lines contains the coordinates, and eventually also the radius of a single sphere

1. x y z (monodisperse system, requires specifying -ro on the command line)
2. id x y z (monodisperse system, requires specifying -ro on the command line)        
3. id x y z radius  (polydisperse system, monodisperse if all radii are equal)
4. samarth-type configuration file (do not specify box dimensions in that case, requires specifying -ro on the command line)

The six values for the box can be either saved in a txt-file (single line, six values xlo xhi ylo yhi zlo zhi separated by blank or commata), or passed over on the command line. The delimiting character can be specified using the -d option.  


## How to run the code. Parameters. Command-line options

      perl GPSD-3D
          -in=<filename>
          -box=<filename> OR -xlo=.. -xhi=.. ylo=.. -yhi=.. -zlo=.. -zhi=..
          [-ro=<positive value>] 
          [-rp=<value>] 
          [-rc=<value>] 
          [-q=<integer>] 
          [-o=<filename>]
          [-np=<integer>]
          [-d=<delimiter>]
          [-more]
          [-info]
          [-quiet]
          [-clean]

**-in=**    name of the file containing the material sphere [coordinates](#format) (for polydisperse systems also the sphere radii).

**-box=**   name of the file containing the box [geometry](#format) (alternatively, the box size can be passed over on the command line using -xlo= .. -zhi=..).

**-ro=** material circle radius *r<sub>o</sub>* (required, if the particle radii are not contained in the input file). If *r<sub>o</sub>* is specified, existing radii in the input file are ignored).

**-rp=** probe particle radius *r<sub>p</sub>*. If not specified, *r<sub>p</sub>=0* is used.

**-rc=** shell thickness  *r<sub>c</sub>*. If not specified, *r<sub>c</sub>=0* is used.

**-q=** positive quality value (optionally). If not specified, *q=10.0* is used. The number of random shots is *q* times the number of material spheres.

**-o=** name of the resulting file containing a list of pore radii. If not specified, the list is saved in a .gpsd-file.

**-np=** number of threads np to be used. If not specified, *n<sub>p</sub>* is set to the available number of threads.

**-d=** delimiter (single character) present in the configuration and box-files, e.g. -d="\ " for a blank, default is -d=","

**-more**: tell the code to add, besides pore radius *r*, the corresponding probe particle center **p** and pore center coordinate **c** to the output file. This file then has 8 columns: line number, **p**, **c**, *r*

**-info**: triggers storing runtime information (cpu times etc) in a separate [info](#info) file.

**-quiet**: prevents GPSD-3D to create stdout.

**-clean**: removes all temporary directories that may have been generated during previous crashes.

### Comments: 

1. The leading 'perl' is only needed, if your system does not automatically recognize GPSD-3D to be a perl file. If GPSD-3D is not found, call it via: perl ./GPSD-3D.
2. If called without any argument, GPSD-3D displays the documentation.
3. GPSD-3D can be called in parallel. 
4. Each GPSD-3D call runs in a unique temporary directory, upon successful completion the temporary directory is removed.
5. The maximum number of threads used by OpenMP is reported during the installation and also if GPSD-3D is called without arguments. For very large systems with, say, more than 1000000 spheres, running at the maximum number of threats must not be an advantage and you should check the speed also using a single processor, using -np=1
6. Negative values *r<sub>p</sub>* and *r<sub>c</sub>* are allowed as long as *r<sub>p</sub>*+*r<sub>c</sub>*+*r<sub>o</sub>* is positive. A [negative coating thickness](#about) effectively reduces the material sphere radius.


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

A number of configurations and corresponding box-files are available from the current respository. They are named .benchmark-x-config and .benchmark-x-box. A test call, using 10 of the available threads, 20000 Monte Carlo trials (*q=10*), for a probe sphere with zero radius, and *N=2000* materials spheres of radius *r*<sub>o</sub>=1.0 is 

        perl ./GPSD-3D -in=.benchmark-7-config -box=.benchmark-7-box -rp=0.0 -ro=1.0 -q=10 -np=10 -fortran

As we did not suppress stdout via -quiet, it should produce the following within a few seconds:

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
        [INFO] created files in .tmp-GPSD-3D-49303 including .parameters.
        [INFO] monodisperse system. The particle radius is taken as 1, shell thickness 0, test particle radius 0.
        [GPSD-3D] Using 20000 shots on 10 threads
        [GPSD-3D] Please stand by ..
        [GPSD-3D]                               reading box ..
        [GPSD-3D]                                          box      24.0000       24.0000       24.0000
        [VORO++]                                  max vertices           53
        [VORO++]                                     max faces           26
        [VORO++]                             max face-vertices           12
        [GPSD-3D]                                    triangles        59622
        [GPSD-3D]                      parallel processes (np)           10
        [GPSD-3D]                         material spheres (N)         2000
        [GPSD-3D]                              UpperPoreRadius       2.8465
        [GPSD-3D]                                           ro       1.0000
        [GPSD-3D]                                           rc       0.0000
        [GPSD-3D]                                           rp       0.0000
        [GPSD-3D]                                 reff = rc+rp       0.0000
        [GPSD-3D]                                 rs = ro+reff       1.0000
        [GPSD-3D]                       max triangle extension       2.8398
        [INFO]                  note: vertices inside material
        [GPSD-3D]                    creating neighbor list ..
        [GPSD-3D]                               neighborlist_M            4             4             4
        [GPSD-3D]                            neighborlist_size       6.0000        6.0000        6.0000
        [GPSD-3D]                           triangles per cell     931.5938
        [GPSD-3D]                      starting Monte Carlo ..
        [GPSD-3D]                            volume V=V(0,-ro)   13824.0000
        [GPSD-3D]                    volume fraction phi(reff)       0.3638
        [GPSD-3D]                                    V(0|reff)    8794.8288
        [GPSD-3D]                              min pore radius       0.0122
        [GPSD-3D]                             mean pore radius       1.5452 +/     0.0033
        [GPSD-3D]                              max pore radius       2.8465
        [GPSD-3D]             created a list {r} of pore radii
        [GPSD-3D]                    shots (use -q to enlarge)        20000
        [GPSD-3D]       cpu+real time spent in overhead [secs]       0.0003        0.0000
        [GPSD-3D]      cpu+real time spent in read_voro [secs]       0.0799        0.2500
        [GPSD-3D]cpu+real time spent in setup_triangles [secs]       0.0004        0.0000
        [GPSD-3D]     cpu+real time spent in MonteCarlo [secs]      14.3359        1.3750
        [GPSD-3D]         cpu+real time per 10000 shots [secs]       7.1680        0.6875
        [GPSD-3D] completed
        [GPSD-3D] created: .benchmark-7-config-ro=1-rp=0-rc=0.gpsd

and the following file (a list of roughly 20000 *r* values) should have been generated (if you do not see it, type: ls -lat): 

        .benchmark-7-config-ro=1-rp=0-rc=0.gpsd

With such list of values at hand, creating the normalized histogram (the pore radius distribution) is straightforward using any software that can bin the values, and visualize a graph. Some quantities derived from the list, such as minimum and maximum pore radius, as well as the mean pore radius including its standard error are mentioned in the above stdout. If you are not satisfied with the name of the resulting file, use the -o option. 

## -info file <a name="info">

If you repeat the above command, now using the -info option

        perl ./GPSD-3D -in=.benchmark-7-config -box=.benchmark-7-box -rp=0.0 -ro=1.0 -q=10 -np=10 -info -fortran

a second file will have been generated (all entries in this file are described in brackets)

        [GPSD-3D] created: .benchmark-7-config-ro=1-rp=0-rc=0.info

        N=2000                             (number of material spheres)
        ro=1.00000                         (material sphere radius)
        rp=0.00000                         (probe sphere radius)
        rc=0.00000                         (shell thickness)
        V=13824.00000                      (box volume)
        triangles=59622                    (total number of triangles)
        shots=20000                        (total number of shots)
        radii=12656                        (total number of pore radius values)
        min_pore_radius=0.01316            (minimum pore radius detected)
        max_pore_radius=2.84655            (maximum pore radius detected)
        mean_pore_radius=1.54575           (mean pore radius)
        stderr_pore_radius=0.00333         (standard error of the mean pore radius)
        phi_reff=0.36720                   (phi(reff))
        cells=64                           (neighbor list cells)
        threads=72                         (number of threads)
        voro_cpu_time=0.11232              (cpu time used to process the voro++ output [secs])
        voro_real_time=0.25000             (real time used to process the voro++ output [secs])
        triangles_setup_cpu_time=0.00046   (cpu time used to setup triangles [secs])
        triangles_setup_real-time=0.00000  (real time used to setup triangles [secs])
        MonteCarlo_cpu_time=29.61485       (cpu time spent during Monte Carlo [secs])
        MonteCarlo_real_time=0.50000       (real time spent during Monte Carlo [secs])
        walltime=0.736114025115967         (total wall time [secs])

## -inf file <a name="inf">

A text-free version of the .info-file is available in the .inf-file. 

## Polydisperse systems: Grid-based <a name="hardcoded">

For the case of polydisperse systems, the GPSD-3D script contains two lines that may be edited by users to increase or reduce the resolution further. The default setting is:  

        $min_delta_grid     = 0.005;      # USER-defined minimum grid spacing (in units of the effective particle radius ro+rc)
        $maxvoxels_grid     = 1000000;    # USER-defined upper limit for number of voxels 

## Benchmarks (fortran90 version) 

Benchmark configurations are available as .benchmark-#-config and .benchmark-#-box files. 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | *S* | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 312 | 30 | 99840 | 9389 | 1.36 | 1.67 | 0 s | 0.92 s | 
|2 | 1080 | 0.1 | 0 | 0 | 92 | 30 | 99360 | 40338 | 31.28 | 31.3 | 0.25 s | 0.73 s | 
|3 | 64 | 0.1 | 0 | 0 | 1562 | 30 | 99968 | 2488 | 1.37 | 1.38 | 0.13 s | 0.43 s | 
|4 | 100000 | 0.1 | 0 | 0 | 1 | 30 | 100000 | 4058940 | 2.47 | 3.86 | 14.88 s | 22.56 s | 
|5 | 139218 | 0.1 | 0 | 0 | 1 | 30 | 100000 | 3917190 | 1.38 | 2.14 | 17 s | 21.62 s | 
|6 | 3200 | 0.1 | 0 | 0 | 31 | 30 | 99200 | 39504 | 1.31 | 1.31 | 1 s | 2.91 s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 30 | 99200 | 39504 | 0.34 | 0.41 | 1 s | 1.36 s | 
|6 | 3200 | 1 | 0.1 | 0 | 31 | 30 | 99200 | 39504 | 0.38 | 0.41 | 1 s | 1.41 s | 
|6 | 3200 | 1 | 0.1 | 0.1 | 31 | 30 | 99200 | 39504 | 0.29 | 0.31 | 1 s | 1.26 s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 3 | 99200 | 39504 | 0.34 | 0.41 | 1 s | 3.39 s | 
|7 | 2000 | 0.1 | 0 | 0 | 50 | 30 | 100000 | 59622 | 2.57 | 3.75 | 0.25 s | 4.82 s | 
|8 | 10000 | 0.1 | 0 | 0 | 10 | 30 | 100000 | 307609 | 3.21 | 4.83 | 1.25 s | 8.87 s | 
|9 | 100000 | 0.1 | 0 | 0 | 1 | 30 | 100000 | 4057056 | 0.38 | 0.4 | 14.5 s | 17.42 s | 
|10 | 200 | 0.1 | 0 | 0 | 500 | 30 | 100000 | 8028 | 1.61 | 2.03 | 0 s | 0.6 s | 

## Benchmarks (c++ version)

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | *S* | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|




## About <a name="about">

Related publication (GPSD-3D): Comput. Phys. Commun. (2023) submitted

Related publication (GPSD-2D): Phys. Rev. E 107 (2023) 015307. Link: http://doi.org/DOI:10.1103/PhysRevE.107.015307

GPSD-3D Code available from: https://github.com/mkmat/CODE-GPSD-3D

GPSD-2D Code available from: https://github.com/mkmat/CODE-GPSD-2D
