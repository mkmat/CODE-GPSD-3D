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

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 312 | 30 | 99840 | 9389 | 1.36 | 1.67 | 0.13 s | 1.13 s | 11.31 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 92 | 30 | 99360 | 40338 | 31.28 | 31.3 | 0.13 s | 0.98 s | 9.85 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 1562 | 30 | 99968 | 2488 | 1.37 | 1.38 | 0.13 s | 0.66 s | 6.56 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 1 | 30 | 100000 | 4058940 | 2.47 | 3.86 | 14.63 s | 22.19 s | 221.94 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 1 | 30 | 100000 | 3917190 | 1.38 | 2.14 | 17 s | 21.94 s | 219.38 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 31 | 30 | 99200 | 39504 | 1.31 | 1.31 | 1 s | 3.23 s | 32.56 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 30 | 99200 | 39504 | 0.34 | 0.41 | 1.13 s | 1.58 s | 15.97 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0 | 31 | 30 | 99200 | 39504 | 0.37 | 0.41 | 1 s | 1.64 s | 16.56 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0.1 | 31 | 30 | 99200 | 39504 | 0.29 | 0.31 | 1 s | 1.48 s | 14.93 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 3 | 99200 | 39504 | 0.34 | 0.41 | 1 s | 3.56 s | 35.9 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 50 | 30 | 100000 | 59622 | 2.57 | 3.75 | 0.25 s | 5.08 s | 50.78 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 10 | 30 | 100000 | 307609 | 3.21 | 4.83 | 1.25 s | 9.08 s | 90.84 &mu;s | 
|9 | 100000 | 0.1 | 0 | 0 | 1 | 30 | 100000 | 4057056 | 0.38 | 0.4 | 14.13 s | 17.35 s | 173.5 &mu;s | 
|10 | 200 | 0.1 | 0 | 0 | 500 | 30 | 100000 | 8028 | 1.61 | 2.03 | 0 s | 0.83 s | 8.3 &mu;s | 


## Benchmarks (c++ version)

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total |  time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 312 | 1 | 99840 | -1 | 1.36 | 0 | 0 s | 288.44 s | 2889.03 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 92 | 1 | 99360 | -1 | 31.29 | 31.3 | 0 s | 9441.82 s | 95026.32 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 1562 | 1 | 99968 | -1 | 1.37 | 1.38 | 0 s | 146.93 s | 1469.81 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 1 | 1 | 100000 | -1 | 2.47 | 3.66 | 0 s | 963.47 s | 9634.66 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 1 | 1 | 100000 | -1 | 1.38 | 2.08 | 0.01 s | 484.15 s | 4841.45 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 31 | 1 | 99200 | -1 | 1.31 | 1.31 | 0 s | 2514.03 s | 25343.04 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 1 | 99200 | -1 | 0.34 | 0.41 | 0 s | 540.02 s | 5443.77 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0 | 31 | 1 | 99200 | -1 | 0.14 | 0.31 | 0 s | 539.67 s | 5440.21 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0.1 | 31 | 1 | 99200 | -1 | 0.08 | 0.21 | 0 s | 282.03 s | 2843.03 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 1 | 99200 | -1 | 0.34 | 0.41 | 0 s | 540.49 s | 5448.48 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 50 | 1 | 100000 | -1 | 2.58 | 3.65 | 0 s | 1912.21 s | 19122.1 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 10 | 1 | 100000 | -1 | 3.22 | 4.8 | 0 s | 3627.17 s | 36271.69 &mu;s | 
|9 | 100000 | 0.1 | 0 | 0 | 1 | 1 | 100000 | -1 | 0.38 | 0.4 | 0 s | 4905.32 s | 49053.22 &mu;s | 
|10 | 200 | 0.1 | 0 | 0 | 500 | 30 | 100000 | -1 | --- | --- | --- | --- | crashed | 

## Benchmarks (nlin version) 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 312 | 1 | 99840 | 0 | 1.24 | 1.67 | 0 s | 1624.64 s | 16272.41 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 92 | 1 | 99360 | 0 | 30.82 | 31.3 | 0 s | 955.19 s | 9613.41 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 1562 | 1 | 99968 | 0 | 1.36 | 1.38 | 0 s | 96.23 s | 962.62 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 1 | 1 | 100000 | 0 | 2.26 | 3.86 | 0 s | 1053.45 s | 10534.46 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 1 | 1 | 100000 | 0 | 1.26 | 2.14 | 0 s | 1818.53 s | 18185.35 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 31 | 1 | 99200 | 0 | 1.31 | 1.31 | 0 s | 2646.88 s | 26682.26 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 1 | 99200 | 0 | 0.34 | 0.41 | 0 s | 2201.02 s | 22187.72 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 50 | 1 | 100000 | 0 | 2.36 | 3.75 | 0 s | 3259.57 s | 32595.74 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 10 | 1 | 100000 | 0 | 2.96 | 4.83 | 0 s | 3896.62 s | 38966.25 &mu;s | 
|9 | 100000 | 0.1 | 0 | 0 | 1 | 1 | 100000 | -1 | --- | --- | --- | --- | crashed | 
|10 | 200 | 0.1 | 0 | 0 | 500 | 1 | 100000 | 0 | 1.46 | 2.03 | 0 s | 894.89 s | 8948.91 &mu;s | 


## FEW SHOTS single processor Benchmarks (fortran90 version) 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 3 | 1 | 960 | 9389 | 1.35 | 1.67 | 0 s | 0.51 s | 534.34 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 1 | 1 | 1000 | 40338 | 31.27 | 31.3 | 0.25 s | 0.57 s | 574.41 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 15 | 1 | 960 | 2488 | 1.37 | 1.38 | 0 s | 0.33 s | 344.23 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 4058940 | 2.47 | 3.86 | 14.63 s | 16.34 s | 16339.69 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 3917190 | 1.38 | 2.14 | 17 s | 18.43 s | 18429.47 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 39504 | 1.31 | 1.31 | 1 s | 1.8 s | 1795.71 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 0 | 1 | 1000 | 39504 | 0.34 | 0.41 | 1 s | 1.39 s | 1389.76 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0 | 0 | 1 | 1000 | 39504 | 0.36 | 0.41 | 1.13 s | 1.4 s | 1400.5 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0.1 | 0 | 1 | 1000 | 39504 | 0.29 | 0.31 | 1 s | 1.35 s | 1349.7 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 1 | 1 | 1000 | 59622 | 2.56 | 3.75 | 0.25 s | 1.78 s | 1780.73 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 307609 | 3.22 | 4.83 | 1.25 s | 3.15 s | 3154.72 &mu;s | 
|9 | 100000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 4057056 | 0.38 | 0.4 | 14.13 s | 15.3 s | 15295.63 &mu;s | 
|10 | 200 | 0.1 | 0 | 0 | 5 | 1 | 1000 | 8028 | 1.61 | 2.03 | 0 s | 0.43 s | 426.67 &mu;s | 


## FEW SHOTS single processor Benchmarks (c++ version) 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 3 | 1 | 960 | -1 | 1.36 | 1.67 | 0 s | 3.02 s | 3145.32 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 1 | 1 | 1000 | -1 | 31.25 | 31.3 | 0 s | 102.97 s | 102965.45 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 15 | 1 | 960 | -1 | 1.37 | 1.38 | 0 s | 1.64 s | 1707.8 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | -1 | 2.48 | 3.35 | 0 s | 23.18 s | 23180.67 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 0 | 1 | 1000 | -1 | 1.39 | 2 | 0.01 s | 17.12 s | 17119.82 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 0 | 1 | 1000 | -1 | 1.31 | 1.31 | 0 s | 26.05 s | 26045.12 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 0 | 1 | 1000 | -1 | 0.35 | 0.41 | 0 s | 5.35 s | 5350.34 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0 | 0 | 1 | 1000 | -1 | 0.15 | 0.31 | 0 s | 5.34 s | 5342.38 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0.1 | 0 | 1 | 1000 | -1 | 0.1 | 0.21 | 0 s | 3.17 s | 3170.17 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 0 | 1 | 1000 | -1 | 0.35 | 0.41 | 0 s | 5.36 s | 5357.05 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 1 | 1 | 1000 | -1 | 2.59 | 3.44 | 0 s | 22.43 s | 22428.1 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | -1 | 3.19 | 4.83 | 0 s | 39.1 s | 39101.76 &mu;s | 
|9 | 100000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | -1 | 0.4 | 0 | 0 s | 42.95 s | 5357.05 &mu;s | 
|10 | 200 | 0.1 | 0 | 0 | 5 | 1 | 1000 | 0 | --- | --- | --- | --- | crashed | 


## FEW SHOTS single processor Benchmarks (nlin version) 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 3 | 1 | 960 | 0 | 1.23 | 1.67 | 0 s | 5.07 s | 5280.63 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 1 | 1 | 1000 | 0 | 30.83 | 31.3 | 0 s | 8.36 s | 8361.14 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 15 | 1 | 960 | 0 | 1.36 | 1.38 | 0 s | 0.58 s | 605.69 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 0 | 2.26 | 3.44 | 0 s | 4.25 s | 4251.54 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 0 | 1.26 | 2.14 | 0 s | 6.1 s | 6096.2 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 0 | 1.31 | 1.31 | 0 s | 18.99 s | 18986.37 &mu;s |
|6 | 3200 | 1 | 0 | 0 | 0 | 1 | 1000 | 0 | 0.35 | 0.41 | 0 s | 20.61 s | 20612.61 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 1 | 1 | 1000 | 0 | 2.38 | 3.75 | 0 s | 12.11 s | 12113.66 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 0 | 1 | 1000 | 0 | 3.02 | 4.82 | 0 s | 19.09 s | 19090.51 &mu;s | 
|10 | 200 | 0.1 | 0 | 0 | 5 | 1 | 1000 | 0 | 1.46 | 2.03 | 0 s | 2.78 s | 2780.57 &mu;s | 

## About <a name="about">

Related publication (GPSD-3D): Comput. Phys. Commun. (2023) submitted

Related publication (GPSD-2D): Phys. Rev. E 107 (2023) 015307. Link: http://doi.org/DOI:10.1103/PhysRevE.107.015307

GPSD-3D Code available from: https://github.com/mkmat/CODE-GPSD-3D

GPSD-2D Code available from: https://github.com/mkmat/CODE-GPSD-2D
