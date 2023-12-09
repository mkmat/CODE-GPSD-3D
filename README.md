# CODE-GPSD-3D
Generalized geometric pore size distribution (G-PSD) for periodic systems composed of spheres. This quantity was defined by Gelb and Gubbins, there are other pore size distributions such as T-PSD, as discussed in detail in the accompanying [article](#about). The spheres constitute the 'material', which is surrounded by 'pore space'. This codes allows to calculate the distribution *P(r;r<sub>p</sub>|r<sub>c</sub>)*  of pore radii *r* for such system as function of the radius *r<sub>p</sub>* of a thought probe particle and a thought coating thickness *r<sub>c</sub>* of the material spheres. 

<img src="images/schematic-GPSD3D-github.png" width="100%">  

For monodisperse systems GPSD-3D uses the advanced (grid-free) voronoi-based algorithm. For polydisperse systems, it uses a constrained nonlinear optimization strategy or alternatively a grid-based algorithm, whose resolution is limited by the amount of available memory. See [below](#polydisperse) for details for polydisperse systems. Lammps-users can call GPSD-3D directly from within their lammps script as shown [here](#lammps).

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

But, if you have voro++ not yet installed, or no fortran compiler, you will end up with an error message that provides you with a download link. In this case you need to do the installation, and call the installation script INSTALL.pl again. Furthermore, if you'd like to have GPSD-3D offering the -nlin option, the INSTALL.pl will ask you to install NLopt. You may choose to install GPSD-3D without the -nlin option. 

## Configuration and box file formats <a name=format>

The input required by GPSD-3D are (i) coordinates: the center positions of *N* monodisperse or polydisperse spheres, that are located inside a rectangular, periodic box, whose corners are specified by (ii) box geometry: 6 values: xlo,xhi,ylo,yhi,zlo,zhi. We offer various file formats, where each of the *N* lines contains the coordinates, and eventually also the radius of a single sphere

1. x y z (monodisperse system, requires specifying -ro on the command line)
2. id x y z (monodisperse system, requires specifying -ro on the command line)        
3. id x y z radius  (polydisperse system, monodisperse if all radii are equal)
4. samarth-type configuration file (do not specify box dimensions in that case, requires specifying -ro on the command line)
5. lammps-users can call GPSD-3D from within their script as shown [here](#lammps)

The six values for the box can be either saved in a txt-file (single line, six values xlo xhi ylo yhi zlo zhi separated by blank or commata), or passed over on the command line. The delimiting character can be specified using the -d option.  


## How to run the code. Parameters. Command-line options

      perl GPSD-3D
          -in=<filename>
          -box=<filename> OR -xlo=.. -xhi=.. ylo=.. -yhi=.. -zlo=.. -zhi=..
          [-ro=<positive value>] 
          [-rp=<value>] 
          [-rc=<value>] 
          [-q=<integer>] 
          [-TPSD] 

          [-more]
          [-info]
          [-np=<integer>]
          [-d=<delimiter>]

          [-grid]
          [-griddelta=..]
          [-gridmax=..]

          [-c++]

          [-nlin]
          [-kmpr=..]
          
          [-o=<filename>]
          [-list]
          
          [-quiet]
          [-clean]

**-in=**    name of the file containing the material sphere [coordinates](#format) (for polydisperse systems also the sphere radii).

**-box=**   name of the file containing the box [geometry](#format) (alternatively, the box size can be passed over on the command line using -xlo= .. -zhi=..).

**-ro=** material circle radius *r<sub>o</sub>* (required, if the particle radii are not contained in the input file). If *r<sub>o</sub>* is specified, existing radii in the input file are ignored).

**-rp=** probe particle radius *r<sub>p</sub>*. If not specified, *r<sub>p</sub>=0* is used.

**-rc=** shell thickness  *r<sub>c</sub>*. If not specified, *r<sub>c</sub>=0* is used.

**-q=** positive quality value (optionally). If not specified, *q=10.0* is used. The number of random shots is *q* times the number of material spheres.

**-TPSD** calculate the TPSDs in addition create a *.tpsd file. 

**-more**: tell the code to add, besides pore radius *r*, the corresponding probe particle center **p** and pore center coordinate **c** to the output file. This file then has 8 columns: line number, **p**, **c**, *r*

**-info**: triggers storing runtime information (cpu times etc) in a separate [info](#info) file.

**-np=** number of threads np to be used. If not specified, *n<sub>p</sub>* is set to the available number of threads.

**-d=** delimiter (single character) present in the configuration and box-files, e.g. -d="\ " for a blank, default is -d=","

**-grid** enforce using the grid-based method

**-griddelta=..** minimum grid spacing (in units of effective particle radius ro+rc) (default: 0.005)

**-gridmax=..** maximum number of grid voxels (default: 1000000)

**-c++** enforce using the c++ version (default: voronoi fortran version)

**-nlin** enforce using the constrained nonlinear optimization [if installed]

**-kmpr=..** specify, if known, the maximum pore radius to speed up the -nlin algorithm

**-o=** name of the resulting file containing a list of pore radii. If not specified, the list is saved in a .gpsd-file.

**-list** use a list of p vectors (stored in p.txt, format: id px py pz) instead of randomly shooting

**-quiet**: prevents GPSD-3D to create stdout.

**-clean**: removes all temporary directories that may have been generated during previous crashes.

### Comments: 

1. The leading 'perl' is only needed, if your system does not automatically recognize GPSD-3D to be a perl file. If GPSD-3D is not found, call it via: perl ./GPSD-3D.
2. If called without any argument, GPSD-3D displays the documentation.
3. GPSD-3D can be called in parallel. 
4. Each GPSD-3D call runs in a unique temporary directory, upon successful completion the temporary directory is removed.
5. The maximum number of threads used by OpenMP is reported during the installation and also if GPSD-3D is called without arguments. For very large systems with, say, more than 1000000 spheres, running at the maximum number of threats must not be an advantage and you should check the speed also using a single processor, using -np=1
6. Negative values *r<sub>p</sub>* and *r<sub>c</sub>* are allowed as long as *r<sub>p</sub>*+*r<sub>c</sub>*+*r<sub>o</sub>* is positive. A [negative coating thickness](#about) effectively reduces the material sphere radius.
7. Windows-users: Your perl may not accept arguments on the command line. If so, see [windows-user](#windows)

## Output

GPSD-3D returns a list of pore radii *r* in a file, either in configfilename.gpsd or outputfilename, if the latter had explicitly been defined on the command line. 

        0.658312
        0.274754
        1.070546
        0.685289
        ...

This list of *r* values (for the chosen values *r<sub>p</sub>* and *r<sub>c</sub>*) gives rise to a distribution of pore radii, so called generalized geometric pore radius distribution *P(r;r<sub>p</sub>|r<sub>c</sub>)*. The bare radius of the particles *r<sub>o</sub>* is usually not mentioned here, as it belongs to the system. For monodisperse systems only the sum or *r<sub>o</sub>+r<sub>c</sub>* matters. For polydisperse systems each spherical particle has its own radius according to the configuration file, and *r<sub>c</sub>* can be used to effectively modify the stored particle radii, without changing the configuration file. If you call GPSD-3D with the -more option, the same file will contain seven columns (no header)

        px py pz cx cy cz r

where *px*,*py*,*pz* are the coordinates of the shot into the void space that gave rise to the center coordinates *cx*, *cy*, *cz* of the pore with radius *r*.

## Test configurations and test runs

A number of configurations and corresponding box-files are available from the current respository. They are named .benchmark-x-config and .benchmark-x-box. A test call, using 10 of the available threads, 20000 Monte Carlo trials (*q=10*), for a probe sphere with zero radius, and *N=2000* materials spheres of radius *r*<sub>o</sub>=1.0 is 

        perl ./GPSD-3D -in=.benchmark-7-config -box=.benchmark-7-box -rp=0.0 -ro=1.0 -q=10 -np=10 

As we did not suppress stdout via -quiet, it should produce the following within a few seconds:

      
            _______________________________________________________________________________________________________________________________
            
            This is GPSD-3D version 1.0
            written Nov 2023 by Martin Kroger, ETH Zurich, www.complexfluids.ethz.ch, Email: mk@mat.ethz.ch
            
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
            [INFO] using ro=1, rc=0, rp=0
            [INFO] created files in .tmp-GPSD-3D-733750 including .parameters.
            [INFO] monodisperse system. The particle radius is taken as 1, shell thickness 0, test particle radius 0.
            [GPSD-3D] using voronoi (fortran) algorithm[GPSD-3D] Using 20000 shots on 10 threads
            [GPSD-3D] Please stand by ..
            [GPSD-3D]                               reading box ..
            [GPSD-3D]                                          box      24.0000       24.0000       24.0000
            [VORO++]                                             N         2000
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
            [GPSD-3D]                  creating T-neighbor list ..
            [GPSD-3D]                              Tneighborlist_M            4             4             4
            [GPSD-3D]                           Tneighborlist_size       6.0000        6.0000        6.0000
            [GPSD-3D]                           triangles per cell     931.5938
            [GPSD-3D]                  creating X-neighbor list ..
            [GPSD-3D]                              Xneighborlist_M           24            24            24
            [GPSD-3D]                           Xneighborlist_size       1.0000        1.0000        1.0000
            [GPSD-3D]                           particles per cell       0.1447
            [GPSD-3D]                      starting Monte Carlo ..
            [GPSD-3D]                            volume V=V(0,-ro)   13824.0000
            [GPSD-3D]                    volume fraction phi(reff)       0.3667
            [GPSD-3D]                                    V(0|reff)    8754.7392
            [GPSD-3D]                              min pore radius       0.0147
            [GPSD-3D]                             mean pore radius       1.5456 +/     0.0033
            [GPSD-3D]                              max pore radius       2.8465
            [GPSD-3D]                    pore radius {r} list size        12666
            [GPSD-3D]                    shots (use -q to enlarge)        20000
            [GPSD-3D]       cpu+real time spent in overhead [secs]       0.0003        0.0000
            [GPSD-3D]      cpu+real time spent in read_voro [secs]       0.0770        0.2500
            [GPSD-3D]cpu+real time spent in setup_triangles [secs]       0.0004        0.0000
            [GPSD-3D]     cpu+real time spent in MonteCarlo [secs]       7.3443        0.7500
            [GPSD-3D]         cpu+real time per 10000 shots [secs]       3.6722        0.3750
            [GPSD-3D] completed GPSD
            [GPSD-3D] created: .benchmark-7-config-ro=1-rp=0-rc=0.gpsd

and the following file (a list of roughly 13000 *r* values) should have been generated (if you do not see it, type: ls -lat): 

        .benchmark-7-config-ro=1-rp=0-rc=0.gpsd

With such list of values at hand, creating the normalized histogram (the pore radius distribution) is straightforward using any software that can bin the values, and visualize a graph. Some quantities derived from the list, such as minimum and maximum pore radius, as well as the mean pore radius including its standard error are mentioned in the above stdout. If you are not satisfied with the name of the resulting file, use the -o option. 

## -info file <a name="info">

If you repeat the above command, now using the -info option

        perl ./GPSD-3D -in=.benchmark-7-config -box=.benchmark-7-box -rp=0.0 -ro=1.0 -q=10 -np=10 -info

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

## Polydisperse systems<a name="polydisperse">

While GPSD-3D can handle polydisperse conifigurations, we do not recommend using it, as the Voronoi-based method cannot be used directly; GPSD-3D falls back using the costrained nonlinear optimization or grid-based method. While the constrained nonlinear optimization is slow and must not produce absolutely correct results, the grid-based method is memory-consuming and slow as well. A future release of GPSD-3D will be able to handle polydisperse systems quickly. 

### Constrained nonlinear optimization

For the case of polydisperse systems, by default the GPSD-3D script employs a constrained nonlinear solver. There is one option that may be used to speed up the solver if an upper limit for the pore radius is already known

        -nlin
        -kmpr=5.3              # USER-defined upper limit of the pore radius (here: 5.3)

### Grid-based <a name="hardcoded">

For the case of polydisperse systems, a grid-based solver is used if the -grid option is given. There are two options, that may be used to increase or reduce the resolution further. The default setting is:  

        -grid 
        -griddelta=0.005       # USER-defined minimum grid spacing (in units of the effective particle radius ro+rc)
        -gridmax=1000000       # USER-defined upper limit for the number of voxels 

## lammps-users<a name="lammps">

There are several ways to use GPSD-3D with lammps. 

### Version 1

Call GPSD-3D directly from within your lammps script. GPSD-3D uses by default the box dimensions provided by lammps, no need to specify them in your GPSD-3D call.

        dump           GPSDdump all custom 10000 lammps.dump id x y z
        run            10000
        shell          perl ./GPSD-3D -in=lammps.dump -rp=.. -ro=.. -rc=.. -q=.. -info

The above example creates a single lammps.gpsd file at the end of your run, after 10000 time steps. You have still the freedom to change box dimensions on the command line via -xlo=.. -xhi=.. etc. Because the -info option has been given, a lammps.info and lammps.inf file are created as well.

### Version 2

If you want to calculate the PSD in regular intervals during a lammps simulation, create a lammps loop. 

        label         loop 
        variable      a loop 15
        dump          GPSDdump all custom 1000 lammps.dump id x y z
        run           1000
        shell         perl ./GPSD-3D -in=lammps.dump -rp=.. -ro=.. -rc=.. -q=.. -info -o=result${a} 
        # here you may call your own script that analyses result.gpsd or result.info 
        undump        GPSDdump
        next          a
        jump          SELF loop      
        label         break
        variable      a delete

where SELF should be replaced by the file name of your lammps script, if you running lammps in parallel. In the above example, the GPSD is calculated 15 times, each 1000 steps, and the resulting pore radii are saved in result1.gpsd, result2.gpsd, ..., result15.gpsd. If you prefer to not save all pore radii, but instead accumulate a histogram, you can replace result${a} by result, write a script that reads the result.gpsd file, calculates and accumulates the histogram, and saves it. Using the -info option, result1.info etc. files are also generated. 

### Version 3 (not recommended)

Store a xyz-formatted file using lammps-commands

        dump           GPSDdump all xyz 1000 lammps.xyz
        dump_modify    GPSDdump header no

Then call GPSD-3D from the command line, after lammps finished. Because the xyz-formmated lammps file does not contain box sizes, you have the freedom to specify them manually. 

        perl ./GPSD-3D -in=lammps.xyz -xlo=... -xhi=... -ylo=... -yhi=... -zlo=... -zhi=... -rp=... -ro=... -q=...

Be aware that the lammps-formatted xyz-file may contain scaled coordinates, so that you have to specify ro, rp, and rc also in scaled form.

## Benchmarks (100000 shots, using -np=30) 

Benchmark configurations are available as .benchmark-#-config and .benchmark-#-box files. 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 312 | 30  | 9389 | 1.36 | 1.67 | 0 s | 0.85 s | 8.5 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 92 | 30  | 40338 | 31.28 | 31.3 | 0.13 s | 0.69 s | 7.0 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 1562 | 30  | 2488 | 1.37 | 1.38 | 0 s | 0.42 s | 4.2 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 1 | 30  | 4058940 | 2.47 | 3.86 | 14.25 s | 20.92 s | 209.3 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 1 | 30  | 3917190 | 1.38 | 2.14 | 16.63 s | 20.39 s | 203.9 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 31 | 30  | 39504 | 1.31 | 1.31 | 1 s | 2.74 s | 27.6 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 30  | 39504 | 0.34 | 0.41 | 1 s | 1.31 s | 13.3 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0 | 31 | 30  | 39504 | 0.37 | 0.41 | 1 s | 1.37 s | 13.8 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0.1 | 31 | 30  | 39504 | 0.29 | 0.31 | 1 s | 1.21 s | 12.2 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 50 | 30  | 59622 | 2.58 | 3.75 | 0.25 s | 4.46 s | 44.6 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 10 | 30  | 307609 | 3.21 | 4.83 | 1.13 s | 8.2 s | 82.0 &mu;s | 
|9 | 100000 | 0.1 | 0 | 0 | 1 | 30  | 4057056 | 0.38 | 0.4 | 13.88 s | 16.73 s | 167.3 &mu;s | 
|10 | 200 | 0.1 | 0 | 0 | 500 | 30  | 8028 | 1.61 | 2.03 | 0.13 s | 0.58 s | 5.8 &mu;s | 
|11 | 64 | 0.1 | 0 | 0 | 1562 | 30  | 768 | 1.63 | 1.63 | 0 s | 0.41 s | 4.1 &mu;s | 
|12 | 192 | 0.1 | 0 | 0 | 520 | 30  | 5880 | 3.64 | 3.91 | 0 s | 0.47 s | 4.7 &mu;s | 


## Benchmarks (100000 shots, using -np=1) 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1 | 320 | 0.1 | 0 | 0 | 312 | 1  | 9389 | 1.36 | 1.67 | 0 s | 21.41 s | 214.5 &mu;s | 
|2 | 1080 | 0.1 | 0 | 0 | 92 | 1  | 40338 | 31.28 | 31.3 | 0.25 s | 12.18 s | 122.6 &mu;s | 
|3 | 64 | 0.1 | 0 | 0 | 1562 | 1  | 2488 | 1.37 | 1.38 | 0 s | 6.3 s | 63.0 &mu;s | 
|4 | 100000 | 0.1 | 0 | 0 | 1 | 1  | 4058940 | 2.47 | 3.86 | 14.25 s | 120.87 s | 1208.7 &mu;s | 
|5 | 139218 | 0.1 | 0 | 0 | 1 | 1  | 3917190 | 1.38 | 2.14 | 16.5 s | 74.7 s | 747.0 &mu;s | 
|6 | 3200 | 0.1 | 0 | 0 | 31 | 1  | 39504 | 1.31 | 1.31 | 1 s | 49.85 s | 502.5 &mu;s | 
|6 | 3200 | 1 | 0 | 0 | 31 | 1  | 39504 | 0.34 | 0.41 | 1 s | 7.52 s | 75.8 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0 | 31 | 1  | 39504 | 0.37 | 0.41 | 1 s | 9.38 s | 94.6 &mu;s | 
|6 | 3200 | 1 | 0.1 | 0.1 | 31 | 1  | 39504 | 0.29 | 0.31 | 1 s | 4.97 s | 50.1 &mu;s | 
|7 | 2000 | 0.1 | 0 | 0 | 50 | 1  | 59622 | 2.58 | 3.75 | 0.25 s | 122.89 s | 1228.9 &mu;s | 
|8 | 10000 | 0.1 | 0 | 0 | 10 | 1  | 307609 | 3.21 | 4.83 | 1.13 s | 164.99 s | 1649.9 &mu;s | 
|9 | 100000 | 0.1 | 0 | 0 | 1 | 1  | 4057056 | 0.39 | 0.4 | 13.88 s | 64.3 s | 643.0 &mu;s | 
|10 | 200 | 0.1 | 0 | 0 | 500 | 1  | 8028 | 1.61 | 2.03 | 0.13 s | 13.08 s | 130.8 &mu;s | 
|11 | 64 | 0.1 | 0 | 0 | 1562 | 1  | 768 | 1.63 | 1.63 | 0.13 s | 0.5 s | 5.0 &mu;s | 
|12 | 192 | 0.1 | 0 | 0 | 520 | 1  | 5880 | 3.64 | 3.91 | 0.13 s | 6.35 s | 63.7 &mu;s | 



## Benchmarks (c++ version)

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total |  time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|


## Benchmarks (nlin version) 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

## Benchmarks (grid version) 

| no | *N* | *r<sub>o</sub>* | *r<sub>p</sub>* | *r<sub>c</sub>* | *q* | *n<sub>p</sub>* | shots | triangles | $\langle r\rangle$ | *r<sub>max</sub>* | voro++ | total | time/shot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

## Windows-users<a name="windows"></a>

To our surprize, we successfully installed voro++ and GPSD-3D under windows 10 using /mingw64/bin/gfortran and /mingw64/bin/g++. You need to have perl installed. After successful installation, we were faced with the problem that our perl under windows did not accept command line arguments. To circumvent this problem, add the full GPSD-3D command to a file, say, windows.pl and execute this file from the command line via: 

        perl ./windows.pl

An example windows.pl file comes with the GPSD-3D distribution. 

## About <a name="about"></a>

Related publication (GPSD-3D): Comput. Phys. Commun. (2023) submitted

Related publication (GPSD-2D): Phys. Rev. E 107 (2023) 015307. Link: http://doi.org/DOI:10.1103/PhysRevE.107.015307

GPSD-3D Code available from: https://github.com/mkmat/CODE-GPSD-3D

GPSD-2D Code available from: https://github.com/mkmat/CODE-GPSD-2D
