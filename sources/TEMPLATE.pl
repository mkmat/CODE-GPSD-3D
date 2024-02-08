#! /usr/bin/perl

# (c) 28 sep 2023 mk@mat.ethz.ch
# GPSD-3D had been generated from TEMPLATE.pl

$installation_dir   = "InstallationDirectory";  # set by INSTALL.pl  

$code               = "$installation_dir/GPSD-3D";
$refcode            = "$installation_dir/GPSD-3D-grid.ex";
$nlincode           = "$installation_dir/GPSD-3D-nlin.ex"; 
$voroparser0        = "$installation_dir/voro-parser-0.ex";
$voroparser         = "$installation_dir/voro-parser.ex";
$nlopt_build        = "$installation_dir/nlopt-build"; 
$converter1         = "$installation_dir/convert-voro-output-to-triangles.pl";
$converter2         = "$installation_dir/convert_samarth_to_config.pl";
$maxnp              = `$installation_dir/.maxnp.ex`+0;
$zeocode            = "$installation_dir/ZEO++/zeoplusplus-master/network"; 
$MESSAGE            = "$installation_dir/MESSAGE.txt";
$pwd                = `pwd`; chomp $pwd; 

if (-d "$nlopt_build") { } else { $option_nlin = " [not installed]"; }; 
if (-d "$code-cpp.ex") { } else { $option_cpp  = " [not installed]"; }; 

# ------------------------------------------------
sub PRINT   { if ($quiet eq "true") { } else { print "$_[0]"; }; };
sub green   { "\033[1;32;40m$_[0]\033[1;37;40m"; };
sub red     { "\033[1;31;40m$_[0]\033[1;37;40m"; };
sub orange  { "\033[1;33;40m$_[0]\033[1;37;40m"; };
sub showerr { PRINT(red("________________________________________________\n\nERROR: $ERROR\n________________________________________________ cd $workdir\n")); exit; };
sub strip   { chomp $_[0]; $_[0]=~s/\s+/ /g; $_[0]=~s/,/ /g; $_[0]=~s/$delimiter/ /g; $_[0]=~s/^ //; $_[0]; };

# ------------------------------------------------
sub USAGE { print green("perl $0 -in=<configfile> [-box=<boxfile>] [-rp=<rp>] [-ro=<ro>] [-rc=<rc>] [-q=<integer>] [-o=<outputfile>] [-more] [-info] [-quiet] [-np=<integer>]\n");
print "or\n";
print green("perl $0 -clean\n");
print<<EOF;
\nwhere 

-in=..      : name of a configuration file (contains as many lines as material particles, three to five columns)
              (A) possible file format: x y z            (monodisperse system)
              (B) possible file format: id x y z         (monodisperse system)
              (C) possible file format: id x y z radius  (polydisperse system, monodisperse if all radii are equal)
              (D) samarth-type configuration file (do not specify -boxfile in this case)
              (E) lammps.dump file (no need to specify a -boxfile in that case)  
              (F) lammps.dump trajectory file (no need to specify a -boxfile in that case, choose time step via the -dumpstep=.. option)
-box=..     : name of a file that contains the box bounds: xlo xhi ylo yhi zlo zhi on a single line
              Alternatively, box bounds can be entered on the command line via -xlo=.. -xhi=.. -ylo=.. -yhi=.. -zlo=.. -zhi=..
-ro=..      : particle radius (required for formats A and B) (if set and if system is monodisperse, overwrites radius in file format C. For polydisperse systems see note 6) 
-rp=..      : probe particle radius (default: rp=0)
-rc=..      : particle coating thickness (default: rc=0) (only taken into account if file format C is used, ro=radius+rc is used)
-q=..       : A positive number (default: 10), determines the quality. The number of shots is q times the number of material spheres
-TPSD       : calculate T-PSDs in addition

-more       : if -more is added, not only the pore radius r (column 1) but also the center of the pore (columns 2-4) is added to outputfile.
-info       : if -info is added, runtime information is collected in files whose names end with .info and .inf (more information, see github repository)
-np=..      : maximum number of parallel processes to be used (between 1 and $maxnp)
-d=..       : delimiter in configuration and box-files, default: blank. Example -d="," for comma-delimited data.

-grid       : enforce using the grid-based method
-griddelta  : minimum grid spacing (in units of effective particle radius ro+rc) (default: 0.005)
-gridmax    : maximum number of grid voxels (default: 1000000)

-c++        :$option_cpp enforce using the c++ version (default: voronoi fortran version)
-nlin       :$option_nlin enforce using the constrained nonlinear optimization
-kmpr       :$option_nlin specify, if known, the maximum pore radius to speed up the -nlin algorithm

-o=..       : name of the resulting file containing a list of pore radii (if not specified, resulting file ends with .gpsd)
-list       : use a list of p vectors (stored in p.txt, format: id px py pz) instead of randomly shooting
-dumpstep=..: if the configuration file is a lammps dump-trajectory file, choose a snapshot at given time step
-quiet      : if -quiet is added, GPSD-3D won't produce stdout.
-clean      : remove temporary directories .tmp-GPSD*, if existing (such directories may exist after a crash or interrupt).


Note: 

1) This software and its documentation is hosted at https://github.com/mkmat/CODE-GPSD-3D
2) GPSD-3D can be called in parallel and sent to background. Each instance runs in a separate temporary directory.
3) All temporary directories can be viewed via: ls -1 .tmp-GPSD*; 
4) Particle coordinates that do not fall into the simulation box are just ignored.
5) The coating thickness can be chosen to be negative, while ro+rc should be semipositive.
6) For polydisperse systems, GPSD-3D uses -nlin by default. See https://github.com/mkmat/CODE-GPSD-3D for more information.
7) Windows-user: Your perl may not accept command line arguments. If so, edit windows.pl and call it via: perl ./windows.pl

Resulting files: 

<configfile>.gpsd         : contains a list of pore radii (monodisperse systems, voronoi- or grid-based G-PSD)
<configfile>.tpsd         : contains a list of pore radii (T-PSDs)
<configfile>.info         : -info file (human readable)
<configfile>.inf          : -info file (numbers only)

Reproduce the example mentioned at our github site using a benchmark configuration available in your installation directory:

perl ./GPSD-3D -in=.benchmark-7-config -box=.benchmark-7-box -rp=0.0 -ro=1.0 -q=10 -np=10 

Examples using other existing benchmarks
EOF
$benchmarks = `ls -1 .bench*config | grep -c .`+0; 
foreach $no (0 .. $benchmarks) { 
   if (-s ".benchmark-$no-config") { 
      PRINT("perl $0 -in=.benchmark-$no-config -box=.benchmark-$no-box -rp=0.0 -ro=1.0 -q=10 [-nlin -info -more -quiet]\n");
   };
}; 
if ($ERROR) { showerr; }; 
exit;
};

# init walltime clock
$walltime_tic=`date +"%s.%N"`+0; 

# defaults
$rc               = 0; 
$ro               = 0; 
$rp               = 0; 
$q                = 10; 
$box              = 0; 
$more             = "false"; 
$fortran          = "true"; 
$cpp              = "false";
$TPSD             = "false"; 
$np               = 1; 
$gridbased        = "false"; 
$nlin             = "false"; 
$zeo              = "false"; 
$griddelta        = 0.005; 
$maxvoxels_grid   = 1000000; 
$delimiter        = " ";
$use_shots_list   = "false";
$known_max_pore_radius = -1; 
$dumpstep         = -1; 

# check for -quiet
$quiet="false"; foreach $arg (@ARGV) { if ($arg eq "-quiet") { $quiet="true"; } }; 

if ("$quiet" eq "false") { print `cat $MESSAGE`; };

# scan arguments
foreach $arg (@ARGV) { $arg=~s/=/==/; ($field,$value)=split(/==/,$arg);
    if ($field eq "-in")            { $config="$value";   PRINT("[INFO] using configuration file $config\n"); 
    } elsif ($field eq "-box")      { $boxfile="$value";  PRINT("[INFO] using box file $boxfile\n"); 
    } elsif ($field eq "-rp")       { $rp=$value+0; 
    } elsif ($field eq "-ro")       { $ro=$value+0; $setro=1; 
    } elsif ($field eq "-rc")       { $rc=$value+0; 
    } elsif ($field eq "-q")        { $q=$value+0; if ($q eq 0) { PRINT("[INFO] q=0 means: no shots\n"); }; 
    } elsif ($field eq "-clean")    { `rm -rf .tmp-GPSD-3D-*; rm -f .UpperPoreRadius .benchmark*.inf .benchmark*.info .benchmark*psd `; PRINT("[INFO] cleaned\n"); exit; 
    } elsif ($field eq "-quiet")    { $quiet="true"; 
    } elsif ($field eq "-more")     { $more="true"; 
    } elsif ($field eq "-info")     { $info=1; 
    } elsif ($field eq "-xlo")      { $xlo=$value+0; $box+=1; 
    } elsif ($field eq "-xhi")      { $xhi=$value+0; $box+=1; 
    } elsif ($field eq "-ylo")      { $ylo=$value+0; $box+=1; 
    } elsif ($field eq "-yhi")      { $yhi=$value+0; $box+=1; 
    } elsif ($field eq "-zlo")      { $zlo=$value+0; $box+=1; 
    } elsif ($field eq "-zhi")      { $zhi=$value+0; $box+=1; 
    } elsif ($field eq "-o")        { $outputfile="$value"; 
    } elsif ($field eq "--v")       { USAGE;
    } elsif ($field eq "-h")        { USAGE; 
    } elsif ($field eq "-np")       { $np=$value; if ($np<1) { $np=1; }; if ($np>$maxnp) { $np=$maxnp; }; $np_set=1; 
    } elsif ($field eq "-TPSD")     { $TPSD="true"; 
    } elsif ($field eq "-d")        { $delimiter=substr("$value ",1);
    } elsif ($field eq "-grid")     { $gridbased="true";
    } elsif ($field eq "-griddelta"){ $griddelta=$value+0; 
    } elsif ($field eq "-gridmax")  { $maxvoxels_grid=$value+0;
    } elsif ($field eq "-c++")      { $cpp="true";
    } elsif ($field eq "-nlin")     { if (-d "$nlopt_build") { $nlin="true"; } else { $nlin="false"; PRINT("[WARNING] the -nlin option is not available"); }; 
    } elsif ($field eq "-list")     { $use_shots_list="true"; $use_shots_file="$value"; 
    } elsif ($field eq "-kmpr")     { $known_max_pore_radius = $value+0; 
    } elsif ($field eq "-zeo")      { $zeo="true"; $np=1; 
    } elsif ($field eq "-fortran")  { $fortran="true"; 
    } elsif ($field eq "-dumpstep") { $dumpstep=$value+0; 
    } else                          { $ERROR="unknown argument $arg"; USAGE; };
};
$CONFIG = $config; 

if (-s "$config") { } else  { $ERROR="missing config file $config"; USAGE; };
if ($TPSD eq "true") { $fortran="true"; $nlin="false"; $cpp="false"; $gridbased="false"; }; 
if ($use_shots_list eq "true") { if (-s "$use_shots_file") { } else { $ERROR="missing shots file $use_shots_file"; USAGE; }; };
if (($fortran eq "true")&&(!$np_set)) { $np=$maxnp; }; 
if ($zeo eq "true") { if ($ro+$rs-$rp eq 0) { } else { $ERROR="zeo++ requires ro+rs=rp"; USAGE; }; };

# use temporary directory and copy files there, erase blanks to make voro++ work
$workdir=".tmp-GPSD-3D-$$"; 

# lammps-dump file
if (-s $config) {
   open(IN,"<$config"); $line=<IN>; $line=strip($line); 
   if ($line eq "ITEM: TIMESTEP") { 
      $firstdump = 1; 
      $dumpdone  = 0;
      if (!$setro) { close(IN); $ERROR="You have to specify a material sphere radius via -ro=..\n"; USAGE; }; 
      $newconfig = ".config.dump";
      open(NEW,">$newconfig");
      PRINT("[LAMMPS] handling lammps-dump file [$newconfig]. Overtaking box dimensions.\n");
      while ((!eof(IN))&&($dumpdone eq 0)) { 
         if ($firstdump eq 0) { $line=<IN>; }; 
         $timestep=<IN>+0; 
         $line=<IN>; 
         $N=<IN>+0; 
         $line=<IN>; 
         $line=<IN>; $line=strip($line); ($xlo,$xhi)=split(/ /,$line);
         $line=<IN>; $line=strip($line); ($ylo,$yhi)=split(/ /,$line);
         $line=<IN>; $line=strip($line); ($zlo,$zhi)=split(/ /,$line);
         $line=<IN>; $line=strip($line); @fields=split(/ /,$line);
         if (($dumpstep eq -1)||($dumpstep eq $timestep)) { 
            PRINT("[INFO] using dump timestep $timestep ($N atoms)\n");
            foreach $ifield (2 .. $#fields) {
               if      ($fields[$ifield] eq "x")  { $fx=$ifield-2;
               } elsif ($fields[$ifield] eq "y")  { $fy=$ifield-2;
               } elsif ($fields[$ifield] eq "z")  { $fz=$ifield-2;
               }; 
            };
            if (!$fx) { $ERROR="Your lammps dump-file has no x column"; USAGE; }; 
            if (!$fy) { $ERROR="Your lammps dump-file has no y column"; USAGE; };
            if (!$fz) { $ERROR="Your lammps dump-file has no z column"; USAGE; };
            foreach $id (1 .. $N) { 
               $line=<IN>; $line=strip($line); @fields=split(/ /,$line);
               print NEW "$fields[$fx] $fields[$fy] $fields[$fz]\n";
               $dumpdone=1;
            }; 
         } else {
            foreach $id (1 .. $N) { $line=<IN>; }; 
            $firstdump = 0; 
         }; 
      }; 
      close(NEW);
      $config=$newconfig;
      if ($dumpdone eq 0) { $ERROR="your dump file does not contain timestep $dumpstep!\n"; USAGE; }; 
      $using_lammps_dump = 1; 
      $box = 6; 
   };
   close(IN);
};

# set xlo, xhi, ... zhi
if ($boxfile) {  
    `rm -rf $workdir; mkdir $workdir`;
    PRINT("[PREPARING] scanning $boxfile\n");
    open(IN,"<$boxfile"); open(OUT,">$workdir/box.txt"); $line=<IN>; $line=strip($line); print OUT $line; close(IN); close(OUT);
    ($xlo,$xhi,$ylo,$yhi,$zlo,$zhi)=split(/ /,$line); 
} elsif ($box eq 6) {
    `rm -rf $workdir; mkdir $workdir`;
    PRINT("[PREPARING] creating boxfile\n");
    open(OUT,">$workdir/box.txt"); print OUT "$xlo $xhi $ylo $yhi $zlo $zhi"; close(OUT); 
} else { 
    open(IN,"<$config"); $line=<IN>; $line=strip($line); @c = split(/ /,$line); $cols=$#c+1; close(IN);
    if ($line=~/headers/) {
        `rm -rf $workdir; mkdir $workdir`;
        PRINT("[PREPARING] calling $converter2\n");
        `cp $config $workdir/config-samarth.csv`; 
        `cd $workdir; perl $converter2 config-samarth.csv 1.0`;
        if (-s "$workdir/config.txt") { } else { $ERROR="$converter2 failed in $workdir"; USAGE; }; 
        $config = "$workdir/config-samarth.txt";
        $boxfile = "$workdir/box-samarth.txt";
        `mv $workdir/box.txt $boxfile`; 
        `mv $workdir/config.txt $config`; 
        open(IN,"<$boxfile"); open(OUT,">$workdir/box.txt"); $line=<IN>; $line=strip($line); print OUT $line; close(IN); close(OUT);
        ($xlo,$xhi,$ylo,$yhi,$zlo,$zhi)=split(/ /,$line);
    } else { 
        $ERROR="box information missing. Specify file name via -box=... or use -xlo=.. -xhi=.. ... -zhi=.."; USAGE; 
    };
};

# create config*.txt files in $workdir, obtain N and set radii
if (-s $config) { 
    # get cols, create files with (grid-based) and without (voro-based) radius column
    open(IN,"<$config"); $line=<IN>; $line=strip($line); @c = split(/ /,$line); $cols=$#c+1; close(IN);
    if ($cols eq 3) {       # x y z
        PRINT("[PREPARING] recognized format (A)\n");
        $monodisperse = 1; $N = 0;
        open(IN,"<$config");
        open(OUT1,">$workdir/config.txt"); $N=0;
        open(OUT2,">$workdir/config-with-radii.txt"); 
        while (!eof(IN)) { $line=<IN>; $line=strip($line);
            if ($line) { $N+=1; print OUT1 "$N $line\n"; print OUT2 "$N $line $ro\n"; }; 
        };
        close(IN); close(OUT1); close(OUT2);  
    } elsif ($cols eq 4) { # id x y z
        PRINT("[PREPARING] recognized format (B)\n");
        $monodisperse = 1; $N = 0;
        open(IN,"<$config");
        open(OUT1,">$workdir/config.txt"); $N=0;
        open(OUT2,">$workdir/config-with-radii.txt");
        while (!eof(IN)) { $line=<IN>; $line=strip($line); @c=split(/ /,$line); 
            if ($line) { $N+=1; print OUT1 "$N $c[1] $c[2] $c[3]\n"; print OUT2 "$N $c[1] $c[2] $c[3] $ro\n"; }; 
        };
        close(IN); close(OUT1); close(OUT2);
    } elsif ($cols eq 5) {  # id x y z radius
        PRINT("[PREPARING] recognized format (C)\n");
        $monodisperse = 1; $N = 0;
        $radius[1] = $c[$#c];
        open(IN,"<$config"); 
        while (!eof(IN)) { $line=<IN>; $line=strip($line); @c = split(/ /,$line);
            if ($line) { $N+=1; $radius[$N]=$c[$#c]; if ($radius[$N] eq $radius[1]) { } else { $monodisperse=0; }; };            
        }; 
        close(IN);
        open(OUT2,">$workdir/config-with-radii.txt"); 
        open(IN,"<$config"); $N=0;
        if ($monodisperse eq 1) { 
            PRINT("[PREPARING] $config contains spheres with equal radius $radius[1] (monodisperse)\n");
            if ($setro eq 1) { } else { $ro=$radius[1]; }; 
            PRINT("[PREPARING] using ro=$ro because -ro was specified\n");
            open(OUT1,">$workdir/config.txt"); 
            while (!eof(IN)) { $line=<IN>; $line=strip($line); @c = split(/ /,$line);
                if ($line) { $N+=1; print OUT1 "$N $c[1] $c[2] $c[3]\n"; print OUT2 "$N $c[1] $c[2] $c[3] $ro\n"; };
            };
            close(OUT1); 
        } else {
            while (!eof(IN)) { $line=<IN>; $line=strip($line); @c = split(/ /,$line);
                if ($line) { $N+=1; print OUT2 "$c[0] $c[1] $c[2] $c[3] $c[4]\n"; };  
            };
        };
        close(IN); close(OUT2);
    } else {
        $ERROR="configuration file must have 3, 4 or 5 columns"; USAGE;
    }; 
    PRINT("[INFO] monodisperse: $monodisperse\n");
} else {
    $ERROR="configuration file not specified via -in=..."; USAGE;
};

PRINT("[INFO] $config contains $N particle coordinates ($cols columns)\n");
PRINT("[INFO] using ro=$ro, rc=$rc, rp=$rp\n");
PRINT("[INFO] created files in $workdir including .parameters.\n");

sub VORO_fortran {
# quality settings
$shots = int(0.5+$q*$N);
open(P,">$workdir/.parameters"); print P<<EOF;
\&list
N                                       = $N          ! number of particles
rp                                      = $rp         ! test particle radius    
ro                                      = $ro         ! particle radius
rc                                      = $rc         ! shell radius
shots                                   = $shots      ! former 1e6
more                                    = .$more.
quiet                                   = .$quiet.
np                                      = $np
TPSD                                    = .$TPSD.
use_shots_list                          = .$use_shots_list.
/
EOF
close(P);
PRINT("[GPSD-3D] Using $shots shots on $np threads\n");
PRINT("[GPSD-3D] Please stand by ..\n");
if (-s ".seed") { `cp .seed $workdir`; }; 
if ($quiet eq "true") {
    `cd $workdir; ($voroparser0; $voroparser) | $code-fortran.ex`;
} else {
    open my $cmd_mk, "cd $workdir; ($voroparser0; $voroparser) | $code-fortran.ex |";
    while (<$cmd_mk>) { print "$_"; };
};
if (-s "$workdir/.seed") { `cp $workdir/.seed ./`; }; 
if (-s "$workdir/_r")       { PRINT("[GPSD-3D] completed GPSD\n"); } else { $ERROR="GPSD-3D crashed. Remainings in $workdir"; showerr; };
if ($outputfile) { } else { $outputfile="$CONFIG-ro=$ro-rp=$rp-rc=$rc"; }; 
if (-s "$workdir/_r-TPSD")  { PRINT("[TPSD-3D] completed TPSD\n"); `mv $workdir/_r-TPSD $outputfile.tpsd`; $TPSDmessage="[TPSD-3D] created: $outputfile.tpsd\n"; };
`mv $workdir/_r $outputfile.gpsd`; 
if ($info eq 1) { `mv $workdir/_inf.txt $outputfile.inf`; }; 
if ($info eq 1) { `mv $workdir/_info.txt $outputfile.info`; };
if ($info eq 1) { `mv $workdir/.UpperPoreRadius ./`; };           # may help to speed up NONLIN
};

sub VORO_cpp {
# quality settings
$shots = int(0.5+$q*$N);
if ($more eq "true") { $add_cpp .= "-more"; }; 
if ($info eq 1)      { $add_cpp .= "-info"; }; 
PRINT("[GPSD-3D] Using $shots shots on $np threads\n");
PRINT("[GPSD-3D] Please stand by ..\n");
if (-s ".seed") { `cp .seed $workdir`; };
if ($quiet eq "true") {
    `cd $workdir; $code-cpp.ex -in=config.txt -ro=$ro -rp=$rp -rc=$rc -xlo=$xlo -xhi=$xhi -ylo=$ylo -yhi=$yhi -zlo=$zlo -zhi=$zhi -q=$q -d=\" \" $add_cpp`;
} else {
    open my $cmd_mk, "cd $workdir; $code-cpp.ex -in=config.txt -ro=$ro -rp=$rp -rc=$rc -xlo=$xlo -xhi=$xhi -ylo=$ylo -yhi=$yhi -zlo=$zlo -zhi=$zhi -q=$q -d=\" \" $add_cpp|";
    while (<$cmd_mk>) { print "$_"; };
};
if (-s "$workdir/.seed") { `cp $workdir/.seed ./`; };
if (-s "$workdir/results.gpsd") { PRINT("[GPSD-3D] completed GPSD\n"); } else { $ERROR="GPSD-3D crashed. Remainings in $workdir"; showerr; };
if ($outputfile) { } else { $outputfile="$CONFIG-ro=$ro-rp=$rp-rc=$rc"; };
`mv $workdir/results.gpsd $outputfile.gpsd`;
if ($info eq 1) { `mv $workdir/results-INFO.gpsd $outputfile.info`; };
if ($info eq 1) { `mv $workdir/results-INFO-numbers.gpsd $outputfile.inf`; };
};

sub NONLIN {
# create parameters-nlin for nonlinear optimization code
$shots = int(0.5+$q*$N);
PRINT("[GPSD-NLIN] Using $shots shots on $np threads\n");
open(P,">$workdir/parameters-nlin.txt"); print P<<EOF;
\&list
N                                       = $N          ! number of particles
rp                                      = $rp         ! test particle radius
ro                                      = $ro         ! particle radius
rc                                      = $rc
shots                                   = $shots
more                                    = .$more.
quiet                                   = .$quiet.
info                                    = .true.
use_radii_from_config                   = .true.
use_shots_list                          = .$use_shots_list.
known_max_pore_radius                   = $known_max_pore_radius
/
EOF
close(P);
PRINT("[GPSD-3D-NLIN] Executing $nlincode, please stand by ..\n");
open my $cmd_mk, "cd $workdir; export LD_LIBRARY_PATH=$nlopt_build; $nlincode |";
while (<$cmd_mk>) { print "$_"; };
if (-s "$workdir/radii.txt") { PRINT("[GPSD-3D-NLIN] completed\n"); } else { $ERROR="GPSD-3D-NLIN  crashed. Remainings in $workdir"; showerr; };
if ($outputfile) { } else { $outputfile="$CONFIG-ro=$ro-rp=$rp-rc=$rc"; };
`mv $workdir/radii.txt $outputfile.gpsd`;
if ($info eq 1) { `mv $workdir/radii.info $outputfile.info`; }; 
if ($info eq 1) { `mv $workdir/radii.inf  $outputfile.inf`; };
};

sub ZEO {
$shots = int(0.5+$q*$N);
PRINT("[GPSD-3D+ZEO+] Using $shots shots on $np threads\n");   
open(C,"<$workdir/config.txt"); 
open(Z,">$workdir/config.cssr"); 
$boxx=$xhi-$xlo; 
$boxy=$yhi-$ylo;
$boxz=$zhi-$zlo;
print Z "$boxx $boxy $boxz\n";
print Z "90   90   90   SPGR =  1 P 1         OPT = 1\n";
print Z "$N 0\n";
print Z "0 GPSD : GPSD\n";
foreach $id (1 .. $N) { $line=<C>; $line=strip($line); 
   ($i,$x,$y,$z)=split(/ /,$line); 
   $x/=$boxx;
   $y/=$boxy;
   $z/=$boxz;
   print Z "$id MK $x $y $z 0 0 0 0 0 0 0 0 0.00\n";  
};
close(C);
close(Z); 
`echo "MK $ro" > $workdir/use.rad`; 
`echo "MK 1.0" > $workdir/use.mass`; 
print "[ZEO++ INFO] Using use.rad:\n"; print `cat $workdir/use.rad`; 
print "[ZEO++ INFO] Using use.mass:\n"; print `cat $workdir/use.mass`; 
print "[ZEO++ INFO] Using config.cssr (first 7 lines, scaled coordinates):\n"; print `head -7 $workdir/config.cssr`;  
print "[ZEO++ INFO] Using zeo command: $zeocode -r use.rad -mass use.mass -psd $rp $rp $shots config.cssr\n";
PRINT("[GPSD-3D+ZEO+] Executing $zeocode, please stand by ..\n");
`cd $workdir; $zeocode -r use.rad -mass use.mass -psd $rp $rp $shots config.cssr 2> log.zeo++`; 
if (-s "$workdir/config.psd_histo") { PRINT("[GPSD-3D+ZEO+] completed\n"); } else { $ERROR="GPSD-3D+ZEO+ crashed. Remainings in $workdir"; showerr; };
if ($outputfile) { } else { $outputfile="$CONFIG-ro=$ro-rp=$rp-rc=$rc"; };
open(Z,"<$workdir/config.psd_histo"); 
open(R,">$outputfile.gpsd"); 
while (!eof(Z)) {
   $line=<Z>; $line=strip($line);
   if ($line =~ /Bin Count/) {
      $radii = 0;
      $mean_pore_radius = 0; 
      $mean_pore_radius_count = 0; 
      $max_pore_radius = 0;
      $min_pore_radius = 1e30;
      $stderr_pore_radius = 0;
      while (!eof(Z)) {
         $line=<Z>; $line=strip($line);
         ($diameter,$count,$cumulative,$derivative)=split(/ /,$line);
         if ($count>0) { 
            $radii += $count; 
            $pore_radius = $diameter/2.0; 
            if ($pore_radius<$min_pore_radius) { $min_pore_radius=$pore_radius; };
            if ($pore_radius>$max_pore_radius) { $max_pore_radius=$pore_radius; };
            foreach $mycount (1 .. $count) { print R "$pore_radius\n"; $mean_pore_radius+=$pore_radius; $mean_pore_radius_count+=1; $stderr_pore_radius+=$pore_radius**2; }; 
         }; 
      };
   };
};
close(R);
close(Z);
if ($mean_pore_radius_count>0) { 
   $mean_pore_radius/=$mean_pore_radius_count; 
   $stderr_pore_radius/=$mean_pore_radius_count;
   $stderr_pore_radius=sqrt($stderr_pore_radius-$mean_pore_radius**2)/sqrt($mean_pore_radius_count); 
}; 
PRINT("[GPSD-3D+ZEO+] min pore radius : $min_pore_radius\n");
PRINT("[GPSD-3D+ZEO+] max pore radius : $max_pore_radius\n");
PRINT("[GPSD-3D+ZEO+] mean pore radius: $mean_pore_radius +/- $stderr_pore_radius\n");
if ($info eq 1) { 
   $Volume = $boxx*$boxy*$boxz; 
   $phi_reff = 1.0-$radii/$shots;
   open(INFO,">$outputfile.info"); 
   print INFO<<EOF;
N=$N
ro=$ro
rp=$rp
rc=$rc
V=$Volume
triangles=-1
shots=$shots
radii=$radii
min_pore_radius=$min_pore_radius
max_pore_radius=$max_pore_radius
mean_pore_radius=$mean_pore_radius
stderr_pore_radius=$stderr_pore_radius
phi_reff=$phi_reff
cells=0
threads=$np           
voro_cpu_time=-1
voro_real_time=-1
triangles_setup_cpu_time=-1
triangles_setup_real-time=-1
MonteCarlo_cpu_time=-1 
MonteCarlo_real_time=-1
EOF
close(INFO);
open(INF,">$outputfile.inf");
   print INF<<EOF;
$N
$ro
$rp
$rc
$Volume
-1
$shots
$radii
$min_pore_radius
$max_pore_radius
$mean_pore_radius
$stderr_pore_radius
$phi_reff
0
$np
-1
-1
-1
-1
-1
-1
EOF
close(INF);
};
}; 

sub GRID {
# create parameters-grid for reference code
$min_delta_grid     = $griddelta*($ro+$rc);
PRINT("[GPSD-3D-GRID] using min_delta_grid = $min_delta_grid\n");
PRINT("[GPSD-3D-GRID] using maxvoxels = $maxvoxels_grid\n");
open(P,">$workdir/.parameters-grid"); print P<<EOF;
\&list
N                                       = $N          ! number of particles
rp                                      = $rp         ! test particle radius   
ro                                      = $ro         ! particle radius
rc                                      = $rc
min_delta                               = $min_delta_grid
maxvoxels                               = $maxvoxels_grid
use_radii_from_config                   = .true.
more                                    = .$more.
quiet                                   = .$quiet.
/
EOF
close(P);
PRINT("[GPSD-3D-GRID] Executing $refcode, please stand by ..\n");
if ($quiet eq "true") {
    `cd $workdir; $refcode`;
} else {
    open my $cmd_mk, "cd $workdir; $refcode |";
    while (<$cmd_mk>) { print "$_"; };
};
if (-s "$workdir/_radii") { PRINT("[GPSD-3D-GRID] completed\n"); } else { $ERROR="GPSD-3D-GRID crashed. Remainings in $workdir"; showerr; }; 
if ($outputfile) { } else { $outputfile="$CONFIG-ro=$ro-rp=$rp-rc=$rc"; }; 
`mv $workdir/_radii $outputfile.gpsd`; 
if ($info eq 1) { `mv $workdir/_info-grid.txt $outputfile.inf`; };
};

sub REMOVE_BLANKS {
    open(IN,"<$_[0]"); open(OUT,">$_[0].tmp"); 
    while (!eof(IN)) { $line=<IN>; $line=strip($line); print OUT "$line\n"; }; close(IN); close(OUT);
    `mv $_[0].tmp $_[0]`; 
};

if ($use_shots_list eq "true") { `cp $use_shots_file $workdir/p.txt`; }; 

if ($monodisperse) { 
    PRINT("[INFO] monodisperse system. The particle radius is taken as $ro, shell thickness $rc, test particle radius $rp.\n"); 
    if ($gridbased eq "true") { 
      PRINT("[GPSD-3D-GRID] using grid-based algorithm");
      GRID; 
    } elsif ($cpp eq "true") {  
      PRINT("[GPSD-3D] using voronoi (c++) algorithm");
      VORO_cpp;
    } elsif ($nlin eq "true") { 
      PRINT("[GPSD-3D-NLIN] using constrained optimization algorithm");
      NONLIN; 
    } elsif ($zeo  eq "true") {
      PRINT("[GPSD-3D+ZEO+] using zeo++");
      ZEO;
    } else {  
      PRINT("[GPSD-3D] using voronoi (fortran) algorithm");
      VORO_fortran; 
    }; 
} else {
    PRINT("[INFO] polydisperse system. The particle radius is overtaken from the configuration, shell thickness $rc, test particle radius $rp.\n"); 
    if ($gridbased eq "true") {
      GRID;
    } else { 
      NONLIN; 
    };  
};

# stop walltime clock
$walltime_toc=`date +"%s.%N"`+0; $walltime=$walltime_toc-$walltime_tic; 
if ($info eq 1) { `echo "$walltime">>$outputfile.inf`; };
if ($info eq 1) { `echo "walltime=$walltime">>$outputfile.info`; };

REMOVE_BLANKS("$outputfile.gpsd");
if ($TPSD eq "true") { REMOVE_BLANKS("$outputfile.tpsd"); }; 
if ($info eq 1)      { REMOVE_BLANKS("$outputfile.info"); PRINT("[GPSD-3D] created: $outputfile.info\n"); };
if ($info eq 1)      { REMOVE_BLANKS("$outputfile.inf");  PRINT("[GPSD-3D] created: $outputfile.inf\n"); };
`rm -rf $workdir`; 
PRINT("[GPSD-3D] created: $outputfile.gpsd\n"); 
PRINT("$TPSDmessage"); 
