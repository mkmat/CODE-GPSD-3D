#! /usr/bin/perl

# (c) 28 sep 2023 mk@mat.ethz.ch
# GPSD-3D had been generated from TEMPLATE.pl

$min_delta_grid     = 0.005;      # USER defines minimum grid spacing here (in units of effective particle radius ro+rc) 
$maxvoxels_grid     = 1000000;    # USER defines upper limit here

$installation_dir   = "InstallationDirectory";  # set by INSTALL.pl  

$code               = "$installation_dir/GPSD-3D.ex";
$refcode            = "$installation_dir/GPSD-3D-grid.ex";
$voroparser         = "$installation_dir/voro-parser.ex";
$converter1         = "$installation_dir/convert-voro-output-to-triangles.pl";
$converter2         = "$installation_dir/convert_samarth_to_config.pl";
$maxnp              = `$installation_dir/.maxnp.ex`+0;
$pwd                = `pwd`; chomp $pwd; 

# ------------------------------------------------
sub PRINT { if ($quiet eq "true") { } else { print "$_[0]"; }; };

# ------------------------------------------------
sub USAGE { print<<EOF;
perl $0 [-in=<configfile>] [-box=<boxfile>] [-rp=<rp>] [-ro=<ro>] [-rc=<rc>] [-q=<integer>] [-o=<outputfile>] [-more] [-info] [-quiet] [-np=<integer>]
or
perl $0 clean

where 

-in=..      : name of a configuration file (contains as many lines as material particles, three to five columns)
              (A) possible file format: x y z            (monodisperse system)
              (B) possible file format: id x y z         (monodisperse system)
              (C) possible file format: id x y z radius  (polydisperse system, monodisperse if all radii are equal)
              (D) samarth-type configuration file (do not specify -boxfile in this case)
-box=..     : name of a file that contains the box bounds: xlo xhi ylo yhi zlo zhi on a single line
              Alternatively, box bounds can be entered on the command line via -xlo=.. -xhi=.. -ylo=.. -yhi=.. -zlo=.. -zhi=..
-ro=..      : particle radius (required for formats A and B) (if set, ignores radius in file format C, all particles have radius ro then)
-rp=..      : probe particle radius (default: rp=0)
-rc=..      : particle coating thickness (only taken into account if file format C is used, ro=radius+rc is used)
-q=..       : A positive number (default: 10), determines the quality. The number of shots is q times the number of material spheres
-o=..       : name of the resulting file containing a list of pore radii (if not specified, resulting file ends with .gpsd)
-more       : if -more is added, not only the pore radius r (column 1) but also the center of the pore (columns 2-4) is added to outputfile.
-info       : if -info is added, runtime information is collected in a file whose name ends with -info (more information, see github repository)
-quiet      : if -quiet is added, GPSD-3D won't produce stdout.
-clean      : remove temporary directories .tmp-GPSD*, if existing (such directories may exist after a crash or interrupt).
-np=..      : maximum number of parallel processes to be used (between 1 and $maxnp)
-TPSD       : calculate T-PSDs in addition

Note: 

1) GPSD-3D can be called in parallel and sent to background. Each instance runs in a separate temporary directory and on a single core.
2) All temporary directories can be viewed via: ls -1 .tmp-GPSD*; 
3) Particle coordinates that do not fall into the simulation box are just ignored.
4) The coating thickness can be chosen to be negative, while ro+rc should be semipositive.

Resulting files: 

<configfile>-radii-GPSD-3D.txt          : contains a list of pore radii (monodisperse systems, voronoi-based)
<configfile>-radii-GPSD-3D-grid.txt     : contains a list of pore radii (polydisperse systems, grid-based)

Example mentioned at our github site: (using a benchmark configuration available in your installation directory)
perl ./GPSD-3D -in=.benchmark-7-config -box=.benchmark-7-box -rp=0.0 -ro=1.0 -q=10 -np=10

Examples using other existing benchmarks
EOF
$benchmarks = `ls -1 .bench*config | grep -c .`+0; 
foreach $no (1 .. $benchmarks) { 
    PRINT("perl $0 -in=.benchmark-$no-config -box=.benchmark-$no-box -rp=0.0 -ro=1.0 -q=10\n");
}; 
if ($ERROR) { 
    PRINT("________________________________________________\n\nERROR: $ERROR\n________________________________________________ cd $workdir\n"); 
}; 
exit;
};

# check for -quiet
if ("-quiet"~~@ARGV) { $quiet="true"; } else { $quiet="false";  print `cat MESSAGE.txt`; }; 

# init walltime clock
$walltime_tic=`date +"%s.%N"`+0; 

# scan arguments
$rc=0; $ro=0; $rp=0; $q=10; $box=0; $more="false"; $TPSD="false"; $np=$maxnp; 
foreach $arg (@ARGV) { $arg=~s/=/==/; ($field,$value)=split(/==/,$arg);
    if ($field eq "-in")            { $config="$value";   PRINT("[INFO] using configuration file $config\n"); 
    } elsif ($field eq "-box")      { $boxfile="$value";  PRINT("[INFO] using box file $boxfile\n"); 
    } elsif ($field eq "-rp")       { $rp=$value+0; 
    } elsif ($field eq "-ro")       { $ro=$value+0; $setro=1; 
    } elsif ($field eq "-rc")       { $rc=$value+0; 
    } elsif ($field eq "-q")        { $q=$value+0; if ($q < 1) { $q=1; }; 
    } elsif ($field eq "-clean")    { `rm -rf .tmp-GPSD-3D-*`; PRINT("[INFO] cleaned\n"); exit; 
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
    } elsif ($field eq "-np")       { $np=$value; if ($np<1) { $np=1; }; if ($np>$maxnp) { $np=$maxnp; }; 
    } elsif ($field eq "-TPSD")     { $TPSD="true"; 
    } else                          { $ERROR="unknown argument $arg"; USAGE; };
};
$CONFIG = $config; 

if (-s "$config") { } else { $ERROR="missing config file $config [$field:$value]"; USAGE; };

# use temporary directory and copy files there, erase blanks to make voro++ work
$workdir=".tmp-GPSD-3D-$$"; 

# set xlo, xhi, ... zhi
if ($boxfile) {  
    `rm -rf $workdir; mkdir $workdir`;
    PRINT("[PREPARING] scanning $boxfile\n");
    open(IN,"<$boxfile"); open(OUT,">$workdir/box.txt"); $line=<IN>; $line=~s/\s+/ /g; $line=~s/^ //; print OUT $line; close(IN); close(OUT);
    chomp $line; ($xlo,$xhi,$ylo,$yhi,$zlo,$zhi)=split(/ /,$line); 
} elsif ($box eq 6) {
    `rm -rf $workdir; mkdir $workdir`;
    PRINT("[PREPARING] creating boxfile\n");
    open(OUT,">$workdir/box.txt"); print OUT "$xlo $xhi $ylo $yhi $zlo $zhi"; close(OUT); 
} else { 
    open(IN,"<$config"); $line=<IN>; $line=~s/\s+/ /g; $line=~s/,/ /g; $line=~s/^ //; @columns = split(/ /,$line); $cols=$#columns+1; close(IN);
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
        open(IN,"<$boxfile"); open(OUT,">$workdir/box.txt"); $line=<IN>; $line=~s/\s+/ /g; $line=~s/^ //; print OUT $line; close(IN); close(OUT);
        chomp $line; ($xlo,$xhi,$ylo,$yhi,$zlo,$zhi)=split(/ /,$line);
    } else { 
        $ERROR="box information missing. Specify file name via -box=... or use -xlo=.. -xhi=.. ... -zhi=.."; USAGE; 
    };
};

# create config*.txt files in $workdir, obtain N and set radii
if (-s $config) { 
    # get cols, create files with (grid-based) and without (voro-based) radius column
    open(IN,"<$config"); $line=<IN>; $line=~s/\s+/ /g; $line=~s/,/ /g; $line=~s/^ //; @columns = split(/ /,$line); $cols=$#columns+1; close(IN);
    if ($cols eq 3) {       # x y z
        PRINT("[PREPARING] recognized format (A)\n");
        $monodisperse = 1; $N = 0;
        open(IN,"<$config");
        open(OUT1,">$workdir/config.txt"); $N=0;
        open(OUT2,">$workdir/config-with-radii.txt"); 
        while (!eof(IN)) { $line=<IN>; chomp $line; $line=~s/\s+/ /g; $line=~s/,/ /g; $line=~s/^ //; 
            if ($line) { $N+=1; print OUT1 "$N $line\n"; print OUT2 "$N $line $ro\n"; }; 
        };
        close(IN); close(OUT1); close(OUT2);  
    } elsif ($cols eq 4) { # id x y z
        PRINT("[PREPARING] recognized format (B)\n");
        $monodisperse = 1; $N = 0;
        open(IN,"<$config");
        open(OUT1,">$workdir/config.txt"); $N=0;
        open(OUT2,">$workdir/config-with-radii.txt");
        while (!eof(IN)) { $line=<IN>; chomp $line; $line=~s/\s+/ /g; $line=~s/,/ /g; $line=~s/^ //; @c=split(/ /,$line); 
            if ($line) { $N+=1; print OUT1 "$line\n"; print OUT2 "$c[0] $c[1] $c[2] $c[3] $ro\n"; }; 
        };
        close(IN); close(OUT1); close(OUT2);
    } elsif ($cols eq 5) {  # id x y z radius
        PRINT("[PREPARING] recognized format (C)\n");
        $monodisperse = 1; $N = 0;
        $radius[1] = $columns[$#columns];
        open(IN,"<$config"); 
        while (!eof(IN)) { $line=<IN>; chomp $line; $line=~s/\s+/ /g; $line=~s/,/ /g; $line=~s/^ //; @columns = split(/ /,$line);
            if ($line) { $N+=1; $radius[$N]=$columns[$#columns]; if ($radius[$N] eq $radius[1]) { } else { $monodisperse=0; }; };            
        }; 
        close(IN);
        open(OUT2,">$workdir/config-with-radii.txt"); 
        open(IN,"<$config"); $N=0;
        if ($monodisperse eq 1) { 
            PRINT("[PREPARING] $config contains spheres with equal radius $radius[1] (monodisperse)\n");
            if ($setro eq 1) { } else { $ro=$radius[1]; }; 
            PRINT("[PREPARING] using ro=$ro because -ro was specified\n");
            open(OUT1,">$workdir/config.txt"); 
            while (!eof(IN)) { $line=<IN>; chomp $line; $line=~s/\s+/ /g; $line=~s/,/ /g; $line=~s/^ //; @columns = split(/ /,$line);
                if ($line) { $N+=1; print OUT1 "$c[0] $c[1] $c[2] $c[3]\n"; print OUT2 "$c[0] $c[1] $c[2] $c[3] $ro\n"; };
            }; 
            close(OUT1); 
        } else {
            while (!eof(IN)) { $line=<IN>; chomp $line; $line=~s/\s+/ /g; $line=~s/,/ /g; $line=~s/^ //; @columns = split(/ /,$line);
                if ($line) { $N+=1; print OUT2 "$c[0] $c[1] $c[2] $c[3] $ro\n"; };
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
PRINT("[INFO] created files in $workdir including .parameters.\n");


sub VORO {
# quality settings
$shots = $q*$N;
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
/
EOF
close(P);
PRINT("[GPSD-3D] Using $shots shots on $np threads\n");
PRINT("[GPSD-3D] Please stand by ..\n");
if ($quiet eq "true") {
    `cd $workdir; ($voroparser; $voroparser) | $code`;
} else {
    open my $cmd_mk, "cd $workdir; ($voroparser; $voroparser) | $code |";
    while (<$cmd_mk>) { print "$_"; };
};
if (-s "$workdir/_r")       { PRINT("[GPSD-3D] completed GPSD\n"); } else { $ERROR="GPSD-3D crashed. Remainings in $workdir"; USAGE; };
if ($outputfile) { } else { $outputfile="$CONFIG-ro=$ro-rp=$rp-rc=$rc"; }; 
if (-s "$workdir/_r-TPSD")  { PRINT("[TPSD-3D] completed TPSD\n"); `mv $workdir/_r-TPSD $outputfile.tpsd`; $TPSDmessage="[TPSD-3D] created: $outputfile.tpsd\n"; };
`mv $workdir/_r $outputfile.gpsd`; 
if ($info eq 1) { `mv $workdir/_summary.txt $outputfile.info`; }; 
};

sub GRID {
# create parameters-grid for reference code
$min_delta_grid     = $min_delta_grid*($ro+$rc);
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
create_positions_plus_radii             = .false.
more                                    = .$more.
quiet                                   = .$quiet.
/
EOF
close(P);
PRINT("[GPSD-3D-GRID] Please stand by ..\n");
if ($quiet eq "true") {
    `cd $workdir; $refcode`;
} else {
    open my $cmd_mk, "cd $workdir; $refcode |";
    while (<$cmd_mk>) { print "$_"; };
};
if (-s "$workdir/_radii") { PRINT("[GPSD-3D-GRID] completed\n"); } else { $ERROR="GPSD-3D-GRID crashed. Remainings in $workdir"; USAGE; }; 
if ($outputfile) { } else { $outputfile="$CONFIG-ro=$ro-rp=$rp-rc=$rc"; }; 
`mv $workdir/_radii $outputfile.gpsd`; 
};

sub REMOVE_BLANKS {
    open(IN,"<$_[0]"); open(OUT,">$_[0].tmp"); 
    while (!eof(IN)) { $line=<IN>; $line=~s/\s+ / /g; $line=~s/^ //; print OUT $line; }; close(IN); close(OUT);
    `mv $_[0].tmp $_[0]`; 
};

if ($monodisperse) { 
    PRINT("[INFO] monodisperse system. The particle radius is taken as $ro, shell thickness $rc, test particle radius $rp.\n"); 
    VORO;
    # GRID; 
} else {
    PRINT("[INFO] polydisperse system. The particle radius is overtaken from the configuration, shell thickness $rc, test particle radius $rp.\n"); 
    GRID;
};

# stop walltime clock
$walltime_toc=`date +"%s.%N"`+0; $walltime=$walltime_toc-$walltime_tic; 
if ($info eq 1) { `echo "$walltime">>$outputfile.info`; };

REMOVE_BLANKS("$outputfile.gpsd");
if ($TPSD eq "true") { REMOVE_BLANKS("$outputfile.tpsd"); }; 
if ($info eq 1)      { REMOVE_BLANKS("$outputfile.info"); print "[GPSD-3D] created: $outputfile.info\n"; };
`rm -rf $workdir`; 
PRINT("[GPSD-3D] created: $outputfile.gpsd\n");
PRINT("$TPSDmessage"); 
