#! /usr/bin/perl

# (c) 28 sep 2023 mk@mat.ethz.ch

# fortran compilers
$compiler[0]            = "ifort";
$compiler_options[0]    = "-O3 -fopenmp";

$compiler[1]            = "gfortran"; 
$compiler_options[1]    = "-Ofast -ffree-line-length-none -fopenmp";

# c++ compiler
$compiler_cpp           = "g++"; 

# codes
$code                   = "GPSD-3D";
$refcode                = "GPSD-3D-grid";
$threadcode             = "get-max-threads";
$voroparser             = "voro-parser";
$converter              = "convert_samarth_to_config.pl";

# voro++
$voro                   = "voro++"; 
$vorolib                = "libvoro++.a"; 

# default voro++ paths (if non-existent, INSTALL will search elsewhere)
$vorosrc_path           = "/usr/local/include/voro++/"; 
$vorolib_path           = "/usr/local/lib/libvoro++.a"; 

$install_dir            = `pwd`; chomp $install_dir; 

@pack                   = ("$code.f90","$refcode.f90","$converter","$voroparser.cpp","$threadcode.f90","MESSAGE.txt","LICENSE.txt","TEMPLATE.pl","INSTALL.pl");

print `cat MESSAGE.txt`; $i="[INSTALL]";

# check if default voro++ directories exist
if (-d "$vorosrc_path") { } else { $vorosrc_path=""; }; 
if (-s "$vorolib_path") { } else { $vorolib_path=""; }; 

# check if installation directory is complete
foreach $p (@pack) { if (-s "$p") { } else { print "incomplete GPSD-3D package, $p missing\n"; exit; }; };

# required compiler
foreach $j (0 .. $#compiler) { $exe=$compiler[$j]; $which=`which $exe`; chomp $which; if (-s "$which") { print "$i found: $which\n"; $use_compiler=$exe; $use_options=$compiler_options[$j]; }; };
if (!$use_compiler) { print "missing fortran90 compiler (gfortran, ifort, ..)"; $m1=1; };

# required system-wide executables and check if they exist or are found
@exeutables=("$voro","$compiler_cpp");   
foreach $exe (@exeutables) { $which=`which $exe`; chomp $which; if (-s "$which") { print "$i found: $which\n"; } else { print "missing executable $exe [which=$which]\n"; $m1=1; }; };

# find voro src
if (!$vorosrc_path) { 
    $vorosrc_path=`locate voro++.hh | head -1`; chomp $vorosrc_path; $vorosrc_path=~s/voro\+\+.hh$//; 
}; 
if (-d "$vorosrc_path") { print "$i found: $vorosrc_path\n"; `rm -f voro++; ln -s $vorosrc_path ./voro++`; } else { print "missing $vorosrc_path\n"; $m1=1; };

# required library
if (!$vorolib_path) { 
    $vorolib_path=`locate $vorolib | head -1`; chomp $vorolib_path; if (-s "$vorolib_path") { print "$i found: $vorolib_path\n";  } else { print "missing $vorolib_path\n"; $m1=1; };
}; 

# required source codes and check if they exist
@sources=("$code.f90","$refcode.f90","$voroparser.cpp","$converter","$threadcode.f90");
foreach $source (@sources) { if (-s "$source") { print "$i found: $source\n"; } else { print "missing source code $source\n"; $m2=1; }; };

# required local executables and check if they exist
@exeutables=("$code.ex","$refcode.ex","$voroparser.ex","$converter",".maxnp.ex"); 
foreach $exe (@exeutables) { if (-s "$exe") { `chmod 700 $exe`; } else { $m3=1; }; };
$m3 = 1;    # enforce compilation

$M0     = "You have successfully installed GPSD-3D.\n";
$M0    .= "Copy GPSD-3D (just this single file) to a directory where you'll need it or where it can be found.\n"; 
$M0    .= "Then call (i) GPSD-3D or (ii) ./GPSD-3D or (iii) perl GPSD-3D to start GPSD-3D (includes documentation)\n"; 
$M1     = "If $voro or $vorolib is missing: Download and install voro++ from https://math.lbl.gov/voro++/download/\n"; 
$M1    .= "If fortran90 compiler is missing: Download and install gfortran from https://fortran-lang.org/en/learn/os_setup/install_gfortran/ or choose another compiler in $0";
$M1    .= "If $compiler_cpp is missing: Download and install g++ or choose another compiler in $0";
$M2     = "You did not download the missing files which are distributed along with the GPSD-3D code\n";

if ($m1 eq 1) { print "$M1\n"; exit; };
if ($m2 eq 1) { print "$M2\n"; exit; };

if (($m3 eq 1)&&(!$m1)&&(!$m2)) {
    print "$i compiling $threadcode.f90 using $use_compiler ..\n"; print `$use_compiler $use_options $threadcode.f90 -o .maxnp.ex`; 
    print "$i compiling $code.f90 ..\n";       print `$use_compiler $use_options $code.f90 -o $code.ex`; 
    print "$i compiling $refcode.f90 ..\n";    print `$use_compiler $use_options $refcode.f90 -o $refcode.ex`;
    print "$i compiling $voroparser.cpp using $compiler_cpp ..\n"; print `$compiler_cpp -O3 -I $vorosrc_path $voroparser.cpp $vorolib_path -o $voroparser.ex`; 
};

# required local executables and check if they exist
@exeutables=("$code.ex","$refcode.ex","$voroparser.ex",".maxnp.ex");
foreach $exe (@exeutables) { if (-s "$exe") { print "$i found: $exe\n"; } else { print "ERROR. $exe did not compile\n"; $m3=1; exit; }; };

# retrieve maximum number of threads
$maxnp = `.maxnp.ex`+0; print "[INFO] available threads: $maxnp\n";

# pack, if this was chosen
if ($ARGV[0] eq "-pack") { print `tar -c -v -f GPSD-3D.tar @pack .bench*config .bench*box`; exit; }; 

open(S,">GPSD-3D"); 
open(T,"<TEMPLATE.pl"); while (!eof(T)) { $line=<T>; 
    $line=~s/InstallationDirectory/$install_dir/; 
    print S $line; 
}; 
close(S); close(T);
print "\n$M0\n"; 

print `cat NOTES.txt`; 


