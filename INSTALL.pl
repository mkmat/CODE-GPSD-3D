#! /usr/bin/perl

# (c) 28 sep 2023 mk@mat.ethz.ch

$compiler           = "gfortran"; 
$compiler_options   = "-Ofast -ffree-line-length-none";
$code               = "GPSD-3D";
$refcode            = "GPSD-3D-grid";
$vorocode           = "voro++"; 
$converter1         = "convert-voro-output-to-triangles.pl";
$converter2         = "convert_samarth_to_config"; 
$install_dir        = `pwd`; chomp $install_dir; 

print `cat MESSAGE.txt`; 

# required system-wide executables and check if they exist
@exeutables=("voro++","$compiler");   
foreach $exe (@exeutables) { $which=`which $exe`; chomp $which; if (-s "$which") { print "found: $exe\n"; } else { print "missing executable $exe [which=$which]\n"; $m1=1; }; };

# required source codes and check if they exist
@sources=("$code.f90","$refcode.f90","$converter1","$converter2");
foreach $source (@sources) { if (-s "$source") { print "found: $source\n"; } else { print "missing source code $source\n"; $m2=1; }; };

# required local executables and check if they exist
@exeutables=("$code.ex","$refcode.ex","$converter1","$converter2"); 
foreach $exe (@exeutables) { if (-s "$exe") { `chmod 700 $exe`; } else { $m3=1; }; };

$M0     = "You have successfully installed GPSD-3D.\n";
$M0    .= "Copy GPSD-3D (just this single file) to a directory where you'll need it or where it can be found.\n"; 
$M0    .= "Then call (i) GPSD-3D or (ii) ./GPSD-3D or (iii) perl GPSD-3D to start GPSD-3D (includes documentation)\n"; 
$M1     = "If voro++ is missing: Download and install voro++ from https://math.lbl.gov/voro++/download/\n"; 
$M1    .= "If gfortran is missing: Download and install gfortran from https://fortran-lang.org/en/learn/os_setup/install_gfortran/";
$M2     = "You did not download the missing files which are distributed along with the GPSD-3D code\n";

if ($m1 eq 1) { print "$M1\n"; exit; };
if ($m2 eq 1) { print "$M2\n"; exit; };

if (($m3 eq 1)&&(!$m1)&&(!$m2)) {
    print "compiling ..\n";
    print `$compiler $compiler_options $code.f90 -o $code.ex`; 
    print `$compiler $compiler_options $refcode.f90 -o $refcode.ex`;
};

# required local executables and check if they exist
@exeutables=("$code.ex","$refcode.ex");
foreach $exe (@exeutables) { if (-s "$exe") { print "found: $exe\n"; } else { print "ERROR. $exe did not compile\n"; $m3=1; exit; }; };

open(S,">GPSD-3D"); 
open(T,"<TEMPLATE.pl"); while (!eof(T)) { $line=<T>; 
    $line=~s/InstallationDirectory/$install_dir/; 
    print S $line; 
}; 
close(S); close(T);
print "\n$M0\n"; 

print `cat NOTES.txt`; 
