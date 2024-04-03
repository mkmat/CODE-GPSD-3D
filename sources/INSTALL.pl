#! /usr/bin/perl

# (c) 13 nov 2023 mk@mat.ethz.ch

# codes
$code                   = "GPSD-3D";
$refcode                = "GPSD-3D-grid";
$nlincode               = "GPSD-3D-nlin";
$zeocode                = "ZEO++/zeoplusplus-master/network";
$threadcode             = "get-max-threads";
$voroparser0            = "voro-parser-0";
$voroparser             = "voro-parser";
$converter              = "convert_samarth_to_config.pl";
$replicater             = "replicate.pl";
$benchmark_tester       = "GPSD-3D-benchmark-tests.pl"; 
$hardcode               = "hard-code-software-location.pl";
$windows_script         = "windows.pl"; 

$install_dir            = `pwd`; chomp $install_dir;

$cpp                    = "false";   # developer
$nlin                   = "true";

sub green   { "\033[1;32;40m$_[0]\033[1;37;40m"; };
sub red     { "\033[1;31;40m$_[0]\033[1;37;40m"; };
sub orange  { "\033[1;33;40m$_[0]\033[1;37;40m"; }; 

$found                  = green("found");
$notfound               = red("not found");
$compiling              = orange("compiling");
$hline                  = "_______________________________________________________________________________________________________________________________\n";
if ($ARGV[0] eq "-static") { $static=1; };
$message_help = orange("Note: If you run into tremendous trouble with the installation, you may read the file software-locations-example.txt\n");
if      (`which locate`) { $locate="locate"; 
} elsif (`which mdfind`) { $locate="mdfind -name "; 
} else                   { $locate="locate";
}; 

# software locations
if (-s "software-locations.txt") { 
   require "./software-locations.txt";
   if (-s "$libnlopt") { $nlin="true"; } else { $nlin="false"; };  
} else {
   print `cat MESSAGE.txt`; 
   print "$hline\nGPSD-3D is looking for required software ..\n$hline\n";
   $hard = "./$hardcode"; 
   require "$hard"; 
   $message1 = "-----------------------------------------------\n"; $message2 = "Afterwards, start $0 again.\n$message1";

   # directories
   #$instdir[$#instdir+1]   = "nlopt-2.7.1";
   #$instdir[$#instdir+1]   = "voro++-0.4.6";

   foreach $ie (0 .. $#exe) { $e=$exe[$ie];
      $found_exe[$ie]=`which $e 2> /dev/null`; chomp $found_exe[$ie];
      if ($found_exe[$ie]) {
         print "$found: \$exe[$ie]\t= $e at $found_exe[$ie]\n";
      } else {
         $found_exe[$ie]=`find ./ -name "$e" -executable -print | head -1`; chomp $found_exe[$ie];
         if ($found_exe[$ie]) {
            $found_exe[$ie] = "$install_dir/$found_exe[$ie]"; 
            print "$found: \$exe[$ie]\t= $e at $found_exe[$ie]\n";
         } else {
            print "$message1";
            print "$notfound: \$exe[$ie]\t= $e\n";
            print "install $e from find $exe_download[$ie] and/or find $e on your disk and\nadd the full absolute directory name as \$exe[$ie] = .. to section (A) of $hard\n";
            print "$message2\n$message_help"; exit;
         };
      };
   };

   foreach $il (0 .. $#lib) { $l=$lib[$il];
      if (-s "$lib[$il]") { 
         $found_lib[$il] = $lib[$il]; 
      } else { 
         $found_lib[$il] = `$locate $l | head -1`; chomp $found_lib[$il];
         if (!$found_lib[$il]) { $found_lib[$il] = `find $install_dir -name $l -print | head -1`; chomp $found_lib[$il]; };
      }
      if ($found_lib[$il]) {
         print "$found: \$lib[$il]\t= $l at $found_lib[$il]\n";
      } else {
         print "$message1";
         print "$notfound: \$lib[$il]\t= $l\n";
         print "install $l from $lib_download[$il] and/or find $l on your disk and\nadd the full absolute directory name as \$lib[$il] = .. to section (B) of $hard\n";
         if ($l eq "libnlopt.so") {
            print "BUT: you have also the freedom to NOT install the constrained nonlinear optimization.\n";
            print "If so, GPSD-3D won't offer the -nlin option.\nWill you try to install $l? [y/n] ";
            $yesno = <>; chomp $yesno;
            if (($yesno eq "n")||($yesno eq "no")) { 
               $nlin = "false"; print "$message2";
            } else {
               print "$message2\n$message_help"; exit;
            };
         };        
      };
   };

   # locate voro++ installation directory
   foreach $il (0 .. $#lib) { $l=$lib[$il];
      if ($l =~ /libvoro\+\+.a/) {
         if (-d "$build[$il]") {
            print "$found: \$build[$il]\t= $build[$il]\n";
         } else {
            $build[$il] = $found_lib[$il]; $build[$il] =~ s/src\/libvoro\+\+.a//;
            if (-d "$build[$il]") {
               print "$found: \$build[$il]\t= $build[$il]\n";
            } else {
               $build[$il] = $found_lib[$il]; $build[$il] =~ s/lib\/libvoro\+\+.a//;
               if (-d "$build[$il]") {
                  print "$found: \$build[$il]\t= $build[$il]\n";
               } else {
                  print "$message1";
                  print "$notfound: \$build[$il]\t= $build[$il]\n";
                  print "Add the full absolute directory name of the voro++ installation directory as \$build[$il]=\"..\" to section (C) of $hard\n";
                  print "$message2\n$message_help"; exit;
               }; 
            };
         };
      };
   };

   # locate voro++.hh
   foreach $il (0 .. $#lib) { $l=$lib[$il];
      if ($l =~ /libvoro\+\+.a/) {
         if (-d "$file[$il]") {
            print "$found: \$file[$il] = $file[$il]\n";
         } else {
            if (-s "$build[$il]/src/$file[$il]") {
               $file[$il] = "$build[$il]/src/$file[$il]";
               print "$found: \$file[$il] = $file[$il]\n";
            } else {
               $file[$il] = `find $build[$il] -name $file[il]  -print | head -1`; chomp $file[$il];
               if (-s "$file[$il]") {
                  print "$found: \$file[$il] = $file[$il]\n";
               } else {
                  print "$message1";
                  print "$notfound: \$file[$il] = $file[$il]\n";
                  print "Add the full absolute file name of voro++.hh as \$file[$il] = \"..\" to section (D) of $hard\n";
                  print "$message2\n$message_help"; exit;
               };
            };
         };
      };
   };

   # locate nlopt build directory
   if ($nlin eq "true") { 
      foreach $il (0 .. $#lib) { $l=$lib[$il];
         if ($l =~ /libnlopt.so/) {
            if (-d "$build[$il]") {
               print "$found: \$build[$il]\t= $build[$il]\n";
            } else {
               $build[$il] = $found_lib[$il]; $build[$il] =~ s/\/libnlopt.so//;
               if (-d "$build[$il]") {
                  print "$found: \$build[$il]\t= $build[$il]\n";
               } else {
                  print "$message1";
                  print "$notfound: $build[$il]\t= $build[$il]\n";
                  print "Add the full absolute directory name of the NLopt installation directory as \$build[$il]=\"..\" to section (C) of $hard\n";
                  print "$message2\n$message_help"; exit;
               };
            };
         };
      };
   };

   # compiler options
   if ($static eq 1) {
      print "[INSTALL] static linking ..\n"; 
      if ($exe[0] eq "g++")      { $c_compiler_options = "-O3 -std=c++11 -I$install_dir -static"; };
      if ($exe[1] eq "ifort")    { $f_compiler_options = "-O3 -fopenmp -static_intel"; };
      if ($exe[1] eq "gfortran") { $f_compiler_options = "-Ofast -ffree-line-length-none -fopenmp -static"; };
   } else {
      if ($exe[0] eq "g++")      { $c_compiler_options = "-O3 -std=c++11 -I$install_dir "; };
      if ($exe[1] eq "ifort")    { $f_compiler_options = "-O3 -fopenmp"; };
      if ($exe[1] eq "gfortran") { $f_compiler_options = "-Ofast -ffree-line-length-none -fopenmp"; };
   }; 

   # now save the information
   $OUT=<<EOF;
# codes 
\$compiler_cpp = "$found_exe[0]";
\$compiler_f90 = "$found_exe[1]";
\$voro         = "$found_exe[2]";

# libraries
\$libvoro      = "$found_lib[0]";
\$libnlopt     = "$found_lib[1]";
\$build_voro   = "$build[0]";
\$build_nlopt  = "$build[1]";
\$vorohh       = "$file[0]";

# compiler options (please review)
\$f_compiler_options = "$f_compiler_options";
\$c_compiler_options = "$c_compiler_options";
EOF
   open(INFO,">software-locations.txt");
   print INFO $OUT;
   close(INFO);
   print green("\ncreated software-locations.txt"),". It contains:\n\n";
   print $OUT;
   `sleep 5`; 
};
require "./software-locations.txt";

open(VORO,">include-voro++.hh"); 
print VORO "#include \"$vorohh\""; 
close(VORO);

@pack     = ("$code.f90","$refcode.f90","$nlincode.f90","$converter","$replicater","$voroparser.cpp","$voroparser0.cpp","$threadcode.f90","MESSAGE.txt","LICENSE.txt","TEMPLATE.pl","INSTALL.pl");
@pack     = (@pack,"$benchmark_tester","$hardcode","$windows_script"); 
@pack     = (@pack,"software-locations-example.txt","README.txt"); 
if ($cpp eq "true") {
 @pack    = (@pack,"$code.cpp","simulation_box.hh","triangle.hh","voronoi_faces.hh","voronoi_particle.hh","coords.hh"); 
};
print `cat MESSAGE.txt`; $i="[INSTALL]";

# check if installation directory is complete
foreach $p (@pack) { if (-s "$p") { } else { print "incomplete GPSD-3D package, $p missing\n"; exit; }; };

# required source codes and check if they exist
@sources=("$code.f90","$refcode.f90","$nlincode.f90","$voroparser.cpp","$voroparser0.cpp","$converter","$replicater","$threadcode.f90");
foreach $source (@sources) { if (-s "$source") { print "$i $found: $source\n"; } else { print "$notfound: missing source code $source\n"; $m2=1; }; };

if ($cpp eq "true") {
   $source = "$code.cpp"; if (-s "$source") { print "$i $found: $source\n"; } else { print "$notfound: missing source code $source\n"; $m2=1; };
};

# required local executables and check if they exist
@exeutables=("$code-fortran.ex","$code-cpp.ex","$refcode.ex","$nlincode.ex","$voroparser.ex","$voroparser0.ex",".maxnp.ex"); 
foreach $exe (@exeutables) { if (-s "$exe") { `chmod 700 $exe`; } else { $m3=1; }; };

$m3 = 1;    # enforce compilation

$M0     = green("$hline\nYou have successfully installed GPSD-3D.\n");
$M0    .= "Copy GPSD-3D (just this single file) to a directory where you'll need it or where it can be found.\n"; 
$M0    .= "Then call (i) GPSD-3D or (ii) ./GPSD-3D or (iii) perl GPSD-3D to start GPSD-3D (includes documentation)\n"; 
$M0    .= "[For windows users: If your perl does not allow for command line arguments, see $windows_script]\n"; 
$M2     = "You did not download the missing files which are distributed along with the GPSD-3D code\n";

if ($m2 eq 1) { print "$M2\n"; exit; };



if (($m3 eq 1)&&(!$m2)) {
    # compile
    print "$i $compiling $compiler_f90 $f_compiler_options $threadcode.f90 -o .maxnp.ex ..\n"; 
     print `$compiler_f90 $f_compiler_options $threadcode.f90 -o .maxnp.ex`; 
    print "$i $compiling $compiler_f90 $f_compiler_options $code.f90 -o $code-fortran.ex ..\n";        
     print `$compiler_f90 $f_compiler_options $code.f90 -o $code-fortran.ex`; 
    print "$i $compiling $compiler_f90 $f_compiler_options $refcode.f90 -o $refcode.ex ..\n";     
     print `$compiler_f90 $f_compiler_options $refcode.f90 -o $refcode.ex`;
    print "$i $compiling $compiler_cpp -O3 -I $build_voro $voroparser.cpp $libvoro -o $voroparser.ex ..\n"; 
     print `$compiler_cpp -O3 -I $build_voro $voroparser.cpp $libvoro -o $voroparser.ex`; 
    print "$i $compiling $compiler_cpp -O3 -I $build_voro $voroparser0.cpp $libvoro -o $voroparser0.ex ..\n"; 
     print `$compiler_cpp -O3 -I $build_voro $voroparser0.cpp $libvoro -o $voroparser0.ex`;
    if ($cpp eq "true") { 
     print "$i $compiling $compiler_cpp $c_compiler_options -I ./ -I $build_voro $code.cpp $libvoro -o $code-cpp.ex ..\n";
     print `$compiler_cpp $c_compiler_options -I ./ -I $build_voro $code.cpp $libvoro -o $code-cpp.ex`;
    }; 
    if ($nlin eq "true") {
      print "$i $compiling $compiler_f90 $f_compiler_options -L $build_nlopt -lnlopt -lm $nlincode.f90 -o $nlincode.ex ..\n";    
      print `$compiler_f90 $f_compiler_options -L $build_nlopt -lnlopt -lm $nlincode.f90 -o $nlincode.ex`; 
    };
};

# required local executables and check if they exist
@exeutables=("$code-fortran.ex","$refcode.ex","$voroparser.ex","$voroparser0.ex",".maxnp.ex","$converter","$replicater","$benchmark_tester","$hardcode");
if ($nlin eq "true") { @exeutables=(@exeutables,"$nlincode.ex"); };
if ($cpp  eq "true") { @exeutables=(@exeutables,"$code-cpp.ex"); }; 
   
foreach $exe (@exeutables) { if (-s "$exe") { print "$i $found: $exe\n"; } else { print "ERROR. $exe did not compile\n"; $m3=1; exit; }; };

# retrieve maximum number of threads
$maxnp = `./.maxnp.ex`+0; print "[INFO]    ",green("available threads: $maxnp"),"\n";

# pack, if this was chosen
if ($ARGV[0] eq "-pack") { print `tar -c -v -f GPSD-3D.tar @pack benchmark`; exit; }; 

open(S,">GPSD-3D"); 
open(T,"<TEMPLATE.pl"); while (!eof(T)) { $line=<T>; 
    $line=~s/InstallationDirectory/$install_dir/; 
    print S $line; 
}; 
close(S); close(T);

# perform some benchmark tests
require "./$benchmark_tester";

print "\n$M0\n"; 
