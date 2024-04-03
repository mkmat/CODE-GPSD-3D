#! /usr/bin/perl

# mk@mat.ethz.ch

sub strip  { chomp $_[0]; $_[0]=~s/^\s+//g; $_[0]=~s/\s+$//; $_[0]=~s/\s+/ /g; $_[0]; }; 

$shots = 100000;
$ro=0.1; $rc=0; 

$target[0][1]  = 1.365; 
$target[0][2]  = 31.283;
$target[0][3]  = 1.37;
$target[0][6]  = 1.3137;
$target[0][7]  = 2.57;
$target[0][10] = 1.61;

$target[1][1]  = 1.329;
$target[1][2]  = 31.28;
$target[1][3]  = 1.37;
$target[1][6]  = 1.3137;
$target[1][7]  = 2.52;
$target[1][10] = 1.58;

@nos  = (1,2,3,6,10); if ($maxnp < 10) { @nos = (1,2); };
@irps = (0,1);          $noirps = ($#nos+1)*($#irps+1); 

print "[CHECKING] ",orange("now reproducing selected benchmark results ..\n");
$count=0;
foreach $no (@nos) {
foreach $irp (@irps) { $rp=$irp*0.1; 
   $count+=1;
   $N=`grep -c . benchmark/benchmark-$no-config`+0; 
   $q = $shots/$N; if ($q>1) { $q=int($shots/$N); }
   $command = "perl ./GPSD-3D -in=benchmark/benchmark-$no-config -box=benchmark/benchmark-$no-box -rp=$rp -ro=$ro -q=$q -o=tmp -info"; 
   print "[CHECKING $count/$noirps] $command; "; 
   `$command`; 
   $line=`grep "mean_pore_radius" tmp.info`; $line=strip($line); @tmp=split(/=/,$line); $mean=$tmp[1]; 
   $line=`grep "stderr_pore_radius" tmp.info`; $line=strip($line); @tmp=split(/=/,$line); $stderr=$tmp[1];
   # print "[benchmark-$no-config] <r> = $mean +/- $stderr\n"; 
   if (abs(1-$target[$irp][$no]/$mean-$target) < 0.01) {
      print green("PASSED")," <r> = $mean +/- $stderr\n";
   } else {
      print red("ERROR <r> = $mean +/- $stderr\n");
   }; 
   `rm -f tmp.gpsd tmp.inf tmp.info`; 
}; 
};

if ($cpp eq "true") {
 print "[CHECK++ ] ",orange("now reproducing selected benchmark results ..\n");
 $count=0;
 foreach $no (@nos) {
 foreach $irp (@irps) { $rp=$irp*0.1;
   $count+=1;
   $N=`grep -c . benchmark/benchmark-$no-config`+0;
   $q = $shots/$N; if ($q>1) { $q=int($shots/$N); }
   $command = "perl ./GPSD-3D -in=benchmark/benchmark-$no-config -box=benchmark/benchmark-$no-box -rp=$rp -ro=$ro -q=$q -o=tmp -info -c++";
   print "[CHECK++  $count/$noirps] $command; ";
   `$command`;
   $line=`grep "mean_pore_radius" tmp.info`; $line=strip($line); @tmp=split(/=/,$line); $mean=$tmp[1];
   $line=`grep "stderr_pore_radius" tmp.info`; $line=strip($line); @tmp=split(/=/,$line); $stderr=$tmp[1];
   # print "[benchmark-$no-config] <r> = $mean +/- $stderr\n";
   if (abs(1-$target[$irp][$no]/$mean-$target) < 0.01) {
      print green("PASSED")," <r> = $mean +/- $stderr\n";
   } else {
      print red("ERROR <r> = $mean +/- $stderr\n");
   };
   `rm -f tmp.gpsd tmp.inf tmp.info`;
 };
 };
};


1;
