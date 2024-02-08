#! /usr/bin/perl

if ($#ARGV eq 1) { } else { print "args: samarth-file usefraction\n"; exit; }; 

$file=$ARGV[0]; if (-s "$file") { } else { print "arg: samarth-file\n"; exit; }; 
$usefraction=$ARGV[1]; print "usefraction = $usefraction\n"; 

open(F,"<$file"); 
$line=<F>; chomp $line; ($dummy,$headers)=split(/=/,$line);
foreach $i (1 .. $headers) {
 $line=<F>; chomp $line;
 if ($line=~/^N=/)      { ($dummy,$N)=split(/=/,$line); }; 
 if ($line=~/^lattice/) { ($dummy,$lattice)=split(/=/,$line); };
 if ($line=~/^folded/)  { ($dummy,$folded)=split(/=/,$line); }; 
 if ($line=~/^phi=/)    { ($dummy,$phi)=split(/=/,$line); }; 
 if ($line=~/^D=/)      { ($dummy,$dim)=split(/=/,$line); };
 if ($line=~/^x0_lo=/)  { ($dummy,$lo[1])=split(/=/,$line); };
 if ($line=~/^x1_lo=/)  { ($dummy,$lo[2])=split(/=/,$line); };
 if ($line=~/^x2_lo=/)  { ($dummy,$lo[3])=split(/=/,$line); };
 if ($line=~/^x0_hi=/)  { ($dummy,$hi[1])=split(/=/,$line); };
 if ($line=~/^x1_hi=/)  { ($dummy,$hi[2])=split(/=/,$line); };
 if ($line=~/^x2_hi=/)  { ($dummy,$hi[3])=split(/=/,$line); };
};
$Voccupied=0; $pi=3.1415927;
open(L,">lattice.txt"); print L $lattice; close(L);
$min_radius = 1e30; $max_radius = -1e30; 
foreach $i (1 .. $N) { $line=<F>; chomp $line; 
 if ($dim eq 3) { 
  ($id,$x[$i],$y[$i],$z[$i],$assignedSeedStatus[$i],$currentSeedStatus[$i],$diameter,$attachments[$i],$rest)=split(/,/,$line); 
  if ($lattice eq 1) { $Voccupied += $diameter**3; };
  if ($lattice eq 0) { $Voccupied += $pi/6*$diameter**3; }; 
 } elsif ($dim eq 2) {
  ($id,$x[$i],$y[$i],$assignedSeedStatus[$i],$currentSeedStatus[$i],$diameter,$attachments[$i],$rest)=split(/,/,$line);
  $z[$i]=0; 
  if ($lattice eq 1) { $Voccupied += $diameter**2; };
  if ($lattice eq 0) { $Voccupied += $pi/4*$diameter**2; };
 }; 
 $radius[$i]=$diameter/2; 
 if ($radius[$i]<$min_radius) { $min_radius = $radius[$i]; }; 
 if ($radius[$i]>$max_radius) { $max_radius = $radius[$i]; }; 
}; 
if (eof(F)) { } else { print "file format error\n"; exit; }; 
close(F);

foreach $k (1 .. 3) { $box[$k]=$hi[$k]-$lo[$k]; }; print "OLD box sizes $box[1] $box[2] $box[3]\n";

$newbox[1] = $usefraction**(1.0/$dim)*$box[1];
$newbox[2] = $usefraction**(1.0/$dim)*$box[2];
$newbox[3] = $usefraction**(1.0/$dim)*$box[3];

$RAND = rand()*($box[1]-$newbox[1]); $lo[1]+=$RAND; $hi[1]=$lo[1]+$newbox[1];
$RAND = rand()*($box[2]-$newbox[2]); $lo[2]+=$RAND; $hi[2]=$lo[2]+$newbox[2];
$RAND = rand()*($box[3]-$newbox[3]); $lo[3]+=$RAND; $hi[3]=$lo[3]+$newbox[3];

print "$N particles in $dim dimensions (original, before usefraction)\n";
if ($dim eq 2) { $lo[3]=-0.001; $hi[3]=0.001; };
foreach $k (1 .. 3) { $box[$k]=$hi[$k]-$lo[$k]; }; print "NEW box sizes $box[1] $box[2] $box[3]\n";
$V=1; foreach $k (1 .. $dim) { $V*=$box[$k]; }; 
$myphi = $Voccupied/$V; 
print "samarth phi = $phi versus my phi = $myphi (FORGET IT)\n";

open(F,">box.txt"); print F "$lo[1] $hi[1] $lo[2] $hi[2] $lo[3] $hi[3]\n"; close(F);
open(F,">config.txt"); 
$newN=0;
if (abs(1.0-$min_radius/$max_radius)<0.001) { 
    $monodisperse = 1; 
    open(MONO,">.mono"); print MONO "$radius[1]"; close(MONO); 
    foreach $i (1 .. $N) { 
        if (($x[$i]>=$lo[1])&&($x[$i]<=$hi[1])&&($y[$i]>=$lo[2])&&($y[$i]<=$hi[2])&&($z[$i]>=$lo[3])&&($z[$i]<=$hi[3])) { 
            $newN+=1; print F "$newN $x[$i] $y[$i] $z[$i]\n"; 
        }; 
    }; 
} else {
    $monodisperse = 0; 
    `rm -f .mono`; 
    foreach $i (1 .. $N) {
        if (($x[$i]>=$lo[1])&&($x[$i]<=$hi[1])&&($y[$i]>=$lo[2])&&($y[$i]<=$hi[2])&&($z[$i]>=$lo[3])&&($z[$i]<=$hi[3])) {
            $newN+=1; print F "$newN $x[$i] $y[$i] $z[$i] $radius[$i]\n";
        };
    };
}; 
close(F);
print "$newN particles in $dim dimensions (using usefraction). Started from $N particles.\n";
$observedfraction = 100*$newN/$N; print "$observedfraction % realized\n";
print "monodisperse: $monodisperse (monodisperse saved without radius)\n"; 
print "created box.txt and config.txt and lattice.txt with $newN particles\n";

