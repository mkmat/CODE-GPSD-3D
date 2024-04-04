#! /usr/bin/perl

if ($#ARGV ge 4) { } else { print<<EOF;
This script creates replicated counterparts of the configuration + box-file
and saves them with a replicate-suffix.

use: perl $0 <configuration-file> <boxfile> replicate-x replicate-y replicate-z [noise]

Example:
perl $0 config.txt box.txt 2 1 1  # creates config.txt.211 and box.txt.211
EOF
exit;
};
if (-s "$ARGV[0]") { } else { print "missing $ARGV[0]\n"; exit; }; 
if (-s "$ARGV[1]") { } else { print "missing $ARGV[1]\n"; exit; };

open(B,"<$ARGV[1]"); $line=<B>; close(B); 
chomp $line; $line=~s/\s++/ /g; $line=~s/^ //; $line=~s/,/ /g; ($xlo,$xhi,$ylo,$yhi,$zlo,$zhi)=split(/ /,$line); 
$boxx=$xhi-$xlo;
$boxy=$yhi-$ylo;
$boxz=$zhi-$zlo;

$xhi+=($ARGV[2]-1)*$boxx;
$yhi+=($ARGV[3]-1)*$boxy;
$zhi+=($ARGV[4]-1)*$boxz;

open(C,"<$ARGV[0]"); @C=<C>; close(C);
$newfile = "$ARGV[0].$ARGV[2]$ARGV[3]$ARGV[4]"; 
open(C,">$newfile"); 

$id=0;
foreach $rx (1 .. $ARGV[2]) {
foreach $ry (1 .. $ARGV[3]) {
foreach $rz (1 .. $ARGV[4]) {
foreach $i (0 .. $#C) {
    $line=$C[$i]; chomp $line; $line=~s/\s++/ /g; $line=~s/^ //; $line=~s/,/ /g; 
    @tmp=split(/ /,$line);
    if ($#tmp eq 3) {
      ($origid,$x,$y,$z)=split(/ /,$line);
    } else {
      ($x,$y,$z)=split(/ /,$line);
    };
    $id+=1;
    $X=$x+($rx-1)*$boxx;
    $Y=$y+($ry-1)*$boxy;
    $Z=$z+($rz-1)*$boxz;
    if ($ARGV[5]) { 
      $X += $ARGV[5]*rand(); if ($X>$xhi) { $X-=($xhi-$xlo); }; 
      $Y += $ARGV[5]*rand(); if ($Y>$yhi) { $Y-=($yhi-$ylo); };
      $Z += $ARGV[5]*rand(); if ($Z>$zhi) { $Z-=($zhi-$zlo); };
    }; 
    print C "$id $X $Y $Z\n";
};
};
};
};
close(C);
print "created: $newfile\n";

$newfile = "$ARGV[1].$ARGV[2]$ARGV[3]$ARGV[4]";
open(B,">$newfile"); print B "$xlo $xhi $ylo $yhi $zlo $zhi"; close(B);
print "created: $newfile\n";

