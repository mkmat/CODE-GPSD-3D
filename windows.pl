# windows-user: Here is an example for a GPSD-3D call. You need to put your command into a file 
# like the present windows.pl, and then call the script from a windows (or vscode etc.) console via:
# perl ./windows.pl 

print `perl ./GPSD-3D -in=.benchmark-13-config -box=.benchmark-13-box -rp=0.0 -ro=1.0 -q=10`;
