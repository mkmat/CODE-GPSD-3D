# Note for windows-user: use / rather than \ in absolute path names

# section (A)
$exe[$#exe+1] = "g++";             $exe_download[$#exe] = "https://gcc.gnu.org/";
$exe[$#exe+1] = "gfortran";        $exe_download[$#exe] = "https://gcc.gnu.org/";
$exe[$#exe+1] = "voro++";          $exe_download[$#exe] = "https://math.lbl.gov/voro++/";

# section (B)
$lib[$#lib+1] = "libvoro++.a";     $lib_download[$#lib] = "https://math.lbl.gov/voro++/";
$lib[$#lib+1] = "libnlopt.so";     $lib_download[$#lib] = "https://pypi.org/project/nlopt/";

# section (C)
$build[$#build+1] = "automatic";
$build[$#build+1] = "automatic";

# section (D)
$file[$#file+1] = "voro++.hh";     # part of the voro++ distribution

# section (E) unused software
$unused[$#unused+1] = "network";   $unused_download[$#exe] = "https://github.com/richardjgowers/zeoplusplus";
