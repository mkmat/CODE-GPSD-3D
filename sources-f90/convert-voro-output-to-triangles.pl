#! /usr/bin/perl

sub USAGE { print<<EOF;
$0 <config.vol>
EOF
exit;
}; 

if ($#ARGV eq -1) { USAGE; }; 
$configvol = $ARGV[0]; if (-s "$configvol") { } else { print "missing arg <config.vol>, missing $configvol\n"; exit; }; 

open(C,"<$configvol"); 
$total_triangles = 0; $total_faces = 0; 
open(TRI,">_triangles.txt");
while (!eof(C)) { 
    $line = <C>; chomp $line; $line=~s/\(//g; $line=~s/\)//g; 
    @tmp=split(/ /,$line); 
    $col=-1; 
    $col+=1; $id    = $tmp[$col]; 
    $col+=1; $x     = $tmp[$col]; 
    $col+=1; $y     = $tmp[$col];
    $col+=1; $z     = $tmp[$col]; 
    $col+=1; $faces = $tmp[$col]; 
    foreach $face (1 .. $faces) {
        $col+=1; ($nx[$face],$ny[$face],$nz[$face])=split(/,/,$tmp[$col]);  # face normal
        # print "n[$face] = $nx[$face] $ny[$face] $nz[$face]\n";
    }; 
    $col+=1; $vertices = $tmp[$col]; 
    foreach $vertex (0 .. $vertices-1) {
        $col+=1; ($rel_vertex_x,$rel_vertex_y,$rel_vertex_z)=split(/,/,$tmp[$col]);      # vertex coordinates
        $vertex_x[$vertex]=$rel_vertex_x+$x; 
        $vertex_y[$vertex]=$rel_vertex_y+$y;
        $vertex_z[$vertex]=$rel_vertex_z+$z;
    }; 
    foreach $face (1 .. $faces) {
        $col+=1; $vertexnos[$face]=$tmp[$col]; 
    }; 
    # print "col=$col from cols=$#tmp [$faces faces and $vertices vertices]\n";

    # create triangles 
    $total_faces += $faces; 
    foreach $face (1 .. $faces) {
        @iBs = split(/,/,$vertexnos[$face]); 
        if ($#iBs<2) { die 'whats that'; }; 
        # get one A (X-projected on face) for all triangles
        $Bx = $vertex_x[$iBs[0]];
        $By = $vertex_y[$iBs[0]];
        $Bz = $vertex_z[$iBs[0]];
        $Cx = $vertex_x[$iBs[1]];
        $Cy = $vertex_y[$iBs[1]];
        $Cz = $vertex_z[$iBs[1]];
        $d  = $Bx*$nx[$face]+$By*$ny[$face]+$Bz*$nz[$face];     # plane n.x = d
        $D  = $x*$nx[$face]+$y*$ny[$face]+$z*$nz[$face]-$d; 
        $Ax = $x-$D*$nx[$face]; 
        $Ay = $y-$D*$ny[$face];
        $Az = $z-$D*$nz[$face];
        $triangles = $#iBs+1; 
        $total_triangles += $triangles; 
        print TRI "$id $triangles\n";
        print TRI "$x $y $z\n";
        print TRI "$Ax $Ay $Az\n";
        foreach $j (0 .. $#iBs) {
            $j1 = ($j % ($#iBs+1));
            $j2 = (($j+1) % ($#iBs+1)); 
            $Bx = $vertex_x[$iBs[$j1]];
            $By = $vertex_y[$iBs[$j1]];
            $Bz = $vertex_z[$iBs[$j1]];
            $Cx = $vertex_x[$iBs[$j2]];
            $Cy = $vertex_y[$iBs[$j2]];
            $Cz = $vertex_z[$iBs[$j2]];      
            # print "[$face $j1 $j2 $#iBs]\n";
            #$normAB = sqrt(($Ax-$Bx)**2+($Ay-$By)**2+($Az-$Bz)**2);
            #$normAC = sqrt(($Ax-$Cx)**2+($Ay-$Cy)**2+($Az-$Cz)**2);
            #$normBC = sqrt(($Cx-$Bx)**2+($Cy-$By)**2+($Cz-$Bz)**2);
            #$Mx = ($normBC*$Ax+$normAC*$Bx+$normAB*$Cx)/($normBC+$normAC+$normAB); 
            #$My = ($normBC*$Ay+$normAC*$By+$normAB*$Cy)/($normBC+$normAC+$normAB);      # M incircle center
            #$Mz = ($normBC*$Az+$normAC*$Bz+$normAB*$Cz)/($normBC+$normAC+$normAB);
            #$normMA = sqrt(($Mx-$Ax)**2+($My-$Ay)**2+($Mz-$Az)**2);
            #$normMB = sqrt(($Mx-$Bx)**2+($My-$By)**2+($Mz-$Bz)**2);
            #$normMC = sqrt(($Mx-$Cx)**2+($My-$Cy)**2+($Mz-$Cz)**2);
            #$maxD = 0; 
            #if ($normMA>$normMB) { $maxD = $normMA; } else { $maxD = $normMB; };        # maxD max distance between incircle center and corner 
            #if ($normMC>$maxD)   { $maxD = $normMC; }; 
            #print TRI "$Ax $Ay $Az $Bx $By $Bz $Cx $Cy $Cz $x $y $z $Mx $My $Mz $maxD\n";
            print TRI "$Bx $By $Bz $Cx $Cy $Cz\n";
        }; 
    };  
    

}; 
close(C); 
close(TRI);

open(TRI,">.triangles-and-faces"); print TRI "$total_triangles $total_faces"; close(TRI); 
