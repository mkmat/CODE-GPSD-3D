module shared
implicit none

    ! (c) mk@mat.ethz.ch 21 oct 2023
    ! Voronoi-based analytic approach
    ! monodisperse systems only
    ! GPSD by default, optionally TPSDs in addition

    integer         :: stdout,np,replicate(3)
    integer(kind=4) :: N,triangles,shots,count_r
    real*8          :: ro,rp,rc,ro_plus_rc,rs,rs2,rs24,lo(3),hi(3),box(3),lo_plus_halfbox(3),volume,volume_fraction
    real*8          :: UpperPoreRadius,min_pore_radius,max_pore_radius,mean_pore_radius,std_pore_radius
    integer         :: Tneighborlist_M(3),Xneighborlist_M(3)
    real*8          :: Tneighborlist_size(3),Xneighborlist_size(3),max_triangle_max_extension
    logical         :: quiet,more,TPSD
    logical         :: vertices_inside_material
    logical         :: use_shots_list

    real*8, dimension(:,:), allocatable  :: coord
    real*8, dimension(:,:), allocatable  :: triangle_o
    real*8, dimension(:,:), allocatable  :: triangle_n
    real*8, dimension(:,:), allocatable  :: triangle_A
    real*8, dimension(:,:), allocatable  :: triangle_B
    real*8, dimension(:,:), allocatable  :: triangle_C
    real*8, dimension(:,:), allocatable  :: triangle_vM
    real*8, dimension(:),   allocatable  :: triangle_Xn
    real*8, dimension(:),   allocatable  :: triangle_Rho
    real*8, dimension(:),   allocatable  :: triangle_chi

    integer, dimension(:,:,:), allocatable  :: Tneighborlist_first
    integer, dimension(:), allocatable      :: Tneighborlist_next
    integer, dimension(:,:,:), allocatable  :: Xneighborlist_first
    integer, dimension(:), allocatable      :: Xneighborlist_next
    integer :: Xi1min,Xi1max,Xi2min,Xi2max,Xi3min,Xi3max

    character(len=100), parameter        :: form100='("[GPSD-3D]",A45,3(I13,1x))'
    character(len=100), parameter        :: form101='("[GPSD-3D]",A45,3(F13.4,1x))'
    character(len=100), parameter        :: form102='("[GPSD-3D]",A45,F13.4,A3,F11.4)'
    character(len=100), parameter        :: form104='("[VORO++] ",A45,3(I13,1x))'
    character(len=100), parameter        :: form200='(7(F12.5,1x),I2)'
    character(len=100), parameter        :: form201='(F12.5)'
    character(len=100), parameter        :: form105='("[GPSD-3D]",A45,3(I3,1x))'
    character(len=100), parameter        :: form106='("[INFO]   ",A45)'

    real    :: cpu_and_real_time(5,2)

    namelist /list/ ro,rp,rc,shots,N,quiet,more,np,TPSD,use_shots_list

    contains 

    subroutine read_parameters
    logical :: EX
        open(2,file='.parameters'); rewind(2)
        read(2,nml=list)
        close(2)
        if (quiet) then
            stdout = 9
            open(stdout,file='log-GPSD-3D.txt'); rewind(stdout)
        else
            stdout = 6
        endif
        if (use_shots_list) then
            inquire(file='p.txt',exist=EX)   
            if (.not.EX) stop 'missing p.txt'
        endif
    return
    end

    subroutine read_box
        write(stdout,form100) 'reading box ..'
        open(2,file='box.txt'); read(2,*) lo(1),hi(1),lo(2),hi(2),lo(3),hi(3); close(2)
        box = hi-lo
        lo_plus_halfbox = lo+box/2.D0
        write(stdout,form101) 'box',box
    return
    end

    subroutine constants
        ro_plus_rc      = ro+rc
        rs              = ro+rc+rp
        rs2             = rs**2
        rs24            = 4.D0*rs2
    return
    end
    
    function distance_fold(vec)
    real*8  :: distance_fold(3),vec(3)
        distance_fold  = vec-anint(vec/box)*box
    return
    end

    function position_fold(vec)
    real*8  :: position_fold(3),vec(3)
         position_fold = vec - lo_plus_halfbox
         position_fold = position_fold - anint(position_fold/box)*box
         position_fold = position_fold + lo_plus_halfbox
    return
    end

    function norm(vec)
    real*8  :: norm,vec(3)
        norm = dsqrt(norm2(vec))
    return
    end

    function norm2(vec)
    real*8 :: norm2,vec(3)
        norm2 = dot_product(vec,vec)
    return
    end

    subroutine cross(vec1,vec2,vec3)
    real*8, intent(in)  :: vec1(3),vec2(3)
    real*8, intent(out) :: vec3(3)
        vec3(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
        vec3(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
        vec3(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)
    return
    end

    subroutine manage_cpu_and_real_clock(no)
    integer :: no,realcount,countrate
    real    :: realtime,cputime
        call system_clock(realcount,countrate)
        call cpu_time(cputime)
        realtime = realcount/dble(countrate)
        if (no.eq.1) then 
            cpu_and_real_time(no,:)=[cputime,realtime]
        else
            cpu_and_real_time(no,:)  =[cputime,realtime]
            cpu_and_real_time(no-1,:)=cpu_and_real_time(no,:)-cpu_and_real_time(no-1,:)
        endif
    return
    end

    subroutine init_seed
    logical EX
    integer oldseed(33)
        inquire(file='.seed',EXIST=EX)
        if (EX) then
            open(2,file='.seed'); rewind(2); read(2,*) oldseed; close(2)
            open(2,file='.oldseed'); rewind(2); write(2,*) oldseed; close(2)
            call random_seed(put=oldseed)
        else
            call random_seed()
        endif
    return
    end

    subroutine finalize_seed
    integer oldseed(33)
        call random_seed(get=oldseed)
        open(2,file='.seed'); rewind(2); write(2,*) oldseed; close(2)
    return
    end

    subroutine read_voro_parser
    ! parser
    integer :: voro_max_vertices,voro_max_faces,voro_max_vertices_for_face
    integer :: voro_vertices,voro_faces
    integer :: itmp(2),j2
    integer, dimension(:), allocatable   :: voro_vertices_for_face
    integer, dimension(:,:), allocatable :: voro_vid
    real*8, dimension(:,:), allocatable  :: voro_vertex
    ! triangles
    integer(kind=4)     :: tri,face,i,j,id
    real*8  :: A(3),B(3),C(3),E(3),X(3),facenormal(3)
    real*8  :: CB(3),CE(3),XA(3)
    real*8  :: normGA,normGB,normGC
    real*8  :: normAX,normBX,normCX
    real*8  :: triangle_max_extension
    ! parser
    real*8  :: vA(3),vB(3),vC(3),vE(3),o(3),ref(3)
    real*8  :: GeomCenter(3)

    ! read once to determine field sizes
    voro_max_vertices           = 0
    voro_max_faces              = 0
    voro_max_vertices_for_face  = 0
    triangles                   = 0

    ! attention: overwriting N as seen by voro++
    read(*,*) N,voro_max_vertices,voro_max_faces,voro_max_vertices_for_face,triangles
    allocate(coord(N,3))

    write(stdout,form104) 'N',N
    write(stdout,form104) 'max vertices',voro_max_vertices
    write(stdout,form104) 'max faces',voro_max_faces
    write(stdout,form104) 'max face-vertices',voro_max_vertices_for_face
    write(stdout,form100) 'triangles',triangles

    allocate(voro_vertices_for_face(voro_max_vertices))
    allocate(voro_vid(voro_max_faces,voro_max_vertices_for_face))
    allocate(voro_vertex(voro_max_vertices,3))

    allocate(triangle_o(triangles,3))   ! o
    allocate(triangle_n(triangles,3))   ! n
    allocate(triangle_A(triangles,3))   ! A  = vA-o
    allocate(triangle_B(triangles,3))   ! B  = vB-vA
    allocate(triangle_C(triangles,3))   ! C  = vC-vA
    allocate(triangle_vM(triangles,3))  ! triangle center
    allocate(triangle_Xn(triangles))    ! X.n with X=x-o
    allocate(triangle_Rho(triangles))   ! rho
    allocate(triangle_chi(triangles))   ! chi

    tri = 0
    UpperPoreRadius = 0.D0
    max_triangle_max_extension = 0.D0
    vertices_inside_material = .false.

    read(*,*) itmp(1)        ! attention: overwriting N as seen by voro++
    itmp(2) = 0
    if (N.ne.itmp(1)) stop 'parser error (called twice?)'
    do id=1,N
         ! I am assuming absolute, folded, vertex coordinates here
        read(*,*,end=111) coord(id,:),voro_vertices,voro_faces,(voro_vertices_for_face(face),face=1,voro_faces)
        read(*,*) (voro_vertex(j,:),j=1,voro_vertices),((voro_vid(face,j),j=1,voro_vertices_for_face(face)),face=1,voro_faces)   
        itmp(2) = itmp(2)+1
        do face=1,voro_faces
            ! determine geometrical face center at vA
            ref = voro_vertex(voro_vid(face,1),:)       ! voro_vertex in global coordinate system
            vA  = 0.D0
            do j=1,voro_vertices_for_face(face)
                vE = distance_fold(voro_vertex(voro_vid(face,j),:)-ref)
                vA = vA + vE 
            enddo
            vA  = position_fold(ref + vA/voro_vertices_for_face(face))
            ! determine face normal, origin, and A = vA-o
            vB  = voro_vertex(voro_vid(face,1),:) 
            vC  = voro_vertex(voro_vid(face,2),:)   ! requires vertices_for_face > 2
            vE  = voro_vertex(voro_vid(face,3),:)
            CB  = distance_fold(vB-vC)
            CE  = distance_fold(vE-vC)
            call cross(CB/norm(CB),CE/norm(CE),facenormal) 
            facenormal = facenormal/norm(facenormal)                
            XA  = distance_fold(coord(id,:)-vA)
            o   = coord(id,:)-dot_product(XA,facenormal)*facenormal   ! origin
            o   = position_fold(o)
            A   = distance_fold(vA-o)      
            X   = distance_fold(coord(id,:)-o)
            facenormal = X/norm(X)                                  ! the sign is used to ensure Pn>0 lateron
            do j=1,voro_vertices_for_face(face)   ! (0 .. $#iBs) {
                j2  = j+1; if (j2.eq.1+voro_vertices_for_face(face)) j2=1     ! sorted via voro++
                vB  = voro_vertex(voro_vid(face,j) ,:)
                vC  = voro_vertex(voro_vid(face,j2),:)
                B   = distance_fold(vB-vA)
                C   = distance_fold(vC-vA)
                tri = tri + 1
                triangle_o(tri,:)  = o
                triangle_A(tri,:)  = A
                triangle_B(tri,:)  = B
                triangle_C(tri,:)  = C
                triangle_n(tri,:)  = facenormal 
                triangle_Xn(tri)   = dot_product(X,facenormal)      ! Xn > 0
                ! geometric center of triangle
                GeomCenter = (B+C)/3.D0     ! relative to vertex A
                triangle_vM(tri,:) = position_fold(vA + GeomCenter)
                ! max distance vertex - geometric center  
                normGA = norm(GeomCenter)
                normGB = norm(GeomCenter-B)
                normGC = norm(GeomCenter-C)
                triangle_max_extension = max(normGA,max(normGB,normGC))
                triangle_chi(tri) = triangle_max_extension
                max_triangle_max_extension = max(max_triangle_max_extension,triangle_max_extension)
                normAX  = norm(A-X)
                normBX  = norm(A+B-X)
                normCX  = norm(A+C-X)
                ! old: rho = |X-vA|-ro-rc
                if (min(normAX,min(normBX,normCX)).lt.ro+rc) vertices_inside_material = .true.
                triangle_Rho(tri)   = max(normAX,max(normBX,normCX))-rs ! (largest R=r+rp we may produce using this triangle)                      
                UpperPoreRadius     = max(UpperPoreRadius,triangle_Rho(tri)+rp)  ! r=R+rp
            enddo
        enddo
    enddo
111 continue
    if (itmp(2).lt.N-100) stop 'voro++ inconsistency'
    N= itmp(2) ! taking care of seldom voro++ inconsistencies
    write(stdout,form100) 'parallel processes (np)',np
    write(stdout,form100) 'material spheres (N)',N
    write(stdout,form101) 'UpperPoreRadius',UpperPoreRadius     ! this value is correct
    write(stdout,form101) 'ro',ro
    write(stdout,form101) 'rc',rc
    write(stdout,form101) 'rp',rp
    write(stdout,form101) 'reff = rc+rp',rc+rp
    write(stdout,form101) 'rs = ro+reff',rs
    write(stdout,form101) 'max triangle extension',max_triangle_max_extension
    if (vertices_inside_material) write(stdout,form106) 'note: vertices inside material'
    if (.not.vertices_inside_material) write(stdout,form106) 'note: all vertices in void space'

    replicate = 1+UpperPoreRadius/(box/2.D0) 
    if (max(replicate(1),max(replicate(2),replicate(3))).gt.1) then
        write(stdout,form100) 'ERROR: modify your configuration first using'
        write(stdout,form105) '... perl replicate <yourconfig> <yourboxfile>',replicate
        write(stdout,form100) '... Then call GPSD-3D again'
        write(stdout,form100) '... using the new files'
        stop
    endif

    if (triangles.ne.tri) stop 'format error'
    
    deallocate(voro_vertices_for_face)
    deallocate(voro_vid)
    deallocate(voro_vertex)
    ! if (.not.TPSD.and..not.more) deallocate(coord)        ! taken out since Xneighbor

    return
    end

    subroutine setup_particles_neighborlist
    integer(kind=4)     :: i
    integer     :: icell(3),ncells
        write(stdout,form100) 'creating X-neighbor list ..'

        Xneighborlist_size   = max(ro+rc,box/100.D0)
        Xneighborlist_M      = max(1,int(box/Xneighborlist_size))
        Xneighborlist_size   = box/Xneighborlist_M
        ncells               = Xneighborlist_M(1)*Xneighborlist_M(2)*Xneighborlist_M(3)

        write(stdout,form100) 'Xneighborlist_M',Xneighborlist_M
        write(stdout,form101) 'Xneighborlist_size',Xneighborlist_size
        write(stdout,form101) 'particles per cell',N/dble(ncells)

        allocate(Xneighborlist_first(Xneighborlist_M(1),Xneighborlist_M(2),Xneighborlist_M(3)))
        allocate(Xneighborlist_next(N))
        Xneighborlist_first = 0
        Xneighborlist_next  = 0
        do i=1,N
            icell = 1+int((coord(i,:)-lo)/Xneighborlist_size)
            icell = max(1,min(Xneighborlist_M,icell))
            Xneighborlist_next(i) = Xneighborlist_first(icell(1),icell(2),icell(3))
            Xneighborlist_first(icell(1),icell(2),icell(3)) = i
        enddo 
        ! 
        Xi1min = -1; Xi1max = 1
        Xi2min = -1; Xi2max = 1
        Xi3min = -1; Xi3max = 1
        if (Xneighborlist_M(1).eq.1) Xi1min=1
        if (Xneighborlist_M(1).eq.2) Xi1min=0
        if (Xneighborlist_M(2).eq.1) Xi2min=1
        if (Xneighborlist_M(2).eq.2) Xi2min=0
        if (Xneighborlist_M(3).eq.1) Xi3min=1
        if (Xneighborlist_M(3).eq.2) Xi3min=0
    return
    end

    function inside_coated_material(pos)
    logical :: inside_coated_material
    real*8  :: pos(3)
    integer     :: icell(3),bead,jcell(3),i1,i2,i3
         icell = 1+int((pos-lo)/Xneighborlist_size)
         icell = max(1,min(Xneighborlist_M,icell)) 
         inside_coated_material = .true.
         do i1=Xi1min,Xi1max
            jcell(1) = icell(1)+i1
            if (jcell(1).eq.Xneighborlist_M(1)+1) jcell(1)=1
            if (jcell(1).eq.0) jcell(1)=Xneighborlist_M(1)
            do i2=Xi2min,Xi2max
               jcell(2) = icell(2)+i2
               if (jcell(2).eq.Xneighborlist_M(2)+1) jcell(2)=1
               if (jcell(2).eq.0) jcell(2)=Xneighborlist_M(2)
               do i3=Xi3min,Xi3max
                  jcell(3) = icell(3)+i3
                  if (jcell(3).eq.Xneighborlist_M(3)+1) jcell(3)=1
                  if (jcell(3).eq.0) jcell(3)=Xneighborlist_M(3)
                  bead = Xneighborlist_first(jcell(1),jcell(2),jcell(3))
                  do while (bead.ne.0)  
                     if (norm(distance_fold(coord(bead,:)-pos)).le.ro+rc) return
                     bead = Xneighborlist_next(bead)
                  enddo
               enddo
            enddo
         enddo
         inside_coated_material = .false.
    return
    end

    subroutine setup_triangles_neighborlist
    integer(kind=4)     :: i
    integer     :: icell(3),ncells
        write(stdout,form100) 'creating T-neighbor list ..'

        Tneighborlist_size   = UpperPoreRadius-rp+max_triangle_max_extension
        Tneighborlist_M      = max(1,int(box/Tneighborlist_size))
        Tneighborlist_size   = box/Tneighborlist_M
        ncells               = Tneighborlist_M(1)*Tneighborlist_M(2)*Tneighborlist_M(3)

        write(stdout,form100) 'Tneighborlist_M',Tneighborlist_M
        write(stdout,form101) 'Tneighborlist_size',Tneighborlist_size
        write(stdout,form101) 'triangles per cell',triangles/dble(ncells)

        allocate(Tneighborlist_first(Tneighborlist_M(1),Tneighborlist_M(2),Tneighborlist_M(3)))
        allocate(Tneighborlist_next(triangles))
        Tneighborlist_first = 0
        Tneighborlist_next  = 0
        do i=1,triangles
            icell = 1+int((triangle_vM(i,:)-lo)/Tneighborlist_size)   
            icell = max(1,min(Tneighborlist_M,icell))
            Tneighborlist_next(i) = Tneighborlist_first(icell(1),icell(2),icell(3))
            Tneighborlist_first(icell(1),icell(2),icell(3)) = i
        enddo
        ! deallocate(triangle_vM)
    return
    end

    subroutine MonteCarlo
    real    :: myrand3(3)
    real*8  :: R
    real*8  :: r_is_above
    real*8  :: vp(3),P(3),Pn,v_center(3),PM(3)
    real*8  :: o(3),A(3),B(3),C(3),X(3),Xn
    real*8  :: facenormal(3)
    integer :: icell(3),jcell(3),k,i1,i2,i3
    real*8  :: tmp_r,tmp_Center(3)
    integer(kind=4) :: shot,i
    integer :: i1min,i2min,i3min
    integer :: i1max,i2max,i3max

    integer :: ruled_out,ruled_out_A,ruled_out_B,ruled_out_C,ruled_out_D

        write(stdout,form100) 'starting Monte Carlo ..'

        min_pore_radius        = 100.0*UpperPoreRadius
        max_pore_radius        = -1.0
        mean_pore_radius       = 0.D0
        std_pore_radius        = 0.D0
        count_r                = 0
        volume                 = box(1)*box(2)*box(3)
        write(stdout,form101) 'volume V=V(0,-ro)',volume

        ruled_out   = 0
        ruled_out_A = 0
        ruled_out_B = 0
        ruled_out_C = 0
        ruled_out_D = 0

        call init_seed

        i1min = -1; i1max = 1
        i2min = -1; i2max = 1
        i3min = -1; i3max = 1
        if (Tneighborlist_M(1).eq.1) i1min=1
        if (Tneighborlist_M(1).eq.2) i1min=0
        if (Tneighborlist_M(2).eq.1) i2min=1
        if (Tneighborlist_M(2).eq.2) i2min=0
        if (Tneighborlist_M(3).eq.1) i3min=1
        if (Tneighborlist_M(3).eq.2) i3min=0

        open(3,file='_r'); rewind(3)
!$omp parallel num_threads ( np) shared ( i1min, i1max, i2min, i2max ) private ( shot, myrand3, vp, P, Pn, o, icell, jcell, k, i, A, facenormal, B, C, tmp_r, tmp_Center, r, v_center, X, Xn, PM, r_is_above )

!$omp do reduction (min : min_pore_radius) reduction (max : max_pore_radius) reduction( + : mean_pore_radius, count_r, std_pore_radius)
        do shot=1,shots
          call random_number(myrand3)
          vp    = lo + myrand3*box
          if (.not.inside_coated_material(vp)) then 
            R = 0.D0
            icell = 1+int((vp-lo)/Tneighborlist_size)
            do i1=i1min,i1max
             jcell(1) = icell(1)+i1
             if (jcell(1).eq.Tneighborlist_M(1)+1) jcell(1)=1
             if (jcell(1).eq.0) jcell(1)=Tneighborlist_M(1)
            do i2=i2min,i2max
             jcell(2) = icell(2)+i2
             if (jcell(2).eq.Tneighborlist_M(2)+1) jcell(2)=1
             if (jcell(2).eq.0) jcell(2)=Tneighborlist_M(2)
            do i3=i3min,i3max
             jcell(3) = icell(3)+i3
             if (jcell(3).eq.Tneighborlist_M(3)+1) jcell(3)=1
             if (jcell(3).eq.0) jcell(3)=Tneighborlist_M(3)
                i = Tneighborlist_first(jcell(1),jcell(2),jcell(3))
                do while (i.ne.0) 
                    if (R.lt.triangle_Rho(i)) then ! otherwise, we cannot improve R using this triangle
                        facenormal  = triangle_n(i,:)
                        o           = triangle_o(i,:)
                        P           = distance_fold(vp-o)
                        Pn          = dot_product(P,facenormal)
                        if (Pn.lt.0.D0) then
                            facenormal = -facenormal
                            Pn         = -Pn
                        endif
                        r_is_above  = Pn+rp  ! Pn+rp 
                        if (r_is_above.le.UpperPoreRadius) then     ! otherwise we won't find any R (ok)
                            PM = distance_fold(vp-triangle_vM(i,:))
                            r_is_above = norm(PM)+rp-triangle_chi(i)
                            if (r_is_above.le.UpperPoreRadius) then ! otherwise we won't find any R (ok)
                                Xn = triangle_Xn(i)
                                X  = Xn*facenormal
                                if (norm2(P-X).ge.rs2) then
                                    A  = triangle_A(i,:)
                                    B  = triangle_B(i,:)
                                    C  = triangle_C(i,:)
                                    call candidate_for_triangle(facenormal,A,B,C,X,Xn,P, tmp_R,tmp_Center)
                                    if (tmp_R.gt.R) then    
                                        R = tmp_R                   ! this is r-rp
                                        v_center = o+tmp_Center     ! actually only needed when more=.true.
                                    endif
                                endif
                            endif
                        endif
                    endif         
                    i = Tneighborlist_next(i)
                enddo
            enddo
            enddo
            enddo

            if (R.gt.0.D0) then
                r = R+rp
                count_r                 = count_r + 1
                if (more) then
                    write(3,form200) vp,position_fold(v_center),r
                else
                    write(3,form201) r
                endif
                min_pore_radius         = min(min_pore_radius,r)
                max_pore_radius         = max(max_pore_radius,r)
                mean_pore_radius        = mean_pore_radius + r
                std_pore_radius         = std_pore_radius  + r**2
            endif
          endif ! not inside material
        enddo ! shot
!$omp end do
!$omp end parallel
        close(3)

        volume_fraction  = 1.D0-dble(count_r)/dble(shots)
        mean_pore_radius = mean_pore_radius / dble(max(1,count_r))
        std_pore_radius  = std_pore_radius / dble(max(1,count_r))
        std_pore_radius  = (std_pore_radius-mean_pore_radius**2)/dsqrt(dble(max(1,count_r)))


        write(stdout,form101) 'volume fraction phi(reff)',volume_fraction
        write(stdout,form101) 'V(0|reff)',(1-volume_fraction)*volume
        write(stdout,form101) 'min pore radius',min_pore_radius
        write(stdout,form102) 'mean pore radius',mean_pore_radius,' +/-',std_pore_radius
        write(stdout,form101) 'max pore radius',max_pore_radius
        write(stdout,form100) 'pore radius {r} list size',count_r
    return
    end

    subroutine MonteCarlo_using_shots_list
    real    :: myrand3(3)
    real*8  :: R
    real*8  :: r_is_above
    real*8  :: vp(3),P(3),Pn,v_center(3),PM(3)
    real*8  :: o(3),A(3),B(3),C(3),X(3),Xn
    real*8  :: facenormal(3)
    integer :: icell(3),jcell(3),k,i1,i2,i3
    real*8  :: tmp_r,tmp_Center(3)
    integer(kind=4) :: shot,i,shots_target
    integer :: i1min,i2min,i3min
    integer :: i1max,i2max,i3max
    integer :: id

    integer :: ruled_out,ruled_out_A,ruled_out_B,ruled_out_C,ruled_out_D

        write(stdout,form100) 'starting Monte Carlo ..'

        min_pore_radius        = 100.0*UpperPoreRadius
        max_pore_radius        = -1.0
        mean_pore_radius       = 0.D0
        std_pore_radius        = 0.D0
        count_r                = 0
        volume                 = box(1)*box(2)*box(3)
        write(stdout,form101) 'volume V=V(0,-ro)',volume

        ruled_out   = 0
        ruled_out_A = 0
        ruled_out_B = 0
        ruled_out_C = 0
        ruled_out_D = 0

        call init_seed

        i1min = -1; i1max = 1
        i2min = -1; i2max = 1
        i3min = -1; i3max = 1
        if (Tneighborlist_M(1).eq.1) i1min=1
        if (Tneighborlist_M(1).eq.2) i1min=0
        if (Tneighborlist_M(2).eq.1) i2min=1
        if (Tneighborlist_M(2).eq.2) i2min=0
        if (Tneighborlist_M(3).eq.1) i3min=1
        if (Tneighborlist_M(3).eq.2) i3min=0

        open(20,file='p.txt'); rewind(20)

        open(3,file='_r'); rewind(3)
        shots_target = shots
        shots = 0
        do while (shots.le.shots_target) 
            R     = 0.D0
            read(20,*,end=1) id,vp     
            shots = shots + 1
            if (.not.inside_coated_material(vp)) then
               icell = 1+int((vp-lo)/Tneighborlist_size)
               do i1=i1min,i1max
                  jcell(1) = icell(1)+i1
                  if (jcell(1).eq.Tneighborlist_M(1)+1) jcell(1)=1
                  if (jcell(1).eq.0) jcell(1)=Tneighborlist_M(1)
                  do i2=i2min,i2max
                     jcell(2) = icell(2)+i2
                     if (jcell(2).eq.Tneighborlist_M(2)+1) jcell(2)=1
                     if (jcell(2).eq.0) jcell(2)=Tneighborlist_M(2)
                     do i3=i3min,i3max
                        jcell(3) = icell(3)+i3
                        if (jcell(3).eq.Tneighborlist_M(3)+1) jcell(3)=1
                        if (jcell(3).eq.0) jcell(3)=Tneighborlist_M(3)
                        i = Tneighborlist_first(jcell(1),jcell(2),jcell(3))
                        do while (i.ne.0)
                           if (R.lt.triangle_Rho(i)) then ! otherwise, we cannot improve R using this triangle
                              facenormal  = triangle_n(i,:)
                              o           = triangle_o(i,:)
                              P           = distance_fold(vp-o)
                              Pn          = dot_product(P,facenormal)
                              if (Pn.lt.0.D0) then
                                 facenormal = -facenormal
                                 Pn         = -Pn
                              endif
                              r_is_above  = Pn+rp  ! Pn+rp
                              if (r_is_above.le.UpperPoreRadius) then     ! otherwise we won't find any R (ok)
                                 PM = distance_fold(vp-triangle_vM(i,:))
                                 r_is_above = norm(PM)+rp-triangle_chi(i)
                                 if (r_is_above.le.UpperPoreRadius) then ! otherwise we won't find any R (ok)
                                    Xn = triangle_Xn(i)
                                    X  = Xn*facenormal
                                    if (norm2(P-X).ge.rs2) then
                                       A  = triangle_A(i,:)
                                       B  = triangle_B(i,:)
                                       C  = triangle_C(i,:)
                                       call candidate_for_triangle(facenormal,A,B,C,X,Xn,P, tmp_R,tmp_Center)
                                       if (tmp_R.gt.R) then
                                          R = tmp_R                   ! this is r-rp
                                          v_center = o+tmp_Center     ! actually only needed when more=.true.
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                           i = Tneighborlist_next(i)
                        enddo
                     enddo
                  enddo
               enddo

               if (R.gt.0.D0) then
                  r = R+rp
                  count_r                 = count_r + 1
                  if (more) then
                    write(3,form200) vp,position_fold(v_center),r
                  else
                    write(3,form201) r
                  endif
                  min_pore_radius         = min(min_pore_radius,r)
                  max_pore_radius         = max(max_pore_radius,r)
                  mean_pore_radius        = mean_pore_radius + r
                  std_pore_radius         = std_pore_radius  + r**2
               endif
            endif
        enddo ! while
1       continue
        close(20)
        close(3)

        volume_fraction  = 1.D0-dble(count_r)/dble(shots)
        mean_pore_radius = mean_pore_radius / dble(max(1,count_r))
        std_pore_radius  = std_pore_radius / dble(max(1,count_r))
        std_pore_radius  = (std_pore_radius-mean_pore_radius**2)/dsqrt(dble(max(1,count_r)))

        write(stdout,form100) 'shots into void space',count_r
        write(stdout,form101) 'volume fraction phi(reff)',volume_fraction
        write(stdout,form101) 'V(0|reff)',(1-volume_fraction)*volume
        write(stdout,form101) 'min pore radius',min_pore_radius
        write(stdout,form102) 'mean pore radius',mean_pore_radius,' +/-',std_pore_radius
        write(stdout,form101) 'max pore radius',max_pore_radius
        write(stdout,form101) 'created a list {r} of pore radii'
    return
    end


    subroutine candidate_for_triangle(facenormal,A,B,C,X,Xn,P, R,Center)
        real*8, parameter       :: mkeps    = 1.D-10
        real*8, intent(in)      :: facenormal(3),A(3),B(3),C(3),X(3),Xn,P(3)
        real*8, intent(out)     :: R,Center(3)
        real*8  :: ex(3),ey(3),ez(3)
        real*8  :: P2,test_Center(3)
        real*8  :: Qsquared,Q1squared,Q2squared
        real*8  :: Ap,Aq,Bp,Bq,Cp,Cq,Pn,Pp,Xn2
        real*8  :: t(2),u(2)
        logical :: useful(2)

        ! introduce Q1 = R1+rs and Q2 = R2+rs and Q1squared = Q1^2

        Qsquared    = rs**2
        Center      = 0.D0

        ! base system ex,ey,ez
        ex  = facenormal                    ! ex = n
        ey  = P-dot_product(P,ex)*ex        ! ey = n_p
        ey  = ey/norm(ey)
        call cross(ex,ey,ez)                ! ez = n_q

        Bp   = dot_product(B,ey)
        Bq   = dot_product(B,ez)
        Cp   = dot_product(C,ey)
        Cq   = dot_product(C,ez)
        Pn   = dot_product(P,ex)          ! py=0
        Pp   = dot_product(P,ey)
        Ap   = dot_product(A,ey)
        Aq   = dot_product(A,ez)
        P2   = Pn**2+Pp**2
        Xn2  = Xn**2

        ! type 1 (corner B)
        test_Center = A+B
        ! Q1squared = norm2(X-test_Center)                Xn**2
        Q1squared = Xn2+(Aq+Bq)**2+(Ap+Bp)**2
        if (Q1squared.ge.Qsquared) then
            Q2squared = (norm(P-test_Center)+rs)**2
            if (Q1squared.ge.Q2squared-mkeps) then
                Qsquared = Q1squared
                Center   = test_Center
            endif
        endif
        ! type 2 (corner C)
        test_Center  = A+C
        ! Q1squared = norm2(X-test_Center)
        Q1squared = Xn2+(Aq+Cq)**2+(Ap+Cp)**2
        if (Q1squared.ge.Qsquared) then
            Q2squared = (norm(P-test_Center)+rs)**2
            if (Q1squared.ge.Q2squared-mkeps) then
                Qsquared = Q1squared
                Center   = test_Center
            endif
        endif

        ! case iii (edge B-C) (u=1-t) note R1=R2
        call make_case_iii(Ap,Aq,Bp,Bq,Cp,Cq,Xn2,Pp,P2, t,useful)
        if (useful(1)) then
            test_Center = A+C+t(1)*(B-C)        ! A + t*B + u*C
            Q1squared = norm2(X-test_Center)
            if (Q1squared.ge.Qsquared) then
                Qsquared = Q1squared
                Center   = test_Center          ! type 3
            endif
        endif
        if (useful(2)) then
            test_Center = A+C+t(2)*(B-C)        ! A + t*B + u*C
            Q1squared = norm2(X-test_Center)
            if (Q1squared.ge.Qsquared) then
                Qsquared = Q1squared
                Center   = test_Center          ! type 4
            endif
        endif

        ! case iv (at Sq=0)  note R1=R2
        if (dabs(Cq).gt.0.D0) then
            call make_case_iv(Ap,Aq,Bp,Bq,Cp,Cq,Xn2,Pp,P2, t,u,useful)
            if (useful(1)) then
                test_Center = A+t(1)*B+u(1)*C
                Q1squared = norm2(X-test_Center)
                if (Q1squared.ge.Qsquared) then
                    Qsquared = Q1squared
                    Center   = test_Center      ! type 5
                endif
            endif
            if (useful(2)) then
                test_Center = A+t(2)*B+u(2)*C
                Q1squared = norm2(X-test_Center)
                if (Q1squared.ge.Qsquared) then
                    Qsquared = Q1squared
                    Center   = test_Center      ! type 6
                endif
            endif
        endif

        R = dsqrt(Qsquared)-rs

    return
    end

    subroutine make_case_iii(Ap,Aq,Bp,Bq,Cp,Cq,Xn2,Pp,P2, t,useful)
    real*8  :: Ups,t(2)
    real*8  :: Ap,Aq,Bp,Bq,Cp,Cq,Xn2,Pp,P2
    real*8  :: a,b,c,d,e
    real*8  :: b2e,ApCp,AqCq,BpCp,BqCq,sqrtUps,dab
    logical :: useful(2)
        ApCp = Ap+Cp
        BpCp = Bp-Cp
        a = 2.D0*ApCp*Pp+rs2+Xn2-P2
        b = 2.D0*BpCp*Pp
        if (a.lt.0.D0.and.a+b.lt.0.D0) then     
            useful = .false.; return
        endif
        AqCq = Aq+Cq
        BqCq = Bq-Cq
        c = rs24*(Xn2+ApCp**2+AqCq**2)        ! c = 4*rs2*(X2+dot(A+C,A+C));
        d = rs24*(BpCp*ApCp+BqCq*AqCq)          ! d = 4*rs2*dot(B-C,A+C);
        e = rs24*(BpCp**2+BqCq**2)              ! e = 4*rs2*dot(B-C,B-C);
        b2e = b**2-e
        Ups = a**2*e-2.D0*a*b*d+d**2+b2e*c;
        if (Ups.ge.0.D0) then
            sqrtUps = dsqrt(Ups)
            dab     = d-a*b
            t(1) = (dab+sqrtUps)/b2e
            t(2) = (dab-sqrtUps)/b2e
            useful = .true.
            if (a+b*t(1).lt.0.D0.or.t(1).lt.0.D0.or.t(1).gt.1.D0) useful(1)=.false.
            if (a+b*t(2).lt.0.D0.or.t(2).lt.0.D0.or.t(2).gt.1.D0) useful(2)=.false.
        else
            useful = .false.
        endif
    return
    end

    subroutine make_case_iv(Ap,Aq,Bp,Bq,Cp,Cq,Xn2,Pp,P2, t,u,useful)
    real*8  :: Ups,t(2),u(2)
    real*8  :: Ap,Aq,Bp,Bq,Cp,Cq,Xn2,Pp,P2
    real*8  :: a,b,c,d,e,f,h
    real*8  :: b2e,AphCp,BpfCp,BqfCq,sqrtUps,dab
    logical :: useful(2)
        h = Aq/Cq
        f = Bq/Cq 
        AphCp = Ap-h*Cp
        BpfCp = Bp-f*Cp
        a = 2.D0*AphCp*Pp+rs2+Xn2-P2
        b = 2.D0*BpfCp*Pp
        if (a.lt.0.D0.and.a+b.lt.0.D0) then     ! 30 oct 2023
            useful = .false.; return
        endif
        BqfCq = Bq-f*Cq
        c = rs24*(Xn2+AphCp**2)
        d = rs24*AphCp*BpfCp   
        e = rs24*( BpfCp**2 + BqfCq**2 ) 
        b2e = b**2-e
        Ups = a**2*e-2.D0*a*b*d+d**2+b2e*c 
        if (Ups.ge.0.D0) then
            sqrtUps = dsqrt(Ups)
            dab     = d-a*b
            t(1) = (dab+sqrtUps)/b2e
            t(2) = (dab-sqrtUps)/b2e
            u    = -h-f*t
            useful = .true.
            if (a+b*t(1).lt.0.D0.or.t(1).lt.0.D0.or.t(1).gt.1.D0.or.u(1).lt.0.D0.or.t(1)+u(1).gt.1.D0) useful(1)=.false.
            if (a+b*t(2).lt.0.D0.or.t(2).lt.0.D0.or.t(2).gt.1.D0.or.u(2).lt.0.D0.or.t(2)+u(2).gt.1.D0) useful(2)=.false.
        else
            useful = .false.
        endif
    return
    end

    subroutine calc_TPSD
    integer     :: M(3),icell(3),jcell(3),i,k,i1,i2,i3,bead
    real*8      :: nsize(3),pos(3),mydist,r
    real        :: myrand3(3)
    integer     :: i1min,i1max,i2min,i2max,i3min,i3max
    integer(kind=4) :: shot
    integer, dimension(:,:,:), allocatable  :: firstbead
    integer, dimension(:), allocatable      :: nextbead
        deallocate(triangle_o)
        deallocate(triangle_n)
        deallocate(triangle_A)
        deallocate(triangle_B)
        deallocate(triangle_C)
        deallocate(triangle_vM)
        deallocate(triangle_Xn)
        deallocate(triangle_Rho)
        deallocate(triangle_chi)
        deallocate(Tneighborlist_first)
        deallocate(Tneighborlist_next)

        write(stdout,form100) 'creating _r-TPSD  ..'
        M     = max(1,int(box/(UpperPoreRadius+ro+rc)))  ! because we need V(0|-ro) as well
        nsize = box/M
        write(stdout,form101) '[TPSD-3D] cutoff',nsize
        write(stdout,form100) '[TPSD-3D] cells',M
        allocate(firstbead(M(1),M(2),M(3)))
        allocate(nextbead(N))
        firstbead = 0
        nextbead = 0
        do i=1,N
            icell = 1 + int((coord(i,:)-lo)/nsize)
            nextbead(i) = firstbead(icell(1),icell(2),icell(3))
            firstbead(icell(1),icell(2),icell(3)) = i
        enddo
        i1min = -1; i1max = 1
        i2min = -1; i2max = 1
        i3min = -1; i3max = 1
        if (nsize(1).eq.1) i1min=1
        if (nsize(1).eq.2) i1min=0
        if (nsize(2).eq.1) i2min=1
        if (nsize(2).eq.2) i2min=0
        if (nsize(3).eq.1) i3min=1
        if (nsize(3).eq.2) i3min=0
        open(3,file='_r-TPSD'); rewind(3)
        do shot=1,shots
            r = UpperPoreRadius+ro+rc
            call random_number(myrand3)
            pos = lo+box*myrand3
            icell = 1 + int((pos-lo)/nsize)
            do i1=i1min,i1max
                jcell(1) = icell(1)+i1
                if (jcell(1).eq.M(1)+1) jcell(1)=1 
                if (jcell(1).eq.0)      jcell(1)=M(1)
            do i2=i2min,i2max
                jcell(2) = icell(2)+i2
                if (jcell(2).eq.M(2)+1) jcell(2)=1
                if (jcell(2).eq.0)      jcell(2)=M(2)
            do i3=i3min,i3max
                jcell(3) = icell(3)+i3
                if (jcell(3).eq.M(3)+1) jcell(3)=1
                if (jcell(3).eq.0)      jcell(3)=M(3)
                bead = firstbead(jcell(1),jcell(2),jcell(3))
                do while (bead.ne.0)
                    mydist = norm(distance_fold(pos-coord(bead,:)))
                    if (mydist.lt.r) r=mydist      
                    bead = nextbead(bead)
                enddo
            enddo   
            enddo
            enddo
            write(3,*) r         ! smallest distance to material particle center
        enddo
        close(3)
    return
    end

end module shared

program GPSD_3D
use omp_lib
use shared
    if (np.gt.OMP_GET_MAX_THREADS()) stop 'np > OMP_GET_MAX_THREADS'

    call manage_cpu_and_real_clock(1)
    call read_parameters
    call constants
    call read_box
    call init_seed
    call manage_cpu_and_real_clock(2); call read_voro_parser
    call manage_cpu_and_real_clock(3); call setup_triangles_neighborlist; call setup_particles_neighborlist
    call manage_cpu_and_real_clock(4); if (use_shots_list) then; call MonteCarlo_using_shots_list; else; call MonteCarlo; endif
    call manage_cpu_and_real_clock(5)
    call finalize_seed

    ! speed etc. information about the run
    write(stdout,form100) 'shots (use -q to enlarge)',shots
    write(stdout,form101) 'cpu+real time spent in overhead [secs]',cpu_and_real_time(1,:)
    write(stdout,form101) 'cpu+real time spent in read_voro [secs]',cpu_and_real_time(2,:)
    write(stdout,form101) 'cpu+real time spent in setup_triangles [secs]',cpu_and_real_time(3,:)
    write(stdout,form101) 'cpu+real time spent in MonteCarlo [secs]',cpu_and_real_time(4,:)
    write(stdout,form101) 'cpu+real time per 10000 shots [secs]',10000.D0*cpu_and_real_time(4,:)/dble(shots)

    ! Tc-PSD and Ts-PSD
    if (TPSD) then
        call calc_TPSD    
        call manage_cpu_and_real_clock(6)
        write(stdout,form101) 'cpu+real time spent in TPSD [secs]',cpu_and_real_time(5,:)
    endif

    open(2,file='_info.txt'); rewind(2)
    write(2,400) 'N=',N
    write(2,401) 'ro=',ro
    write(2,401) 'rp=',rp
    write(2,401) 'rc=',rc
    write(2,401) 'V=',volume
    write(2,400) 'triangles=',triangles
    write(2,400) 'shots=',shots
    write(2,400) 'radii=',count_r
    write(2,401) 'min_pore_radius=',min_pore_radius
    write(2,401) 'max_pore_radius=',UpperPoreRadius
    write(2,401) 'mean_pore_radius=',mean_pore_radius
    write(2,401) 'stderr_pore_radius=',std_pore_radius
    write(2,401) 'phi_reff=',volume_fraction  ! phi(reff)
    write(2,400) 'cells=',Tneighborlist_M(1)*Tneighborlist_M(2)*Tneighborlist_M(3)    ! cells
    write(2,400) 'threads=',np
    write(2,401) 'voro_cpu_time=',cpu_and_real_time(2,1)
    write(2,401) 'voro_real_time=',cpu_and_real_time(2,2)
    write(2,401) 'triangles_setup_cpu_time=',cpu_and_real_time(3,1)
    write(2,401) 'triangles_setup_real-time=',cpu_and_real_time(3,2)
    write(2,401) 'MonteCarlo_cpu_time=',cpu_and_real_time(4,1)
    write(2,401) 'MonteCarlo_real_time=',cpu_and_real_time(4,2)
    close(2)

    open(2,file='_inf.txt'); rewind(2)
    write(2,*) N
    write(2,*) ro
    write(2,*) rp
    write(2,*) rc
    write(2,*) volume
    write(2,*) triangles
    write(2,*) shots
    write(2,*) count_r
    write(2,*) min_pore_radius
    write(2,*) UpperPoreRadius
    write(2,*) mean_pore_radius
    write(2,*) std_pore_radius
    write(2,*) volume_fraction  ! phi(reff)
    write(2,*) Tneighborlist_M(1)*Tneighborlist_M(2)*Tneighborlist_M(3)    ! cells
    write(2,*) np
    write(2,*) cpu_and_real_time(2,1)
    write(2,*) cpu_and_real_time(2,2)
    write(2,*) cpu_and_real_time(3,1)
    write(2,*) cpu_and_real_time(3,2)
    write(2,*) cpu_and_real_time(4,1)
    write(2,*) cpu_and_real_time(4,2)
    close(2)

    open(2,file='.UpperPoreRadius'); rewind(2)
    write(2,*) UpperPoreRadius
    close(2)

400 format(A,I12)
401 format(A,F12.5)
stop
end
