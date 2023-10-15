module shared
implicit none

    integer         :: stdout,np
    integer(kind=4) :: N,triangles,shots,count_r
    real*8          :: ro,rp,rc,ro_plus_rc,rs,rs2,rs24,lo(3),hi(3),box(3),lo_plus_halfbox(3),volume,volume_fraction
    real*8          :: UpperPoreRadius,min_pore_radius,max_pore_radius,mean_pore_radius,std_pore_radius
    integer         :: neighborlist_M(3)
    real*8          :: neighborlist_size(3),max_triangle_max_extension
    logical         :: quiet,more

    real*8, dimension(:,:), allocatable  :: triangle_o
    real*8, dimension(:,:), allocatable  :: triangle_n
    real*8, dimension(:,:), allocatable  :: triangle_A
    real*8, dimension(:,:), allocatable  :: triangle_B
    real*8, dimension(:,:), allocatable  :: triangle_C
    real*8, dimension(:,:), allocatable  :: triangle_vM
    real*8, dimension(:),   allocatable  :: triangle_Xn
    real*8, dimension(:),   allocatable  :: triangle_rho

    integer, dimension(:,:,:), allocatable  :: neighborlist_firstT
    integer, dimension(:), allocatable      :: neighborlist_nextT

    character(len=100), parameter        :: form100='("[GPSD-3D]",A45,3(I13,1x))'
    character(len=100), parameter        :: form101='("[GPSD-3D]",A45,3(F13.4,1x))'
    character(len=100), parameter        :: form102='("[GPSD-3D]",A45,F13.4,A3,F11.4)'
    character(len=100), parameter        :: form104='("[VORO++] ",A45,3(I13,1x))'
    character(len=100), parameter        :: form200='(4(F12.5,1x),I2)'
    character(len=100), parameter        :: form201='(F12.5)'

    real cpu_and_real_time(5,2)

    namelist /list/ ro,rp,rc,shots,N,quiet,more,np

    contains 

    subroutine read_parameters
        open(2,file='.parameters'); rewind(2)
        read(2,nml=list)
        close(2)
        if (quiet) then
            stdout = 9
            open(stdout,file='log-GPSD-3D.txt'); rewind(stdout)
        else
            stdout = 6
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
    
!    function dot(vec1,vec2)
!    real*8 :: dot,vec1(3),vec2(3)
!        dot = vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
!    return
!    end

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
        ! print *,'DEBUG ',realcount,countrate,realcount/dble(countrate)
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
    integer oldseed(20)
        inquire(file='_seed',EXIST=EX)
        if (EX) then
            open(2,file='_seed'); rewind(2); read(2,*) oldseed; close(2)
            open(2,file='_oldseed'); rewind(2); write(2,*) oldseed; close(2)
            call random_seed(put=oldseed)
        else
            call random_seed()
        endif
    return
    end

    subroutine finalize_seed
    integer oldseed(20)
        call random_seed(get=oldseed)
        open(2,file='_seed'); rewind(2); write(2,*) oldseed; close(2)
    return
    end

    subroutine read_voro_parser
    ! parser
    integer :: voro_max_vertices,voro_max_faces,voro_max_vertices_for_face
    integer :: voro_vertices,voro_faces
    integer :: itmp(2000),j2
    integer, dimension(:), allocatable   :: voro_vertices_for_face
    integer, dimension(:,:), allocatable :: voro_vid
    real*8, dimension(:,:), allocatable  :: voro_vertex
    ! triangles
    integer(kind=4)     :: tri,face,i,j,id
    real*8  :: coord(3)
    real*8  :: A(3),B(3),C(3),E(3),X(3),facenormal(3)
    real*8  :: CB(3),CE(3),XA(3)
    real*8  :: normGA,normGB,normGC
    real*8  :: normBX,normCX,pTmax
    real*8  :: triangle_max_extension
    ! parser
    real*8  :: vA(3),vB(3),vC(3),vE(3),o(3),ref(3)
    real*8  :: GeomCenter(3)

    ! read once to determine field sizes
    voro_max_vertices           = 0
    voro_max_faces              = 0
    voro_max_vertices_for_face  = 0
    triangles                   = 0

    read(*,*) N        ! attention: overwriting N as seen by voro++
    do id=1,N
        read(*,*) coord,voro_vertices,voro_faces,(itmp(face),face=1,voro_faces)
        read(*,*) 
        voro_max_vertices = max(voro_max_vertices,voro_vertices)
        voro_max_faces    = max(voro_max_faces,voro_faces)
        do face=1,voro_faces
            voro_max_vertices_for_face = max(voro_max_vertices_for_face,itmp(face))
            triangles = triangles + itmp(face)
        enddo
    enddo

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
    allocate(triangle_rho(triangles))   ! rho

    tri = 0
    UpperPoreRadius = 0.D0
    max_triangle_max_extension = 0.D0

    read(*,*) itmp(1)        ! attention: overwriting N as seen by voro++
    if (N.ne.itmp(1)) stop 'parser error (called twice?)'
    do id=1,N
         ! I am assuming absolute, folded, vertex coordinates here
        read(*,*) coord,voro_vertices,voro_faces,(voro_vertices_for_face(face),face=1,voro_faces)
        read(*,*) (voro_vertex(j,:),j=1,voro_vertices),((voro_vid(face,j),j=1,voro_vertices_for_face(face)),face=1,voro_faces)   
        do face=1,voro_faces
            ! determine geometrical face center at vA
            ref = voro_vertex(voro_vid(face,1),:)
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
            facenormal = facenormal/norm(facenormal)                ! sign does not matter here
            XA  = distance_fold(coord-vA)
            o   = coord-dot_product(XA,facenormal)*facenormal   ! origin
            o   = position_fold(o)
            A   = distance_fold(vA-o)      
            X   = distance_fold(coord-o)
            facenormal = X/norm(X)                                  ! as the sign is used elsewhere
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
                triangle_Xn(tri)   = dot_product(X,facenormal)
                ! geometric center
                GeomCenter = (B+C)/3.D0     ! relative to vertex A
                triangle_vM(tri,:) = position_fold(vA + GeomCenter)
                ! max distance vertex - geometric center
                normGA = norm(GeomCenter)
                normGB = norm(GeomCenter-B)
                normGC = norm(GeomCenter-C)
                triangle_max_extension = (2.D0/3.D0)*max(normGA,max(normGB,normGC))
                max_triangle_max_extension = max(max_triangle_max_extension,triangle_max_extension)
                normBX  = norm(A+B-X)
                normCX  = norm(A+C-X)
                triangle_rho(tri) = max(normBX,normCX)-ro_plus_rc                    
                UpperPoreRadius = max(UpperPoreRadius,triangle_rho(tri))
            enddo
        enddo
    enddo
    write(stdout,form100) 'parallel processes (np)',np
    write(stdout,form100) 'material spheres (N)',N
    write(stdout,form101) 'UpperPoreRadius',UpperPoreRadius     ! this value is correct
    write(stdout,form101) 'ro',ro
    write(stdout,form101) 'rc',rc
    write(stdout,form101) 'rp',rp
    write(stdout,form101) 'reff = rc+rp',rc+rp
    write(stdout,form101) 'rs = ro+reff',rs
    write(stdout,form101) 'max triangle extension',max_triangle_max_extension
    if (triangles.ne.tri) stop 'format error'
    
    deallocate(voro_vertices_for_face)
    deallocate(voro_vid)
    deallocate(voro_vertex)

    return
    end

    subroutine setup_triangles_neighborlist
    integer(kind=4)     :: i
    integer     :: icell(3),ncells
        write(stdout,form100) 'creating neighbor list ..'

        neighborlist_size   = UpperPoreRadius-rp+max_triangle_max_extension
        neighborlist_M      = max(1,int(box/neighborlist_size))
        neighborlist_size   = box/neighborlist_M
        ncells              = neighborlist_M(1)*neighborlist_M(2)*neighborlist_M(3)

        write(stdout,form100) 'neighborlist_M',neighborlist_M
        write(stdout,form101) 'neighborlist_size',neighborlist_size
        write(stdout,form101) 'triangles per cell',triangles/dble(ncells)

        allocate(neighborlist_firstT(neighborlist_M(1),neighborlist_M(2),neighborlist_M(3)))
        allocate(neighborlist_nextT(triangles))
        neighborlist_firstT = 0
        neighborlist_nextT  = 0
        do i=1,triangles
            icell = 1+int((triangle_vM(i,:)-lo)/neighborlist_size)   
            icell = max(1,min(neighborlist_M,icell))
            neighborlist_nextT(i) = neighborlist_firstT(icell(1),icell(2),icell(3))
            neighborlist_firstT(icell(1),icell(2),icell(3)) = i
        enddo
        deallocate(triangle_vM)
    return
    end

    subroutine MonteCarlo
    real    :: myrand3(3)
    real*8  :: R
    real*8  :: pTmax,pTmin_Upper
    real*8  :: vp(3),P(3),Pn,center(3)
    real*8  :: o(3),A(3),B(3),C(3),X(3),Xn
    real*8  :: facenormal(3)
    integer :: icell(3),jcell(3),k,i1,i2,i3
    real*8  :: tmp_r,tmp_center(3)
    integer(kind=4) :: shot,i
    integer :: i1min,i2min,i3min
    integer :: i1max,i2max,i3max

        write(stdout,form100) 'starting Monte Carlo ..'

        min_pore_radius        = 100.0*UpperPoreRadius
        max_pore_radius        = -1.0
        mean_pore_radius       = 0.D0
        std_pore_radius        = 0.D0
        count_r                = 0
        pTmin_Upper            = UpperPoreRadius+ro_plus_rc-rp
        volume                 = box(1)*box(2)*box(3)

        call init_seed

        i1min = -1; i1max = 1
        i2min = -1; i2max = 1
        i3min = -1; i3max = 1
        if (neighborlist_M(1).eq.1) i1min=1
        if (neighborlist_M(1).eq.2) i1min=0
        if (neighborlist_M(2).eq.1) i2min=1
        if (neighborlist_M(2).eq.2) i2min=0
        if (neighborlist_M(3).eq.1) i3min=1
        if (neighborlist_M(3).eq.2) i3min=0

        open(3,file='_r'); rewind(3)
!$omp parallel num_threads ( np) shared ( i1min, i1max, i2min, i2max, neighborlist_M, neighborlist_size, triangle_rho, triangle_B, triangle_C, triangle_A ,triangle_n, triangle_Xn, neighborlist_nextT, lo, box, more ) private ( shot, myrand3, vp, P, Pn, o, icell, jcell, k, i, pTmax, A, facenormal, B, C, tmp_r, tmp_center, r, center, X, Xn )

!$omp do reduction (min : min_pore_radius) reduction (max : max_pore_radius) reduction( + : mean_pore_radius, count_r, std_pore_radius)
        do shot=1,shots
            call random_number(myrand3)
            R     = 0.D0
            vp    = lo + myrand3*box
            icell = 1+int((vp-lo)/neighborlist_size)
            do i1=i1min,i1max
            do i2=i2min,i2max
            do i3=i3min,i3max
                jcell(1) = icell(1)+i1
                jcell(2) = icell(2)+i2
                jcell(3) = icell(3)+i3
                do k=1,3
                    if (jcell(k).eq.neighborlist_M(k)+1) then
                        jcell(k)=1
                    elseif (jcell(k).eq.0) then
                        jcell(k)=neighborlist_M(k)
                    endif
                enddo
                i = neighborlist_firstT(jcell(1),jcell(2),jcell(3))
                do while (i.ne.0) 
                    pTmax = triangle_rho(i)
                    if (pTmax.ge.R) then ! otherwise, we cannot improve R using this triangle
                        facenormal  = triangle_n(i,:)
                        o           = triangle_o(i,:)
                        P           = distance_fold(vp-o)
                        Pn          = dot_product(P,facenormal)
                        if (dabs(Pn).le.pTmin_Upper) then      ! EDIT
                            Xn = triangle_Xn(i)
                            X  = Xn*facenormal
                            if (norm2(P-X).ge.rs2) then
                                A  = triangle_A(i,:)
                                B  = triangle_B(i,:)
                                C  = triangle_C(i,:)
                                call candidate_for_triangle(facenormal,o,A,B,C,X,Xn,P, tmp_R,tmp_center)
                                if (tmp_R.gt.R) then    
                                    R = tmp_R       ! this is r-rp
                                    center = tmp_center
                                endif
                            endif
                        endif
                    endif         
                    i = neighborlist_nextT(i)
                enddo
            enddo
            enddo
            enddo

            if (R.gt.0.D0) then
                r = R+rp
                if (more) then
                    write(3,form200) r,position_fold(o+center)
                else
                    write(3,form201) r
                endif
                min_pore_radius         = min(min_pore_radius,r)
                max_pore_radius         = max(max_pore_radius,r)
                mean_pore_radius        = mean_pore_radius + r
                std_pore_radius         = std_pore_radius  + r**2
                count_r                 = count_r + 1
            endif
        enddo ! shot
!$omp end do
!$omp end parallel
        close(3)

        volume_fraction  = 1.D0-dble(count_r)/dble(shots)
        mean_pore_radius = mean_pore_radius / dble(count_r)
        std_pore_radius  = std_pore_radius / dble(count_r)
        std_pore_radius  = (std_pore_radius-mean_pore_radius**2)/dsqrt(dble(count_r))


        write(stdout,form101) 'volume fraction phi(reff)',volume_fraction
        write(stdout,form101) 'V(0|reff)',(1-volume_fraction)*volume
        write(stdout,form101) 'min pore radius',min_pore_radius
        write(stdout,form102) 'mean pore radius',mean_pore_radius,'+/-',std_pore_radius
        write(stdout,form101) 'max pore radius',max_pore_radius
        write(stdout,form101) 'created a list {r} of pore radii'
    return
    end

    subroutine candidate_for_triangle(facenormal,o,A,B,C,X,Xn,P, R,center)
        real*8, parameter       :: mkeps    = 1.D-10
        real*8, intent(in)      :: o(3),facenormal(3),A(3),B(3),C(3),X(3),Xn,P(3)
        real*8, intent(out)     :: R,center(3)
        real*8  :: ex(3),ey(3),ez(3)
        real*8  :: P2,test_center(3)
        real*8  :: R1,R2
        real*8  :: Ap,Aq,Bp,Bq,Cp,Cq,Pn,Pp
        real*8  :: t(2),u(2)

        r           = 0.D0
        center      = 0.D0

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

        ! type 1 (corner B)
        test_center = A+B
        R1 = norm(X-test_center)-rs
        if (R1.ge.r) then
            R2 = norm(P-test_center)
            if (R1.ge.R2-mkeps) then
                R        = R1
                center   = test_center
            endif
        endif
        ! type 2 (corner C)
        test_center  = A+C
        R1 = norm(X-test_center)-rs
        if (R1.ge.r) then
            R2  = norm(P-test_center)
            if (R1.ge.R2-mkeps) then
                R        = R1
                center   = test_center
            endif
        endif

        ! abbreviations
        P2 = Pn**2+Pp**2

        ! case iii (edge B-C)
        call make_case_iii(Ap,Aq,Bp,Bq,Cp,Cq,Xn,Pp,P2, t,u)
        if (t(1).ge.0.D0.and.t(1).le.1.D0) then
            test_center = A+t(1)*B+u(1)*C
            R1      = norm(X-test_center)-rs
            R2 = norm(P-test_center)                    
            if (R1.ge.r.and.R1.ge.R2) then
                R        = R1
                center   = test_center          ! type 3
            endif
        endif
        if (t(2).ge.0.D0.and.t(2).le.1.D0) then
            test_center = A+t(2)*B+u(2)*C
            R1      = norm(X-test_center)-rs
            R2 = norm(P-test_center)
            if (R1.ge.r.and.R1.ge.R2) then
                R        = R1
                center   = test_center          ! type 4
            endif
        endif

        ! case iv (at Sq=0) 
        if (dabs(Cq).gt.0.D0) then
            call make_case_iv(Ap,Aq,Bp,Bq,Cp,Cq,Xn,Pp,P2, t,u)
            if (t(1).ge.0.D0.and.t(1).le.1.D0.and.u(1).ge.0.D0.and.u(1).le.1.D0.and.t(1)+u(1).le.1.D0) then
                test_center = A+t(1)*B+u(1)*C
                R1      = norm(X-test_center)-rs
                R2 = norm(P-test_center)
                if (R1.ge.r.and.R1.ge.R2) then
                    R        = R1
                    center   = test_center      ! type 5
                endif
            endif
            if (t(2).ge.0.D0.and.t(2).le.1.D0.and.u(2).ge.0.D0.and.u(2).le.1.D0.and.t(2)+u(2).le.1.D0) then
                test_center = A+t(2)*B+u(2)*C
                R1      = norm(X-test_center)-rs
                R2 = norm(P-test_center)
                if (R1.ge.r.and.R1.ge.R2) then
                    R        = R1
                    center   = test_center      ! type 6
                endif
            endif
        endif

    return
    end

    subroutine make_case_iii(Ap,Aq,Bp,Bq,Cp,Cq,Xn,Pp,P2, t,u)
    real*8  :: Ups,t(2),u(2)
    real*8  :: Ap,Aq,Bp,Bq,Cp,Cq,Xn,Pp,P2
    real*8  :: a,b,c,d,e
    real*8  :: b2e,ApCp,AqCq,BpCp,BqCq,sqrtUps
        ApCp = Ap+Cp
        AqCq = Aq+Cq
        BpCp = Bp-Cp
        BqCq = Bq-Cq
        a = 2.D0*ApCp*Pp+rs2+Xn**2-P2;
        b = 2.D0*BpCp*Pp;
        c = rs24*(Xn**2+ApCp**2+AqCq**2)        ! c = 4*rs2*(X2+dot(A+C,A+C));
        d = rs24*(BpCp*ApCp+BqCq*AqCq)          ! d = 4*rs2*dot(B-C,A+C);
        e = rs24*(BpCp**2+BqCq**2)              ! e = 4*rs2*dot(B-C,B-C);
        b2e = b**2-e
        Ups = a**2*e-2.D0*a*b*d+d**2+b2e*c;
        if (Ups.ge.0.D0) then
            sqrtUps = dsqrt(Ups)
            t(1) = (d-a*b+sqrtUps)/b2e
            t(2) = (d-a*b-sqrtUps)/b2e
            u    = 1.D0-t
        else
            t = -1.D0
            u = -1.D0
        endif
    return
    end

    subroutine make_case_iv(Ap,Aq,Bp,Bq,Cp,Cq,Xn,Pp,P2, t,u)
    real*8  :: Ups,t(2),u(2)
    real*8  :: Ap,Aq,Bp,Bq,Cp,Cq,Xn,Pp,P2
    real*8  :: a,b,c,d,e,f,h
    real*8  :: b2e,AphCp,BpfCp,BqfCq,sqrtUps
        h = Aq/Cq
        f = Bq/Cq 
        AphCp = Ap-h*Cp
        BpfCp = Bp-f*Cp
        BqfCq = Bq-f*Cq
        a = 2.D0*AphCp*Pp+rs2+Xn**2-P2
        b = 2.D0*BpfCp*Pp;
        c = rs24*(Xn**2+AphCp**2);
        d = rs24*AphCp*BpfCp   
        e = rs24*( BpfCp**2 + BqfCq**2 ) 
        b2e = b**2-e
        Ups = a**2*e-2.D0*a*b*d+d**2+b2e*c 
        if (Ups.ge.0.D0) then
            sqrtUps = dsqrt(Ups)
            t(1) = (d-a*b+sqrtUps)/b2e
            t(2) = (d-a*b-sqrtUps)/b2e
            u    = -h-f*t
        else
            t = -1.D0
            u = -1.D0
        endif
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
    call manage_cpu_and_real_clock(3); call setup_triangles_neighborlist
    call manage_cpu_and_real_clock(4); call MonteCarlo
    call manage_cpu_and_real_clock(5)
    call finalize_seed

    ! speed etc. information about the run
    write(stdout,form100) 'shots (use -q to enlarge)',shots
    write(stdout,form101) 'cpu+real time spent in overhead [secs]',cpu_and_real_time(1,:)
    write(stdout,form101) 'cpu+real time spent in read_voro [secs]',cpu_and_real_time(2,:)
    write(stdout,form101) 'cpu+real time spent in setup_triangles [secs]',cpu_and_real_time(3,:)
    write(stdout,form101) 'cpu+real time spent in MonteCarlo [secs]',cpu_and_real_time(4,:)
    write(stdout,form101) 'cpu+real time per 10000 shots [secs]',10000.D0*cpu_and_real_time(4,:)/dble(shots)


    open(2,file='_summary.txt'); rewind(2)
    write(2,*) N
    write(2,*) ro
    write(2,*) rp
    write(2,*) rc
    write(2,*) volume
    write(2,*) triangles
    write(2,*) shots
    write(2,*) count_r
    write(2,*) UpperPoreRadius
    write(2,*) min_pore_radius
    write(2,*) mean_pore_radius
    write(2,*) max_pore_radius
    write(2,*) std_pore_radius
    write(2,*) volume_fraction  ! phi(reff)
    write(2,*) neighborlist_M(1)*neighborlist_M(2)*neighborlist_M(3)    ! cells
    write(2,*) np
    write(2,*) cpu_and_real_time(2,1)
    write(2,*) cpu_and_real_time(2,2)
    write(2,*) cpu_and_real_time(3,1)
    write(2,*) cpu_and_real_time(3,2)
    write(2,*) cpu_and_real_time(4,1)
    write(2,*) cpu_and_real_time(4,2)
    close(2)

stop
end
