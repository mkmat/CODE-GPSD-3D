module shared
implicit none

    integer         :: stdout
    integer(kind=4) :: N,faces,triangles,shots
    real*8          :: ro,rp,rc,ro_plus_rc,rs,rs2,lo(3),hi(3),box(3),lo_plus_halfbox(3)
    real*8          :: UpperPoreRadius
    integer         :: neighborlist_M(3)
    real*8          :: neighborlist_size(3)
    real*8          :: max_triangle_max_extension
    logical         :: quiet,more
    real*8, dimension(:,:), allocatable  :: triangle_A
    real*8, dimension(:,:), allocatable  :: triangle_B
    real*8, dimension(:,:), allocatable  :: triangle_C
    real*8, dimension(:,:), allocatable  :: triangle_X
    real*8, dimension(:,:), allocatable  :: triangle_M
    real*8, dimension(:,:), allocatable  :: triangle_n
    real*8, dimension(:),   allocatable  :: triangle_rho

    integer, dimension(:,:,:), allocatable  :: neighborlist_firstT
    integer, dimension(:), allocatable      :: neighborlist_nextT

    character(len=100), parameter        :: form100='("[GPSD-3D]",A40,3(I10,1x))'
    character(len=100), parameter        :: form101='("[GPSD-3D]",A40,3(F10.3,1x))'
    character(len=100), parameter        :: form102='("[GPSD-3D]",A40,F10.3,A,F10.3)'
    character(len=100), parameter        :: form200='(4(F12.4,1x))'

    namelist /list/ ro,rp,rc,shots,N,triangles,faces,quiet,more

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
        write(stdout,form101) 'box',box
    return
    end

    subroutine constants
        lo_plus_halfbox = lo+box/2.D0
        ro_plus_rc      = ro+rc
        rs              = ro+rc+rp
        rs2             = rs**2
    return
    end
    
!    function dot(vec1,vec2)
!    real*8 :: dot,vec1(3),vec2(3)
!        dot = vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
!    return
!    end

    subroutine fold(vec)
    real*8  :: vec(3)
        vec  = vec-anint(vec/box)*box
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

    subroutine read_voro_output
    integer(kind=4)     :: tri,face,i,j,id
    real*8  :: coord(3)
    real*8  :: A(3),B(3),C(3)
    real*8  :: normAB,normAC,normBC
    real*8  :: normAM,normBM,normCM
    real*8  :: normBX,normCX,normpM,pTmin,pTmax
    real*8  :: triangle_max_extension

        allocate(triangle_A(triangles,3))
        allocate(triangle_B(triangles,3))
        allocate(triangle_C(triangles,3))
        allocate(triangle_X(triangles,3))
        allocate(triangle_M(triangles,3))
        allocate(triangle_n(triangles,3))
        allocate(triangle_rho(triangles))

        write(stdout,form100) 'reading _triangles.txt'
        open(2,file='_triangles.txt'); rewind(2)
        UpperPoreRadius = 0.D0
        write(stdout,form100) 'reading and processing triangles ..'
        tri = 0
        max_triangle_max_extension = 0.D0
        do face=1,faces
            read(2,*) id,i
            read(2,*) coord
            read(2,*) A
            do j=1,i
                tri = tri + 1
                read(2,*) B,C
                triangle_A(tri,:) = A
                triangle_B(tri,:) = B
                triangle_C(tri,:) = C
                triangle_X(tri,:) = coord
                triangle_n(tri,:) = (coord-A)/norm(coord-A)
                normAB = norm(A-B)
                normAC = norm(A-C)
                normBC = norm(B-C)
                triangle_M(tri,:) = (normAB*C+normAC*B+normBC*A)/(normAB+normAC+normBC)     ! triangle center
                normBX = norm(B-coord)
                normCX = norm(C-coord)
                triangle_rho(tri) = max(normBX,normCX)                                     ! >= r+ro_plus_rc 
                UpperPoreRadius = max(UpperPoreRadius,triangle_rho(tri)-ro_plus_rc)
                normAM = norm(A-triangle_M(tri,:))
                normBM = norm(B-triangle_M(tri,:))
                normCM = norm(C-triangle_M(tri,:))
                triangle_max_extension = max(max(normAM,normBM),normCM)                     ! >= distance from M to any triangle corner
                max_triangle_max_extension = max(max_triangle_max_extension,triangle_max_extension) 
            enddo
        enddo
        close(2)
        write(stdout,form100) 'N',N
        write(stdout,form100) 'triangles',triangles
        write(stdout,form101) 'UpperPoreRadius',UpperPoreRadius
        write(stdout,form101) 'ro',ro
        write(stdout,form101) 'rc',rc
        write(stdout,form101) 'rp',rp
        write(stdout,form101) 'reff = rc+rp',rc+rp
        ! write(stdout,form101) 'max_triangle_max_extension',max_triangle_max_extension
        if (triangles.ne.tri) stop 'format error'

        ! clear histogram
        ! binsize = UpperPoreRadius/bins
        ! PoreRadiusDistribution = 0.D0
    return
    end

    subroutine setup_triangles_neighborlist
    integer(kind=4)     :: i
    integer     :: icell(3)
        write(stdout,form100) 'creating neighbor list ..'

        neighborlist_size   = UpperPoreRadius-rp+max_triangle_max_extension
        neighborlist_M      = max(1,int(box/neighborlist_size))
        neighborlist_size   = box/neighborlist_M

        write(stdout,form100) 'neighborlist_M',neighborlist_M
        write(stdout,form101) 'neighborlist_size',neighborlist_size

        allocate(neighborlist_firstT(neighborlist_M(1),neighborlist_M(2),neighborlist_M(3)))
        allocate(neighborlist_nextT(triangles))
        neighborlist_firstT = 0
        neighborlist_nextT  = 0
        do i=1,triangles
            triangle_M(i,:) = triangle_M(i,:) - anint((triangle_M(i,:)-lo_plus_halfbox)/box)*box + lo_plus_halfbox  ! special fold
            icell = 1+int((triangle_M(i,:)-lo)/neighborlist_size)   
            icell = max(1,min(neighborlist_M,icell))
            neighborlist_nextT(i) = neighborlist_firstT(icell(1),icell(2),icell(3))
            neighborlist_firstT(icell(1),icell(2),icell(3)) = i
        enddo
        deallocate(triangle_M)
    return
    end

    subroutine MonteCarlo
    real    :: tic,toc,myrand3(3)
    real*8  :: r,min_pore_radius,max_pore_radius,mean_pore_radius,std_pore_radius,volume,volume_fraction
    real*8  :: pTmin,pTmax,pTmin_Upper
    real*8  :: p(3),pA(3),center(3)
    real*8  :: A(3),B(3),C(3),X(3)
    real*8  :: facenormal(3)
    integer :: icell(3),jcell(3),ii,i1,i2,i3
    integer :: tmp_type
    real*8  :: tmp_r,tmp_rprime,tmp_center(3)
    integer(kind=4) :: shot,i,count_mean_pore_radius
    integer :: i1min,i2min,i3min
    integer :: i1max,i2max,i3max

        write(stdout,form100) 'start MC ..'
        call cpu_time(tic)

        min_pore_radius        = 100.0*UpperPoreRadius
        max_pore_radius        = -1.0
        mean_pore_radius       = 0.D0
        std_pore_radius        = 0.D0
        count_mean_pore_radius = 0
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
        do shot=1,shots
            call random_number(myrand3)
            r     = 0.D0
            p     = lo + myrand3*box
            icell = 1+int((p-lo)/neighborlist_size)
            do i1=i1min,i1max
            do i2=i2min,i2max
            do i3=i3min,i3max
                jcell(1) = icell(1)+i1
                jcell(2) = icell(2)+i2
                jcell(3) = icell(3)+i3
                do ii=1,3
                    if (jcell(ii).eq.neighborlist_M(ii)+1) then
                        jcell(ii)=1
                    elseif (jcell(ii).eq.0) then
                        jcell(ii)=neighborlist_M(ii)
                    endif
                enddo
                i = neighborlist_firstT(jcell(1),jcell(2),jcell(3))
                do while (i.ne.0) 
                    pTmax = triangle_rho(i)-ro_plus_rc
                    if (pTmax.ge.r) then ! otherwise, we cannot improve r using this triangle
                        A           = triangle_A(i,:)
                        facenormal  = triangle_n(i,:)
                        pA          = p-A
                        pTmin       = norm(pA-facenormal*dot_product(facenormal,pA))
                        if (pTmin.le.pTmin_Upper) then   
                            B = triangle_B(i,:)
                            C = triangle_C(i,:)
                            X = triangle_X(i,:)
                            call candidate_for_triangle(facenormal,A,B,C,X,p, tmp_r,tmp_center,tmp_type,tmp_rprime)
                            if (tmp_r.ge.r) then    
                                r = tmp_r
                                center = tmp_center
                            endif
                        endif
                    endif         
                    i = neighborlist_nextT(i)
                enddo
            enddo
            enddo
            enddo

            if (r.gt.0.D0) then
                if (more) then
                    write(3,form200) r,center
                else
                    write(3,form200) r
                endif
                min_pore_radius         = min(min_pore_radius,r)
                max_pore_radius         = max(max_pore_radius,r)
                mean_pore_radius        = mean_pore_radius + r
                std_pore_radius         = std_pore_radius  + r**2
                count_mean_pore_radius  = count_mean_pore_radius + 1
            endif
        enddo ! shot
        close(3)
        call cpu_time(toc)

        volume_fraction  = 1.D0-dble(count_mean_pore_radius)/dble(shots)
        mean_pore_radius = mean_pore_radius / dble(count_mean_pore_radius)
        std_pore_radius  = std_pore_radius / dble(count_mean_pore_radius)
        std_pore_radius  = (std_pore_radius-mean_pore_radius**2)/dsqrt(dble(count_mean_pore_radius))


        write(stdout,form101) 'cpu per 1000 shots [secs]',1000.D0*(toc-tic)/dble(shots)
        write(stdout,form101) 'volume fraction phi(reff)',volume_fraction
        write(stdout,form101) 'V(0|rc)',(1-volume_fraction)*volume
        write(stdout,form101) 'min pore radius',min_pore_radius
        write(stdout,form102) 'mean pore radius',mean_pore_radius,' +/- ',std_pore_radius
        write(stdout,form101) 'max pore radius',max_pore_radius
        write(stdout,form101) 'created a list {r} of pore radii'
    return
    end

    subroutine candidate_for_triangle(facenormal,A,B,C,X,p, r,center,type,rprime)
        real*8, parameter       :: mkeps    = 1.D-10
        real*8, intent(in)      :: facenormal(3),A(3),B(3),C(3),X(3),p(3)
        real*8, intent(out)     :: r,center(3),rprime
        integer, intent(out)    :: type
        real*8  :: myshift(3)
        real*8  :: ex(3),ey(3),ez(3)
        real*8  :: BB(3),CC(3),XX(3)
        real*8  :: B2,C2
        real*8  :: largeR,largeRprime,largeCenter(3)
        real*8  :: pp(3),ppp(3),p2
        real*8  :: Bx,By,Bz,Cx,Cy,Cz
        real*8  :: Bx2,By2,Cx2,Cy2,BCx,BCy,BCx2,BCy2,BC2
        real*8  :: px,pz,px2,pz2,z,z2,s,normp,Ups,fUps,xp,xm,Q,t1p,t1m,t2p,t2m
        real*8  :: o1,o2,o3

        ! center = xp*ex+y*ey gives center in old coordinate system
        ! center = [xp y 0] gives center in new coordinate system

        ! overlap check
        if (norm2(p-X).lt.rs2) then; type=0; return; endif  

        type    = -1
        r       = 0.D0
        rprime  = 0.D0
        center  = 0.D0

        ! translate all to new origin at A
        myshift = A
        BB = B-A; call fold(BB)
        CC = C-A; call fold(CC)
        XX = X-A; call fold(XX)
        pp = p-A; call fold(pp)

        B2 = norm2(BB)
        C2 = norm2(CC)

        ! type 1 and 2 (triangle corners), return immediately if test passed for longer side 
        if (B2.ge.C2) then
            largeR      = norm(XX-BB)-ro_plus_rc        ! rho-ro-rc
            largeRprime = norm(pp-BB)+rp
            if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) then
                r       = largeR
                rprime  = largeRprime
                center  = B
                type    = 1
                return
            endif
            largeR      = norm(XX-CC)-ro_plus_rc
            largeRprime = norm(pp-CC)+rp
            if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) then
                r       = largeR
                rprime  = largeRprime
                center  = C
                type    = 2
            endif
        else
            largeR      = norm(XX-CC)-ro_plus_rc
            largeRprime = norm(pp-CC)+rp
            if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) then
                r       = largeR
                rprime  = largeRprime
                center  = C
                type    = 1
                return
            endif
            largeR      = norm(XX-BB)-ro_plus_rc
            largeRprime = norm(pp-BB)+rp
            if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) then
                r       = largeR
                rprime  = largeRprime
                center  = B
                type    = 2
            endif
        endif

        ! base system ex,ey,ez
        ez  = XX/norm(XX)
        ppp = pp - dot_product(pp,ez)*ez
        ex  = ppp/norm(ppp) 
        call cross(ez,ex,ey) 

        ! convert all to new base system
        ! A = [dot_product(A,ex),dot_product(A,ey),dot_product(A,ez)]; ! (0,0,0) 
        Bx   = dot_product(BB,ex)
        By   = dot_product(BB,ey)
        Bz   = dot_product(BB,ez)
        Cx   = dot_product(CC,ex)
        Cy   = dot_product(CC,ey)
        Cz   = dot_product(CC,ez)
        px   = dot_product(pp,ex)          ! py=0
        pz   = dot_product(pp,ez)          
        z    = dot_product(XX,ez)          ! x=y=0

        ! abbreviations
        Bx2  = Bx**2
        By2  = By**2
        Cx2  = Cx**2
        Cy2  = Cy**2
        px2  = px**2
        pz2  = pz**2
        p2   = px2+pz2
        z2   = z**2
        s    = px2-rs2
        BCx  = Cx-Bx
        BCy  = Cy-By
        BCx2 = BCx**2
        BCy2 = BCy**2
        BC2  = BCx2+BCy2
        normp = dsqrt(p2)
        fUps = (rs-normp-z)*(rs+normp-z)*(rs-normp+z)*(rs+normp+z)

        ! type 3 grazing contact using y=0
        Q   = z2**2+2.D0*(s-pz2)*z2+(s+pz2)**2 
        if (Q.ge.0.D0) then
            o1  = px*(s+pz2-z2)
            o2  = rs*dsqrt(Q)
            o3  = 2.D0*s 
            xp  = (o1+o2)/o3
            xm  = (o1-o2)/o3
            o1  = Cy/(Bx*Cy-Cx*By)
            o2  = By/(Cx*By-Bx*Cy)
            t1p = xp*o1
            t2p = xp*o2
            t1m = xm*o1
            t2m = xm*o2 
            if (t1p.ge.0.D0.and.t1p.le.1D0.and.t2p.gt.0.D0.and.t2p.le.1.D0-t1p) then    ! inside triangle
                largeCenter = [xp,0.D0,0.D0] 
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime)
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,3,r,rprime,center,type)
            endif
            if (t1m.ge.0D0.and.t1m.le.1.D0.and.t2m.ge.0.D0.and.t2m.le.1.D0-t1m) then    ! inside triangle
                largeCenter = [xm,0.D0,0.D0]
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime)
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,3,r,rprime,center,type)
            endif
        endif

        ! type 4 triangle edge AC
        Ups = 4.D0*Cx2*px2*z2+C2*fUps
        if (Ups.ge.0.D0) then
            o1  = Cx*px*(rs2-p2+z2)
            o2  = rs*dsqrt(Ups)
            o3  = 2.D0*(C2*rs2-Cx2*px2)
            t1p = (o1+o2)/o3
            t1m = (o1-o2)/o3
            if (t1p.ge.0.D0.and.t1p.le.1.D0) then
                largeCenter = [t1p*Cx,t1p*Cy,0.D0] 
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime) 
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,4,r,rprime,center,type)
            endif
            if (t1m.ge.0.D0.and.t1m.le.1.D0) then
                largeCenter = [t1m*Cx,t1m*Cy,0.D0] 
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime)       ! DEBUG THIS COMMAND WAS MISSING
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,4,r,rprime,center,type)
            endif
        endif

        ! type 5 triangle edge AB
        Ups = 4.D0*Bx2*px2*z2+B2*fUps
        if (Ups.ge.0.D0) then
            o1  = Bx*px*(rs2-p2+z2)
            o2  = rs*dsqrt(Ups)
            o3  = 2.D0*(B2*rs2-Bx2*px2)
            t1p = (o1+o2)/o3
            t1m = (o1-o2)/o3
            if (t1p.ge.0.D0.and.t1p.le.1.D0) then
                largeCenter = [t1p*Bx,t1p*By,0.D0]
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime)
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,5,r,rprime,center,type)
            endif
            if (t1m.ge.0.D0.and.t1m.le.1.D0) then
                largeCenter = [t1m*Bx,t1m*By,0.D0]
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime)
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,5,r,rprime,center,type)
            endif
        endif

        ! type 6 triangle edge BC
        Ups = (8*BCy*By*rs2 + 8*BCx*Bx*(rs-px)*(rs+px)-4*BCx*px*(rs2-p2+z2))**2+16.D0*(BCy**2*rs2+BCx**2*(rs2-px2))*(-4.D0*By2*rs2+(rs2-p2)**2+4.D0*Bx2*(-rs2+px2)-2*(rs2+p2)*z2+z2**2+4.D0*Bx*px*(rs2-p2+z2))
        if (Ups.ge.0.D0) then
            o1  = rs2*(-2.D0*BCy*By+BCx*(-2.D0*Bx+px))+BCx*px*(2.D0*Bx*px+z2-p2)
            o2  = dsqrt(Ups)/4.D0
            o3  = 2.D0*BC2*rs2-2*BCx2*px2
            t1p = (o1+o2)/o3
            t1m = (o1-o2)/o3
            if (t1p.ge.0.D0.and.t1p.le.1.D0) then
                largeCenter = [Bx+t1p*BCx,By+t1p*BCy,0.D0] 
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime)
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,6,r,rprime,center,type)
            endif
            if (t1m.ge.0.D0.and.t1m.le.1.D0) then
                largeCenter = [Bx+t1m*BCx,By+t1m*BCy,0.D0] 
                call largeRs(z,px,pz,largeCenter, largeR,largeRprime)
                if (largeR.ge.largeRprime-mkeps.and.largeR.ge.r) call add(largeR,largeRprime,largeCenter,6,r,rprime,center,type)
            endif
        endif

        ! convert center to original coordinate
        center = myshift + center(1)*ex + center(2)*ey + center(3)*ez
    return
    end

    subroutine add(largeR,largeRprime,largeCenter,in,r,rprime,center,type)
    real*8, intent(in)   :: largeR,largeRprime,largeCenter(3)
    integer, intent(in)  :: in
    real*8, intent(out)  :: r,rprime,center(3)
    integer, intent(out) :: type
        r       = largeR
        rprime  = largeRprime
        center  = largeCenter
        type    = in
    return
    end

    subroutine largeRs(z,px,pz,largeCenter, largeR,largeRprime)
    real*8, intent(in)  :: z,px,pz,largeCenter(3)
    real*8, intent(out) :: largeR,largeRprime
    real*8 :: norm
    external norm
        largeR      = dsqrt(largeCenter(1)**2+largeCenter(2)**2+z**2)-ro_plus_rc    ! X-Center
        largeRprime = dsqrt((px-largeCenter(1))**2+largeCenter(2)**2+pz**2)+rp      ! p-Center
    return
    end

end module shared

program GPSD_3D
use shared
    real :: cpu_read_voro_output
    real :: cpu_setup_triangles
    real :: cpu_MonteCarlo

    call read_parameters
    call constants
    call read_box
    call cpu_time(tic); call read_voro_output; call cpu_time(toc); cpu_read_voro_output=toc-tic
    call cpu_time(tic); call setup_triangles_neighborlist; ; call cpu_time(toc); cpu_setup_triangles=toc-tic
    call init_seed
    call cpu_time(tic); call MonteCarlo; call cpu_time(toc); call cpu_time(toc); cpu_MonteCarlo=toc-tic
    call finalize_seed

    write(stdout,form100) 'shots (use -q to enlarge)',shots
    write(stdout,form101) 'time spent in read_voro_output [secs]',cpu_read_voro_output
    write(stdout,form101) 'time spent in setup_triangles [secs]',cpu_setup_triangles
    write(stdout,form101) 'time spent in MonteCarlo [secs]',cpu_MonteCarlo

stop
end
