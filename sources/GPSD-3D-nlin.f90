module shared

! 1 dec 2023 mk@mat.ethz.ch employing the libnlopt.so library

! routines available in this code: 

! subroutine f(R1, ndim, c, grad, need_gradient, f_data) defines the nonlinear optimization target
! subroutine fc(val, ndim, c, grad, need_gradient, p) defines the nonlinear optimization problem
! subroutines init_seed, finalize_seed handle the random number seed
! subroutine check_result(p,c) is not in use
! subroutine manage_cpu_and_real_clock(no,cpu_and_real_time) measures and reports cpu time

! functions avaiable in this code:
! function smallest_coated_surface_distance_allspheres_point(p) calculates smallest point-surface distances

   integer, parameter :: dim = 3
   integer  :: N
   real*8, dimension(:,:), allocatable :: x
   real*8, dimension(:), allocatable   :: radius
   real*8   :: box(dim)
   real*8   :: rc,rp,max_radius
   real*8   :: xlo,xhi,ylo,yhi,zlo,zhi

   integer, dimension(:,:,:), allocatable    :: firstbead
   integer, dimension(:), allocatable        :: nextbead
   real*8   :: neighborlist_size(3)
   integer  :: neighborlist_M(3),bin(3),binstart(3),binend(3)

   contains

   function smallest_coated_surface_distance_allspheres_point(p)
   implicit none
   real*8   :: smallest_coated_surface_distance_allspheres_point
   real*8   :: p(3),px(3),distance
   integer  :: bin(3)
   integer  :: i1,i2,i3,j1,j2,j3,id
      smallest_coated_surface_distance_allspheres_point = 1.D30
      bin = 1+int((p-[xlo,ylo,zlo])/neighborlist_size)
      bin = max(1,min(neighborlist_M,bin))
      do i1=binstart(1),binend(1)
         j1=bin(1)+i1
         if (j1.eq.0) j1=neighborlist_M(1)
         if (j1.eq.neighborlist_M(1)+1) j1=1
         do i2=binstart(2),binend(2)
            j2=bin(2)+i2
            if (j2.eq.0) j2=neighborlist_M(2)
            if (j2.eq.neighborlist_M(2)+1) j2=1
            do i3=binstart(3),binend(3)
               j3=bin(3)+i3
               if (j3.eq.0) j3=neighborlist_M(3)
               if (j3.eq.neighborlist_M(3)+1) j3=1
               id = firstbead(j1,j2,j3)
               do while (id.ne.0)
                  px = p-x(id,:)
                  px = px-box*anint(px/box)
                  distance = dsqrt(dot_product(px,px))-radius(id)-rc
                  if (distance.lt.smallest_coated_surface_distance_allspheres_point) then
                     smallest_coated_surface_distance_allspheres_point = distance
                  endif
                  id = nextbead(id)
               enddo
            enddo
         enddo
      enddo
   return
   end
   

end module shared

     program GPSD_nlin
     use shared
     implicit none

! (c) 13 nov 2023 mk@mat.ethz.ch

! GPSD via constrained nonlinear optimization

! requires nlopt-2.7.1 or higher
! https://pypi.org/project/nlopt/

! export LD_LIBRARY_PATH=/home/mkroeger/PC/vsandra/NONLINEAR-OPTIMIZATION/nlopt-2.7.1/build;
! gfortran -Ofast -ffree-line-length-none -L nlopt-2.7.1/build -lnlopt -lm GPSD-3D-nlin.f90 -o GPSD-3D-nlin.ex

     external f, fc
     real*8, parameter  :: upper_limit_maxtime = 1.D0
     real*8, parameter  :: pi = acos(-1.0)
     real*8, parameter  :: pi2 = 2.D0*pi
     real*8  :: lb(dim),ub(dim)
     integer*8 :: opt,local_opt
     real*8  :: c(dim), R,px(dim),distance,cx(dim)
     integer :: ires,j,id,count_no_space
     real*8  :: f_data,p(dim)
     logical :: EX,quiet,info,more,use_radii_from_config,use_shots_list,c_outside_rs,bad
     real*8  :: ro
     real*8  :: phi,z,phi_reff
     real    :: myrand3(3)
     real    :: myrand2(2)
     real    :: cpu_and_real_time(2,2)
     real    :: cputime
     real*8  :: maxtime,tic,toc
     real*8  :: min_pore_radius,max_pore_radius,mean_pore_radius,stderr_pore_radius,known_max_pore_radius
     integer(kind=4) :: shots,shots_target,material_shots,pore_radii,bads

     integer, parameter :: NLOPT_LN_SBPLX=29
     integer, parameter :: NLOPT_AUGLAG=36

     namelist /list/ N,ro,rc,rp,shots,quiet,info,more,use_radii_from_config,use_shots_list,known_max_pore_radius
    
     call init_seed
     call manage_cpu_and_real_clock(1,cpu_and_real_time)

     f_data = 0.D0

     ! read parameters
     inquire(file='parameters-nlin.txt',exist=EX)
     if (.not.EX) stop 'missing parameters-nlin.txt'
     known_max_pore_radius = -1.D0
     open(2,file='parameters-nlin.txt'); rewind(2)
     read(2,nml=list)
     close(2)
     if (.not.quiet) then
         print 100,'N',N
         print 101,'ro',ro
         print 101,'rc',rc
         print 101,'rp',rp
     endif

     ! read box
     inquire(file='box.txt',exist=EX)
     if (.not.EX) stop 'missing box.txt'
     open(2,file='box.txt'); rewind(2)
     read(2,*) xlo,xhi,ylo,yhi,zlo,zhi
     close(2)
     box = [xhi-xlo,yhi-ylo,zhi-zlo]
     if (.not.quiet) print 101,'box',box

     if (known_max_pore_radius.lt.0.D0) then
         known_max_pore_radius = max(box(1),max(box(2),box(3)))
     endif
     if (.not.quiet) print 101,'known max pore radius',known_max_pore_radius

     ! read configuration
     allocate(x(N,dim))
     allocate(radius(N))
     max_radius = 0.D0
     if (use_radii_from_config) then
      inquire(file='config-with-radii.txt',exist=EX)
      if (.not.EX) stop 'missing config-with-radii.txt'
      open(2,file='config-with-radii.txt'); rewind(2)
      do id=1,N
         read(2,*) j,x(id,:),radius(id)
         max_radius = max(max_radius,radius(id))
      enddo
     else
      inquire(file='config.txt',exist=EX)
      if (.not.EX) stop 'missing config.txt'
      open(2,file='config.txt'); rewind(2)
      do id=1,N
         read(2,*) j,x(id,:)
         radius(id) = ro
         max_radius = max(max_radius,radius(id))
      enddo
     endif
     close(2)

     ! setup neighbor list
     neighborlist_size = min(box,0*bin+known_max_pore_radius+max_radius+rc)
     neighborlist_M    = max(1,int(box/neighborlist_size))
     neighborlist_size = box/neighborlist_M
     allocate(firstbead(neighborlist_M(1),neighborlist_M(2),neighborlist_M(3)))
     allocate(nextbead(N))
     firstbead = 0
     nextbead = 0
     do id=1,N
       bin = 1+int((x(id,:)-[xlo,ylo,zlo])/neighborlist_size) 
       bin = max(1,min(neighborlist_M,bin))
       nextbead(id) = firstbead(bin(1),bin(2),bin(3))
       firstbead(bin(1),bin(2),bin(3)) = id
     enddo
     binstart = -1
     binend   = 1
     do id=1,3
      if (neighborlist_M(id).eq.1) then
         binstart(id) = 0
         binend(id) = 0
      elseif (neighborlist_M(id).eq.2) then
         binstart(id) = 0
         binend(id) = 1
      endif
     enddo
     if (.not.quiet) then
        print 100,'neighborlist_M',neighborlist_M
        print 101,'neighborlist_size',neighborlist_size
        print 100,'binstart',binstart
        print 100,'binend',binend
     endif
      

     open(2,file='radii.txt'); rewind(2)
     material_shots = 0
     pore_radii = 0
     min_pore_radius = 1.D30
     max_pore_radius = 0.D0
     mean_pore_radius = 0.D0
     stderr_pore_radius = 0.D0

     if (use_shots_list) then
         open(20,file='p.txt'); rewind(20)
     endif

     ! reset counters
     shots_target = shots
     shots = 0
     bads  = 0
     call cpu_time(tic)

     ! initialize maxtime
     maxtime = max(0.5D0,dble(N)/dble(1e5))          ! default, will be overwritten automatically
     if (.not.quiet) print 101,'maxtime (init)',maxtime

     do while (shots.lt.shots_target) 

         if (use_shots_list) then
            read(20,*,end=1) id,p
         else
            call random_number(myrand3)
            p = [xlo,ylo,zlo]+myrand3*box
         endif

         ! skip and count if inside coated material
         if (smallest_coated_surface_distance_allspheres_point(p).le.0.D0) then
            shots = shots + 1
            material_shots = material_shots + 1
         else
            call nlo_create(local_opt, NLOPT_LN_SBPLX, dim)    ! new
            call nlo_create(opt, NLOPT_AUGLAG, dim)     ! new
            call nlo_set_local_optimizer(ires, opt, local_opt) 

            lb = [xlo,ylo,zlo]
            call nlo_set_lower_bounds(ires, opt, lb)
            ub = [xhi,yhi,zhi]
            call nlo_set_upper_bounds(ires, opt, ub)
         
            call nlo_set_max_objective(ires, opt, f, f_data)
            call nlo_add_inequality_constraint(ires, opt, fc, p, 1.D-8)
            call nlo_set_xtol_rel(ires, opt, 1.D-4)  
            call nlo_set_maxtime(ires, opt, maxtime) 
         
            ! choose c somewhere outside
            c_outside_rs = .false.
            count_no_space = 0
            do while (.not.c_outside_rs) 
               count_no_space = count_no_space + 1
               call random_number(myrand2)
               phi = pi2*myrand2(1)
               z   = -1.D0+2.D0*myrand2(2)
               c = p + 1.001*rp*[cos(phi)*dsqrt(1.D0-z**2),sin(phi)*dsqrt(1.D0-z**2),z]
               c_outside_rs = .true.
               if (smallest_coated_surface_distance_allspheres_point(c).le.rp) c_outside_rs=.false.
               if (count_no_space.eq.100) then  ! give up
                  shots = shots + 1
                  goto 2
               endif
            enddo
         
            ! calling the nonlinear optimizer
            call nlo_optimize(ires, opt, c, R)
         
            bad = .false.
            if (ires.le.0.or.ires.ge.5) bad=.true.
            if (ires.eq.6) then
               maxtime = min(1.1*maxtime,upper_limit_maxtime)
               if (maxtime.lt.upper_limit_maxtime) print 101,'WARNING large maxtime',maxtime
            endif

            if (bad) then
               bads = bads + 1
            else
               shots = shots + 1
               pore_radii = pore_radii + 1
               r = R+rp
               if (more) then
                  write(2,'(I10,7(F14.6,1x))') pore_radii,p,c,r
               else
                  write(2,*) r
               endif
               if (r.lt.min_pore_radius) then
                  min_pore_radius = r
               elseif (r.gt.max_pore_radius) then
                  max_pore_radius = r
               endif
               mean_pore_radius = mean_pore_radius + r
               stderr_pore_radius = stderr_pore_radius + r**2
            endif
2        continue
         endif

         if (shots.eq.100) then
            call cpu_time(toc)
            maxtime = min(upper_limit_maxtime,2*(toc-tic)/dble(max(1,pore_radii)))
            print 101,'tic toc',tic,toc
            print 100,'shots so far',shots
            print 100,'pore radii so far',pore_radii
            print 101,'nlin maxtime',maxtime
         endif

      enddo 
1    continue
     close(2)
     close(20)
     
     phi_reff = dble(material_shots)/dble(shots)
     mean_pore_radius = mean_pore_radius/dble(max(1,pore_radii))
     stderr_pore_radius = stderr_pore_radius/dble(max(1,pore_radii))
     stderr_pore_radius = dsqrt(stderr_pore_radius-mean_pore_radius**2)/dsqrt(dble(max(1,pore_radii)))
     if (.not.quiet) then
         print 100,'material shots',material_shots
         print 100,'pore radii',pore_radii
         print 100,'bads (ignored)',bads 
         print 101,'phi_reff',phi_reff
     endif
     
     call nlo_destroy(opt)
     call finalize_seed
     call manage_cpu_and_real_clock(2,cpu_and_real_time)

     if (.not.quiet) then
         print 101,'maxtime (final)',maxtime
         print 100,'shots',shots
         print 100,'material shots',material_shots
         print 100,'pore radii',pore_radii
         print 100,'bads',bads
         print 101,'phi_reff',phi_reff
         print 101,'min_pore_radius',min_pore_radius
         print 101,'mean_pore_radius',mean_pore_radius
         print 101,'max_pore_radius',max_pore_radius
         print 101,'stderr_pore_radius',stderr_pore_radius
         print 101,'cpu time',cpu_and_real_time(1,1)
         print 101,'real time',cpu_and_real_time(1,2)
     endif

     if (info) then
      ! see https://github.com/mkmat/CODE-GPSD-3D

      open(2,file='radii.info'); rewind(2)
      write(2,400) 'N=',N
      write(2,401) 'ro=',ro
      write(2,401) 'rp=',rp
      write(2,401) 'rc=',rc
      write(2,401) 'V=',box(1)*box(2)*box(3)
      write(2,400) 'triangles=',0
      write(2,400) 'shots=',shots
      write(2,400) 'radii=',pore_radii
      write(2,401) 'min_pore_radius=',min_pore_radius
      write(2,401) 'max_pore_radius=',max_pore_radius
      write(2,401) 'mean_pore_radius=',mean_pore_radius
      write(2,401) 'stderr_pore_radius=',stderr_pore_radius
      write(2,401) 'phi_reff=',phi_reff
      write(2,400) 'cells=',1
      write(2,400) 'threads=',1
      write(2,401) 'voro_cpu_time=',0.D0
      write(2,401) 'voro_real_time=',0.D0
      write(2,401) 'triangles_setup_cpu_time=',0.D0
      write(2,401) 'triangles_setup_real-time=',0.D0
      write(2,401) 'MonteCarlo_cpu_time=',cpu_and_real_time(1,1)
      write(2,401) 'MonteCarlo_real_time=',cpu_and_real_time(1,2)
      close(2)

      open(2,file='radii.inf'); rewind(2)
      write(2,*) N
      write(2,*) ro
      write(2,*) rc
      write(2,*) rp
      write(2,*) box(1)*box(2)*box(3)
      write(2,*) 0
      write(2,*) shots
      write(2,*) pore_radii
      write(2,*) min_pore_radius
      write(2,*) max_pore_radius
      write(2,*) mean_pore_radius
      write(2,*) stderr_pore_radius
      write(2,*) phi_reff
      write(2,*) 1
      write(2,*) 1
      write(2,*) 0.D0
      write(2,*) 0.D0
      write(2,*) 0.D0
      write(2,*) 0.D0
      write(2,*) cpu_and_real_time(1,1)
      write(2,*) cpu_and_real_time(1,2)
      close(2)
     endif
     
100 format(A40,3(I10,1x))
101 format(A40,3(F10.4,1x))
400 format(A,I12)
401 format(A,F12.5)
     stop
     end

     subroutine f(R1, ndim, c, grad, need_gradient, f_data)
     use shared
     implicit none
     integer :: ndim
     real*8  :: c(ndim), grad(ndim)
     real*8  :: f_data,R1
     integer :: need_gradient
         if (need_gradient.ne.0) stop 'ERROR. Nonlinear optimization failed. Contact the author.'
         R1 = smallest_coated_surface_distance_allspheres_point(c)-rp
     end

     subroutine fc(val, ndim, c, grad, need_gradient, p)
     use shared
     implicit none
     integer :: ndim
     integer :: need_gradient
     real*8  :: val, c(ndim), grad(ndim), p(dim)
     real*8  :: R1,R2
     real*8  :: pc(dim)
      if (need_gradient.ne.0) stop 'ERROR. Nonlinear optimization failed. contact the author.'
      R1 = smallest_coated_surface_distance_allspheres_point(c)-rp
      pc = p-c
      pc = pc-box*anint(pc/box)
      R2 = dsqrt(dot_product(pc,pc))
      val = R2-R1     ! R2-R1 <= 0 
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
 

    subroutine check_result(p,c) ! not in use
    use shared
    implicit none
    real*8   :: p(dim),c(dim)
    integer  :: id
    real*8   :: R1,R2,singleR1,cp(dim),cx(dim)
       cp = c-p
       cp = cp-box*anint(cp/box)
       R2 = dsqrt(dot_product(cp,cp))
       ! get R1
       R1 = 1.D30
       do id=1,N
          cx = c-x(id,:); cx=cx-box*anint(cx/box)
          singleR1 = dsqrt(dot_product(cx,cx)) - radius(id)-rc-rp
          if (singleR1.lt.R1) then
             R1 = singleR1
          endif
       enddo
       if (R1.lt.0.D0)        stop 'R1 < 0 ERROR. Contact the author'
       if (R1-R2.lt.-1.D-6)   stop 'R1 > R2 ERROR. Contact the author'
    return
    end

    subroutine manage_cpu_and_real_clock(no,cpu_and_real_time)
    use shared
    integer :: no,realcount,countrate
    real    :: realtime,cputime
    real    :: cpu_and_real_time(2,2)
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

