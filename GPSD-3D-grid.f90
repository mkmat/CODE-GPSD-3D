module shared
implicit none

! (c) 27 sep 2023 mk@mat.ethz.ch
! grid-based GPSD 

! routines offered by this code: 

! subroutine read_parameters  ! reads simulation parameters
! subroutine read_config      ! reads the configuration files
! subroutine allocate_grids   ! allocates memory
! subroutine print_something  ! returns some numbers to stdout
! subroutine init_seed        ! initialize seed
! subroutine finalize_seed    ! finalize seed
! subroutine make_radial_list ! creates a list of voxel contained in a sphere
! recursive subroutine quicksort(left,right) ! sort routine
! subroutine exchange(i,j)    ! exchanges entries i and  j
! subroutine fill_grid        ! fills the voxel grid using material sphere data
! subroutine manage_cpu_and_real_clock(no)
! subroutine create_maps      ! prevents mod in other places
! subroutine create_ID2_bottleneck ! creates the Euclidean distance map
! subroutine lookup_ID2       ! looks up the Euclidean distance map
! subroutine get_jp           ! creates a counter that walks spiral-wise around a voxel
! subroutine create_info      ! summarizes the results and creates files


    integer     :: stdout
    integer     :: maxGradius
    integer     :: N,M(3),i,j,k,bin(3),i1,i2,i3,j1,j2,j3
    integer(kind=8) :: maxdata,maxvoxels,voxels,voxels1,jp
    real*8      :: lo(3),hi(3),box(3),boxhalf(3),max_length,vol,pore_radius,min_pore_radius,max_pore_radius,mean_pore_radius,out(3)
    real*8      :: delta,min_delta,rp,rc,ro,volume_fraction
    integer, dimension(:,:), allocatable    :: data 
    integer, dimension(:,:,:), allocatable  :: iD2
    logical, dimension(:,:,:), allocatable  :: grid
    real*8, dimension(:,:), allocatable     :: x
    real*8, dimension(:), allocatable       :: radius
    integer, dimension(:), allocatable      :: map1,map2,map3
    character(len=150), parameter  :: form100='("[GPSD-3D-grid]",A40,3(I12,1x))'
    character(len=150), parameter  :: form101='("[GPSD-3D-grid]",A40,3(F12.3,1x))'
    real    :: cpu_and_real_time(7,2)
    logical :: use_radii_from_config
    logical :: quiet,more

    ! ------------------------------------------------ parameters
    namelist /list/ N,min_delta,maxvoxels,rp,rc,ro,use_radii_from_config,more,quiet
    ! ------------------------------------------------

    contains 

    subroutine read_parameters
    logical     :: EX
        use_radii_from_config = .true.          ! default
        more = .false.   ! default
        inquire(file='.parameters-grid',exist=EX)
        if (.not.EX) stop 'missing .parameters-grid'
        open(2,file='.parameters-grid'); rewind(2)
        read(2,nml=list)
        close(2)
        if (quiet) then
            stdout = 9
            open(stdout,file='log-GPSD-3D-grid.txt'); rewind(9)
        else
            stdout = 6
        endif
    return
    end

    subroutine read_config
    integer :: id,dummy
        write(stdout,form100) 'reading config-with-radii.txt and box.txt'
        allocate(x(N,3))
        allocate(radius(N))
        open(2,file='config-with-radii.txt'); rewind(2)
        if (use_radii_from_config) then
            ro = 0.D0
            do id=1,N; read(2,*) dummy,x(id,:),radius(id); enddo
            do id=1,N; ro = ro + radius(id)/dble(N); enddo
        else
            do id=1,N; read(2,*) dummy,x(id,:); radius(id)=ro; enddo 
        endif
        close(2)
        open(2,file='box.txt'); rewind(2)
        read(2,*) lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)
        close(2)
        box     = hi-lo
        boxhalf = box/2.D0
        vol     = box(1)*box(2)*box(3)
        ! fold required by voro++
        do id=1,N
            x(id,:)=x(id,:)-lo-boxhalf
            x(id,:)=x(id,:)-anint(x(id,:)/box)*box
            x(id,:)=x(id,:)+lo+boxhalf
        enddo
        ! overwrite config.txt
        ! open(2,file='config.txt'); rewind(2)
        ! do id=1,N; write(2,'(I6,1x,4(F14.5,1x))') id,x(id,:),radius(id); enddo
        ! close(2)
        ! 
        ! open(2,file='config-no-radius.txt'); rewind(2)
        ! do id=1,N; write(2,'(I6,1x,3(F14.5,1x))') id,x(id,:); enddo
        ! close(2)
    return
    end

    subroutine allocate_grids
        write(stdout,form100) 'allocating grids'
        delta       = (vol/maxvoxels)**(1.0/3.0)
        delta       = max(min_delta,delta)
        M           = 1+int(box/delta)
        voxels      = M(1)*M(2)*M(3)
        maxGradius  = max(M(1),max(M(2),M(3)))/2      ! 23 sep added /2 here
        if (maxGradius.eq.0) stop 'increase -griddata or decrease -griddelta'
        maxdata     = (2*maxGradius+1)**3
        allocate(iD2(M(1),M(2),M(3)))
        allocate(grid(M(1),M(2),M(3)))
        allocate(data(maxdata,0:3))
    return
    end

    subroutine print_something
        write(stdout,form100) 'N',N
        write(stdout,form101) 'rp',rp
        write(stdout,form101) 'rc',rc
        if (use_radii_from_config) then
            write(stdout,form101) '<ro> [from file]',ro
        else
            write(stdout,form101) 'ro [set]',ro
        endif
        write(stdout,form100) 'maxvoxels',maxvoxels
        write(stdout,form100) 'voxels',voxels
        write(stdout,form101) 'min_delta',min_delta
        write(stdout,form101) 'delta',delta
        write(stdout,form101) 'box',box
        write(stdout,form100) 'M',M
        write(stdout,form101) 'vol',vol
        write(stdout,form100) 'maxGradius [delta units]',maxGradius
    return
    end

    subroutine init_seed
    logical EX
    integer, allocatable :: oldseed(:)
    integer n
        inquire(file='.seed',EXIST=EX)
        call random_seed(size=n)
        allocate(oldseed(n))
        if (EX) then
            open(2,file='.seed'); rewind(2)
            read(2,*) oldseed
            close(2)
            call random_seed(put=oldseed)
        else
            call random_seed()
        endif
     return
     end

     subroutine finalize_seed
     integer, allocatable :: oldseed(:)
     integer n
        call random_seed(size=n)
        allocate(oldseed(n))
        call random_seed(get=oldseed)
        open(2,file='.seed'); rewind(2)
        write(2,*) oldseed
        close(2)
     return
     end

    subroutine make_radial_list
	integer(kind=8) left,right
	integer(kind=8) j
	integer i0,i1,i2,i3
        ! create data list which contains
        ! pixel values in increasing distance order
        ! generate date: distance to (0,0,0) 
	    j = 0
	    do i1=-maxGradius,maxGradius
	    do i2=-maxGradius,maxGradius
	    do i3=-maxGradius,maxGradius
	        j  = j + 1
	        i0 = i1**2 + i2**2 + i3**2
	        data(j,:) = [i0,i1,i2,i3]
	    enddo
	    enddo
	    enddo
        ! sort
	    left  = 1
	    right = maxdata
	    call quicksort(left,right)
	return
	end

	recursive subroutine quicksort(left,right)
	integer(kind=8) left,right,divisor
	    do while (right.gt.left)
	        divisor = divide(left,right)
	        if (right-divisor.gt.divisor-left) then
	            call quicksort(left,divisor-1); left = divisor + 1
	        else
	            call quicksort(divisor+1,right); right = divisor - 1
	        endif
	    enddo
	return
	end

	function divide(left,right)
	integer(kind=8) :: divide,pivot
    integer(kind=8) :: left,right,i,j
	    i = left; j = right - 1; pivot = data(right,0)
	    do while (i.lt.j) 
	        do while (data(i,0).le.pivot.and.i.lt.right); i=i+1; enddo
	        do while (data(j,0).ge.pivot.and.j.gt.left);  j=j-1; enddo
	        if (i.lt.j) then
	            call exchange(i,j)
	        endif
	    enddo
	    if (data(i,0).gt.pivot) call exchange(i,right)
   	    divide = i
	return
	end

	subroutine exchange(i,j)
	integer(kind=8) i,j,tmp(0:3)
	    tmp       = data(i,:)
        data(i,:) = data(j,:)
        data(j,:) = tmp
	return
	end

    subroutine fill_grid
    integer :: limit
        write(stdout,form100) 'create void grid ..'
        grid=.false.
        do i=1,N
            limit = nint((radius(i)+rc)/delta)**2
            do k=1,3; bin(k)=max(1,min(M(k),1+int((x(i,k)-lo(k))/delta))); enddo
            grid(bin(1),bin(2),bin(3))=.true.
            do j=1,maxdata
                j1 = map1(bin(1) + data(j,1))
                j2 = map2(bin(2) + data(j,2))
                j3 = map3(bin(3) + data(j,3))
                if (data(j,0).gt.limit) goto 1
                grid(j1,j2,j3)=.true.
            enddo
1           continue
        enddo
        deallocate(x)
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


    subroutine create_maps  ! just to prevent mod in other places
    integer :: j
        allocate(map1(-2*M(1):3*M(1)))
        allocate(map2(-2*M(2):3*M(2)))
        allocate(map3(-2*M(3):3*M(3)))
        do j=-2*M(1),3*M(1); map1(j) = mod(j+2*M(1),M(1))+1; enddo 
        do j=-2*M(2),3*M(2); map2(j) = mod(j+2*M(2),M(2))+1; enddo
        do j=-2*M(3),3*M(3); map3(j) = mod(j+2*M(3),M(3))+1; enddo
    return
    end

    subroutine create_ID2_bottleneck
    integer         :: i1,i2,i3
    integer(kind=8) :: jmax
        write(stdout,form100) 'create iD2 on void space ..'
        iD2 = 0
        voxels1 = 0
        do i1=1,M(1)
        do i2=1,M(2)
        do i3=1,M(3)
        if (.not.grid(i1,i2,i3)) then
            ! find distance > rp to nearest occupied grid
            jmax = 1
            do j=jp,maxdata
                j1 = map1(i1 + data(j,1))
                j2 = map2(i2 + data(j,2))
                j3 = map3(i3 + data(j,3))
                if (grid(j1,j2,j3)) goto 2
                jmax=j
            enddo
2           continue
            ! then fill the sphere with this jmax value
            ! if it exceeds the jp threshold (see get_jp)
            if (jmax.ge.jp) then
                do j=1,jmax
                    j1 = map1(i1 + data(j,1))
                    j2 = map2(i2 + data(j,2))
                    j3 = map3(i3 + data(j,3))
                    iD2(j1,j2,j3)=max(iD2(j1,j2,j3),jmax)
                enddo
            endif
        else
            voxels1=voxels1+1
        endif
        enddo
        enddo
        enddo
        do i1=1,M(1)
        do i2=1,M(2)
        do i3=1,M(3)
            if (iD2(i1,i2,i3).gt.0) iD2(i1,i2,i3)=data(iD2(i1,i2,i3),0)     ! squared in delta units
        enddo
        enddo
        enddo
        volume_fraction = dble(voxels1)/dble(voxels)
        write(stdout,form100) '1-voxels',voxels1
        write(stdout,form101) 'volume fraction',volume_fraction
        open(98,file='_volume_fraction'); write(98,*) volume_fraction; close(98)
    return
    end

    subroutine lookup_ID2
    integer :: i1,i2,i3
    integer(kind=4) :: count
        write(stdout,form100) 'lookup iD2 +/- delta ..'
        min_pore_radius=1.D30
        max_pore_radius=0.D0
        mean_pore_radius=0.D0
        open(99,file='_radii'); rewind(99)
        count = 0
        do i1=1,M(1)
        do i2=1,M(2)
        do i3=1,M(3)
            if (.not.grid(i1,i2,i3)) then
            if (id2(i1,i2,i3).ne.0) then
                pore_radius = dsqrt(iD2(i1,i2,i3)*delta**2)
                min_pore_radius=min(min_pore_radius,pore_radius)
                max_pore_radius=max(max_pore_radius,pore_radius)
                mean_pore_radius=mean_pore_radius+pore_radius
                count = count + 1
                if (more) then
                    out(1)=lo(1)+(i1-0.5)*delta
                    out(2)=lo(2)+(i2-0.5)*delta
                    out(3)=lo(3)+(i3-0.5)*delta
                    do k=1,3; out(k)=out(k)-anint((out(k)-(lo(k)+hi(k))/2.D0)/box(k))*box(k); enddo
                    write(99,'(4(F14.6,1x))') pore_radius,out
                else
                    write(99,*) pore_radius
                endif
            endif
            endif
        enddo
        enddo
        enddo
        mean_pore_radius=mean_pore_radius/dble(count)
        write(stdout,form101) 'min pore radius',min_pore_radius
        write(stdout,form101) 'mean pore radius',mean_pore_radius 
        write(stdout,form101) 'max pore radius',max_pore_radius
        close(99)
        open(99,file='_max_pore_radius'); rewind(99); write(99,*) max_pore_radius; close(99)
        write(stdout,form100) 'created _radii (crude result)'
        write(stdout,form100) 'created _volume_fraction'
        write(stdout,form100) 'created _max_pore_radius'
        write(stdout,form100) 'NOT creating _radii_fine'
        if (more) write(stdout,form100) 'created _radii'
        ! call system('./voxels-poreradii_2021_postprocessing.ex')
        ! call system('mv _radii_fine _radii')
    return
    end

    subroutine get_jp
    integer :: ratio
        ratio = (rp/delta)**2
        jp = 1
        do while (data(jp,0).lt.ratio)
            jp = jp + 1
        enddo
        write(stdout,form100) 'jp [d-unit]',jp
        write(stdout,form101) 'grid rp',delta*dsqrt(dble(data(jp,0)))
        write(stdout,form100) 'maxdata',maxdata
    return
    end

    subroutine create_info
    integer     :: stage
        write(stdout,form101) 'overhead [cpu+real secs]',cpu_and_real_time(1,:)
        write(stdout,form101) 'maps [cpu+real secs]',cpu_and_real_time(2,:)
        write(stdout,form101) 'radial lists [cpu+real secs]',cpu_and_real_time(3,:)
        write(stdout,form101) 'fill grid [cpu+real secs]',cpu_and_real_time(4,:)
        write(stdout,form101) 'create ID2 [cpu+real secs]',cpu_and_real_time(5,:)
        write(stdout,form101) 'lookup ID2 [cpu+real secs]',cpu_and_real_time(6,:)
        open(22,file='_info-grid.txt'); rewind(22)
        write(22,*) N
        write(22,*) ro
        write(22,*) rp
        write(22,*) rc
        write(22,*) box(1)*box(2)*box(3)
        write(22,*) M(1)*M(2)*M(3)
        write(22,*) mean_pore_radius
        write(22,*) max_pore_radius
        write(22,*) volume_fraction  ! phi(reff)
        write(22,*) delta
        do stage=1,6
        write(22,*) cpu_and_real_time(stage,1)
        write(22,*) cpu_and_real_time(stage,2)
        enddo
        close(22)
    return
    end

end module shared

program GPSD_3D_grid
use shared

    call manage_cpu_and_real_clock(1)
    call init_seed
    call read_parameters
	call read_config
    call allocate_grids
    call manage_cpu_and_real_clock(2)
    call create_maps
    call print_something

	write(stdout,form100) 'creating radial list ONCE'
    call manage_cpu_and_real_clock(3)
	call make_radial_list
    call get_jp
    call manage_cpu_and_real_clock(4)
    call fill_grid
    call manage_cpu_and_real_clock(5)
    call create_ID2_bottleneck
    call manage_cpu_and_real_clock(6)
    call lookup_ID2
    call manage_cpu_and_real_clock(7)
    call create_info
	call finalize_seed

100	format(A50,3(I10,1x))
101	format(A50,3(F10.4,1x))
	stop
	end
