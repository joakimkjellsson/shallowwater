!!
!! A simple ocean model
!! Based on one-layer, linearised shallow water equations
!! Written in Fortran 90
!!
!!
!! Compile with typing "make" after you tuned the Makefile
!!
!!
!! AUTHORS: Joakim Kjellsson and Laurent Brodeau
!!
!!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!! HISTORY
!!
!! - sept. 2010, Brodeau: added netcdf and namelist support
!!
!! - 2010, original code by Joakim Kjellsson
!!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!
PROGRAM SHALLOW
  !!
  !! This is a preprocessor directive
  !! When compiling, the compiler will first read these and evaluate the 
  !! if statement. If with_netcdf is an argument to the compiler then 
  !! this USE statement will be included in the file. 
  !! That way, we can control whether to compile the model with 
  !! netcdf or not. 
#if defined with_netcdf
  USE mod_write_nc
#endif
  !!
  !!
  !! If a variable is not decleared (see below), then Fortran will 
  !! assume variables starting with i,j,k,l,m,n are integers
  !! and all others are floats. 
  !! This is often not desired, so all Fortran code should say 
  !! IMPLICIT none 
  !! which will cause the program to stop if it encounters an 
  !! undeclared variable rather than assuming its form. 
  IMPLICIT none 
  !!
  !! Now we have to declare variables
  !! Variables can be floats, integers, logicals (true/false)
  !! characters (strings) etc. 
  !!
  !!
  !! Starting with parameters
  !!
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !! Parameters which are given a value in the namelist:
  !! - setting a default value anyway
  !! => change these parameters in the namelist file, NOT HERE!
  !! => you can save your different namelist files...
  !!
  REAL ::          &
       & D0    = 1000.,       &        !average depth
       & f0    = 1e-4,        &        !Coriolis constant
       & g     = 9.81,        &        !gravity
       & gamma = 0.1,         &        !Asselin coefficient
       & Lx    = 5.*1e+7,     &        !width of domain [m]
       & Ly    = 5.*1e+7,     &
       & Ldx   = 1.*1e+6,     &
       & Ldy   = 1.*1e+6,     &
       & tm    = 60*60*24*100.,&       !length of run [s]
       & alfa  = 0.,          &        !slope of bottom   (dD/dy)
       & beta  = 0.                    !slope of Coriolis (df/dy)
  !!
  INTEGER :: &
       & imt   = 101, &         !number of x-points
       & jmt   = 101, &         !number of y-points
       & nsubcycles = 100
  !!
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!
  !! Now we declare parameters which are not modified by the namelist
  LOGICAL :: lexist
  !!
  !! Size of grid cells and time step is determined by the domain size,
  !! number of time steps etc.
  REAL :: dx, dy, dt
  !!
  !! We need to define some physical constants too
  !! Most models will define pi, universal gas constant, etc
  REAL :: pi = 3.1415
  !!
  !!
  !!
  ! Define arrays
  !! u = zonal velocity 
  !! v = meridional velocity 
  !! h = surface height
  !! Also define temporary versions of them
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: u, v, h, u_tmp, v_tmp, h_tmp
  !! spongeu, spongev, spongeh = 
  !! relaxes all variables to zero along boundaries to avoid reflecting waves
  !! f = Coriolis
  !! D = depth
  REAL, DIMENSION(:,:),   ALLOCATABLE :: spongeu, spongev,spongeh,f,D
  !! x, y and time arrays
  REAL, DIMENSION(:),     ALLOCATABLE :: vx, vy, vtime
  !!
  !! Indices for loops
  INTEGER :: nt, i, j, k, l, m, n, ic, im, ip, jc, jp, jm, nc, nm, np
  !!
  !! Tendencies of u,v,h each step, i.e. the new step is found by e.g.
  !! h(t+1) = h(t) + dh * dt (Euler step)
  !! or 
  !! h(t+1) = h(t-1) + dh * 2*dt (Leap frog)
  REAL    :: du,dv, dh
  !!
  !! Defining namelist section "nsetup":
  NAMELIST /nsetup/ D0, f0, beta, alfa, g, gamma, imt, jmt, Lx, Ly, Ldx, Ldy, tm, nsubcycles
  !!
  !!
  !!
  !! Declarations are done, starting the program
  !!
  !! First testing if the namelist file is present in the current directory:
  INQUIRE(FILE='namelist', EXIST=lexist )
  IF ( .NOT. lexist ) THEN
     PRINT *, 'ERROR: file "namelist" not found!'; STOP
  END IF
  !!
  !! Opening and reading the namelist
  !! This will overwrite the defaults from the declarations above
  OPEN( UNIT=11, FILE='namelist', FORM='FORMATTED', STATUS='OLD' )
  READ(11,nsetup)
  CLOSE(11)
  !!
  !!
  !! Compute grid cell sizes
  dx = Lx/REAL(imt-1)  !delta x
  dy = Ly/REAL(jmt-1)  !delta y
  !!
  !!
  !! Time step length is determined by grid box size,
  !! a CFL-number, and phase speed (c^2 = g*D)
  !! We choose CFL = 0.1
  !! but larger numbers may be ok. 
  dt = 0.1*min(dx,dy)/sqrt(g*D0)
  !! We add +1 to make sure that the model runs to at least tm
  nt = int((int(tm/dt)+1)/nsubcycles)+1  
  !!
  !!
  !! ====================================================
  !! === Start ===
  !! ====================================================
  print*,' === Shallow water model === '
  print*,' dx = ',dx
  print*,' dy = ',dy
  print*,' dt = ',dt
  !!
  !!
  !! Allocating arrays
  !! Now we allocate all variables in memory
  !!
  !! Note: u,v,h are the arrays where store results each time step
  !!       u_tmp etc are arrays in which we store non-saved steps
  !!       A leap frog scheme only needs three time levels, so 
  !!       arrays only need to have 3 steps
  !!
  ALLOCATE ( u(imt,jmt-1,nt), v(imt-1,jmt,nt), h(imt-1,jmt-1,nt), &
       &     vx(imt), vy(jmt), vtime(nt), &
       &     u_tmp(imt,jmt-1,3), v_tmp(imt-1,jmt,3), h_tmp(imt-1,jmt-1,3),&
       &     spongeu(imt,jmt-1), spongev(imt-1,jmt), spongeh(imt-1,jmt-1),&
       &     f(imt-1,jmt-1), D(imt-1,jmt-1) )
  !!
  !! Fill coordinates vectors:
  DO ic=1, imt
     vx(ic) = (ic - 1)*dx
  END DO
  DO jc=1, jmt
     vy(jc) = (jc - 1)*dy
  END DO
  DO n=1, nt
     vtime(n) = (n - 1)*dt*float(nsubcycles)
  END DO
  print*,' tmin = ',vtime(1)
  print*,' tmax (days) = ',vtime(nt)/(3600.0*24.0)
  !!
  !!
  ! ====================================================
  ! === Reset matrices and variables ===
  ! ====================================================
  !!
  !! It is good practice to set all arrays to zero before starting
  du = 0.
  dv = 0.
  dh = 0.
  u = 0.
  v = 0.
  h = 0.
  u_tmp = 0.
  v_tmp = 0.
  h_tmp = 0.
  D = 0.
  f = 0.
  !!
  !! ====================================================
  !! === Define topography and beta plane ===
  !! ====================================================
  !!
  D = D0
  f = f0
  do j=1,jmt-1
       f(:,j) = f0 + beta*(vy(j)-Ly/2)
       D(:,j) = D0 + alfa*(vy(j))
  enddo
  print*,' fmax, fmin = ',maxval(f),minval(f)
  print*,' Dmax, Dmin = ',maxval(D),minval(D)
  !!
  !! ====================================================
  !! === Initial condition ===
  !! ====================================================
  !!
  ! Initial Gaussian disturbance
  DO i=1,imt-1
       DO j=1,jmt-1
            h(i,j,1) = 2.*exp(-( ((vx(i)-Lx/2)/Ldx)**2 + ((vy(j)-Ly/2)/Ldy)**2 ))
       ENDDO
  ENDDO
  
  ! In geostrophic balance
  DO ic=2,imt-2
       DO jc=2,JMT-2
            ip = ic+1 ; im = ic-1
            jp = jc+1 ; jm = jc-1
            
            u(ic,jc,1) = -g/(0.5*(f(ic,jc)+f(im,jc)))*0.25/dy*&
            (h(ic,jp,1)+h(im,jp,1)-h(ic,jm,1)-h(im,jm,1))
            v(ic,jc,1) =  g/(0.5*(f(ic,jc)+f(ic,jm)))*0.25/dx*&
            (h(ip,jc,1)+h(ip,jm,1)-h(im,jc,1)-h(im,jm,1))
       ENDDO
  ENDDO
  
  !!
  !!
  !! ====================================================
  !! === Sponge zone ===
  !! ====================================================
  !!
  spongeu = 0.
  spongev = 0.
  spongeh = 0.
  
  !! Define a sponge zone which is 1 in most of the domain
  !! but approaches zero near boundaries. 
  !! Then we multiply the velocities etc with the sponge
  
  ! Only sponge at y boundaries
  DO j=1,11
       spongev(:,j) =       0.5+0.5*cos(pi*(j-1)/10.)
       spongeh(:,j) =       0.5+0.5*cos(pi*(j-1)/10.)
       spongev(:,JMT-j+1) = 0.5+0.5*cos(pi*(j-1)/10.)
       spongeh(:,JMT-j) =   0.5+0.5*cos(pi*(j-1)/10.)
       spongeu(:,j) =       0.5+0.5*cos(pi*(j-1)/10.)
       spongeu(:,JMT-j) =   0.5+0.5*cos(pi*(j-1)/10.)
  ENDDO
  
  spongeu = dble(1) - spongeu
  spongev = dble(1) - spongev
  spongeh = dble(1) - spongeh
  
  !! ====================================================
  !! === Solving with Euler forward on a C-grid ===
  !! ====================================================
  !!
  !!
  !! Loop from 2nd index to 2nd last index
  DO jc=2,jmt-1
     
     jm = jc-1 ; jp = jc+1
     
     DO ic=1,imt-1
        
        im = ic-1 ; ip = ic+1
        !! Periodic boundaries in the zonal direction
        !! i = 0 is equal to last index
        if(im == 0)   im = IMT-1
        if(ip == IMT) ip = 1
        
        !! Reset tendencies
        du = 0.
        dv = 0.
        dh = 0.
        
        ! Height gradient
        du = du - g/dx*( h(ic,jc,1)-h(im,jc,1) )
        dv = dv - g/dy*( h(ic,jc,1)-h(ic,jm,1) )
        ! Coriolis
        du = du + (f(ic,jc)+f(im,jc))/8.*&
                  ( v(ic,jc,1)+v(im,jc,1)+v(ic,jp,1)+v(im,jp,1) )
        dv = dv - (f(ic,jc)+f(ic,jm))/8.*&
                  ( u(ic,jc,1)+u(ic,jm,1)+u(ip,jc,1)+u(ip,jm,1) )
        ! Convergence/divergence
        dh = dh - D(ic,jc)*((u(ip,jc,1)-u(ic,jc,1))/dx + (v(ic,jp,1)-v(ic,jc,1))/dy)&
	              - 0.25/dx*(u(ic,jc,1)+u(ip,jc,1))*(D(ip,jc)-D(im,jc))
	    ! Special case at northern boundary. Southern boundary solved further down. 
	    if(jc.lt.jmt-1) dh = dh - 0.25/dy*(v(ic,jc,1)+v(ic,jp,1))*(D(ic,jp)-D(ic,jm))
        ! Time step
        u(ic,jc,2) = (u(ic,jc,1) + du*dt)*spongeu(ic,jc)
        v(ic,jc,2) = (v(ic,jc,1) + dv*dt)*spongev(ic,jc)
        h(ic,jc,2) = (h(ic,jc,1) + dh*dt)*spongeh(ic,jc)
     END DO
  END DO
  
  !! === With periodic boundary conditions in x,
  !! === there is only one boundary left
  
  ! Calculate h(2:imt-1,1) at southern boundary
  jc = 1 ; jp = 2
  
  DO ic=1,imt-1
     
     im = ic-1 ; ip = ic+1
     if(im == 0) im = IMT-1
     if(ip == IMT) ip = 1
     
     du = 0.
     dh = 0.
     
     ! Height gradient
     du = du - g/dx*( h(ic,jc,1)-h(im,jc,1) )
     ! Coriolis
     du = du + (f(ic,jc)+f(im,jc))/8.*&
               ( v(ic,jc,1)+v(im,jc,1)+v(ic,jp,1)+v(im,jp,1) )
     ! Convergence/divergence
     dh = dh - D(ic,jc)*((u(ip,jc,1)-u(ic,jc,1))/dx + (v(ic,jp,1)-v(ic,jc,1))/dy)&
	           - 0.25/dx*(u(ic,jc,1)+u(ip,jc,1))*(D(ip,jc)-D(im,jc))
     ! Time step
     u(ic,jc,2) = (u(ic,jc,1) + du*dt)*spongeu(ic,jc)
     h(ic,jc,2) = (h(ic,jc,1) + dh*dt)*spongeh(ic,jc)
     
  END DO
  
  
  ! ====================================================
  
  print*,'Done with first Euler step'
  
  u_tmp(:,:,1) = u(:,:,1)
  u_tmp(:,:,2) = u(:,:,2)
  v_tmp(:,:,1) = v(:,:,1)
  v_tmp(:,:,2) = v(:,:,2)
  h_tmp(:,:,1) = h(:,:,1)
  h_tmp(:,:,2) = h(:,:,2)
  
  ! ====================================================
  ! === Solving with Leap frog on a C-grid ===
  ! ====================================================
  
  ! Time loops
  DO n=3,nt
     
     ! We only save data every full step
     ! Then we evolve the model between full steps
     ! Another way to do it would be to run only full steps
     ! but write output every 10 steps or so. 
     DO m=1,nsubcycles
        
        ! Space loops
        DO jc=2,jmt-1
           
           jp = jc+1 ; jm = jc-1
           
           DO ic=1,imt-1
              
              im = ic-1 ; ip = ic+1
              if(im == 0) im = IMT-1
              if(ip == IMT) ip = 1
              
              du = 0.
              dv = 0.
              dh = 0.
              
              ! Height gradient
              du = du - g/dx*( h_tmp(ic,jc,2)-h_tmp(im,jc,2) )
              dv = dv - g/dy*( h_tmp(ic,jc,2)-h_tmp(ic,jm,2) )
              ! Coriolis
              du = du + (f(ic,jc)+f(im,jc))/8.*&
                              ( v_tmp(ic,jc,2)+v_tmp(im,jc,2)&
                               +v_tmp(ic,jp,2)+v_tmp(im,jp,2) )
              dv = dv - (f(ic,jc)+f(ic,jm))/8.*&
                              ( u_tmp(ic,jc,2)+u_tmp(ip,jc,2)&
                               +u_tmp(ic,jm,2)+u_tmp(ip,jm,2) )          
              ! Convergence/divergence
              dh = dh - D(ic,jc)*&
                         ((u_tmp(ip,jc,2)-u_tmp(ic,jc,2))/dx+&
                          (v_tmp(ic,jp,2)-v_tmp(ic,jc,2))/dy)&
                - 0.25/dx*(u_tmp(ic,jc,2)+u_tmp(ip,jc,2))*(D(ip,jc)-D(im,jc))
              if(jc.lt.jmt-1) dh = dh - 0.25/dy*&
                          (v_tmp(ic,jc,2)+v_tmp(ic,jp,2))*(D(ic,jp)-D(ic,jm))
              ! Time step
              u_tmp(ic,jc,3) = u_tmp(ic,jc,1) + du*2.*dt
              v_tmp(ic,jc,3) = v_tmp(ic,jc,1) + dv*2.*dt
              h_tmp(ic,jc,3) = h_tmp(ic,jc,1) + dh*2.*dt
              ! Asselin filter
              u_tmp(ic,jc,2) = u_tmp(ic,jc,2) + gamma * &
                   (u_tmp(ic,jc,1)+u_tmp(ic,jc,3)-2*u_tmp(ic,jc,2))
              v_tmp(ic,jc,2) = v_tmp(ic,jc,2) + gamma * &
                   (v_tmp(ic,jc,1)+v_tmp(ic,jc,3)-2*v_tmp(ic,jc,2))
              h_tmp(ic,jc,2) = h_tmp(ic,jc,2) + gamma * &
                   (h_tmp(ic,jc,1)+h_tmp(ic,jc,3)-2*h_tmp(ic,jc,2))
              ! Sponge
              u_tmp(ic,jc,3) = u_tmp(ic,jc,3)*spongeu(ic,jc)
              v_tmp(ic,jc,3) = v_tmp(ic,jc,3)*spongev(ic,jc)
              h_tmp(ic,jc,3) = h_tmp(ic,jc,3)*spongeh(ic,jc)
           END DO
        END DO
        
        ! === Boundary ===
                
        ! === Calculate h(2:imt-1,1) ===
        
        jc = 1 ; jp = 2
        
        DO ic=1,imt-1
           
           im = ic-1 ; ip = ic+1
           if(im == 0) im = IMT-1
           if(ip == IMT) ip = 1
           
           du = 0.
           dh = 0.
           
           ! Height gradient
           du = du - g/dx*( h_tmp(ic,jc,2)-h_tmp(im,jc,2) )
           ! Coriolis
           du = du + (f(ic,jc)+f(im,jc))/8.*&
                           ( v_tmp(im,jp,2)+v_tmp(ic,jp,2)&
                            +v_tmp(im,jc,2)+v_tmp(ic,jc,2) )                
           ! Convergence/divergence
           dh = dh - D(ic,jc)*&
                            ((u_tmp(ip,jc,2)-u_tmp(ic,jc,2))/dx+&
                             (v_tmp(ic,jp,2)-v_tmp(ic,jc,2))/dy)&
                   - 0.25/dx*(u_tmp(ic,jc,2)+u_tmp(ip,jc,2))*(D(ip,jc)-D(im,jc))
           ! Time step
           u_tmp(ic,jc,3) = u_tmp(ic,jc,1) + du*2.*dt
           h_tmp(ic,jc,3) = h_tmp(ic,jc,1) + dh*2.*dt
           ! Asselin filter
           u_tmp(ic,jc,2) = u_tmp(ic,jc,2) + gamma * &
                (u_tmp(ic,jc,1)+u_tmp(ic,jc,3)-2*u_tmp(ic,jc,2))    
           h_tmp(ic,jc,2) = h_tmp(ic,jc,2) + gamma * &
                (h_tmp(ic,jc,1)+h_tmp(ic,jc,3)-2*h_tmp(ic,jc,2))
           ! Sponge
           u_tmp(ic,jc,3) = u_tmp(ic,jc,3)*spongeu(ic,jc)
           h_tmp(ic,jc,3) = h_tmp(ic,jc,3)*spongeh(ic,jc)
        END DO
        
        
        ! Re-arrage matrices
        ! Move level 2 to level 1
        !  and level 3 to level 2
        ! Then reset level 3
        u_tmp(:,:,1) = u_tmp(:,:,2)
        v_tmp(:,:,1) = v_tmp(:,:,2)
        h_tmp(:,:,1) = h_tmp(:,:,2)
        u_tmp(:,:,2) = u_tmp(:,:,3)
        v_tmp(:,:,2) = v_tmp(:,:,3)
        h_tmp(:,:,2) = h_tmp(:,:,3)
        u_tmp(:,:,3) = 0.
        v_tmp(:,:,3) = 0.
        h_tmp(:,:,3) = 0.
        
     END DO
     
     ! Save to u,v,h matrices
     u(:,:,n) = u_tmp(:,:,2)
     v(:,:,n) = v_tmp(:,:,2)
     h(:,:,n) = h_tmp(:,:,2)
     print*,'Done with n = ',n,' of ',nt,' timesteps'
  END DO
  
  
  ! Write data to a file

#if defined with_netcdf
  !!
  CALL WRITE_NC(imt,jmt-1,nt, vx, (vy(1:jmt-1)+dy/2), vtime, u, 'u_test.nc', 'u')
  CALL WRITE_NC(imt-1,jmt,nt, (vx(1:imt-1)+dx/2), vy, vtime, v, 'v_test.nc', 'v')
  CALL WRITE_NC(imt-1,jmt-1,nt, (vx(1:imt-1)+dx/2),&
               (vy(1:jmt-1)+dy/2), vtime, h, 'h_test.nc', 'h')
  !!
#else
  !! File must be opened, written to, then closed.
  !! It must be given an id (unit), and the format must be specified.
  !! 'unformatted' means that the file will be binary.
  open(unit=111,file='u_field.bin',form='unformatted')
  write(unit=111) imt,jmt-1,nt,dx,dy,dt*dble(nsubcycles),&
                  vx,(vy(1:jmt-1)+dy/2), vtime,u
  close(unit=111)
  !!
  open(unit=222,file='v_field.bin',form='unformatted')
  write(unit=222) imt-1,jmt,nt,dx,dy,dt*dble(nsubcycles),&
                  (vx(1:imt-1)+dx/2),vy, vtime,v
  close(unit=222)
  !!
  open(unit=333,file='h_field.bin',form='unformatted')
  write(unit=333) imt-1,jmt-1,nt,dx,dy,dt*dble(nsubcycles),&
                  (vx(1:imt-1)+dx/2),(vy(1:jmt-1)+dy/2), vtime,h
  close(unit=333)
  !!
#endif
  !!
  !!

  

  ! === THE END ===
  
  
  !!
END PROGRAM SHALLOW
