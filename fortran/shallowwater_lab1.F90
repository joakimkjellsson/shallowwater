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
#if defined with_netcdf
  USE mod_write_nc
#endif
  !!
  !!
  IMPLICIT none ! Perhaps the most important line in the code
  !!
  !! First, we have to declare variables
  !! Starting with parameters
  !!
  !!
  !!
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !! Parameters which are given a value in the namelist:
  !! - setting a default value anyway
  !! => change these parameters in the namelist file, NOT HERE!
  !! => you can save your different namelist files...
  !!
  REAL ::          &
       & D     = 4000.,       &        !average depth
       & f0    = 1e-4,        &        !Coriolis constant
       & g     = 9.81,        &        !gravity
       & gamma = 0.1,         &        !Asselin coefficient
       & Lx    = 5.*1e+7,     &        !width of domain [m]
       & Ly    = 5.*1e+7,     &
       & tm    = 60*60*24*10.          !length of run [s]
  !!
  INTEGER :: &
       & imt   = 101, &         !number of x-points
       & jmt   = 101, &         !number of y-points
       & nsubcycles = 50
  !!
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!
  LOGICAL :: lexist
  !!
  REAL :: dx, dy, dt
  !!
  !!
  !!
  ! Define matrices
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: u, v, h, u_tmp, v_tmp, h_tmp
  REAL, DIMENSION(:),     ALLOCATABLE :: vx, vy, vtime
  !!
  !! Indices for loops
  INTEGER :: nt, i, j, k, l, m, n, ic, im, ip, jc, jp, jm, nc, nm, np
  REAL    :: du,dv, dh
  !!
  !! Defining namelist section "nsetup":
  NAMELIST /nsetup/ D, f0, g, gamma, imt, jmt, Lx, Ly, tm, nsubcycles
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
  !! Opening and reading the namelist:
  OPEN( UNIT=11, FILE='namelist', FORM='FORMATTED', STATUS='OLD' )
  READ(11,nsetup)
  CLOSE(11)
  !!
  !!
  dx = Lx/REAL(imt-1)  !delta x
  dy = Ly/REAL(jmt-1)  !delta y
  !!
  !!
  !! Time step length is determined by grid box size,
  !! a CFL-number, and phase speed (c^2 = g*D)
  dt = 0.1*min(dx,dy)/sqrt(g*D)
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
  !! Allocating arrays:
  ALLOCATE ( u(imt,jmt-1,nt), v(imt-1,jmt,nt), h(imt-1,jmt-1,nt), &
       &     vx(imt), vy(jmt), vtime(nt), &
       &     u_tmp(imt,jmt-1,3), v_tmp(imt-1,jmt,3), h_tmp(imt-1,jmt-1,3) )
  !!
  !! Building coordinates vectors:
  DO ic=1, imt
     vx(ic) = (ic - 1)*dx
  END DO
  DO jc=1, jmt
     vy(jc) = (jc - 1)*dy
  END DO
  DO n=1, nt
     vtime(n) = (n - 1)*dt
  END DO
  !!
  !!
  ! ====================================================
  ! === Reset matrices and variables ===
  ! ====================================================
  !!
  du = 0.
  dv = 0.
  dh = 0.
  u = 0.
  v = 0.
  h = 0.
  u_tmp = 0.
  v_tmp = 0.
  h_tmp = 0.
  !!
  !! ====================================================
  !! === Initial condition ===
  !! ====================================================
  
  !!
  !! ====================================================
  !! === Sponge zone ===
  !! ====================================================
  !!
  !! ====================================================
  !! === Solving with Euler forward on a C-grid ===
  !! ====================================================
  !!
  !!
  !!
  DO jc=2,jmt-1
     
     jm = jc-1 ; jp = jc+1
     
     DO ic=2,imt-1
        
        im = ic-1 ; ip = ic+1
        
        du = 0.
        dv = 0.
        dh = 0.
        
        ! Height gradient
        du = du - g/dx*( h(ic,jc,1)-h(im,jc,1) )
        dv = dv - g/dy*( h(ic,jc,1)-h(ic,jm,1) )
        ! Coriolis
        du = du + f0/4.*( v(ic,jc,1)+v(im,jc,1)+v(ic,jp,1)+v(im,jp,1) )
        dv = dv - f0/4.*( u(ic,jc,1)+u(ic,jm,1)+u(ip,jc,1)+u(ip,jm,1) )
        ! Convergence/divergence
        dh = dh - D*((u(ip,jc,1)-u(ic,jc,1))/dx + (v(ic,jp,1)-v(ic,jc,1))/dy)
        ! Time step
        u(ic,jc,2) = u(ic,jc,1) + du*dt
        v(ic,jc,2) = v(ic,jc,1) + dv*dt
        h(ic,jc,2) = h(ic,jc,1) + dh*dt
	
     END DO
  END DO
  
  ! Calculate h(1,2:jmt-1)
  ic = 1 ; ip = 2
  
  DO jc=2,jmt-1
     
     jp = jc+1 ; jm = jc-1	
     
     dh = 0.
     dv = 0.
     
     ! Height gradient
     dv = dv - g/dy*( h(ic,jc,1)-h(ic,jm,1) )
     ! Coriolis
     dv = dv - f0/4.*( u(ic,jc,1)+u(ip,jc,1)+u(ic,jm,1)+u(ip,jm,1) )
     ! Convergence/divergence
     dh = dh - D*((u(ip,jc,1)-u(ic,jc,1))/dx + (v(ic,jp,1)-v(ic,jc,1))/dy)
     ! Time step
     v(ic,jc,2) = v(ic,jc,1) + dv*dt
     h(ic,jc,2) = h(ic,jc,1) + dh*dt
     
  END DO
  
  ! Calculate h(2:imt-1,1)
  jc = 1 ; jp = 2
  
  DO ic=2,imt-1
     
     im = ic-1 ; ip = ic+1
     
     du = 0.
     dh = 0.
     
     ! Height gradient
     du = du - g/dx*( h(ic,jc,1)-h(im,jc,1) )
     ! Coriolis
     du = du + f0/4.*( v(ic,jc,1)+v(im,jc,1)+v(ic,jp,1)+v(im,jp,1) )
     ! Convergence/divergence
     dh = dh - D*((u(ip,jc,1)-u(ic,jc,1))/dx + (v(ic,jp,1)-v(ic,jc,1))/dy)
     ! Time step
     u(ic,jc,2) = u(ic,jc,1) + du*dt
     h(ic,jc,2) = h(ic,jc,1) + dh*dt
     
  END DO
  
  ! Calculate h(1,1)
  jc = 1 ; jp = 2
  ic = 1 ; ip = 2
  
  dh = 0.
  
  ! Convergence/divergence
  dh = dh - D*((u(ip,jc,1)-u(ic,jc,1))/dx + (v(ic,jp,1)-v(ic,jc,1))/dy)
  ! Time step
  h(ic,jc,2) = h(ic,jc,1) + dh*dt
  
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
     
     DO m=1,nsubcycles
        
        ! Space loops
        DO jc=2,jmt-1
           
           jp = jc+1 ; jm = jc-1
           
           DO ic=2,imt-1
              
              im = ic-1 ; ip = ic+1
              
              du = 0.
              dv = 0.
              dh = 0.
              
              ! Height gradient
              du = du - g/dx*( h_tmp(ic,jc,2)-h_tmp(im,jc,2) )
              dv = dv - g/dy*( h_tmp(ic,jc,2)-h_tmp(ic,jm,2) )
              ! Coriolis
              du = du + f0/4.*( v_tmp(ic,jc,2)+v_tmp(im,jc,2)&
                               +v_tmp(ic,jp,2)+v_tmp(im,jp,2) )
              dv = dv - f0/4.*( u_tmp(ic,jc,2)+u_tmp(ip,jc,2)&
                               +u_tmp(ic,jm,2)+u_tmp(ip,jm,2) )          
              ! Convergence/divergence
              dh = dh - D*((u_tmp(ip,jc,2)-u_tmp(ic,jc,2))/dx&
                          +(v_tmp(ic,jp,2)-v_tmp(ic,jc,2))/dy)
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
              
           END DO
        END DO
        
        ! === Boundaries ===
        
        ! Calculate h(1,2:jmt-1)
        ic = 1 ; ip = 2
        
        DO jc=2,jmt-1
           
           jm = jc-1 ; jp = jc+1
           
           dv = 0.
           dh = 0.
           
           ! Height gradient
           dv = dv - g/dy*( h_tmp(ic,jc,2)-h_tmp(ic,jm,2) )
           ! Coriolis
           dv = dv - f0/4.*( u_tmp(ic,jc,2)+u_tmp(ip,jc,2)&
                            +u_tmp(ic,jm,2)+u_tmp(ip,jm,2) )
           ! Convergence/divergence
           dh = dh - D*((u_tmp(ip,jc,2)-u_tmp(ic,jc,2))/dx&
                       +(v_tmp(ic,jp,2)-v_tmp(ic,jc,2))/dy)
           ! Time step
           v_tmp(ic,jc,3) = v_tmp(ic,jc,1) + dv*2.*dt
           h_tmp(ic,jc,3) = h_tmp(ic,jc,1) + dh*2.*dt
           ! Asselin filter
           v_tmp(ic,jc,2) = v_tmp(ic,jc,2) + gamma * &
                (v_tmp(ic,jc,1)+v_tmp(ic,jc,3)-2*v_tmp(ic,jc,2))
           h_tmp(ic,jc,2) = h_tmp(ic,jc,2) + gamma * &
                (h_tmp(ic,jc,1)+h_tmp(ic,jc,3)-2*h_tmp(ic,jc,2))
           
        END DO
        
        ! === Calculate h(2:imt-1,1) ===
        
        jc = 1 ; jp = 2
        
        DO ic=2,imt-1
           
           im = ic-1 ; ip = ic+1
           
           du = 0.
           dh = 0.
           
           ! Height gradient
           du = du - g/dx*( h_tmp(ic,jc,2)-h_tmp(im,jc,2) )
           ! Coriolis
           du = du + f0/4.*( v_tmp(im,jp,2)+v_tmp(ic,jp,2)&
                            +v_tmp(im,jc,2)+v_tmp(ic,jc,2) )                
           ! Convergence/divergence
           dh = dh - D*((u_tmp(ip,jc,2)-u_tmp(ic,jc,2))/dx&
                       +(v_tmp(ic,jp,2)-v_tmp(ic,jc,2))/dy)
           ! Time step
           u_tmp(ic,jc,3) = u_tmp(ic,jc,1) + du*2.*dt
           h_tmp(ic,jc,3) = h_tmp(ic,jc,1) + dh*2.*dt
           ! Asselin filter
           u_tmp(ic,jc,2) = u_tmp(ic,jc,2) + gamma * &
                (u_tmp(ic,jc,1)+u_tmp(ic,jc,3)-2*u_tmp(ic,jc,2))    
           h_tmp(ic,jc,2) = h_tmp(ic,jc,2) + gamma * &
                (h_tmp(ic,jc,1)+h_tmp(ic,jc,3)-2*h_tmp(ic,jc,2))
           
        END DO
        
        ! === Calculate h(1,1)    
        ic = 1 ; ip = 2
        jc = 1 ; jp = 2
        
        dh = 0.
        
        ! Convergence/divergence
        dh = dh - D*((u_tmp(ip,jc,2)-u_tmp(ic,jc,2))/dx&
                    +(v_tmp(ic,jp,2)-v_tmp(ic,jc,2))/dy)
        ! Time step   
        h_tmp(ic,jc,3) = h_tmp(ic,jc,1) + dh*2.*dt
        ! Asselin filter
        h_tmp(ic,jc,2) = h_tmp(ic,jc,2) + gamma * &
             (h_tmp(ic,jc,1)+h_tmp(ic,jc,3)-2*h_tmp(ic,jc,2))
        
        ! Re-arrage matrices
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
  print*,imt,jmt-1,nt,vx(1),vy(1)+dy/2
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
