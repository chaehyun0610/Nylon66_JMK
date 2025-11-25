      PROGRAM POST

      IMPLICIT NONE


      INTEGER :: i, im, ia, ntime, nmol, natom, ndata, dt
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION (:,:,:) :: x, y, z
      INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: ix, iy, iz
      INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: id, mol, type1

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION (:,:,:) :: ux, uy, uz
      REAL,ALLOCATABLE,DIMENSION (:) :: lx, ly, lz
      REAL,ALLOCATABLE,DIMENSION (:) :: xlo, ylo, zlo
      REAL,ALLOCATABLE,DIMENSION (:) :: xhi, yhi, zhi

      INTEGER :: j, k, anum
      INTEGER :: tau, nt, stepcount

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Rx,Ry,Rz,R_mag
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: uRx,uRy,uRz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sumx,sumy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: sumz,sumxy,sumxyz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xx,yy,zz,xy,xyz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tcfxx,tcfyy
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tcfzz,tcfxy,tcf
      
      CHARACTER(len=50) :: txt
      CHARACTER(len=3) :: File_im

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: theta
      DOUBLE PRECISION,PARAMETER :: pi=3.1415926535897d0 
  
      INTEGER :: ib1, ibin1, ib2, ibin2, ib3, ibin3, ndist
      DOUBLE PRECISION :: delr1, delr2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: histo1, histo2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: histo3
 
      INTEGER,PARAMETER :: atype = 9
      REAL(8) :: mass(atype)

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: x_cm,y_cm,z_cm
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: totmass
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Rgsq, Rg
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Rgx,Rgy,Rgz
       

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: PT

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: sum_Rg,sum_Rete
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: sum_Rgxx,sum_Rgxy,
     & sum_Rgyy,sum_Rgzz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: sum_Retexx,sum_Retexy
     & ,sum_Reteyy,sum_Retezz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: avg_Rgxx,avg_Rgxy,
     & avg_Rgyy,avg_Rgzz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: avg_Retexx,avg_Retexy
     & ,avg_Reteyy,avg_Retezz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: avg_Rg,avg_Rete


      ntime = 300
      nmol = 300
      natom = 254
      dt = 5

      ndist = 100
      delr1 = 0.25d0
      delr2 = 1.0d0
  
c      mass(1) = 14
      mass(1) = 12.01115  ! c
      mass(2) = 12.01115  ! c1
      mass(3) = 12.01115  ! c3h
      mass(4) = 12.01115  ! cp
      mass(5) =  1.00797  ! hc
      mass(6) =  1.00797  ! ho
      mass(7) = 15.99940  ! o3e
      mass(8) = 15.99940  ! oc
      mass(9) = 15.99940  ! oh

 
      ALLOCATE (histo1(0:ndist), histo2(-180:180),histo3(0:ndist))
      ALLOCATE (theta(ntime,nmol), R_mag(ntime,nmol))
  
      ALLOCATE (x_cm(ntime,nmol), y_cm(ntime,nmol),z_cm(ntime,nmol))
      ALLOCATE (totmass(ntime,nmol),Rgsq(ntime,nmol),Rg(ntime,nmol))
      ALLOCATE (Rgx(ntime,nmol),Rgy(ntime,nmol),Rgz(ntime,nmol))

      ALLOCATE (sum_Rg(ntime), sum_Rete(ntime))
      ALLOCATE (sum_Rgxx(ntime), sum_Rgxy(ntime), sum_Rgyy(ntime), 
     & sum_Rgzz(ntime))
      ALLOCATE (sum_Retexx(ntime), sum_Retexy(ntime), sum_Reteyy(ntime),
     & sum_Retezz(ntime))
      ALLOCATE (avg_Rg(ntime), avg_Rete(ntime))
      ALLOCATE (avg_Rgxx(ntime), avg_Rgxy(ntime), avg_Rgyy(ntime), 
     & avg_Rgzz(ntime))
      ALLOCATE (avg_Retexx(ntime), avg_Retexy(ntime), avg_Reteyy(ntime),
     &  avg_Retezz(ntime))

      ALLOCATE (PT(ntime,6,2))
  
      ALLOCATE(x(ntime,nmol,natom),y(ntime,nmol,natom),
     & z(ntime,nmol,natom))
      ALLOCATE(ux(ntime,nmol,natom),uy(ntime,nmol,natom),
     & uz(ntime,nmol,natom))
      ALLOCATE(ix(ntime,nmol,natom),iy(ntime,nmol,natom),
     & iz(ntime,nmol,natom))
      ALLOCATE(id(ntime,nmol,natom),mol(ntime,nmol,natom),
     & type1(ntime,nmol,natom))
      ALLOCATE(xlo(ntime),ylo(ntime),zlo(ntime))
      ALLOCATE(xhi(ntime),yhi(ntime),zhi(ntime))
      ALLOCATE(lx(ntime),ly(ntime),lz(ntime))
        
      ALLOCATE (Rx(ntime,nmol),Ry(ntime,nmol),Rz(ntime,nmol))
      ALLOCATE (uRx(ntime,nmol),uRy(ntime,nmol),uRz(ntime,nmol))
      ALLOCATE (sumx(nmol,ntime),sumxy(nmol,ntime))
      ALLOCATE (sumy(nmol,ntime),sumz(nmol,ntime),sumxyz(nmol,ntime))
      ALLOCATE (xx(nmol,ntime),xy(nmol,ntime))
      ALLOCATE (yy(nmol,ntime),zz(nmol,ntime),xyz(nmol,ntime))
      ALLOCATE (tcfxx(ntime),tcfxy(ntime))
      ALLOCATE (tcfyy(ntime),tcfzz(ntime),tcf(ntime))


      R_mag = 0 ; theta = 0.d0; histo1 = 0; histo2 = 0; histo3 = 0
      sum_Rg = 0.d0; sum_Rete = 0.d0; avg_Rg = 0.d0; avg_Rete = 0.d0
      sum_Rgxx = 0.d0; sum_Rgxy = 0.d0; sum_Rgyy = 0.d0; sum_Rgzz = 0.d0
      sum_Retexx = 0.d0; sum_Retexy = 0.d0; sum_Reteyy = 0.d0; 
      sum_Retezz = 0.d0
      avg_Rgxx = 0.d0; avg_Rgxy = 0.d0; avg_Rgyy = 0.d0; avg_Rgzz = 0.d0
      avg_Retexx = 0.d0; avg_Retexy = 0.d0; avg_Reteyy = 0.d0; 
      avg_Retezz = 0.d0

  
      x=0.0; y=0.0; z=0.0; ix=0; iy=0; iz=0; id=0; mol=0; type1=0
      xlo=0.0; xhi=0.0; ylo=0.0; yhi=0.0; zlo=0.0; zhi=0.0
      ux=0.d0; uy=0.d0; uz=0.d0; lx=0.0; ly=0.0; lz=0.0

      Rx=0.d0; Ry=0.d0; Rz=0.d0; uRx=0.d0; uRy=0.d0; uRz=0.d0
      sumx=0.d0; sumy=0.d0; sumxy=0.d0; sumz=0.d0; sumxyz=0.d0
      xx=0.d0; xy=0.d0; yy=0.d0; zz=0.d0; xyz=0.d0
      tcfxx=0.d0; tcfxy=0.d0; tcfyy=0.d0; tcfzz=0.d0; tcf=0.d0


      OPEN (1, file="trj.lammpstrj", form='formatted')
      OPEN (2, file='trj_image.xyz', position='append',
     & form='formatted')

c      OPEN(3,file='Rete.dat',position='append',form='formatted')
      OPEN(4,file='Rete_PDF_0.25.dat',position='append',
     &form='formatted')

      OPEN(5,file='theta.dat',position='append',form='formatted')
      OPEN(6,file='theta_PDF_0.25.dat',position='append',
     &form='formatted')
        
      DO i = 1, ntime

      stepcount = 0; anum = 0
 
        READ (1,*) txt
        READ (1,*) stepcount
        PRINT *, stepcount
        READ (1,*) txt
        READ (1,*) anum
        READ (1,*) txt
        READ (1,*) xlo(i), xhi(i)
        READ (1,*) ylo(i), yhi(i)
        READ (1,*) zlo(i), zhi(i)
        READ (1,*) txt

        lx(i) = xhi(i) - xlo(i)
        ly(i) = yhi(i) - ylo(i)
        lz(i) = zhi(i) - zlo(i)

        DO im = 1, nmol
          DO ia = 1, natom
            READ(1,*) id(i,im,ia), mol(i,im,ia), type1(i,im,ia),
     & x(i,im,ia),y(i,im,ia),z(i,im,ia), 
     & ix(i,im,ia), iy(i,im,ia), iz(i,im,ia)
          END DO
        END DO
      END DO

      DO i = 1, ntime
        WRITE(2,*) anum
        WRITE(2,*)
        DO im = 1, nmol
          DO ia = 1, natom
            ux(i,im,ia) = x(i,im,ia) + ix(i,im,ia)*lx(i)
            uy(i,im,ia) = y(i,im,ia) + iy(i,im,ia)*ly(i)
            uz(i,im,ia) = z(i,im,ia) + iz(i,im,ia)*lz(i)
  
            WRITE(2,*) 'C ', ux(i,im,ia), uy(i,im,ia), uz(i,im,ia)
          END DO      
        END DO
      END DO


!!!   MOVIE for single chain   !!!
      DO i = ntime/2, ntime

        DO im = 1, nmol, 10

          IF (im .lt. 10) THEN
            WRITE(File_im(1:1),'(I1)') 0
            WRITE(File_im(2:2),'(I1)') 0
            WRITE(File_im(3:3),'(I1)') im
          ELSE IF (im .lt. 100) THEN
            WRITE(File_im(1:1),'(I1)') 0
            WRITE(File_im(2:3),'(I2)') im
          ELSE IF (im .lt. 1000) THEN
            WRITE(File_im(1:3),'(I3)') im
          END IF


          OPEN(1000+im,file='single'//'_'//File_im//'.xyz',
     & position='append',form='formatted')

          WRITE(1000+im,*) natom
          WRITE(1000+im,*) " "

          DO ia = 1, natom
            WRITE(1000+im,10) type1(i,im,ia) 
     & ,x(i,im,ia),y(i,im,ia),z(i,im,ia) 
          END DO 

          CLOSE(1000+im)
 
        END DO

      END DO


10    format (I2,1x,3(f12.6,1x))



!     Calculate Rete vector !

      DO i = 1, ntime
         DO im = 1, nmol
            Rx(i,im) = ux(i,im,1)-ux(i,im,natom)
            Ry(i,im) = uy(i,im,1)-uy(i,im,natom)
            Rz(i,im) = uz(i,im,1)-uz(i,im,natom)
            R_mag(i,im) = sqrt(Rx(i,im)*Rx(i,im)+Ry(i,im)*Ry(i,im)+
     &                   Rz(i,im)*Rz(i,im))
 
          IF (i .gt. ntime*0.5) THEN
           theta(i,im) = acos(Rx(i,im)/R_mag(i,im))*(180/pi) 
           WRITE (5,*) i, im, theta(i,im)
c           WRITE (3,11)i,im,R_mag(i,im)**2,Rx(i,im),Ry(i,im),Rz(i,im)
         END IF
 
         uRx(i,im) = Rx(i,im) / R_mag(i,im)
         uRy(i,im) = Ry(i,im) / R_mag(i,im)
         uRz(i,im) = Rz(i,im) / R_mag(i,im)
 
         END DO
      END DO
 
      DO i= 1, ntime 
        IF (i .gt. (ntime*0.5)) THEN
          DO im = 1, nmol

            ibin1 = NINT (R_mag(i,im) / delr1)
             histo1(ibin1) = histo1(ibin1) + 1.d0
  
            ibin2 = NINT (theta(i,im) / delr2)
            histo2(ibin2) = histo2(ibin2) + 1.d0

          END DO
        END IF
      END DO
  
      ! PDF for Rete !
      DO ib1 = 0, ndist
        WRITE(4,*)ib1*delr1,histo1(ib1),histo1(ib1)/(nmol*
     &(ntime*0.5))
      END DO

      ! PDF for Orientation angle !
      DO ib2 = 0, 180
        WRITE(6,*)ib2*delr2,histo2(ib2),histo2(ib2)/(nmol*
     &(ntime*0.5))
      END DO 
  
!     Calculation of time correlation function  !

      DO i = 1, nmol
        tau = -1
        DO ia = 1, ntime
          tau = tau + 1
          nt = ntime - tau
          DO k = 1, nt
            sumx(i,ia) = sumx(i,ia) + uRx(k,i)*uRx(k+tau,i)
            sumy(i,ia) = sumy(i,ia) + uRy(k,i)*uRy(k+tau,i)
            sumz(i,ia) = sumz(i,ia) + uRz(k,i)*uRz(k+tau,i)
            sumxy(i,ia) = sumxy(i,ia) + uRx(k,i)*uRy(k+tau,i)
            sumxyz(i,ia) = sumxyz(i,ia) +(uRx(k,i)*uRx(k+tau,i) +
     & uRy(k,i)*uRy(k+tau,i) + uRz(k,i)*uRz(k+tau,i))

          END DO
          xx(i,ia) = sumx(i,ia)/dble(nt)
          xy(i,ia) = sumxy(i,ia)/dble(nt)
          yy(i,ia) = sumy(i,ia)/dble(nt)
          zz(i,ia) = sumz(i,ia)/dble(nt)
          xyz(i,ia) = sumxyz(i,ia)/dble(nt)
        END DO
      END DO

      DO ia = 1, ntime
        DO i = 1, nmol
          tcfxx(ia) = tcfxx(ia) + xx(i,ia)
          tcfxy(ia) = tcfxy(ia) + xy(i,ia)
          tcfyy(ia) = tcfyy(ia) + yy(i,ia)
          tcfzz(ia) = tcfzz(ia) + zz(i,ia)
          tcf(ia) = tcf(ia) + xyz(i,ia)
        END DO
        tcfxx(ia) = tcfxx(ia)/dble(nmol)
        tcfxy(ia) = tcfxy(ia)/dble(nmol)
        tcfyy(ia) = tcfyy(ia)/dble(nmol)
        tcfzz(ia) = tcfzz(ia)/dble(nmol)
        tcf(ia) = tcf(ia)/dble(nmol)
      END DO


      OPEN (7,FILE='tcf_uRete.dat')
      DO i = 1, ntime*0.8
         WRITE(7,*) i, i*dt, tcf(i)
      END DO

      OPEN (8,FILE='tcf_xyz_uRete.dat')
      DO i = 1, ntime
         WRITE(8,12) i,tcf(i),tcfxx(i),tcfxy(i),tcfyy(i),tcfzz(i)
      END DO

12    format(I6,1x,4(f12.8,1x))
20    format(4(f20.10,2x))
30    format(I5,1x,f10.2,1x,f20.10)

11      format(I6,1x,I3,1x,4(f12.6,1x))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! CALCULATE THE RADIUS OF GYRATION !!!!!


c      OPEN(9, file='Rg_mol.dat', position='append', form='formatted')
      OPEN(10, file='Rg_PDF_0.25.dat', position='append',
     & form='formatted')


      DO i = 1, ntime

      IF (i .gt. 0.5*ntime) THEN

      DO im = 1, nmol 
      
      x_cm = 0.0; y_cm = 0.0; z_cm = 0.0; totmass = 0.0

        DO ia = 1, natom
c          totmass(i,im) = totmass(i,im) + mass(1)
c          x_cm(i,im) = x_cm(i,im) + mass(1) * ux(i,im,ia)
c          y_cm(i,im) = y_cm(i,im) + mass(1) * uy(i,im,ia)
c          z_cm(i,im) = z_cm(i,im) + mass(1) * uz(i,im,ia)

          x_cm(i,im) = x_cm(i,im) + mass(type1(i,im,ia)) * ux(i,im,ia)
          y_cm(i,im) = y_cm(i,im) + mass(type1(i,im,ia)) * uy(i,im,ia)
          z_cm(i,im) = z_cm(i,im) + mass(type1(i,im,ia)) * uz(i,im,ia)
        END DO
        x_cm(i,im) = x_cm(i,im) / totmass(i,im)
        y_cm(i,im) = y_cm(i,im) / totmass(i,im)
        z_cm(i,im) = z_cm(i,im) / totmass(i,im)

        Rgsq = 0.0

        DO ia = 1, natom

          Rgx(i,im) = Rgx(i,im) + (ux(i,im,ia) - x_cm(i,im))
          Rgy(i,im) = Rgy(i,im) + (uy(i,im,ia) - y_cm(i,im))
          Rgz(i,im) = Rgz(i,im) + (uz(i,im,ia) - z_cm(i,im))

          Rgsq(i,im) = Rgsq(i,im) + (Rgx(i,im)**2 + Rgy(i,im)**2 + 
     & Rgz(i,im)**2)
        
        END DO

        Rg(i,im) = SQRT(Rgsq(i,im) / natom)

c        PRINT *, Rg(i,im)
  
        WRITE (11,*) i, im, Rg(i,im)
         
        ibin3 = NINT (Rg(i,im) / delr1)
        histo3(ibin3) = histo3(ibin3) + 1.d0
                
 
      END DO   ! for nmol
 
      END IF
 
      END DO   ! for ntime

      DO ib3 = 0, ndist
        WRITE(10,*)ib3*delr1,histo3(ib3),histo3(ib3)/(nmol*
     &(ntime*0.5))
      END DO
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      OPEN(11, file='avg_R.dat', position='append', form='formatted')
      OPEN(12, file='GT.dat', position='append', form='formatted')
      OPEN(13, file='CT.dat', position='append', form='formatted')
        
      DO i = 1, ntime
        
        IF (i .gt. 0.5*ntime) THEN
        
        DO im = 1, nmol

        sum_Rg(i) = sum_Rg(i) + Rg(i,im)
        sum_Rgxx(i) = sum_Rgxx(i) + Rgx(i,im)**2
        sum_Rgxy(i) = sum_Rgxy(i) + Rgx(i,im)*Rgy(i,im)
        sum_Rgyy(i) = sum_Rgyy(i) + Rgy(i,im)**2
        sum_Rgzz(i) = sum_Rgzz(i) + Rgz(i,im)**2

        sum_Rete(i) = sum_Rete(i) + R_mag(i,im)  
        sum_Retexx(i) = sum_Retexx(i) + Rx(i,im)**2
        sum_Retexy(i) = sum_Retexy(i) + Rx(i,im)*Ry(i,im)
        sum_Reteyy(i) = sum_Reteyy(i) + Ry(i,im)**2
        sum_Retezz(i) = sum_Retezz(i) + Rz(i,im)**2

        END DO
        
        avg_Rg(i) = sum_Rg(i)/nmol
        avg_Rgxx(i) = sum_Rgxx(i)/nmol
        avg_Rgxy(i) = sum_Rgxy(i)/nmol
        avg_Rgyy(i) = sum_Rgyy(i)/nmol
        avg_Rgzz(i) = sum_Rgzz(i)/nmol

        avg_Rete(i) = sum_Rete(i)/nmol
        avg_Retexx(i) = sum_Retexx(i)/nmol
        avg_Retexy(i) = sum_Retexy(i)/nmol
        avg_Reteyy(i) = sum_Reteyy(i)/nmol
        avg_Retezz(i) = sum_Retezz(i)/nmol
        
        WRITE (11,*) i, avg_Rg(i), avg_Rete(i)

        WRITE (12,*) i, avg_Rgxx(i),avg_Rgxy(i),avg_Rgyy(i),avg_Rgzz(i)
        WRITE (13,*) i, avg_Retexx(i), avg_Retexy(i), avg_Reteyy(i), 
     & avg_Retezz(i)

        END IF
        
      END DO

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! REARRANGE THE PRESSURE TENSOR !!!!!

        OPEN(unit=14,file='press.tensor',status='unknown') !,form='formatted')
        OPEN(unit=15, file='PT_tsr.dat', position='append')

        ! Skip initial header lines
        READ(14,*) txt
        READ(14,*) txt
        READ(14,*) txt

        DO i = 1, 3000  !ntime
          READ(14,*) txt
          DO j = 1, 6
            READ(14,*) PT(i,j,1), PT(i,j,2)
          END DO
          WRITE(15,'(6(F12.5,1X))') (PT(i,j,2), j=1,6)
        END DO

        CLOSE(14)
        CLOSE(15)

  
      STOP

      END PROGRAM POST

