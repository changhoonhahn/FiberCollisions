c ifort -arch pn4 -O4 -convert big_endian -o FFTil4_aniso_be.exe FFTil4_aniso.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw   
c ifort -arch pn4 -O4 -o FFTil4_aniso_Lag.exe FFTil4_aniso.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw   
      parameter(Nparm=2*10**8)
      parameter(pi=3.141592654,tpi=2.*pi)
      real,allocatable :: r(:,:)
      real rm(6),xscale,rlow,Rx,Ry,Rz,akx,aky,akz,phys_nyq
      integer Npar,Lx,Ly,Lz,Ngroup,i 
      integer*8 planf
      character Lxstr*200,Lystr*200,Lzstr*200
      character infile*200
      character mockfile*200
      character paramfile*200
      include '/home/users/hahn/powercode/fftw_f77.i'
     
      write(*,*) 'Particle Positions file :'
      call getarg(1,infile)
      Ngroup=1
      write(*,*)'allocate memory'
      allocate(r(3,Nparm))

      write(*,*)'reading file'
      call input(Nparm,r,infile,xscale,rlow,Npar,Rx,Ry,Rz,Ngroup)
      write(*,*)'xscale=',xscale
      write(*,*)'Npar=',Npar
      write(*,*)'Rx,Ry,Rz=',Rx,Ry,Rz
      write(*,*)'rlow=',rlow
      call getarg(2,mockfile)

      open(unit=8,file=mockfile,status='unknown',form='formatted')
      do i=1,Npar
        write(8,'(3(F,2x))') r(1,i), r(2,i), r(3,i)
      enddo
      close(8)
      
      call getarg(3,paramfile)

      open(unit=9,file=paramfile,status='unknown',form='formatted')
      write(9,'(F,2x,I,2x,4(F,2x))') xscale,Npar,Rx,Ry,Rz,rlow
      close(9)

      stop
      end
cc*******************************************************************
      subroutine input(N,r,startfile,xscale,rlow,Np,Rx,Ry,Rz,Ngroup)
cc*******************************************************************
      real r(3,N),xscale,rlow,Rx,Ry,Rz,r1,r2,r3
      integer Np,idum,ipar,Ngroup
      character startfile*(*)

      rlow=0.
      
      open(11,file=startfile,status='old',form='unformatted')
cc     This is for Fastfood format:      
      read(11) idum, Np,idum,idum,idum
      read(11) Rx,idum,idum,idum,idum,idum,idum,idum,idum
      read(11) idum

!c      read(11) Np, Rx !,Ry,Rz
cc     This is for std format:      
c      read(11) Np
      Rx=2400. !Oriana
c      Rx=1000. !Carmen
cc      Rx=7680.  !MICE
      xscale=Rx
      Ry=Rx
      Rz=Rx
      if (Np.gt.N) then
         write(*,*)'Number of particles:',Np
         write(*,*)'too large, increase Nparm'
         stop
      endif

cc     This is for Fastfood format:      
      read(11)(r(1,j),j=1,Np)
      read(11)(r(2,j),j=1,Np)
      read(11)(r(3,j),j=1,Np)

cc     This is for std format:      
c      read(11)((r(i,j),i=1,3),j=1,Np)
      close(11)
c      write(*,*)'file read'


c!      Np=0
c!      open(11,file=startfile,status='old',form='formatted')
c!      do i=1,N
cc         read(11,*,end=12)r(1,i),r(2,i),r(3,i) !,idum !
c         read(11,*,end=12)ipar,i1,i2,i3,idum,idum,idum !mice halos
c         if(ipar.gt. Ngroup) then
c            Np=Np+1
c            r(1,Np)=float(i1)
c            r(2,Np)=float(i2)
c            r(3,Np)=float(i3)
c         endif
c!         read(11,*,end=12)idum,idum,r1,r2,r3 !,idum,idum,idum !Lagrange halos
         !oriana1 for some reason has no velocities e.g.
         ! /export/chichipio1/kcc274/Oriana_Gaussian/Oriana_1/LagrangeHalo_FoFpt156_20_36_z0.dat
c!            Np=Np+1
c!            r(1,Np)=r1
c!            r(2,Np)=r2
c!            r(3,Np)=r3
c!      enddo
c! 12   continue  
c!      close(11)
         write(*,*)'Number of particles:',Np

      rup=rlow+xscale
      rupr=rup*(1.-1.e-6)
      do i=1,Np !make sure there is no problem later with interpolation
         r(1,i)=max(rlow,min(rupr,r(1,i)))
         r(2,i)=max(rlow,min(rupr,r(2,i)))
         r(3,i)=max(rlow,min(rupr,r(3,i)))
      end do
      return
      end
cc*******************************************************************
