c  ifort -fast -o bisp_fast_bin_fftw2_quad.exe bisp_fast_bin_fftw2_quad.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw
      implicit none
      include '/usr/local/src/fftw-2.1.5/fortran/fftw_f77.i'
      integer nside,nsideD,Dim,nn(3),nmax,ncheck,istep
      integer*8 planb,planf
c* INPUT (Check) ************************************************************
c      parameter(nside=400,nmax=50)
c      parameter(nside=240,nmax=40)
      parameter(nside=375,nmax=40)
c****************************************************************************
      parameter(Dim=3,nsideD=nside**Dim) 
      data nn/nside,nside,nside/
      integer i, j, l, m, n, nk(0:2*nside),nbk(0:2*nside+1),iflag,irsd
      integer iseed, indx(nsideD), id, jd, ld, nmodes, ndim,k,ix,iy,iz
c      real mapk(nsideD,nmax),m1(2,nsideD), map1(nsideD), map2(nsideD)  
      real pow(nmax),I10,I12,I22,I13,I23,I33,P0,alpha
      real pow2(nmax), I12d,I13d,alphaSTD, poww(nmax)
      real, allocatable :: mapk(:,:),m1(:,:),map1(:),map2(:)
      real, allocatable :: m1xx(:,:),map1xx(:),map2xx(:)
      real, allocatable :: m1w(:,:),map1w(:),map2w(:)
c      real normk(nsideD,nmax), norm1(2,nsideD)
      real di, dj, dl, eightpi2, bi, step, bi2, biw
      real dist, bisp(nmax,nmax,nmax), q(nmax,nmax,nmax)
      real bisp2(nmax,nmax,nmax), q2(nmax,nmax,nmax)
      real bispw(nmax,nmax,nmax) 
      real*8 coun(nmax,nmax,nmax), avg, sum, sum2, avgw
      real rand,fac,Fr,Fi,Pt,Rr,Ri
      complex, allocatable :: dclr1(:,:,:),dclr2(:,:,:)
      complex, allocatable :: dclr3(:,:,:),dclr4(:,:,:)
      complex, allocatable :: dclr5(:,:,:),dclr6(:,:,:)
      integer nobj,ip,aexp,p,dpnp1,dpn,ibuf,lm,mcode,Ncut,ng,Ncuts
      character*255 filebisp, filecoef,filecounts,rfile,gfile,homedir
      character*200 iflagstr,irsdstr
      common/discrete/I10,I12,I22,I13,I23,I33,P0,alpha
      common/discrete2/I12d,I13d,alphaSTD

c* INPUT (Check) ************************************************************
      filecounts='~/Code/Fortran/counts2quad_n360_nmax40_ncut3_s3'
c      filecounts='~/Code/Fortran/counts2quad_n240_nmax40_ncut2_s2'
c      filecounts='~/Code/Fortran/counts2quad_n360_nmax60_ncut2_s2'
      Ncut=3
      istep=3
      ncuts=ncut/int(istep)
c****************************************************************************
      
      call fftw3d_f77_create_plan(planf,nside,nside,nside,
     $ FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)

      step=real(istep)
      eightpi2 = 78.95683521
     
      allocate(map1(nside**Dim),map2(nside**Dim))
      allocate(map1xx(nside**Dim),map2xx(nside**Dim))

c      write(*,*) 'Nbody/periodic (1) or data/cut-sky (2)?'
c      read(*,*)iflag
      call getarg(1,iflagstr)
      read(iflagstr,*)iflag 

      if (iflag.eq.1) then

c         write(*,*) 'Fourier file :'
c         read(*,'(a)') filecoef
c         write(*,*) 'Bispectrum file :'
c         read(*,'(a)') filebisp
         call getarg(2,filecoef)
         call getarg(3,filebisp)
         call getarg(4,irsdstr) !redshift-space direction (1:x,2:y,3:z)
         read(irsdstr,*)irsd 

         allocate(dclr1(nside/2+1,nside,nside))
         allocate(dclr2(nside/2+1,nside,nside))
         write(*,*)'memory allocated'
c         call inputNB(nside,Dim,filecoef,map1,map2,dclr1)
         call inputNB(nside,Dim,filecoef,map1,map2,dclr1,dclr2,
     $                map1xx,map2xx,irsd)
         deallocate(dclr1,dclr2)
         write(*,*)'input done and memory deallocated'

      elseif (iflag.eq.2) then

C         write(*,*) 'Random Fourier file :'
C         read(*,'(a)') rfile
C         write(*,*) 'Data Fourier file :'
C         read(*,'(a)') gfile
C         write(*,*) 'Bispectrum file :'
C         read(*,'(a)') filebisp         
         call getarg(2,rfile)
         call getarg(3,gfile)
         call getarg(4,filebisp)

         allocate(dclr1(nside/2+1,nside,nside))
         allocate(dclr2(nside/2+1,nside,nside))
         allocate(dclr3(nside/2+1,nside,nside))
         allocate(dclr4(nside/2+1,nside,nside))
         allocate(dclr5(nside/2+1,nside,nside))
         allocate(dclr6(nside/2+1,nside,nside))
         allocate(map1w(nside**Dim),map2w(nside**Dim))
         write(*,*)'memory allocated'
         call inputDA(nside,Dim,gfile,rfile,map1,map2,dclr1,dclr2,
     $        map1xx,map2xx,dclr3,dclr4,map1w,map2w,dclr5,dclr6)
         deallocate(dclr1,dclr2,dclr3,dclr4,dclr5,dclr6)
         write(*,*)'input done and memory deallocated'

      endif
      do l=1,nmax
         do j=1,nmax
            do i=1,nmax
               bisp(i,j,l) = 0.
               bisp2(i,j,l) = 0.
               if (iflag.eq.2) then
                  bispw(i,j,l) = 0.
               endif
            enddo
         enddo
      enddo
      
      do i=0,2*nside
         nk(i)=0
      enddo

      write(*,*) 'Find modes of amplitude |k|'

      do i = 0, nside -1
         do j = 0, nside  -1
            do l = 0, nside  -1
               di = float(min(i,nside-i))
               dj = float(min(j,nside-j))
               dl = float(min(l,nside-l))
               dist = sqrt(di*di + dj*dj +dl*dl)
               nk(int(dist/step+0.5)) = nk(int(dist/step+0.5)) + 1 
            enddo 
         enddo 
      enddo 
      
      nbk(0) = 0 
      do i = 0, 2*nside
         nbk(i+1) = nbk(i) + nk(i)
         nk(i) = 0 
      enddo 
      
      write(*,*) 'Save coordinates'
  
      m = 0
      do i = 0, nside -1 
         do j = 0, nside -1
            do l = 0, nside -1
               di = float(min(i,nside-i))
               dj = float(min(j,nside-j))
               dl = float(min(l,nside-l))
               dist = sqrt(di*di + dj*dj +dl*dl)
               nk(int(dist/step+0.5)) = nk(int(dist/step+0.5)) + 1  
               n = nbk(int(dist/step+0.5)) + nk(int(dist/step+0.5))
               m = m+ 1
               indx(n) = m 
            enddo 
         enddo 
      enddo 

      write(*,*) 'Calculate k maps'

      ndim=3

      allocate(m1(2,nsideD),mapk(nsideD,nmax),m1w(2,nsideD))
      do i = ncuts, nmax
         do n = 1, nsideD 
            m1(1,n) = 0.
            m1(2,n) = 0.
            if (iflag.eq.2) then
               m1w(1,n) = 0.
               m1w(2,n) = 0.
            endif   
         enddo 
         nmodes = 0
         do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = map1(indx(n))
            m1(2,indx(n)) = map2(indx(n))
            if (iflag.eq.2) then
               m1w(1,indx(n)) = map1w(indx(n))
               m1w(2,indx(n)) = map2w(indx(n))
            endif   
            nmodes =nmodes +1 
         enddo 
         call fftwnd_f77_one(planf,m1,0)
         avg = 0.d0 
         if (iflag.eq.2) then
            call fftwnd_f77_one(planf,m1w,0)
            avgw = 0.d0 
         endif   
         do n = 1, nsideD
            avg = avg + dble(m1(1,n))*dble(m1(1,n))
            if (iflag.eq.2) then
               avgw = avgw + dble(m1w(1,n))*dble(m1(1,n))
            endif
            mapk(n,i) = m1(1,n)
         enddo
         pow(i)=real(avg)/float(nsided)/float(nmodes)
         if (iflag.eq.2) then 
            pow(i)=(pow(i)-(1.+alpha)*I12)/I22 !correct discreteness
            poww(i)=real(avgw)/float(nsided)/float(nmodes)
         endif
      enddo
      
      deallocate(m1,map1)
      if (iflag.eq.2) then
         deallocate(m1w,map1w)
      endif   

      write(*,*) 'Read counts'

      open(unit=2,status='old',form='unformatted',file=filecounts)
      read(2) coun
      close(2)

      write(*,*) 'counts were read'
      allocate(m1xx(2,nsideD))

      write(*,*) 'Compute bisp monopole and quadrupole!'
      
      
      do l = ncuts, nmax !keep ffts in outer loop
         do n = 1, nsideD 
            m1xx(1,n) = 0.        
            m1xx(2,n) = 0.        
         enddo
         nmodes = 0
         do n = nbk(l)+1, nbk(l+1)
            m1xx(1,indx(n)) = map1xx(indx(n))
            m1xx(2,indx(n)) = map2xx(indx(n))
            nmodes = nmodes + 1 
         enddo 
         call fftwnd_f77_one(planf,m1xx,0)
         avg = 0.d0 
         do n = 1, nsideD
            avg = avg + dble(m1xx(1,n))*dble(mapk(n,l))
         enddo
         pow2(l)=real(avg)/float(nsided)/float(nmodes) !quadrupole power
         if (iflag.eq.2) then 
            pow2(l)=(pow2(l))/I22 !correct discreteness
         endif

         do j = ncuts, l
            do i = max(ncuts,l-j), j
               sum = 0.d0 
               sum2 = 0.d0 
               do n = 1, nsideD
                  sum = sum 
     $            + dble(mapk(n,i))*dble(mapk(n,j))*dble(mapk(n,l))
                  sum2 = sum2 
     $            + dble(mapk(n,i))*dble(mapk(n,j))*dble(m1xx(1,n))
               enddo
               bi=real(sum/coun(i,j,l))
               bi2=real(sum2/coun(i,j,l))
               if (iflag.eq.2) then !shot noise correction
                  biw=(bi-(poww(i)+poww(j)+poww(l)) +
     $               2.*(I13d-alpha**2*I13))/I33
                  bi=(bi-(pow(i)+pow(j)+pow(l))*I23 -
     $               (1.-alpha**2)*I13)/I33
                  bi2=(bi2-pow2(l)*I23)/I33
               endif
               bisp(i,j,l)=bi
               if (iflag.eq.2) then
                  bispw(i,j,l)=biw
               endif                  
               q(i,j,l)=bi/(pow(i)*pow(j)+pow(j)*pow(l)
     $                 +pow(l)*pow(i))
               bisp2(i,j,l)=bi2
               q2(i,j,l)=bi2/(pow(i)*pow(j)+pow(j)*pow(l)
     $                 +pow(l)*pow(i))
             enddo
         enddo 
      enddo          
      
      open(unit=7,file=filebisp,status='unknown',form='formatted')
      write(*,*) 'output'
      do l = ncuts, nmax
         do j = ncuts, l
            do i = ncuts,j
               fac=1.
               if(coun(i,j,l).ne.0.d0) then 
                  if(j.eq.l .and. i.eq.j) fac=6.
                  if(i.eq.j .and. j.ne.l) fac=2.
                  if(i.eq.l .and. l.ne.j) fac=2.
                  if(j.eq.l .and. l.ne.i) fac=2.
                  coun(i,j,l)=coun(i,j,l)/dble(fac*float(nsided))
                  if (iflag.eq.1) then 
                  write(7,1000) int(step)*l,int(step)*j,int(step)*i,
     $                 pow(l),pow(j),pow(i),bisp(i,j,l),q(i,j,l),
     $                 pow2(l),pow2(j),pow2(i),bisp2(i,j,l),q2(i,j,l),
     $                 real(coun(i,j,l))
                  else
                  write(7,1000) int(step)*l,int(step)*j,int(step)*i,
     $                 pow(l),pow(j),pow(i),bisp(i,j,l),q(i,j,l),
     $                 pow2(l),pow2(j),pow2(i),bisp2(i,j,l),q2(i,j,l),
     $                 real(coun(i,j,l)),bispw(i,j,l)                  
                  endif
               endif
            enddo
         enddo
      enddo 
      close(7)

 1000 format(3I4, 12e13.5)
 2000 stop 
      end 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inputNB(nside,Dim,filecoef,map1,map2,dclr1,dclr2,
     $                map1xx,map2xx,irsd)
      implicit none
      integer nside,Dim,n,i,j,l,id,jd,ld,Lx,Ly,Lz,Npar
      integer irsd,ix,ikx,iy,iky,iz,ikz
      real akx,aky,akz,phys_nyq,rk,amu
      real map1(nside**Dim), map2(nside**Dim) 
      real map1xx(nside**Dim), map2xx(nside**Dim) 
      complex dclr1(nside/2+1,nside,nside)
      complex dclr2(nside/2+1,nside,nside)
      character*255 filecoef

      open(unit=1,status='old',file=filecoef,form='unformatted')
      read(1) Lx,Ly,Lz,Npar,akx,aky,akz,phys_nyq
      Lx=nside
      Ly=Lx
      Lz=Lx
      read(1)(((dclr1(i,j,l),i=1,Lx/2+1),j=1,Ly),l=1,Lz) !delta
      close(1)
      write(*,*)'data read'

      do 100 iz=1,Lz !build quadrupole
         ikz=mod(iz+Lz/2-2,Lz)-Lz/2+1
         do 100 iy=1,Ly
            iky=mod(iy+Ly/2-2,Ly)-Ly/2+1
            do 100 ix=1,Lx/2+1
               ikx=mod(ix+Lx/2-2,Lx)-Lx/2+1
               rk=sqrt(float(ikx**2+iky**2+ikz**2))
               if(rk.gt.0.)then
                  if (irsd.eq.3) then 
                     amu=float(ikz)/rk
                  elseif (irsd.eq.2) then
                     amu=float(iky)/rk
                  elseif (irsd.eq.1) then
                     amu=float(ikx)/rk !unit vectors
                  else
                     write(*,*)'choose RSD direction'
                     stop
                  endif   
                  dclr2(ix,iy,iz)=(7.5*amu**2-2.5)*dclr1(ix,iy,iz) !weight by 5*Leg2
               end if
 100  continue
      write(*,*)'quadrupole built'


      n = 0 
      do i = 1, nside 
         do j = 1, nside 
            do l = 1, nside 
               n = n+1
               if (i .le. nside/2+1) then 
                  map1(n) = real(dclr1(i,j,l))
                  map2(n) = aimag(dclr1(i,j,l))
                  map1xx(n) = real(dclr2(i,j,l))
                  map2xx(n) = aimag(dclr2(i,j,l))
               else 
                  id = mod(nside-i+1,nside)+1
                  jd = mod(nside-j+1,nside)+1
                  ld = mod(nside-l+1,nside)+1
                  map1(n) = real(dclr1(id,jd,ld))
                  map2(n) = -aimag(dclr1(id,jd,ld))
                  map1xx(n) = real(dclr2(id,jd,ld))
                  map2xx(n) = -aimag(dclr2(id,jd,ld))
               endif
               if (mod(i-1,nside/2)+mod(j-1,nside/2)+
     $              mod(l-1,nside/2).eq.0) then
                  map2(n) = 0.
                  map2xx(n) = 0.
               endif  
            enddo
         enddo 
      enddo
c      map2(1)=0.
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inputDA(nside,Dim,gfile,rfile,map1,map2,dcg,dcr,
     $           map1xx,map2xx,dcgxx,dcrxx,map1w,map2w,dcgw,dcrw)
      implicit none
      integer nside,Dim,n,i,j,l,id,jd,ld,Lx,Ly,Lz,Nr,Ng
      real akx,aky,akz,phys_nyq,I10,I12,I22,I13,I23,I33,P0,alpha,Ngsys
      real map1(nside**Dim), map2(nside**Dim) ,Nrsyscomp,Nrsystot
      real map1xx(nside**Dim), map2xx(nside**Dim) 
      real map1w(nside**Dim), map2w(nside**Dim) 
      complex dcg(nside/2+1,nside,nside),dcgxx(nside/2+1,nside,nside)
      complex dcr(nside/2+1,nside,nside),dcrxx(nside/2+1,nside,nside)
      complex dcrw(nside/2+1,nside,nside),dcgw(nside/2+1,nside,nside)
      character*255 gfile,rfile
      integer Ngrid,I12F,Ngrid2
      real Ngsyscomp,xminD,xmaxD,yminD,ymaxD,zminD,zmaxD,Ngsystot
      real Nrsys,thetaobs,phiobs,xminR,xmaxR,yminR,ymaxR,zminR,zmaxR
      real I12d,I12xx,I12yy,I12zz,I12xy,I12yz,I12zx,P0m,I13d,alphaSTD
      common/discrete/I10,I12,I22,I13,I23,I33,P0,alpha
      common/discrete2/I12d,I13d,alphaSTD
      
      open(unit=4,file=gfile,status='old',form='unformatted')
      read(4)Ngrid
      if (Ngrid.ne.nside) then
         write(*,*)'Ngrid does not match nside',Ngrid,nside
         stop
      endif      
      read(4)dcg
cc      read(4)I10,I12d,I22,I13,I23,I33,I12F !as in v2 of FFT
      read(4)I10,I12d,I22,I13d,I23,I33 !v3 of FFT
      read(4)P0m,Ng,Ngsys,Ngsyscomp !,Ngsystot
      read(4)xminD,xmaxD,yminD,ymaxD,zminD,zmaxD
      read(4)dcgxx
      read(4)dcgw !for now ignore
      close(4)

      open(unit=4,file=rfile,status='old',form='unformatted')
      read(4)Ngrid2
      if (Ngrid2.ne.Ngrid) then
         write(*,*)'Ngrids do not match',Ngrid,Ngrid2
         stop
      endif      
      read(4)dcr
      read(4)I10,I12,I22,I13,I23,I33
      read(4)P0,Nr,Nrsys,Nrsyscomp,Nrsystot
      read(4)thetaobs,phiobs
      read(4)xminR,xmaxR,yminR,ymaxR,zminR,zmaxR
      read(4)dcrxx
      read(4)I12xx,I12yy,I12zz,I12xy,I12yz,I12zx
      read(4)dcrw !for now ignore
      close(4)

      if (P0m.ne.P0) then
         write(*,*)'P0s do not match',P0m,P0
         stop
      endif

      if (P0m.ne.P0) then
         write(*,*)'P0s do not match',P0m,P0
         stop
      endif
      if (xminD.lt.xminR) then 
         write(*,*)'WARNING: xminD=',xminD,'<xminR=',xminR
      endif
      if (yminD.lt.yminR) then 
         write(*,*)'WARNING: yminD=',yminD,'<yminR=',yminR
      endif
      if (zminD.lt.zminR) then 
         write(*,*)'WARNING: zminD=',zminD,'<zminR=',zminR
      endif
      if (xmaxD.gt.xmaxR) then 
         write(*,*)'WARNING: xmaxD=',xmaxD,'>xmaxR=',xmaxR
      endif
      if (ymaxD.gt.ymaxR) then 
         write(*,*)'WARNING: ymaxD=',ymaxD,'>ymaxR=',ymaxR
      endif
      if (zmaxD.gt.zmaxR) then 
         write(*,*)'WARNING: zmaxD=',zmaxD,'>zmaxR=',zmaxR
      endif
      alphaSTD=float(Ng)/float(Nr) !now scale random integrals by alpha
      write(*,*)'alphaSTD is',alphaSTD
      alpha=Ngsys/Nrsys !now scale random integrals by alphasys
      write(*,*)'alpha sys is',alpha
      alpha=Ngsystot/Nrsystot !now scale random integrals by everything
      write(*,*)'alpha tot is',alpha
      alpha=Ngsyscomp/Nrsyscomp !now scale random integrals by everything that is not FKP
      write(*,*)'alpha syscomp is',alpha
      I10=I10*alpha
      I12=I12*alpha
      I22=I22*alpha !STD
      I13=I13*alpha
      I23=I23*alpha
      I33=I33*alpha !STD
      I12xx=I12xx*alpha
      I12yy=I12yy*alpha
      I12zz=I12zz*alpha
      I12xy=I12xy*alpha
      I12yz=I12yz*alpha
      I12zx=I12zx*alpha
      

      do l=1,nside
         do j=1,nside
            do i=1,nside/2+1
               dcg(i,j,l)=dcg(i,j,l)-alpha*dcr(i,j,l)
               dcgxx(i,j,l)=dcgxx(i,j,l)-alpha*dcrxx(i,j,l)
               dcgw(i,j,l)=dcgw(i,j,l)+alpha**2*dcrw(i,j,l)
            enddo
         enddo
      enddo
      
 
      n = 0 
      do i = 1, nside 
         do j = 1, nside 
            do l = 1, nside 
               n = n+1
               if (i .le. nside/2+1) then 
                  map1(n) = real(dcg(i,j,l))
                  map2(n) = aimag(dcg(i,j,l))
                  map1xx(n) = real(dcgxx(i,j,l))
                  map2xx(n) = aimag(dcgxx(i,j,l))
                  map1w(n) = real(dcgw(i,j,l))
                  map2w(n) = aimag(dcgw(i,j,l))
               else 
                  id = mod(nside-i+1,nside)+1
                  jd = mod(nside-j+1,nside)+1
                  ld = mod(nside-l+1,nside)+1
                  map1(n) = real(dcg(id,jd,ld))
                  map2(n) = -aimag(dcg(id,jd,ld))
                  map1xx(n) = real(dcgxx(id,jd,ld))
                  map2xx(n) = -aimag(dcgxx(id,jd,ld))
                  map1w(n) = real(dcgw(id,jd,ld))
                  map2w(n) = -aimag(dcgw(id,jd,ld))
               endif
               if (mod(i-1,nside/2)+mod(j-1,nside/2)+
     $              mod(l-1,nside/2).eq.0) then
                  map2(n) = 0.
                  map2xx(n) = 0.
                  map2w(n) = 0.
               endif  
            enddo
         enddo 
      enddo
c      map2(1)=0.
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
