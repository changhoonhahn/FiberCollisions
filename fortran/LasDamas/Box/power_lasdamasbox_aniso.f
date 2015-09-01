c     power spectrum multipoles
c  ifort -arch pn4 -O4 -o power3s_aniso.exe power3s_aniso.f -L/usr/local/fftw_intel/lib -lsrfftw -lsfftw
      parameter(Nbinm=500)
      complex,allocatable :: dclr(:,:,:)
      complex ct
      real mode,rev,imv,amp,pha,phi1,cot1,sit1
      real avgk(Nbinm),avgP(Nbinm),coga,Le2,Le4
      real RP(Nbinm),cp,sp,cc,thetaobs, phiobs,co(Nbinm)
      real avgP4(Nbinm),avgP2(Nbinm),Delta(Nbinm)
      integer a,b,c,Ngrid,Nbins,iobs,obs,ng2
      real kfun,akx,aky,akz,phys_nyq,rkx,rky,rkz
      integer Lx,Ly,Lz,Npar
      character*255 filepower,filecoef
      character*200 nbinstr
      pi=3.1415926

      write(*,*) 'Fourier file :'
      call getarg(1,filecoef)
      write(*,*) filecoef
      write(*,*) 'Power Spectrum file :'
      call getarg(2,filepower) 
      write(*,*) filepower
      write(*,*) 'Nbins :'
      call getarg(3,nbinstr) 
      read(nbinstr,*) Nbins
      if (Nbins.gt.Nbinm) then
         write(*,*)Nbins,'larger than',Nbinm
         write(*,*)'too large Nbins'
         stop
      endif
      iobs=1
      obs=mod(iobs,3)
      if (obs.eq.1) then
         thetaobs= 0.5*pi !0.837*pi
         phiobs=0. !0.595*pi
      elseif (obs.eq.2) then
         thetaobs=0.5*pi !0.527*pi
         phiobs=0.5*pi !1.53*pi
      elseif (obs.eq.0) then
         thetaobs=0. !0.235*pi
         phiobs=0. !0.783*pi
      endif

c      write(*,*) thetaobs/pi,phiobs/pi
c      write(*,*) 'INPUT theta, phi (in units of Pi) for the observer:'
c      read(*,*) thetaobs, phiobs

      do 10 i=1,Nbins
         avgk(i)=0.
         avgP(i)=0.
         avgP2(i)=0.
         avgP4(i)=0.
         RP(i)=0.
         co(i)=0.
 10   continue
      
c******************************************************************
      open(unit=2,status='old',file=filecoef,form='unformatted')
      read(2) Lx,Ly,Lz,Npar,akx,aky,akz,phys_nyq
      write(*,*)'allocate memory'
      allocate(dclr(Lx/2+1,Ly,Lz))
      read(2) (((dclr(ix,iy,iz),ix=1,Lx/2+1),iy=1,Ly),iz=1,Lz)
      close(2)
c******************************************************************
      do 100 iz=1,Lz
         icz=mod(Lz+1-iz,Lz)+1
         rkz=akz*float(mod(iz+Lz/2-2,Lz)-Lz/2+1) 
         do 100 iy=1,Ly
            icy=mod(Ly+1-iy,Ly)+1
            rky=aky*float(mod(iy+Ly/2-2,Ly)-Ly/2+1)
c            do 100 ix=1,Lx/2+1
c               rkx=akx*float(ix-1)
            do 100 ix=1,Lx
               icx=mod(Lx+1-ix,Lx)+1
               rkx=akx*float(mod(ix+Lx/2-2,Lx)-Lx/2+1)
               rk=sqrt(rkx**2+rky**2+rkz**2)
               imk=nint(Nbins*rk/phys_nyq)  
               if(imk.le.Nbins .and. imk.ne.0)then
                  cot1=rkz/rk
                  sit1=sqrt(1.-cot1*cot1)
                  if (sit1.gt.0.) then
                     cp=rkx/(rk*sit1)
                     sp=rky/(rk*sit1)
                     cc=sin(phiobs)*sp+cos(phiobs)*cp
                  else
                     cc=0.
                  endif
                  coga=cos(thetaobs)*cot1+sin(thetaobs)*sit1*cc
c        coga=cot1 !use z-direction as redshift mapping
                  Le2=-5.e-1+1.5e0*coga**2
                  Le4=3.75e-1-3.75e0*coga**2+4.375e0*coga**4
                  co(imk)=co(imk)+1.
                  if (ix.le.Lx/2+1) then 
                     ct=dclr(ix,iy,iz)
                  else !use cc
                     ct=dclr(icx,icy,icz)
                  endif
                  pk=(cabs(ct))**2
                  avgk(imk)=avgk(imk)+rk
                  avgP(imk)=avgP(imk)+pk
                  avgP2(imk)=avgP2(imk)+pk*5.*Le2
                  avgP4(imk)=avgP4(imk)+pk*9.*Le4
               end if
 100        continue
c******************************************************************
      open(4,file=filepower,status='unknown',form='formatted')
c      write(4,*)Npar
      do 110 Ibin=1,Nbins
         if(co(Ibin).gt.0.)then
            avgk(Ibin)=avgk(Ibin)/co(Ibin)  
            avgP(Ibin)=(avgP(Ibin)/co(Ibin)-1./(9.44451e-5))
            avgP2(Ibin)=avgP2(Ibin)/co(Ibin)/(akx*aky*akz)
            avgP4(Ibin)=avgP4(Ibin)/co(Ibin)/(akx*aky*akz)
c            avgP(Ibin)=avgP(Ibin)/co(Ibin)/(akx*aky*akz)
c            avgP2(Ibin)=avgP2(Ibin)/co(Ibin)/(akx*aky*akz)
c            avgP4(Ibin)=avgP4(Ibin)/co(Ibin)/(akx*aky*akz)
c            RP(Ibin)=avgP2(Ibin)/avgP(Ibin)
c            Delta(Ibin)=4*3.1415926*avgk(Ibin)**3*avgP(Ibin)
         write(4,1015) avgk(Ibin),avgP(Ibin),
     &   avgP2(Ibin),avgP4(Ibin),co(Ibin)
      end if
 110  continue
      close(4)

 1015 format(2x,5e16.6)
      stop
      end
c******************************************************************

