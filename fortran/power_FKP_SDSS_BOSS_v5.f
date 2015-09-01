! ifort -fast -o power_FKP_SDSS_BOSS_v5.exe power_FKP_SDSS_BOSS_v5.f
! use: power_FKP_SDSS_BOSS_v5.exe Fourier_R Fourier_D powerfile Lbox Nbins 
      implicit none ! as v4 but removes phiobs,thetaobs and Iijkl
      integer Ngrid,ix,iy,iz,Nbins,nyq,ikx,iky,ikz,imk,i,Ibin,Nr,Ng
      integer icz,icy,icx,Ngrid2
      complex, allocatable :: dcg(:,:,:),dcr(:,:,:)
      complex, allocatable :: dcgw(:,:,:),dcrw(:,:,:)
      complex, allocatable :: dcgxx(:,:,:),dcrxx(:,:,:)
      complex, allocatable :: dcgwxx(:,:,:),dcrwxx(:,:,:)
      complex, allocatable :: dcgxxxx(:,:,:),dcrxxxx(:,:,:)
      real*8, allocatable :: avgk(:),avgPg0(:),avgPr0(:),co(:)
      real*8, allocatable :: avgPg4fast(:),avgPg4fft15(:)
      real*8, allocatable :: avgPg2fast(:),avgPr2fast(:),avgPg0T(:)
      real*8, allocatable :: avgRWr0(:),avgIWr0(:)
      real*8, allocatable :: avgP00w(:),avgP02w(:),avgP20w(:)
      character filecoef*200,filecoefr*200,filepower*200,Rboxstr*200
      character Nbinstr*200
      real akfun,I10,I12,I22,I13,I23,I33,P0,alpha,P0m,rk
      real I12xx,I12yy,I12zz,I12xy,I12yz,I12zx, alphaSTD
      real I10d,I12d,I22d,I13d,I23d,I33d,cnorm
      real a12xx,a12yy,a12zz,a12xy,a12yz,a12zx,p4
      real pk,Ngsys,pq,thetaobs,phiobs,Nrsys,Ngsyscomp,Ngsystot
      real xminD,xmaxD,yminD,ymaxD,zminD,zmaxD
      real xminR,xmaxR,yminR,ymaxR,zminR,zmaxR
      real kxh,kyh,kzh,Nrsyscomp,Nrsystot,px,pp,py,pz
      complex ct,ctQ,ctx,ctP,cty
      
      call getarg(1,filecoefr) !Random Fourier file
      call getarg(2,filecoef) !Data Fourier file
      call getarg(3,filepower) !power spectrum file
      call getarg(4,Rboxstr) !box size (Mpc/h)
      read(Rboxstr,*)akfun
      call getarg(5,Nbinstr) !Number of bins
      read(Nbinstr,*)Nbins

      open(unit=4,file=filecoef,status='old',form='unformatted')
      read(4)Ngrid
      
      allocate(dcg(Ngrid/2+1,Ngrid,Ngrid),dcr(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcgw(Ngrid/2+1,Ngrid,Ngrid),dcrw(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcgwxx(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcrwxx(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcgxx(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcrxx(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcgxxxx(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcrxxxx(Ngrid/2+1,Ngrid,Ngrid))

      read(4)dcg
      read(4)I10d,I12d,I22d,I13d,I23d,I33d
      read(4)P0m,Ng,Ngsys,Ngsyscomp,Ngsystot
      read(4)xminD,xmaxD,yminD,ymaxD,zminD,zmaxD
      read(4)dcgxx ! 5 delta2
      read(4)dcgw
      read(4)dcgwxx
      read(4)dcgxxxx !9 delta4
      close(4)
      open(unit=4,file=filecoefr,status='old',form='unformatted')
      read(4)Ngrid2
      if (Ngrid2.ne.Ngrid) then
         write(*,*)'Ngrids do not match',Ngrid,Ngrid2
         stop
      endif      
      read(4)dcr
      read(4)I10,I12,I22,I13,I23,I33
      read(4)P0,Nr,Nrsys,Nrsyscomp,Nrsystot
      read(4)xminR,xmaxR,yminR,ymaxR,zminR,zmaxR
      read(4)dcrxx
      read(4)dcrw
      read(4)dcrwxx
      read(4)dcrxxxx
      close(4)
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
      I22=I22*alphaSTD
      I13=I13*alpha
      I23=I23*alpha
      I33=I33*alpha

      if (mod(Ngrid,2).eq.0) then
         nyq=float(Ngrid/2)
      else
         nyq=float((Ngrid-1)/2)
      endif   
      write(*,*)'random, data shot noise=',I12/I22,I12d/I22
      allocate(avgk(Nbins),avgPg0(Nbins),avgPr0(Nbins),co(Nbins))
      allocate(avgPg2fast(Nbins),avgPr2fast(Nbins))
      allocate(avgP00w(Nbins),avgP02w(Nbins),avgP20w(Nbins))
      allocate(avgPg4fft15(Nbins),avgPg0T(Nbins),avgPg4fast(Nbins))
      allocate(avgRWr0(Nbins),avgIWr0(Nbins))
      do 10 i=1,Nbins
         avgk(i)=0.d0
         avgPg0(i)=0.d0
         avgPg0T(i)=0.d0
         avgP00w(i)=0.d0
         avgP02w(i)=0.d0
         avgP20w(i)=0.d0
         avgPr0(i)=0.d0
         co(i)=0.d0
         avgPg2fast(i)=0.d0
         avgPr2fast(i)=0.d0
         avgPg4fft15(i)=0.d0
         avgPg4fast(i)=0.d0
         avgRWr0(i)=0.d0
         avgIWr0(i)=0.d0
 10   continue
      do iz=1,Ngrid
         do iy=1,Ngrid
            do ix=1,Ngrid/2+1
               dcg(ix,iy,iz)=dcg(ix,iy,iz) -alpha*dcr(ix,iy,iz)
               dcgxx(ix,iy,iz)=dcgxx(ix,iy,iz)-alpha*dcrxx(ix,iy,iz)
               dcgw(ix,iy,iz)=dcgw(ix,iy,iz)+alpha**2*dcrw(ix,iy,iz)
            dcgwxx(ix,iy,iz)=dcgwxx(ix,iy,iz)-alpha*dcrwxx(ix,iy,iz)
            dcgxxxx(ix,iy,iz)=dcgxxxx(ix,iy,iz)-alpha*dcrxxxx(ix,iy,iz)
            enddo
         enddo
      enddo
      cnorm=alpha*dcr(1,1,1)
      write(*,*)'window normalization=',cnorm
      
      do 100 iz=1,Ngrid
         ikz=mod(iz+Ngrid/2-2,Ngrid)-Ngrid/2+1
         icz=mod(Ngrid+1-iz,Ngrid)+1
         do 100 iy=1,Ngrid
            iky=mod(iy+Ngrid/2-2,Ngrid)-Ngrid/2+1
            icy=mod(Ngrid+1-iy,Ngrid)+1
            do 100 ix=1,Ngrid
               ikx=mod(ix+Ngrid/2-2,Ngrid)-Ngrid/2+1
               icx=mod(Ngrid+1-ix,Ngrid)+1
               rk=sqrt(float(ikx**2+iky**2+ikz**2))
               imk=nint(Nbins*rk/nyq)
               if(imk.le.Nbins .and. imk.ne.0)then
                  co(imk)=co(imk)+1.d0
                  avgk(imk)=avgk(imk)+dble(rk)                  
                  if (ix.le.Ngrid/2+1) then 
                     ct=dcg(ix,iy,iz)
                     ctQ=dcgxx(ix,iy,iz)
                     ctx=dcgw(ix,iy,iz)
                     ctP=dcgxxxx(ix,iy,iz)
                     cty=dcgwxx(ix,iy,iz)
                  else !use cc (kx changes sign!)
                     ct=dcg(icx,icy,icz)
                     ctQ=dcgxx(icx,icy,icz)
                     ctx=dcgw(icx,icy,icz)
                     ctP=dcgxxxx(icx,icy,icz)
                     cty=dcgwxx(icx,icy,icz)
                  endif
                  pk=(cabs(ct))**2
                  pq=ctQ*conjg(ct) !converting to real part gives the same
                  px=ctx*conjg(ct)
                  pp=ctP*conjg(ct)
                  p4=(cabs(ctQ))**2
                  py=cty*conjg(ct)
                  pz=ctx*conjg(ctQ)
                  avgPg0(imk)      = avgPg0(imk)      + dble(pk)
                  avgPg2fast(imk)  = avgPg2fast(imk)  + dble(pq)
                  avgP00w(imk)     = avgP00w(imk)     + dble(px)
                  avgP02w(imk)     = avgP02w(imk)     + dble(py)
                  avgP20w(imk)     = avgP20w(imk)     + dble(pz)
                  avgPg4fft15(imk) = avgPg4fft15(imk) + dble(pp)
                  avgPg4fast(imk)  = avgPg4fast(imk)  + dble(p4)
                  if (ix.le.Ngrid/2+1) then 
                     ct=alpha*dcr(ix,iy,iz)
                     ctQ=alpha*dcrxx(ix,iy,iz) 
                  else !use cc
                     ct=alpha*dcr(icx,icy,icz)
                     ctQ=alpha*dcrxx(icx,icy,icz)
                  endif
                  pk=(cabs(ct))**2
                  pq=ctQ*conjg(ct) !converting to real part gives the same
                  avgPr0(imk)     = avgPr0(imk)     + dble(pk)
                  avgPr2fast(imk) = avgPr2fast(imk) + dble(pq) 
                  avgRWr0(imk)= avgRWr0(imk)+ dble(real(ct))**2
                  avgIWr0(imk)= avgIWr0(imk)+ dble(imag(ct))**2
               end if
 100  continue


      akfun=6.28319/akfun
      open(4,file=filepower,status='unknown',form='formatted')
      do 110 Ibin=1,Nbins
         if(co(Ibin).gt.0.)then
            avgPg4fft15(Ibin)=avgPg4fft15(Ibin)/co(Ibin)/dble(I22) 
            avgPg4fast(Ibin)=avgPg4fast(Ibin)/co(Ibin)/dble(I22) 
            avgPg2fast(Ibin)=avgPg2fast(Ibin)/co(Ibin)/dble(I22) 
            avgk(Ibin)=avgk(Ibin)/co(Ibin)*dble(akfun)
            avgPg4fast(Ibin)= 0.7d0*avgPg4fast(Ibin)-avgPg2fast(Ibin)
     &       -3.5d0*avgPg0(Ibin)/co(Ibin)/dble(I22)
            avgP00w(Ibin)=avgP00w(Ibin)/dble(I23)*dble(I22)/avgPg0(Ibin) !in units of expected
            avgPg0T(Ibin)=(avgPg0(Ibin)/co(Ibin)-
     &      (dble(I12d)+dble(alpha*I12)) )/dble(I22) !"True" shot noise
            avgPg0(Ibin)=(avgPg0(Ibin)/co(Ibin)-
     &      dble(1.+alpha)*dble(I12))/dble(I22) !expected shot noise
            avgPr0(Ibin)=(avgPr0(Ibin)/co(Ibin)-
     &      dble(alpha)*dble(I12))/dble(I22) !monopole window sq
            avgPr0(Ibin)=avgPr0(Ibin)*dble(I22)/dble(cnorm**2)
            avgPr2fast(Ibin)=avgPr2fast(Ibin)/co(Ibin)/dble(I22) !quadrupole window sq
            avgRWr0(Ibin)=avgRWr0(Ibin)/co(Ibin)/dble(cnorm)**2
            avgIWr0(Ibin)=avgIWr0(Ibin)/co(Ibin)/dble(cnorm)**2
            avgPr2fast(Ibin)=avgPr2fast(Ibin)*dble(I22)/dble(cnorm**2)
            avgP02w(Ibin)=avgP02w(Ibin)/co(Ibin)/dble(I23)/
     &       avgPg2fast(Ibin) !in units of expected
            avgP20w(Ibin)=avgP20w(Ibin)/co(Ibin)/dble(I23)/
     &       avgPg2fast(Ibin) !in units of expected
      write(4,1015) avgk(Ibin),avgPg0T(Ibin),avgPg2fast(Ibin),
     & avgPg4fast(Ibin),avgPg4fft15(Ibin),avgPg0(Ibin),avgPr0(Ibin),
     & avgPr2fast(Ibin),avgRWr0(Ibin),avgIWr0(Ibin),avgP00w(Ibin),
     & avgP02w(Ibin),avgP20w(Ibin),co(Ibin)
      end if
 110  continue
      close(4)
 1015 format(2x,17e16.6)
      
      stop
      end
