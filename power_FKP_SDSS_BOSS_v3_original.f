! ifort -fast -o power_FKP_SDSS_BOSS_v3.exe power_FKP_SDSS_BOSS_v3.f
      implicit none ! change since v2: includes new weighted field 
      integer Ngrid,ix,iy,iz,Nbins,nyq,ikx,iky,ikz,imk,i,Ibin,Nr,Ng
      integer icz,icy,icx,Ngrid2
      complex, allocatable :: dcg(:,:,:),dcr(:,:,:)
      complex, allocatable :: dcgw(:,:,:),dcrw(:,:,:)
      complex, allocatable :: dcgxx(:,:,:),dcrxx(:,:,:)
      real*8, allocatable :: avgk(:),avgPg(:),avgPr(:),co(:),dk(:)
      real*8, allocatable :: avgPx(:)
      real*8, allocatable :: avgPg2Y(:),avgPr2Y(:),avgI12Q(:),avgPgT(:)
      character filecoef*200,filecoefr*200,filepower*200,Rboxstr*200
      character Nbinstr*200
      real akfun,I10,I12,I22,I13,I23,I33,P0,alpha,P0m,rk
      real I12xx,I12yy,I12zz,I12xy,I12yz,I12zx, alphaSTD
      real I10d,I12d,I22d,I13d,I23d,I33d
      real pk,Ngsys,pq,thetaobs,phiobs,Nrsys,Ngsyscomp,Ngsystot,px
      real xminD,xmaxD,yminD,ymaxD,zminD,zmaxD
      real xminR,xmaxR,yminR,ymaxR,zminR,zmaxR
      real kxh,kyh,kzh,Nrsyscomp,Nrsystot
      complex ct,ctQ,ctx
      
c      write(*,*) 'Random Fourier file :'
c      read(*,'(a)') filecoefr
      call getarg(1,filecoefr)
c      write(*,*) 'LSS/Mock Fourier file :'
c      read(*,'(a)') filecoef
      call getarg(2,filecoef)
c      write(*,*) 'INPUT Power Spectrum file :'
c      read(*,'(a)') filepower
      call getarg(3,filepower)
c      write(*,*)'Survey box size (Mpc/h)'
c      read(*,*) akfun      
      call getarg(4,Rboxstr)
      read(Rboxstr,*)akfun
c      write(*,*)'Number of bins'
      call getarg(5,Nbinstr)
      read(Nbinstr,*)Nbins
    
      open(unit=4,file=filecoef,status='old',form='unformatted')
      read(4)Ngrid
      
c      Nbins=Ngrid/2
      allocate(dcg(Ngrid/2+1,Ngrid,Ngrid),dcr(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcgw(Ngrid/2+1,Ngrid,Ngrid),dcrw(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcgxx(Ngrid/2+1,Ngrid,Ngrid))
      allocate(dcrxx(Ngrid/2+1,Ngrid,Ngrid))

      read(4)dcg
      read(4)I10d,I12d,I22d,I13d,I23d,I33d
      read(4)P0m,Ng,Ngsys,Ngsyscomp,Ngsystot
      read(4)xminD,xmaxD,yminD,ymaxD,zminD,zmaxD
      read(4)dcgxx
      read(4)dcgw
      close(4)
c      write(*,*)P0m,Ng,Ngsys
c      write(*,*)I10,I12,I22,I13,I23,I33
      open(unit=4,file=filecoefr,status='old',form='unformatted')
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
      read(4)dcrw
c      write(*,*)P0,Nr
c      write(*,*)I10,I12,I22,I13,I23,I33
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
      write(*,*)'galaxies',Ng,Ngsys,Ngsystot,Ngsyscomp
      write(*,*)'randoms',Nr,Nrsys,Nrsystot,Nrsyscomp
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
      I12xx=I12xx*alpha
      I12yy=I12yy*alpha
      I12zz=I12zz*alpha
      I12xy=I12xy*alpha
      I12yz=I12yz*alpha
      I12zx=I12zx*alpha


      nyq=float(Ngrid/2)
      write(*,*)'random, data shot noise=',I12/I22,I12d/I22
      allocate(avgk(Nbins),avgPg(Nbins),avgPr(Nbins),co(Nbins))
      allocate(dk(Nbins),avgPx(Nbins),avgPg2Y(Nbins),avgPr2Y(Nbins))
      allocate(avgI12Q(Nbins),avgPgT(Nbins))
      do 10 i=1,Nbins
         avgk(i)=0.d0
         avgPg(i)=0.d0
         avgPgT(i)=0.d0
         avgPx(i)=0.d0
         avgPr(i)=0.d0
         co(i)=0.d0
         avgPg2Y(i)=0.d0
         avgPr2Y(i)=0.d0
         avgI12Q(i)=0.d0
 10   continue
      do iz=1,Ngrid
         do iy=1,Ngrid
            do ix=1,Ngrid/2+1
               dcg(ix,iy,iz)=dcg(ix,iy,iz) -alpha*dcr(ix,iy,iz)
               dcgxx(ix,iy,iz)=dcgxx(ix,iy,iz)-alpha*dcrxx(ix,iy,iz)
               dcgw(ix,iy,iz)=dcgw(ix,iy,iz)+alpha**2*dcrw(ix,iy,iz)
            enddo
         enddo
      enddo

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
                  kxh=float(ikx)/rk !unit vectors
                  kyh=float(iky)/rk
                  kzh=float(ikz)/rk

                  co(imk)=co(imk)+1.d0
                  avgk(imk)=avgk(imk)+dble(rk)                  
                  if (ix.le.Ngrid/2+1) then 
                     ct=dcg(ix,iy,iz)
                     ctQ=dcgxx(ix,iy,iz)
                     ctx=dcgw(ix,iy,iz)
                  else !use cc (kx changes sign!)
                     ct=dcg(icx,icy,icz)
                     ctQ=dcgxx(icx,icy,icz)
                     ctx=dcgw(icx,icy,icz)
                  endif
                  pk=(cabs(ct))**2
                  pq=ctQ*conjg(ct) !converting to real part gives the same
                  px=ctx*conjg(ct)
                  avgPg(imk)=avgPg(imk)+dble(pk)
                  avgPx(imk)=avgPx(imk)+dble(px)
                  avgPg2Y(imk)=avgPg2Y(imk)+dble(pq)
                  avgI12Q(imk)=avgI12Q(imk)
     &             +dble(kxh**2*I12xx+kyh**2*I12yy+kzh**2*I12zz)
     &             +2.d0*dble(kxh*kyh*I12xy+kyh*kzh*I12yz+kzh*kxh*I12zx)
                  if (ix.le.Ngrid/2+1) then 
                     ct=alpha*dcr(ix,iy,iz)
                     ctQ=alpha*dcrxx(ix,iy,iz)
                  else !use cc
                     ct=alpha*dcr(icx,icy,icz)
                     ctQ=alpha*dcrxx(icx,icy,icz)
                  endif
                  pk=(cabs(ct))**2
                  pq=ctQ*conjg(ct) !converting to real part gives the same
                  avgPr(imk)=avgPr(imk)+dble(pk)
                  avgPr2Y(imk)=avgPr2Y(imk)+dble(pq) 
               end if
 100  continue


      akfun=6.28319/akfun
      open(4,file=filepower,status='unknown',form='formatted')
      do 110 Ibin=1,Nbins
         if(co(Ibin).gt.0.d0)then
            avgk(Ibin)=avgk(Ibin)/co(Ibin)*dble(akfun)
            avgPgT(Ibin)=(avgPg(Ibin)/co(Ibin)-     !"True" shot noise corrected power
     &      (dble(I12d)+dble(alpha*I12)) )/dble(I22) 
c            write(*,*)avgPg(Ibin)/co(Ibin),avgPr(Ibin)/co(Ibin),
c     &       dble(alpha)*dble(I12) 
            avgPg(Ibin)=(avgPg(Ibin)/co(Ibin)-  ! <shot noise> corrected power
     &      dble(1.+alpha)*dble(I12))/dble(I22) 
            avgPx(Ibin)=avgPx(Ibin)/co(Ibin)/dble(I23)/avgPg(Ibin) !in units of expected
            avgPr(Ibin)=(avgPr(Ibin)/co(Ibin)-
     &      dble(alpha)*dble(I12))/dble(I22) 
            dk(Ibin)=avgPg(Ibin)*avgk(Ibin)**3 /19.7392d0 
            avgPg2Y(Ibin)=avgPg2Y(Ibin)/co(Ibin)/dble(I22) 
            avgPr2Y(Ibin)=avgPr2Y(Ibin)/co(Ibin)/dble(I22)
            avgI12Q(Ibin)=7.5d0*avgI12Q(Ibin)/co(Ibin)-2.5d0*dble(I12)
            avgI12Q(Ibin)=dble(1.+alpha)*avgI12Q(Ibin)/dble(I22)
            avgPg2Y(Ibin)=avgPg2Y(Ibin)-avgI12Q(Ibin) !shot noise correction
      write(4,1015) avgk(Ibin),avgPg(Ibin),avgPg2Y(Ibin),avgPr(Ibin),
     &   avgPr2Y(Ibin),avgPgT(Ibin),avgPx(Ibin),co(Ibin)
      end if
 110  continue
      close(4)
 1015 format(2x,17e16.6)
      
      stop
      end
