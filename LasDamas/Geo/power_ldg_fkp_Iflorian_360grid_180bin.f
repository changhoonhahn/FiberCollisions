c PowerSpectrum Code for LasDamasGeo using Florian's normalization and shot noise term calculation 
c alpha = N_(gal,sys)/N_(random)
c I values comptued using Beutler et al. 2013
      implicit none 
      integer Ngrid,ix,iy,iz,Nbins,nyq,iky,ikz,imk,i,Ibin,Ng,Nr
      parameter(Ngrid=360,Nbins=180)
      complex, allocatable :: dcg(:,:,:),dcr(:,:,:)
c      complex dcg(Ngrid/2+1,Ngrid,Ngrid),dcr(Ngrid/2+1,Ngrid,Ngrid)
      real avgk(Nbins),avgPg(Nbins),avgPr(Nbins),co(Nbins),rk,dk(Nbins)
      real avgPg2(Nbins),avgPr2(Nbins),avgPg4(Nbins),avgPr4(Nbins)
      character filecoef*200,filecoefr*200,filepower*200
      character lssfname*200,randfname*200,powerfname*200,sscale*200
      real akfun,P0,alpha,P0m,wsys, wsysr
      real gI10,gI12,gI22,gI13,gI23,gI33
      real rI10,rI12,rI22,rI13,rI23,rI33
      real cot1,coga,Le2,Le4,pk
      complex ct
      
      call getarg(1,lssfname)
      call getarg(2,randfname)
      call getarg(3,powerfname)
      call getarg(4,sscale)
      read(sscale,*) akfun
      allocate(dcg(Ngrid/2+1,Ngrid,Ngrid),dcr(Ngrid/2+1,Ngrid,Ngrid))

      open(unit=4,file=lssfname,status='old',form='unformatted')
      read(4)dcg
      read(4)gI10,gI12,gI22,gI13,gI23,gI33
      read(4)P0m,Ng,wsys
      close(4)
      write(*,*) 'gIs',gI10,gI12,gI22,gI13,gI23,gI33
      open(unit=4,file=randfname,status='old',form='unformatted')
      read(4)dcr
      read(4)rI10,rI12,rI22,rI13,rI23,rI33
      read(4)P0,Nr,wsysr
      close(4)
      write(*,*) 'rIs',rI10,rI12,rI22,rI13,rI23,rI33
      if (P0m.ne.P0) then
         write(*,*)'P0s do not match'
         stop
      endif
      WRITE(*,*) "Ngrid=",Ngrid,"Box=",akfun,"P0=",P0
      WRITE(*,*) "Ng=",Ng,"Ng,sys=",wsys,"Nr=",Nr,"Nr,sys=",wsysr
      alpha=wsys/float(Nr) !now scale random integrals by alpha
      nyq=float(Ngrid/2)
      do 10 i=1,Nbins
         avgk(i)=0.
         avgPg(i)=0.
         avgPr(i)=0.
         avgPg2(i)=0.
         avgPr2(i)=0.
         avgPg4(i)=0.
         avgPr4(i)=0.
         co(i)=0.
 10   continue
      do iz=1,Ngrid
         do iy=1,Ngrid
            do ix=1,Ngrid/2+1
               dcg(ix,iy,iz)=dcg(ix,iy,iz)-alpha*dcr(ix,iy,iz)
            enddo
         enddo
      enddo
      do 100 iz=1,Ngrid
         ikz=mod(iz+Ngrid/2-2,Ngrid)-Ngrid/2+1
         do 100 iy=1,Ngrid
            iky=mod(iy+Ngrid/2-2,Ngrid)-Ngrid/2+1
            do 100 ix=1,Ngrid/2+1
               rk=sqrt(real((ix-1)**2+iky**2+ikz**2))
               imk=nint(Nbins*rk/nyq)
               if(imk.le.Nbins .and. imk.ne.0)then
                  cot1=real(ikz)/rk
                  coga=cot1 !use z-direction as redshift mapping
                  Le2=-5.e-1+1.5e0*coga**2
                  Le4=3.75e-1-3.75e0*coga**2+4.375e0*coga**4
                  co(imk)=co(imk)+1.
                  avgk(imk)=avgk(imk)+rk                  
                  ct=dcg(ix,iy,iz)
                  pk=(cabs(ct))**2
                  avgPg(imk)=avgPg(imk)+pk
                  avgPg2(imk)=avgPg2(imk)+pk*5.*Le2
                  avgPg4(imk)=avgPg4(imk)+pk*9.*Le4                  
                  ct=alpha*dcr(ix,iy,iz)
                  pk=(cabs(ct))**2
                  avgPr(imk)=avgPr(imk)+pk
                  avgPr2(imk)=avgPr2(imk)+pk*5.*Le2
                  avgPr4(imk)=avgPr4(imk)+pk*9.*Le4                  
               end if
 100  continue
      akfun=6.28319/akfun
      write(*,*) 'gI12',gI12
      write(*,*) 'alpha*rI12',alpha*rI12
      write(*,*) 'alpha*rI12 - gI12',alpha*rI12-gI12
      write(*,*) 'alpha*rI22 - gI22',alpha*rI22-gI22
      write(*,*) (1.+alpha)*alpha*rI12/(alpha*rI22)-
     &   (gI12+alpha**2*rI12)/gI22
      open(4,file=powerfname,status='unknown',form='formatted')
      do 110 Ibin=1,Nbins
         if(co(Ibin).gt.0.)then
            avgk(Ibin)=avgk(Ibin)/co(Ibin)*akfun
c            avgPg(Ibin)=(avgPg(Ibin)/co(Ibin)-(1.+alpha)*I12)/I22
            avgPg(Ibin)=(avgPg(Ibin)/co(Ibin)-(gI12+alpha**2*rI12))
     &   /(alpha*rI22) 
            avgPr(Ibin)=avgPr(Ibin)/co(Ibin) 
            avgPg2(Ibin)=avgPg2(Ibin)/co(Ibin)/(alpha*rI22)
            avgPr2(Ibin)=avgPr2(Ibin)/co(Ibin)/(alpha*rI22)
            avgPg4(Ibin)=avgPg4(Ibin)/co(Ibin)/(alpha*rI22)
            avgPr4(Ibin)=avgPr4(Ibin)/co(Ibin)/(alpha*rI22)
            dk(Ibin)=avgPg(Ibin)*avgk(Ibin)**3 /19.7392 
      write(4,1015) avgk(Ibin),avgPg(Ibin),avgPg2(Ibin),avgPg4(Ibin),
     &   avgPr(Ibin),avgPr2(Ibin),avgPr4(Ibin),dk(Ibin),co(Ibin)
      end if
 110  continue
      close(4)
 1015 format(2x,9e16.6)
      
      stop
      end
