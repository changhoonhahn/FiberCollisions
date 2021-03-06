      implicit none 
      integer Ngrid,ix,iy,iz,Nbins,nyq,iky,ikz,imk,i,Ibin,Ng,Nr,ikx
      integer icx,icy,icz
      parameter(Ngrid=960,Nbins=480)
      complex, allocatable :: dcg(:,:,:),dcr(:,:,:)
c      complex dcg(Ngrid/2+1,Ngrid,Ngrid),dcr(Ngrid/2+1,Ngrid,Ngrid)
      real avgk(Nbins),avgPg(Nbins),avgPr(Nbins),co(Nbins),rk,dk(Nbins)
      real avgPg2(Nbins),avgPr2(Nbins),avgPg4(Nbins),avgPr4(Nbins)
      character filecoef*200,filecoefr*200,filepower*200
      character randfft*200,lssfft*200,powername*200,akfunstr*200
      real akfun,I10,I12,I22,I13,I23,I33,P0,alpha,P0m,rwsys
      real cot1,coga,Le2,Le4,pk,wsys
      complex ct
      
      write(*,*) 'Random Fourier file :'
      CALL GETARG(1,randfft)
      write(*,*) 'LSS/Mock Fourier file :'
      CALL GETARG(2,lssfft)
      write(*,*) 'INPUT Power Spectrum file :'
      CALL GETARG(3,powername)
      write(*,*)'Survey scale (Mpc/h)'
      CALL GETARG(4,akfunstr)
      READ(akfunstr,*) akfun
      allocate(dcg(Ngrid/2+1,Ngrid,Ngrid),dcr(Ngrid/2+1,Ngrid,Ngrid))

      open(unit=4,file=lssfft,status='old',form='unformatted')
      read(4)dcg
      read(4)P0m,Ng,wsys 
      close(4)
      open(unit=4,file=randfft,status='old',form='unformatted')
      read(4)dcr
      read(4)I10,I12,I22,I13,I23,I33
      read(4)P0,Nr,rwsys
      close(4)
      if (P0m.ne.P0) then
         write(*,*)'P0s do not match'
         stop
      endif
      !alpha=float(Ng)/float(Nr) !now scale random integrals by alpha
      alpha=wsys/rwsys
      WRITE(*,*) 'Ngalsys=',wsys,'Ngal=',Ng,'ratio=',wsys/float(Ng)
      WRITE(*,*) 'Alpha=',alpha
      I10=I10*alpha
      I12=I12*alpha
      I22=I22*alpha
      I13=I13*alpha
      I23=I23*alpha
      I33=I33*alpha
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
         icz=mod(Ngrid+1-iz,Ngrid)+1
         do 100 iy=1,Ngrid
            iky=mod(iy+Ngrid/2-2,Ngrid)-Ngrid/2+1
            icy=mod(Ngrid+1-iy,Ngrid)+1
            do 100 ix=1,Ngrid
               ikx=mod(ix+Ngrid/2-2,Ngrid)-Ngrid/2+1
               icx=mod(Ngrid+1-ix,Ngrid)+1
               rk=sqrt(real(ikx**2+iky**2+ikz**2))
               imk=nint(Nbins*rk/nyq)
               if(imk.le.Nbins .and. imk.ne.0)then
                  cot1=real(ikz)/rk
                  coga=cot1 !use z-direction as redshift mapping
                  Le2=-5.e-1+1.5e0*coga**2
                  Le4=3.75e-1-3.75e0*coga**2+4.375e0*coga**4
                  co(imk)=co(imk)+1.
                  avgk(imk)=avgk(imk)+rk                  
                  if (ix.le.Ngrid/2+1) then 
                     ct=dcg(ix,iy,iz)
                  else !use cc
                     ct=dcg(icx,icy,icz)
                  endif
                  pk=(cabs(ct))**2
                  avgPg(imk)=avgPg(imk)+pk
                  avgPg2(imk)=avgPg2(imk)+pk*5.*Le2
                  avgPg4(imk)=avgPg4(imk)+pk*9.*Le4                  
                  if (ix.le.Ngrid/2+1) then 
                     ct=alpha*dcr(ix,iy,iz)
                  else !use cc
                     ct=alpha*dcr(icx,icy,icz)
                  endif
                  pk=(cabs(ct))**2
                  avgPr(imk)=avgPr(imk)+pk
                  avgPr2(imk)=avgPr2(imk)+pk*5.*Le2
                  avgPr4(imk)=avgPr4(imk)+pk*9.*Le4                  
               end if
 100  continue
      akfun=6.28319/akfun
      open(4,file=powername,status='unknown',form='formatted')
      do 110 Ibin=1,Nbins
         if(co(Ibin).gt.0.)then
            avgk(Ibin)=avgk(Ibin)/co(Ibin)*akfun
            avgPg(Ibin)=(avgPg(Ibin)/co(Ibin)-(1.+alpha)*I12)/I22
            avgPr(Ibin)=avgPr(Ibin)/co(Ibin) 
            avgPg2(Ibin)=avgPg2(Ibin)/co(Ibin)/I22
            avgPr2(Ibin)=avgPr2(Ibin)/co(Ibin)/I22 
            avgPg4(Ibin)=avgPg4(Ibin)/co(Ibin)/I22
            avgPr4(Ibin)=avgPr4(Ibin)/co(Ibin)/I22 
            dk(Ibin)=avgPg(Ibin)*avgk(Ibin)**3 /19.7392 
      write(4,1015) avgk(Ibin),avgPg(Ibin),avgPg2(Ibin),avgPg4(Ibin),
     &   avgPr(Ibin),avgPr2(Ibin),avgPr4(Ibin),dk(Ibin),co(Ibin)
      end if
 110  continue
      close(4)
 1015 format(2x,9e16.6)
      
      stop
      end
