! ifort  -fast -o FFT_FKP_BOSS_cic_il4_v5.exe FFT_FKP_BOSS_cic_il4_v5.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw -lm
! as v4 but builds arrays in turn to save memory (8 dc arrays at most are held in memory for IL4 data)
! use: FFT_FKP_BOSS_cic_il4_v5.exe idata Lbox Lm interpol idatran P0 ifc icomp inputfile {izbin} outputfile
      implicit none  
      integer Nran,i,iwr,Ngal,Nmax,n,kx,ky,kz,Lm,Ngrid,ix,iy,iz,Nsel
      integer Ng,Nr,iflag,ic,Nbin,interpol,Nariel,idata,iveto,icosmo
      integer*8 planf,plan_real
      real pi,wfc,wrf,xmin,xmax,ymin,ymax,zmin,zmax,cspeed
!      parameter(Nsel=181,Nmax=8*10**7,Nbin=201,Nsel2=10)
      parameter(Nsel=181,Nmax=3*10**8,Nbin=201)
      parameter(cspeed=299800.0,pi=3.141592654)
      integer grid,ifc,icomp,j,k,l,izbin
      dimension grid(3)
      integer, allocatable :: ig(:),ir(:)      
      real zbin(Nbin),dbin(Nbin),sec3(Nbin),zt,dum,gfrac,az2
      real cz,Om0,OL0,chi,nbar,Rbox,x0,y0,z0,wnoz,wcp
      real, allocatable :: nbg(:),nbr(:),rg(:,:),rr(:,:)
      real, allocatable :: wg(:),wr(:)
      real selfun(Nsel),z(Nsel),sec(Nsel),az,ra,dec,rad,numden,zend
      real w, nbb, xcm, ycm, zcm, rcm, bias, wboss, veto
      real thetaobs,phiobs
      real alpha,P0,nb,weight,ar,akf,Fr,Fi,Gr,Gi,wsys,wred,comp
      real*8 I10,I12,I22,I13,I23,I33
      real*8 compavg,Nrsys,Nrsyscomp,Ngsys,Ngsyscomp,Ngsystot,Nrsystot
      real rsq,xr,yr,zr
      real*8 I12xx,I12yy,I12zz,I12xy,I12yz,I12zx
      real kdotr,vol,xscale,rlow,rm(2)
      integer ikz,icz,iky,icy,ikx,icx
      real rk,kxh,kyh,kzh
      complex, allocatable :: dcg(:,:,:),dcr(:,:,:)
      complex, allocatable :: dcgxx(:,:,:),dcrxx(:,:,:)
      complex, allocatable :: dcgyy(:,:,:),dcryy(:,:,:)
      complex, allocatable :: dcgzz(:,:,:),dcrzz(:,:,:)
      complex, allocatable :: dcgxy(:,:,:),dcrxy(:,:,:)
      complex, allocatable :: dcgyz(:,:,:),dcryz(:,:,:)
      complex, allocatable :: dcgzx(:,:,:),dcrzx(:,:,:)
 
      complex, allocatable :: dcgw(:,:,:),dcrw(:,:,:)
      complex, allocatable :: dcgwxx(:,:,:),dcrwxx(:,:,:)
      complex, allocatable :: dcgwyy(:,:,:),dcrwyy(:,:,:)
      complex, allocatable :: dcgwzz(:,:,:),dcrwzz(:,:,:)
      complex, allocatable :: dcgwxy(:,:,:),dcrwxy(:,:,:)
      complex, allocatable :: dcgwyz(:,:,:),dcrwyz(:,:,:)
      complex, allocatable :: dcgwzx(:,:,:),dcrwzx(:,:,:)

      complex, allocatable :: dcgxxxx(:,:,:),dcrxxxx(:,:,:)
      complex, allocatable :: dcgyyyy(:,:,:),dcryyyy(:,:,:)
      complex, allocatable :: dcgzzzz(:,:,:),dcrzzzz(:,:,:)
      complex, allocatable :: dcgxxxy(:,:,:),dcrxxxy(:,:,:)
      complex, allocatable :: dcgxxxz(:,:,:),dcrxxxz(:,:,:)
      complex, allocatable :: dcgyyyx(:,:,:),dcryyyx(:,:,:)
      complex, allocatable :: dcgyyyz(:,:,:),dcryyyz(:,:,:)
      complex, allocatable :: dcgzzzx(:,:,:),dcrzzzx(:,:,:)
      complex, allocatable :: dcgzzzy(:,:,:),dcrzzzy(:,:,:)
      complex, allocatable :: dcgxxyy(:,:,:),dcrxxyy(:,:,:)
      complex, allocatable :: dcgxxzz(:,:,:),dcrxxzz(:,:,:)
      complex, allocatable :: dcgyyzz(:,:,:),dcryyzz(:,:,:)
      complex, allocatable :: dcgxxyz(:,:,:),dcrxxyz(:,:,:)
      complex, allocatable :: dcgyyxz(:,:,:),dcryyxz(:,:,:)
      complex, allocatable :: dcgzzxy(:,:,:),dcrzzxy(:,:,:)

      complex, allocatable :: dcgwxxxx(:,:,:),dcrwxxxx(:,:,:)
      complex, allocatable :: dcgwyyyy(:,:,:),dcrwyyyy(:,:,:)
      complex, allocatable :: dcgwzzzz(:,:,:),dcrwzzzz(:,:,:)
      complex, allocatable :: dcgwxxxy(:,:,:),dcrwxxxy(:,:,:)
      complex, allocatable :: dcgwxxxz(:,:,:),dcrwxxxz(:,:,:)
      complex, allocatable :: dcgwyyyx(:,:,:),dcrwyyyx(:,:,:)
      complex, allocatable :: dcgwyyyz(:,:,:),dcrwyyyz(:,:,:)
      complex, allocatable :: dcgwzzzx(:,:,:),dcrwzzzx(:,:,:)
      complex, allocatable :: dcgwzzzy(:,:,:),dcrwzzzy(:,:,:)
      complex, allocatable :: dcgwxxyy(:,:,:),dcrwxxyy(:,:,:)
      complex, allocatable :: dcgwxxzz(:,:,:),dcrwxxzz(:,:,:)
      complex, allocatable :: dcgwyyzz(:,:,:),dcrwyyzz(:,:,:)
      complex, allocatable :: dcgwxxyz(:,:,:),dcrwxxyz(:,:,:)
      complex, allocatable :: dcgwyyxz(:,:,:),dcrwyyxz(:,:,:)
      complex, allocatable :: dcgwzzxy(:,:,:),dcrwzzxy(:,:,:)

      character lssfile*200,randomfile*200,filecoef*200
      character dummy*200,fname*200,outname*200,Omstr*200
      character Rboxstr*200,Ngridstr*200,interpolstr*200,iflagstr*200
      character P0str*200,typestr*200,lssinfofile*200,cosmostr*200
      character*200 icompstr,ifcstr,izbinstr
      common /interpol/z,selfun,sec
      common /interp3/dbin,zbin,sec3
      common /radint/Om0,OL0
      external nbar,chi,PutIntoBox,assign2,fcomb
      include '/usr/local/src/fftw-2.1.5/fortran/fftw_f77.i'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! idata= 1 (BOSS data), 2(LasDamas), 3-4 (QPM mocks N-S)
      ! idata= 5 (LOWZ data), 6(PATCHY), 7 (Nseries)
      call getarg(1,typestr)
      read(typestr,*)idata 
    
      call getarg(2,cosmostr)
      read(cosmostr,*)icosmo
      call cosmology(idata,icosmo,Om0,OL0) !fiducial cosmology(0), survey cosmology(1)
      call nbtable(idata,Nsel,z,selfun,sec)
      call radialdist(idata,Nsel,z,Nbin,zbin,dbin,sec3)

      call getarg(3,Rboxstr)
      read(Rboxstr,*)Rbox !box side (survey origin will be at box center)
      xscale=RBox !actual side of Fourier Box
      Rbox=0.5*RBox !"radius of the box"
      rlow=-Rbox !lowest coordinate in any direction

      call getarg(4,Ngridstr)
      read(Ngridstr,*)Ngrid !FFT Grid Size
      Lm=Ngrid
      rm(1)=float(Lm)/xscale
      rm(2)=1.-rlow*rm(1)

      call getarg(5,interpolstr)
      read(interpolstr,*)interpol !2nd order cic (2) or 4th-order interlaced(4)?
      if (interpol.eq.2) then !CIC needs RtoC transform
         call rfftw3d_f77_create_plan(plan_real,Lm,Lm,Lm,
     &   FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE + FFTW_IN_PLACE)
      else
         grid(1) = Lm
         grid(2) = Lm
         grid(3) = Lm
         call fftwnd_f77_create_plan(planf,3,grid,FFTW_BACKWARD,
     $   FFTW_ESTIMATE + FFTW_IN_PLACE)
      endif

      call getarg(6,iflagstr) 
      read(iflagstr,*)iflag !mock (0) or random mock(1)?
      call getarg(7,P0str)
      read(P0str,*)P0 !FKP weight

      if (iflag.eq.0) then ! run on mock
      
         call getarg(8,lssfile) !Mock Survey File
         allocate(rg(3,Nmax),nbg(Nmax),ig(Nmax),wg(Nmax))
         call readmock(idata,lssfile,Nmax,nbg,wg,rg,P0,
     $        Ngal,Ngsys,Ngsyscomp,Ngsystot,compavg,
     $        xmin,xmax,ymin,ymax,zmin,zmax)
         
         call PutIntoBox(Ngal,rg,Rbox,ig,Ng,Nmax)

         gfrac=100. *float(Ng)/float(Ngal)
         write(*,*)'Number of Galaxies in Box=',Ng,gfrac,'percent'
         Ngsys=Ngsys*dble(Ng)/dble(Ngal) !in reality we have to do the Ngsys sum again!
         Ngsyscomp=Ngsyscomp*dble(Ng)/dble(Ngal) !in reality we have to do the Ngsys sum again!
         Ngsystot=Ngsystot*dble(Ng)/dble(Ngal) !in reality we have to do the Ngsys sum again!
         ! but it does not matter as we always go for 100% galaxies inside FFT box
         write(*,*)'upweighted-Galaxies in Box=',Ngsys
         write(*,*)'comp+upweighted-Galaxies in Box=',Ngsyscomp
         write(*,*)'comp+upweighted+FKP-Galaxies in Box=',Ngsystot
         compavg=compavg/dble(Ngal)
         write(*,*)'average comp=',compavg

         call computeIij(iflag,Nmax,nbg,wg,ig,Ngal,
     $                   I10,I12,I22,I13,I23,I33)

         if (interpol.eq.2) then !CIC
            allocate(dcg(Lm/2+1,Lm,Lm))
            allocate(dcgxx(Lm/2+1,Lm,Lm),dcgyy(Lm/2+1,Lm,Lm))
            allocate(dcgzz(Lm/2+1,Lm,Lm),dcgxy(Lm/2+1,Lm,Lm))
            allocate(dcgyz(Lm/2+1,Lm,Lm),dcgzx(Lm/2+1,Lm,Lm))

            allocate(dcgw(Lm/2+1,Lm,Lm))
            allocate(dcgwxx(Lm/2+1,Lm,Lm),dcgwyy(Lm/2+1,Lm,Lm))
            allocate(dcgwzz(Lm/2+1,Lm,Lm),dcgwxy(Lm/2+1,Lm,Lm))
            allocate(dcgwyz(Lm/2+1,Lm,Lm),dcgwzx(Lm/2+1,Lm,Lm))

            allocate(dcgxxxx(Lm/2+1,Lm,Lm))
            allocate(dcgyyyy(Lm/2+1,Lm,Lm),dcgzzzz(Lm/2+1,Lm,Lm))
            allocate(dcgxxxy(Lm/2+1,Lm,Lm),dcgxxxz(Lm/2+1,Lm,Lm))
            allocate(dcgyyyx(Lm/2+1,Lm,Lm),dcgyyyz(Lm/2+1,Lm,Lm))
            allocate(dcgzzzx(Lm/2+1,Lm,Lm),dcgzzzy(Lm/2+1,Lm,Lm))
            allocate(dcgxxyy(Lm/2+1,Lm,Lm),dcgxxzz(Lm/2+1,Lm,Lm))
            allocate(dcgyyzz(Lm/2+1,Lm,Lm),dcgxxyz(Lm/2+1,Lm,Lm))
            allocate(dcgyyxz(Lm/2+1,Lm,Lm),dcgzzxy(Lm/2+1,Lm,Lm))

            call assign_CIC(Ngal,rg,rm,Lm,dcg,wg,ig,0,0,0,0) 

            call assign_CIC(Ngal,rg,rm,Lm,dcgxx,wg,ig,1,1,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyy,wg,ig,2,2,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzz,wg,ig,3,3,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxy,wg,ig,1,2,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyz,wg,ig,2,3,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzx,wg,ig,3,1,0,0) 

            call assign_CIC(Ngal,rg,rm,Lm,dcgxxxx,wg,ig,1,1,1,1) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyyyy,wg,ig,2,2,2,2) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzzzz,wg,ig,3,3,3,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxxxy,wg,ig,1,1,1,2) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxxxz,wg,ig,1,1,1,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyyyx,wg,ig,2,2,2,1) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyyyz,wg,ig,2,2,2,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzzzx,wg,ig,3,3,3,1) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzzzy,wg,ig,3,3,3,2) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxxyy,wg,ig,1,1,2,2) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxxzz,wg,ig,1,1,3,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyyzz,wg,ig,2,2,3,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxxyz,wg,ig,1,1,2,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyyxz,wg,ig,2,2,1,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzzxy,wg,ig,3,3,1,2) 

            do i=1,Ngal
               wg(i)=wg(i)**2
            enddo   
            call assign_CIC(Ngal,rg,rm,Lm,dcgw,wg,ig,0,0,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgwxx,wg,ig,1,1,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgwyy,wg,ig,2,2,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgwzz,wg,ig,3,3,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgwxy,wg,ig,1,2,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgwyz,wg,ig,2,3,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgwzx,wg,ig,3,1,0,0) 

            call rfftwnd_f77_one_real_to_complex(plan_real,dcg,dcg)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgxx,dcgxx)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgyy,dcgyy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgzz,dcgzz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgxy,dcgxy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgyz,dcgyz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgzx,dcgzx)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgw,dcgw)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcgwxx,dcgwxx)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcgwyy,dcgwyy)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcgwzz,dcgwzz)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcgwxy,dcgwxy)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcgwyz,dcgwyz)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcgwzx,dcgwzx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgxxxx,dcgxxxx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgyyyy,dcgyyyy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgzzzz,dcgzzzz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgxxxy,dcgxxxy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgxxxz,dcgxxxz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgyyyx,dcgyyyx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgyyyz,dcgyyyz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgzzzx,dcgzzzx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgzzzy,dcgzzzy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgxxyy,dcgxxyy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgxxzz,dcgxxzz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgyyzz,dcgyyzz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgxxyz,dcgxxyz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgyyxz,dcgyyxz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcgzzxy,dcgzzxy)
            call correct(Lm,Lm,Lm,dcg) 
            call correct(Lm,Lm,Lm,dcgxx) 
            call correct(Lm,Lm,Lm,dcgyy) 
            call correct(Lm,Lm,Lm,dcgzz) 
            call correct(Lm,Lm,Lm,dcgxy) 
            call correct(Lm,Lm,Lm,dcgyz) 
            call correct(Lm,Lm,Lm,dcgzx) 
            call correct(Lm,Lm,Lm,dcgw) 
            call correct(Lm,Lm,Lm,dcgwxx) 
            call correct(Lm,Lm,Lm,dcgwyy) 
            call correct(Lm,Lm,Lm,dcgwzz) 
            call correct(Lm,Lm,Lm,dcgwxy) 
            call correct(Lm,Lm,Lm,dcgwyz) 
            call correct(Lm,Lm,Lm,dcgwzx) 
            call correct(Lm,Lm,Lm,dcgxxxx) 
            call correct(Lm,Lm,Lm,dcgyyyy) 
            call correct(Lm,Lm,Lm,dcgzzzz) 
            call correct(Lm,Lm,Lm,dcgxxxy) 
            call correct(Lm,Lm,Lm,dcgxxxz) 
            call correct(Lm,Lm,Lm,dcgyyyx) 
            call correct(Lm,Lm,Lm,dcgyyyz) 
            call correct(Lm,Lm,Lm,dcgzzzx) 
            call correct(Lm,Lm,Lm,dcgzzzy) 
            call correct(Lm,Lm,Lm,dcgxxyy) 
            call correct(Lm,Lm,Lm,dcgxxzz) 
            call correct(Lm,Lm,Lm,dcgyyzz) 
            call correct(Lm,Lm,Lm,dcgxxyz) 
            call correct(Lm,Lm,Lm,dcgyyxz) 
            call correct(Lm,Lm,Lm,dcgzzxy) 

            do 98 iz=1,Ngrid !build quadrupole
             ikz=mod(iz+Ngrid/2-2,Ngrid)-Ngrid/2+1
             do 98 iy=1,Ngrid
              iky=mod(iy+Ngrid/2-2,Ngrid)-Ngrid/2+1
              do 98 ix=1,Ngrid/2+1
               ikx=mod(ix+Ngrid/2-2,Ngrid)-Ngrid/2+1
               rk=sqrt(float(ikx**2+iky**2+ikz**2))
               if(rk.gt.0.)then
                  kxh=float(ikx)/rk !unit vectors
                  kyh=float(iky)/rk
                  kzh=float(ikz)/rk
                     dcgxx(ix,iy,iz)=7.5*(dcgxx(ix,iy,iz)*kxh**2 
     &                  +dcgyy(ix,iy,iz)*kyh**2
     &                  +dcgzz(ix,iy,iz)*kzh**2 
     &                  +2.*dcgxy(ix,iy,iz)*kxh*kyh
     &                  +2.*dcgyz(ix,iy,iz)*kyh*kzh
     &                  +2.*dcgzx(ix,iy,iz)*kzh*kxh)
     &                  -2.5*dcg(ix,iy,iz)  !quadrupole field: 5 delta2
                     dcgwxx(ix,iy,iz)=7.5*(dcgwxx(ix,iy,iz)*kxh**2 
     &                  +dcgwyy(ix,iy,iz)*kyh**2
     &                  +dcgwzz(ix,iy,iz)*kzh**2 
     &                  +2.*dcgwxy(ix,iy,iz)*kxh*kyh
     &                  +2.*dcgwyz(ix,iy,iz)*kyh*kzh
     &                  +2.*dcgwzx(ix,iy,iz)*kzh*kxh)
     &                  -2.5*dcgw(ix,iy,iz)  !quadrupole field: 5 delta2
               end if
 98         continue
      
            deallocate(dcgyy,dcgzz,dcgxy,dcgyz,dcgzx)
            deallocate(dcgwyy,dcgwzz,dcgwxy,dcgwyz,dcgwzx)

            do 99 iz=1,Ngrid !build hexadecapole
             ikz=mod(iz+Ngrid/2-2,Ngrid)-Ngrid/2+1
             do 99 iy=1,Ngrid
              iky=mod(iy+Ngrid/2-2,Ngrid)-Ngrid/2+1
              do 99 ix=1,Ngrid/2+1
               ikx=mod(ix+Ngrid/2-2,Ngrid)-Ngrid/2+1
               rk=sqrt(float(ikx**2+iky**2+ikz**2))
               if(rk.gt.0.)then
                  kxh=float(ikx)/rk !unit vectors
                  kyh=float(iky)/rk
                  kzh=float(ikz)/rk
                  dcgxxxx(ix,iy,iz)=35.*9./8.*(dcgxxxx(ix,iy,iz)*kxh**4 
     &              +dcgyyyy(ix,iy,iz)*kyh**4+dcgzzzz(ix,iy,iz)*kzh**4
     & +4.*dcgxxxy(ix,iy,iz)*kxh**3*kyh+4.*dcgxxxz(ix,iy,iz)*kxh**3*kzh 
     & +4.*dcgyyyx(ix,iy,iz)*kyh**3*kxh+4.*dcgyyyz(ix,iy,iz)*kyh**3*kzh 
     & +4.*dcgzzzx(ix,iy,iz)*kzh**3*kxh+4.*dcgzzzy(ix,iy,iz)*kzh**3*kyh 
     & +6.*dcgxxyy(ix,iy,iz)*kxh**2*kyh**2
     & +6.*dcgxxzz(ix,iy,iz)*kxh**2*kzh**2
     & +6.*dcgyyzz(ix,iy,iz)*kyh**2*kzh**2
     & +12.*dcgxxyz(ix,iy,iz)*kxh**2*kyh*kzh
     & +12.*dcgyyxz(ix,iy,iz)*kyh**2*kxh*kzh
     & +12.*dcgzzxy(ix,iy,iz)*kzh**2*kxh*kyh)
     & -9./2.*dcgxx(ix,iy,iz) -63./8.*dcg(ix,iy,iz)  !hexadecapole field: 9 delta4
               end if
 99         continue
      
            deallocate(dcgyyyy,dcgzzzz,dcgxxxy,dcgxxxz,dcgyyyx,dcgyyyz)
            deallocate(dcgzzzx,dcgzzzy,dcgxxyy,dcgxxzz,dcgyyzz,dcgxxyz)
            deallocate(dcgyyxz,dcgzzxy)

         else !4th-order interlaced

            allocate(dcg(Lm,Lm,Lm)) !do delta_g
            call assign(Ngal,rg,rm,Lm,dcg,P0,nbg,ig,wg,0,0,0,0)
            call fftwnd_f77_one(planf,dcg,dcg)      
            call fcomb(Lm,dcg,Ng)

            allocate(dcgxx(Lm,Lm,Lm),dcgyy(Lm,Lm,Lm),dcgzz(Lm,Lm,Lm)) !do 5 delta2_g (1st part)
            call FiveDelta2g_1(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg, 
     $      dcgxx,dcgyy,dcgzz,planf)
            deallocate(dcgyy,dcgzz) !done 5 delta2_g (1st part)

            allocate(dcgxy(Lm,Lm,Lm),dcgyz(Lm,Lm,Lm),dcgzx(Lm,Lm,Lm)) !do 5 delta2_g (2nd part)
            call FiveDelta2g_2(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg, 
     $      dcg,dcgxx,dcgxy,dcgyz,dcgzx,planf)
            deallocate(dcgxy,dcgyz,dcgzx) !done 5 delta2_g (2nd part)

            allocate(dcgxxxx(Lm,Lm,Lm),dcgyyyy(Lm,Lm,Lm)) !do 9 delta4_g (1st part)
            allocate(dcgzzzz(Lm,Lm,Lm),dcgxxxy(Lm,Lm,Lm))
            allocate(dcgxxxz(Lm,Lm,Lm))
            call NineDelta4g_1(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,
     $           dcgxxxx,dcgyyyy,dcgzzzz,dcgxxxy,dcgxxxz,planf)
            deallocate(dcgyyyy,dcgzzzz,dcgxxxy,dcgxxxz) !done 9 delta4_g (1st part)

            allocate(dcgyyyx(Lm,Lm,Lm))!do 9 delta4_g (2nd part)
            allocate(dcgyyyz(Lm,Lm,Lm),dcgzzzx(Lm,Lm,Lm))
            allocate(dcgzzzy(Lm,Lm,Lm),dcgxxyy(Lm,Lm,Lm))
            call NineDelta4g_2(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,
     $           dcgxxxx,dcgyyyx,dcgyyyz,dcgzzzx,dcgzzzy,dcgxxyy,planf)
            deallocate(dcgyyyx,dcgyyyz,dcgzzzx,dcgzzzy,dcgxxyy)!done 9 delta4_g (2nd part)

            allocate(dcgxxzz(Lm,Lm,Lm)) !do 9 delta4_g (3rd part)
            allocate(dcgyyzz(Lm,Lm,Lm),dcgxxyz(Lm,Lm,Lm))
            allocate(dcgyyxz(Lm,Lm,Lm),dcgzzxy(Lm,Lm,Lm))
            call NineDelta4g_3(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,dcg,
     $      dcgxx,dcgxxxx,dcgxxzz,dcgyyzz,dcgxxyz,dcgyyxz,dcgzzxy,planf)
            deallocate(dcgxxzz,dcgyyzz,dcgxxyz,dcgyyxz,dcgzzxy)!done 9 delta4_g (3rd part)

            do i=1,Ngal !build reweighted fields now
               wg(i)=wg(i)**2
            enddo   
            allocate(dcgw(Lm,Lm,Lm)) !do deltaw_g
            call assign(Ngal,rg,rm,Lm,dcgw,P0,nbg,ig,wg,0,0,0,0)
            call fftwnd_f77_one(planf,dcgw,dcgw)      
            call fcomb(Lm,dcgw,Ng)

            allocate(dcgwxx(Lm,Lm,Lm)) !do 5 delta2w_g (1st part)
            allocate(dcgwyy(Lm,Lm,Lm),dcgwzz(Lm,Lm,Lm))
            call FiveDelta2g_1(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg, 
     $      dcgwxx,dcgwyy,dcgwzz,planf)
            deallocate(dcgwyy,dcgwzz) !done 5 delta2w_g (1st part)

            allocate(dcgwxy(Lm,Lm,Lm)) !do 5 delta2w_g (2nd part)
            allocate(dcgwyz(Lm,Lm,Lm),dcgwzx(Lm,Lm,Lm))
            call FiveDelta2g_2(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg, 
     $      dcgw,dcgwxx,dcgwxy,dcgwyz,dcgwzx,planf)
            deallocate(dcgwxy,dcgwyz,dcgwzx) !done 5 delta2w_g (2nd part)

            allocate(dcgwxxxx(Lm,Lm,Lm),dcgwyyyy(Lm,Lm,Lm)) !do 9 delta4w_g (1st part)
            allocate(dcgwzzzz(Lm,Lm,Lm),dcgwxxxy(Lm,Lm,Lm))
            allocate(dcgwxxxz(Lm,Lm,Lm))
            call NineDelta4g_1(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,
     $           dcgwxxxx,dcgwyyyy,dcgwzzzz,dcgwxxxy,dcgwxxxz,planf)
            deallocate(dcgwyyyy,dcgwzzzz,dcgwxxxy,dcgwxxxz) !done 9 delta4w_g (1st part)

            allocate(dcgwyyyx(Lm,Lm,Lm))!do 9 delta4w_g (2nd part)
            allocate(dcgwyyyz(Lm,Lm,Lm),dcgwzzzx(Lm,Lm,Lm))
            allocate(dcgwzzzy(Lm,Lm,Lm),dcgwxxyy(Lm,Lm,Lm))
            call NineDelta4g_2(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,
     $           dcgwxxxx,dcgwyyyx,dcgwyyyz,dcgwzzzx,dcgwzzzy,
     $           dcgwxxyy,planf)
            deallocate(dcgwyyyx,dcgwyyyz,dcgwzzzx,dcgwzzzy,dcgwxxyy)!done 9 delta4w_g (2nd part)

            allocate(dcgwxxzz(Lm,Lm,Lm)) !do 9 delta4w_g (3rd part)
            allocate(dcgwyyzz(Lm,Lm,Lm),dcgwxxyz(Lm,Lm,Lm))
            allocate(dcgwyyxz(Lm,Lm,Lm),dcgwzzxy(Lm,Lm,Lm))
            call NineDelta4g_3(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,dcgw,
     $      dcgwxx,dcgwxxxx,dcgwxxzz,dcgwyyzz,dcgwxxyz,dcgwyyxz,
     $      dcgwzzxy,planf)
            deallocate(dcgwxxzz,dcgwyyzz,dcgwxxyz,dcgwyyxz,dcgwzzxy)!done 9 delta4w_g (3rd part)

         endif
         
         call getarg(9,filecoef) !Fourier file
         open(unit=4,file=filecoef,status='unknown',form='unformatted')
         write(4)Lm
         write(4)(((dcg(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)real(I10),real(I12),real(I22),real(I13),real(I23),
     &   real(I33) 
         write(4)P0,Ng,real(Ngsys),real(Ngsyscomp), real(Ngsystot)
         write(4)xmin,xmax,ymin,ymax,zmin,zmax
         write(4)(((dcgxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcgw(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcgwxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcgxxxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcgwxxxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         close(4)

       elseif (iflag.eq.1) then ! compute discretness integrals and FFT random mock

         call getarg(8,randomfile) !Random Survey File
         allocate(rr(3,Nmax),nbr(Nmax),ir(Nmax),wr(Nmax))

         call readrand(idata,randomfile,Nmax,nbr,wr,rr,P0,
     $           Nran,Nrsys,Nrsyscomp,Nrsystot,compavg,izbin,
     $           xmin,xmax,ymin,ymax,zmin,zmax,thetaobs,phiobs)
      
         call PutIntoBox(Nran,rr,Rbox,ir,Nr,Nmax)
         gfrac=100. *float(Nr)/float(Nran)
         write(*,*)'Number of Randoms in Box=',Nr,gfrac,'percent'
         Nrsys=Nrsys*dble(Nr)/dble(Nran) !scale in case not 100 % inside mask
         Nrsyscomp=Nrsyscomp*dble(Nr)/dble(Nran) !scale in case not 100 % inside mask
         Nrsystot=Nrsystot*dble(Nr)/dble(Nran) !scale in case not 100 % inside mask
         write(*,*)'upweighted randoms in Box=',Nrsys
         write(*,*)'comp+upweighted-randoms in Box=',Nrsyscomp
         write(*,*)'comp+upweighted+FKP-randoms in Box=',Nrsystot
         compavg=compavg/dble(Nran)
         write(*,*)'average comp=',compavg

         call computeIij(iflag,Nmax,nbr,wr,ir,Nran,
     $                   I10,I12,I22,I13,I23,I33)

         if (interpol.eq.2) then !CIC
            allocate(dcr(Lm/2+1,Lm,Lm),dcrw(Lm/2+1,Lm,Lm))
            allocate(dcrxx(Lm/2+1,Lm,Lm),dcryy(Lm/2+1,Lm,Lm))
            allocate(dcrzz(Lm/2+1,Lm,Lm),dcrxy(Lm/2+1,Lm,Lm))
            allocate(dcryz(Lm/2+1,Lm,Lm),dcrzx(Lm/2+1,Lm,Lm))
            allocate(dcrwxx(Lm/2+1,Lm,Lm),dcrwyy(Lm/2+1,Lm,Lm))
            allocate(dcrwzz(Lm/2+1,Lm,Lm),dcrwxy(Lm/2+1,Lm,Lm))
            allocate(dcrwyz(Lm/2+1,Lm,Lm),dcrwzx(Lm/2+1,Lm,Lm))
            allocate(dcrxxxx(Lm/2+1,Lm,Lm),dcryyyy(Lm/2+1,Lm,Lm))
            allocate(dcrzzzz(Lm/2+1,Lm,Lm),dcrxxxy(Lm/2+1,Lm,Lm))
            allocate(dcrxxxz(Lm/2+1,Lm,Lm),dcryyyx(Lm/2+1,Lm,Lm))
            allocate(dcryyyz(Lm/2+1,Lm,Lm),dcrzzzx(Lm/2+1,Lm,Lm))
            allocate(dcrzzzy(Lm/2+1,Lm,Lm),dcrxxyy(Lm/2+1,Lm,Lm))
            allocate(dcrxxzz(Lm/2+1,Lm,Lm),dcryyzz(Lm/2+1,Lm,Lm))
            allocate(dcrxxyz(Lm/2+1,Lm,Lm),dcryyxz(Lm/2+1,Lm,Lm))
            allocate(dcrzzxy(Lm/2+1,Lm,Lm))
            call assign_CIC(Nran,rr,rm,Lm,dcr,wr,ir,0,0,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxx,wr,ir,1,1,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcryy,wr,ir,2,2,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzz,wr,ir,3,3,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxy,wr,ir,1,2,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcryz,wr,ir,2,3,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzx,wr,ir,3,1,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxxxx,wr,ir,1,1,1,1) 
            call assign_CIC(Nran,rr,rm,Lm,dcryyyy,wr,ir,2,2,2,2) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzzzz,wr,ir,3,3,3,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxxxy,wr,ir,1,1,1,2) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxxxz,wr,ir,1,1,1,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcryyyx,wr,ir,2,2,2,1) 
            call assign_CIC(Nran,rr,rm,Lm,dcryyyz,wr,ir,2,2,2,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzzzx,wr,ir,3,3,3,1) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzzzy,wr,ir,3,3,3,2) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxxyy,wr,ir,1,1,2,2) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxxzz,wr,ir,1,1,3,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcryyzz,wr,ir,2,2,3,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxxyz,wr,ir,1,1,2,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcryyxz,wr,ir,2,2,1,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzzxy,wr,ir,3,3,1,2) 
            do i=1,Nran
               wr(i)=wr(i)**2
            enddo   
            call assign_CIC(Nran,rr,rm,Lm,dcrw,wr,ir,0,0,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrwxx,wr,ir,1,1,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrwyy,wr,ir,2,2,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrwzz,wr,ir,3,3,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrwxy,wr,ir,1,2,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrwyz,wr,ir,2,3,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrwzx,wr,ir,3,1,0,0) 
            call rfftwnd_f77_one_real_to_complex(plan_real,dcr,dcr)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrxx,dcrxx)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcryy,dcryy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrzz,dcrzz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrxy,dcrxy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcryz,dcryz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrzx,dcrzx)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrw,dcrw)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcrwxx,dcrwxx)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcrwyy,dcrwyy)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcrwzz,dcrwzz)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcrwxy,dcrwxy)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcrwyz,dcrwyz)
          call rfftwnd_f77_one_real_to_complex(plan_real,dcrwzx,dcrwzx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrxxxx,dcrxxxx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcryyyy,dcryyyy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrzzzz,dcrzzzz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrxxxy,dcrxxxy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrxxxz,dcrxxxz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcryyyx,dcryyyx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcryyyz,dcryyyz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrzzzx,dcrzzzx)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrzzzy,dcrzzzy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrxxyy,dcrxxyy)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrxxzz,dcrxxzz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcryyzz,dcryyzz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrxxyz,dcrxxyz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcryyxz,dcryyxz)
        call rfftwnd_f77_one_real_to_complex(plan_real,dcrzzxy,dcrzzxy)
            call correct(Lm,Lm,Lm,dcr) 
            call correct(Lm,Lm,Lm,dcrxx) 
            call correct(Lm,Lm,Lm,dcryy) 
            call correct(Lm,Lm,Lm,dcrzz) 
            call correct(Lm,Lm,Lm,dcrxy) 
            call correct(Lm,Lm,Lm,dcryz) 
            call correct(Lm,Lm,Lm,dcrzx) 
            call correct(Lm,Lm,Lm,dcrw) 
            call correct(Lm,Lm,Lm,dcrwxx) 
            call correct(Lm,Lm,Lm,dcrwyy) 
            call correct(Lm,Lm,Lm,dcrwzz) 
            call correct(Lm,Lm,Lm,dcrwxy) 
            call correct(Lm,Lm,Lm,dcrwyz) 
            call correct(Lm,Lm,Lm,dcrwzx) 
            call correct(Lm,Lm,Lm,dcrxxxx) 
            call correct(Lm,Lm,Lm,dcryyyy) 
            call correct(Lm,Lm,Lm,dcrzzzz) 
            call correct(Lm,Lm,Lm,dcrxxxy) 
            call correct(Lm,Lm,Lm,dcrxxxz) 
            call correct(Lm,Lm,Lm,dcryyyx) 
            call correct(Lm,Lm,Lm,dcryyyz) 
            call correct(Lm,Lm,Lm,dcrzzzx) 
            call correct(Lm,Lm,Lm,dcrzzzy) 
            call correct(Lm,Lm,Lm,dcrxxyy) 
            call correct(Lm,Lm,Lm,dcrxxzz) 
            call correct(Lm,Lm,Lm,dcryyzz) 
            call correct(Lm,Lm,Lm,dcrxxyz) 
            call correct(Lm,Lm,Lm,dcryyxz) 
            call correct(Lm,Lm,Lm,dcrzzxy) 

         else !4th-order interlaced

            allocate(dcr(Lm,Lm,Lm)) !do delta_g
            call assign(Nran,rr,rm,Lm,dcr,P0,nbr,ir,wr,0,0,0,0)
            call fftwnd_f77_one(planf,dcr,dcr)      
            call fcomb(Lm,dcr,Nr)

            allocate(dcrxx(Lm,Lm,Lm),dcryy(Lm,Lm,Lm),dcrzz(Lm,Lm,Lm)) !do 5 delta2_g (1st part)
            call FiveDelta2g_1(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr, 
     $      dcrxx,dcryy,dcrzz,planf)
            deallocate(dcryy,dcrzz) !done 5 delta2_g (1st part)

            allocate(dcrxy(Lm,Lm,Lm),dcryz(Lm,Lm,Lm),dcrzx(Lm,Lm,Lm)) !do 5 delta2_g (2nd part)
            call FiveDelta2g_2(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr, 
     $      dcr,dcrxx,dcrxy,dcryz,dcrzx,planf)
            deallocate(dcrxy,dcryz,dcrzx) !done 5 delta2_g (2nd part)

            allocate(dcrxxxx(Lm,Lm,Lm),dcryyyy(Lm,Lm,Lm)) !do 9 delta4_g (1st part)
            allocate(dcrzzzz(Lm,Lm,Lm),dcrxxxy(Lm,Lm,Lm))
            allocate(dcrxxxz(Lm,Lm,Lm))
            call NineDelta4g_1(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr,
     $           dcrxxxx,dcryyyy,dcrzzzz,dcrxxxy,dcrxxxz,planf)
            deallocate(dcryyyy,dcrzzzz,dcrxxxy,dcrxxxz) !done 9 delta4_g (1st part)

            allocate(dcryyyx(Lm,Lm,Lm))!do 9 delta4_g (2nd part)
            allocate(dcryyyz(Lm,Lm,Lm),dcrzzzx(Lm,Lm,Lm))
            allocate(dcrzzzy(Lm,Lm,Lm),dcrxxyy(Lm,Lm,Lm))
            call NineDelta4g_2(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr,
     $           dcrxxxx,dcryyyx,dcryyyz,dcrzzzx,dcrzzzy,dcrxxyy,planf)
            deallocate(dcryyyx,dcryyyz,dcrzzzx,dcrzzzy,dcrxxyy)!done 9 delta4_g (2nd part)

            allocate(dcrxxzz(Lm,Lm,Lm)) !do 9 delta4_g (3rd part)
            allocate(dcryyzz(Lm,Lm,Lm),dcrxxyz(Lm,Lm,Lm))
            allocate(dcryyxz(Lm,Lm,Lm),dcrzzxy(Lm,Lm,Lm))
            call NineDelta4g_3(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr,dcr,
     $      dcrxx,dcrxxxx,dcrxxzz,dcryyzz,dcrxxyz,dcryyxz,dcrzzxy,planf)
            deallocate(dcrxxzz,dcryyzz,dcrxxyz,dcryyxz,dcrzzxy)!done 9 delta4_g (3rd part)

            do i=1,Nran !build reweighted fields now
               wr(i)=wr(i)**2
            enddo   
            allocate(dcrw(Lm,Lm,Lm)) !do deltaw_g
            call assign(Nran,rr,rm,Lm,dcrw,P0,nbr,ir,wr,0,0,0,0)
            call fftwnd_f77_one(planf,dcrw,dcrw)      
            call fcomb(Lm,dcrw,Nr)

            allocate(dcrwxx(Lm,Lm,Lm)) !do 5 delta2w_g (1st part)
            allocate(dcrwyy(Lm,Lm,Lm),dcrwzz(Lm,Lm,Lm))
            call FiveDelta2g_1(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr, 
     $      dcrwxx,dcrwyy,dcrwzz,planf)
            deallocate(dcrwyy,dcrwzz) !done 5 delta2w_g (1st part)

            allocate(dcrwxy(Lm,Lm,Lm)) !do 5 delta2w_g (2nd part)
            allocate(dcrwyz(Lm,Lm,Lm),dcrwzx(Lm,Lm,Lm))
            call FiveDelta2g_2(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr, 
     $      dcrw,dcrwxx,dcrwxy,dcrwyz,dcrwzx,planf)
            deallocate(dcrwxy,dcrwyz,dcrwzx) !done 5 delta2w_g (2nd part)

            allocate(dcrwxxxx(Lm,Lm,Lm),dcrwyyyy(Lm,Lm,Lm)) !do 9 delta4_g (1st part)
            allocate(dcrwzzzz(Lm,Lm,Lm),dcrwxxxy(Lm,Lm,Lm))
            allocate(dcrwxxxz(Lm,Lm,Lm))
            call NineDelta4g_1(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr,
     $           dcrwxxxx,dcrwyyyy,dcrwzzzz,dcrwxxxy,dcrwxxxz,planf)
            deallocate(dcrwyyyy,dcrwzzzz,dcrwxxxy,dcrwxxxz) !done 9 delta4_g (1st part)

            allocate(dcrwyyyx(Lm,Lm,Lm))!do 9 delta4_g (2nd part)
            allocate(dcrwyyyz(Lm,Lm,Lm),dcrwzzzx(Lm,Lm,Lm))
            allocate(dcrwzzzy(Lm,Lm,Lm),dcrwxxyy(Lm,Lm,Lm))
            call NineDelta4g_2(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr,
     $           dcrwxxxx,dcrwyyyx,dcrwyyyz,dcrwzzzx,dcrwzzzy,dcrwxxyy,
     $           planf)
            deallocate(dcrwyyyx,dcrwyyyz,dcrwzzzx,dcrwzzzy,dcrwxxyy)!done 9 delta4_g (2nd part)

            allocate(dcrwxxzz(Lm,Lm,Lm)) !do 9 delta4_g (3rd part)
            allocate(dcrwyyzz(Lm,Lm,Lm),dcrwxxyz(Lm,Lm,Lm))
            allocate(dcrwyyxz(Lm,Lm,Lm),dcrwzzxy(Lm,Lm,Lm))
            call NineDelta4g_3(Nran,Nr,Nmax,rr,rm,Lm,P0,nbr,ir,wr,dcrw,
     $      dcrwxx,dcrwxxxx,dcrwxxzz,dcrwyyzz,dcrwxxyz,dcrwyyxz,
     $      dcrwzzxy,planf)
            deallocate(dcrwxxzz,dcrwyyzz,dcrwxxyz,dcrwyyxz,dcrwzzxy)!done 9 delta4_g (3rd part)

         endif


         call getarg(9,filecoef) !Fourier file
         open(unit=4,file=filecoef,status='unknown',form='unformatted')
         write(4)Lm
         write(4)(((dcr(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)real(I10),real(I12),real(I22),real(I13),real(I23),
     &   real(I33)
         write(4)P0,Nr,real(Nrsys),real(Nrsyscomp),real(Nrsystot)
         write(4)xmin,xmax,ymin,ymax,zmin,zmax
         write(4)(((dcrxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcrw(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcrwxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcrxxxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcrwxxxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         close(4)
         
      endif
      

 1025 format(2x,6e14.6)
 123  stop
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cosmology(idata,icosmo,Om0,OL0)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit none
      real Om0,OL0
      integer idata,icosmo
      ! fiducial cosmology (OmM=0.31 h=0.676 Ol=0.69 Obh2=0.022)
      if (icosmo.eq.0) then 
         Om0=0.31
      else
         if (idata.eq.1) then !CMASS sample
            Om0=0.274
         elseif (idata.eq.2) then !LasDamas
            Om0=0.25
         elseif (idata.eq.3) then !QPM north
!            Om0=0.29 !dr12c mocks
            Om0=0.31 !dr12d mocks
         elseif (idata.eq.4) then !QPM south
!            Om0=0.29 !dr12c mocks
            Om0=0.31 !dr12d mocks
         elseif (idata.eq.5) then !LOWZ sample
            Om0=0.3
         elseif (idata.eq.6) then !PATCHY
            Om0=0.307115
         elseif (idata.eq.7) then !PTHALOS CMASS north
            Om0=0.274
         elseif (idata.eq.8) then !PTHALOS CMASS south
            Om0=0.274
         else
            write(*,*)'specify which dataset you want!'
            stop
         endif   
      endif 
      OL0=1.-Om0 !assume flatness
      return
      end
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine nbtable(idata,Nsel,z,selfun,sec)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !build nbar table from Anderson files for use in QPM/PTH mocks
      implicit none
      integer i,idata,Nsel,Nsel2
      real selfun(Nsel),z(Nsel),sec(Nsel),dum
      character dummy*200,selfunfile*200,dir*200
      
      if (idata.eq.7) then 
         dir='/mount/riachuelo1/hahn/data/Nseries/'
         selfunfile=dir(1:len_trim(dir))//
     $    'nbar-nseries-fibcoll.dat'
      endif
      if (idata.eq.7) then
         open(unit=4,file=selfunfile,status='old',form='formatted')
         do i=1,Nsel
            read(4,*)z(i),dum,dum,selfun(i)
         enddo 
         close(4)
         call spline(z,selfun,Nsel,3e30,3e30,sec)
      endif 

      return
      end
            
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine radialdist(idata,Nsel,z,Nbin,zbin,dbin,sec3)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit none
      integer idata,Nsel,ic,Nbin
      real zend,z(Nsel),zt,zbin(Nbin),dbin(Nbin),sec3(Nbin)
      real chi
      external chi
      
      zend=1.1 !up to what z we build radial distances
      if (idata.lt.7) then
         if (zend.le.z(Nsel)) then
            write(*,*)'increase zend'
            stop
         endif
      endif
      do ic=1,Nbin !build comoving distance vs redshift relation
         zt=zend*float(ic-1)/float(Nbin-1)
         zbin(ic)=zt
         dbin(ic)=chi(zt)
      enddo
      call spline(dbin,zbin,Nbin,3e30,3e30,sec3)

      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine readmock(idata,lssfile,Nmax,nbg,wg,rg,P0,
     $           Ngal,Ngsys,Ngsyscomp,Ngsystot,compavg,
     $           xmin,xmax,ymin,ymax,zmin,zmax) !for data mocks
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit none
      integer i,idata,Ngal,Nariel,Nmax
      real xmin,xmax,ymin,ymax,zmin,zmax,xcm,ycm,zcm,nbar,chi,pi
      real ra,dec,az,comp,nbb,wsys,wred,cspeed,dum,az2,wnoz,bias,wfc
      real wboss,veto,rad,P0,nbg(Nmax),wg(Nmax),rg(3,Nmax)
      parameter(cspeed=299800.0,pi=3.141592654)
      real*8 Ngsys,Ngsyscomp,Ngsystot,compavg
      character*200 lssfile,lssinfofile,dummy  
      external nbar,chi

      open(unit=4,file=lssfile,status='old',form='formatted')
      Ngal=0 !true Ngal (=Ng) will get determined later after survey is put into a box 
      Ngsys=0.d0
      Ngsyscomp=0.d0
      Ngsystot=0.d0
      compavg=0.d0
      xmin=1000.
      xmax=-1000.
      ymin=1000.
      ymax=-1000.
      zmin=1000.
      zmax=-1000.
      xcm=0.
      ycm=0.
      zcm=0.
      if (idata.eq.1) then 
         read(4,*)Nariel !ariel has Nobjects in first line:
      elseif (idata.eq.7) then 
         read(4,'(a)')dummy !skip comment line
      endif
      do i=1,Nmax
         if (idata.eq.7) then !Nseries
            read(4,*,end=13)ra,dec,az,wfc,comp
            wsys=1.0
            wred=wfc
            nbb=nbar(az)
            nbg(i)=nbb*comp
         endif
            
         compavg=compavg+dble(comp)
         wg(i)=wsys*wred/(1.+nbg(i)*P0/comp)/comp !all weights
         Ngal=Ngal+1
         Ngsys=Ngsys+dble(wsys*wred)
         Ngsyscomp=Ngsyscomp+dble(wsys*wred/comp)
         Ngsystot=Ngsystot+dble(wg(i))
         ra=ra*(pi/180.)
         dec=dec*(pi/180.)
         rad=chi(az)
         rg(1,i)=rad*cos(dec)*cos(ra)
         rg(2,i)=rad*cos(dec)*sin(ra)
         rg(3,i)=rad*sin(dec)
         xmin=min(xmin,rg(1,i))
         xmax=max(xmax,rg(1,i))
         ymin=min(ymin,rg(2,i))
         ymax=max(ymax,rg(2,i))
         zmin=min(zmin,rg(3,i))
         zmax=max(zmax,rg(3,i))
         xcm=xcm+rg(1,i)
         ycm=ycm+rg(2,i)
         zcm=zcm+rg(3,i)
      enddo
 13   continue
      close(4)
      write(*,*)xmin,'<= xD <=',xmax
      write(*,*)ymin,'<= yD <=',ymax
      write(*,*)zmin,'<= zD <=',zmax
      xcm=xcm/float(Ngal)
      ycm=ycm/float(Ngal)
      zcm=zcm/float(Ngal)
c         x0=0.5*(xmin+xmax)
c         y0=0.5*(ymin+ymax)
c         z0=0.5*(zmin+zmax)
c         write(*,*)'CoM is at',xcm,ycm,zcm
c         write(*,*)'origin at',x0,y0,z0
c         write(*,*)'enter origin of coordinates (from randoms)'
c         read(*,*)x0,y0,z0
c         do i=1,Ngal !move origin
c            rg(1,i)=rg(1,i)-x0
c            rg(2,i)=rg(2,i)-y0
c            rg(3,i)=rg(3,i)-z0
c         enddo

      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine readrand(idata,randomfile,Nmax,nbr,wr,rr,
     $           P0,Nran,Nrsys,Nrsyscomp,Nrsystot,compavg,izbin,
     $           xmin,xmax,ymin,ymax,zmin,zmax,thetaobs,phiobs) !for random mocks
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit none
      integer i,idata,ifc,icomp,Nran,Nariel,Nmax,iveto,izbin
      real xmin,xmax,ymin,ymax,zmin,zmax,xcm,ycm,zcm,nbar,chi,pi
      real ra,dec,az,comp,nbb,wsys,wred,cspeed,dum,az2,wnoz,wcp,bias
      real wboss,veto,rad,P0,nbr(Nmax),wr(Nmax),rr(3,Nmax)
      real rcm,thetaobs,phiobs
      parameter(cspeed=299800.0,pi=3.141592654)
      real*8 Nrsys,Nrsyscomp,Nrsystot,compavg
      character*200 randomfile,dummy  
      external nbar,chi

         open(unit=4,file=randomfile,status='old',form='formatted')
         Nran=0 !Ngal will get determined later after survey is put into a box (Nr)
         Nrsys=0.d0
         Nrsyscomp=0.d0
         Nrsystot=0.d0
         compavg=0.d0
         xmin=1000.
         xmax=-1000.
         ymin=1000.
         ymax=-1000.
         zmin=1000.
         zmax=-1000.
         xcm=0.
         ycm=0.
         zcm=0.
         if (idata.eq.1) then !BOSS
            read(4,*)Nariel !ariel has Nobjects in first line
         elseif (idata.eq.7) then
            read(4,'(a)')dum !skip comment line
         endif 
         do i=1,Nmax
            if (idata.eq.7) then ! Nseries
               read(4,*,end=15)ra,dec,az,comp
               wsys=1.
               wred=1.
               nbb=nbar(az)
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            endif
            
            Nran=Nran+1
            if (Nran.eq.Nmax) then
               write(*,*)'increase Nmax',Nran
               stop
            endif   

            compavg=compavg+dble(comp)
            wr(i) =wsys*wred/(1.+nbr(i)*P0/comp)/comp !all weights
            Nrsys=Nrsys+dble(wsys*wred)
            Nrsyscomp=Nrsyscomp+dble(wsys*wred/comp)
            Nrsystot=Nrsystot+dble(wr(i))

            ra=ra*(pi/180.)
            dec=dec*(pi/180.)
            rad=chi(az)
            rr(1,i)=rad*cos(dec)*cos(ra)
            rr(2,i)=rad*cos(dec)*sin(ra)
            rr(3,i)=rad*sin(dec)
            xmin=min(xmin,rr(1,i))
            xmax=max(xmax,rr(1,i))
            ymin=min(ymin,rr(2,i))
            ymax=max(ymax,rr(2,i))
            zmin=min(zmin,rr(3,i))
            zmax=max(zmax,rr(3,i))
            xcm=xcm+rr(1,i)
            ycm=ycm+rr(2,i)
            zcm=zcm+rr(3,i)
         enddo
 15      continue
         close(4)
         write(*,*)xmin,'<= xR <=',xmax
         write(*,*)ymin,'<= yR <=',ymax
         write(*,*)zmin,'<= zR <=',zmax
         xcm=xcm/float(Nran)
         ycm=ycm/float(Nran)
         zcm=zcm/float(Nran)
c         x0=0.5*(xmin+xmax)
c         y0=0.5*(ymin+ymax)
c         z0=0.5*(zmin+zmax)
         write(*,*)'CoM is at',xcm,ycm,zcm !save for LOS direction
c         write(*,*)'origin at',x0,y0,z0
c         do i=1,Nran !move origin
c            rr(1,i)=rr(1,i)-x0
c            rr(2,i)=rr(2,i)-y0
c            rr(3,i)=rr(3,i)-z0
c         enddo
C         rcm=sqrt(xcm**2+ycm**2+zcm**2)
C         thetaobs=acos(zcm/rcm)
C         phiobs=acos(xcm/sqrt(xcm**2+ycm**2))
C         write(*,*)'observer to CoM angles:',thetaobs,phiobs

      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine computeIij(iflag,Nmax,nbg,wg,ig,Ngal,
     $                      I10,I12,I22,I13,I23,I33)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit none
      integer Nmax,ig(Nmax),i,iflag,Ngal
      real*8 I10,I12,I22,I13,I23,I33
      real nbg(Nmax),wg(Nmax),nb,weight
         I10=0.d0
         I12=0.d0
         I22=0.d0
         I13=0.d0
         I23=0.d0
         I33=0.d0
         do i=1,Ngal 
            nb=nbg(i)  
            weight=wg(i)
            if (ig(i).eq.1) then 
               I10=I10+1.d0  
               I12=I12+dble(weight**2)  
               I22=I22+dble(nb*weight**2)
               I13=I13+dble(weight**3) 
               I23=I23+dble(nb*weight**3 ) 
               I33=I33+dble(nb**2 *weight**3) 
            endif
         enddo
         if (iflag.eq.0) then !data
            write(*,*)'these are normalization integrals from data'
         else !random
            write(*,*)'these are normalization integrals from random'
            write(*,*)'to be scaled by alpha in power spectrum code'
         endif
         write(*,*)'I10=',I10
         write(*,*)'I12=',I12
         write(*,*)'I22=',I22
         write(*,*)'I13=',I13
         write(*,*)'I23=',I23
         write(*,*)'I33=',I33
         
         return
         end

       subroutine FiveDelta2g_1(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg, 
     $      dcgxx,dcgyy,dcgzz,planf)
       implicit none
       integer Ngal,Lm,Nmax,ig(Nmax),Ng,ix,iy,iz,ikx,iky,ikz
       integer*8 planf
       real rg(3,Nmax),rm,P0,nbg(Nmax),wg(Nmax),kxh,kyh,kzh,rk
       complex dcgxx(Lm,Lm,Lm),dcgyy(Lm,Lm,Lm)
       complex dcgzz(Lm,Lm,Lm)
       
            call assign(Ngal,rg,rm,Lm,dcgxx,P0,nbg,ig,wg,1,1,0,0)
            call assign(Ngal,rg,rm,Lm,dcgyy,P0,nbg,ig,wg,2,2,0,0)
            call assign(Ngal,rg,rm,Lm,dcgzz,P0,nbg,ig,wg,3,3,0,0)
            call fftwnd_f77_one(planf,dcgxx,dcgxx)      
            call fftwnd_f77_one(planf,dcgyy,dcgyy)      
            call fftwnd_f77_one(planf,dcgzz,dcgzz)      
            call fcomb(Lm,dcgxx,Ng)
            call fcomb(Lm,dcgyy,Ng)
            call fcomb(Lm,dcgzz,Ng)
            do 104 iz=1,Lm !build quadrupole
               ikz=mod(iz+Lm/2-2,Lm)-Lm/2+1
               do 104 iy=1,Lm
                  iky=mod(iy+Lm/2-2,Lm)-Lm/2+1
                  do 104 ix=1,Lm/2+1
                     ikx=mod(ix+Lm/2-2,Lm)-Lm/2+1
                     rk=sqrt(float(ikx**2+iky**2+ikz**2))
                     if(rk.gt.0.)then
                        kxh=float(ikx)/rk !unit vectors
                        kyh=float(iky)/rk
                        kzh=float(ikz)/rk
                        dcgxx(ix,iy,iz)=7.5*(dcgxx(ix,iy,iz)*kxh**2 
     &                  +dcgyy(ix,iy,iz)*kyh**2
     &                  +dcgzz(ix,iy,iz)*kzh**2)
                     end if
 104        continue

       return
       end

       subroutine FiveDelta2g_2(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg, 
     $      dcg,dcgxx,dcgxy,dcgyz,dcgzx,planf)
       implicit none
       integer Ngal,Lm,Nmax,ig(Nmax),Ng,ix,iy,iz,ikx,iky,ikz
       integer*8 planf
       real rg(3,Nmax),rm,P0,nbg(Nmax),wg(Nmax),kxh,kyh,kzh,rk
       complex dcg(Lm,Lm,Lm),dcgxx(Lm,Lm,Lm),dcgxy(Lm,Lm,Lm)
       complex dcgyz(Lm,Lm,Lm),dcgzx(Lm,Lm,Lm)
       
            call assign(Ngal,rg,rm,Lm,dcgxy,P0,nbg,ig,wg,1,2,0,0)
            call assign(Ngal,rg,rm,Lm,dcgyz,P0,nbg,ig,wg,2,3,0,0)
            call assign(Ngal,rg,rm,Lm,dcgzx,P0,nbg,ig,wg,3,1,0,0)
            call fftwnd_f77_one(planf,dcgxy,dcgxy)      
            call fftwnd_f77_one(planf,dcgyz,dcgyz)      
            call fftwnd_f77_one(planf,dcgzx,dcgzx)      
            call fcomb(Lm,dcgxy,Ng)
            call fcomb(Lm,dcgyz,Ng)
            call fcomb(Lm,dcgzx,Ng)
            do 105 iz=1,Lm !build quadrupole
               ikz=mod(iz+Lm/2-2,Lm)-Lm/2+1
               do 105 iy=1,Lm
                  iky=mod(iy+Lm/2-2,Lm)-Lm/2+1
                  do 105 ix=1,Lm/2+1
                     ikx=mod(ix+Lm/2-2,Lm)-Lm/2+1
                     rk=sqrt(float(ikx**2+iky**2+ikz**2))
                     if(rk.gt.0.)then
                        kxh=float(ikx)/rk !unit vectors
                        kyh=float(iky)/rk
                        kzh=float(ikz)/rk
                        dcgxx(ix,iy,iz)=dcgxx(ix,iy,iz) + 7.5*( 
     &                   2.*dcgxy(ix,iy,iz)*kxh*kyh
     &                  +2.*dcgyz(ix,iy,iz)*kyh*kzh
     &                  +2.*dcgzx(ix,iy,iz)*kzh*kxh)
     &                  -2.5*dcg(ix,iy,iz)   
                     end if
 105        continue

       return
       end

       subroutine FiveDelta2g(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,
     $           dcg,dcgxx,dcgyy,dcgzz,dcgxy,dcgyz,dcgzx,planf)
       implicit none
       integer Ngal,Lm,Nmax,ig(Nmax),Ng,ix,iy,iz,ikx,iky,ikz
       integer*8 planf
       real rg(3,Nmax),rm,P0,nbg(Nmax),wg(Nmax),kxh,kyh,kzh,rk
       complex dcg(Lm,Lm,Lm)
       complex dcgxx(Lm,Lm,Lm),dcgyy(Lm,Lm,Lm),dcgzz(Lm,Lm,Lm)
       complex dcgxy(Lm,Lm,Lm),dcgyz(Lm,Lm,Lm),dcgzx(Lm,Lm,Lm)
       
       call assign(Ngal,rg,rm,Lm,dcgxx,P0,nbg,ig,wg,1,1,0,0)
       call assign(Ngal,rg,rm,Lm,dcgyy,P0,nbg,ig,wg,2,2,0,0)
       call assign(Ngal,rg,rm,Lm,dcgzz,P0,nbg,ig,wg,3,3,0,0)
       call assign(Ngal,rg,rm,Lm,dcgxy,P0,nbg,ig,wg,1,2,0,0)
       call assign(Ngal,rg,rm,Lm,dcgyz,P0,nbg,ig,wg,2,3,0,0)
       call assign(Ngal,rg,rm,Lm,dcgzx,P0,nbg,ig,wg,3,1,0,0)
       call fftwnd_f77_one(planf,dcgxx,dcgxx)      
       call fftwnd_f77_one(planf,dcgyy,dcgyy)      
       call fftwnd_f77_one(planf,dcgzz,dcgzz)      
       call fftwnd_f77_one(planf,dcgxy,dcgxy)      
       call fftwnd_f77_one(planf,dcgyz,dcgyz)      
       call fftwnd_f77_one(planf,dcgzx,dcgzx)      
       call fcomb(Lm,dcgxx,Ng)
       call fcomb(Lm,dcgyy,Ng)
       call fcomb(Lm,dcgzz,Ng)
       call fcomb(Lm,dcgxy,Ng)
       call fcomb(Lm,dcgyz,Ng)
       call fcomb(Lm,dcgzx,Ng)
       do 100 iz=1,Lm !build quadrupole
          ikz=mod(iz+Lm/2-2,Lm)-Lm/2+1
             do 100 iy=1,Lm
                iky=mod(iy+Lm/2-2,Lm)-Lm/2+1
                do 100 ix=1,Lm/2+1
                   ikx=mod(ix+Lm/2-2,Lm)-Lm/2+1
                   rk=sqrt(float(ikx**2+iky**2+ikz**2))
                   if(rk.gt.0.)then
                      kxh=float(ikx)/rk !unit vectors
                      kyh=float(iky)/rk
                      kzh=float(ikz)/rk
                      dcgxx(ix,iy,iz)=7.5*(dcgxx(ix,iy,iz)*kxh**2 
     &                +dcgyy(ix,iy,iz)*kyh**2+dcgzz(ix,iy,iz)*kzh**2 
     &          +2.*dcgxy(ix,iy,iz)*kxh*kyh+2.*dcgyz(ix,iy,iz)*kyh*kzh
     &               +2.*dcgzx(ix,iy,iz)*kzh*kxh)-2.5*dcg(ix,iy,iz) 
                     end if
 100   continue

       return
       end
       
       subroutine NineDelta4g_1(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,
     $           dcgxxxx,dcgyyyy,dcgzzzz,dcgxxxy,dcgxxxz,planf)
       implicit none
       integer Ngal,Lm,Nmax,ig(Nmax),Ng,ix,iy,iz,ikx,iky,ikz
       integer*8 planf
       real rg(3,Nmax),rm,P0,nbg(Nmax),wg(Nmax),kxh,kyh,kzh,rk
       complex dcgxxxx(Lm,Lm,Lm),dcgyyyy(Lm,Lm,Lm),dcgzzzz(Lm,Lm,Lm)
       complex dcgxxxy(Lm,Lm,Lm),dcgxxxz(Lm,Lm,Lm) 
       
            call assign(Ngal,rg,rm,Lm,dcgxxxx,P0,nbg,ig,wg,1,1,1,1) 
            call assign(Ngal,rg,rm,Lm,dcgyyyy,P0,nbg,ig,wg,2,2,2,2) 
            call assign(Ngal,rg,rm,Lm,dcgzzzz,P0,nbg,ig,wg,3,3,3,3) 
            call assign(Ngal,rg,rm,Lm,dcgxxxy,P0,nbg,ig,wg,1,1,1,2) 
            call assign(Ngal,rg,rm,Lm,dcgxxxz,P0,nbg,ig,wg,1,1,1,3) 
            call fftwnd_f77_one(planf,dcgxxxx,dcgxxxx)
            call fftwnd_f77_one(planf,dcgyyyy,dcgyyyy)
            call fftwnd_f77_one(planf,dcgzzzz,dcgzzzz)
            call fftwnd_f77_one(planf,dcgxxxy,dcgxxxy)
            call fftwnd_f77_one(planf,dcgxxxz,dcgxxxz)
            call fcomb(Lm,dcgxxxx,Ng)
            call fcomb(Lm,dcgyyyy,Ng)
            call fcomb(Lm,dcgzzzz,Ng)
            call fcomb(Lm,dcgxxxy,Ng)
            call fcomb(Lm,dcgxxxz,Ng)
            do 101 iz=1,Lm !build hexadecapole (part 2)
               ikz=mod(iz+Lm/2-2,Lm)-Lm/2+1
               do 101 iy=1,Lm
                  iky=mod(iy+Lm/2-2,Lm)-Lm/2+1
                  do 101 ix=1,Lm/2+1
                     ikx=mod(ix+Lm/2-2,Lm)-Lm/2+1
                     rk=sqrt(float(ikx**2+iky**2+ikz**2))
                     if(rk.gt.0.)then
                        kxh=float(ikx)/rk !unit vectors
                        kyh=float(iky)/rk
                        kzh=float(ikz)/rk
                 dcgxxxx(ix,iy,iz)=35.*9./8.*(dcgxxxx(ix,iy,iz)*kxh**4 
     &              +dcgyyyy(ix,iy,iz)*kyh**4+dcgzzzz(ix,iy,iz)*kzh**4
     & +4.*dcgxxxy(ix,iy,iz)*kxh**3*kyh+4.*dcgxxxz(ix,iy,iz)*kxh**3*kzh)
                     end if
 101        continue

       return
       end

       subroutine NineDelta4g_2(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,
     $           dcgxxxx,dcgyyyx,dcgyyyz,dcgzzzx,dcgzzzy,dcgxxyy,planf)
       implicit none
       integer Ngal,Lm,Nmax,ig(Nmax),Ng,ix,iy,iz,ikx,iky,ikz
       integer*8 planf
       real rg(3,Nmax),rm,P0,nbg(Nmax),wg(Nmax),kxh,kyh,kzh,rk
       complex dcgyyyx(Lm,Lm,Lm),dcgyyyz(Lm,Lm,Lm),dcgzzzx(Lm,Lm,Lm)
       complex dcgzzzy(Lm,Lm,Lm),dcgxxyy(Lm,Lm,Lm),dcgxxxx(Lm,Lm,Lm) 
       
            call assign(Ngal,rg,rm,Lm,dcgyyyx,P0,nbg,ig,wg,2,2,2,1) 
            call assign(Ngal,rg,rm,Lm,dcgyyyz,P0,nbg,ig,wg,2,2,2,3) 
            call assign(Ngal,rg,rm,Lm,dcgzzzx,P0,nbg,ig,wg,3,3,3,1) 
            call assign(Ngal,rg,rm,Lm,dcgzzzy,P0,nbg,ig,wg,3,3,3,2) 
            call assign(Ngal,rg,rm,Lm,dcgxxyy,P0,nbg,ig,wg,1,1,2,2) 
            call fftwnd_f77_one(planf,dcgyyyx,dcgyyyx)
            call fftwnd_f77_one(planf,dcgyyyz,dcgyyyz)
            call fftwnd_f77_one(planf,dcgzzzx,dcgzzzx)
            call fftwnd_f77_one(planf,dcgzzzy,dcgzzzy)
            call fftwnd_f77_one(planf,dcgxxyy,dcgxxyy)
            call fcomb(Lm,dcgyyyx,Ng)
            call fcomb(Lm,dcgyyyz,Ng)
            call fcomb(Lm,dcgzzzx,Ng)
            call fcomb(Lm,dcgzzzy,Ng)
            call fcomb(Lm,dcgxxyy,Ng)
            do 102 iz=1,Lm !build hexadecapole (part 2)
               ikz=mod(iz+Lm/2-2,Lm)-Lm/2+1
               do 102 iy=1,Lm
                  iky=mod(iy+Lm/2-2,Lm)-Lm/2+1
                  do 102 ix=1,Lm/2+1
                     ikx=mod(ix+Lm/2-2,Lm)-Lm/2+1
                     rk=sqrt(float(ikx**2+iky**2+ikz**2))
                     if(rk.gt.0.)then
                        kxh=float(ikx)/rk !unit vectors
                        kyh=float(iky)/rk
                        kzh=float(ikz)/rk
                 dcgxxxx(ix,iy,iz)=dcgxxxx(ix,iy,iz) + 35.*9./8.*(
     &  4.*dcgyyyx(ix,iy,iz)*kyh**3*kxh+4.*dcgyyyz(ix,iy,iz)*kyh**3*kzh 
     & +4.*dcgzzzx(ix,iy,iz)*kzh**3*kxh+4.*dcgzzzy(ix,iy,iz)*kzh**3*kyh 
     & +6.*dcgxxyy(ix,iy,iz)*kxh**2*kyh**2)
                     end if
 102        continue

       return
       end

       subroutine NineDelta4g_3(Ngal,Ng,Nmax,rg,rm,Lm,P0,nbg,ig,wg,dcg,
     $      dcgxx,dcgxxxx,dcgxxzz,dcgyyzz,dcgxxyz,dcgyyxz,dcgzzxy,planf)
       implicit none
       integer Ngal,Lm,Nmax,ig(Nmax),Ng,ix,iy,iz,ikx,iky,ikz
       integer*8 planf
       real rg(3,Nmax),rm,P0,nbg(Nmax),wg(Nmax),kxh,kyh,kzh,rk
       complex dcg(Lm,Lm,Lm),dcgxx(Lm,Lm,Lm)
       complex dcgxxxx(Lm,Lm,Lm),dcgxxzz(Lm,Lm,Lm),dcgyyzz(Lm,Lm,Lm)
       complex dcgxxyz(Lm,Lm,Lm),dcgyyxz(Lm,Lm,Lm),dcgzzxy(Lm,Lm,Lm) 
       
            call assign(Ngal,rg,rm,Lm,dcgxxzz,P0,nbg,ig,wg,1,1,3,3) 
            call assign(Ngal,rg,rm,Lm,dcgyyzz,P0,nbg,ig,wg,2,2,3,3) 
            call assign(Ngal,rg,rm,Lm,dcgxxyz,P0,nbg,ig,wg,1,1,2,3) 
            call assign(Ngal,rg,rm,Lm,dcgyyxz,P0,nbg,ig,wg,2,2,1,3) 
            call assign(Ngal,rg,rm,Lm,dcgzzxy,P0,nbg,ig,wg,3,3,1,2) 
            call fftwnd_f77_one(planf,dcgxxzz,dcgxxzz)
            call fftwnd_f77_one(planf,dcgyyzz,dcgyyzz)
            call fftwnd_f77_one(planf,dcgxxyz,dcgxxyz)
            call fftwnd_f77_one(planf,dcgyyxz,dcgyyxz)
            call fftwnd_f77_one(planf,dcgzzxy,dcgzzxy)
            call fcomb(Lm,dcgxxzz,Ng)
            call fcomb(Lm,dcgyyzz,Ng)
            call fcomb(Lm,dcgxxyz,Ng)
            call fcomb(Lm,dcgyyxz,Ng)
            call fcomb(Lm,dcgzzxy,Ng)
            do 103 iz=1,Lm !build hexadecapole (part 2)
               ikz=mod(iz+Lm/2-2,Lm)-Lm/2+1
               do 103 iy=1,Lm
                  iky=mod(iy+Lm/2-2,Lm)-Lm/2+1
                  do 103 ix=1,Lm/2+1
                     ikx=mod(ix+Lm/2-2,Lm)-Lm/2+1
                     rk=sqrt(float(ikx**2+iky**2+ikz**2))
                     if(rk.gt.0.)then
                        kxh=float(ikx)/rk !unit vectors
                        kyh=float(iky)/rk
                        kzh=float(ikz)/rk
                 dcgxxxx(ix,iy,iz)=dcgxxxx(ix,iy,iz) + 35.*9./8.*(
     &                          6.*dcgxxzz(ix,iy,iz)*kxh**2*kzh**2
     &                         +6.*dcgyyzz(ix,iy,iz)*kyh**2*kzh**2
     &                         +12.*dcgxxyz(ix,iy,iz)*kxh**2*kyh*kzh
     &                         +12.*dcgyyxz(ix,iy,iz)*kyh**2*kxh*kzh
     &                         +12.*dcgzzxy(ix,iy,iz)*kzh**2*kxh*kyh)
     &                   -9./2.*dcgxx(ix,iy,iz) -63./8.*dcg(ix,iy,iz)
                     end if
 103        continue

       return
       end



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function nbar(QQ) !nbar(z)   for QPM mocks                  
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      implicit none
      integer Nsel
      parameter(Nsel=181)      
      real z(Nsel),selfun(Nsel),sec(Nsel),qq,az
      common /interpol/z,selfun,sec
      az=QQ
      call splint(z,selfun,sec,Nsel,az,nbar)
      RETURN
      END
cc*******************************************************************
      subroutine assign_CIC(Npar,r,rm,Lm,dtl,wg,ig,ia,ib,ic,id)
cc*******************************************************************
      implicit none
      integer Lm,ixp,iyp,izp,ixa,iya,iza,i,ix,iy,iz,ia,ib,Npar,ic,id
      real r(3,Npar),rm(2),wg(Npar),we
      real dtl(Lm+2,Lm,Lm)
      real dx,dy,dz,rx,ry,rz,rnorm,rvec(3)
      integer ig(Npar)
c
      do 11 iz=1,Lm
       do 11 iy=1,Lm
        do 11 ix=1,Lm+2
11        dtl(ix,iy,iz)=0.

      do i=1,Npar
      if (ig(i).eq.1) then
         rx=rm(1)*r(1,i)+rm(2)
         ry=rm(1)*r(2,i)+rm(2)
         rz=rm(1)*r(3,i)+rm(2)
         ixp=int(rx)
         iyp=int(ry)
         izp=int(rz)
         dx=rx-real(ixp)
         dy=ry-real(iyp)
         dz=rz-real(izp)
         ixa=mod(ixp,Lm)+1
         iya=mod(iyp,Lm)+1
         iza=mod(izp,Lm)+1

C            we=wg(i)
         if (ia.eq.0 .and. ib.eq.0 .and. ic.eq.0 .and. id.eq.0) then !FFT delta
            we=wg(i)
         elseif (ic.eq.0 .and. id.eq.0) then !FFT Qij
            rnorm=r(1,i)**2+r(2,i)**2+r(3,i)**2
            we=wg(i)*r(ia,i)*r(ib,i)/rnorm
         else !FFT Qijkl
            rnorm=r(1,i)**2+r(2,i)**2+r(3,i)**2
            we=wg(i)*r(ia,i)*r(ib,i)*r(ic,i)*r(id,i)/rnorm**2
         endif

         dtl(ixa,iya,iza) = dtl(ixa,iya,iza)+dx*dy*dz *we
         dtl(ixa,iya,izp) = dtl(ixa,iya,izp)+dx*dy*(1.-dz) *we

         dtl(ixp,iya,iza) = dtl(ixp,iya,iza)+(1.-dx)*dy*dz *we
         dtl(ixp,iya,izp) = dtl(ixp,iya,izp)+(1.-dx)*dy*(1.-dz) *we

         dtl(ixa,iyp,iza) = dtl(ixa,iyp,iza)+dx*(1.-dy)*dz *we
         dtl(ixa,iyp,izp) = dtl(ixa,iyp,izp)+dx*(1.-dy)*(1.-dz) *we

         dtl(ixp,iyp,iza) = dtl(ixp,iyp,iza)+(1.-dx)*(1.-dy)*dz *we
         dtl(ixp,iyp,izp) = dtl(ixp,iyp,izp)+(1.-dx)*(1.-dy)*(1.-dz)
     $    *we
      endif
      enddo

C      if (ia.gt.0 .and. ib.gt.0 .and. ic.eq.0 .and. id.eq.0) then  !FFT Qij
C        do iz=1,Lm
C         rvec(3)=-float(Lm)/2.+float(iz-1)
C         do iy=1,Lm
C           rvec(2)=-float(Lm)/2.+float(iy-1)
C           do ix=1,Lm
C             rvec(1)=-float(Lm)/2.+float(ix-1)
C             rnorm=rvec(1)**2+rvec(2)**2+rvec(3)**2
C             if (rnorm.gt.0.) then 
C               dtl(ix,iy,iz)=dtl(ix,iy,iz)*rvec(ia)*rvec(ib)/rnorm
C             endif  
C           enddo
C         enddo
C        enddo 
C       endif
C
C      if (ia.gt.0 .and. ib.gt.0 .and. ic.gt.0 .and. id.gt.0) then  !FFT Qijkl
C        do iz=1,Lm
C         rvec(3)=-float(Lm)/2.+float(iz-1)
C         do iy=1,Lm
C           rvec(2)=-float(Lm)/2.+float(iy-1)
C           do ix=1,Lm
C             rvec(1)=-float(Lm)/2.+float(ix-1)
C             rnorm=rvec(1)**2+rvec(2)**2+rvec(3)**2
C             if (rnorm.gt.0.) then 
C               dtl(ix,iy,iz)=dtl(ix,iy,iz)*rvec(ia)*rvec(ib)*rvec(ic)
C     $                       *rvec(id)/rnorm**2
C             endif  
C           enddo
C         enddo
C        enddo 
C       endif

      return
      end
cc*******************************************************************
      subroutine correct(Lx,Ly,Lz,dtl)
cc*******************************************************************
      complex dtl(Lx/2+1,Ly,Lz)
      real rkz,rky,rkx,tpiLx,tpiLy,tpiLz,Wkz,Wky,Wkx,cf,cfac
      integer icz,icy,icx,iflag
      iflag=2

      tpi=6.283185307
      cf=1.
      tpiLx=tpi/float(Lx)
      tpiLy=tpi/float(Ly)
      tpiLz=tpi/float(Lz)
      do 300 iz=1,Lz/2+1
         icz=mod(Lz-iz+1,Lz)+1
         rkz=tpiLz*float(iz-1)
         Wkz=1.
         if(rkz.ne.0.)Wkz=(sin(rkz/2.)/(rkz/2.))**iflag
         do 300 iy=1,Ly/2+1
            icy=mod(Ly-iy+1,Ly)+1
            rky=tpiLy*float(iy-1)
            Wky=1.
            if(rky.ne.0.)Wky=(sin(rky/2.)/(rky/2.))**iflag
            do 300 ix=1,Lx/2+1
               rkx=tpiLx*float(ix-1)
               Wkx=1.
               if(rkx.ne.0.)Wkx=(sin(rkx/2.)/(rkx/2.))**iflag
               cfac=cf/(Wkx*Wky*Wkz)
               dtl(ix,iy,iz)=dtl(ix,iy,iz)*cfac
               if(iz.ne.icz) dtl(ix,iy,icz)=dtl(ix,iy,icz)*cfac
               if(iy.ne.icy) dtl(ix,icy,iz)=dtl(ix,icy,iz)*cfac
               if(iz.ne.icz .and. iy.ne.icy) then
                  dtl(ix,icy,icz)=dtl(ix,icy,icz)*cfac
               endif
 300              continue
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Subroutine PutIntoBox(Ng,rg,Rbox,ig,Ng2,Nmax)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit none
      integer Nmax
      integer Ng,Ng2,j,i,ig(Nmax)
      real rg(3,Nmax),Rbox,acheck
      j=0
      do i=1,Ng
         acheck=abs(rg(1,i))+abs(rg(2,i))+abs(rg(3,i))
         if (acheck.gt.0.) then 
            if (abs(rg(1,i)).lt.Rbox .and. abs(rg(2,i)).lt.Rbox .and. 
     $           abs(rg(3,i)).lt.Rbox) then !put into box
               j=j+1
               ig(i)=1
            else
               ig(i)=0
            endif
         else
            ig(i)=0
         endif
      enddo
      Ng2=j
      RETURN
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine assign(N,r,rm,L,dtl,P0,nbg,ig,w,ia,ib,ic,id) !with FKP weighing
cc*******************************************************************
      implicit none
      integer N,ig(N),ix,iy,iz,L,i,ixm1,iym1,izm1,ixm2,ixp1,ixp2
      real dtl(2*L,L,L),r(3,N),rm(2),we,ar,nb,nbg(N),w(N)
      real rca,rcb,rx,ry,rz,tx,ty,tz,hx,hx2,hxm2,hxm1,hxp2,hxp1
      integer iym2,iyp1,iyp2,izm2,izp1,izp2,nxm1,nym1,nzm1
      real hy,hy2,hym2,hym1,hyp2,hyp1,hz,hz2,hzm1,hzm2,hzp1,hzp2
      real gx,gx2,gxm2,gxm1,gxp2,gxp1,gy,gy2,gym2,gym1,gyp2,gyp1
      integer nym2,nyp1,nyp2,nzm2,nzp1,nzp2,nxm2,nxp1,nxp2
      real gz,gz2,gzm2,gzm1,gzp2,gzp1,P0,rnorm,rvec(3)
      integer ia,ib,ic,id
c
      do 1 iz=1,L
       do 1 iy=1,L
        do 1 ix=1,2*L
1        dtl(ix,iy,iz)=0.
c
      rca=rm(1)
      rcb=rm(2)
c
      do 2 i=1,N
      if (ig(i).eq.1) then

       ! here goes the FKP weight and any other per-galaxy quantity
         if (ia.eq.0 .and. ib.eq.0 .and. ic.eq.0 .and. id.eq.0) then !FFT delta
            we=w(i)
         elseif (ic.eq.0 .and. id.eq.0) then !FFT Qij
            rnorm=r(1,i)**2+r(2,i)**2+r(3,i)**2
            we=w(i)*r(ia,i)*r(ib,i)/rnorm
         else !FFT Qijkl
            rnorm=r(1,i)**2+r(2,i)**2+r(3,i)**2
            we=w(i)*r(ia,i)*r(ib,i)*r(ic,i)*r(id,i)/rnorm**2
         endif

       rx=rca*r(1,i)+rcb
       ry=rca*r(2,i)+rcb
       rz=rca*r(3,i)+rcb
       tx=rx+0.5
       ty=ry+0.5
       tz=rz+0.5
       ixm1=int(rx)
       iym1=int(ry)
       izm1=int(rz)
       ixm2=2*mod(ixm1-2+L,L)+1
       ixp1=2*mod(ixm1,L)+1
       ixp2=2*mod(ixm1+1,L)+1
       hx=rx-ixm1
       ixm1=2*ixm1-1
       hx2=hx*hx
       hxm2=(1.-hx)**3
       hxm1=4.+(3.*hx-6.)*hx2
       hxp2=hx2*hx
       hxp1=6.-hxm2-hxm1-hxp2
c
       iym2=mod(iym1-2+L,L)+1
       iyp1=mod(iym1,L)+1
       iyp2=mod(iym1+1,L)+1
       hy=ry-iym1
       hy2=hy*hy
       hym2=(1.-hy)**3
       hym1=4.+(3.*hy-6.)*hy2
       hyp2=hy2*hy
       hyp1=6.-hym2-hym1-hyp2
c
       izm2=mod(izm1-2+L,L)+1
       izp1=mod(izm1,L)+1
       izp2=mod(izm1+1,L)+1
       hz=rz-izm1
       hz2=hz*hz
       hzm2=(1.-hz)**3
       hzm1=4.+(3.*hz-6.)*hz2
       hzp2=hz2*hz
       hzp1=6.-hzm2-hzm1-hzp2
c
       nxm1=int(tx)
       nym1=int(ty)
       nzm1=int(tz)
c
       gx=tx-nxm1
       nxm1=mod(nxm1-1,L)+1
       nxm2=2*mod(nxm1-2+L,L)+2
       nxp1=2*mod(nxm1,L)+2
       nxp2=2*mod(nxm1+1,L)+2
       nxm1=2*nxm1
       gx2=gx*gx
       gxm2=(1.-gx)**3
       gxm1=4.+(3.*gx-6.)*gx2
       gxp2=gx2*gx
       gxp1=6.-gxm2-gxm1-gxp2
c
       gy=ty-nym1
       nym1=mod(nym1-1,L)+1
       nym2=mod(nym1-2+L,L)+1
       nyp1=mod(nym1,L)+1
       nyp2=mod(nym1+1,L)+1
       gy2=gy*gy
       gym2=(1.-gy)**3
       gym1=4.+(3.*gy-6.)*gy2
       gyp2=gy2*gy
       gyp1=6.-gym2-gym1-gyp2
c
       gz=tz-nzm1
       nzm1=mod(nzm1-1,L)+1
       nzm2=mod(nzm1-2+L,L)+1
       nzp1=mod(nzm1,L)+1
       nzp2=mod(nzm1+1,L)+1
       gz2=gz*gz
       gzm2=(1.-gz)**3
       gzm1=4.+(3.*gz-6.)*gz2
       gzp2=gz2*gz
       gzp1=6.-gzm2-gzm1-gzp2
c
       dtl(ixm2,iym2,izm2)   = dtl(ixm2,iym2,izm2)+ hxm2*hym2 *hzm2*we
       dtl(ixm1,iym2,izm2)   = dtl(ixm1,iym2,izm2)+ hxm1*hym2 *hzm2*we
       dtl(ixp1,iym2,izm2)   = dtl(ixp1,iym2,izm2)+ hxp1*hym2 *hzm2*we
       dtl(ixp2,iym2,izm2)   = dtl(ixp2,iym2,izm2)+ hxp2*hym2 *hzm2*we
       dtl(ixm2,iym1,izm2)   = dtl(ixm2,iym1,izm2)+ hxm2*hym1 *hzm2*we
       dtl(ixm1,iym1,izm2)   = dtl(ixm1,iym1,izm2)+ hxm1*hym1 *hzm2*we
       dtl(ixp1,iym1,izm2)   = dtl(ixp1,iym1,izm2)+ hxp1*hym1 *hzm2*we
       dtl(ixp2,iym1,izm2)   = dtl(ixp2,iym1,izm2)+ hxp2*hym1 *hzm2*we
       dtl(ixm2,iyp1,izm2)   = dtl(ixm2,iyp1,izm2)+ hxm2*hyp1 *hzm2*we
       dtl(ixm1,iyp1,izm2)   = dtl(ixm1,iyp1,izm2)+ hxm1*hyp1 *hzm2*we
       dtl(ixp1,iyp1,izm2)   = dtl(ixp1,iyp1,izm2)+ hxp1*hyp1 *hzm2*we
       dtl(ixp2,iyp1,izm2)   = dtl(ixp2,iyp1,izm2)+ hxp2*hyp1 *hzm2*we
       dtl(ixm2,iyp2,izm2)   = dtl(ixm2,iyp2,izm2)+ hxm2*hyp2 *hzm2*we
       dtl(ixm1,iyp2,izm2)   = dtl(ixm1,iyp2,izm2)+ hxm1*hyp2 *hzm2*we
       dtl(ixp1,iyp2,izm2)   = dtl(ixp1,iyp2,izm2)+ hxp1*hyp2 *hzm2*we
       dtl(ixp2,iyp2,izm2)   = dtl(ixp2,iyp2,izm2)+ hxp2*hyp2 *hzm2*we
       dtl(ixm2,iym2,izm1)   = dtl(ixm2,iym2,izm1)+ hxm2*hym2 *hzm1*we
       dtl(ixm1,iym2,izm1)   = dtl(ixm1,iym2,izm1)+ hxm1*hym2 *hzm1*we
       dtl(ixp1,iym2,izm1)   = dtl(ixp1,iym2,izm1)+ hxp1*hym2 *hzm1*we
       dtl(ixp2,iym2,izm1)   = dtl(ixp2,iym2,izm1)+ hxp2*hym2 *hzm1*we
       dtl(ixm2,iym1,izm1)   = dtl(ixm2,iym1,izm1)+ hxm2*hym1 *hzm1*we
       dtl(ixm1,iym1,izm1)   = dtl(ixm1,iym1,izm1)+ hxm1*hym1 *hzm1*we
       dtl(ixp1,iym1,izm1)   = dtl(ixp1,iym1,izm1)+ hxp1*hym1 *hzm1*we
       dtl(ixp2,iym1,izm1)   = dtl(ixp2,iym1,izm1)+ hxp2*hym1 *hzm1*we
       dtl(ixm2,iyp1,izm1)   = dtl(ixm2,iyp1,izm1)+ hxm2*hyp1 *hzm1*we
       dtl(ixm1,iyp1,izm1)   = dtl(ixm1,iyp1,izm1)+ hxm1*hyp1 *hzm1*we
       dtl(ixp1,iyp1,izm1)   = dtl(ixp1,iyp1,izm1)+ hxp1*hyp1 *hzm1*we
       dtl(ixp2,iyp1,izm1)   = dtl(ixp2,iyp1,izm1)+ hxp2*hyp1 *hzm1*we
       dtl(ixm2,iyp2,izm1)   = dtl(ixm2,iyp2,izm1)+ hxm2*hyp2 *hzm1*we
       dtl(ixm1,iyp2,izm1)   = dtl(ixm1,iyp2,izm1)+ hxm1*hyp2 *hzm1*we
       dtl(ixp1,iyp2,izm1)   = dtl(ixp1,iyp2,izm1)+ hxp1*hyp2 *hzm1*we
       dtl(ixp2,iyp2,izm1)   = dtl(ixp2,iyp2,izm1)+ hxp2*hyp2 *hzm1*we
       dtl(ixm2,iym2,izp1)   = dtl(ixm2,iym2,izp1)+ hxm2*hym2 *hzp1*we
       dtl(ixm1,iym2,izp1)   = dtl(ixm1,iym2,izp1)+ hxm1*hym2 *hzp1*we
       dtl(ixp1,iym2,izp1)   = dtl(ixp1,iym2,izp1)+ hxp1*hym2 *hzp1*we
       dtl(ixp2,iym2,izp1)   = dtl(ixp2,iym2,izp1)+ hxp2*hym2 *hzp1*we
       dtl(ixm2,iym1,izp1)   = dtl(ixm2,iym1,izp1)+ hxm2*hym1 *hzp1*we
       dtl(ixm1,iym1,izp1)   = dtl(ixm1,iym1,izp1)+ hxm1*hym1 *hzp1*we
       dtl(ixp1,iym1,izp1)   = dtl(ixp1,iym1,izp1)+ hxp1*hym1 *hzp1*we
       dtl(ixp2,iym1,izp1)   = dtl(ixp2,iym1,izp1)+ hxp2*hym1 *hzp1*we
       dtl(ixm2,iyp1,izp1)   = dtl(ixm2,iyp1,izp1)+ hxm2*hyp1 *hzp1*we
       dtl(ixm1,iyp1,izp1)   = dtl(ixm1,iyp1,izp1)+ hxm1*hyp1 *hzp1*we
       dtl(ixp1,iyp1,izp1)   = dtl(ixp1,iyp1,izp1)+ hxp1*hyp1 *hzp1*we
       dtl(ixp2,iyp1,izp1)   = dtl(ixp2,iyp1,izp1)+ hxp2*hyp1 *hzp1*we
       dtl(ixm2,iyp2,izp1)   = dtl(ixm2,iyp2,izp1)+ hxm2*hyp2 *hzp1*we
       dtl(ixm1,iyp2,izp1)   = dtl(ixm1,iyp2,izp1)+ hxm1*hyp2 *hzp1*we
       dtl(ixp1,iyp2,izp1)   = dtl(ixp1,iyp2,izp1)+ hxp1*hyp2 *hzp1*we
       dtl(ixp2,iyp2,izp1)   = dtl(ixp2,iyp2,izp1)+ hxp2*hyp2 *hzp1*we
       dtl(ixm2,iym2,izp2)   = dtl(ixm2,iym2,izp2)+ hxm2*hym2 *hzp2*we
       dtl(ixm1,iym2,izp2)   = dtl(ixm1,iym2,izp2)+ hxm1*hym2 *hzp2*we
       dtl(ixp1,iym2,izp2)   = dtl(ixp1,iym2,izp2)+ hxp1*hym2 *hzp2*we
       dtl(ixp2,iym2,izp2)   = dtl(ixp2,iym2,izp2)+ hxp2*hym2 *hzp2*we
       dtl(ixm2,iym1,izp2)   = dtl(ixm2,iym1,izp2)+ hxm2*hym1 *hzp2*we
       dtl(ixm1,iym1,izp2)   = dtl(ixm1,iym1,izp2)+ hxm1*hym1 *hzp2*we
       dtl(ixp1,iym1,izp2)   = dtl(ixp1,iym1,izp2)+ hxp1*hym1 *hzp2*we
       dtl(ixp2,iym1,izp2)   = dtl(ixp2,iym1,izp2)+ hxp2*hym1 *hzp2*we
       dtl(ixm2,iyp1,izp2)   = dtl(ixm2,iyp1,izp2)+ hxm2*hyp1 *hzp2*we
       dtl(ixm1,iyp1,izp2)   = dtl(ixm1,iyp1,izp2)+ hxm1*hyp1 *hzp2*we
       dtl(ixp1,iyp1,izp2)   = dtl(ixp1,iyp1,izp2)+ hxp1*hyp1 *hzp2*we
       dtl(ixp2,iyp1,izp2)   = dtl(ixp2,iyp1,izp2)+ hxp2*hyp1 *hzp2*we
       dtl(ixm2,iyp2,izp2)   = dtl(ixm2,iyp2,izp2)+ hxm2*hyp2 *hzp2*we
       dtl(ixm1,iyp2,izp2)   = dtl(ixm1,iyp2,izp2)+ hxm1*hyp2 *hzp2*we
       dtl(ixp1,iyp2,izp2)   = dtl(ixp1,iyp2,izp2)+ hxp1*hyp2 *hzp2*we
       dtl(ixp2,iyp2,izp2)   = dtl(ixp2,iyp2,izp2)+ hxp2*hyp2 *hzp2*we
c
       dtl(nxm2,nym2,nzm2)   = dtl(nxm2,nym2,nzm2)+ gxm2*gym2 *gzm2*we
       dtl(nxm1,nym2,nzm2)   = dtl(nxm1,nym2,nzm2)+ gxm1*gym2 *gzm2*we
       dtl(nxp1,nym2,nzm2)   = dtl(nxp1,nym2,nzm2)+ gxp1*gym2 *gzm2*we
       dtl(nxp2,nym2,nzm2)   = dtl(nxp2,nym2,nzm2)+ gxp2*gym2 *gzm2*we
       dtl(nxm2,nym1,nzm2)   = dtl(nxm2,nym1,nzm2)+ gxm2*gym1 *gzm2*we
       dtl(nxm1,nym1,nzm2)   = dtl(nxm1,nym1,nzm2)+ gxm1*gym1 *gzm2*we
       dtl(nxp1,nym1,nzm2)   = dtl(nxp1,nym1,nzm2)+ gxp1*gym1 *gzm2*we
       dtl(nxp2,nym1,nzm2)   = dtl(nxp2,nym1,nzm2)+ gxp2*gym1 *gzm2*we
       dtl(nxm2,nyp1,nzm2)   = dtl(nxm2,nyp1,nzm2)+ gxm2*gyp1 *gzm2*we
       dtl(nxm1,nyp1,nzm2)   = dtl(nxm1,nyp1,nzm2)+ gxm1*gyp1 *gzm2*we
       dtl(nxp1,nyp1,nzm2)   = dtl(nxp1,nyp1,nzm2)+ gxp1*gyp1 *gzm2*we
       dtl(nxp2,nyp1,nzm2)   = dtl(nxp2,nyp1,nzm2)+ gxp2*gyp1 *gzm2*we
       dtl(nxm2,nyp2,nzm2)   = dtl(nxm2,nyp2,nzm2)+ gxm2*gyp2 *gzm2*we
       dtl(nxm1,nyp2,nzm2)   = dtl(nxm1,nyp2,nzm2)+ gxm1*gyp2 *gzm2*we
       dtl(nxp1,nyp2,nzm2)   = dtl(nxp1,nyp2,nzm2)+ gxp1*gyp2 *gzm2*we
       dtl(nxp2,nyp2,nzm2)   = dtl(nxp2,nyp2,nzm2)+ gxp2*gyp2 *gzm2*we
       dtl(nxm2,nym2,nzm1)   = dtl(nxm2,nym2,nzm1)+ gxm2*gym2 *gzm1*we
       dtl(nxm1,nym2,nzm1)   = dtl(nxm1,nym2,nzm1)+ gxm1*gym2 *gzm1*we
       dtl(nxp1,nym2,nzm1)   = dtl(nxp1,nym2,nzm1)+ gxp1*gym2 *gzm1*we
       dtl(nxp2,nym2,nzm1)   = dtl(nxp2,nym2,nzm1)+ gxp2*gym2 *gzm1*we
       dtl(nxm2,nym1,nzm1)   = dtl(nxm2,nym1,nzm1)+ gxm2*gym1 *gzm1*we
       dtl(nxm1,nym1,nzm1)   = dtl(nxm1,nym1,nzm1)+ gxm1*gym1 *gzm1*we
       dtl(nxp1,nym1,nzm1)   = dtl(nxp1,nym1,nzm1)+ gxp1*gym1 *gzm1*we
       dtl(nxp2,nym1,nzm1)   = dtl(nxp2,nym1,nzm1)+ gxp2*gym1 *gzm1*we
       dtl(nxm2,nyp1,nzm1)   = dtl(nxm2,nyp1,nzm1)+ gxm2*gyp1 *gzm1*we
       dtl(nxm1,nyp1,nzm1)   = dtl(nxm1,nyp1,nzm1)+ gxm1*gyp1 *gzm1*we
       dtl(nxp1,nyp1,nzm1)   = dtl(nxp1,nyp1,nzm1)+ gxp1*gyp1 *gzm1*we
       dtl(nxp2,nyp1,nzm1)   = dtl(nxp2,nyp1,nzm1)+ gxp2*gyp1 *gzm1*we
       dtl(nxm2,nyp2,nzm1)   = dtl(nxm2,nyp2,nzm1)+ gxm2*gyp2 *gzm1*we
       dtl(nxm1,nyp2,nzm1)   = dtl(nxm1,nyp2,nzm1)+ gxm1*gyp2 *gzm1*we
       dtl(nxp1,nyp2,nzm1)   = dtl(nxp1,nyp2,nzm1)+ gxp1*gyp2 *gzm1*we
       dtl(nxp2,nyp2,nzm1)   = dtl(nxp2,nyp2,nzm1)+ gxp2*gyp2 *gzm1*we
       dtl(nxm2,nym2,nzp1)   = dtl(nxm2,nym2,nzp1)+ gxm2*gym2 *gzp1*we
       dtl(nxm1,nym2,nzp1)   = dtl(nxm1,nym2,nzp1)+ gxm1*gym2 *gzp1*we
       dtl(nxp1,nym2,nzp1)   = dtl(nxp1,nym2,nzp1)+ gxp1*gym2 *gzp1*we
       dtl(nxp2,nym2,nzp1)   = dtl(nxp2,nym2,nzp1)+ gxp2*gym2 *gzp1*we
       dtl(nxm2,nym1,nzp1)   = dtl(nxm2,nym1,nzp1)+ gxm2*gym1 *gzp1*we
       dtl(nxm1,nym1,nzp1)   = dtl(nxm1,nym1,nzp1)+ gxm1*gym1 *gzp1*we
       dtl(nxp1,nym1,nzp1)   = dtl(nxp1,nym1,nzp1)+ gxp1*gym1 *gzp1*we
       dtl(nxp2,nym1,nzp1)   = dtl(nxp2,nym1,nzp1)+ gxp2*gym1 *gzp1*we
       dtl(nxm2,nyp1,nzp1)   = dtl(nxm2,nyp1,nzp1)+ gxm2*gyp1 *gzp1*we
       dtl(nxm1,nyp1,nzp1)   = dtl(nxm1,nyp1,nzp1)+ gxm1*gyp1 *gzp1*we
       dtl(nxp1,nyp1,nzp1)   = dtl(nxp1,nyp1,nzp1)+ gxp1*gyp1 *gzp1*we
       dtl(nxp2,nyp1,nzp1)   = dtl(nxp2,nyp1,nzp1)+ gxp2*gyp1 *gzp1*we
       dtl(nxm2,nyp2,nzp1)   = dtl(nxm2,nyp2,nzp1)+ gxm2*gyp2 *gzp1*we
       dtl(nxm1,nyp2,nzp1)   = dtl(nxm1,nyp2,nzp1)+ gxm1*gyp2 *gzp1*we
       dtl(nxp1,nyp2,nzp1)   = dtl(nxp1,nyp2,nzp1)+ gxp1*gyp2 *gzp1*we
       dtl(nxp2,nyp2,nzp1)   = dtl(nxp2,nyp2,nzp1)+ gxp2*gyp2 *gzp1*we
       dtl(nxm2,nym2,nzp2)   = dtl(nxm2,nym2,nzp2)+ gxm2*gym2 *gzp2*we
       dtl(nxm1,nym2,nzp2)   = dtl(nxm1,nym2,nzp2)+ gxm1*gym2 *gzp2*we
       dtl(nxp1,nym2,nzp2)   = dtl(nxp1,nym2,nzp2)+ gxp1*gym2 *gzp2*we
       dtl(nxp2,nym2,nzp2)   = dtl(nxp2,nym2,nzp2)+ gxp2*gym2 *gzp2*we
       dtl(nxm2,nym1,nzp2)   = dtl(nxm2,nym1,nzp2)+ gxm2*gym1 *gzp2*we
       dtl(nxm1,nym1,nzp2)   = dtl(nxm1,nym1,nzp2)+ gxm1*gym1 *gzp2*we
       dtl(nxp1,nym1,nzp2)   = dtl(nxp1,nym1,nzp2)+ gxp1*gym1 *gzp2*we
       dtl(nxp2,nym1,nzp2)   = dtl(nxp2,nym1,nzp2)+ gxp2*gym1 *gzp2*we
       dtl(nxm2,nyp1,nzp2)   = dtl(nxm2,nyp1,nzp2)+ gxm2*gyp1 *gzp2*we
       dtl(nxm1,nyp1,nzp2)   = dtl(nxm1,nyp1,nzp2)+ gxm1*gyp1 *gzp2*we
       dtl(nxp1,nyp1,nzp2)   = dtl(nxp1,nyp1,nzp2)+ gxp1*gyp1 *gzp2*we
       dtl(nxp2,nyp1,nzp2)   = dtl(nxp2,nyp1,nzp2)+ gxp2*gyp1 *gzp2*we
       dtl(nxm2,nyp2,nzp2)   = dtl(nxm2,nyp2,nzp2)+ gxm2*gyp2 *gzp2*we
       dtl(nxm1,nyp2,nzp2)   = dtl(nxm1,nyp2,nzp2)+ gxm1*gyp2 *gzp2*we
       dtl(nxp1,nyp2,nzp2)   = dtl(nxp1,nyp2,nzp2)+ gxp1*gyp2 *gzp2*we
       dtl(nxp2,nyp2,nzp2)   = dtl(nxp2,nyp2,nzp2)+ gxp2*gyp2 *gzp2*we
       endif
2     continue
c
c      write(*,*)rca,rcb,P0

C       if (ia.gt.0 .and. ib.gt.0) then  !this is a bit risky given the details of interlacing
C        do iz=1,L
C         rvec(3)=-float(L)/2.+float(iz-1)
C         do iy=1,L
C           rvec(2)=-float(L)/2.+float(iy-1)
C           do ix=1,L
C             rvec(1)=-float(L)/2.+float(ix-1)
C             rnorm=rvec(1)**2+rvec(2)**2+rvec(3)**2
C             if (rnorm.gt.0.) then 
C               dtl(ix,iy,iz)=dtl(ix,iy,iz)*rvec(ia)*rvec(ib)/rnorm
C             endif  
C           enddo
C         enddo
C        enddo 
C        do iz=1,L
C         rvec(3)=-float(L)/2.+float(iz-1)+0.5
C         do iy=1,L
C           rvec(2)=-float(L)/2.+float(iy-1)+0.5
C           do ix=L+1,2*L
C             rvec(1)=-float(L)/2.+float(ix-1-L)+0.5
C             rnorm=rvec(1)**2+rvec(2)**2+rvec(3)**2
C             if (rnorm.gt.0.) then 
C               dtl(ix,iy,iz)=dtl(ix,iy,iz)*rvec(ia)*rvec(ib)/rnorm
C             endif  
C           enddo
C         enddo
C        enddo 
C       endif
C
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function zdis(ar) !interpolation redshift(distance)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      parameter(Nbin=151)
      common /interp3/dbin,zbin,sec3
      real dbin(Nbin),zbin(Nbin),sec3(Nbin)
      call splint(dbin,zbin,sec3,Nbin,ar,zdis)
      RETURN
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function chi(x) !radial distance in Mpc/h as a function of z
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real*8 qmin,qmax,ans,rdi,epsabs,epsrel,abserr
      parameter (limit=1000)
      integer neval,ier,iord(limit),last
      real*8 alist(limit),blist(limit),elist(limit),rlist(limit)
      external rdi,dqage
      common/radint/Om0,OL0
      qmin=0.d0
      qmax=dble(x)
      epsabs=0.d0
      epsrel=1.d-2                                                            
      call dqage(rdi,qmin,qmax,epsabs,epsrel,30,limit,ans,abserr,
     $ neval,ier,alist,blist,rlist,elist,iord,last)
      chi=real(ans)
      RETURN
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real*8 function rdi(z) !radial distance integrand
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      common/radint/Om0,OL0
      real Om0,OL0
      real*8 z
      rdi=3000.d0/dsqrt(OL0+(1.d0-Om0-OL0)*(1.d0+z)**2+Om0*(1.d0+z)**3)
      return
      end
cc*******************************************************************
      subroutine fcomb(L,dcl,N)
cc*******************************************************************
      implicit none
      integer L,Lnyq,ix,iy,iz,icx,icy,icz,N
      real cf,rkx,rky,rkz,wkx,wky,wkz,cfac
      complex dcl(L,L,L)
      real*8 tpi,tpiL,piL
      complex*16 rec,xrec,yrec,zrec
      complex c1,ci,c000,c001,c010,c011,cma,cmb,cmc,cmd
      tpi=6.283185307d0
c
      cf=1./(6.**3*4.) !*float(N)) not needed for FKP
      Lnyq=L/2+1
      tpiL=tpi/float(L)
      piL=-tpiL/2.
      rec=cmplx(dcos(piL),dsin(piL))
      c1=cmplx(1.,0.)
      ci=cmplx(0.,1.)
      zrec=c1
      do 301 iz=1,Lnyq
       icz=mod(L-iz+1,L)+1
       rkz=tpiL*(iz-1)
       Wkz=1.
       if(rkz.ne.0.)Wkz=(sin(rkz/2.)/(rkz/2.))**4
       yrec=c1
       do 302 iy=1,Lnyq
        icy=mod(L-iy+1,L)+1
        rky=tpiL*(iy-1)
        Wky=1.
        if(rky.ne.0.)Wky=(sin(rky/2.)/(rky/2.))**4
        xrec=c1
        do 303 ix=1,Lnyq
         icx=mod(L-ix+1,L)+1
         rkx=tpiL*(ix-1)
         Wkx=1.
         if(rkx.ne.0.)Wkx=(sin(rkx/2.)/(rkx/2.))**4
         cfac=cf/(Wkx*Wky*Wkz)
c
         cma=ci*xrec*yrec*zrec
         cmb=ci*xrec*yrec*conjg(zrec)
         cmc=ci*xrec*conjg(yrec)*zrec
         cmd=ci*xrec*conjg(yrec*zrec)
c
         c000=dcl(ix,iy ,iz )*(c1-cma)+conjg(dcl(icx,icy,icz))*(c1+cma)
         c001=dcl(ix,iy ,icz)*(c1-cmb)+conjg(dcl(icx,icy,iz ))*(c1+cmb)
         c010=dcl(ix,icy,iz )*(c1-cmc)+conjg(dcl(icx,iy ,icz))*(c1+cmc)
         c011=dcl(ix,icy,icz)*(c1-cmd)+conjg(dcl(icx,iy ,iz ))*(c1+cmd)
c
c
         dcl(ix,iy ,iz )=c000*cfac
         dcl(ix,iy ,icz)=c001*cfac
         dcl(ix,icy,iz )=c010*cfac
         dcl(ix,icy,icz)=c011*cfac
         dcl(icx,iy ,iz )=conjg(dcl(ix,icy,icz))
         dcl(icx,iy ,icz)=conjg(dcl(ix,icy,iz ))
         dcl(icx,icy,iz )=conjg(dcl(ix,iy ,icz))
         dcl(icx,icy,icz)=conjg(dcl(ix,iy ,iz ))
c
         xrec=xrec*rec
303     continue
        yrec=yrec*rec
302    continue
       zrec=zrec*rec
301   continue
c
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      include '/home/users/hahn/powercode/dqage.f'
      include '/home/users/hahn/powercode/d1mach.f'
      include '/home/users/hahn/powercode/dqawfe.f'
      include '/home/users/hahn/powercode/spline.f'
      include '/home/users/hahn/powercode/splint.f'
      include '/home/users/hahn/powercode/gabqx.f'
