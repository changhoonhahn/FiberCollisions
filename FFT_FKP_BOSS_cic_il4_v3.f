! ifort  -fast -o FFT_FKP_BOSS_cic_il4_v3.exe FFT_FKP_BOSS_cic_il4_v3.f -L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw -lm
! as v2 but with different comp weights and shot noise subtraction
! also optional comp weighting and for QPM optional fiber collisions
      implicit none  !as FFT_FKP_BOSS_cic_il4 but outputs quadrupole field directly (instead of 6 files)
      integer Nsel,Nran,i,iwr,Ngal,Nmax,n,kx,ky,kz,Lm,Ngrid,ix,iy,iz,j,k
      integer Ng,Nr,iflag,ic,Nbin, l,interpol,Nariel,idata, Nsel2, iveto
      integer*8 planf,plan_real
      real pi,cspeed,wfc,wrf,xmin,xmax,ymin,ymax,zmin,zmax,nbar2
      parameter(Nsel=181,Nmax=3*10**8,Nbin=201,Nsel2=10)
      parameter(pi=3.141592654)
      integer grid, ifc, icomp
      dimension grid(3)
      parameter(cspeed=299800.0)
      integer, allocatable :: ig(:),ir(:)      
      real zbin(Nbin),dbin(Nbin),sec3(Nbin),zt,dum,gfrac,az2
      real cz,Om0,OL0,chi,nbar,Rbox,x0,y0,z0,wnoz,wcp
      real, allocatable :: nbg(:),nbr(:),rg(:,:),rr(:,:)
      real, allocatable :: wg(:),wr(:),ws(:)
      real*8, allocatable :: compg(:),compr(:)
      real selfun(Nsel),z(Nsel),sec(Nsel),az,ra,dec,rad,numden,zend
      real w, nbb, xcm, ycm, zcm, rcm, bias, wboss, veto
      real thetaobs,phiobs, selfun2(Nsel2),z2(Nsel2),sec2(Nsel2)
      real alpha,P0,nb,weight,ar,akf,Fr,Fi,Gr,Gi,wsys,wred,comp,wfkp
      real*8 I10,I12,I22,I13,I23,I33,I12F
      real*8 compavg,Nrsys,Nrsyscomp,Ngsys,Ngsyscomp,Ngsystot,Nrsystot
      real rsq,xr,yr,zr,xg,yg,zg
      real*8 I12xx,I12yy,I12zz,I12xy,I12yz,I12zx
      real kdotr,vol,xscale,rlow,rm(2)
      integer ikz,icz,iky,icy,ikx,icx
      real rk,kxh,kyh,kzh
      complex, allocatable :: dcg(:,:,:),dcr(:,:,:)
      complex, allocatable :: dcgw(:,:,:),dcrw(:,:,:)
      complex, allocatable :: dcgxx(:,:,:),dcrxx(:,:,:)
      complex, allocatable :: dcgyy(:,:,:),dcryy(:,:,:)
      complex, allocatable :: dcgzz(:,:,:),dcrzz(:,:,:)
      complex, allocatable :: dcgxy(:,:,:),dcrxy(:,:,:)
      complex, allocatable :: dcgyz(:,:,:),dcryz(:,:,:)
      complex, allocatable :: dcgzx(:,:,:),dcrzx(:,:,:)
      character selfunfile*200,lssfile*200,randomfile*200,filecoef*200
      character dummy*200,fname*200,outname*200,dir*200,Omstr*200
      character Rboxstr*200,Ngridstr*200,interpolstr*200,iflagstr*200
      character P0str*200,typestr*200,lssinfofile*200
      character*200 icompstr,ifcstr
      common /interpol/z,selfun,sec
      common /interpol2/z2,selfun2,sec2
      common /interp3/dbin,zbin,sec3
      common /radint/Om0,OL0
      external nbar,chi,PutIntoBox,assign2,fcomb,nbar2
      include '/usr/local/src/fftw-2.1.5/fortran/fftw_f77.i'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,*)'Omega_m at z=0: Om0'
c      read(*,*)Om0
c      OL0=1.-Om0
c      call getarg(1,Omstr)
c      read(Omstr,*)Om0
      ! idata= 1 (BOSS data), 2(LasDamas), 3-4 (QPM mocks)
      call getarg(1,typestr)
      read(typestr,*)idata 
      
      ! fiducial cosmology (OmM=0.31 h=0.676 Ol=0.69 Obh2=0.022)
      if (idata.eq.1) then !CMASS sample
         Om0=0.274
      elseif (idata.eq.2) then !LasDamas
         Om0=0.25
      elseif (idata.eq.3) then !QPM north
!         Om0=0.29 !dr12c mocks
         Om0=0.31 !dr12d mocks
      elseif (idata.eq.4) then !QPM south
!         Om0=0.29 !dr12c mocks
         Om0=0.31 !dr12d mocks
      elseif (idata.eq.5) then !LOWZ sample
         Om0=0.3
      elseif (idata.eq.6) then !PATCHY
         Om0=0.307115
      elseif (idata.eq.7) then !PTHALOS CMASS north
         Om0=0.274
      elseif (idata.eq.8) then !PTHALOS CMASS south
         Om0=0.274
      elseif (idata.eq.9) then ! TILING MOCK 
         Om0=0.274
      elseif (idata.eq.10) then !Nseries
         Om0=0.31 ! fiducial cosmology
      elseif (idata.eq.11) then ! Downsampled LasDamas
         Om0=0.25 
      else
         write(*,*)'specify which dataset you want!'
         stop
      endif   
      OL0=1.-Om0 !assume flatness

      !build nbar table from Anderson files for use in QPM/PTH mocks
      if (idata.eq.4) then !QPMsouth
         dir='/mount/riachuelo2/rs123/BOSS/QPM/cmass/'
         selfunfile=dir(1:len_trim(dir))//'cmass-dr12v1-sgc.zsel' !dr12c
         selfunfile=dir(1:len_trim(dir))//
     $    'nbar-cmass-dr12v4-S-Reid-om0p31.dat' !dr12d 
      elseif (idata.eq.3) then !QPMnorth
         dir='/mount/riachuelo2/rs123/BOSS/QPM/cmass/'
         selfunfile=dir(1:len_trim(dir))//'cmass-dr12v1-ngc.zsel' !dr12c
         selfunfile=dir(1:len_trim(dir))//
     $    'nbar-cmass-dr12v4-N-Reid-om0p31.dat' !dr12d
      elseif (idata.eq.7) then !PTHnorth
         dir='/mount/riachuelo2/rs123/BOSS/PTHalos/'
         selfunfile=dir(1:len_trim(dir))//
     $    'nzfit_dr11_vm22_north.txt'
      elseif (idata.eq.8) then !PTHsouth
         dir='/mount/riachuelo2/rs123/BOSS/PTHalos/'
         selfunfile=dir(1:len_trim(dir))//
     $    'nzfit_dr11_vm22_south.txt'
      elseif (idata.eq.10) then ! Nseries
         selfunfile='/mount/riachuelo1/hahn/data/Nseries'//
     $    '/nbar-nseries-fibcoll.dat'
      elseif (idata.eq.11) then 
         selfunfile='/mount/riachuelo1/hahn/data/LasDamas/Geo'//
     $    '/nbar-lasdamasgeo.down_nz.dat'
      else !just to have z(Nsel)            
         dir='/mount/riachuelo2/rs123/BOSS/QPM/cmass/'
         ! directy has been changed in order to remove header
         selfunfile=dir(1:len_trim(dir))//'cmass-dr12v1-ngc.zsel'
      endif
      if (idata.eq.4 .or. idata.eq.3) then 
         open(unit=4,file=selfunfile,status='old',form='formatted')
         do i=1,2 !skip *3* comment lines
            read(4,'(a)')dummy
         enddo
         do i=1,Nsel
            read(4,*)z(i),dum,dum,selfun(i),dum,dum,dum
         enddo   
         close(4)
         call spline(z,selfun,Nsel,3e30,3e30,sec)
      elseif (idata.lt.7 .or. idata.eq.9) then 
         open(unit=4,file=selfunfile,status='old',form='formatted')
         do i=1,3 !skip 2 comment lines
            read(4,'(a)')dummy
         enddo
         do i=1,Nsel
            read(4,*)z(i),selfun(i)
         enddo   
         close(4)
         call spline(z,selfun,Nsel,3e30,3e30,sec)
      elseif (idata.eq.10 .or. idata.eq.11) then 
         ! Nseries fiducial cosmology or Downsampled LasDamas Geo
         open(unit=4,file=selfunfile,status='old',form='formatted')
         do i=1,Nsel
            read(4,*)z(i),dum,dum,selfun(i)
         enddo   
         close(4)
         call spline(z,selfun,Nsel,3e30,3e30,sec)
      else
         open(unit=4,file=selfunfile,status='old',form='formatted')
         do i=1,Nsel2
            read(4,*)z2(i),selfun2(i)
         enddo   
         close(4)
         call spline(z2,selfun2,Nsel2,3e30,3e30,sec2)
      endif      
            
      zend=1.1 !up to what z we build radial distances
      if (zend.le.z(Nsel)) then
         write(*,*)'increase zend'
         stop
      endif
      do ic=1,Nbin !build comoving distance vs redshift relation
         zt=zend*float(ic-1)/float(Nbin-1)
         zbin(ic)=zt
         dbin(ic)=chi(zt)
      enddo
      call spline(dbin,zbin,Nbin,3e30,3e30,sec3)

c      write(*,*)'box side (survey origin will be at box center)'
c      read(*,*)Rbox 
      call getarg(2,Rboxstr)
      read(Rboxstr,*)Rbox
      xscale=RBox !actual side of Fourier Box
      Rbox=0.5*RBox !"radius of the box"
      rlow=-Rbox !lowest coordinate in any direction

c      write(*,*)'FFT Grid Size (Ngrid)'
c      read(*,*)Ngrid
      call getarg(3,Ngridstr)
      read(Ngridstr,*)Ngrid
      grid(1) = Ngrid
      grid(2) = Ngrid
      grid(3) = Ngrid

c      write(*,*)'2nd order cic (2) or 4th-order interlaced(4)?'
c      read(*,*)interpol
      call getarg(4,interpolstr)
      read(interpolstr,*)interpol
      if (interpol.eq.2) then !CIC needs RtoC transform
         call rfftw3d_f77_create_plan(plan_real,Ngrid,Ngrid,Ngrid,
     &   FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE + FFTW_IN_PLACE)
      else
         call fftwnd_f77_create_plan(planf,3,grid,FFTW_BACKWARD,
     $   FFTW_ESTIMATE + FFTW_IN_PLACE)
      endif
      Lm=Ngrid
      rm(1)=float(Lm)/xscale
      rm(2)=1.-rlow*rm(1)

      
c      write(*,*)'mock (0) or random mock(1)?'
c      read(*,*)iflag
      call getarg(5,iflagstr) 
      read(iflagstr,*)iflag
c      write(*,*)'FKP weight P0?'
c      read(*,*)P0
      call getarg(6,P0str)
      read(P0str,*)P0

      call getarg(7,ifcstr) 
      read(ifcstr,*)ifc
      call getarg(8,icompstr) 
      read(icompstr,*)icomp

      if (iflag.eq.0) then ! run on mock
      
c         write(*,*)'Mock Survey File'
c         read(*,'(a)')lssfile
         call getarg(9,lssfile)
c         fname='/mount/chichipio2/hahn/data/'//lssfile 
         allocate(rg(3,Nmax),nbg(Nmax),ig(Nmax),wg(Nmax),compg(Nmax))
         allocate(ws(Nmax))
         write(*,*)lssfile
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
            read(4,*)Nariel !ariel has Nobjects in first line
         elseif (idata.eq.3 .or. idata.eq.4) then !QPM info mock files
            lssinfofile=lssfile(1:len_trim(lssfile))//'.info'
            !open(unit=5,file=lssinfofile,status='old',form='formatted')
            !do i=1,3 !skip 3 comment lines
            !   read(5,'(a)')dummy
            !enddo
         elseif(idata.eq.7 .or. idata.eq.8) then
            read(4,'(a)')dummy !skip comment line
         endif
         do i=1,Nmax
            if (idata.eq.1) then !BOSS
               read(4,*,end=13)ra,dec,az,comp,nbb,wsys,wred
               nbg(i)=nbb*comp
            elseif (idata.eq.2) then !LasDamas
               read(4,*,end=13)ra,dec,az,wred
               az=az!/cspeed
               nbb=0.0000944233
               wsys=1.
               comp=1.
               nbg(i)=nbb*comp
            elseif (idata.eq.11) then ! Down sampled LasDamas
               read(4,*,end=13)ra,dec,az,wred
               az=az
               nbb=nbar(az)
               wsys=1.
               comp=1.
               nbg(i)=nbb*comp
            elseif (idata.eq.3 .or. idata.eq.4 .or. idata.eq.10) then ! QPM and Nseries
 33            read(4,*,end=13)ra,dec,az,wred,comp
               !read(5,*,end=13)dum,comp,dum,az2,dum,dum,dum
               !if (abs(az2/az-1.).gt.1.e-5) then
               !   write(*,*)'problem matching info to std mock file'
               !   write(*,*)az,az2
               !   stop
               !endif
               ! IFC NOT NECESSARY SINCE I FEED IN ALREADY FIBERCOLLIDED
               ! DATA
               !if (ifc.eq.1) then !impose fiber collisions + veto 
               !   if (wred.eq.0.) then !collision, read another
               !      goto 33
               !   endif
               !else   
               !   wred=1. !take it, and set unit weight
               !endif
               nbb=nbar(az)
               wsys=1
               nbg(i)=nbb*comp
            elseif (idata.eq.5) then
               read(4,*,end=13)ra,dec,az,wsys,wnoz,wcp,nbb,comp
               wred=wnoz+wcp-1.
               nbg(i)=nbb*comp
            elseif (idata.eq.6) then
c               read(4,*,end=13)ra,dec,az,nbb
               read(4,*,end=13)ra,dec,az,dum,nbb,comp,bias
               wsys=1.
               wred=1.
               nbg(i)=nbb*comp
            elseif (idata.eq.7 .or. idata.eq.8) then !PTH
 11            read(4,*,end=13)ra,dec,az,dum,wboss,wcp,wnoz,
     $          dum,dum,dum,dum,veto
               if (wboss*wcp*wnoz*veto.eq.0.) then !read another entry
                  goto 11
               endif   
               wsys=1.
               wred=wnoz+wcp-1.
               comp=1.
               nbb=nbar2(az)
               nbg(i)=nbb*comp
            elseif (idata.eq.9) then 
               read(4,*,end=13)ra,dec,az,nbb,wred
               wsys=1.
               comp=1.
               nbg(i)=nbb
            endif
            
            if (icomp.eq.1) then ! COMP=1 analysis
               comp=1.
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
 13      continue
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
         I10=0.d0
         I12=0.d0
         I12F=0.d0
         I22=0.d0
         I13=0.d0
         I23=0.d0
         I33=0.d0

         do i=1,Ngal
            nb=nbg(i)  
            weight=wg(i)
            If (ig(i).eq.1) then
               I10=I10+1.d0  
               I12=I12+dble(weight**2)  
               I22=I22+dble(nb*weight**2)
               I13=I13+dble(weight**3) 
               I23=I23+dble(nb*weight**3 ) 
               I33=I33+dble(nb**2 *weight**3) 
            endif
         enddo
         write(*,*)'these are normalization integrals from data'
         write(*,*)'I10=',I10
         write(*,*)'I12=',I12
         write(*,*)'I22=',I22
         write(*,*)'I13=',I13
         write(*,*)'I23=',I23
         write(*,*)'I33=',I33
       

         if (interpol.eq.2) then !CIC
            allocate(dcg(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcgw(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcgxx(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcgyy(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcgzz(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcgxy(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcgyz(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcgzx(Ngrid/2+1,Ngrid,Ngrid))
            call assign_CIC(Ngal,rg,rm,Lm,dcg,wg,ig,0,0) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxx,wg,ig,1,1) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyy,wg,ig,2,2) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzz,wg,ig,3,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgxy,wg,ig,1,2) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgyz,wg,ig,2,3) 
            call assign_CIC(Ngal,rg,rm,Lm,dcgzx,wg,ig,3,1) 
            do i=1,Ngal
               wg(i)=wg(i)**2
            enddo   
            call assign_CIC(Ngal,rg,rm,Lm,dcgw,wg,ig,0,0) 
c            write(*,*) 'assign done!'
c            write(*,*)'doing fft'
            call rfftwnd_f77_one_real_to_complex(plan_real,dcg,dcg)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgw,dcgw)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgxx,dcgxx)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgyy,dcgyy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgzz,dcgzz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgxy,dcgxy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgyz,dcgyz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcgzx,dcgzx)
c            write(*,*) 'FFT done!'
            call correct(Lm,Lm,Lm,dcg) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcgw) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcgxx) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcgyy) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcgzz) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcgxy) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcgyz) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcgzx) !correct density for interpolation
c            write(*,*) 'correction for interpolation done!'
         else !4th-order interlaced
            allocate(dcg(Ngrid,Ngrid,Ngrid))
            allocate(dcgw(Ngrid,Ngrid,Ngrid))
            allocate(dcgxx(Ngrid,Ngrid,Ngrid))
            allocate(dcgyy(Ngrid,Ngrid,Ngrid))
            allocate(dcgzz(Ngrid,Ngrid,Ngrid))
            allocate(dcgxy(Ngrid,Ngrid,Ngrid))
            allocate(dcgyz(Ngrid,Ngrid,Ngrid))
            allocate(dcgzx(Ngrid,Ngrid,Ngrid))
            call assign(Ngal,rg,rm,Lm,dcgxx,P0,nbg,ig,wg,1,1)
            call assign(Ngal,rg,rm,Lm,dcgyy,P0,nbg,ig,wg,2,2)
            call assign(Ngal,rg,rm,Lm,dcgzz,P0,nbg,ig,wg,3,3)
            call assign(Ngal,rg,rm,Lm,dcgxy,P0,nbg,ig,wg,1,2)
            call assign(Ngal,rg,rm,Lm,dcgyz,P0,nbg,ig,wg,2,3)
            call assign(Ngal,rg,rm,Lm,dcgzx,P0,nbg,ig,wg,3,1)
            call assign(Ngal,rg,rm,Lm,dcg,P0,nbg,ig,wg,0,0)
            do i=1,Ngal
               wg(i)=wg(i)**2
            enddo   
            call assign(Ngal,rg,rm,Lm,dcgw,P0,nbg,ig,wg,0,0)
c            write(*,*) 'assign done!'
c            write(*,*)'doing fft'
            call fftwnd_f77_one(planf,dcg,dcg)      
            call fftwnd_f77_one(planf,dcgw,dcgw)      
            call fftwnd_f77_one(planf,dcgxx,dcgxx)      
            call fftwnd_f77_one(planf,dcgyy,dcgyy)      
            call fftwnd_f77_one(planf,dcgzz,dcgzz)      
            call fftwnd_f77_one(planf,dcgxy,dcgxy)      
            call fftwnd_f77_one(planf,dcgyz,dcgyz)      
            call fftwnd_f77_one(planf,dcgzx,dcgzx)      
c            write(*,*) 'FFT done!'
            call fcomb(Lm,dcg,Ng)
            call fcomb(Lm,dcgw,Ng)
            call fcomb(Lm,dcgxx,Ng)
            call fcomb(Lm,dcgyy,Ng)
            call fcomb(Lm,dcgzz,Ng)
            call fcomb(Lm,dcgxy,Ng)
            call fcomb(Lm,dcgyz,Ng)
            call fcomb(Lm,dcgzx,Ng)
c            write(*,*) 'correction for interpolation done!'
         endif
         
         
      do 100 iz=1,Ngrid !build quadrupole
         ikz=mod(iz+Ngrid/2-2,Ngrid)-Ngrid/2+1
         do 100 iy=1,Ngrid
            iky=mod(iy+Ngrid/2-2,Ngrid)-Ngrid/2+1
            do 100 ix=1,Ngrid/2+1
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
     &                  -2.5*dcg(ix,iy,iz)  !quadrupole field
               end if
 100  continue


c         write(*,*) 'Fourier file :'
c         read(*,'(a)') filecoef
         call getarg(10,filecoef)
         open(unit=4,file=filecoef,status='unknown',form='unformatted')
         write(4)Lm
         write(4)(((dcg(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)real(I10),real(I12),real(I22),real(I13),real(I23),
     &   real(I33) !,real(I12F)
         write(4)P0,Ng,real(Ngsys),real(Ngsyscomp), real(Ngsystot)
         write(4)xmin,xmax,ymin,ymax,zmin,zmax
         write(4)(((dcgxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)(((dcgw(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         close(4)

       elseif (iflag.eq.1) then ! compute discretness integrals and FFT random mock
c         write(*,*)'Random Survey File'
c         read(*,'(a)')randomfile
         call getarg(9,randomfile)
         write(*,*)'Randomfile',randomfile
c         fname='/mount/chichipio2/hahn/data/'//randomfile
         allocate(rr(3,Nmax),nbr(Nmax),ir(Nmax),wr(Nmax),compr(Nmax))
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
         endif   
         do i=1,Nmax
            if (idata.eq.1) then !BOSS
c               read(4,*,end=15)ra,dec,az,nbb
               read(4,*,end=15)ra,dec,az,comp,nbb,wsys,wred
               if (wsys.ne.1. .or. wred.ne.1.) then
                  write(*,*)'randoms have bad systot weights',wsys,wred
                  stop
               endif   
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            elseif (idata.eq.2) then !LasDamas   
               read(4,*,end=15)ra,dec,az
               az=az/cspeed         ! ldg random has cz instead of z
               nbb=0.0000944233
               wsys=1.
               wred=1.
               comp=1.
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            elseif (idata.eq.11) then !Downsampled LasDamas   
               read(4,*,end=15)ra,dec,az
               nbb=nbar(az)
               wsys=1.
               wred=1.
               comp=1.
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            elseif (idata.eq.3 .or. idata.eq.4 .or. idata.eq.10) then !QPM and Nseries
 17            read(4,*,end=15)ra,dec,az,comp
               !17            read(4,*,end=15)ra,dec,az,comp,iveto
               !if (ifc.eq.1) then ! fiber colls + veto mask
               !   if (iveto.eq.1) then  
               !      goto 17 !in veto mask, read another entry
               !   endif   
               !endif 
               nbb=nbar(az)
               wsys=1.
               wred=1.
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            elseif (idata.eq.5) then
               read(4,*,end=15)ra,dec,az,nbb,comp
               wsys=1.
               wred=1.
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            elseif (idata.eq.6) then
               read(4,*,end=15)ra,dec,az,nbb,comp
               wsys=1.
               wred=1.
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            elseif (idata.eq.7 .or. idata.eq.8) then !PTH
 12            read(4,*,end=15)ra,dec,az,dum,wboss,wcp,wnoz,veto
               if (wboss*wcp*wnoz*veto.eq.0.) then !read another entry
                  goto 12
               endif   
               wsys=1.
               wred=wnoz+wcp-1.
               comp=1.
               nbb=nbar2(az)
               nbr(i)=nbb*comp ! number density as given in randoms (comp weighted)
            elseif (idata.eq.9) then 
               read(4,*,end=15)ra,dec,az,nbb
               wsys=1.
               wred=1.
               comp=1.
               nbr(i)=nbb
            endif
            
            Nran=Nran+1
            if (Nran.eq.Nmax) then
               write(*,*)'increase Nmax'
               stop
            endif   

            if (icomp.eq.1) then ! COMP=1 analysis
               comp=1.
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
         rcm=sqrt(xcm**2+ycm**2+zcm**2)
         thetaobs=acos(zcm/rcm)
         phiobs=acos(xcm/sqrt(xcm**2+ycm**2))
         write(*,*)'observer to CoM angles:',thetaobs,phiobs
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

         I10=0.d0
         I12=0.d0
         I22=0.d0
         I13=0.d0
         I23=0.d0
         I33=0.d0
         I12xx=0.d0
         I12yy=0.d0
         I12zz=0.d0
         I12xy=0.d0
         I12yz=0.d0
         I12zx=0.d0

         do i=1,Nran
            nb=nbr(i)  
            weight=wr(i)
            xr=rr(1,i)
            yr=rr(2,i)
            zr=rr(3,i)
            rsq=xr**2+yr**2+zr**2
            If (ir(i).eq.1) then
               I10=I10+1.d0 
               I12=I12+dble(weight**2)  
               I22=I22+dble(nb*weight**2)  
               I13=I13+dble(weight**3)  
               I23=I23+dble(nb*weight**3 )  
               I33=I33+dble(nb**2 *weight**3)  
               I12xx=I12xx+dble(weight**2*xr**2/rsq) 
               I12yy=I12yy+dble(weight**2*yr**2/rsq) 
               I12zz=I12zz+dble(weight**2*zr**2/rsq) 
               I12xy=I12xy+dble(weight**2*xr*yr/rsq) 
               I12yz=I12yz+dble(weight**2*yr*zr/rsq) 
               I12zx=I12zx+dble(weight**2*zr*xr/rsq) 
            endif
         enddo
         write(*,*)'(random) Inm to be scaled by alpha later'
         write(*,*)'I10=',I10
         write(*,*)'I12=',I12
         write(*,*)'I22=',I22
         write(*,*)'I13=',I13
         write(*,*)'I23=',I23
         write(*,*)'I33=',I33
         write(*,*)'I12xx/I12=',I12xx/I12
         write(*,*)'I12yy/I12=',I12yy/I12
         write(*,*)'I12zz/I12=',I12zz/I12
         write(*,*)'I12xy/I12=',I12xy/I12
         write(*,*)'I12yz/I12=',I12yz/I12
         write(*,*)'I12zx/I12=',I12zx/I12

         if (interpol.eq.2) then !CIC
            allocate(dcr(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcrw(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcrxx(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcryy(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcrzz(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcrxy(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcryz(Ngrid/2+1,Ngrid,Ngrid))
            allocate(dcrzx(Ngrid/2+1,Ngrid,Ngrid))
            call assign_CIC(Nran,rr,rm,Lm,dcr,wr,ir,0,0) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxx,wr,ir,1,1) 
            call assign_CIC(Nran,rr,rm,Lm,dcryy,wr,ir,2,2) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzz,wr,ir,3,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcrxy,wr,ir,1,2) 
            call assign_CIC(Nran,rr,rm,Lm,dcryz,wr,ir,2,3) 
            call assign_CIC(Nran,rr,rm,Lm,dcrzx,wr,ir,3,1) 
            do i=1,Nran
               wr(i)=wr(i)**2
            enddo   
            call assign_CIC(Nran,rr,rm,Lm,dcrw,wr,ir,0,0) 
c            write(*,*) 'assign done!'
c            write(*,*)'doing fft'
            call rfftwnd_f77_one_real_to_complex(plan_real,dcr,dcr)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrw,dcrw)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrxx,dcrxx)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcryy,dcryy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrzz,dcrzz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrxy,dcrxy)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcryz,dcryz)
            call rfftwnd_f77_one_real_to_complex(plan_real,dcrzx,dcrzx)
c            write(*,*) 'FFT done!'
            call correct(Lm,Lm,Lm,dcr) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcrw) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcrxx) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcryy) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcrzz) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcrxy) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcryz) !correct density for interpolation
            call correct(Lm,Lm,Lm,dcrzx) !correct density for interpolation
c            write(*,*) 'correction for interpolation done!'
         else !4th-order interlaced
            allocate(dcr(Ngrid,Ngrid,Ngrid))
            allocate(dcrw(Ngrid,Ngrid,Ngrid))
            allocate(dcrxx(Ngrid,Ngrid,Ngrid))
            allocate(dcryy(Ngrid,Ngrid,Ngrid))
            allocate(dcrzz(Ngrid,Ngrid,Ngrid))
            allocate(dcrxy(Ngrid,Ngrid,Ngrid))
            allocate(dcryz(Ngrid,Ngrid,Ngrid))
            allocate(dcrzx(Ngrid,Ngrid,Ngrid))
            call assign(Nran,rr,rm,Lm,dcrxx,P0,nbr,ir,wr,1,1)
            call assign(Nran,rr,rm,Lm,dcryy,P0,nbr,ir,wr,2,2)
            call assign(Nran,rr,rm,Lm,dcrzz,P0,nbr,ir,wr,3,3)
            call assign(Nran,rr,rm,Lm,dcrxy,P0,nbr,ir,wr,1,2)
            call assign(Nran,rr,rm,Lm,dcryz,P0,nbr,ir,wr,2,3)
            call assign(Nran,rr,rm,Lm,dcrzx,P0,nbr,ir,wr,3,1)
            call assign(Nran,rr,rm,Lm,dcr,P0,nbr,ir,wr,0,0)
            do i=1,Nran
               wr(i)=wr(i)**2
            enddo   
            call assign(Nran,rr,rm,Lm,dcrw,P0,nbr,ir,wr,0,0)
c            write(*,*) 'assign done!'
c            write(*,*)'doing fft'
            call fftwnd_f77_one(planf,dcr,dcr)      
            call fftwnd_f77_one(planf,dcrw,dcrw)      
            call fftwnd_f77_one(planf,dcrxx,dcrxx)      
            call fftwnd_f77_one(planf,dcryy,dcryy)      
            call fftwnd_f77_one(planf,dcrzz,dcrzz)      
            call fftwnd_f77_one(planf,dcrxy,dcrxy)      
            call fftwnd_f77_one(planf,dcryz,dcryz)      
            call fftwnd_f77_one(planf,dcrzx,dcrzx)      
c            write(*,*) 'FFT done!'
            call fcomb(Lm,dcr,Nr)
            call fcomb(Lm,dcrw,Nr)
            call fcomb(Lm,dcrxx,Nr)
            call fcomb(Lm,dcryy,Nr)
            call fcomb(Lm,dcrzz,Nr)
            call fcomb(Lm,dcrxy,Nr)
            call fcomb(Lm,dcryz,Nr)
            call fcomb(Lm,dcrzx,Nr)
c            write(*,*) 'correction for interpolation done!'
         endif

      do 101 iz=1,Ngrid !build quadrupole
         ikz=mod(iz+Ngrid/2-2,Ngrid)-Ngrid/2+1
         do 101 iy=1,Ngrid
            iky=mod(iy+Ngrid/2-2,Ngrid)-Ngrid/2+1
            do 101 ix=1,Ngrid/2+1
               ikx=mod(ix+Ngrid/2-2,Ngrid)-Ngrid/2+1
               rk=sqrt(float(ikx**2+iky**2+ikz**2))
               if(rk.gt.0.)then
                  kxh=float(ikx)/rk !unit vectors
                  kyh=float(iky)/rk
                  kzh=float(ikz)/rk
                     dcrxx(ix,iy,iz)=7.5*(dcrxx(ix,iy,iz)*kxh**2 
     &                  +dcryy(ix,iy,iz)*kyh**2
     &                  +dcrzz(ix,iy,iz)*kzh**2 
     &                  +2.*dcrxy(ix,iy,iz)*kxh*kyh
     &                  +2.*dcryz(ix,iy,iz)*kyh*kzh
     &                  +2.*dcrzx(ix,iy,iz)*kzh*kxh)
     &                  -2.5*dcr(ix,iy,iz)  !quadrupole field
               end if
 101  continue


c         write(*,*) 'Fourier file :'
c         read(*,'(a)') filecoef
         call getarg(10,filecoef)
c         outname='/mount/chichipio2/hahn/FFT/'//filecoef
         open(unit=4,file=filecoef,status='unknown',form='unformatted')
         write(4)Lm
         write(4)(((dcr(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)real(I10),real(I12),real(I22),real(I13),real(I23),
     &   real(I33)
         write(4)P0,Nr,real(Nrsys),real(Nrsyscomp),real(Nrsystot)
         write(4)thetaobs,phiobs
         write(4)xmin,xmax,ymin,ymax,zmin,zmax
         write(4)(((dcrxx(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         write(4)real(I12xx),real(I12yy),real(I12zz),real(I12xy),
     &   real(I12yz),real(I12zx)
         write(4)(((dcrw(ix,iy,iz),ix=1,Lm/2+1),iy=1,Lm),iz=1,Lm)
         close(4)
         
      endif
      

 1025 format(2x,6e14.6)
 123  stop
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                                                                                                                                                        
      REAL function nbar2(QQ) !nbar(z) for PTHalos mocks                                                                                                                                                                                                                                                                                         
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      implicit none
      integer Nsel2
      parameter(Nsel2=10)      
      real z2(Nsel2),selfun2(Nsel2),sec2(Nsel2),qq,az
      common /interpol2/z2,selfun2,sec2
      az=QQ
      call splint(z2,selfun2,sec2,Nsel2,az,nbar2)
      RETURN
      END
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
      subroutine assign_CIC(Npar,r,rm,Lm,dtl,wg,ig,ia,ib)
cc*******************************************************************
      implicit none
      integer Lm,ixp,iyp,izp,ixa,iya,iza,i,ix,iy,iz,ia,ib,Npar
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

         if (ia.eq.0 .and. ib.eq.0) then !FFT delta
            we=wg(i)
         else !FFT Qij
            rnorm=r(1,i)**2+r(2,i)**2+r(3,i)**2
            we=wg(i)*r(ia,i)*r(ib,i)/rnorm
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

C      if (ia.gt.0 .and. ib.gt.0) then  !FFT Qij
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
      subroutine assign(N,r,rm,L,dtl,P0,nbg,ig,w,ia,ib) !with FKP weighing
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
      integer ia,ib
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
       if (ia.eq.0 .and. ib.eq.0) then !FFT delta
          we=w(i)
       else !FFT Qij
          rnorm=r(1,i)**2+r(2,i)**2+r(3,i)**2
          we=w(i)*r(ia,i)*r(ib,i)/rnorm
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
