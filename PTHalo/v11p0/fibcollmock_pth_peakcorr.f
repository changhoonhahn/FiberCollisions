c Applies peak corrections to PTHalo Mock Galaxy Catalogs
c Generates artificial galaxy dLOS away from galaxy with w_cp > 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c code input: sigma_peak, nbar file, galaxy file, outputfile
c galaxy file input: ra,dec,z,ipoly,wboss,wcp,wrf,redtru,flag,m1,m2,veto
c output file output: ra,dec,z,ipoly,wboss,wcp,wrf,redtru,flag,m1,m2,veto
      implicit none
      integer i,j,k,ii,jj,kk,n,nn
      integer Ngal,Nmax
      integer Nnz,Ntail,wNgal,cpcount,cppeakcount,cptailcount
      integer Ng,l,ipoly,wb,wcp,wred,flag
      real pi,cspeed,Om0,OL0,redtru,m1,m2,veto
      REAL zt,zlo,zhi,garb1,garb2,garb3
      real cpz,cpnbarz
      REAL sigma
      parameter(Nmax=2*10**8,pi=3.141592654)
      parameter(Om0=0.27,OL0=0.73)
      parameter(cspeed=299800.0)
      real gfrac,pr,wsys,wwsys
      real, allocatable :: wg(:),wr(:)
      REAL, ALLOCATABLE :: out_ra(:),out_dec(:),out_az(:)
      REAL, ALLOCATABLE :: out_wboss(:),out_wred(:),out_wcp(:)
      real, allocatable :: out_ipoly(:),out_redtru(:),out_flag(:)
      real, allocatable :: out_m1(:),out_m2(:),out_veto(:)
      real, allocatable :: wwg(:)
      REAL, ALLOCATABLE :: z(:),dm(:),dmsec(:)
      REAL az,ra,dec,rad,numden
      real ran1,ran2,ran3,pz
      real chi,peakpofr,tailnbar,dmtoz
      character lssfile*200,peakfile*200,nbarfile*200,tailnbarfile*200
      CHARACTER outputfile*200,sigmastr*200
      external chi,peakpofr,tailnbar,dmtoz

      Nnz=0 
      ALLOCATE(dm(Nmax),z(Nmax),dmsec(Nmax))
      call getarg(1, sigmastr) 
      read(sigmastr, *) sigma
      CALL getarg(2,nbarfile)
      open(unit=3,file=nbarfile,status='old',form='formatted')
      do l=1,Nmax
            read(3,*,end=12) zt,zlo,zhi,numden,garb1,garb2,garb3
            z(l)=zt
            dm(l)=chi(zt)
            Nnz=Nnz+1 
      enddo
 12   continue
      close(3)
      call spline(dm,z,Nnz,3e30,3e30,dmsec)
      
      CALL getarg(3,lssfile)
      allocate(wg(Nmax),wwg(Nmax))
      allocate(out_ra(Nmax),out_dec(Nmax),out_az(Nmax),out_ipoly(Nmax))
      allocate(out_wboss(Nmax),out_wred(Nmax),out_wcp(Nmax))
      allocate(out_redtru(Nmax),out_flag(Nmax))
      allocate(out_m1(Nmax),out_m2(Nmax),out_veto(Nmax))
      open(unit=4,file=lssfile,status='old',form='formatted')

      Ngal=0 
      wNgal=0
      wsys=0.0
      wwsys=0.0
      cpcount=0
      cppeakcount=0
      cptailcount=0
      CALL RANDOM_SEED
      DO i=1,Nmax
        READ(4,*,END=13)ra,dec,az,ipoly,wb,wcp,wred,redtru,flag,m1
     &      ,m2,veto
            IF (wb.gt.0 .and. wcp.gt.0 .and. wred.gt.0 .and. veto.gt.0) 
     &      THEN
                wwg(i)=float(wb)*(float(wcp)+float(wred)-1.0)
            ELSE
                wwg(i)=0.0
            ENDIF
            wNgal=wNgal+1
            wwsys=wwsys+wwg(i)

            IF (wcp.eq.1 .and. wb.gt.0 .and. wred.gt.0 .and. veto.gt.0) 
     &      THEN 
                out_ra(Ngal+1)=ra
                out_dec(Ngal+1)=dec
                out_az(Ngal+1)=az
                out_ipoly(Ngal+1)=ipoly
                out_wboss(Ngal+1)=wb
                out_wred(Ngal+1)=wred
                out_wcp(Ngal+1)=wcp
                out_redtru(Ngal+1)=redtru
                out_flag(Ngal+1)=flag
                out_m1(Ngal+1)=m1
                out_m2(Ngal+1)=m2
                out_veto(Ngal+1)=veto
                rad=chi(az)
                wg(Ngal+1)=float(wb)*(float(wcp)+float(wred)-1.0)
                wsys=wsys+wg(Ngal+1)
                Ngal=Ngal+1
            ELSEIF (wcp.gt.1 .and. wb.gt.0 .and. wred.gt.0 .and. 
     &      veto.gt.0) THEN
                cpcount=cpcount+1
                out_ra(Ngal+1)=ra
                out_dec(Ngal+1)=dec
                out_az(Ngal+1)=az
                out_ipoly(Ngal+1)=ipoly
                out_wboss(Ngal+1)=float(wb)
                out_wred(Ngal+1)=float(wred)
                out_wcp(Ngal+1)=1.0
                out_redtru(Ngal+1)=redtru
                out_flag(Ngal+1)=flag
                out_m1(Ngal+1)=m1
                out_m2(Ngal+1)=m2
                out_veto(Ngal+1)=veto
                rad=chi(az)
                wg(Ngal+1)=float(wb)*float(wred)
                wsys=wsys+wg(Ngal+1)
                Ngal=Ngal+1
                
                DO WHILE (wcp .gt. 1)
                    CALL RANDOM_NUMBER(ran1) 
                    CALL RANDOM_NUMBER(ran2)
                    CALL RANDOM_NUMBER(ran3)
                    cppeakcount=cppeakcount+1
                    ran2=(-3.0+ran2*6.0)*sigma
                    pr=peakpofr(ran2,sigma)
                    DO WHILE (pr .le. ran1)
                        CALL RANDOM_NUMBER(ran1)
                        CALL RANDOM_NUMBER(ran2)
                        ran2=(-3.0+ran2*6.0)*sigma
                        pr=peakpofr(ran2,sigma)
                    ENDDO
                    rad=rad+ran2
                    
                    wcp=wcp-1 
                    wg(Ngal+1)=float(wb)
                    wsys=wsys+wg(Ngal+1)
                    out_ra(Ngal+1)=ra
                    out_dec(Ngal+1)=dec
                    out_az(Ngal+1)=dmtoz(rad,Nnz,dm,z,dmsec)
                    out_ipoly(Ngal+1)=ipoly
                    out_wboss(Ngal+1)=float(wb)
                    out_wred(Ngal+1)=1.0
                    out_wcp(Ngal+1)=1.0
                    out_redtru(Ngal+1)=redtru
                    out_flag(Ngal+1)=flag
                    out_m1(Ngal+1)=m1
                    out_m2(Ngal+1)=m2
                    out_veto(Ngal+1)=veto
                    Ngal=Ngal+1
                ENDDO 
            ELSE 
                out_ra(Ngal+1)=ra
                out_dec(Ngal+1)=dec
                out_az(Ngal+1)=az
                out_ipoly(Ngal+1)=ipoly
                out_wboss(Ngal+1)=wb
                out_wred(Ngal+1)=wred
                out_wcp(Ngal+1)=wcp
                out_redtru(Ngal+1)=redtru
                out_flag(Ngal+1)=flag
                out_m1(Ngal+1)=m1
                out_m2(Ngal+1)=m2
                out_veto(Ngal+1)=veto
                wg(Ngal+1)=0.0
                wsys=wsys+wg(Ngal+1)
                Ngal=Ngal+1
            ENDIF
         enddo
 13      continue
         close(4)
         WRITE(*,*) 'Total Number of Close Pairs=',cpcount
         WRITE(*,*) 'Peak=',cppeakcount,'Tail=',cptailcount
         WRITE(*,*) 'Normal Wsys=',wwsys,'Normal Ngal=',wNgal
         WRITE(*,*) 'Corrected Wsys=',wsys,'Corrected Ngal=',Ngal
         WRITE(*,*) 'Ngal,sys=',wsys/float(Ngal)
        
         call getarg(4,outputfile)
         
         OPEN(unit=8,file=outputfile,status='unknown',form='formatted') 
         DO i=1,Ngal
             WRITE(8,'(12(E,2x))') out_ra(i),out_dec(i),out_az(i),
     &          out_ipoly(i),out_wboss(i),out_wcp(i),out_wred(i),
     &          out_redtru(i),out_flag(i),out_m1(i),out_m2(i),
     &          out_veto(i)
         ENDDO 
         CLOSE(8)
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function peakpofr(QQQ,sigsig) !peakpofr(r)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      REAL sigma,sigsig
      REAL peakar,QQQ
      peakar=QQQ
      sigma=sigsig
      peakpofr=EXP(-1.0*ABS(peakar)/sigma) 
      RETURN
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function dmtoz(QQQ,N,dm,z,dmsec)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      INTEGER N
      REAL dm(N),z(N),dmsec(N)
      REAL dmself,dmar,QQQ
      dmar=QQQ
      CALL splint(dm,z,dmsec,N,dmar,dmself)
      dmtoz=dmself
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
      real Om0,OL0
      parameter (Om0=0.27,OL0=0.73)  
      real*8 z
      rdi=3000.d0/dsqrt(OL0+(1.d0-Om0-OL0)*(1.d0+z)**2+Om0*(1.d0+z)**3)
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      include '/home/users/hahn/powercode/dqage.f'
      include '/home/users/hahn/powercode/d1mach.f'
      include '/home/users/hahn/powercode/dqawfe.f'
      include '/home/users/hahn/powercode/spline.f'
      include '/home/users/hahn/powercode/splint.f'
      include '/home/users/hahn/powercode/gabqx.f'
