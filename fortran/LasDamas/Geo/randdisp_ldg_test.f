c Applies peak corrections to fibercollided Mock Galaxy Catalogs
c Generates artificial galaxy dLOS away from galaxy with weight > 1 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c code input: sigma_peak, nbar file, galaxy file, outputfile
c galaxy file input: ra,dec,z,wcp
      implicit none
      integer i,j,k,ii,jj,kk,n,nn
      integer Ngal,Nmax
      integer Nnz,Ntail,wNgal,cpcount,cppeakcount,cptailcount
      integer Ng,l,ipoly,wb,wcp,wred,flag
      real pi,cspeed,Om0,OL0,redtru,m1,m2,veto
      REAL zt,zlo,zhi,garb1,garb2,garb3
      real cpz,cpnbarz
      REAL sigma,peakfrac
      parameter(Nmax=2*10**8,pi=3.141592654)
      parameter(Om0=0.27,OL0=0.73)
      parameter(cspeed=299800.0)
      real gfrac,pr,wsys,wwsys
      real, allocatable :: wg(:),wr(:)
      REAL, ALLOCATABLE :: out_ra(:),out_dec(:),out_az(:)
      REAL, ALLOCATABLE :: out_wboss(:),out_wred(:),out_wcp(:)
      real, allocatable :: out_veto(:)
      REAL, ALLOCATABLE :: z(:),dm(:),dmsec(:)
      REAL, ALLOCATABLE :: tailz(:),tailnbarz(:),tailsec(:)
      REAL az,ra,dec,rad,numden
      real ran1,ran2,ran3,pz
      real chi,peakpofr,tailnbar,dmtoz
      character lssfile*200,peakfile*200,nbarfile*200,tailnbarfile*200
      CHARACTER outputfile*200,sigmastr*200,peakfracstr*200
      external chi,peakpofr,tailnbar,dmtoz

      Nnz=0
      allocate(dm(Nmax),z(Nmax),dmsec(Nmax))
      call getarg(1,nbarfile)
      open(unit=3,file=nbarfile,status='old',form='formatted')
      do l=1,Nmax
        read(3,*,end=12) zt,zlo,zhi,numden,garb1,garb2,garb3
        z(l)=zt
        dm(l)=chi(zt)
        Nnz=Nnz+1
      enddo
  12  continue
      close(3)
      call spline(dm,z,Nnz,3e30,3e30,dmsec)
      
      CALL getarg(2,lssfile)
      allocate(wg(Nmax))
      allocate(out_ra(Nmax),out_dec(Nmax),out_az(Nmax),out_wcp(Nmax))
      open(unit=4,file=lssfile,status='old',form='formatted')
      write(*,*) lssfile

      Ngal=0 
      wNgal=0
      wsys=0.0
      wwsys=0.0
      cpcount=0
      cppeakcount=0
      cptailcount=0
      CALL RANDOM_SEED
      DO i=1,Nmax
        READ(4,*,END=13)ra,dec,az
            wg(i)=1.0 

            az=az/cspeed
            CALL RANDOM_NUMBER(ran1) 

            out_ra(i)=ra
            out_dec(i)=dec
            out_wcp(i)=1.0
            IF (ran1 .GT. 0.03) THEN    ! hardcoded fraction of randomly displaced galaxies 
                out_az(i)=az 
            ELSE
                CALL RANDOM_NUMBER(ran2) 
                out_az(i)=0.16+ran2*0.28        ! randomly displace
            ENDIF 
            wsys=wsys+wg(i)
            Ngal=Ngal+1
         enddo
 13      continue
         close(4)
         WRITE(*,*) 'Ngal,sys=',wsys/float(Ngal)
        
         CALL getarg(3,outputfile)
         
         OPEN(unit=8,file=outputfile,status='unknown',form='formatted') 
         DO i=1,Ngal
             WRITE(8,'(4(E,2x))') out_ra(i),out_dec(i),
     &          out_az(i),out_wcp(i)
         ENDDO 
         CLOSE(8)
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function peakpofr(QQQ, sigsig) !peakpofr(r)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      REAL sigma, sigsig
      REAL peakar,QQQ
      peakar=QQQ
      sigma=sigsig 
      peakpofr=EXP(-1.0*ABS(peakar)/sigma) 
      RETURN
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function tailnbar(QQQ)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      REAL tailar, QQQ
      tailar=QQQ
      tailnbar=9.44451e-5
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
