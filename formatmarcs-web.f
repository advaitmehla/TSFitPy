      PROGRAM storemarcs
* 2012-04    Fixed silly output data (large mu & neg: anconv, Pe, Pg)
* 2011-07-12 Changed Pg -> Pg-Pe in the krz models
* 2010-07-06 changed name field delimiter from : to _ for Windows users /BE
* 2009-02-02 corrected .krz abundance scale
* will create gzipped .mod .flx .krz and .opa files for models in file "donames"
* these are put into the temporary directory 'base/data' for checking before publication
* ifort -convert big_endian -save -O3 -fpe3 -ftz -axT -o storemarcs storemarcs.f
* run by storemarcs < donames
* 2008-01-11 store also the right column mass (RHOX) scale in the .mod files
* NOTE: cancel check of file names around line 415 if file names are 'tweaked'
*
      parameter(ndp=100)
      parameter (nmolmax=30)
      CHARACTER modname*24,filepub*72,filemod*86,fileflx*86,filekrz*86
      character filewav*100
      character fileopa*86,filebin*88,cpopa*180,string*256
      character allmodname*73,dofile*69,metallicity*5,alpha*5
      character zipcom*91,dag*10,dirname*43,abufile*57,xyz*25
      dimension tau(ndp),t(ndp),z(ndp),ptot(ndp),prad(ndp)
      dimension pg(ndp),pturb(ndp),pe(ndp),ro(ndp),xkapr(ndp)
      dimension rr(ndp),taus(ndp),xlr(30)
      dimension xlb(155000),w(155000),fluxme(155000),fluxmecont(155000)
      real emu(ndp),cp(ndp),cv(ndp),agrad(ndp),q(ndp),u(ndp),v(ndp)
      real anconv(ndp),hnic(ndp),psum(ndp),other(ndp),rhox(ndp)
      real TEFF,FLUX,PALFA,PNY,PY,PBETA,inputabund(92),xite,rmass,cc
      real xabund,yabund,zabund,xcalc,ycalc,zcalc
      real*8 dluminosity
      integer micro,nstore
      logical sph,erroneous
      common /struct/ teff,tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /termodyn/ emu,cp,cv,agrad,q,u,v,anconv,hnic
      common /pressure/ presmo(nmolmax,ndp),ptio(ndp),nmol
      common /spectr/ nlb,xlb,w,fluxme,fluxmecont
      common /binstruc/bPPR(NDP),bPPT(NDP),bPP(NDP),bGG(NDP),
     & bZZ(NDP),bDD(NDP),
     & bVV(NDP),bFFC(NDP),bPPE(NDP),bTT(NDP),
     & bTAULN(NDP),NbTAU,IbTER,erad(ndp)
      logical lwarning
      data lwarning / .false. /
*
      print*,'MARCS file to convert'
      read(5,'(a)') allmodname
      print*,'mass? (0.0 for PP)'
      read(5,*) rmass
      print*,'microturb? (km/s)'
      read(5,*) xite
      print*,'metallicity? (log)'
      read(5,*) metallicity
      print*,'alpha? (log)'
      read(5,*) alpha
      filepub=allmodname
      sph=.true.
      OPEN(UNIT=10,FILE=allmodname,STATUS='OLD',FORM='UNFORMATTED',
     &     RECL=152600)
      call READMO(10,NDEPTH,GG,RADIUS,modname,
     &            FLUX,PALFA,PNY,PY,PBETA,dluminosity,inputabund,dag)
      close(10)
      if(radius.le.1.0) then
        sph=.false.
        rmass=0.0
      endif
* compute RHOX the integrated column mass from the surface to here
ccc* store only 56 depth points
ccc      nstore=56
      nstore=ndepth
      if(ndepth.lt.nstore) then
        stop 'main: ndepth < 56 !!!'
      endif
      do k=1,ndepth
        rhox(k)=ptot(k)/gg
      enddo
      write(filemod,1010) trim(filepub)
 1010 format('./',a,'.mod')
      write(fileflx,1020) trim(filepub)
 1020 format('./',a,'.flx')
      write(filekrz,1030) trim(filepub)
 1030 format('./',a,'.krz')
      write(fileopa,1040) trim(filepub)
 1040 format('./',a,'.opa')
      write(filebin,1050) trim(filepub)
 1050 format('./',a,'.bin')
      write(filewav,1051) trim(filepub)
 1051 format('./',a,'.wav')
*
* write the model file
      close(20)
      OPEN(UNIT=20,FILE=filemod,form='formatted',STATUS='UNKNOWN')
      write(*,'(a)') trim(filepub)
      write(20,'(a)') trim(filepub)
      write(*,1170)  TEFF,dag(1:8),FLUX,GG,xite
      write(20,1170) TEFF,dag(1:8),FLUX,GG,xite
 1170   format(  f7.0,'      Teff [K].         ',
     &                      'Last iteration; yyyymmdd=',a8,/,
     &         1p,e12.4,0p,' Flux [erg/cm2/s]',/,
     &         1p,e12.4,0p,' Surface gravity [cm/s2]',/,
     &         f5.1,'        Microturbulence parameter [km/s]')
      if(sph) then
* spherical
        write(*, 1171) rmass
        write(20,1171) rmass
 1171   format(f5.1,'        Mass [Msun]')
      else
* plane-parallel
        write(*, 1172) rmass
        write(20,1172) rmass
 1172   format(f5.1,'        No mass for plane-parallel models')
      endif
      write(*, 1173)  metallicity,alpha
      write(20,1173)  metallicity,alpha
 1173 format(x,a5,x,a5,  ' Metallicity [Fe/H] and [alpha/Fe]')
      if(sph) then
* spherical
        write(*, 1174)  RADIUS,dluminosity
        write(20,1174)  RADIUS,dluminosity
 1174   format(1p,e12.4,0p,' Radius [cm] at Tau(Rosseland)=1.0',/,
     &         f12.5,      ' Luminosity [Lsun]')
      else
* plane-parallel
        write(*, 1175)  RADIUS,dluminosity
        write(20,1175)  RADIUS,dluminosity
 1175   format(1p,e12.4,0p,' 1 cm radius for plane-parallel models',/,
     &         1p,e12.4,0p,' Luminosity [Lsun] FOR A RADIUS OF 1 cm!')
      endif
      write(*, 1176) PALFA,PNY,PY,PBETA
      write(20,1176) PALFA,PNY,PY,PBETA
 1176 format(f6.2,f5.2,f6.3,f5.2,
     &       ' are the convection parameters: alpha, nu, y and beta')
      call sxyz(92,inputabund,xcalc,ycalc,zcalc)
        write(xyz,2222) xcalc,ycalc,zcalc
2222    format(f6.4,1x,f6.4,1x,f9.7)
        write(*, 3177) xyz
        write(20,3177) xyz
 3177   format(x,a25,' are X, Y and Z')
c      if(cc.eq.-0.13) then
c        write(*, 1177) xyz
c        write(20,1177) xyz
c 1177   format(x,a25,' are X, Y and Z, 12C/13C=20')
c      else if(cc.eq.-0.38) then
c        write(*, 1178) xyz
c        write(20,1178) xyz
c 1178   format(x,a25,' are X, Y and Z, 12C/13C=4')
c      else if(cc.eq.+0.00 .or. cc.eq.+0.08) then
c        write(*, 1179) xyz
c        write(20,1179) xyz
c 1179   format(x,a25,' are X, Y and Z, 12C/13C=89 (=solar)')
c      else
c        stop 'carbon abundance not expected'
c      endif
      write(20,1180) inputabund
 1180 format('Logarithmic chemical number abundances, H always 12.00',
     &        10(/,10f7.2))
      if(dluminosity.ge.1.e6) then
        lwarning=.true.
      endif
cc      read(xyz,*) xabund,yabund,zabund
cc      if(abs(xabund-xcalc)/xcalc.gt.0.000015) then
cc        print 1710,xabund,xcalc
cc 1710   format('WARNING: X abundances may differ, read, calc=',2f8.5)
cc        stop 'storemarcs: erroneous model'
cc      endif
cc      if(abs(yabund-ycalc)/ycalc.gt.0.000045) then
cc        print 1720,yabund,ycalc
cc 1720   format('WARNING: Y abundances may differ, read, calc=',2f8.5)
cc        stop 'storemarcs: erroneous model'
cc      endif
cc      if(abs(zabund-zcalc)/zcalc.gt.0.005) then
cc        print 1730,zabund,zcalc
cc 1730   format('WARNING: Z abundances may differ, read, calc=',2f14.10)
cc        stop 'storemarcs: erroneous model'
cc      endif
*
* Fix common ugly "unphysicals" in the extreme points **2012-04-03*******
      if(pe(1).lt.0.0) pe(1)=-pe(1)                                     *
      if(pg(1).lt.0.0) pg(1)=-pg(1)                                     *
      if(pe(nstor).lt.0.0) pe(nstor)=-pe(nstor)                                  *
      if(pg(nstor).lt.0.0) pg(nstor)=-pg(nstor)                                  *
      if(ro(nstor).lt.0.0) ro(nstor)=-ro(nstor)                                  *
      if(v(nstor).lt.0.0) v(nstor)=-v(nstor)                                     *
      if(rhox(nstor).lt.0.0) rhox(nstor)=-rhox(nstor)                            *
*************************************************************************
* structure
      write(20,1230) nstore
 1230 format(i4,' Number of depth points')
      write(20,1232)
 1232 format('Model structure')
      write(20,1234)
 1234 format(' k lgTauR  lgTau5    Depth     T  ',
     & '      Pe          Pg         Prad       Pturb')
      do k=1,nstore
        if(pe(k).lt.0.0) then
          write(89,8810) filepub(1:len_trim(filepub))
 8810     format('Negative Pe in ',a)
          erroneous=.true.
        endif
        if(pg(k).lt.0.0) then
          write(89,8820) filepub(1:len_trim(filepub))
 8820     format('Negative Pg in ',a)
          erroneous=.true.
        endif
        write(20,200) k,alog10(tau(k)),alog10(taus(k)),z(k),t(k),
     &                pe(k),pg(k),prad(k),pturb(k)
  200   format(i3,f6.2,f8.4,1p,e11.3,0p,f8.1,1p,5e12.4,0p)
* Format changed 2012-04-03
      enddo
* thermodynamics
      write(20,1235)
 1235 format(' k lgTauR    KappaRoss   Density   Mu      Vconv',
     &       '   Fconv/F      RHOX')
      do k=1,nstore
        if(emu(k).gt.2.50) then
          write(88,8830) emu(k),k,filepub(1:len_trim(filepub))
 8830     format('Unreasonable mu:',f6.1,' k=',i3,' put=2.5 in ',a)
          emu(k)=2.50
        endif
        if(anconv(k).lt.0.0) then
          write(88,8840) anconv(k),k,filepub(1:len_trim(filepub))
 8840     format('Negative anconv:',f7.3,' k=',i3,' put=0.0 in ',a)
          anconv(k)=0.0
        endif
        write(20,210) k,alog10(tau(k)),xkapr(k),ro(k),emu(k),
     &                v(k),anconv(k),rhox(k)
  210   format(i3,f6.2,1p,2e12.4,0p,f6.3,1p,e11.3,0p,f8.5,1p,e14.6,0p)
* Format changed 2012-04-03
      enddo
* molecular partial pressures
      write(20,1233)
 1233 format('Assorted logarithmic partial pressures')
      write(20,1236)
 1236 format(' k  lgPgas   H I    H-     H2     H2+    H2O',
     & '    OH     CH     CO     CN     C2  ')
      do k=1,nstore
        psum(k)=pe(k)+hnic(k)+ptio(k)
        do j=1,16
          psum(k)=psum(k)+presmo(j,k)
        enddo
        do j=18,nmol
          psum(k)=psum(k)+presmo(j,k)
        enddo
        write(20,220) k,alog10(pg(k)),alog10(hnic(k)),
     &                (alog10(presmo(j,k)),j=1,9)
  220   format(i3,f7.3,11f7.2)
      enddo
      write(20,1237)
 1237 format(' k    N2     O2     NO     NH     TiO',
     & '   C2H2    HCN    C2H    HS     SiH    C3H')
      do k=1,nstore
        write(20,222) k,
     &                (alog10(presmo(j,k)),j=10,13),alog10(ptio(k)),
     &                (alog10(presmo(j,k)),j=14,16),
     &                (alog10(presmo(j,k)),j=18,20)
  222   format(i3,15f7.2)
      enddo
      write(20,1238)
 1238 format(' k    C3     CS     SiC   SiC2    NS  ',
     & '   SiN    SiO    SO     S2     SiS   Other')
      do k=1,nstore
        other(k)=pg(k)-psum(k)
        if(other(k).gt.0.0) then
          write(20,224) k,
     &         (alog10(presmo(j,k)),j=21,30),alog10(other(k))
  224     format(i3,12f7.2,f6.3)
        else
          write(20,225) k,(alog10(presmo(j,k)),j=21,30)
  225     format(i3,10f7.2,'   -99.')
        endif
      enddo
      close(20)
      if(erroneous) then
        string='/bin/rm '//filemod
        call system(string)
        write(89,*) filemod(1:len_trim(filemod)),' not written'
        erroneous=.false.
      else
* gzip the model file
cc        zipcom='gzip '//filemod
cc        call system(zipcom)
      endif
*
* just copy the binary file
*     string='cp '//allmodname//' '//filebin
*     call system(string)
      print *,'binary copying disabled, line 325'
*
* check wavelength set:
*
cc      close(26)
      OPEN(UNIT=26,FILE=filewav,form='formatted',
     &     status='unknown')
      do i=1,nlb
        write(26,*) xlb(i)
      enddo
      close(26)
      wfirst=xlb(1)
      wlast=xlb(nlb)
*     print *,abufile,'*****************************************'
*
c      do i=1,nlb
c        if(xlb(i).ge.1299.9) then
c          n=i
c          goto 70
c        endif
c      enddo
c      stop 'Starting wavelength not found'
c   70 continue
c      if(abs(xlb(n)/wfirst-1.0).ge.1.e-5) then
c        print *,'Wavelengths slip:',xlb(n),wfirst
c        stop
c      endif
c      if(abs(xlb(nlb)/wlast-1.0).ge.1.e-5) then
c        print *,'Wavelengths slip:',xlb(nlb),wlast
c        stop
c      endif
      print*,'wavelengths between ',wfirst,' and ',wlast
*
******** Just constructing the wavelength file (once and for all)*******
*     OPEN(UNIT=25,FILE='base/wavelengths.wav',form='formatted',
*    &     STATUS='UNKNOWN')
*     do i=n,nlb
*       if(xlb(i).lt.10000.) then
*         write(25,251) xlb(i)
* 251     format(f9.4)
*       else if(xlb(i).lt.100000.) then
*         write(25,252) xlb(i)
* 252     format(f9.3)
*       else if(xlb(i).lt.1000000.) then
*         write(25,253) xlb(i)
* 253     format(f9.2)
*       else
*         write(25,254) xlb(i)
* 254     format(f9.1)
*       endif
*     enddo
*     close(25)
************************************************************************
*
* write the flux data file:
      OPEN(UNIT=30,FILE=fileflx,form='formatted',STATUS='UNKNOWN')
      do i=1,nlb
        write(30,300) fluxme(i)
  300   format(1p,e11.5,0p)
      enddo
      close(30)
* gzip the flux file
c      write(zipcom,400) fileflx
c  400 format('gzip ',a86)
c      call system(zipcom)
* make and gzip a "Kurucz" type file 
      call marcs35_2krz(allmodname,filekrz,sph,xite,inputabund,nstore)
c      write(zipcom,400) filekrz
c      call system(zipcom)
* copy the opacity data (babsma) file:
      if(modname(6:6).eq.'-') then
* this is a spherical model with negative log g
        mlen=21
      else
        if(modname(15:15).eq.'m') then
* this is a spherical model with positive log g
          mlen=20
        else
* this is a plane parallel model (always with positive log g)
          mlen=16
        endif
      endif
      if(xite.ge.9.6) mlen=mlen+1
 
      if(lwarning) then
        stop 'WARNING model format: Luminosity >= 1.e6'
        stop 'storemarcs: erroneous format'
      endif
 
*
  99  continue
      print *,'Problems: files fort.88 (small) & fort.89 (fatal)'
      stop
      END
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
C
      SUBROUTINE READMO(IARCH,JTAU,G,RADIUS,modname,
     &           FLUX,PALFA,PNY,PY,PBETA,dluminosity,inputabund,dag)
C        THIS ROUTINE READS ONE MODEL
C             ( All features taken from listmo )
      PARAMETER (NDP=100)
C
      DIMENSION ABUND(16),TKORRM(NDP),FCORR(NDP),TAU(NDP),TAUS(NDP),
     *T(NDP),PE(NDP),PG(NDP),PRAD(NDP),PTURB(NDP),XKAPR(NDP),RO(NDP),
     *CP(NDP),CV(NDP),AGRAD(NDP),Q(NDP),U(NDP),V(NDP),ANCONV(NDP),
     *PRESMO(30,NDP),RR(NDP),Z(NDP),EMU(NDP),HNIC(NDP)
     *,NJ(16),XLR(30),IEL(16),ANJON(16,5),PART(16,5),PROV(50,20+1),
     *ABSKA(50),SPRIDA(50),XLB(155000),FLUXME(155000),
     & FLUXMEcont(155000),
     * ABNAME(50),SOURCE(50),PTOT(NDP)
*     dimension FCONV(NDP),FLUMAG(155000),PEP(16)
      DIMENSION W(155000),UW(12),BW(21),VW(25)
      CHARACTER*10 DAG,NAME,NAMEP,KLOCK
      CHARACTER*8 ABNAME,SOURCE
      character*24 realname,modname
*     DIMENSION WAVFLX(10)
      dimension PTIO(NDP)
      real*8 dluminosity
      real abSc,abTi,abV,abMn,abCo
      real inputabund(92)
      common /binstruc/ dummy(11*ndp+2),erad(ndp)
      common /struct/ teff,tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /termodyn/ emu,cp,cv,agrad,q,u,v,anconv,hnic
      common /pressure/ presmo,ptio,nmol
      common /spectr/ nlb,xlb,w,fluxme,fluxmecont
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/,BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/,VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      DATA NAME/'LOCAL'/,NAMEP/'PARSONS'/
      DATA A,B/.34785485,.65214515/
      IREAD=5
C
      REWIND IARCH
      READ(IARCH) DUMMY,erad
      READ(IARCH) INORD,DAG,KLOCK
      READ(IARCH) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
*    &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            realname,
     &            ITMAX,NEL,(ABUND(I),I=1,NEL),abSc,abTi,abV,abMn,abCo
      modname=realname
      GLOG=ALOG10(G)
      FNORD=0.1*INORD
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX
      DO 2 I=1,NEL
    2 ABUND(I)=ALOG10(ABUND(I))+12.
      READ(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      READ(IARCH)JTAU,(TKORRM(I),I=1,JTAU),(FCORR(K),K=1,JTAU)
      NTPO=0
      DO 3 K=1,JTAU
        READ(IARCH) KR,TAU(K),TAUS(K),Z(K),T(K),PE(K),PG(K),PRAD(K),
     &              PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &              AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &              (PRESMO(J,K),J=1,NMOL),ptio(k)
        TAUK=ALOG10(TAU(K))+10.01
        KTAU=TAUK
        IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
        IF(KTAU.EQ.10) K0=K
        NTPO=NTPO+1
   31   CONTINUE
    3 CONTINUE
      Z0=Z(K0)
      DO 5 I=1,JTAU
        Z(I)=Z(I)-Z0
        i1=min(i+1,jtau)
        PTOT(I)=PG(I)+PRAD(I)+0.5*(pturb(i)+pturb(i1))
    5 CONTINUE
***
      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
      DO 22 KTAU=1,NTPO
      DO 20 IE=1,NEL
      NJP=NJ(IE)
      READ(IARCH) KR,TAUI,TI,PEI,IEL(IE),ABUND(IE),
     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
   20 CONTINUE
      DO 21 KLAM=1,NLP
      READ(IARCH) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE
   22 continue
      READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C CONVERT TO 'PHYSICAL' FLUXES
      DO 24 J=1,NLB
   24 FLUXME(J)=3.14159*FLUXME(J)

      dluminosity=0.
      do 25 j=1,nlb
       dluminosity=dluminosity+fluxme(j)*w(j)
25    continue
      dluminosity=dluminosity*4.*3.14159*radius**2/3.82d33
      ddddd=real(dluminosity)
*     print*,'luminosity: ',dluminosity*3.82d33,' erg/s  = ',ddddd,
*    &  ' solar luminosities'
*** Read also continuum flux if MARCS version 10.8 or later
      READ(IARCH,end=999) (FLUXMEcont(J),J=1,NLB)
C CONVERT TO 'PHYSICAL' FLUXES
      DO 26 J=1,NLB
   26 FLUXMEcont(J)=3.14159*FLUXMEcont(J)
      read(iarch,end=999) inputabund
      RETURN
  999 continue
*     write(*,345)
* 345 format(" OLD marcs model, no continuous fluxes or allabundances saved")
***
      RETURN
         END
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
*
      subroutine rename(file,modname,filepub,metallicity,alpha,cc,
     &                  comp,dirname)
* create model name
      implicit none
      character file*256,modname*24,string2*17,string3*13,filepub*72
      character*6 z,a,c,n,o,r,s,zz
      character dirname*42,comp*2
      character g*5,pp*3,xi*3,gsign*1,m*4,xite*3,metallicity*5,alpha*5
      integer ixi,it,ixite
      real rg,cc
*     do i=1,20000
        read(file,'(51x,a24)') modname
        read(file,1010) z,a,c,n,o,r,s,pp,xi,it,string2
 1010   format(7a6,x,a1,3x,a3,x,i4,x,a17)
        read(c,'(x,f5.2)') cc
        dirname=file(1:42)
        read(z,'(x,a5)') metallicity
        read(a,'(x,a5)') alpha
        read(xi,'(x,i2)') ixi
        read(string2,'(a1)') gsign
        if(gsign.eq.'-') then
          read(string2,1020) rg,string3
 1020     format(f4.1,a13)
          write(g,1025) rg
 1025     format('g',f4.1)
        else
          read(string2,1030) rg,string3
 1030     format(f3.1,a13)
          write(g,1035) rg
 1035     format('g+',f3.1)
        endif
        read(string3,1040) zz
 1040   format(a6)
        if(zz.ne.z) then
          print *,file
          stop 'metallicities do not correspond'
        endif
        if(pp.eq.'s') then
          read(string3,1050) m,xite
 1050     format(6x,a4,a3)
        else if(pp.eq.'p') then
          read(string3,1060) xite
 1060     format(6x,a3)
          m='m0.0'
        else
          print *,file
          stop 'pp not ppl or sph!'
        endif
        read(xite,'(x,i2)') ixite
        if(ixite.ne.ixi) then
          print *,file
          print *,ixite,ixi,xite,pp
          stop 'microturb values do not correspond'
        endif
* produce file name:
        write(filepub,2010) pp,it,g,m,xi,comp,z,a,c,n,o,r,s
 2010   format(a1,i4,'_',a5,'_',a4,'_',a3,'_',a2,7('_',a6))
*     enddo
      return
      end
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/sph.t10/2500g-3.0z+0.00m0.5t10
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/ppl.t00/2500g3.5z+0.00t0
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/ppl.t00/2500g4.0z+0.00t0
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/ppl.t00/2500g4.5z+0.00t0
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/ppl.t00/2500g5.0z+0.00t0
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/ppl.t00/2600g3.0z+0.00t0
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/ppl.t00/2600g3.5z+0.00t0
*z+0.00a+0.00c+0.00n+0.00o+0.00r+0.00s+0.00/ppl.t00/2600g4.0z+0.00t0
**...v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
*s2500_g-3.0_m0.5_t10_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
*p2500_g+3.5_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
*p2500_g+4.0_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
*p2500_g+4.5_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
*p2500_g+5.0_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
*p2600_g+3.0_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
*p2600_g+3.5_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
*p2600_g+4.0_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00
      subroutine namemod(dofile,allmodname,comp)
* create model name
      implicit none
      real g,m,feh,alpha,carbon,oxygen
      integer it,iturb
      character P*1,cturb*2,z*5,a*5,c*5,n*5,o*5,r*5,s*5,comp*2,mm*3
      character allmodname*73,dofile*69
*     print *,'subr namemod, dofile=',dofile,'|'
      read(dofile,4010) P,it,g,m,cturb,z,a,c,n,o,r,s
 4010 format(a1,i4,2x,f4.1,2x,f3.1,2x,a2,7(2x,a5))
* determine and label by composition class ******************************
      read(z,'(f5.2)') feh                                              *
      read(a,'(f5.2)') alpha                                            *
      read(c,'(f5.2)') carbon                                           *
      read(o,'(f5.2)') oxygen                                           *
      if((feh.ge.0.00  .and. alpha.eq.0.0 .and. carbon.eq.0.0) .or.     *
     &   (feh.eq.-0.25 .and. alpha.eq.0.1 .and. carbon.eq.0.0) .or.     *
     &   (feh.eq.-0.50 .and. alpha.eq.0.2 .and. carbon.eq.0.0) .or.     *
     &   (feh.eq.-0.75 .and. alpha.eq.0.3 .and. carbon.eq.0.0) .or.     *
     &   (feh.le.-1.00 .and. alpha.eq.0.4 .and. carbon.eq.0.0)) then    *
* this is a standard composition model                                  *
        comp='st'                                                       *
      else if(feh.lt.0.0 .and. alpha.eq.0.0) then                       *
* this is an alpha-poor composition model                               *
        comp='ap'                                                       *
      else if(carbon.eq.-0.13) then                                     *
* this is a moderately CN-cycled composition model                      *
        comp='mc'                                                       *
      else if(carbon.eq.-0.38) then                                     *
* this is a heavily CN-cycled composition model                         *
        comp='hc'                                                       *
      else if(carbon.eq.0.08) then                                      *
* this is a Grevesse&Sauval '98 solar metallicity model:                *
* [Fe/H]=+0.05 [alpha/Fe]=+0.11 [C/Fe]=+0.08 [N/Fe]=+0.09 [O/Fe]=+0.12  *
        comp='gs'                                                       *
      else if(feh.ge.-0.75 .and. alpha.eq.0.4 .and. carbon.eq.0.0) then *
* this is an alpha-enhanced model                                       *
        comp='ae'                                                       *
      else if(alpha.eq.-0.4 .and. carbon.eq.0.0) then                   *
* this is an negative alpha model                                       *
        comp='an'                                                       *
      else if(alpha.eq.+0.6.and.carbon.eq.0.0.and.oxygen.eq.0.6) then   *
* this is a alpha +0.6 model                                            *
        comp='a6'                                                       *
      else if(alpha.eq.+0.8.and.carbon.eq.0.0.and.oxygen.eq.0.8) then   *
* this is a alpha +0.6 model                                            *
        comp='a8'                                                       *
      else                                                              *
        print *,'[Fe/H], [alpha/Fe], [C/Fe]:',feh,alpha,carbon          *
        stop 'unnamed composition?'                                     *
      endif                                                             *
*************************************************************************
      read(cturb,'(i2)') iturb
      if(P.eq.'p') then
        if(g.lt.0.0) then
          if(iturb.ge.10) then
            write(allmodname,1040) z,a,c,n,o,r,s,cturb,it,g,z,iturb
 1040       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/ppl.t',a2,'/',i4,'g',f4.1,'z',a5,'t',i2)
          else
* iturb<10
            write(allmodname,1050) z,a,c,n,o,r,s,cturb,it,g,z,iturb
 1050       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/ppl.t',a2,'/',i4,'g',f4.1,'z',a5,'t',i1)
          endif
        else
* g>=0.0
          if(iturb.ge.10) then
            write(allmodname,1060) z,a,c,n,o,r,s,cturb,it,g,z,iturb
 1060       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/ppl.t',a2,'/',i4,'g',f3.1,'z',a5,'t',i2)
          else
* iturb<10
            write(allmodname,1070) z,a,c,n,o,r,s,cturb,it,g,z,iturb
 1070       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/ppl.t',a2,'/',i4,'g',f3.1,'z',a5,'t',i1)
          endif
        endif
      else
* P='s'
        if(m.ge.0.0 .and. m.lt.9.95) then
          write(mm,'(f3.1)') m
        else if(m.le.99.1) then
          write(mm,'(f3.0)') m
        else
          print *,'mass=',m
          stop 'subr. namemod: mass out of range'
        endif
        if(g.lt.0.0) then
          if(iturb.ge.10) then
            write(allmodname,2040) z,a,c,n,o,r,s,cturb,it,g,z,mm,iturb
 2040       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/sph.t',a2,'/',i4,'g',f4.1,'z',a5,'m',a3,'t',i2)
          else
* iturb<10
            write(allmodname,2050) z,a,c,n,o,r,s,cturb,it,g,z,mm,iturb
 2050       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/sph.t',a2,'/',i4,'g',f4.1,'z',a5,'m',a3,'t',i1)
          endif
        else
* g>=0.0
          if(iturb.ge.10) then
            write(allmodname,2060) z,a,c,n,o,r,s,cturb,it,g,z,mm,iturb
 2060       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/sph.t',a2,'/',i4,'g',f3.1,'z',a5,'m',a3,'t',i2)
          else
* iturb<10
            write(allmodname,2070) z,a,c,n,o,r,s,cturb,it,g,z,mm,iturb
 2070       format('z',a5,'a',a5,'c',a5,'n',a5,'o',a5,'r',a5,'s',a5,
     &             '/sph.t',a2,'/',i4,'g',f3.1,'z',a5,'m',a3,'t',i1)
          endif
        endif
      endif
      return
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine marcs35_2krz(file,file2,sph,xite,inputabund,nstore)
* convert MARCS35 to "Kurucz" format
      parameter(NDPSIZ=100)
*     CHARACTER*117 ADUM
      CHARACTER FILE*73,FILE2*86
      LOGICAL OLD_FORMAT,sph
      dimension tau(NDPSIZ),t(NDPSIZ),z(NDPSIZ),ptot(NDPSIZ),
     & prad(NDPSIZ),pg(NDPSIZ),pturb(NDPSIZ),pe(NDPSIZ),ro(NDPSIZ),
*    & xmass(NDPSIZ),
     & xkapr(NDPSIZ)
      dimension rr(NDPSIZ),taus(NDPSIZ),xlr(30)
      real xite,inputabund(92)
      integer nstore
*     dimension gradp(NDPSIZ),gravity(NDPSIZ),pcheck(NDPSIZ),
*    & geff(NDPSIZ),gradptur(NDPSIZ),dp(NDPSIZ),coldens(NDPSIZ)
      dimension xlb(155000),w(155000),fluxme(155000)
      common /struct/ teff,tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /spectrkrz/ nlb,xlb,w,fluxme
      common /pressurekrz/ presmo(30,NDPSIZ),ptio(NDPSIZ)
      common /binstruc/bPPR(NDPSIZ),bPPT(NDPSIZ),bPP(NDPSIZ),
     & bGG(NDPSIZ),bZZ(NDPSIZ),bDD(NDPSIZ),
     & bVV(NDPSIZ),bFFC(NDPSIZ),bPPE(NDPSIZ),bTT(NDPSIZ),
     & bTAULN(NDPSIZ),NbTAU,IbTER,erad(NDPSIZ)
      common /ABUNDANCE/ABUND(16),abSc,abTi,abV,abMn,abCo,NEL
      dimension JEL(16)
      dimension dabund(99),abunds(99)
     
      data dabund/
c         H      He     Li     Be     B      C      N      O
     &  0.911, 0.089,-10.88,-10.89, -9.44, -3.48, -3.99, -3.11,
c         F      Ne     Na     Mg     Al     Si     P      S
     &  -7.48, -3.95, -5.71, -4.46, -5.57, -4.49, -6.59, -4.83,
c         Cl     Ar     K      Ca     Sc     Ti     V      Cr
     &  -6.54, -5.48, -6.82, -5.68, -8.94, -7.05, -8.04, -6.37,
c         Mn     Fe     Co     Ni     Cu     Zn     Ga     Ge
     &  -6.65, -4.37, -7.12, -5.79, -7.83, -7.44, -9.16, -8.63,
c         As     Se     Br     Kr     Rb     Sr     Y      Zr
     &  -9.67, -8.69, -9.41, -8.81, -9.44, -9.14, -9.80, -9.54,
c         Nb     Mo     Tc     Ru     Rh     Pd     Ag     Cd
     & -10.62,-10.12,-20.00,-10.20,-10.92,-10.35,-11.10,-10.18,
c         In     Sn     Sb     Te     I      Xe     Cs     Ba
     & -10.58,-10.04,-11.04, -9.80,-10.53, -9.81,-10.92, -9.91,
c         La     Ce     Pr     Nd     Pm     Sm     Eu     Gd
     & -10.82,-10.49,-11.33,-10.54,-20.00,-11.04,-11.53,-10.92,
c         Tb     Dy     Ho     Er     Tm     Yb     Lu     Hf
     & -11.94,-10.94,-11.78,-11.11,-12.04,-10.96,-11.28,-11.16,
c         Ta     W      Re     Os     Ir     Pt     Au     Hg
     & -11.91,-10.93,-11.77,-10.59,-10.69,-10.24,-11.03,-10.95,
c         Tl     Pb     Bi     Po     At     Rn     Fr     Ra
     & -11.14,-10.19,-11.33,-20.00,-20.00,-20.00,-20.00,-20.00,
c         Ac     Th     Pa     U      Np     Pu     Am     Cm
     & -20.00,-11.92,-20.00,-12.51,-20.00,-20.00,-20.00,-20.00,
c         Bk     Cf     Es
     & -20.00,-20.00,-20.00/
      DATA JEL/1,2,6,7,8,10,11,12,13,14,16,19,20,24,26,28/
*     EXP10(X0)=EXP(X0*2.30258509299405D0)
C
*     CALL GETARG(1,FILE)
*     IF(FILE.EQ.' ') THEN
*       WRITE(*,'('' Enter MARCS35 filename ? ... '')')
*       READ(*,'(A)',END=1,ERR=1) FILE
*     END IF
C
C  Open input file and read the model
C
      OPEN(UNIT=10,FILE=FILE,STATUS='OLD',FORM='UNFORMATTED',
     &     RECL=152600,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        WRITE(*,*) 'marcs35_2krz can''t open MARCS35 data file'
        STOP
      END IF
      CALL readmarcs(10,NDEPTH,GG,RADIUS,OLD_FORMAT)
C
C  OPEN OUTPUT FILE
C
        OPEN(21,FILE=FILE2,STATUS='UNKNOWN',FORM='FORMATTED')
C
      write(21,2000) xite
2000  format('TITLE MARCS 35 model for MLT convection',
     &       ' testing. Vturb ',f4.1,' km/s')
      if(sph) then
        write(21,2001) TEFF,LOG10(GG),radius
2001    format('T EFF=',F6.0,' GRAV=',F4.1,
     & '  MODEL TYPE= 3 WLSTD= 5000 SPHERICAL, RADIUS=',1pe10.3,' cm')
      else
        write(21,2004) TEFF,LOG10(GG)
2004    format('T EFF=',F6.0,' GRAV=',F4.1,
     &         '  MODEL TYPE= 0 WLSTD= 5000 PLANE-PARALLEL.')
      endif
      write(21,2002)
2002  format(' 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0',
     &       ' - OPACITY SWITCHES')
* change abundance scale to ".krz" format
      abusum=0.0
      do ie=1,92
        if(inputabund(ie).gt.-19.9) then
          abusum=abusum+10.**inputabund(ie)
        endif
      enddo
      do ie=1,2
        abunds(ie)=10.**inputabund(ie)/abusum
      enddo
      do ie=3,92
        if(inputabund(ie).gt.-19.9) then
          abunds(ie)=inputabund(ie)-alog10(abusum)
        else
          abunds(ie)=-20.0
        endif
      enddo
      do ie=93,99
        abunds(ie)=-20.0
      enddo
      write (21,'(2f7.3,8f7.2/(10F7.2))') (ABUNDS(IE),IE=1,90)
ccc* store only 56 depth points
ccc      nstore=56
      write (21,'(9F7.2,I3)') (ABUNDS(IE),IE=91,99),nstore
      if(ndepth.lt.nstore) stop 'subr marcs35_2krz: ndepth < 56 !!!'
      do k=1,nstore
******************************modified to correspond to MARCS defs
*       k1=max(1,k-1)
*       if (k.eq.1) then
*         rhox=(pg(1)+prad(1))/GG
*       else
*         rhox=rhox+(z(k)-z(k1))*(ro(k)+ro(k1))*0.5
*       endif
********KE found this error:**corrected to column mass 070531 BE
        rhox=ptot(k)/gg
        if(sph) then
******************************modified to correspond to MARCS defs
*         write(21,2005) rhox,t(k),pe(k)/1.380622e-16/t(k),
*    &                   pg(k)/1.380622e-16/t(k),ro(k),rr(k)-radius
*******************************************************************
          write(21,2005) rhox,t(k),pe(k)/1.380622e-16/t(k),
     &                   (pg(k)-pe(k))/1.380622e-16/t(k),ro(k),0.-z(k)
* !! Pg in krz is gas pressure minus electron pressure !!
2005      format(1PE16.9,',',0PF9.1,',',4(1PE12.5,','))
        else
          write(21,2003) rhox,t(k),pe(k)/1.380622e-16/t(k),
     &                   (pg(k)-pe(k))/1.380622e-16/t(k),ro(k)
2003      format(1PE16.9,',',0PF9.1,',',3(1PE12.5,','))
        endif
      enddo
      close(21)
*
*     open(22,FILE='flux',STATUS='UNKNOWN',FORM='UNFORMATTED')
*     write(22) NLB
*     write(22) (XLB(J),J=1,NLB)
*     write(22) (FLUXME(J),J=1,NLB)
*     close(22)
*
*  1  STOP
      return
      END
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      SUBROUTINE readmarcs(IARCH,JTAU,G,RADIUS,OLD_FORMAT)
C        THIS ROUTINE READS ONE MODEL, TO GET INFO ON PRAD
C             ( All features taken from listmo )
      PARAMETER (NDPSIZ=100)
C
      DIMENSION ABUND(16),TKORRM(NDPSIZ),FCORR(NDPSIZ),TAU(NDPSIZ),
     * TAUS(NDPSIZ),T(NDPSIZ),PE(NDPSIZ),PG(NDPSIZ),PRAD(NDPSIZ),
     * PTURB(NDPSIZ),XKAPR(NDPSIZ),RO(NDPSIZ),CP(NDPSIZ),CV(NDPSIZ),
     * AGRAD(NDPSIZ),Q(NDPSIZ),U(NDPSIZ),V(NDPSIZ),ANCONV(NDPSIZ),
     * PRESMO(30,NDPSIZ),RR(NDPSIZ),Z(NDPSIZ),
     * EMU(NDPSIZ),HNIC(NDPSIZ),NJ(16),XLR(30),IEL(16),ANJON(16,5),
     * PART(16,5),PROV(50,20+1),ABSKA(50),SPRIDA(50),XLB(155000),
     * FLUXME(155000),ABNAME(50),SOURCE(50),PTOT(NDPSIZ)
*      DIMENSION FCONV(NDPSIZ),FLUMAG(155000),PEP(16)
      DIMENSION W(155000),UW(12),BW(21),VW(25)
      CHARACTER*10 DAG,NAME,NAMEP,KLOCK
      CHARACTER*8 ABNAME,SOURCE
*     DIMENSION WAVFLX(10)
      LOGICAL OLD_FORMAT
      dimension PTIO(NDPSIZ)
      real*8 dluminosity
      real abSc,abTi,abV,abMn,abCo
      common /binstruc/ dummy(11*NDPSIZ+2),erad(NDPSIZ)
      common /struct/ teff,tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /spectrkrz/ nlb,xlb,w,fluxme
      common /pressurekrz/ presmo,ptio
      common /ABUNDANCE/ABUND,abSc,abTi,abV,abMn,abCo,NEL
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/,BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/,VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      DATA NAME/'LOCAL'/,NAMEP/'PARSONS'/
      DATA A,B/.34785485,.65214515/
      IREAD=5
      NDP=100
      OLD_FORMAT=.FALSE.
C
  1   REWIND IARCH
      IF(OLD_FORMAT) THEN
        READ(IARCH,ERR=26,END=27) (DUMMY(I),I=1,NDP*11+2)
        READ(IARCH,ERR=26,END=27) INORD,DAG,KLOCK
        READ(IARCH,ERR=26,END=27) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,
     &              ILINE,ISTRAL,MIHAL,
     &              IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              ITMAX,NEL,(ABUND(I),I=1,NEL)
        READ(IARCH,ERR=26,END=27) RADIUS
        RADIUS=1.
        if(RADIUS.le.10.) RADIUS=6.9599e10
      ELSE
        READ(IARCH,ERR=26,END=27) DUMMY,erad
        READ(IARCH,ERR=26,END=27) INORD,DAG,KLOCK
        READ(IARCH,ERR=26,END=27) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,
     &              ILINE,ISTRAL,MIHAL,
     &              IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              ITMAX,NEL,(ABUND(I),I=1,NEL),abSc,abTi,abV,abMn,abCo
        READ(IARCH,ERR=26,END=27) JTAU,NCORE,DIFLOG,TAUM,
     &                            RADIUS,(RR(K),K=1,JTAU)
*       if(RADIUS.le.10.) RADIUS=6.9599e10
      ENDIF      
      GLOG=ALOG10(G)
      FNORD=0.1*INORD
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX
c      DO 2 I=1,NEL
c    2 ABUND(I)=ALOG10(ABUND(I))+12.
c      abSc=alog10(abSc)+12.
c      abTi=alog10(abTi)+12.
c      abV =alog10(abV) +12.
c      abMn=alog10(abMn)+12.
c      abCo=alog10(abCo)+12.
c
      READ(IARCH,ERR=26,END=27) JTAU,(TKORRM(I),I=1,JTAU),
     &                          (FCORR(K),K=1,JTAU)
      NTPO=0
      IF(OLD_FORMAT) THEN
        NMOL=13
        DO 32 K=1,JTAU
          READ(IARCH,ERR=26,END=27) KR,TAU(K),TAUS(K),Z(K),
     &               T(K),PE(K),PG(K),PRAD(K),
     &               PTURB(K),XKAPR(K)
          READ(IARCH,ERR=26,END=27) RO(K),EMU(K),CP(K),CV(K),
     &               AGRAD(K),Q(K),U(K),V(K),ANCONV(K)
          READ(IARCH,ERR=26,END=27) XMOLH
          READ(IARCH,ERR=26,END=27) HNIC(K),(PRESMO(J,K),J=1,NMOL)
          TAUK=ALOG10(TAU(K))+10.01
          KTAU=TAUK
          IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
          IF(KTAU.EQ.10) K0=K
          NTPO=NTPO+1
c          write(*,'(I3,6E12.5)') 
c     &         K,TAU(K),TAUS(K),T(K),PE(K),PG(K),ptio(k)
   31     CONTINUE
   32   CONTINUE
      ELSE
        DO 34 K=1,JTAU
          READ(IARCH,ERR=26,END=27) KR,TAU(K),TAUS(K),Z(K),
     &                T(K),PE(K),PG(K),PRAD(K),
     &                PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &                AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &                (PRESMO(J,K),J=1,NMOL),ptio(k)
          TAUK=ALOG10(TAU(K))+10.01
          KTAU=TAUK
          IF(ABS(TAUK-KTAU).GT.0.02) GO TO 33
          IF(KTAU.EQ.10) K0=K
          NTPO=NTPO+1
c          write(*,'(I3,6E12.5)') 
c     &         K,TAU(K),TAUS(K),T(K),PE(K),PG(K),ptio(k)
   33     CONTINUE
   34   CONTINUE
      ENDIF
      Z0=Z(K0)
      DO 5 I=1,JTAU
        Z(I)=Z(I)-Z0
        i1=min(i+1,jtau)
        PTOT(I)=PG(I)+PRAD(I)+0.5*(pturb(i)+pturb(i1))
    5 CONTINUE
***
      IF(OLD_FORMAT) THEN
        READ(IARCH,ERR=26,END=27) (NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
        NPROV=28
      ELSE
        READ(IARCH,ERR=26,END=27) (NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     &   ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
      ENDIF
      DO 22 KTAU=1,NTPO
      DO 20 IE=1,NEL
      NJP=NJ(IE)
      READ(IARCH,ERR=26,END=27) KR,TAUI,TI,PEI,IEL(IE),ABUND(IE),
     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
   20 CONTINUE
      DO 21 KLAM=1,NLP
      READ(IARCH,ERR=26,END=27) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE
   22 continue
      READ(IARCH,ERR=26,END=27) NLB,(XLB(J),FLUXME(J),J=1,NLB),
     &                          (W(J),J=1,NLB)
C CONVERT TO 'PHYSICAL' FLUXES
      DO 24 J=1,NLB
   24 FLUXME(J)=3.14159*FLUXME(J)

      dluminosity=0.
      do 25 j=1,nlb
       dluminosity=dluminosity+fluxme(j)*w(j)
25    continue
      dluminosity=dluminosity*4.*3.14159*radius**2/3.82d33
      ddddd=real(dluminosity)
*     write(*,200) dluminosity*3.82d33,ddddd
*200  format(' Luminosity: ',1PE12.6,' erg/s  = ',0PF7.3,
*    &       ' solar luminosities')
*     write(*,*) 'Radius is',radius/6.9599d10,' solar radii'
***
      RETURN
C
C IO ERROR
C
26    IF(OLD_FORMAT) THEN
        WRITE(*,*) 'ERROR reading input file'
        STOP
      ELSE
        OLD_FORMAT=.TRUE.
        NDP=40
        GO TO 1
      ENDIF
C
C IO ERROR
C
27    IF(OLD_FORMAT) THEN
        WRITE(*,*) 'EOF reading input file'
        STOP
      ELSE
        OLD_FORMAT=.TRUE.
        NDP=40
        GO TO 1
      ENDIF

      END
************************************************************************/*******
*     program testxyz
* compute mass fraction of H, He and metals
*     implicit none
*     real abund(92),x,y,z
*     data abund /
*    & 12.00, 10.93,  1.10,  1.40,  2.55,  8.40,  7.80,  8.65,  4.56,   | Photosphere Asplund
*    &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.5,   6.40,
*    &  5.12,  6.36,  3.17,  5.02,  4.00,  5.67,  5.39,  7.50,  4.92,
*    &  6.25,  4.21,  4.60,  2.88,  3.41,  2.37,  3.41,  2.63,  3.21,
*    &  2.60,  2.97,  2.24,  2.60,  1.42,  1.92, -9.99,  1.84,  1.12,
*    &  1.69,  0.94,  1.77,  1.66,  2.0,   1.0,   2.24,  1.51,  2.27,
*    &  1.13,  2.13,  1.17,  1.58,  0.71,  1.50, -9.99,  1.01,  0.51,
*    &  1.12, -0.1,   1.14,  0.26,  0.93,  0.00,  1.08,  0.06,  0.88,
*    & -0.13,  1.11,  0.28,  1.45,  1.35,  1.8,   1.01,  1.13,  0.9,
*    &  1.95,  0.71, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99,  0.09,
*    & -9.99, -0.50 /
*     call sxyz(92,abund,x,y,z)
*     print *,x,y,z
*     print *,'this should have been about: 0.73735, 0.24923, 0.01342'
*     end
*************************************************************************
      subroutine sxyz(n,abund4,xx,yy,zz)
* compute number fraction and mass fraction of elements
* program solarabund.f does this too
      implicit none
      integer i,n
      character*2 lel(92)
      real abund4(92),xx,yy,zz
      real*8 abund(92),weight(92),sumass,weighthamu
      real*8 w,rf,fract,fmass,weightamu,abumetals,wmetals,x,y,z,sumab
*
      data weight /
     &1.0079,4.0026, 6.941, 9.012, 10.81, 12.01, 14.01, 15.99, 19.00,
     & 20.18, 22.99, 24.30, 26.98, 28.09, 30.97, 32.06, 35.45, 39.95,
     & 39.10, 40.08, 44.96, 47.90, 50.94, 52.00, 54.94, 55.85, 58.93,
     & 58.71, 63.55, 65.37, 69.72, 72.59, 74.92, 78.96, 79.90, 83.80,
     & 85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 98.91,101.07,102.91,
     &106.4 ,107.87,112.40,114.82,118.69,121.75,127.60,126.90,131.30,
     &132.91,137.34,138.91,140.12,140.91,144.24,146.  ,150.4 ,151.96,
     &157.25,158.93,162.50,164.93,167.26,168.93,170.04,174.97,178.49,
     &180.95,183.85,186.2 ,190.2 ,192.2 ,195.09,196.97,200.59,204.37,
     &207.19,208.98,210.  ,210.  ,222.  ,223.  ,226.03,227.  ,232.04,
     &230.04,238.03 /
      data lel /
     & 'H ' , 'He' , 'Li' , 'Be' , 'B ' , 'C ' , 'N ' , 'O ' , 'F ' , 
     & 'Ne' , 'Na' , 'Mg' , 'Al' , 'Si' , 'P ' , 'S ' , 'Cl' , 'Ar' , 
     & 'K ' , 'Ca' , 'Sc' , 'Ti' , 'V ' , 'Cr' , 'Mn' , 'Fe' , 'Co' , 
     & 'Ni' , 'Cu' , 'Zn' , 'Ga' , 'Ge' , 'As' , 'Se' , 'Br' , 'Kr' , 
     & 'Rb' , 'Sr' , 'Y ' , 'Zr' , 'Nb' , 'Mo' , 'Tc' , 'Ru' , 'Rh' , 
     & 'Pd' , 'Ag' , 'Cd' , 'In' , 'Sn' , 'Sb' , 'Te' , 'I ' , 'Xe' , 
     & 'Cs' , 'Ba' , 'La' , 'Ce' , 'Pr' , 'Nd' , 'Pm' , 'Sm' , 'Eu' , 
     & 'Gd' , 'Tb' , 'Dy' , 'Ho' , 'Er' , 'Tm' , 'Yb' , 'Lu' , 'Hf' , 
     & 'Ta' , 'W ' , 'Re' , 'Os' , 'Ir' , 'Pt' , 'Au' , 'Hg' , 'Tl' , 
     & 'Pb' , 'Bi' , 'Po' , 'At' , 'Rn' , 'Fr' , 'Ra' , 'Ac' , 'Th' , 
     & 'Pa' , 'U ' /
* sum all masses * abundances
      if(n.ne.92) stop 'subr SXYZ not 92 abundance numbers'
      w=0.0d0
      rf=0.0d0
      do i=92,3,-1
        abund(i)=dble(abund4(i))
        if(abund(i).gt.-9.) then
          rf=rf + 10.d0**abund(i)
          w=w + weight(i)*10.d0**abund(i)
        endif
      enddo
      abumetals=dlog10(rf)
*     print *,'abumetals=',abumetals
      wmetals=w
      abund(2)=dble(abund4(2))
      rf=rf + 10.d0**abund(2)
      w=w + weight(2)*10.d0**abund(2)
      abund(1)=dble(abund4(1))
      rf=rf + 10.d0**abund(1)
      w=w + weight(1)*10.d0**abund(1)
*     print *,'Z=',wmetals/w
*     print *,'Y=',weight(2)*10.d0**abund(2)/w
*     print *,'X=',weight(1)*10.d0**abund(1)/w
*
      weightamu=w/rf
* weightamu is the mean weight of an atom (in amu)
*     print 1000
      sumab=0.0
      sumass=0.0
      do 20 i=1,30
        if(abund(i).gt.-3.) then
          fract=10.**abund(i)/rf
          fmass=10.**abund(i)*weight(i)/w
          if(i.eq.1) then
            X=fmass
            weightHamu=weightamu/fract
          else if(i.eq.2) then
            Y=fmass
          endif
        else
          fract=0.0
          fmass=0.0
        endif
        if(i.gt.2) then
          sumab=sumab+fract
          sumass=sumass+fmass
*         print 1010,i,lel(i),weight(i),abund(i),fract,fmass
*1010     format(1x,i2,3x,a2,f10.4,10x,f5.2,8x,1pe10.2,6x,e10.2)
        else
*         print 1020,i,lel(i),weight(i),abund(i),fract,fmass
*1020     format(1x,i2,3x,a2,f10.4,10x,f5.2,8x,f10.6,6x,f10.6)
        endif
*1000   format('Element Weight(amu)  log abundance   Number fraction',
*    &         ' Mass fraction')
   20 continue
      Z=sumass
*     print *,'======================================================'
*     write(*,1030) abumetals
*1030 format(' log. abundance of metal atoms    =   ',f8.5)
*     write(*,1040) sumab
*1040 format(' metal (Z>2) atom number fraction =   ',f8.5)
*     write(*,1050) sumass
*1050 format(' metal (Z>2) atom mass fraction   =   ',f8.5)
*     print *,
*    &'mean mass of a metal atom        =',sumass/sumab*weightamu,' amu'
*     print *,'mean mass per atom               =',weightamu,' amu'
*     print *,'mean mass per hydrogen atom      =',weightHamu,' amu'
*     write(*,1060) x,y,z,x+y+z
*1060 format(' X=',f8.5,'   Y=',f8.5,'   Z=',f8.5,'   X+Y+Z=',f8.5)
      xx=real(x)
      yy=real(y)
      zz=real(z)
      return
      end
