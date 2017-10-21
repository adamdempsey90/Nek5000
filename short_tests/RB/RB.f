c-----------------------------------------------------------------------
      subroutine RB !THIS SETS BASIC PARAMETERS FOR THE PROBLEM
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux
      external kappafunc
      real invtrapz
      integer n,nintervals
      real omega, lowerlim
      real gam
c rotation parameters below
      omega = param(3)
      gam= param(4) !in degrees from vertical for rotation axis
      gam=(gam/180.0D0)*4.0D0*datan(1.0D0)
      cgam=2.0D0*omega*cos(gam)
      sgam=2.0D0*omega*sin(gam)
      n = nx1*ny1*nz1*nelv
      zmin = -1.0D0
      zmax = 2.0D0
c      zmin = glmin(zm1,n)
c      zmax = glmax(zm1,n)

      nintervals = 1000
      lowerlim = 0.0D0
      bctemp = param(5)
      bcflux = 1.0D0
c      if(istep.eq.0) write(*,*) omega,cgam,sgam,bctemp,bcflux
      return
      end
c-----------------------------------------------------------------------
      subroutine my_full_restart
      include 'SIZE'
      include 'TOTAL'

      character*80 s80(3)

      call blank(s80,3*80)
      s80(1) ='rs6RB0.f00004'
      s80(2) ='rs6RB0.f00005'
      s80(3) ='rs6RB0.f00006'

      if(param(18).ne.0) call full_restart(s80,3)  ! Will overload 4-6 onto steps 0-2


      iosave = int(param(17))           ! Trigger save based on iostep
      call full_restart_save(iosave)

      return
      end
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg) !SETS variable viscosity, density. I've not used this
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux
      integer e,f,eg 
      real N1,N2,N3
      real kappafunc
c variable conductivity called if param(30)>0
      if(ifield.eq.1) then ! MOMENTUM EQUATION 
        udiff=param(2); !viscosity
        utrans=param(1); !density
      elseif(ifield.eq.2) then ! HEAT EQUATION
         ! kappa profile 

         udiff = kappafunc(z,param(9),param(6),
     $              param(79),param(5),-1.0D0)

        utrans=param(7); !density* C_p
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,eg) !RHS of momentum equation -- acceleration terms (not force)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,f,eg
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux
c      e = gllel(eg)
c Coriolis force for rotation about axis inclined to z by angle gam
c also buoyancy term in z direction
      ffx=uy*cgam 
      ffy=-ux*cgam+uz*sgam
      ffz= temp-uy*sgam
      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg) !RHS of heat equation
      include 'SIZE'
      include 'TOTAL'      
      include 'NEKUSE'
      integer e,f,eg
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux
      real fac1,fac2, N1, N2, N3, totfac,fac3 

      qvol = -uz*param(5)   !for uz * del_ad

c check the profile
C      if(nid .eq. 0) then !root process
C         open(52,file="qvol.dat",access="append")
C         write(52,*) z,qvol 
C      endif
C      close(52)
      return
      end
c-----------------------------------------------------------------------
      subroutine userchk ! called once per step. I output history file data here
      include 'SIZE'
      include 'TOTAL'

      integer histstep, n, ifla, zavgdumpstep,indfortype
      integer volstep

      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux

      parameter(lt=lx1*ly1*lz1*lelv)  
      parameter(nelz=8) ! Number of elements in z direction 
      common /scrns/ vort(lt,3),w1(lt),w2(lt)

      zavgdumpstep = int(param(80))
      volstep = int(param(81)) !output history points every histstep timesteps
      histstep = int(param(52)) !output history points every histstep timesteps

      call my_full_restart



      if(istep.eq.0) then
        call RB !set up parameters for the problem
      endif



      if(mod(istep,zavgdumpstep). eq.0) then

        call zavg_dump
     
       endif

      if(mod(istep,volstep) .eq. 0) then

          call volavg_dump(0.0D0, 1.0D0)

      endif
      if(mod(istep,histstep) .eq. 0) then

          call hpts()

      endif


c output vorticity into separate .fld files
      if (istep.gt.0.and.mod(istep,iostep).eq.0) then
        call output_vort;
      endif


c output history file data of mean quantities every histstep timesteps

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg) !Boundary conditions
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux
        if(iside.eq.5) then 
            temp = bctemp
        else
            flux = -bcflux
        endif
      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg) !initial conditions
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux
      integer nintervals
      real zcut 
      real scut,scut2,aval,dval
      real z1,z2

      zcut = param(77)
      dval = -1.0D0*param(5)
      aval = param(79)
      nval = -1.0D0*param(75)


      scut = dval*(zcut-1.0D0-log((1+aval*zcut)/(1+aval))/aval)
      scut2 = scut + nval*(1.0D0+2.0D0*zcut)

       z1 = -1.0D0*zcut
       z2 = 1.0D0 + zcut
       if (z.lt.z1) then
          temp = dval*(-1.0D0*(z+1.0D0)
     $         -log((1.0D0-aval*z)/(1+aval))/aval)
        else
           if(z.gt.z2) then
             temp = scut2 + 
     $       dval*(-1.0D0*(z-z2) +
     $          log((1+aval*(z-1.0D0))/(1+aval*zcut))/aval)
           else
               temp = scut + nval*(z + zcut)
           endif
       endif

        

       temp =  temp - param(5)*z

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat ! can modify mesh here
      include 'SIZE'
      include 'TOTAL'
      integer e,f
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2 ! can deform mesh here
      include 'SIZE'
      include 'TOTAL'
      integer e,f
      param(66) = 6 !output file type
      param(67) = 6 !input file type
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3 !not sure what this does...
      include 'SIZE'
      include 'TOTAL'
      return
      end
c-----------------------------------------------------------------------
      subroutine volavg_dump(lowb, upb)
        include 'SIZE'
        include 'TOTAL'
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux

      real kappafunc
      real lowb, upb, work
      real vcut

        integer n 
       
         real  dxux(lx1,ly1,lz1,lelv),
     $     dxuz(lx1,ly1,lz1,lelv),dyux(lx1,ly1,lz1,lelv),
     $     dyuz(lx1,ly1,lz1,lelv),dzux(lx1,ly1,lz1,lelv),
     $     dzuz(lx1,ly1,lz1,lelv),dxuy(lx1,ly1,lz1,lelv),
     $     dyuy(lx1,ly1,lz1,lelv),dzuy(lx1,ly1,lz1,lelv),
     $     dxth(lx1,ly1,lz1,lelv),dyth(lx1,ly1,lz1,lelv),
     $     dzth(lx1,ly1,lz1,lelv),fakevec(lx1,ly1,lz1,lelv)

      common /csij/
     $     sij(lx1*ly1*lz1,6,lelv), ! 6 components
     $     ur(lx1*ly1*lz1), us(lx1*ly1*lz1), ut(lx1*ly1*lz1),
     $     vr(lx1*ly1*lz1), vs(lx1*ly1*lz1), vt(lx1*ly1*lz1),
     $     wr(lx1*ly1*lz1), ws(lx1*ly1*lz1), wt(lx1*ly1*lz1)
     
        real regionvol, totalvol
        real uz2val,
     $        up2val,
     $        KEval,
     $        vortx,
     $        vorty,
     $        vortz,
     $        enstval,
     $        visc1val,
     $        visc2val,
     $        visc3val,
     $        kapval,
     $        thermval,
     $        dtdzval,
     $        Fkval,
     $        Fcval,
     $        Feval
        real uz2_avgU,
     $        up2_avgU,
     $        KE_avgU,
     $        enst_avgU,
     $        visc1_avgU,
     $        visc2_avgU,
     $        visc3_avgU,
     $        therm_avgU,
     $        dtdz_avgU,
     $        Fk_avgU,
     $        Fc_avgU,
     $        Fe_avgU,
     $        area_avgU,
     $        uz2_avgD,
     $        up2_avgD,
     $        KE_avgD,
     $        enst_avgD,
     $        visc1_avgD,
     $        visc2_avgD,
     $        visc3_avgD,
     $        therm_avgD,
     $        dtdz_avgD,
     $        Fk_avgD,
     $        Fc_avgD,
     $        Fe_avgD,
     $        area_avgD

        real uz2_avgU_r,
     $        up2_avgU_r,
     $        KE_avgU_r,
     $        enst_avgU_r,
     $        visc1_avgU_r,
     $        visc2_avgU_r,
     $        visc3_avgU_r,
     $        therm_avgU_r,
     $        dtdz_avgU_r,
     $        Fk_avgU_r,
     $        Fc_avgU_r,
     $        Fe_avgU_r,
     $        area_avgU_r,
     $        uz2_avgD_r,
     $        up2_avgD_r,
     $        KE_avgD_r,
     $        enst_avgD_r,
     $        visc1_avgD_r,
     $        visc2_avgD_r,
     $        visc3_avgD_r,
     $        therm_avgD_r,
     $        dtdz_avgD_r,
     $        Fk_avgD_r,
     $        Fc_avgD_r,
     $        Fe_avgD_r,
     $        area_avgD_r


        n = nx1*ny1*nz1*nelv
        vcut = 0.0D0
       call gradm1(dxux,dyux,dzux,vx) ! 1st derivative
       call gradm1(dxuy,dyuy,dzuy,vy)
       call gradm1(dxuz,dyuz,dzuz,vz)
       call gradm1(dxth,dyth,dzth,t)


         uz2_avgU = 0.0D0
         up2_avgU = 0.0D0
         KE_avgU = 0.0D0
         enst_avgU = 0.0D0
         visc1_avgU = 0.0D0
         visc2_avgU = 0.0D0
         visc3_avgU = 0.0D0
         therm_avgU = 0.0D0
         dtdz_avgU = 0.0D0
         Fk_avgU = 0.0D0
         Fc_avgU = 0.0D0
         Fe_avgU = 0.0D0
         area_avgU = 0.0D0

         uz2_avgD = 0.0D0
         up2_avgD = 0.0D0
         KE_avgD = 0.0D0
         enst_avgD = 0.0D0
         visc1_avgD = 0.0D0
         visc2_avgD = 0.0D0
         visc3_avgD = 0.0D0
         therm_avgD = 0.0D0
         dtdz_avgD = 0.0D0
         Fk_avgD = 0.0D0
         Fc_avgD = 0.0D0
         Fe_avgD = 0.0D0
         area_avgD = 0.0D0
         
         uz2_avgU_r = 0.0D0
         up2_avgU_r = 0.0D0
         KE_avgU_r = 0.0D0
         enst_avgU_r = 0.0D0
         visc1_avgU_r = 0.0D0
         visc2_avgU_r = 0.0D0
         visc3_avgU_r = 0.0D0
         therm_avgU_r = 0.0D0
         dtdz_avgU_r = 0.0D0
         Fk_avgU_r = 0.0D0
         Fc_avgU_r = 0.0D0
         Fe_avgU_r = 0.0D0
         area_avgU_r = 0.0D0

         uz2_avgD_r = 0.0D0
         up2_avgD_r = 0.0D0
         KE_avgD_r = 0.0D0
         enst_avgD_r = 0.0D0
         visc1_avgD_r = 0.0D0
         visc2_avgD_r = 0.0D0
         visc3_avgD_r = 0.0D0
         therm_avgD_r = 0.0D0
         dtdz_avgD_r = 0.0D0
         Fk_avgD_r = 0.0D0
         Fc_avgD_r = 0.0D0
         Fe_avgD_r = 0.0D0
         area_avgD_r = 0.0D0
         
        totalvol = 0.0D0
        regionvol = 0.0D0
        work = 0.0D0
    
        do i=1,n
         uz2val = vz(i,1,1,1)**2
         up2val = vx(i,1,1,1)**2 + vy(i,1,1,1)**2
         KEval = 0.5D0*(uz2val + up2val)
         
         vortx = dyuz(i,1,1,1) - dzuy(i,1,1,1)
         vorty = dzux(i,1,1,1) - dxuz(i,1,1,1)
         vortz = dxuy(i,1,1,1) - dyux(i,1,1,1)
        
         enstval = 0.5D0*(vortx**2 + vorty**2 + vortz**2) 

         visc1val = 2.0D0*param(2)*enstval

c        visc2val = 0.5D0*param(2)*(
c     $        sij(1,1,1)(i,1,1,1)**2       
c     $   +    sij(1,2,1)(i,1,1,1)**2       
c     $   +    sij(1,3,1)(i,1,1,1)**2       
c     $   +    2.0D0*sij(1,4,1)(i,1,1,1)**2       
c     $   +    2.0D0*sij(1,5,1)(i,1,1,1)**2       
c     $   +    2.0D0*sij(1,6,1)(i,1,1,1)**2) 
         
         visc3val = param(2)* (
     $    dxux(i,1,1,1)**2 + dyux(i,1,1,1)**2 + dzux(i,1,1,1)**2
     $   +   dxuy(i,1,1,1)**2 + dyuy(i,1,1,1)**2 + dzuy(i,1,1,1)**2
     $   +   dxuz(i,1,1,1)**2 + dyuz(i,1,1,1)**2 + dzuz(i,1,1,1)**2 )

         kapval = kappafunc(zm1(i,1,1,1),param(9),param(6),
     $              param(79),param(5),-1.0D0)
         thermval = -kapval*(
     $      dxth(i,1,1,1)**2 + dyth(i,1,1,1)**2 + dzth(i,1,1,1)**2)
         dtdzval = dzth(i,1,1,1)

         Fkval = -kapval*dzth(i,1,1,1)
         Fcval = vz(i,1,1,1)*t(i,1,1,1,1)
         Feval = KEval*vz(i,1,1,1)

         totalvol = totalvol + bm1(i,1,1,1)
         if(uzval.ge.vcut) then
           uz2_avgU = uz2_avgU + bm1(i,1,1,1)*uz2val
           up2_avgU = up2_avgU + bm1(i,1,1,1)*up2val
           KE_avgU = KE_avgU + bm1(i,1,1,1)*KEval
           enst_avgU = enst_avgU + bm1(i,1,1,1)*enstval
           visc1_avgU = visc1_avgU + bm1(i,1,1,1)*visc1val
           visc2_avgU = visc2_avgU + bm1(i,1,1,1)*visc2val
           visc3_avgU = visc3_avgU + bm1(i,1,1,1)*visc3val
           therm_avgU = therm_avgU + bm1(i,1,1,1)*thermval
           dtdz_avgU = dtdz_avgU + bm1(i,1,1,1)*dtdzval
           Fk_avgU = Fk_avgU + bm1(i,1,1,1)*Fkval
           Fc_avgU = Fc_avgU + bm1(i,1,1,1)*Fcval
           Fe_avgU = Fe_avgU + bm1(i,1,1,1)*Feval
           area_avgU = area_avgU + bm1(i,1,1,1)*bm1(i,1,1,1)

         
           if((zm1(i,1,1,1).ge.lowb).and.(zm1(i,1,1,1).le.upb)) then 
              uz2_avgU_r = uz2_avgU_r + bm1(i,1,1,1)*uz2val
              up2_avgU_r = up2_avgU_r + bm1(i,1,1,1)*up2val
              KE_avgU_r = KE_avgU_r + bm1(i,1,1,1)*KEval
              enst_avgU_r = enst_avgU_r + bm1(i,1,1,1)*enstval
              visc1_avgU_r = visc1_avgU_r + bm1(i,1,1,1)*visc1val
              visc2_avgU_r = visc2_avgU_r + bm1(i,1,1,1)*visc2val
              visc3_avgU_r = visc3_avgU_r + bm1(i,1,1,1)*visc3val
              therm_avgU_r = therm_avgU_r + bm1(i,1,1,1)*thermval
              dtdz_avgU_r = dtdz_avgU_r + bm1(i,1,1,1)*dtdzval
              Fk_avgU_r = Fk_avgU_r + bm1(i,1,1,1)*Fkval
              Fc_avgU_r = Fc_avgU_r + bm1(i,1,1,1)*Fcval
              Fe_avgU_r = Fe_avgU_r + bm1(i,1,1,1)*Feval
              area_avgU_r = area_avgU_r + bm1(i,1,1,1)*bm1(i,1,1,1)
              regionvol = regionvol + bm1(i,1,1,1)
           endif
         else
           uz2_avgD = uz2_avgD + bm1(i,1,1,1)*uz2val
           up2_avgD = up2_avgD + bm1(i,1,1,1)*up2val
           KE_avgD = KE_avgD + bm1(i,1,1,1)*KEval
           enst_avgD = enst_avgD + bm1(i,1,1,1)*enstval
           visc1_avgD = visc1_avgD + bm1(i,1,1,1)*visc1val
           visc2_avgD = visc2_avgD + bm1(i,1,1,1)*visc2val
           visc3_avgD = visc3_avgD + bm1(i,1,1,1)*visc3val
           therm_avgD = therm_avgD + bm1(i,1,1,1)*thermval
           dtdz_avgD = dtdz_avgD + bm1(i,1,1,1)*dtdzval
           Fk_avgD = Fk_avgD + bm1(i,1,1,1)*Fkval
           Fc_avgD = Fc_avgD + bm1(i,1,1,1)*Fcval
           Fe_avgD = Fe_avgD + bm1(i,1,1,1)*Feval
           area_avgD = area_avgD + bm1(i,1,1,1)*bm1(i,1,1,1)

           if((zm1(i,1,1,1).ge.lowb).and.(zm1(i,1,1,1).le.upb)) then 
              uz2_avgD_r = uz2_avgD_r + bm1(i,1,1,1)*uz2val
              up2_avgD_r = up2_avgD_r + bm1(i,1,1,1)*up2val
              KE_avgD_r = KE_avgD_r + bm1(i,1,1,1)*KEval
              enst_avgD_r = enst_avgD_r + bm1(i,1,1,1)*enstval
              visc1_avgD_r = visc1_avgD_r + bm1(i,1,1,1)*visc1val
              visc2_avgD_r = visc2_avgD_r + bm1(i,1,1,1)*visc2val
              visc3_avgD_r = visc3_avgD_r + bm1(i,1,1,1)*visc3val
              therm_avgD_r = therm_avgD_r + bm1(i,1,1,1)*thermval
              dtdz_avgD_r = dtdz_avgD_r + bm1(i,1,1,1)*dtdzval
              Fk_avgD_r = Fk_avgD_r + bm1(i,1,1,1)*Fkval
              Fc_avgD_r = Fc_avgD_r + bm1(i,1,1,1)*Fcval
              Fe_avgD_r = Fe_avgD_r + bm1(i,1,1,1)*Feval
              area_avgD_r = area_avgD_r + bm1(i,1,1,1)*bm1(i,1,1,1)
              regionvol = regionvol + bm1(i,1,1,1)
           endif
         endif


        enddo
        
        call gop(uz2_avgU,work,'+  ',1)
        call gop(up2_avgU,work,'+  ',1)
        call gop(KE_avgU,work,'+  ',1)
        call gop(enst_avgU,work,'+  ',1)
        call gop(visc1_avgU,work,'+  ',1)
        call gop(visc2_avgU,work,'+  ',1)
        call gop(visc3_avgU,work,'+  ',1)
        call gop(therm_avgU,work,'+  ',1)
        call gop(dtdz_avgU,work,'+  ',1)
        call gop(Fk_avgU,work,'+  ',1)
        call gop(Fc_avgU,work,'+  ',1)
        call gop(Fe_avgU,work,'+  ',1)
        call gop(area_avgU,work,'+  ',1)
        call gop(uz2_avgD,work,'+  ',1)
        call gop(up2_avgD,work,'+  ',1)
        call gop(KE_avgD,work,'+  ',1)
        call gop(enst_avgD,work,'+  ',1)
        call gop(visc1_avgD,work,'+  ',1)
        call gop(visc2_avgD,work,'+  ',1)
        call gop(visc3_avgD,work,'+  ',1)
        call gop(therm_avgD,work,'+  ',1)
        call gop(dtdz_avgD,work,'+  ',1)
        call gop(Fk_avgD,work,'+  ',1)
        call gop(Fc_avgD,work,'+  ',1)
        call gop(Fe_avgD,work,'+  ',1)
        call gop(area_avgD,work,'+  ',1)
        call gop(totalvol,work,'+  ',1)

        call gop(uz2_avgU_r,work,'+  ',1)
        call gop(up2_avgU_r,work,'+  ',1)
        call gop(KE_avgU_r,work,'+  ',1)
        call gop(enst_avgU_r,work,'+  ',1)
        call gop(visc1_avgU_r,work,'+  ',1)
        call gop(visc2_avgU_r,work,'+  ',1)
        call gop(visc3_avgU_r,work,'+  ',1)
        call gop(therm_avgU_r,work,'+  ',1)
        call gop(dtdz_avgU_r,work,'+  ',1)
        call gop(Fk_avgU_r,work,'+  ',1)
        call gop(Fc_avgU_r,work,'+  ',1)
        call gop(Fe_avgU_r,work,'+  ',1)
        call gop(area_avgU_r,work,'+  ',1)
        call gop(uz2_avgD_r,work,'+  ',1)
        call gop(up2_avgD_r,work,'+  ',1)
        call gop(KE_avgD_r,work,'+  ',1)
        call gop(enst_avgD_r,work,'+  ',1)
        call gop(visc1_avgD_r,work,'+  ',1)
        call gop(visc2_avgD_r,work,'+  ',1)
        call gop(visc3_avgD_r,work,'+  ',1)
        call gop(therm_avgD_r,work,'+  ',1)
        call gop(dtdz_avgD_r,work,'+  ',1)
        call gop(Fk_avgD_r,work,'+  ',1)
        call gop(Fc_avgD_r,work,'+  ',1)
        call gop(Fe_avgD_r,work,'+  ',1)
        call gop(area_avgD_r,work,'+  ',1)
        call gop(regionvol,work,'+  ',1)


        uz2_avgU = uz2_avgU / totalvol
        up2_avgU = up2_avgU / totalvol
        KE_avgU = KE_avgU / totalvol
        enst_avgU = enst_avgU / totalvol
        visc1_avgU = visc1_avgU / totalvol
        visc2_avgU = visc2_avgU / totalvol
        visc3_avgU = visc3_avgU / totalvol
        therm_avgU = therm_avgU / totalvol
        dtdz_avgU = dtdz_avgU / totalvol
        Fk_avgU = Fk_avgU / totalvol
        Fc_avgU = Fc_avgU / totalvol
        Fe_avgU = Fe_avgU / totalvol
        area_avgU = area_avgU / totalvol

        uz2_avgD = uz2_avgD / totalvol
        up2_avgD = up2_avgD / totalvol
        KE_avgD = KE_avgD / totalvol
        enst_avgD = enst_avgD / totalvol
        visc1_avgD = visc1_avgD / totalvol
        visc2_avgD = visc2_avgD / totalvol
        visc3_avgD = visc3_avgD / totalvol
        therm_avgD = therm_avgD / totalvol
        dtdz_avgD = dtdz_avgD / totalvol
        Fk_avgD = Fk_avgD / totalvol
        Fc_avgD = Fc_avgD / totalvol
        Fe_avgD = Fe_avgD / totalvol
        area_avgD = area_avgD / totalvol

        uz2_avgU_r = uz2_avgU_r / regionvol
        up2_avgU_r = up2_avgU_r / regionvol
        KE_avgU_r = KE_avgU_r / regionvol
        enst_avgU_r = enst_avgU_r / regionvol
        visc1_avgU_r = visc1_avgU_r / regionvol
        visc2_avgU_r = visc2_avgU_r / regionvol
        visc3_avgU_r = visc3_avgU_r / regionvol
        therm_avgU_r = therm_avgU_r / regionvol
        dtdz_avgU_r = dtdz_avgU_r / regionvol
        Fk_avgU_r = Fk_avgU_r / regionvol
        Fc_avgU_r = Fc_avgU_r / regionvol
        Fe_avgU_r = Fe_avgU_r / regionvol
        area_avgU_r = area_avgU_r / regionvol
        uz2_avgD_r = uz2_avgD_r / regionvol
        up2_avgD_r = up2_avgD_r / regionvol
        KE_avgD_r = KE_avgD_r / regionvol
        enst_avgD_r = enst_avgD_r / regionvol
        visc1_avgD_r = visc1_avgD_r / regionvol
        visc2_avgD_r = visc2_avgD_r / regionvol
        visc3_avgD_r = visc3_avgD_r / regionvol
        therm_avgD_r = therm_avgD_r / regionvol
        dtdz_avgD_r = dtdz_avgD_r / regionvol
        Fk_avgD_r = Fk_avgD_r / regionvol
        Fc_avgD_r = Fc_avgD_r / regionvol
        Fe_avgD_r = Fe_avgD_r / regionvol
        area_avgD_r = area_avgD_r / regionvol


        if(nid .eq. 0) then !root process
         open(53,file="volvg_up.dat",access="append")
          write(53,'(1p20E15.7)') 
     $        time,
     $        uz2_avgU,
     $        up2_avgU,
     $        KE_avgU,
     $        enst_avgU,
     $        visc1_avgU,
     $        visc2_avgU,
     $        visc3_avgU,
     $        therm_avgU,
     $        dtdz_avgU,
     $        Fk_avgU,
     $        Fc_avgU,
     $        Fe_avgU,
     $        area_avgU
          close(53)
         open(53,file="volvg_down.dat",access="append")
          write(53,'(1p20E15.7)') 
     $        time,
     $        uz2_avgD,
     $        up2_avgD,
     $        KE_avgD,
     $        enst_avgD,
     $        visc1_avgD,
     $        visc2_avgD,
     $        visc3_avgD,
     $        therm_avgD,
     $        dtdz_avgD,
     $        Fk_avgD,
     $        Fc_avgD,
     $        Fe_avgD,
     $        area_avgD
          close(53)

         open(52,file="volavgR_up.dat",access="append")
          write(52,'(1p20E15.7)') 
     $        time,
     $        uz2_avgU_r,
     $        up2_avgU_r,
     $        KE_avgU_r,
     $        enst_avgU_r,
     $        visc1_avgU_r,
     $        visc2_avgU_r,
     $        visc3_avgU_r,
     $        therm_avgU_r,
     $        dtdz_avgU_r,
     $        Fk_avgU_r,
     $        Fc_avgU_r,
     $        Fe_avgU_r,
     $        area_avgU_r
         close(52)
         open(52,file="volavgR_down.dat",access="append")
          write(52,'(1p20E15.7)') 
     $        time,
     $        uz2_avgD_r,
     $        up2_avgD_r,
     $        KE_avgD_r,
     $        enst_avgD_r,
     $        visc1_avgD_r,
     $        visc2_avgD_r,
     $        visc3_avgD_r,
     $        therm_avgD_r,
     $        dtdz_avgD_r,
     $        Fk_avgD_r,
     $        Fc_avgD_r,
     $        Fe_avgD_r,
     $        area_avgD_r
         close(52)
          
        endif
        return
        end
c-----------------------------------------------------------------------
      subroutine zavg_dump
       include 'SIZE'
       include 'TOTAL'
      integer eg,ex,ey,ez,f,n,e
      real cgam,sgam,zmin,zmax,bctemp,bcflux
      common /mypar/ cgam,sgam,zmin,zmax,bctemp,bcflux
      integer iflag,nelxy
      real vcut
      parameter(lt=lx1*ly1*lz1*lelv)  
      parameter(nelz=8) ! Number of elements in z direction 
      real kappafunc 
      real uxval,
     $     uyval,
     $     uzval,
     $     ux2val,
     $     uy2val,
     $     uz2val,
     $     KEval,
     $     dzTval,
     $     Tval,
     $     Fcval,
     $     Fkval,
     $     Feval,
     $     kapval,
     $     vortval,
     $     vortxval,
     $     vortyval,
     $     visc3val,
     $     thermval,
     $     Pval,
     $     dzPval,
     $     Pfluxval,
     $     TdzPval,
     $     trip1val,
     $     trip2val,
     $     rey12val,
     $     rey13val,
     $     rey23val,
     $     vfluxval,
     $     kefluxval,
     $     tfluxval,
     $     thermfluxval, 
     $     Fcflux1val,
     $     Fcflux2val,
     $     visc1val,
     $     visc2val,
     $     Aval,
     $     zval

      real uxbarU(lz1,nelz),
     $     uybarU(lz1,nelz),
     $     uzbarU(lz1,nelz),
     $     ux2barU(lz1,nelz),
     $     uy2barU(lz1,nelz),
     $     uz2barU(lz1,nelz),
     $     KEbarU(lz1,nelz),
     $     dzTbarU(lz1,nelz),
     $     TbarU(lz1,nelz),
     $     T2barU(lz1,nelz),
     $     FcbarU(lz1,nelz),
     $     FkbarU(lz1,nelz),
     $     FebarU(lz1,nelz),
     $     kappabarU(lz1,nelz),
     $     vortbarU(lz1,nelz),
     $     vortxbarU(lz1,nelz),
     $     vortybarU(lz1,nelz),
     $     visc3barU(lz1,nelz),
     $     thermbarU(lz1,nelz),
     $     PbarU(lz1,nelz),
     $     dzPbarU(lz1,nelz),
     $     PfluxbarU(lz1,nelz),
     $     TdzPbarU(lz1,nelz),
     $     trip1barU(lz1,nelz),
     $     trip2barU(lz1,nelz),
     $     rey12barU(lz1,nelz),
     $     rey13barU(lz1,nelz),
     $     rey23barU(lz1,nelz),
     $     vfluxbarU(lz1,nelz),
     $     kefluxbarU(lz1,nelz),
     $     tfluxbarU(lz1,nelz),
     $     thermfluxbarU(lz1,nelz), 
     $     Fcflux1barU(lz1,nelz),
     $     Fcflux2barU(lz1,nelz),
     $     visc1barU(lz1,nelz),
     $     visc2barU(lz1,nelz),
     $     areabarU(lz1,nelz),
     $     uxbarD(lz1,nelz),
     $     uybarD(lz1,nelz),
     $     uzbarD(lz1,nelz),
     $     ux2barD(lz1,nelz),
     $     uy2barD(lz1,nelz),
     $     uz2barD(lz1,nelz),
     $     KEbarD(lz1,nelz),
     $     dzTbarD(lz1,nelz),
     $     TbarD(lz1,nelz),
     $     T2barD(lz1,nelz),
     $     FcbarD(lz1,nelz),
     $     FkbarD(lz1,nelz),
     $     FebarD(lz1,nelz),
     $     kappabarD(lz1,nelz),
     $     vortbarD(lz1,nelz),
     $     vortxbarD(lz1,nelz),
     $     vortybarD(lz1,nelz),
     $     visc3barD(lz1,nelz),
     $     thermbarD(lz1,nelz),
     $     PbarD(lz1,nelz),
     $     dzPbarD(lz1,nelz),
     $     PfluxbarD(lz1,nelz),
     $     TdzPbarD(lz1,nelz),
     $     trip1barD(lz1,nelz),
     $     trip2barD(lz1,nelz),
     $     rey12barD(lz1,nelz),
     $     rey13barD(lz1,nelz),
     $     rey23barD(lz1,nelz),
     $     vfluxbarD(lz1,nelz),
     $     kefluxbarD(lz1,nelz),
     $     tfluxbarD(lz1,nelz),
     $     thermfluxbarD(lz1,nelz), 
     $     Fcflux1barD(lz1,nelz),
     $     Fcflux2barD(lz1,nelz),
     $     visc1barD(lz1,nelz),
     $     visc2barD(lz1,nelz),
     $     areabarD(lz1,nelz),
     $     wght(lz1,nelz),
     $     zbar(lz1,nelz),
     $     work(lz1,nelz)
      
 
        real     dxux(lx1,ly1,lz1,lelv),
     $     dxuz(lx1,ly1,lz1,lelv),dyux(lx1,ly1,lz1,lelv),
     $     dyuz(lx1,ly1,lz1,lelv),dzux(lx1,ly1,lz1,lelv),
     $     dzuz(lx1,ly1,lz1,lelv),dxuy(lx1,ly1,lz1,lelv),
     $     dyuy(lx1,ly1,lz1,lelv),dzuy(lx1,ly1,lz1,lelv),
     $     dxth(lx1,ly1,lz1,lelv),dyth(lx1,ly1,lz1,lelv),
     $     dzth(lx1,ly1,lz1,lelv),px(lx1,ly1,lz1,lelv),
     $     py(lx1,ly1,lz1,lelv),pz(lx1,ly1,lz1,lelv),
     $     pm1(lx1,ly1,lz1,lelv),fakevec(lx1,ly1,lz1,lelv)

      common /csij/
     $     sij(lx1*ly1*lz1,6,lelv), ! 6 components
     $     ur(lx1*ly1*lz1), us(lx1*ly1*lz1), ut(lx1*ly1*lz1),
     $     vr(lx1*ly1*lz1), vs(lx1*ly1*lz1), vt(lx1*ly1*lz1),
     $     wr(lx1*ly1*lz1), ws(lx1*ly1*lz1), wt(lx1*ly1*lz1)

       vcut = 0.0D0 
       n = nx1*ny1*nz1*nelv
       call gradm1(dxux,dyux,dzux,vx) ! 1st derivative
       call gradm1(dxuy,dyuy,dzuy,vy)
       call gradm1(dxuz,dyuz,dzuz,vz)
       call gradm1(dxth,dyth,dzth,t)

        if (lx2.ne.lx1) then
           call mappr(pm1,pr,px,py)
           call gradm1(px,py,pz,pm1)
        else
           call gradm1(px,py,pz,pm1)
        endif

        call rzero(uxbarU,nz1*nelz) 
        call rzero(uybarU,nz1*nelz) 
        call rzero(uzbarU,nz1*nelz) !sets to zero
        call rzero(ux2barU,nz1*nelz) 
        call rzero(uy2barU,nz1*nelz) 
        call rzero(uz2barU,nz1*nelz) 
        call rzero(KEbarU,nz1*nelz) 
        call rzero(dzTbarU,nz1*nelz) 
        call rzero(TbarU,nz1*nelz) 
        call rzero(T2barU,nz1*nelz) 
        call rzero(FkbarU,nz1*nelz) 
        call rzero(FcbarU,nz1*nelz) 
        call rzero(FebarU,nz1*nelz) 
        call rzero(kappabarU,nz1*nelz) 
        call rzero(vortbarU,nz1*nelz) 
        call rzero(areabarU,nz1*nelz)
        call rzero(vortxbarU,nz1*nelz)
        call rzero(vortybarU,nz1*nelz)
        call rzero(visc3barU,nz1*nelz)
        call rzero(thermbarU,nz1*nelz)
        call rzero(PbarU,nz1*nelz)
        call rzero(dzPbarU,nz1*nelz)
        call rzero(PfluxbarU,nz1*nelz)
        call rzero(TdzPbarU,nz1*nelz)
        call rzero(trip1barU,nz1*nelz)
        call rzero(trip2barU,nz1*nelz)
        call rzero(rey12barU,nz1*nelz)
        call rzero(rey13barU,nz1*nelz)
        call rzero(rey23barU,nz1*nelz)
        call rzero(vfluxbarU,nz1*nelz)
        call rzero(kefluxbarU,nz1*nelz)
        call rzero(tfluxbarU,nz1*nelz)
        call rzero(thermfluxbarU,nz1*nelz) 
        call rzero(Fcflux1barU,nz1*nelz)
        call rzero(Fcflux2barU,nz1*nelz)
        call rzero(visc1barU,nz1*nelz)
        call rzero(visc2barU,nz1*nelz)

        call rzero(uxbarD,nz1*nelz) 
        call rzero(uybarD,nz1*nelz) 
        call rzero(uzbarD,nz1*nelz) !sets to zero
        call rzero(ux2barD,nz1*nelz) 
        call rzero(uy2barD,nz1*nelz) 
        call rzero(uz2barD,nz1*nelz) 
        call rzero(KEbarD,nz1*nelz) 
        call rzero(dzTbarD,nz1*nelz) 
        call rzero(TbarD,nz1*nelz) 
        call rzero(T2barD,nz1*nelz) 
        call rzero(FkbarD,nz1*nelz) 
        call rzero(FcbarD,nz1*nelz) 
        call rzero(FebarD,nz1*nelz) 
        call rzero(kappabarD,nz1*nelz) 
        call rzero(vortbarD,nz1*nelz) 
        call rzero(areabarD,nz1*nelz)
        call rzero(vortxbarD,nz1*nelz)
        call rzero(vortybarD,nz1*nelz)
        call rzero(visc3barD,nz1*nelz)
        call rzero(thermbarD,nz1*nelz)
        call rzero(PbarD,nz1*nelz)
        call rzero(dzPbarD,nz1*nelz)
        call rzero(PfluxbarD,nz1*nelz)
        call rzero(TdzPbarD,nz1*nelz)
        call rzero(trip1barD,nz1*nelz)
        call rzero(trip2barD,nz1*nelz)
        call rzero(rey12barD,nz1*nelz)
        call rzero(rey13barD,nz1*nelz)
        call rzero(rey23barD,nz1*nelz)
        call rzero(vfluxbarD,nz1*nelz)
        call rzero(kefluxbarD,nz1*nelz)
        call rzero(tfluxbarD,nz1*nelz)
        call rzero(thermfluxbarD,nz1*nelz) 
        call rzero(Fcflux1barD,nz1*nelz)
        call rzero(Fcflux2barD,nz1*nelz)
        call rzero(visc1barD,nz1*nelz)
        call rzero(visc2barD,nz1*nelz)

        call rzero(wght,nz1*nelz)
        call rzero(zbar,nz1*nelz)
        
        nelxy = nelgv/nelz ! Number of elements in xy-plane
        f = 5  ! Pick face 5 to evaluate surface Jacobian (lower face)
        do e=1,nelv
           eg = lglel(e)  !local->global element index
           call get_exyz(ex,ey,ez,eg,nelxy,1,nelz)
           do k=1,nz1
           do i=1,nx1*ny1
              Aval = area(i,1,f,e)
              wght(k,ez) = wght(k,ez)+Aval
              zbar(k,ez) = zbar(k,ez) + Aval*zval
              zval = zm1(i,1,k,e)
              uxval = vx(i,1,k,e)
              uyval = vy(i,1,k,e)
              uzval = vz(i,1,k,e)
              ux2val = uxval**2
              uy2val = uyval**2
              uz2val = uzval**2
              KEval = 0.5D0*(ux2val + uy2val + uz2val)
              dzTval = dzth(i,1,k,e)
              Tval = t(i,1,k,e,1)
              kapval = kappafunc(zval,param(9),param(6),
     $              param(79),param(5),-1.0D0)
             vortval = dxuy(i,1,k,e) - dyux(i,1,k,e)

              vortxval = dyuz(i,1,k,e) - dzuy(i,1,k,e)
              vortyval = dzux(i,1,k,e) - dxuz(i,1,k,e)  
             visc3val = -param(2)* (
     $       dxux(i,1,k,e)**2 + dyux(i,1,k,e)**2 + dzux(i,1,k,e)**2
     $   +   dxuy(i,1,k,e)**2 + dyuy(i,1,k,e)**2 + dzuy(i,1,k,e)**2
     $   +   dxuz(i,1,k,e)**2 + dyuz(i,1,k,e)**2 + dzuz(i,1,k,e)**2 )

             thermval = -kapval*(
     $       dxth(i,1,k,e)**2 + dyth(i,1,k,e)**2 + dzth(i,1,k,e)**2)

              Pval = pm1(i,1,k,e)
              dzPval = pz(i,1,k,e) 
              Pfluxval = uzval*pval
              TdzPval = Tval*dzPval 
              trip1val = Tval*uzval*uzval
              trip2val = Tval*Tval*uzval

              rey12val = uxval*uyval
              rey13val = uxval*uzval
              rey23val = uyval*uzval

              vfluxval = param(2)*dzuz(i,1,k,e) 
              kefluxval = param(2)*( dzuz(i,1,k,e)*uzval +
     $        + dzux(i,1,k,e)*uxval + dzuy(i,1,k,e)*uyval)
              tfluxval = kapval*dzTval 
              thermfluxval = kapval*Tval*dzTval 
              Fcflux1val = param(2)*Tval*dzuz(i,1,k,e)
              Fcflux2val = kapval*uzval*dzTval
              visc1val =  dxuz(i,1,k,e)*dxth(i,1,k,e)
     $        + dyuz(i,1,k,e)*dyth(i,1,k,e) 
     $        + dzuz(i,1,k,e)*dzth(i,1,k,e)

              visc2val = kapval*visc1val
              visc1val = param(2)*visc1val




             if(uzval.ge.vcut) then
                uxbarU(k,ez) = uxbarU(k,ez) + Aval*uxval
                uybarU(k,ez) = uybarU(k,ez) + Aval*uyval
                ux2barU(k,ez) = ux2barU(k,ez) + Aval*ux2val
                uzbarU(k,ez) = uzbarU(k,ez) + Aval*uzval
                uy2barU(k,ez) = uy2barU(k,ez) + Aval*uy2val
                uz2barU(k,ez) = uz2barU(k,ez) + Aval*uz2val
                KEbarU(k,ez) = KEbarU(k,ez) + Aval*KEval
                dzTbarU(k,ez) = dzTbarU(k,ez) + Aval*dzTval
                TbarU(k,ez) = TbarU(k,ez) + Aval*Tval
                T2barU(k,ez) = T2barU(k,ez) + Aval*Tval*Tval
                FkbarU(k,ez) = FkbarU(k,ez) - Aval*kapval*dzTval
                FcbarU(k,ez) = FcbarU(k,ez) + Aval*Tval*uzval
                FebarU(k,ez) = FebarU(k,ez) + Aval*KEval*uzval
                kappabarU(k,ez) = kappabarU(k,ez) + Aval*kapval
                vortbarU(k,ez) = vortbarU(k,ez) + Aval*vortval
                vortxbarU(k,ez) = vortxbarU(k,ez) + Aval*vortxval
                vortybarU(k,ez) = vortybarU(k,ez)+ Aval*vortyval
                visc3barU(k,ez) =  visc3barU(k,ez)+ Aval*visc3val
                thermbarU(k,ez) = thermbarU(k,ez) + Aval*thermval
                PbarU(k,ez) = PbarU(k,ez) + Aval*Pval
                dzPbarU(k,ez) =  dzPbarU(k,ez)+ Aval*dzPval
                PfluxbarU(k,ez) = PfluxbarU(k,ez) + Aval*Pfluxval
                TdzPbarU(k,ez) = TdzPbarU(k,ez) + Aval*TdzPval
                trip1barU(k,ez) = trip1barU(k,ez) + Aval*trip1val
                trip2barU(k,ez) = trip2barU(k,ez) + Aval*trip2val
                rey12barU(k,ez) = rey12barU(k,ez)+ Aval*rey12val
                rey13barU(k,ez) = rey13barU(k,ez) + Aval*rey13val
                rey23barU(k,ez) = rey23barU(k,ez) + Aval*rey23val
                vfluxbarU(k,ez) = vfluxbarU(k,ez) + Aval*vfluxval
                kefluxbarU(k,ez) = kefluxbarU(k,ez) + Aval*kefluxval
                tfluxbarU(k,ez) = tfluxbarU(k,ez) + Aval*tfluxval
                thermfluxbarU(k,ez) =  thermfluxbarU(k,ez) + Aval*thermfluxval
                Fcflux1barU(k,ez) = Fcflux1barU(k,ez) + Aval*Fcflux1val
                Fcflux2barU(k,ez) = Fcflux2barU(k,ez) + Aval*Fcflux2val
                visc1barU(k,ez) = visc1barU(k,ez)+ Aval*visc1val
                visc2barU(k,ez) = visc2barU(k,ez) + Aval*visc2val
                areabarU(k,ez) = areabarU(k,ez) + Aval*Aval
              else
                uxbarD(k,ez) = uxbarD(k,ez) + Aval*uxval
                uybarD(k,ez) = uybarD(k,ez) + Aval*uyval
                ux2barD(k,ez) = ux2barD(k,ez) + Aval*ux2val
                uzbarD(k,ez) = uzbarD(k,ez) + Aval*uzval
                uy2barD(k,ez) = uy2barD(k,ez) + Aval*uy2val
                uz2barD(k,ez) = uz2barD(k,ez) + Aval*uz2val
                KEbarD(k,ez) = KEbarD(k,ez) + Aval*KEval
                dzTbarD(k,ez) = dzTbarD(k,ez) + Aval*dzTval
                TbarD(k,ez) = TbarD(k,ez) + Aval*Tval
                T2barD(k,ez) = T2barD(k,ez) + Aval*Tval*Tval
                FkbarD(k,ez) = FkbarD(k,ez) - Aval*kapval*dzTval
                FcbarD(k,ez) = FcbarD(k,ez) + Aval*Tval*uzval
                FebarD(k,ez) = FebarD(k,ez) + Aval*KEval*uzval
                kappabarD(k,ez) = kappabarD(k,ez) + Aval*kapval
                vortbarD(k,ez) = vortbarD(k,ez) + Aval*vortval
                vortxbarD(k,ez) = vortxbarD(k,ez) + Aval*vortxval
                vortybarD(k,ez) = vortybarD(k,ez)+ Aval*vortyval
                visc3barD(k,ez) =  visc3barD(k,ez)+ Aval*visc3val
                thermbarD(k,ez) = thermbarD(k,ez) + Aval*thermval
                PbarD(k,ez) = PbarD(k,ez) + Aval*Pval
                dzPbarD(k,ez) =  dzPbarD(k,ez)+ Aval*dzPval
                PfluxbarD(k,ez) = PfluxbarD(k,ez) + Aval*Pfluxval
                TdzPbarD(k,ez) = TdzPbarD(k,ez) + Aval*TdzPval
                trip1barD(k,ez) = trip1barD(k,ez) + Aval*trip1val
                trip2barD(k,ez) = trip2barD(k,ez) + Aval*trip2val
                rey12barD(k,ez) = rey12barD(k,ez)+ Aval*rey12val
                rey13barD(k,ez) = rey13barD(k,ez) + Aval*rey13val
                rey23barD(k,ez) = rey23barD(k,ez) + Aval*rey23val
                vfluxbarD(k,ez) = vfluxbarD(k,ez) + Aval*vfluxval
                kefluxbarD(k,ez) = kefluxbarD(k,ez) + Aval*kefluxval
                tfluxbarD(k,ez) = tfluxbarD(k,ez) + Aval*tfluxval
                thermfluxbarD(k,ez) =  thermfluxbarD(k,ez) + Aval*thermfluxval
                Fcflux1barD(k,ez) = Fcflux1barD(k,ez) + Aval*Fcflux1val
                Fcflux2barD(k,ez) = Fcflux2barD(k,ez) + Aval*Fcflux2val
                visc1barD(k,ez) = visc1barD(k,ez)+ Aval*visc1val
                visc2barD(k,ez) = visc2barD(k,ez) + Aval*visc2val
                areabarD(k,ez) = areabarD(k,ez) + Aval*Aval
            endif
           enddo
           enddo
        enddo

        call gop(wght,work,'+  ',nz1*nelz) !2 spaces after + sign needed in both
        call gop(zbar,work,'+  ',nz1*nelz) !note that these are NOT evenly spaced due to GLL point distribution
        call gop(uxbarU,work,'+  ',nz1*nelz) 
        call gop(uybarU,work,'+  ',nz1*nelz) 
        call gop(uzbarU,work,'+  ',nz1*nelz) !gather over processors (MPI_AllReduce )
        call gop(ux2barU,work,'+  ',nz1*nelz) 
        call gop(uy2barU,work,'+  ',nz1*nelz) 
        call gop(uz2barU,work,'+  ',nz1*nelz) 
        call gop(KEbarU,work,'+  ',nz1*nelz) 
        call gop(dzTbarU,work,'+  ',nz1*nelz) 
        call gop(TbarU,work,'+  ',nz1*nelz) 
        call gop(T2barU,work,'+  ',nz1*nelz) 
        call gop(FkbarU,work,'+  ',nz1*nelz) 
        call gop(FcbarU,work,'+  ',nz1*nelz) 
        call gop(FebarU,work,'+  ',nz1*nelz) 
        call gop(kappabarU,work,'+  ',nz1*nelz) 
        call gop(vortbarU,work,'+  ',nz1*nelz) 
        call gop(areabarU,work,'+  ',nz1*nelz) 
        call gop(vortxbarU,work,'+  ',nz1*nelz)
        call gop(vortybarU,work,'+  ',nz1*nelz)
        call gop(visc3barU,work,'+  ',nz1*nelz)
        call gop(thermbarU,work,'+  ',nz1*nelz)
        call gop(PbarU,work,'+  ',nz1*nelz)
        call gop(dzPbarU,work,'+  ',nz1*nelz)
        call gop(PfluxbarU,work,'+  ',nz1*nelz)
        call gop(TdzPbarU,work,'+  ',nz1*nelz)
        call gop(trip1barU,work,'+  ',nz1*nelz)
        call gop(trip2barU,work,'+  ',nz1*nelz)
        call gop(rey12barU,work,'+  ',nz1*nelz)
        call gop(rey13barU,work,'+  ',nz1*nelz)
        call gop(rey23barU,work,'+  ',nz1*nelz)
        call gop(vfluxbarU,work,'+  ',nz1*nelz)
        call gop(kefluxbarU,work,'+  ',nz1*nelz)
        call gop(tfluxbarU,work,'+  ',nz1*nelz)
        call gop(thermfluxbarU,work,'+  ',nz1*nelz)
        call gop(Fcflux1barU,work,'+  ',nz1*nelz)
        call gop(Fcflux2barU,work,'+  ',nz1*nelz)
        call gop(visc1barU,work,'+  ',nz1*nelz)
        call gop(visc2barU,work,'+  ',nz1*nelz)

        call gop(uxbarD,work,'+  ',nz1*nelz) 
        call gop(uybarD,work,'+  ',nz1*nelz) 
        call gop(uzbarD,work,'+  ',nz1*nelz) !gather over processors (MPI_AllReduce )
        call gop(ux2barD,work,'+  ',nz1*nelz) 
        call gop(uy2barD,work,'+  ',nz1*nelz) 
        call gop(uz2barD,work,'+  ',nz1*nelz) 
        call gop(KEbarD,work,'+  ',nz1*nelz) 
        call gop(dzTbarD,work,'+  ',nz1*nelz) 
        call gop(TbarD,work,'+  ',nz1*nelz) 
        call gop(T2barD,work,'+  ',nz1*nelz) 
        call gop(FkbarD,work,'+  ',nz1*nelz) 
        call gop(FcbarD,work,'+  ',nz1*nelz) 
        call gop(FebarD,work,'+  ',nz1*nelz) 
        call gop(kappabarD,work,'+  ',nz1*nelz) 
        call gop(vortbarD,work,'+  ',nz1*nelz) 
        call gop(areabarD,work,'+  ',nz1*nelz) 
        call gop(vortxbarD,work,'+  ',nz1*nelz)
        call gop(vortybarD,work,'+  ',nz1*nelz)
        call gop(visc3barD,work,'+  ',nz1*nelz)
        call gop(thermbarD,work,'+  ',nz1*nelz)
        call gop(PbarD,work,'+  ',nz1*nelz)
        call gop(dzPbarD,work,'+  ',nz1*nelz)
        call gop(PfluxbarD,work,'+  ',nz1*nelz)
        call gop(TdzPbarD,work,'+  ',nz1*nelz)
        call gop(trip1barD,work,'+  ',nz1*nelz)
        call gop(trip2barD,work,'+  ',nz1*nelz)
        call gop(rey12barD,work,'+  ',nz1*nelz)
        call gop(rey13barD,work,'+  ',nz1*nelz)
        call gop(rey23barD,work,'+  ',nz1*nelz)
        call gop(vfluxbarD,work,'+  ',nz1*nelz)
        call gop(kefluxbarD,work,'+  ',nz1*nelz)
        call gop(tfluxbarD,work,'+  ',nz1*nelz)
        call gop(thermfluxbarD,work,'+  ',nz1*nelz)
        call gop(Fcflux1barD,work,'+  ',nz1*nelz)
        call gop(Fcflux2barD,work,'+  ',nz1*nelz)
        call gop(visc1barD,work,'+  ',nz1*nelz)
        call gop(visc2barD,work,'+  ',nz1*nelz)

        do i=1,nz1*nelz
           uxbarU(i,1)=uxbarU(i,1)/wght(i,1) ! average
           uybarU(i,1)=uybarU(i,1)/wght(i,1) 
           uzbarU(i,1)=uzbarU(i,1)/wght(i,1) 
           ux2barU(i,1)=ux2barU(i,1)/wght(i,1) 
           uy2barU(i,1)=uy2barU(i,1)/wght(i,1) 
           uz2barU(i,1)=uz2barU(i,1)/wght(i,1) 
           KEbarU(i,1)=KEbarU(i,1)/wght(i,1) 
           dzTbarU(i,1)=dzTbarU(i,1)/wght(i,1) 
           TbarU(i,1)=TbarU(i,1)/wght(i,1) 
           T2barU(i,1)=T2barU(i,1)/wght(i,1) 
           FkbarU(i,1)=FkbarU(i,1)/wght(i,1) 
           FcbarU(i,1)=FcbarU(i,1)/wght(i,1) 
           FebarU(i,1)=FebarU(i,1)/wght(i,1) 
           kappabarU(i,1)=kappabarU(i,1)/wght(i,1) 
           vortbarU(i,1) = vortbarU(i,1)/wght(i,1) 
           areabarU(i,1) = areabarU(i,1)/wght(i,1) 
           vortxbarU(i,1) = vortxbarU(i,1) /wght(i,1)
           vortybarU(i,1) = vortybarU(i,1)/wght(i,1)
           visc3barU(i,1) =  visc3barU(i,1)/wght(i,1)
           thermbarU(i,1) = thermbarU(i,1) /wght(i,1)
           PbarU(i,1) = PbarU(i,1) /wght(i,1)
           dzPbarU(i,1) =  dzPbarU(i,1)/wght(i,1)
           PfluxbarU(i,1) = PfluxbarU(i,1) /wght(i,1)
           TdzPbarU(i,1) = TdzPbarU(i,1) /wght(i,1)
           trip1barU(i,1) = trip1barU(i,1) /wght(i,1)
           trip2barU(i,1) = trip2barU(i,1) /wght(i,1)
           rey12barU(i,1) = rey12barU(i,1)/wght(i,1)
           rey13barU(i,1) = rey13barU(i,1) /wght(i,1)
           rey23barU(i,1) = rey23barU(i,1) /wght(i,1)
           vfluxbarU(i,1) = vfluxbarU(i,1) /wght(i,1)
           kefluxbarU(i,1) = kefluxbarU(i,1) /wght(i,1)
           tfluxbarU(i,1) = tfluxbarU(i,1) /wght(i,1)
           thermfluxbarU(i,1) =  thermfluxbarU(i,1) /wght(i,1)
           Fcflux1barU(i,1) = Fcflux1barU(i,1) /wght(i,1)
           Fcflux2barU(i,1) = Fcflux2barU(i,1) /wght(i,1)
           visc1barU(i,1) = visc1barU(i,1)/wght(i,1)
           visc2barU(i,1) = visc2barU(i,1) /wght(i,1)

           uxbarD(i,1)=uxbarD(i,1)/wght(i,1) ! average
           uybarD(i,1)=uybarD(i,1)/wght(i,1) 
           uzbarD(i,1)=uzbarD(i,1)/wght(i,1) 
           ux2barD(i,1)=ux2barD(i,1)/wght(i,1) 
           uy2barD(i,1)=uy2barD(i,1)/wght(i,1) 
           uz2barD(i,1)=uz2barD(i,1)/wght(i,1) 
           KEbarD(i,1)=KEbarD(i,1)/wght(i,1) 
           dzTbarD(i,1)=dzTbarD(i,1)/wght(i,1) 
           TbarD(i,1)=TbarD(i,1)/wght(i,1) 
           T2barD(i,1)=T2barD(i,1)/wght(i,1) 
           FkbarD(i,1)=FkbarD(i,1)/wght(i,1) 
           FcbarD(i,1)=FcbarD(i,1)/wght(i,1) 
           FebarD(i,1)=FebarD(i,1)/wght(i,1) 
           kappabarD(i,1)=kappabarD(i,1)/wght(i,1) 
           vortbarD(i,1) = vortbarD(i,1)/wght(i,1) 
           areabarD(i,1) = areabarD(i,1)/wght(i,1) 
           vortxbarD(i,1) = vortxbarD(i,1) /wght(i,1)
           vortybarD(i,1) = vortybarD(i,1)/wght(i,1)
           visc3barD(i,1) =  visc3barD(i,1)/wght(i,1)
           thermbarD(i,1) = thermbarD(i,1) /wght(i,1)
           PbarD(i,1) = PbarD(i,1) /wght(i,1)
           dzPbarD(i,1) =  dzPbarD(i,1)/wght(i,1)
           PfluxbarD(i,1) = PfluxbarD(i,1) /wght(i,1)
           TdzPbarD(i,1) = TdzPbarD(i,1) /wght(i,1)
           trip1barD(i,1) = trip1barD(i,1) /wght(i,1)
           trip2barD(i,1) = trip2barD(i,1) /wght(i,1)
           rey12barD(i,1) = rey12barD(i,1)/wght(i,1)
           rey13barD(i,1) = rey13barD(i,1) /wght(i,1)
           rey23barD(i,1) = rey23barD(i,1) /wght(i,1)
           vfluxbarD(i,1) = vfluxbarD(i,1) /wght(i,1)
           kefluxbarD(i,1) = kefluxbarD(i,1) /wght(i,1)
           tfluxbarD(i,1) = tfluxbarD(i,1) /wght(i,1)
           thermfluxbarD(i,1) =  thermfluxbarD(i,1) /wght(i,1)
           Fcflux1barD(i,1) = Fcflux1barD(i,1) /wght(i,1)
           Fcflux2barD(i,1) = Fcflux2barD(i,1) /wght(i,1)
           visc1barD(i,1) = visc1barD(i,1)/wght(i,1)
           visc2barD(i,1) = visc2barD(i,1) /wght(i,1)
           zbar(i,1)=zbar(i,1)/wght(i,1)
        enddo

        if(nid .eq. 0) then !root process
         open(51,file="zavg_up.dat",access="append")
          do i=1,nz1*nelz
           write(51,'(1p40E15.7)') time,
     $      zbar(i,1),
     $      uxbarU(i,1),
     $      uybarU(i,1),
     $      uzbarU(i,1),
     $      ux2barU(i,1),
     $      uy2barU(i,1),
     $      uz2barU(i,1),
     $      KEbarU(i,1),
     $      dzTbarU(i,1),
     $      TbarU(i,1),
     $      T2barU(i,1),
     $      FkbarU(i,1),
     $      FcbarU(i,1),
     $      FebarU(i,1),
     $      kappabarU(i,1),
     $      vortbarU(i,1),
     $      areabarU(i,1),
     $      vortxbarU(i,1),
     $      vortybarU(i,1),
     $      visc3barU(i,1),
     $      thermbarU(i,1),
     $      PbarU(i,1),
     $      dzPbarU(i,1),
     $      PfluxbarU(i,1),
     $      TdzPbarU(i,1),
     $      trip1barU(i,1),
     $      trip2barU(i,1),
     $      rey12barU(i,1),
     $      rey13barU(i,1),
     $      rey23barU(i,1),
     $      vfluxbarU(i,1),
     $      kefluxbarU(i,1),
     $      tfluxbarU(i,1),
     $      thermfluxbarU(i,1),
     $      Fcflux1barU(i,1),
     $      Fcflux2barU(i,1),
     $      visc1barU(i,1),
     $      visc2barU(i,1)
          enddo
        close(51)
         open(51,file="zavg_down.dat",access="append")
          do i=1,nz1*nelz
           write(51,'(1p40E15.7)') time,
     $      zbar(i,1),
     $      uxbarD(i,1),
     $      uybarD(i,1),
     $      uzbarD(i,1),
     $      ux2barD(i,1),
     $      uy2barD(i,1),
     $      uz2barD(i,1),
     $      KEbarD(i,1),
     $      dzTbarD(i,1),
     $      TbarD(i,1),
     $      T2barD(i,1),
     $      FkbarD(i,1),
     $      FcbarD(i,1),
     $      FebarD(i,1),
     $      kappabarD(i,1),
     $      vortbarD(i,1),
     $      areabarD(i,1),
     $      vortxbarD(i,1),
     $      vortybarD(i,1),
     $      visc3barD(i,1),
     $      thermbarD(i,1),
     $      PbarD(i,1),
     $      dzPbarD(i,1),
     $      PfluxbarD(i,1),
     $      TdzPbarD(i,1),
     $      trip1barD(i,1),
     $      trip2barD(i,1),
     $      rey12barD(i,1),
     $      rey13barD(i,1),
     $      rey23barD(i,1),
     $      vfluxbarD(i,1),
     $      kefluxbarD(i,1),
     $      tfluxbarD(i,1),
     $      thermfluxbarD(i,1),
     $      Fcflux1barD(i,1),
     $      Fcflux2barD(i,1),
     $      visc1barD(i,1),
     $      visc2barD(i,1)
          enddo
        close(51)
        endif
      return
      end
c-----------------------------------------------------------------------
      subroutine output_vort
       include 'SIZE'
       include 'TOTAL'
      integer eg,ex,ey,ez,f,n,e
      integer iflag,nelxy
      parameter(lt=lx1*ly1*lz1*lelv)  
      parameter(nelz=8) ! Number of elements in z direction 
      real Tbar(lz1,nelz), wght(lz1,nelz),work(lz1,nelz)
      real Tf(lx1,ly1,lz1,lelv)
      common /scrns/ vort(lt,3),w1(lt),w2(lt)
      
 

       n = nx1*ny1*nz1*nelv

        call rzero(wght,nz1*nelz)
        call rzero(Tbar,nz1*nelz)
        
        nelxy = nelgv/nelz ! Number of elements in xy-plane
        f = 5  ! Pick face 5 to evaluate surface Jacobian (lower face)
        do e=1,nelv
           eg = lglel(e)  !local->global element index
           call get_exyz(ex,ey,ez,eg,nelxy,1,nelz)
           do k=1,nz1
           do i=1,nx1*ny1
              wght(k,ez) = wght(k,ez) +area(i,1,f,e)
              Tbar(k,ez) = Tbar(k,ez) + area(i,1,f,e)*t(i,1,k,e,1)
           enddo
           enddo
        enddo

        call gop(wght,work,'+  ',nz1*nelz) !2 spaces after + sign needed in both
        call gop(Tbar,work,'+  ',nz1*nelz) 

        do i=1,nz1*nelz
           Tbar(i,1)=Tbar(i,1)/wght(i,1) 
        enddo

        do e=1,nelv
           eg = lglel(e)  !local->global element index
           call get_exyz(ex,ey,ez,eg,nelxy,1,nelz)
           do k=1,nz1
           do i=1,nx1*ny1
              Tf(i,1,k,e) = t(i,1,k,e,1) - Tbar(k,ez)
           enddo
           enddo
        enddo


c        call copy(u,t,n)
c        call opcopy(u,v,w,vx,vy,vz)
c       call no_z_profile(vz)

       call comp_vort3(vort,w1,w2,vx,vy,vz)

       call outpost(vort(1,3),vort(1,2),vort(1,1),pr,Tf,'vor')
c       call outpost(vx,vy,vz,pr,t,'vor')


      return
      end
c-----------------------------------------------------------------------
       real function kappafunc_power(z,k0, n)
       real z, k0, n, res
        

            res = k0*(
     $            (1.0D0 + z)**(-n)
     $          + (2.0D0-z)**(-n)
     $            )/(1.0D0 + 2.0D0**(-n))

        kappafunc_power = res
        return
        end
c----------------------------------------------------------------------
       real function kappafunc_tanh(z,k0, k1, delta) 
c             tanh kappa profile
c             CZ between 0 and 1
       real z, delta, k0, k1, res

        res = LOG(k1) + 0.5D0*LOG(k1/k0)
     $   * (TANH((z-1.0D0)/delta-2.0D0) - TANH(z/delta+2.0D0))
        res = EXP(res)
c        write(*,*) 'K', z, k0, k1, delta, res
        kappafunc_tanh = res
        return
        end
c----------------------------------------------------------------------
       real function kappafunc_const(z,k0) 
c             constant kappa profile
       real z, k0, res

        kappafunc_const = k0
        return
        end
c----------------------------------------------------------------------
c Solar kappa profile 
       real function kappafunc(z,fmin,dval,aval,k0,zmval)  
c        fmin = param(9)
c        k0  = param(5)
c        dval = param(6)
c        aval = param(79)
c        zmval = min z
        real z,fmin,aval,res,dval,zmval,k0
        real zcval,xval,x1val
        real xmval,xm1val


        zcval = (1.0D0-fmin)/aval
        xval = (z - zcval)/dval
        xmval = (zmval-zcval)/dval
        x1val = (1.0D0 - z -zcval)/dval 
        xm1val = (1.0D0-zmval-zcval)/dval


        fac1 = 1.0D0 + exp(xval)
        fac1 = fac1 / (1.0D0 + exp(xmval))

        fac1 = fac1 * (1.0D0 + exp(-x1val))
        fac1 = fac1 / (1.0D0 + exp(-xm1val))

        res = 1.0D0 - aval*zmval - aval*(z-zmval)
c        if(nid.eq.0) then
c            write(*,*) z,zcval,xval,xmval,res,fac1
c        endif

        res = res + aval*dval*log(fac1)

c        if(nid.eq.0) then
c            write(*,*) z,fmin,dval,aval,k0,zmval,res
c        endif

        kappafunc = res / k0
        return
       end
c----------------------------------------------------------------------
       real function invkappafunc(z,k0,dl)  
        real kappafunc
        invkappafunc = 1.0D0/kappafunc(z,k0,dl)
        return
       end
c----------------------------------------------------------------------
c AMD function to do simple trapezoidal integration over the function
c func
      real function trapz(func,a,b,n, arg1,arg2)
      real a, b
      real dx, dx2, xcurr
      real f
      real res
      integer n
      real arg1, arg2, arg3 
      real func

       dx = (b-a)/n
       dx2 = dx/2.0D0
    
       xcurr = a
        
       f = func(xcurr,arg1,arg2)


       res = dx2*f

       xcurr = xcurr + dx


       do i=1,n-2
         f = func(xcurr,arg1,arg2)
         res = res + dx*f
         xcurr = xcurr + dx
       enddo
       f =  func(xcurr,arg1,arg2)
       res = res + dx2*f 
       trapz = res
       return
      END
c-----------------------------------------------------------------------
      real function invtrapz(func,a,b,n,arg1,arg2,arg3)
      real a, b, res, f
      real dx, dx2, xcurr
      real arg1, arg2, arg3
      integer n
      real func


       dx = (b-a)/n
       dx2 = dx/2.0D0
    
       xcurr = a
        
       f = func(xcurr, arg1, arg2,arg3)


       res = dx2/f

       xcurr = xcurr + dx


       do i=1,n-2
         f =  func(xcurr, arg1, arg2, arg3)
         res = res + dx/f
         xcurr = xcurr + dx
       enddo
       f = func(xcurr, arg1, arg2, arg3)
       res = res + dx2/f 
       invtrapz = res
       return
      END
c-----------------------------------------------------------------------

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
