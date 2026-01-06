!----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
      !################################################################
      !Internal Variables
      !################################################################
      integer*4 i 
      real*8 rpi, rmass, dmass, vmag, rep 
      real*8 fqs, fqsx, fqsy,fqsz, fbx, fby, fbz  
      !################################################################
      !External Variables
      !################################################################
      real*8 rmu,rhof,rg
      common /parameters/ rmu,rhof,rg

c      real*8 rhop,dp
c      common /dpinfo/ dp,rhop

      real*8 mup,sigma
      common /tension/ mup,sigma      
!
      real*8 rlz,rrz
      common /domainsize/ rlz,rrz
!
      real*8 dpmin
      common /dpm/ dpmin 
!
      real*8 y_out(PPICLF_LRS         , PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP + 1 , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out
      !################################################################
      !Langevin Variables
      !################################################################
      real*8 lan_c0,lan_cy,lan_cs,lan_S,lan_cg
      real*8 lan_gij,lan_cw,lan_e,f_w
      common /langevin/ lan_c0,lan_cy,lan_cs,lan_S,
     >              lan_cg,lan_gij,lan_cw,lan_e,f_w
      real*8 ppiclf_ran
      !################################################################
      !Droplet Breakup Variables
      !################################################################
      ! Sam - TAB
      real*8 cd_tab, ck_tab, cf_tab, c_tab, cb_tab, r0, y0, ydot0, 
     >       We, td, omega, k_tab, v_child, cv_tab,dp_maj,
     >       nvec1, nvec2, mvec1, mvec2, mvec3, theta_child, rvelx,
     >       rvely, rvelz, F_cutoff, F_parent, P_chi, Q_chi, bound_chi, 
     >       chisq, CMR, rmm, sigma_g, rhat, r32hat, rphat, rchat,
     >       D32, Re32, A_part, V_part, e_corr, r_minor, rmbar,
     >       x_tab, f_interp, CD_disk,CD_sph, f_diskx, f_disky, f_diskz,
     >       vol_children, h_child, LCMR, Lsigma_g, C_spl, Ftab,
     >       sigma_time, tb_mean_est, tb_tab, y_tab, ydot_tab,
     >       yddot_tab, sigma_y, yc_tab, SPL, f1_tabs, f2_tabs,
     >       s1_tabs, g1_tabs, g2_tabs, g3_tabs, k_tabs, omega_tilde,
     >       alpha_tabs

      real*8, dimension(PPICLF_LRP, PPICLF_LPART) :: add_rprop
      real*8, dimension(PPICLF_LRS, PPICLF_LPART) :: add_y
      real*8, dimension(PPICLF_LRS+6, PPICLF_LPART) :: log_spawn

      real*8 r32, r_child,  vol_child, r_cutoff, v_cutoff, 
     >     v_below,
     >     v_belowi, dF_below, F_belowi, uid

      data r_cutoff /5.0d-6/
      integer*4 n_child, n_added, 
     >          flipsign,iprop,j,k,ispl,n_child_i,ilog
      integer max_iter, status_chi, ierr, tindex
      logical add_child
      character(50) :: fname
      !Print Data to file
      integer :: breakup_counter
      logical I_EXIST
      character(len=25) str
      integer*4 f_dump
      save      f_dump
      data      f_dump /1/
      !################################################################
c      integer*4 icalld
c      save      icalld
c      data      icalld /0/
c      if (icalld.eq.3) icalld= 0
c      icalld = icalld + 1
      open(1114, file='breakup_counts.txt')
      !#######################################################################
      !Variable Type Declaration End
      !#######################################################################
      ! initializations & Constants Declaration

      tindex = nint(ppiclf_time/ppiclf_dt)+1   !uid generation
      n_added = 0
      breakup_counter = 0

      rpi   = 4.0*atan(1.0)
      !rmu   = 1.8E-5
      ! TAB Constants Start
      cf_tab = 1.0/3.0
      ck_tab = 8.0
      cd_tab = 5.0
      cb_tab = 0.5
      cv_tab = 1.0
      c_tab = cf_tab / (cb_tab * ck_tab)
      k_tab = 10.0/3.0
      ! stochastic TAB model
      f1_tabs = 5.37
      f2_tabs = 0.450
      s1_tabs = 0.206
      g1_tabs = 8.06
      g2_tabs = 0.512
      g3_tabs = 2.31
      alpha_tabs = 2.58
      k_tabs = 0.17           ! value based on Flock et al. 2012
      sigma_g = exp(0.31)     ! Jain et al. 2015
      Lsigma_g = log(sigma_g) ! Secondary breakup of a drop at moderate Weber numbers

! evaluate ydot
!#######################################################################
!Breakup - Balachandar  !https://doi.org/10.1063/5.0131815
!O'Rourke and Amsden, 1987  ! http://dx.doi.org/10.1016/j.apm.2012.02.015
!#######################################################################
      do i=1,ppiclf_npart

c          SPL   = ppiclf_rprop(PPICLF_R_JSPL,i)
          rvelx = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
          rvely = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
          rvelz = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
          r0    = ppiclf_rprop(PPICLF_R_JR0TAB ,i)

          vmag  = sqrt(rvelx**2 + rvely**2 + rvelz**2)

          !Define unit normal vector to relative velocity (3d) such that n.v=0;|n|=1,n3=0
          nvec2 = 1/sqrt((rvely/rvelx)**2 + 1)
          nvec1 = -nvec2*rvely/rvelx
          !Define another unit normal vector orthogonal to nvec and relative velocity: mvec = v cross n
          mvec1 = -nvec2*rvelz
          mvec2 =  nvec1*rvelz
          mvec3 =  nvec2*rvelx - nvec1*rvely
          
          !Calculate new deformation
c          y0    = ppiclf_rprop(PPICLF_R_JYTAB,   i)
c          ydot0 = ppiclf_rprop(PPICLF_R_JYDOTTAB,i)

          vmag  = sqrt(rvelx**2 + rvely**2 + rvelz**2)
          ppiclf_rprop(PPICLF_R_JWEB,i) = rhof*(vmag**2)*r0/sigma
          We    = ppiclf_rprop(PPICLF_R_JWEB,i)

          ppiclf_rprop(PPICLF_R_JTD,i) =
     >         2*ppiclf_rprop(PPICLF_R_JRHOP,i)*r0**2/(cd_tab*mup)
          td    = ppiclf_rprop(PPICLF_R_JTD,i)

          ppiclf_rprop(PPICLF_R_JOMEGA ,i) = sqrt(ck_tab*sigma/
     >             (ppiclf_rprop(PPICLF_R_JRHOP,i)*r0**3) - 1/td**2)
          omega = ppiclf_rprop(PPICLF_R_JOMEGA ,i)

          y0    = 0.0
          ydot0 = 0.0
          
          ppiclf_rprop(PPICLF_R_JYTAB,i) = c_tab*We 
     >           + exp(-ppiclf_dt/td)
     >           * (1/omega*((y0-c_tab*We)/td+ydot0)
     >           * sin(omega * ppiclf_dt)
     >           + (y0 - c_tab*We)
     >           * cos(omega*ppiclf_dt))

          ppiclf_rprop(PPICLF_R_JYDOTTAB,i) = 
     >    (c_tab*We-ppiclf_rprop(PPICLF_R_JYTAB,i))/td
     >   + omega*exp(-ppiclf_dt/td)
     >   * (1/omega*((y0 - c_tab*We)/td+ydot0)*cos(omega*ppiclf_dt) 
     >   -                     (y0 - c_tab*We)*sin(omega*ppiclf_dt))

         !Update diameter for CD calculation
          x_tab = cb_tab*r0*ppiclf_rprop(PPICLF_R_JYTAB,i)
          ppiclf_rprop(PPICLF_R_JDMAJ,i) = 
     >           2.0*(ppiclf_rprop(PPICLF_R_JR0TAB,i) + x_tab)
         
          ppiclf_rprop(PPICLF_R_JFTAB,i) = ppiclf_ran(2)
c          Ftab = ppiclf_ran(2)!ppiclf_rprop(PPICLF_R_JFTAB,i)
c          ppiclf_rprop(PPICLF_R_JFTAB,i) = cf_tab*rhof*vmag**2
c     >                      *ppiclf_rprop(PPICLF_R_JVOLP,i)/r0
          Ftab = ppiclf_rprop(PPICLF_R_JFTAB,i)
          
          !Calculate standard deviation in y
          omega_tilde = omega*td
          sigma_y = f1_tabs*We**(-f2_tabs)*
     >            exp(-s1_tabs*2.0*rpi/omega_tilde) 
     >          - g1_tabs*We**(-g2_tabs) + g3_tabs

          !phi(sigma_y) & y_critical
          sigma_y = k_tabs*log(exp(alpha_tabs*sigma_y)+1)/alpha_tabs

          if (sigma_y .gt. 0.0) then
           call cdfnor(2, Ftab, 1.0-Ftab, yc_tab, 1.0d-5,
     >                 sigma_y, status_chi, bound_chi)
          else
           yc_tab = 1
          endif

          ppiclf_rprop(PPICLF_R_JSIGMA,i) = sigma_y
          ppiclf_rprop(PPICLF_R_JYCTAB,i) = yc_tab

          !Check for breakup (y>y_c)
          if (ppiclf_rprop(PPICLF_R_JYTAB ,i) .gt.
     >        ppiclf_rprop(PPICLF_R_JYCTAB,i)) then

          ppiclf_rprop(PPICLF_R_JRR,i) = r0 / (1.0 + 8.0*k_tab/20.0
     >           +ppiclf_rprop(PPICLF_R_JRHOP,i)*r0**3.0/sigma 
     >           *ppiclf_rprop(PPICLF_R_JYDOTTAB,i)**2.0
     >           *(6.0*k_tab-5.0)/120.0)

          r32 = ppiclf_rprop(PPICLF_R_JRR,i)

          r32hat = r32/r0
          rphat = 1
          rchat = r_cutoff / r0

          ! log normal distribution conversion count median non-dim radius
          CMR = r32hat / EXP(2.5*LOG(sigma_g)**2)
          LCMR = LOG(CMR)
          ! mass mean non-dim radius corresponds to 3rd
          ! moment distribution and can be used to
          ! find cumulative mass below a value of rhat
          rmm = CMR * EXP(3.5*LOG(sigma_g)**2)

          ! radius of average mass
          rmbar = r0*CMR*EXP(1.5*LOG(sigma_g)**2)

          ! find mass below cutoff
          call cdfnor(1, F_cutoff, Q_chi, LOG(rchat), LOG(rmm),
     >                         Lsigma_g, status_chi, bound_chi)

          call cdfnor(1, F_parent, Q_chi, LOG(rphat), LOG(rmm),
     >                         Lsigma_g, status_chi, bound_chi)

          V_part= ppiclf_rprop(PPICLF_R_JVOLP,i)*(1-F_cutoff/F_parent)

          ! find F values for cutoff and parent
          call cdfnor(1, F_cutoff, Q_chi, LOG(rchat), LOG(CMR),
     >                           Lsigma_g, status_chi, bound_chi)

          ! this works better for cases where most of the volume is
          ! below the cutoff
c         F_cutoff = 0

         call cdfnor(1, F_parent, Q_chi, LOG(rphat), LOG(CMR),
     >                         Lsigma_g, status_chi, bound_chi)

          vol_children = 0.0
          n_child = IDINT((r0/r32)**3)  ! computational particle

          n_child_i = 0
            do j=1,n_child
              if ((ppiclf_npart + n_added) .gt. PPICLF_LPART) then
                error stop "The TAB model tried to add too many
     >            particles, try increasing PPICLF_LPART" 
              end if

              P_chi = ppiclf_ran(2)*(F_parent - F_cutoff) + F_cutoff
              Q_chi = 1 - P_chi

              call cdfnor(2, P_chi, Q_chi, r_child, LOG(CMR),
     >                       Lsigma_g, status_chi, bound_chi)
              r_child = EXP(r_child)*r0

              if (r_child .gt. r_cutoff) then
                n_added   = n_added + 1   ! physical particle
                n_child_i = n_child_i + 1

                vol_child    = 4.0/3.0*rpi*r_child**(3.0)
                vol_children = vol_children + vol_child

                do iprop=1,PPICLF_LRS
                add_y(iprop,n_added) = ppiclf_y(iprop,i)
                end do
                
                do iprop=1,PPICLF_LRP
                add_rprop(iprop, n_added) = ppiclf_rprop(iprop, i)
                end do
                add_rprop(PPICLF_R_JDMAJ   ,n_added) = 2.0*r_child
                add_rprop(PPICLF_R_JVOLP   ,n_added) = vol_child
                add_rprop(PPICLF_R_JR0TAB  ,n_added) = r_child
                add_rprop(PPICLF_R_JYTAB   ,n_added) = 0.0
                add_rprop(PPICLF_R_JYDOTTAB,n_added) = 0.0
                add_rprop(PPICLF_R_JSPL    ,n_added) = 1.0
                add_rprop(PPICLF_R_JFTAB   ,n_added) = 0.0
           
                ! randomly select angle going from nvec to mvec
                v_child = cb_tab*cv_tab*r0*
     >                    ppiclf_rprop(PPICLF_R_JYDOTTAB,i)

                theta_child = 2.0*rpi*ppiclf_ran(2)!theta_child

                ! radial velocity from TAB model
                add_y(PPICLF_JVX,n_added) = add_y(PPICLF_JVX,n_added) +
     >                                 v_child*nvec1*COS(theta_child) +
     >                                 v_child*mvec1*SIN(theta_child)

                add_y(PPICLF_JVY,n_added) = add_y(PPICLF_JVY,n_added) +
     >                                 v_child*nvec2*COS(theta_child) +
     >                                 v_child*mvec2*SIN(theta_child)
                
                ! nvec3=0 by definition
                add_y(PPICLF_JVZ,j) = add_y(PPICLF_JVZ,j) +
     >                                v_child*mvec3*SIN(theta_child)

                ! generate unique ID using Hopcroft and Ullman pairing function
                !-------------------------------
                ! Don't trust this for anything larger than very small
                ! test cases because it tends to overflow
                !-------------------------------
                ! first combine tindex, ppiclf_nid
                uid = 0.5*(tindex+ppiclf_nid-2)*(tindex+ppiclf_nid-1) 
     >                + tindex
                ! second combine the previous result with i
                uid = 0.5*(uid+i-2)*(uid+i-1) + uid
                ! third combine with j
                uid = 0.5*(uid+j-2)*(uid+j-1) + uid
                add_rprop(PPICLF_R_JUID,n_added) = uid
              endif ! r_cutoff
            enddo  ! n-child

           ! correct SPL to conserve mass
            !ppiclf_rprop(PPICLF_R_JSPL,i) = 0.9
            !ppiclf_rprop(PPICLF_R_JPARENT,i) = 1
            C_spl = V_part*ppiclf_rprop(PPICLF_R_JSPL,i) / vol_children

            do j=1,n_child_i
              k = n_added - n_child_i + j
              add_rprop(PPICLF_R_JSPL,k) =
     >          C_spl*add_rprop(PPICLF_R_JSPL,k)

              do ilog=1,PPICLF_LRS
                log_spawn(ilog, k) = add_y(ilog, k)
              enddo
              log_spawn(PPICLF_LRS+1, k) = add_rprop(PPICLF_R_JR0TAB,k)
              log_spawn(PPICLF_LRS+2, k) = add_rprop(PPICLF_R_JSPL,k)
              log_spawn(PPICLF_LRS+3, k) = ppiclf_time
              log_spawn(PPICLF_LRS+4, k) = add_rprop(PPICLF_R_JUID,k)
              log_spawn(PPICLF_LRS+5, k) = ppiclf_rprop(PPICLF_R_JUID,i)
              log_spawn(PPICLF_LRS+6, k) = add_rprop(PPICLF_R_JPARENT,k)
            enddo

            call ppiclf_solve_MarkForRemoval(i)

          endif ! y >= y_cr
c          !call ppiclf_solve_MarkForRemoval(i)
c          !call ppiclf_solve_RemoveParticle(i)

      enddo ! ppiclf_npart

c      call ppiclf_solve_AddParticles(n_added, add_y, add_rprop)   

! End Droplet Breakup Model
!########################################################################
!########################################################################
      do i=1,ppiclf_npart

         dp_maj= ppiclf_rprop(PPICLF_R_JDMAJ, i)
         r0    = ppiclf_rprop(PPICLF_R_JR0TAB,i)
         vmag  = sqrt(
     >   (ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i))**2
     >  +(ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i))**2
     >  +(ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i))**2)

c         ppiclf_rprop(PPICLF_R_JWeb,i) = rhof*vmag**2*r0/sigma
         rep = rhof*vmag*2*r0/rmu
         CD_sph = 0.36+5.48*rep**(-0.573) + 24/rep
         fqs = 0.5*CD_sph*rpi*(2*r0)**2.0/4.0*rhof*vmag
         !fqs = 3.d0*rpi*rmu*r0*(1.d0+0.15d0*rep**0.687d0)

         fqsx  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUX,i)-ppiclf_y(PPICLF_JVX,i))
         fqsy  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUY,i)-ppiclf_y(PPICLF_JVY,i))
         fqsz  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUZ,i)-ppiclf_y(PPICLF_JVZ,i))
         
         ! Deformation drag
         r_minor = r0**3.0 / (dp_maj/2.0)**2.0
         if ((2*r_minor .ne. dp_maj) .and. (vmag .gt. 0.0)) then
           e_corr = 1.0 - r_minor**2.0/(dp_maj/2.0)**2.0
           A_part = 2.0*rpi*(dp_maj/2.0)**2.0 + 
     >           rpi*r_minor**2.0/e_corr*log((1.0+e_corr)/(1.0-e_corr))
           D32 = 6.0*ppiclf_rprop(PPICLF_R_JVOLP,i)/A_part

           Re32 = rhof*vmag*D32 / rmu
           if (Re32 .gt. 0.2) then
           CD_disk = 1.1 + 64.0 / (rpi*Re32)
           f_interp = 1.0 - (1.0/(dp_maj/D32)**3.0)**2.0

           f_diskx = 0.5*CD_disk*rpi*(2*r0)**2.0/4.0*rhof*vmag
     >       *(ppiclf_rprop(PPICLF_R_JUX,i)-ppiclf_y(PPICLF_JVX,i))

           f_disky = 0.5*CD_disk*rpi*(2*r0)**2.0/4.0*rhof*vmag
     >       *(ppiclf_rprop(PPICLF_R_JUY,i)-ppiclf_y(PPICLF_JVY,i))

           f_diskz = 0.5*CD_disk*rpi*(2*r0)**2.0/4.0*rhof*vmag
     >       *(ppiclf_rprop(PPICLF_R_JUZ,i)-ppiclf_y(PPICLF_JVZ,i))

           fqsx = f_interp*f_diskx + (1.0 - f_interp)*fqsx
           fqsy = f_interp*f_disky + (1.0 - f_interp)*fqsy
           fqsz = f_interp*f_diskz + (1.0 - f_interp)*fqsz
           endif
         endif

         ! gravity
         rmass =  ppiclf_rprop(PPICLF_R_JVOLP,i) 
     >           *ppiclf_rprop(PPICLF_R_JRHOP,i)
         dmass =  ppiclf_rprop(PPICLF_R_JVOLP,i) 
     >          *(ppiclf_rprop(PPICLF_R_JRHOP,i)-rhof)

         fbx  = 0.d0
         fby  = dmass*rg
         fbz  = 0.d0

         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ, i) = ppiclf_y(PPICLF_JVZ,i)

         ppiclf_ydot(PPICLF_JVX,i) = (fqsx+fbx)/rmass
         ppiclf_ydot(PPICLF_JVY,i) = (fqsy+fby)/rmass
         ppiclf_ydot(PPICLF_JVZ,i) = (fqsz+fbz)/rmass

      enddo      !ppiclf_npart
!########################################################################
! end evaluate ydot
!########################################################################
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
!
! Input:
!
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
!
! Output:
!
      real*8 map  (PPICLF_LRP_PRO)
!
! Internal:
!
      real*8 dp_norm

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i
      integer*4 j
      real*8 yi    (PPICLF_LRS)
      real*8 rpropi(PPICLF_LRP)
      real*8 yj    (PPICLF_LRS)
      real*8 rpropj(PPICLF_LRP)
!
! Internal:
!
      real*8 ksp,erest
      common /ucollision/ ksp,erest
#ifdef PPICLC
      BIND(C, name="ucollision") :: /ucollision/ ! c binding
#endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SaveRemoved(i)
!
      implicit none
!
      include "PPICLF"
!
      real*8 y_out(PPICLF_LRS         , PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP + 1 , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART)          ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out

! Input:
!
      integer*4 i,j
      real*8 rp,rpi,alpha
!
      npart_out = npart_out + 1

!     SAVE PARTICLE POSITION WHEN OUT OF DOMIAN
      do j=1,PPICLF_LRS
        y_out(j,npart_out) = ppiclf_y(j,i)
      enddo

      do j=1,PPICLF_LRP
        rprop_out(j,npart_out) = ppiclf_rprop(j,i)
      enddo

      rprop_out(PPICLF_LRP+1,npart_out) = ppiclf_time
      
      do j=1,3
        iprop_out(j,npart_out) = ppiclf_iprop(j+4,i)
      enddo
      
      iprop_out(4,npart_out) = ppiclf_cycle
!
      iprop_out(5,npart_out) = 13

      return
      end

!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteParticleOutletVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp((PPICLF_LRP+1)*PPICLF_LPART)
     >       ,rout_lip(5      *PPICLF_LPART)
      character*6 filein
      character*15 vtufile
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth,prevs(2,ppiclf_np)
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 iint, nnp, nxx, npt_total, jx, jy, jz, if_sz, isize,
     >          iadd, if_pos, if_sln, if_lrp, if_lip, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, i, j, ie, nps, nglob, nkey, ndum,
     >          icount_pos, icount_sln, icount_lrp, icount_lip, iorank,
     >          ierr, ivtu_size
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
!
! External:
!
      real*8 y_out(PPICLF_LRS           , PPICLF_LPART) ! Normal ordering
      real*8 rprop_out(PPICLF_LRP+1     , PPICLF_LPART) ! Normal ordering
      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
      integer*4 npart_out
      common /outlet/ y_out,rprop_out,iprop_out,npart_out
!
!
!
      call ppiclf_printsi(' *Begin WriteParticleOutletVTU$'
     >                   ,ppiclf_cycle)

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = npart_out

      npt_total = ppiclf_iglsum(nxx,1)

      jx    = 1
      jy    = 2
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

c      if_sz = len(filein1)
c      if (if_sz .lt. 3) then
c         filein = 'par'
c      else 
         write(filein,'(A6)') filein1
c      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize = 4

      iadd = 0
      if_pos = 3*isize*npt_total
      if_sln = 1*isize*npt_total
      if_lrp = 1*isize*npt_total
      if_lip = 1*isize*npt_total

      ic_pos = iadd
      ic_sln = iadd
      ic_lrp = iadd
      ic_lip = iadd
 
      if (nxx.ne.0) then
      do i=1,nxx
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(y_out(jx,i))
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(y_out(jy,i))
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
            rout_pos(ic_pos) = sngl(y_out(jz,i))
         else
            rout_pos(ic_pos) = 0.0
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,nxx
         ic_sln = ic_sln + 1
         rout_sln(ic_sln) = sngl(y_out(j,i))
      enddo
      enddo
      do j=1,PPICLF_LRP+1
      do i=1,nxx
         ic_lrp = ic_lrp + 1
         rout_lrp(ic_lrp) = sngl(rprop_out(j,i))
      enddo
      enddo
      do j=1,5
      do i=1,nxx
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = real(iprop_out(j,i))
      enddo
      enddo
      endif ! nxx.ne.0

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,nnp
         prevs(1,i) = i-1
         prevs(2,i) = nxx
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call pfgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call pfgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nglob,nkey)

      stride_len = 0
      if (ppiclf_nid .ne. 0) then
      do i=1,ppiclf_nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A6,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then !---------------------------

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') npt_total
      write(vtu,'(A)',advance='yes') '" NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'


      do ie=1,PPICLF_LRS
         write(prostr,'(A1,I2.2)') "y",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,PPICLF_LRP+1
         write(prostr,'(A4,I2.2)') "rprop",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,5
         write(prostr,'(A3,I2.2)') "tag",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      write(vtu,'(A)',advance='yes') '   </PointData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif ! ------------ nid .eq. 0 -------------------------

      call ppiclf_bcast(ivtu_size, isize)

      ! byte-displacements
      idisp_pos = ivtu_size + isize*(3*stride_len + 1)

      ! how much to write
      icount_pos = 3*nxx
      icount_sln = 1*nxx
      icount_lrp = 1*nxx
      icount_lip = 1*nxx

      iorank = -1

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
        call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
        call ppiclf_byte_set_view(idisp_pos,pth)
        call ppiclf_byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
        call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_sln
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*npart_out + 1
   
         ! write
c         if (nxx.ne.0) then
           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
           call ppiclf_byte_set_view(idisp_sln,pth)
           call ppiclf_byte_write_mpi(rout_sln(j),icount_sln,iorank,pth
     >                             ,ierr)
           call ppiclf_byte_close_mpi(pth,ierr)
c        endif ! nxx.ne.0
      enddo

      do i=1,PPICLF_LRP+1
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lrp
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*npart_out + 1
   
         ! write
c        if (nxx.ne.0) then
           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
           call ppiclf_byte_set_view(idisp_lrp,pth)
           call ppiclf_byte_write_mpi(rout_lrp(j),icount_lrp,iorank,pth
     >                             ,ierr)
           call ppiclf_byte_close_mpi(pth,ierr)
c         endif
      enddo

      do i=1,5
         idisp_lip = ivtu_size + isize*(3*npt_total
     >                         + PPICLF_LRS*npt_total
     >                         + (PPICLF_LRP+1)*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + 1 + i)
         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lip
           close(vtu)
         endif

         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*npart_out + 1
   
         ! write
c         if (nxx.ne.0) then
           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
           call ppiclf_byte_set_view(idisp_lip,pth)
           call ppiclf_byte_write_mpi(rout_lip(j),icount_lip,iorank,pth
     >                             ,ierr)
           call ppiclf_byte_close_mpi(pth,ierr)
c         endif
      enddo

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi(' *End WriteParticleOutletVTU$',ppiclf_cycle)

      return
      end

!-----------------------------------------------------------------------
      REAL*8 FUNCTION ppiclf_ran(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate
c between 0.0 and 1.0 (exclusive of the endpoint values).
c Call with idum a negative integer to initialize; thereafter, do not alter
c idum between successive deviates in a sequence. RNMX should approximate the
c largest floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum1=max(-idum,1)
         idum2=idum1
         do j=NTAB+8,1,-1
            k=idum1/IQ1
            idum1=IA1*(idum1-k*IQ1)-k*IR1
            if (idum1.lt.0) idum1=idum1+IM1
            if (j.le.NTAB) iv(j)=idum1
         enddo
         iy=iv(1)
      endif
      k=idum1/IQ1
      idum1=IA1*(idum1-k*IQ1)-k*IR1
      if (idum1.lt.0) idum1=idum1+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum1
      if(iy.lt.1)iy=iy+IMM1
      ppiclf_ran=min(AM*iy,RNMX)
      return
      END

c-----------------------------------------------------------------------
      real*8 function fun_erfinv(x)

      implicit none
      real*8 x,a
      real*8 :: pi=4.*atan(1.d0)

      if (x.eq.0.d0)then
         fun_erfinv=0.d0
      else
         a=8.0*(pi-3.0)/(3.0*pi*(4.0-pi))
         fun_erfinv=x/abs(x)*sqrt(sqrt((2.0/pi/a+0.5*log(1.0-x**2))**2
     $             -log(1.0-x**2)/a)-(2.0/pi/a+0.5*log(1.0-x**2)))
      endif

      return
      end

c**********************************************************************      
      SUBROUTINE cdfnor(which,p,q,x,mean,sd,status,bound)
!
C     Cumulative Distribution Function Normal distribution
C     Calculates any one parameter of the normal
C     distribution given values for the others.
C
C     .. Parameters ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,mean,p,q,sd,x
      INTEGER status,which
C     .. Local Scalars ..
      DOUBLE PRECISION z,pq
C     .. External Functions ..
      DOUBLE PRECISION dinvnr,spmpar
      EXTERNAL dinvnr,spmpar
C     .. External Subroutines ..
      EXTERNAL cumnor
C     .. Executable Statements ..
C     Check arguments
C
      status = 0
      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 4.0D0
   20 status = -1
      RETURN

   30 IF (which.EQ.1) GO TO 70
C
C     P
C
      IF (.NOT. ((p.LE.0.0D0).OR. (p.GT.1.0D0))) GO TO 60
      IF (.NOT. (p.LE.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = 1.0D0
   50 status = -2
      RETURN

   60 CONTINUE
   70 IF (which.EQ.1) GO TO 110
C
C     Q
C
      IF (.NOT. ((q.LE.0.0D0).OR. (q.GT.1.0D0))) GO TO 100
      IF (.NOT. (q.LE.0.0D0)) GO TO 80
      bound = 0.0D0
      GO TO 90

   80 bound = 1.0D0
   90 status = -3
      RETURN

  100 CONTINUE
  110 IF (which.EQ.1) GO TO 150
C
C     P + Q
C
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 140
      IF (.NOT. (pq.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      GO TO 130

  120 bound = 1.0D0
  130 status = 3
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 170
C
C     SD
C
      IF (.NOT. (sd.LE.0.0D0)) GO TO 160
      bound = 0.0D0
      status = -6
      RETURN

  160 CONTINUE
C
C     Calculate ANSWERS
C
  170 IF ((1).EQ. (which)) THEN
C
C     Computing P
C
          z = (x-mean)/sd
          CALL cumnor(z,p,q)

      ELSE IF ((2).EQ. (which)) THEN
C
C     Computing X
C
          z = dinvnr(p,q)
          x = sd*z + mean

      ELSE IF ((3).EQ. (which)) THEN
C
C     Computing the MEAN
C
          z = dinvnr(p,q)
          mean = x - sd*z

      ELSE IF ((4).EQ. (which)) THEN
C
C     Computing SD
C
          z = dinvnr(p,q)
          sd = (x-mean)/z
      END IF

      RETURN

      END
C-----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION spmpar(i)
C
C     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
C     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
C
C        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
C        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
C        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER i
C     .. Local Scalars ..
      DOUBLE PRECISION b,binv,bm1,one,w,z
      INTEGER emax,emin,ibeta,m
C     .. External Functions ..
      INTEGER ipmpar
      EXTERNAL ipmpar
C     .. Intrinsic Functions ..
      INTRINSIC dble
C     .. Executable Statements ..
C
      IF (i.GT.1) GO TO 10
      b = ipmpar(4)
      m = ipmpar(8)
      spmpar = b** (1-m)
      RETURN
C
   10 IF (i.GT.2) GO TO 20
      b = ipmpar(4)
      emin = ipmpar(9)
      one = dble(1)
      binv = one/b
      w = b** (emin+2)
      spmpar = ((w*binv)*binv)*binv
      RETURN
C
   20 ibeta = ipmpar(4)
      m = ipmpar(8)
      emax = ipmpar(10)
C
      b = ibeta
      bm1 = ibeta - 1
      one = dble(1)
      z = b** (m-1)
      w = ((z-one)*b+bm1)/ (b*z)
C
      z = b** (emax-2)
      spmpar = ((w*z)*b)*b

      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION dinvnr(p,q)
C
C                              Function
C     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
C     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
C
C     P --> The probability whose normal deviate is sought.
C                    P is DOUBLE PRECISION
C     Q --> 1-P
C                    P is DOUBLE PRECISION
C
C**********************************************************************
C     .. Parameters ..
      INTEGER maxit
      PARAMETER (maxit=100)
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0D-13)
      DOUBLE PRECISION r2pi
      PARAMETER (r2pi=0.3989422804014326D0)
      DOUBLE PRECISION nhalf
      PARAMETER (nhalf=-0.5D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION p,q
C     .. Local Scalars ..
      DOUBLE PRECISION strtx,xcur,cum,ccum,pp,dx
      INTEGER i
      LOGICAL qporq
C     .. External Functions ..
      DOUBLE PRECISION stvaln
      EXTERNAL stvaln
C     .. External Subroutines ..
      EXTERNAL cumnor
C     .. Statement Functions ..
      DOUBLE PRECISION dennor,x

      dennor(x) = r2pi*exp(nhalf*x*x)
C     FIND MINIMUM OF P AND Q
C
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 10
      pp = p
      GO TO 20

   10 pp = q
C
C     INITIALIZATION STEP
C
   20 strtx = stvaln(pp)
      xcur = strtx
C
C     NEWTON INTERATIONS
C
      DO 30,i = 1,maxit
          CALL cumnor(xcur,cum,ccum)
          dx = (cum-pp)/dennor(xcur)
          xcur = xcur - dx
          IF (abs(dx/xcur).LT.eps) GO TO 40
   30 CONTINUE
      dinvnr = strtx
C
C     IF WE GET HERE, NEWTON HAS FAILED
C
      IF (.NOT.qporq) dinvnr = -dinvnr
      RETURN
C
C     IF WE GET HERE, NEWTON HAS SUCCEDED
C
   40 dinvnr = xcur
      IF (.NOT.qporq) dinvnr = -dinvnr

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE cumnor(arg,result,ccum)

      INTEGER i
      DOUBLE PRECISION a,arg,b,c,d,del,eps,half,p,one,q,result,sixten,
     +                 temp,sqrpi,thrsh,root32,x,xden,xnum,y,xsq,zero,
     +                 min,ccum
      DIMENSION a(5),b(4),c(9),d(8),p(6),q(5)
C------------------------------------------------------------------
C  External Function
C------------------------------------------------------------------
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C------------------------------------------------------------------
C  Mathematical constants
C
C  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
C  THRSH is the argument for which anorm = 0.75.
C------------------------------------------------------------------
      DATA one,half,zero,sixten/1.0D0,0.5D0,0.0D0,1.60D0/,
     +     sqrpi/3.9894228040143267794D-1/,thrsh/0.66291D0/,
     +     root32/5.656854248D0/
C------------------------------------------------------------------
C  Coefficients for approximation in first interval
C------------------------------------------------------------------
      DATA a/2.2352520354606839287D00,1.6102823106855587881D02,
     +     1.0676894854603709582D03,1.8154981253343561249D04,
     +     6.5682337918207449113D-2/
      DATA b/4.7202581904688241870D01,9.7609855173777669322D02,
     +     1.0260932208618978205D04,4.5507789335026729956D04/
C------------------------------------------------------------------
C  Coefficients for approximation in second interval
C------------------------------------------------------------------
      DATA c/3.9894151208813466764D-1,8.8831497943883759412D00,
     +     9.3506656132177855979D01,5.9727027639480026226D02,
     +     2.4945375852903726711D03,6.8481904505362823326D03,
     +     1.1602651437647350124D04,9.8427148383839780218D03,
     +     1.0765576773720192317D-8/
      DATA d/2.2266688044328115691D01,2.3538790178262499861D02,
     +     1.5193775994075548050D03,6.4855582982667607550D03,
     +     1.8615571640885098091D04,3.4900952721145977266D04,
     +     3.8912003286093271411D04,1.9685429676859990727D04/
C------------------------------------------------------------------
C  Coefficients for approximation in third interval
C------------------------------------------------------------------
      DATA p/2.1589853405795699D-1,1.274011611602473639D-1,
     +     2.2235277870649807D-2,1.421619193227893466D-3,
     +     2.9112874951168792D-5,2.307344176494017303D-2/
      DATA q/1.28426009614491121D00,4.68238212480865118D-1,
     +     6.59881378689285515D-2,3.78239633202758244D-3,
     +     7.29751555083966205D-5/
C------------------------------------------------------------------
C  Machine dependent constants
C------------------------------------------------------------------
      eps = spmpar(1)*0.5D0
      min = spmpar(2)
C------------------------------------------------------------------
      x = arg
      y = abs(x)
      IF (y.LE.thrsh) THEN
C------------------------------------------------------------------
C  Evaluate  anorm  for  |X| <= 0.66291
C------------------------------------------------------------------
          xsq = zero
          IF (y.GT.eps) xsq = x*x
          xnum = a(5)*xsq
          xden = xsq
          DO 10 i = 1,3
              xnum = (xnum+a(i))*xsq
              xden = (xden+b(i))*xsq
   10     CONTINUE
          result = x* (xnum+a(4))/ (xden+b(4))
          temp = result
          result = half + temp
          ccum = half - temp
C------------------------------------------------------------------
C  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
C------------------------------------------------------------------
      ELSE IF (y.LE.root32) THEN
          xnum = c(9)*y
          xden = y
          DO 20 i = 1,7
              xnum = (xnum+c(i))*y
              xden = (xden+d(i))*y
   20     CONTINUE
          result = (xnum+c(8))/ (xden+d(8))
          xsq = aint(y*sixten)/sixten
          del = (y-xsq)* (y+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF
C------------------------------------------------------------------
C  Evaluate  anorm  for |X| > sqrt(32)
C------------------------------------------------------------------
      ELSE
          result = zero
          xsq = one/ (x*x)
          xnum = p(6)*xsq
          xden = xsq
          DO 30 i = 1,4
              xnum = (xnum+p(i))*xsq
              xden = (xden+q(i))*xsq
   30     CONTINUE
          result = xsq* (xnum+p(5))/ (xden+q(5))
          result = (sqrpi-result)/y
          xsq = aint(x*sixten)/sixten
          del = (x-xsq)* (x+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF

      END IF

      IF (result.LT.min) result = 0.0D0
      IF (ccum.LT.min) ccum = 0.0D0

      RETURN
      END
**********************************************************************
      DOUBLE PRECISION FUNCTION stvaln(p)

C     .. Scalar Arguments ..
      DOUBLE PRECISION p
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION sign,y,z
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION xden(5),xnum(5)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl
      EXTERNAL devlpl
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,log,sqrt
C     ..
C     .. Data statements ..
      DATA xnum/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,
     +     -0.204231210245D-1,-0.453642210148D-4/
      DATA xden/0.993484626060D-1,0.588581570495D0,0.531103462366D0,
     +     0.103537752850D0,0.38560700634D-2/
C     ..
C     .. Executable Statements ..
      IF (.NOT. (p.LE.0.5D0)) GO TO 10
      sign = -1.0D0
      z = p
      GO TO 20

   10 sign = 1.0D0
      z = 1.0D0 - p
   20 y = sqrt(-2.0D0*log(z))
      stvaln = y + devlpl(xnum,5,y)/devlpl(xden,5,y)
      stvaln = sign*stvaln
      RETURN

      END      

      INTEGER FUNCTION ipmpar(i)
C
C     .. Scalar Arguments ..
      INTEGER i
C     ..
C     .. Local Arrays ..
      INTEGER imach(10)
C     ..
      DATA IMACH( 1) /     2 /
      DATA IMACH( 2) /    31 /
      DATA IMACH( 3) / 2147483647 /
      DATA IMACH( 4) /     2 /
      DATA IMACH( 5) /    24 /
      DATA IMACH( 6) /  -125 /
      DATA IMACH( 7) /   128 /
      DATA IMACH( 8) /    53 /
      DATA IMACH( 9) / -1021 /
      DATA IMACH(10) /  1024 /
C
C
      ipmpar = imach(i)

      RETURN
      END
**********************************************************************
      DOUBLE PRECISION FUNCTION devlpl(a,n,x)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION a(n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION term
      INTEGER i
C     ..
C     .. Executable Statements ..
      term = a(n)
      DO 10,i = n - 1,1,-1
          term = a(i) + term*x
   10 CONTINUE
      devlpl = term

      RETURN
      END
      
