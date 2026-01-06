!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rpi,rp,vmag,rep,rmass,dmass,fbx,fby,fbz,fqs,fqsx,fqsy,fqsz
      real*8 rr
      integer*4 i
!
! External:
!
      real*8 rmu,rhof,rg
      common /parameters/ rmu,rhof,rg

      real*8 dpmin
      common /dpm/ dpmin

      real*8 rlx,rrx,rly,rry,rlz,rrz
      common /domainsize/ rlx,rrx,rly,rry,rlz,rrz

      real*8 ppiclf_ran

      rpi  = 4.0*atan(1.0)

! evaluate ydot
!-----------------------------------------------------------------------
      do i=1,ppiclf_npart

         if (ppiclf_rprop(PPICLF_R_JDP,i).le.dpmin) then ! Small particles

          ppiclf_ydot(PPICLF_JX ,i) = ppiclf_rprop(PPICLF_R_JUX,i)
          ppiclf_ydot(PPICLF_JY ,i) =(ppiclf_rprop(PPICLF_R_JUY,i)
     >                              + ppiclf_rprop(PPICLF_R_JVS,i))
          ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_rprop(PPICLF_R_JUZ,i)
  
          ppiclf_ydot(PPICLF_JVX,i) = 0.0
          ppiclf_ydot(PPICLF_JVY,i) = 0.0
          ppiclf_ydot(PPICLF_JVZ,i) = 0.0
         else
c        Particle motion with force models
         vmag  = sqrt((ppiclf_rprop(PPICLF_R_JUX,i)
     >                -ppiclf_y(PPICLF_JVX,i))**2
     >               +(ppiclf_rprop(PPICLF_R_JUY,i)
     >                -ppiclf_y(PPICLF_JVY,i))**2
     >               +(ppiclf_rprop(PPICLF_R_JUZ,i)
     >                -ppiclf_y(PPICLF_JVZ,i))**2)

         ppiclf_rprop(PPICLF_R_JREP,i) = 
     >           rhof*vmag*ppiclf_rprop(PPICLF_R_JDP,i)/rmu
         rep = ppiclf_rprop(PPICLF_R_JREP,i)

         ! gravity
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)
         dmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *(ppiclf_rprop(PPICLF_R_JRHOP,i)-rhof)


         fbx  = 0.d0
         fby  = rmass*rg
         fbz  = 0.d0

         ! quasi-steady
         ppiclf_rprop(PPICLF_R_JFQS,i) = 3.d0*rpi*rmu
     >   *ppiclf_rprop(PPICLF_R_JDP,i)*(1.d0+0.15d0*(rep**0.687d0))
         fqs = ppiclf_rprop(PPICLF_R_JFQS,i)

         fqsx  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUX,i)-ppiclf_y(PPICLF_JVX,i))
         fqsy  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUY,i)-ppiclf_y(PPICLF_JVY,i))
         fqsz  = fqs
     >   *(ppiclf_rprop(PPICLF_R_JUZ,i)-ppiclf_y(PPICLF_JVZ,i))


         ! set ydot for all PPICLF_SLN number of equations
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_y(PPICLF_JVZ,i)

         ppiclf_ydot(PPICLF_JVX,i) = (fbx+fqsx)/rmass
         ppiclf_ydot(PPICLF_JVY,i) = (fby+fqsy)/rmass
         ppiclf_ydot(PPICLF_JVZ,i) = (fbz+fqsz)/rmass

         
         endif

      enddo 
! evaluate ydot


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
!

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
      
c!-----------------------------------------------------------------------
c      subroutine ppiclf_solve_SaveRemoved(i)
c!
c      implicit none
c!
c      include "PPICLF"
c!
c      real*8 y_out(PPICLF_LRS         , PPICLF_LPART) ! Normal ordering
c      real*8 rprop_out(PPICLF_LRP + 1 , PPICLF_LPART) ! Normal ordering
c      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
c      integer*4 npart_out
c      common /outlet/ y_out,rprop_out,iprop_out,npart_out
c
c
c
c! Input:
c!
c      integer*4 i,j
c      real*8 rp
c!
c      npart_out = npart_out + 1
c
c!     SAVE PARTICLE POSITION WHEN OUT OF DOMIAN
c      do j=1,PPICLF_LRS
c        y_out(j,npart_out) = ppiclf_y(j,i)
c      enddo
c
c      do j=1,PPICLF_LRP
c        rprop_out(j,npart_out) = ppiclf_rprop(j,i)
c      enddo
c
c      rprop_out(PPICLF_LRP+1,npart_out) = ppiclf_time
c      
c      do j=1,3
c        iprop_out(j,npart_out) = ppiclf_iprop(j+4,i)
c      enddo
c      
c      iprop_out(4,npart_out) = ppiclf_cycle
c!     
c
c      iprop_out(5,npart_out) = 13
c
c      return
c      end
c
!-----------------------------------------------------------------------
c      subroutine ppiclf_io_WriteParticleOutletVTU(filein1)
c!
c      implicit none
c!
c      include "PPICLF"
c      include 'mpif.h'
c!
c! Input:
c!
c      character (len = *) filein1
c!
c! Internal:
c!
c      real*4  rout_pos(3      *PPICLF_LPART) 
c     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
c     >       ,rout_lrp((PPICLF_LRP+1)*PPICLF_LPART)
c     >       ,rout_lip(5      *PPICLF_LPART)
c      character*6 filein
c      character*15 vtufile
c      character*6  prostr
c      integer*4 icalld1
c      save      icalld1
c      data      icalld1 /0/
c      integer*4 vtu,pth,prevs(2,ppiclf_np)
c      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
c      integer*4 iint, nnp, nxx, npt_total, jx, jy, jz, if_sz, isize,
c     >          iadd, if_pos, if_sln, if_lrp, if_lip, ic_pos, ic_sln,
c     >          ic_lrp, ic_lip, i, j, ie, nps, nglob, nkey, ndum,
c     >          icount_pos, icount_sln, icount_lrp, icount_lip, iorank,
c     >          ierr, ivtu_size
c      integer*4 ppiclf_iglsum
c      external ppiclf_iglsum
c!
c! External:
c!
c      real*8 y_out(PPICLF_LRS           , PPICLF_LPART) ! Normal ordering
c      real*8 rprop_out(PPICLF_LRP+1     , PPICLF_LPART) ! Normal ordering
c      integer*4 iprop_out(5  , PPICLF_LPART) ! Normal ordering
c      integer*4 npart_out
c      common /outlet/ y_out,rprop_out,iprop_out,npart_out
c!
c!
c!
c      call ppiclf_printsi(' *Begin WriteParticleOutletVTU$'
c     >                   ,ppiclf_cycle)
c
c      icalld1 = icalld1+1
c
c      nnp   = ppiclf_np
c      nxx   = npart_out
c
c      npt_total = ppiclf_iglsum(nxx,1)
c
c      jx    = 1
c      jy    = 2
c      jz    = 1
c      if (ppiclf_ndim .eq. 3)
c     >jz    = 3
c
cc      if_sz = len(filein1)
cc      if (if_sz .lt. 3) then
cc         filein = 'par'
cc      else 
c         write(filein,'(A6)') filein1
cc      endif
c
c! --------------------------------------------------
c! COPY PARTICLES TO OUTPUT ARRAY
c! --------------------------------------------------
c
c      isize = 4
c
c      iadd = 0
c      if_pos = 3*isize*npt_total
c      if_sln = 1*isize*npt_total
c      if_lrp = 1*isize*npt_total
c      if_lip = 1*isize*npt_total
c
c      ic_pos = iadd
c      ic_sln = iadd
c      ic_lrp = iadd
c      ic_lip = iadd
c 
c      if (nxx.ne.0) then
c      do i=1,nxx
c         ic_pos = ic_pos + 1
c         rout_pos(ic_pos) = sngl(y_out(jx,i))
c         ic_pos = ic_pos + 1
c         rout_pos(ic_pos) = sngl(y_out(jy,i))
c         ic_pos = ic_pos + 1
c         if (ppiclf_ndim .eq. 3) then
c            rout_pos(ic_pos) = sngl(y_out(jz,i))
c         else
c            rout_pos(ic_pos) = 0.0
c         endif
c      enddo
c      do j=1,PPICLF_LRS
c      do i=1,nxx
c         ic_sln = ic_sln + 1
c         rout_sln(ic_sln) = sngl(y_out(j,i))
c      enddo
c      enddo
c      do j=1,PPICLF_LRP+1
c      do i=1,nxx
c         ic_lrp = ic_lrp + 1
c         rout_lrp(ic_lrp) = sngl(rprop_out(j,i))
c      enddo
c      enddo
c      do j=1,5
c      do i=1,nxx
c         ic_lip = ic_lip + 1
c         rout_lip(ic_lip) = real(iprop_out(j,i))
c      enddo
c      enddo
c      endif ! nxx.ne.0
c
c! --------------------------------------------------
c! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
c! --------------------------------------------------
c      do i=1,nnp
c         prevs(1,i) = i-1
c         prevs(2,i) = nxx
c      enddo
c
c      nps   = 1 ! index of new proc for doing stuff
c      nglob = 1 ! unique key to sort by
c      nkey  = 1 ! number of keys (just 1 here)
c      ndum = 2
c      call pfgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
c     >                 ndum,nnp,nnp,nps)
c      call pfgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
c     >                 ndum,nnp,nglob,nkey)
c
c      stride_len = 0
c      if (ppiclf_nid .ne. 0) then
c      do i=1,ppiclf_nid
c         stride_len = stride_len + prevs(2,i)
c      enddo
c      endif
c
c! ----------------------------------------------------
c! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
c! ----------------------------------------------------
c      write(vtufile,'(A6,I5.5,A4)') filein,icalld1,'.vtu'
c
c      if (ppiclf_nid .eq. 0) then !---------------------------
c
c      vtu=867+ppiclf_nid
c      open(unit=vtu,file=vtufile,status='replace')
c
c! ------------
c! FRONT MATTER
c! ------------
c      write(vtu,'(A)',advance='no') '<VTKFile '
c      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
c      write(vtu,'(A)',advance='no') 'version="1.0" '
c      if (ppiclf_iendian .eq. 0) then
c         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
c      elseif (ppiclf_iendian .eq. 1) then
c         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
c      endif
c
c      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'
c
c      write(vtu,'(A)',advance='yes') '  <FieldData>' 
c      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
c      write(vtu,'(A)',advance='no') 'type="Float32" '
c      write(vtu,'(A)',advance='no') 'Name="TIME" '
c      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
c      write(vtu,'(A)',advance='no') 'format="ascii"> '
c      write(vtu,'(E14.7)',advance='no') ppiclf_time
c      write(vtu,'(A)',advance='yes') ' </DataArray> '
c
c      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
c      write(vtu,'(A)',advance='no') 'type="Int32" '
c      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
c      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
c      write(vtu,'(A)',advance='no') 'format="ascii"> '
c      write(vtu,'(I0)',advance='no') ppiclf_cycle
c      write(vtu,'(A)',advance='yes') ' </DataArray> '
c
c      write(vtu,'(A)',advance='yes') '  </FieldData>'
c      write(vtu,'(A)',advance='no') '  <Piece '
c      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
c      write(vtu,'(I0)',advance='no') npt_total
c      write(vtu,'(A)',advance='yes') '" NumberOfCells="0"> '
c
c! -----------
c! COORDINATES 
c! -----------
c      iint = 0
c      write(vtu,'(A)',advance='yes') '   <Points>'
c      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
c      iint = iint + 3*isize*npt_total + isize
c      write(vtu,'(A)',advance='yes') '   </Points>'
c
c! ----
c! DATA 
c! ----
c      write(vtu,'(A)',advance='yes') '   <PointData>'
c
c
c      do ie=1,PPICLF_LRS
c         write(prostr,'(A1,I2.2)') "y",ie
c         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
c         iint = iint + 1*isize*npt_total + isize
c      enddo
c
c      do ie=1,PPICLF_LRP+1
c         write(prostr,'(A4,I2.2)') "rprop",ie
c         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
c         iint = iint + 1*isize*npt_total + isize
c      enddo
c
c      do ie=1,5
c         write(prostr,'(A3,I2.2)') "tag",ie
c         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
c         iint = iint + 1*isize*npt_total + isize
c      enddo
c
c      write(vtu,'(A)',advance='yes') '   </PointData> '
c
c! ----------
c! END MATTER
c! ----------
c      write(vtu,'(A)',advance='yes') '   <Cells> '
c      write(vtu,'(A)',advance='no')  '    <DataArray '
c      write(vtu,'(A)',advance='no') 'type="Int32" '
c      write(vtu,'(A)',advance='no') 'Name="connectivity" '
c      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
c      write(vtu,'(A)',advance='no') '    <DataArray '
c      write(vtu,'(A)',advance='no') 'type="Int32" '
c      write(vtu,'(A)',advance='no') 'Name="offsets" '
c      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
c      write(vtu,'(A)',advance='no') '    <DataArray '
c      write(vtu,'(A)',advance='no') 'type="Int32" '
c      write(vtu,'(A)',advance='no') 'Name="types" '
c      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
c      write(vtu,'(A)',advance='yes') '   </Cells> '
c      write(vtu,'(A)',advance='yes') '  </Piece> '
c      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '
c
c! -----------
c! APPEND DATA  
c! -----------
c      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
c      close(vtu)
c
c      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
c     >    ,position='append')
c      write(vtu) '_'
c      close(vtu)
c
c      inquire(file=vtufile,size=ivtu_size)
c      endif ! ------------ nid .eq. 0 -------------------------
c
c      call ppiclf_bcast(ivtu_size, isize)
c
c      ! byte-displacements
c      idisp_pos = ivtu_size + isize*(3*stride_len + 1)
c
c      ! how much to write
c      icount_pos = 3*nxx
c      icount_sln = 1*nxx
c      icount_lrp = 1*nxx
c      icount_lip = 1*nxx
c
c      iorank = -1
c
c      ! integer write
c      if (ppiclf_nid .eq. 0) then
c        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
c     >      ,position='append')
c        write(vtu) if_pos
c        close(vtu)
c      endif
c
c      call mpi_barrier(ppiclf_comm,ierr)
c
c      ! write
c        call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
c        call ppiclf_byte_set_view(idisp_pos,pth)
c        call ppiclf_byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
c        call ppiclf_byte_close_mpi(pth,ierr)
c
c      call mpi_barrier(ppiclf_comm,ierr)
c
c      do i=1,PPICLF_LRS
c         idisp_sln = ivtu_size + isize*(3*npt_total 
c     >                         + (i-1)*npt_total
c     >                         + (1)*stride_len
c     >                         + 1 + i)
c
c         ! integer write
c         if (ppiclf_nid .eq. 0) then
c           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
c     >         ,position='append')
c           write(vtu) if_sln
c           close(vtu)
c         endif
c   
c         call mpi_barrier(ppiclf_comm,ierr)
c
c         j = (i-1)*npart_out + 1
c   
c         ! write
cc         if (nxx.ne.0) then
c           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
c           call ppiclf_byte_set_view(idisp_sln,pth)
c           call ppiclf_byte_write_mpi(rout_sln(j),icount_sln,iorank,pth
c     >                             ,ierr)
c           call ppiclf_byte_close_mpi(pth,ierr)
cc        endif ! nxx.ne.0
c      enddo
c
c      do i=1,PPICLF_LRP+1
c         idisp_lrp = ivtu_size + isize*(3*npt_total  
c     >                         + PPICLF_LRS*npt_total
c     >                         + (i-1)*npt_total
c     >                         + (1)*stride_len
c     >                         + 1 + PPICLF_LRS + i)
c
c         ! integer write
c         if (ppiclf_nid .eq. 0) then
c           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
c     >         ,position='append')
c           write(vtu) if_lrp
c           close(vtu)
c         endif
c   
c         call mpi_barrier(ppiclf_comm,ierr)
c
c         j = (i-1)*npart_out + 1
c   
c         ! write
cc        if (nxx.ne.0) then
c           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
c           call ppiclf_byte_set_view(idisp_lrp,pth)
c           call ppiclf_byte_write_mpi(rout_lrp(j),icount_lrp,iorank,pth
c     >                             ,ierr)
c           call ppiclf_byte_close_mpi(pth,ierr)
cc         endif
c      enddo
c
c      do i=1,5
c         idisp_lip = ivtu_size + isize*(3*npt_total
c     >                         + PPICLF_LRS*npt_total
c     >                         + (PPICLF_LRP+1)*npt_total
c     >                         + (i-1)*npt_total
c     >                         + (1)*stride_len
c     >                         + 1 + PPICLF_LRS + PPICLF_LRP + 1 + i)
c         ! integer write
c         if (ppiclf_nid .eq. 0) then
c           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
c     >         ,position='append')
c           write(vtu) if_lip
c           close(vtu)
c         endif
c
c         call mpi_barrier(ppiclf_comm,ierr)
c
c         j = (i-1)*npart_out + 1
c   
c         ! write
cc         if (nxx.ne.0) then
c           call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
c           call ppiclf_byte_set_view(idisp_lip,pth)
c           call ppiclf_byte_write_mpi(rout_lip(j),icount_lip,iorank,pth
c     >                             ,ierr)
c           call ppiclf_byte_close_mpi(pth,ierr)
cc         endif
c      enddo
c
c      if (ppiclf_nid .eq. 0) then
c      vtu=867+ppiclf_nid
c      open(unit=vtu,file=vtufile,status='old',position='append')
c
c      write(vtu,'(A)',advance='yes') '</AppendedData>'
c      write(vtu,'(A)',advance='yes') '</VTKFile>'
c
c      close(vtu)
c      endif
c
c      call ppiclf_printsi(' *End WriteParticleOutletVTU$',ppiclf_cycle)
c
c      return
c      end

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
