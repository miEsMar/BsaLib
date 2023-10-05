!! This file is part of BSA Library.
!! Copyright (C) 2023  Michele Esposito Marzino 
!!
!! BSA Library is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BSA Library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BSA Library.  If not, see <https://www.gnu.org/licenses/>.
submodule(BsaLib) BsaLib_MesherImpl

   use BsaLib_Data
   use BsaLib_MPolicy
   use BsaLib_IO,          only: unit_dump_bfm_
   use BsaLib_MPoint,      only: MPoint_t, MPoint
   use BsaLib_MRectZone,   only: MRectZone_t, MRectZone
   use BsaLib_MTriangZone, only: MTriangZone_t, MTriangZone
   use BsaLib_Functions,   only: prefetchSVDWorkDim_  &
      , NFREQS, NNODES, NNODESL, NLIBS, NLIBSL        &
      , NMODES, NMODES_EFF, MODES                     &
      , NPSDEL, NTCOMPS, NDIRS, TCOMPS, DIRS          &
      , MSHR_SVD_INFO, MSHR_SVD_LWORK, MSHR_SVD_WORK
   implicit none (type, external)

   ! to have a local instance to be referenced
   integer(bsa_int_t) :: NM__, NM_EFF__
   character(len = *), parameter :: bfm_dump_file_name_ = 'dumpfile'
   
   ! BUG: let the user choose how many modes it allows to be covered max.
   integer(int32), parameter :: N_RES_PEAK_IN_BKG_ZONE_DIV_FCT_ = 4

   interface getEquivalentLooperIterator
      module procedure getEquivalentLooperIterator_char
      module procedure getEquivalentLooperIterator_real
   end interface


contains




   module subroutine mainMesher_(m3mf_msh, m3mr_msh)
      use BsaLib_Functions, only: cleanSVDWorkInfo_
      real(bsa_real_t), target, allocatable :: m3mf_msh(:), m3mr_msh(:)
      integer(int32) :: istat
      character(len = 256) :: emsg
      logical :: lflag

#ifdef _BSA_DEBUG
      write(unit_debug_, *) ' @BsaMesherImpl::mainMesher_() : Init BSA-Mesher main...'
#endif

#ifdef _BSA_USE_CACHED_POD_DATA
      print '(1x, a, a/)', NOTEMSG, 'Using version with POD caching.' 
#endif

      NM__     = struct_data%modal_%nm_
      NM_EFF__ = struct_data%modal_%nm_eff_

      msh_bfmpts_pre_  = 0
      msh_bfmpts_post_ = 0
      msh_brmpts_post_ = 0

      allocate(m3mf_msh(dimM_bisp_), stat=istat, errmsg=emsg)
      if (istat /= 0) call allocKOMsg('m3mf_msh', istat, emsg)
      m3mf_msh = 0._bsa_real_t
      m3mf_msh_ptr_ => m3mf_msh

      call PreMesh()
      lflag = settings%i_only_diag_ == 0 .and. settings%i_use_svd_ == 1
#ifndef _BSA_USE_CACHED_POD_DATA
      if (lflag .and. test_no_bfm_mlr_) call cleanSVDWorkInfo_()
#endif

#ifdef _BSA_CHECK_NOD_COH_SVD
      goto 998
#endif


      ! ! post-meshing -> BRM
      allocate(m3mr_msh(dimM_bisp_), stat=istat, errmsg=emsg)
      if (istat /= 0) call allocKOMsg('m3mr_msh', istat, emsg)
      m3mr_msh = 0._bsa_real_t
      m3mr_msh_ptr_ => m3mr_msh
      call Mesh()


      ! cleanup
      998 if (lflag .and. .not. test_no_bfm_mlr_) call cleanSVDWorkInfo_()
      inquire(unit = unit_dump_bfm_, opened = lflag)
      if (lflag) close(unit_dump_bfm_)


      print *
      print '(1x, a, a)', &
         INFOMSG, ' Resume of total n. of points in meshing procedure:'
#ifndef _BSA_USE_CACHED_POD_DATA
      print '(1x, a, a, i0)', MSGCONT, ' PRE-MESH   (BFM) : ', msh_bfmpts_pre_
#endif
      print '(1x, a, a, i0)', MSGCONT, ' POST-MESH  (BFM) : ', msh_bfmpts_post_
      print '(1x, a, a, i0)', MSGCONT, ' POST-MESH  (BRM) : ', msh_brmpts_post_

      
#ifdef _BSA_DEBUG
      write(unit_debug_, *) ' @BsaMesherImpl::mainMesher_() : Init BSA-Mesher main -- ok.'
#endif
   end subroutine mainMesher_




   !> Computes Pre-meshing phase for BFM, 
   !> to be stored (dumped) for bRM meshing actual computation
   !>
   !> NOTE: here, we have the same meshing technique
   !>       no matter the triple mode combination i-j-k.
   !>
   !> TODO: implement also PSDFM computation!
   !>
   !> IMPLEMENTATION DETAILS:
   !> ***********************
   !>
   !> Regarding the bispectrum computation optimisation,
   !> we have since long time already understood it is a complex subject.
   !> Why? A priori 5D matrices (which might be reshaped 3D)
   !> 
   !> What we also know is that, from a specific point of view, 
   !> to obtain the Bispectrum of modal response for a given structure,
   !> the most expensive operation is the obtainment of the
   !> bispectrum of the relative modal FORCE, from the know input,
   !> which we know being the PSDs of the wind turbulent components u,v,w.
   !> Why? Because for each triplet (combination) of modes, at each pair of 
   !> frequencies, the value of the bispectrum of the (triplet of) modal force
   !> is given by the 'projection' of the (FULL) bispectrum of nodal forces
   !> in the modal base (for that specific triplet of modes).
   !> So, if this is the most expensive phase, what could be a solution?
   !> Minimising its impact by minimising the count of times it needs to be performed.
   !>
   !> Could we do it?
   !>
   !> We could think so, since by looking at the "classic" shape of a 
   !> bispectrum of a modal response, still referring to the applied load
   !> (turbulent wind loading), we see that it shows one main peak centered
   !> at the origin (0,0) of the frequency space, plus some little 'crests'
   !> on top of the three main axes (w1=0, w2=0, w1=w2).
   !> We could prove this also analytically.
   !> What about the rest? Almost flat, where one could even suppose constness.
   !> Still, in a "classical" approach values (of BFM) are computed exactly, 
   !> no matter the location in the frequency space w1-w2. It is now clear that,
   !> most of these (EXPENSIVE) operations are useless. Because in those regions
   !> outside the 4 aforementioned we already know that values are almost identical
   !> towards zero (meaning of some orders of magnitude smaller than values at the
   !> central peak for example).
   !> So, the idea here is computing, for what concerns the bispectra computation
   !> at the modal forces stage (for all possible triplets of modes), only what we know
   !> A PRIORI being useful, playing an important role in the estimation of correspondent
   !> value of the bispectra of the relative modal response.
   !> 
   !> What is it? Well, what we said before. The centered main peak, plus
   !> those crests along the three main lines in the Cartesian plane.
   !> 
   !> But how can this be more effective?
   !> 
   subroutine PreMesh()
      real(real64), parameter :: cst_sqrt2d2 = sqrt(2._real64) / 2._real64      
      integer(int32) :: NLims, iost
      
      logical :: iun_open
      real(bsa_real_t) :: base_i, base_j, max_ext
      real(bsa_real_t) :: deltaI_S2_2, deltaI_2_S2_2
      real(bsa_real_t) :: df_I_ref, df_J_ref
      real(bsa_real_t), allocatable :: limits(:)
      type(MPolicy_t), allocatable, target :: policies(:)

      type(MRectZone_t) :: bkgz
      character(len = :), allocatable :: zone_title


#ifdef _BSA_DEBUG
      write(unit_debug_, *) ' @BsaMesherImpl::PreMesh() : Init BSA-Mesher pre meshing phase...'
#endif

      ! check that unit dedicated to dump is not in use already, then, open it
      inquire(unit=unit_dump_bfm_, opened=iun_open)
      do while (iun_open)
         unit_dump_bfm_ = unit_dump_bfm_ + 1
         inquire(unit=unit_dump_bfm_, opened=iun_open)
      enddo
      open(&
         unit=unit_dump_bfm_,      &
         file=bfm_dump_file_name_, &
         form=IO_FORM_UNFORMATTED, &
         access=IO_ACCESS_STREAM,  &
         iostat=iost)
      if (iost /= 0) call bsa_Abort('Cannot open dump file'//bfm_dump_file_name_)
      rewind(unit_dump_bfm_) ! make sure its at init byte

      ! NOTE: iost =  0.
      !       Write it twice to make space for info that
      !       will be OVERRIDEN once we have it
      !       namely "dimM_bisp_" and "MAX ZONE n. of points" 
      write(unit_dump_bfm_) iost
      write(unit_dump_bfm_) iost
      write(unit_dump_bfm_) iost
      write(unit_dump_bfm_) iost



      write(*, '(1x, a)') '-----------------------------------------------------------'
      write(*, '(1x, a)') '--------------------    PRE - MESH     --------------------'
      write(*, '(1x, a)') '-----------------------------------------------------------'


      ! NOTE: bkg_peak_width_ is already given in [Hz]
      ! BUG:  is really a 3x3 matrix ??
      if (do_restrict_bkgpeak_) then
         bkg_peakw_ = maxval(struct_data%bkg_peak_width_(:, 1))
      else
         bkg_peakw_ = maxval(struct_data%bkg_peak_width_(:, :))
      endif
      base_i = bkg_peakw_ * settings%bkg_area_extension_
      base_j = base_i

      deltaI_S2_2   = base_i * cst_sqrt2d2
      deltaI_2_S2_2 = deltaI_S2_2 / 2._bsa_real_t

      max_ext = getMaxSpaceExtension_() ! Get max point in space to reach


      ! NOTE: use static module variable msh_ZoneLimsInterestModes 
      !       so that we can reuse in post-meshing phase. 
      call prefetchZoneLimits_(base_i / 2, limits, policies, NLims, msh_ZoneLimsInterestModes)

      if (settings%i_only_diag_ == 0 .and. settings%i_use_svd_ == 1) call prefetchSVDWorkDim_()


#ifdef _BSA_CHECK_NOD_COH_SVD
      block
         double precision :: tmpv(1, NNODESL)

         call wd%getFull2DNodCorrMat(NNODES, nod_corr_full_)
         if (allocated(nod_corr_full_)) then
            allocate(nod_corr_EVLs_(NNODESL))
            allocate(nod_corr_EVTs_(NNODESL, NNODESL))

            nod_corr_full_ = nod_corr_full_(struct_data%n_load_, struct_data%n_load_)
            nod_corr_EVTs_ = nod_corr_full_

            call dgesvd(&
               'O' &           ! min(M,N) columns of U are overwritten on array A (saves memory)
               , 'N' &           ! no rows of V are computed
               , NNODESL    &    ! n. of rows M
               , NNODESL    &    ! n. of cols N
               , nod_corr_EVTs_ &! A matrix (overwritten with left-singular vectors)
               , NNODESL    &
               , nod_corr_EVLs_ &! singular values
               , tmpv       &    ! U
               , 1          & 
               , tmpv       &    ! VT
               , 1          &
               , MSHR_SVD_WORK  &
               , MSHR_SVD_LWORK &
               , MSHR_SVD_INFO  &
            )
            if (MSHR_SVD_INFO /= 0) call bsa_Abort("Error while computing SVD of nodal correlation")
            
            write(4397, *) NNODESL
            do iost = 1, NNODESL
               write(4397, *) nod_corr_full_(:, iost)
            enddo
            
            write(4398, *) NNODESL
            write(4398, *) nod_corr_EVLs_
            do iost = 1, NNODESL
               write(4398, *) nod_corr_EVTs_(:, iost)
            enddo
         endif
      end block
#endif

      
      !===================================================================================
      ! BKG peak
      !
      ! NOTE: we keep it in memory, since it will serve as reference
      !       point for other nearby zones correct identification.
      call timer%init()
      zone_title = 'BKG center peak'
      bkgz       = MRectZone(0._bsa_real_t, zone_title)
      if (settings%i_bisp_sym_ == BSA_SPATIAL_SYM_HALF) then
         base_i = base_i / 2._bsa_real_t
         call bkgz%define(MPoint(0._bsa_real_t, - base_i), 'i', base_i, base_j)
      else
         call bkgz%define(MPoint_t(), 'c', base_i, base_j)
      endif
      call bkgz%setPolicy(MPolicy_PEAK)


      ! NOTE: 0 denotes that interest modes are to be inferenced from index 1.
      !       In fact, there are 3 scenarios.
      !          1.  next zone is pre-peak, and next peak interest modes' start from 1.
      !              BKG does not include any resonant peak.
      !          2.  next zone is pre-peak, and next peak interest modes' DO NOT start from 1.
      !              BKG does include some resonant peaks (from 1-less-index or next peak zone)
      !          3.  next zone is peak.
      !              BKG does include this, plus all previous resonant peaks.
      call bkgz%setInterestModeIndexPtr(0)
      

      iost = settings%bkg_base_rfmnt_
      ! if (.not. test_no_bfm_mlr_) then
      !    iost = iost - 1  ! get n.  of segments
      !    iost = iost / min(bkgz%policy_%interp_bfm_I_fct_, bkgz%policy_%interp_bfm_J_fct_)
      !    iost = iost + 1  ! get back actualised n. of points
      ! endif
      if (settings%i_bisp_sym_ == BSA_SPATIAL_SYM_HALF) then
         call bkgz%setRefinements(((iost-1)/2) + 1, iost)
      else
         call bkgz%setRefinements(iost, iost)
      endif

      ! NOTE: if HALF symmetry, we need to check for the correct n. of ref. points.
      if (settings%i_bisp_sym_ == BSA_SPATIAL_SYM_HALF) then
         iost = (bkgz%nj_ - 1) / 2 + 1
         if (bkgz%ni_ /= iost) call bkgz%setRefinements(iost, bkgz%nj_, .true.)
      endif

      ! backup reference deltas
      df_I_ref = bkgz%deltaf_I_
      df_J_ref = bkgz%deltaf_J_
      if (.not. abs(df_I_ref - df_J_ref) < MACHINE_PRECISION) &
         call bsa_Abort("Freq deltas of BKG peak zone differ. Should be the same.")

#ifdef _BSA_EXPORT_POD_TRUNC_INFO

# ifdef _OPENMP
#  define __export_POD_trunc_id__  omp_get_thread_num()+1
# else
#  define __export_POD_trunc_id__  1
# endif

      allocate(do_export_POD_trunc_(16))
      do_export_POD_trunc_    = .false.
      do_export_POD_trunc_(1) = .true.  ! <-- NO OMP parall here regardless.
#endif
      
      call bkgz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
      call logger_debug%logZonePremeshingTotTime(&
         zone_title, timer%time(), msh_bfmpts_pre_, .true.)
#endif
         
      if (.not. allocated(limits)) goto 998  ! NOTE: BKG zone covers them all, bad..


#define __return_debug__ goto 998
#ifdef _BSA_CHECK_NOD_COH_SVD
      __return_debug__
#endif


      ! __return_debug__


      ! ALL OTHER ZONES (IF NOT BKG COVERS EVERYTHING)

      write(*, '(1x, a, a, /, 10(" ", f10.4))') &
         INFOMSG, '  Limits frontiers:', limits

      block

         !> Number of main directions ['NORTH', 'EAST ', 'SOUTH', 'WEST ']
         integer(int32), parameter :: N_DIRS_FULL = 4_int32
         integer(int32), parameter :: N_DIRS_HALF = 3_int32
         integer(int32) :: n_dirs_

         !> Main directions labels 
         !> ['NORTH', 'EAST ', 'SOUTH', 'WEST ']
         character(len = 5), parameter  :: DIRS_LABELS(4) = ['NORTH', 'EAST ', 'SOUTH', 'WEST ']

         !> Main diag-crests directions labels
         !> ['NORTH-EAST', 'SOUTH-EAST', 'SOUTH-WEST', 'NORTH-WEST']
         character(len = 10), parameter :: DIRS_DIAG_LABELS(4) = &
            ['NORTH-EAST', 'SOUTH-EAST', 'SOUTH-WEST', 'NORTH-WEST']
         
         !> Set of pair rotations, per main direction
         real(real64), parameter :: ROTATIONS(6) = &
            [CST_PIt3d2, 0._real64, CST_PId2, CST_PIGREC, CST_PIt3d2, 0._real64]

         !> Which base is actually passed, depending on N-E-S-W direction
         character(len = 1), parameter :: COORDS_DIR_CH(4) = ['j', 'i', 'j', 'i']

         !> Signs to apply at limits, depending on N-E-S-W direction
         integer(int32), parameter :: LIM_SIGN_DIRS(4) = [1, 1, -1, -1]

         integer(int32), parameter :: LEFT_RZ_SIGNS(4)  = [1, -1, -1, 1]

         real(real64), parameter :: DF_FCT_CST = 0._real64

         integer(int32)   :: i, NLimsP1, idir, iim, ilim, nim, idir_t2
         integer(int32)   :: id_im_last, ilim_init_ = 0, n_bfm_pts_pre_
         real(bsa_real_t) :: lim, rtmp
         real(bsa_real_t) :: df_I, df_J
         integer(int32)   :: sign_dir
         real(bsa_real_t) :: maxF, base_i_
         real(bsa_real_t), allocatable :: maxext_sym_(:), bases_i_(:)

         logical :: warn_zone_over_limits = .false.

         type(MPolicy_t), pointer :: policy_ptr => null()
         type(MPolicy_t)          :: pol

         character(len = 1) :: coord_dir

         ! allocatables depending on Nlim +1
         real(bsa_real_t), allocatable   :: rots(:)
         real(bsa_real_t), allocatable   :: deltas(:, :)
         integer(bsa_int_t), allocatable :: refmts(:, :), inter_modes_(:)
         integer(bsa_int_t), allocatable :: int_modes_(:)
         real(bsa_real_t) :: rval

         !> general rectangular zone
         type(MRectZone_t) :: rz
         type(MPoint_t)    :: basePts(4), ptI, ptE

         integer(int32) :: N_THREADS_MIN_


         ! 
         NLimsP1 = NLims + 1
         allocate(rots(NLimsP1))
         rots   = 0._bsa_real_t
         allocate(deltas(2, NLimsP1))
         deltas = 0._bsa_real_t
         allocate(refmts(2, NLimsP1))
         refmts = 0_bsa_int_t
         allocate(inter_modes_(NLimsP1))
         inter_modes_ = 0

         
         maxF       = limits(NLims) + deltaI_2_S2_2
         id_im_last = size(msh_ZoneLimsInterestModes) - 1
         
         pol = MPolicy_PRE_PEAK_2

         ! extend limits and policies to covering after last peak
         limits(NLimsP1)   = maxF
         policies(NLimsP1) = pol

#ifdef _BSA_DEBUG
         write(*, *) '  Interest modes : '
         write(*, *) msh_ZoneLimsInterestModes
#endif

         ! array of base points from BKG peak reference
         basePts = [&
            bkgz%getAPoint(), &
            bkgz%Ept_,        &
            bkgz%getBPoint(), & 
            bkgz%Ipt_         &
         ]

         allocate(bases_i_(4))
         if (settings%i_bisp_sym_ == BSA_SPATIAL_SYM_HALF) then
            n_dirs_ = N_DIRS_HALF
            bases_i_(1) = base_i
            bases_i_(2) = base_j
            bases_i_(3) = base_i
         
            N_THREADS_MIN_ = 1
         else
            n_dirs_     = N_DIRS_FULL
            bases_i_(:) = base_i

            N_THREADS_MIN_ = 2
         endif


         !$ if (allocated(zone_title)) deallocate(zone_title)


         !$omp parallel do &
         !$omp   default(firstprivate),        &
         !$omp   private(idir, zone_title),    &
         !$omp   shared(DIRS_LABELS, ROTATIONS, NLimsP1, COORDS_DIR_CH, LIM_SIGN_DIRS &
         !$omp          , policies, limits, df_I_ref, df_J_ref, msh_ZoneLimsInterestModes &
         !$omp          , refmts, deltas, inter_modes_, basePts, base_i, bases_i_ &
         !$omp          , struct_data, wd, settings, logger_debug &
         !$omp          , id_im_last, maxF, NLims, getBFM_msh, pol &
         !$omp          , NFREQS, NNODES, NNODESL, NLIBS, NLIBSL &
         !$omp          , NMODES, NMODES_EFF, MODES &
         !$omp          , NPSDEL, NTCOMPS, NDIRS, TCOMPS, DIRS          &
         !$omp          , MSHR_SVD_INFO, MSHR_SVD_LWORK, MSHR_SVD_WORK  &
         !$omp          , do_export_POD_trunc_   &
         !$omp          , msh_NZones, msh_bfmpts_pre_, msh_max_zone_NPts, m3mf_msh_ptr_), &
         !$omp   num_threads(n_dirs_)
         do idir = 1, n_dirs_

            n_bfm_pts_pre_ = 0

            call timer%init()
            zone_title = 'Crest zone at  '//DIRS_LABELS(idir)

            ! init rect zone
            call rz%zoneName(zone_title)
            call rz%setRotation(ROTATIONS(idir + 1))

            ! get base init point
            ptI = basePts(idir)

            coord_dir = COORDS_DIR_CH(idir)
            sign_dir  = LIM_SIGN_DIRS(idir)

            base_i_ = bases_i_(idir)

            iim = 1
            do ilim = 1, NLims
               
               lim = limits(ilim)
               policy_ptr => policies(ilim)

               ! BUG: fix this frequencies definition!
               df_I = df_I_ref * policy_ptr%delta_fI_fct_
               df_J = df_J_ref * policy_ptr%delta_fJ_fct_

               call rz%setPolicy(policy_ptr)
               
               ! get this zone interest modes
               if (idir == 2) inter_modes_(ilim) = iim
               if (ilim == 1 .or. policies(ilim+1) == MPolicy_PEAK) then
                  nim = msh_ZoneLimsInterestModes(iim)
                  if (nim < 0) then ! pre-peak
                     iim = iim + 1
                     nim = msh_ZoneLimsInterestModes(iim)
                  endif
                  int_modes_ = msh_ZoneLimsInterestModes(iim + 1  :  iim + nim)

#ifdef _BSA_DEBUG
                  !$omp critical
                  if (idir == 1) print *, int_modes_
                  !$omp end critical
#endif
               
                  ! set current zone interest modes pointer (before update)
                  call rz%setInterestModeIndexPtr(iim)
                  iim  = iim + nim + 1
               endif

               call rz%defineFromEndPtCoordAndBase(&
                  ptI, sign_dir * lim, coord_dir, base_i_, 'i', df_I, df_J)


               ! BUG: optimise this!!
               select case (idir)
                  case (1, 3)
                     call ptI%move(0._bsa_real_t, sign_dir * lim - ptI%freqJ())  ! move along Y
                  case (2, 4)
                     call ptI%move(sign_dir * lim - ptI%freqI(), 0._bsa_real_t)  ! move along X
               end select


! #ifdef _BSA_DEBUG
!                if (idir == 1) then
!                   write(*, '( 1x, a, i5, l, g16.5, "  ->  ", *(" ", i0) )', advance='yes') &
!                      ' ilim,  isPeak,  rval,  int_modes_,  policies = ', &
!                      ilim, policy_ptr == MPolicy_PEAK, rval, int_modes_ &
!                      , policy_ptr%delta_fI_fct_, policy_ptr%delta_fJ_fct_
!                endif
! #endif

! #ifdef _OPENMP
               ! NOTE: store refinements for later use
               if (idir == 2) then
                  refmts(1:2, ilim) = rz%refinements()
                  deltas(1, ilim)   = rz%deltaf_I_
                  deltas(2, ilim)   = rz%deltaf_J_
               endif
! #endif

#ifdef _BSA_EXPORT_POD_TRUNC_INFO
               if (idir == 1 .or. idir == 3) do_export_POD_trunc_(__export_POD_trunc_id__) = .true.
#endif
               call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
               n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif

            enddo ! limits

            df_I = df_I_ref * pol%delta_fI_fct_
            df_J = df_J_ref * pol%delta_fJ_fct_

            call rz%setPolicy(pol)
            call rz%setInterestModeIndexPtr(id_im_last)  ! this zone has only Last mode of intrest
            call rz%defineFromEndPtCoordAndBase(&
               ptI, sign_dir * maxF, coord_dir, base_i_, 'i', df_I, df_J)

! #ifdef _OPENMP
            if (idir == 2) then
               refmts(1:2, NLimsP1)  = rz%refinements()
               deltas(1, NLimsP1)    = rz%deltaf_I_
               deltas(2, NLimsP1)    = rz%deltaf_J_
            endif
! #endif

#ifdef _BSA_EXPORT_POD_TRUNC_INFO
            if (idir == 1 .or. idir == 3) do_export_POD_trunc_(__export_POD_trunc_id__) = .true.
#endif
            call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
            n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif

#ifndef _BSA_USE_CACHED_POD_DATA
            !$omp critical
            call logger_debug%logZonePremeshingTotTime(&
               zone_title, timer%time(), n_bfm_pts_pre_, .true.)
            !$omp end critical
#endif
            
         enddo ! n dirs
         !$omp end parallel do
      
         deallocate(bases_i_)

         ! BUG: might be removed, code duplication for little CPU improvement..
#ifndef _OPENMP
         refmts(1:2, NLimsP1)  = rz%refinements()
         deltas(1, NLimsP1)    = rz%deltaf_I_
         deltas(2, NLimsP1)    = rz%deltaf_J_
#endif
         inter_modes_(NLimsP1) = id_im_last


         print '(1x, a, a, i0, a/)', &
            INFOMSG, 'Done with   ', msh_NZones, '  pre meshing zones.'


#ifdef _BSA_EXPORT_POD_TRUNC_INFO
         ! From now on, no need for this anymore.
         do_export_POD_trunc_(:) = .false.
#endif



         ! __return_debug__




         if (ipre_mesh_type == BSA_PREMESH_TYPE_DIAG_CREST_NO) then


            block
               character(len = 1), allocatable :: bases_ch(:)
               real(bsa_real_t) :: init_freq_, rbase_
               real(bsa_real_t) :: main_rz_rot_, left_rz_rot_, right_rz_rot_
               real(bsa_real_t) :: delta_main_rz_, rlimit_
               character(len = 40) :: z_name_
               character(len = 1)  :: left_known_coord_, right_known_coord_
               integer(int32)      :: left_sign_, right_sign_, main_refs_
               integer(int32)      :: idirP1_

               if (allocated(zone_title)) deallocate(zone_title)
               if (allocated(rots))       deallocate(rots)
               allocate(bases_ch(N_DIRS_FULL))
               bases_ch = ' '

               ! reset correct base points (per direction)
               ptI = basePts(1)   ! backup 4th quadrant pt
               basePts(1:3) = basePts(2:4)
               basePts(4)   = ptI

               ! computing base of first opening rect zone
               init_freq_ = basePts(1)%freqI()
               
               main_refs_ = maxval(refmts(:, 1))
               bases_ch   = getEquivalentLooperIterator(N_DIRS_FULL, 'ij')

               if (settings%i_bisp_sym_ == BSA_SPATIAL_SYM_HALF) n_dirs_ = n_dirs_ - 1  ! ==2

               !$omp parallel do &
               !$omp   default(firstprivate), &
               !$omp   shared(ROTATIONS, LIM_SIGN_DIRS, main_refs_, bases_ch, inter_modes_ &
               !$omp          , basePts, policies, deltas, LEFT_RZ_SIGNS, DIRS_DIAG_LABELS &
               !$omp          , struct_data, wd, settings, logger_debug, limits &
               !$omp          , NFREQS, NNODES, NNODESL, NLIBS, NLIBSL &
               !$omp          , NMODES, NMODES_EFF, MODES &
               !$omp          , NPSDEL, NTCOMPS, NDIRS, TCOMPS, DIRS &
               !$omp          , MSHR_SVD_INFO, MSHR_SVD_LWORK, MSHR_SVD_WORK  &
               !$omp          , msh_NZones, msh_bfmpts_pre_, msh_max_zone_NPts, m3mf_msh_ptr_), &
               !$omp   num_threads(n_dirs_)
               do idir = 1, n_dirs_

                  call timer%init()
                  idirP1_ = idir + 1

                  ! treat initial opening rect zone
                  ! NOTE: treat it singularly, since its complete definition 
                  !       will serve as base for definition of later rect zones.
                  write(unit=z_name_, fmt='(a, i0, 3a)') &
                     'Zones in quadrant n.  ', idir, '  (', DIRS_DIAG_LABELS(idir), ')'
                  call rz%zoneName(z_name_)

                  left_rz_rot_  = ROTATIONS(idir)
                  main_rz_rot_  = ROTATIONS(idirP1_)
                  right_rz_rot_ = ROTATIONS(idirP1_ + 1)   ! idir + 2
                  call rz%setRotation(main_rz_rot_)
                  
                  pol = policies(1)
                  call rz%setPolicy(pol)

                  ptI    = basePts(idir)
                  rbase_ = limits(1) - init_freq_
                  if (rbase_ <= 0._bsa_real_t) call bsa_Abort("Negative base.")

                  ! BUG: better understand how to correctly define these zones..
                  call rz%define(ptI, 'i', rbase_, rbase_)
                  call rz%setRefinements(main_refs_, main_refs_)
                  iim = 1
                  call rz%setInterestModeIndexPtr(iim)
! #ifdef _BSA_DEBUG
!                   print '(1x, a, 2i5, 4g12.5)', '  @idir, ilim,  ptI,  ptE  =  ', &
!                      idir, 1, rz%Ipt_%freqI(), rz%Ipt_%freqJ(), rz%Ept_%freqI(), rz%Ept_%freqJ()
! #endif
                  call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                  n_bfm_pts_pre_ = rz%zoneTotNPts()
#endif

                  ptE = rz%Ept_   ! saving end point
                  left_sign_  = LEFT_RZ_SIGNS(idir)
                  right_sign_ = LIM_SIGN_DIRS(idir)

                  left_known_coord_  = bases_ch(idir)
                  right_known_coord_ = bases_ch(5 - idir)
                  
                  
                  do ilim = 2, NLimsP1

                     call rz%setInterestModeIndexPtr(inter_modes_(ilim))

                     rlimit_ = limits(ilim)
                     rbase_  = rlimit_ - limits(ilim - 1)
                     ptI     = ptE

                     !
                     ! main rect zone
                     !

                     ! NOTE: this policy is shared by the other zones as well
                     call rz%setPolicy(policies(ilim))
                     
                     delta_main_rz_ = minval(deltas(:, ilim))
                     call rz%setRotation(main_rz_rot_)
                     call rz%defineFromEndPtCoordAndBase(&
                        ptI, left_sign_ * rlimit_, 'j', &
                        rbase_, left_known_coord_, delta_main_rz_, delta_main_rz_)
! #ifdef _BSA_DEBUG
!                      print '(1x, a, 2i5, 4g12.5)', '  @idir, ilim,  ptI,  ptE  =  ', &
!                         idir, ilim, rz%Ipt_%freqI(), rz%Ipt_%freqJ(), rz%Ept_%freqI(), rz%Ept_%freqJ()
! #endif
                     call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                     n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif

                     ptE = rz%Ept_   ! saving end point (becomes init point at next iter)

                     delta_main_rz_ = deltas(2, ilim)

                     !
                     ! side (left) rect zone (blue)
                     !
                     call rz%setRotation(left_rz_rot_)
                     call rz%defineFromEndPtCoordAndBase(&
                        ptI, left_sign_ * init_freq_, left_known_coord_, &
                        rbase_, 'i', delta_main_rz_, delta_main_rz_)
! #ifdef _BSA_DEBUG
!                      print '(1x, a, 2i5, 4g12.5)', '  @idir, ilim,  ptI,  ptE  =  ', &
!                         idir, ilim, rz%Ipt_%freqI(), rz%Ipt_%freqJ(), rz%Ept_%freqI(), rz%Ept_%freqJ()
! #endif
                     call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                     n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif

                     !
                     ! side (right) rect zone (pink)
                     !
                     call rz%setRotation(right_rz_rot_)
                     call rz%defineFromEndPtCoordAndBase(&
                        ptI, right_sign_ * init_freq_, right_known_coord_, &
                        rbase_, 'j', delta_main_rz_, delta_main_rz_)
! #ifdef _BSA_DEBUG
!                      print '(1x, a, 2i5, 4g12.5)', '  @idir, ilim,  ptI,  ptE  =  ', &
!                         idir, ilim, rz%Ipt_%freqI(), rz%Ipt_%freqJ(), rz%Ept_%freqI(), rz%Ept_%freqJ()
! #endif
                     call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                     n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif

                  enddo ! ilim = 2, NLimsP1

#ifndef _BSA_USE_CACHED_POD_DATA
						!$omp critical
                  call logger_debug%logZonePremeshingTotTime(&
                     z_name_, timer%time(), n_bfm_pts_pre_, .true.)
                  !$omp end critical
#endif
               enddo ! idir
               !$omp end parallel do

               print '(1x, a, a, i0, a/)', &
                  INFOMSG, 'Done with   ', msh_NZones, '  pre meshing zones.'
            end block





         elseif (ipre_mesh_type == BSA_PREMESH_TYPE_DIAG_CREST_YES) then



            !===================================================
            ! RECT PADDING ZONE 1 (NE & SW)
            !
            pol  = MPolicy_PAD_ZONE_INTERN
            df_I = df_I_ref * pol%delta_fI_fct_
            df_J = df_J_ref * pol%delta_fJ_fct_

            ! BUG: "0" == all modes, not optimal at all
            call rz%setInterestModeIndexPtr(0)

#ifdef _OPENMP
            if (allocated(zone_title)) deallocate(zone_title)
#endif

            !$omp parallel do &
            !$omp   default(firstprivate), &
            !$omp   private(zone_title),   &
            !$omp   shared(ROTATIONS, LIM_SIGN_DIRS &
            !$omp          , maxF, basePts, df_I, df_J, pol &
            !$omp          , struct_data, wd, settings, logger_debug &
            !$omp          , NFREQS, NNODES, NNODESL, NLIBS, NLIBSL &
            !$omp          , NMODES, NMODES_EFF, MODES &
            !$omp          , NPSDEL, NTCOMPS, NDIRS, TCOMPS, DIRS &
            !$omp          , MSHR_SVD_INFO, MSHR_SVD_LWORK, MSHR_SVD_WORK  &
            !$omp          , msh_NZones, msh_bfmpts_pre_, msh_max_zone_NPts, m3mf_msh_ptr_), &
            !$omp   num_threads(N_THREADS_MIN_)
            do idir = 1, N_THREADS_MIN_

               call timer%init()

               idir_t2 = idir * 2

               if (idir == 1) then
                  zone_title = 'Internal rect Padding (NORTH-EAST)'
               else
                  zone_title = 'Internal rect Padding (SOUTH-WEST)'
               endif
               call rz%zoneName(zone_title)
               call rz%setRotation(ROTATIONS(idir_t2))
               call rz%setPolicy(pol)

               sign_dir = LIM_SIGN_DIRS(idir_t2)
               rtmp     = sign_dir * maxF
               call rz%defineFromDeltas(&
                  basePts(idir_t2), 'i', df_I, df_J, rtmp, rtmp, force=.true.)
               call rz%compute()

#ifndef _BSA_USE_CACHED_POD_DATA
               !$omp critical
               call logger_debug%logZonePremeshingTotTime(&
                  zone_title, timer%time(), rz%zoneTotNPts(), .true.)
               !$omp end critical
#endif
               
            enddo ! idir
            !$omp end parallel do

            print '(1x, a, a, i0, a/)', &
               INFOMSG, 'Done with   ', msh_NZones, '  pre meshing zones.'


            !===================================================
            ! SE & NW diagonal crests
            !

            ! reset base points (use old variable)
            basePts(2) = basePts(1)
            basePts(1) = basePts(3)

            block
               type(MTriangZone_t) :: tz
               type(MPoint_t)      :: ptA, ptB
               real(bsa_real_t)    :: lim_I, lim_J, tmprot, tmprots(2), tmpdelta

               integer(int32) :: rftmp(2), ni, nj, idrot

               ! reset rotations for inclined crests
               tmprots(1) = CST_PId2   + CST_PId4
               tmprots(2) = CST_PIt3d2 + CST_PId4


#ifdef _OPENMP
               if (allocated(zone_title)) deallocate(zone_title)
               if (allocated(rots))       deallocate(rots)
               allocate(rots(NLimsP1 * 2))
#endif


               !$omp parallel do &
               !$omp   default(firstprivate), &
               !$omp   private(zone_title),   &
               !$omp   shared(ROTATIONS, LIM_SIGN_DIRS, DIRS_DIAG_LABELS &
               !$omp          , msh_ZoneLimsInterestModes, bkgz, deltaI_S2_2 &
               !$omp          , maxF, basePts, tmprots, ipre_mesh_mode &
               !$omp          , NLimsP1, limits, refmts, policies, deltas &
               !$omp          , base_i, id_im_last, NLims &
               !$omp          , struct_data, wd, settings, logger_debug &
               !$omp          , NFREQS, NNODES, NNODESL, NLIBS, NLIBSL &
               !$omp          , NMODES, NMODES_EFF, MODES &
               !$omp          , NPSDEL, NTCOMPS, NDIRS, TCOMPS, DIRS &
               !$omp          , MSHR_SVD_INFO, MSHR_SVD_LWORK, MSHR_SVD_WORK  &
               !$omp          , msh_NZones, msh_bfmpts_pre_, msh_max_zone_NPts, m3mf_msh_ptr_), &
               !$omp   num_threads(N_THREADS_MIN_)
               do idir = 1, N_THREADS_MIN_

                  n_bfm_pts_pre_ = 0
                  idir_t2        = idir * 2
                  
                  call timer%init()
                  zone_title = 'Diagonal crest  -  '//DIRS_DIAG_LABELS(idir_t2)

                  iim  = 1  ! reset pointer for interest modes
                  
                  ! caching
                  sign_dir = LIM_SIGN_DIRS(idir_t2) ! BUG: remove, just for check
                  tmpdelta = deltaI_S2_2 * sign_dir

                  ! opening triang zone
                  pol = MPolicy_PEAK
                  tz  = MTriangZone(zone_title)

                  lim_I = basePts(idir)%freqI()
                  lim_J = basePts(idir)%freqJ()
                  ptI   = MPoint(lim_I + tmpdelta,  lim_J)
                  ptA   = MPoint(lim_I,  lim_J - tmpdelta)

                  ! TODO: check these lines
                  rftmp = bkgz%refinements()

                  ! NOTE: take same interest modes as first zone limit
                  call tz%setInterestModeIndexPtr(1)

                  call tz%setRefinements(rftmp(1), rftmp(2))
                  call tz%setPolicy(pol)
                  call tz%define(basePts(idir), ptI, ptA)
                  call tz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                  n_bfm_pts_pre_ = n_bfm_pts_pre_ + tz%zoneTotNPts()
#endif

                  tmprot = tmprots(idir)  ! main rotation for current direction


                  if (ipre_mesh_mode == BSA_PREMESH_MODE_BASE) then

                     
                     ! IMPLEMENT & VERIFY
                     print '(1x, a, a)', &
                        ERRMSG, '"BASE"  pre mesh mode not yet implemented.'
                     call bsa_Abort()


                  elseif (ipre_mesh_mode == BSA_PREMESH_MODE_ZONE_REFINED) then


                     ! Find all covered limits by initial triang-zone
                     if (ilim_init_ == 0) then
                        ilim_init_ = 1
                        do while(abs(ptI%freqI()) > limits(ilim_init_) .and. ilim_init_ <= NLims)
                           ilim_init_ = ilim_init_ + 1
                        enddo
                        !$omp critical
                        print '(/ 1x, a, a, i0, a)', &
                           WARNMSG, 'Init triang zone at diag-crest covers   ', &
                              ilim_init_ - 1, '   limit(s).'
                        !$omp end critical
                        warn_zone_over_limits = ilim_init_ == NLims + 1
                     endif


                     call rz%zoneName(zone_title)


                     ! Renew rotations Looper Iterator (since now inclined zones)
                     ! NOTE: *2 because we index it twice in the same loop!
                     ! BUG: adapt to actual n. of limits done.
                     ! BUG: aviod creation of temporary
                     rots = getEquivalentLooperIterator(NLimsP1 * 2, &
                        [tmprot - tmprots(1), tmprot + tmprots(1)])
                     idrot = 1


                     ! Loop over n of (remained) limits 
                     ! NOTE: +1 for that little padding after last peak
                     do ilim = ilim_init_, NLimsP1


                        lim = limits(ilim)

                        ! =================================
                        ! Main rect zone

                        ! NOTE: same refinements as previous crests
                        call rz%setRefinements(bkgz%ni_, refmts(2, ilim))

                        ! set pointer to interest modes
                        call rz%setInterestModeIndexPtr(iim)
                        call rz%setRotation(tmprot)
                        call rz%setPolicy(policies(ilim))
                        call rz%defineFromEndPtCoordAndBase(&
                           ptI, - lim * LIM_SIGN_DIRS(idir_t2), 'j', base_i, 'i', called=.false.)
                        call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                        n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif

                        ! NOTE: save since RZ will be overridden
                        ptA = rz%getAPoint()
                        ptB = rz%getBPoint()
                        ptE = rz%Ept_
                        if (policies(ilim) == MPolicy_PEAK) then
                           pol = MPolicy_CREST
                        else
                           pol = MPolicy_BASIN
                        endif



                        ! =================================
                        ! north-east (south-west)

                        ! 1. triang leveling zone

                        call tz%setInterestModeIndexPtr(iim)

                        ! update interest modes pointer for next iteration
                        nim = msh_ZoneLimsInterestModes(iim)
                        if (nim < 0) then ! pre-peak
                           iim = iim + 1
                        else! peak
                           iim = iim + nim + 1
                        endif

                        ptI = MPoint(lim * LIM_SIGN_DIRS(idir_t2), rz%Ipt_%freqJ())
                        call tz%setRefinements(rz%nj_, rz%nj_)
                        call tz%setPolicy(pol)
                        call tz%define(ptI, ptA, rz%Ipt_)
                        call tz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                        n_bfm_pts_pre_ = n_bfm_pts_pre_ + tz%zoneTotNPts()
#endif
                        rtmp = tz%baseI()

                        ! 2. rect closing zone
                        call rz%setRotation(rots(idrot))
                        idrot = idrot + 1
                        call rz%setPolicy(pol)
                        call rz%defineFromEndPtCoordAndBase(&
                           rz%Ipt_, lim_J, 'j', rtmp, 'i', &
                           deltas(2, ilim) * pol%delta_fI_fct_, deltas(2, ilim) * pol%delta_fJ_fct_)
                        call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                        n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif



                        ! =================================
                        ! south-west (north-east)

                        ! 1. triang leveling zone
                        ptI = MPoint(ptB%freqI(), - lim * LIM_SIGN_DIRS(idir_t2))
                        call tz%define(ptI, ptB, ptE)
                        call tz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                        n_bfm_pts_pre_ = n_bfm_pts_pre_ + tz%zoneTotNPts()
#endif

                        ! 2. rect closing zone
                        !
                        ! NOTE: the +1 should reproduce the next() method
                        !       called recursively inside the ilim loop..
                        call rz%setRotation(rots(idrot))
                        idrot = idrot + 1
                        call rz%defineFromEndPtCoordAndBase(&
                           tz%Cpt_, lim_I, 'i', rtmp, 'i', &
                           deltas(2, ilim) * pol%delta_fI_fct_, deltas(2, ilim) * pol%delta_fJ_fct_)
                        
                        call rz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                        n_bfm_pts_pre_ = n_bfm_pts_pre_ + rz%zoneTotNPts()
#endif

                        
                        ! save for starting next iteration
                        ptI = ptA

                     enddo ! N limits + 1


                     ! closing triangle
                     if (.not. warn_zone_over_limits) then

                        ! BUG: improve
                        pol = MPolicy_PEAK

                        call tz%setInterestModeIndexPtr(id_im_last)

                        call tz%setPolicy(pol)

                        ni = settings%bkg_base_rfmnt_ * pol%delta_fI_fct_
                        nj = settings%bkg_base_rfmnt_ * pol%delta_fJ_fct_
                        call tz%setRefinements(ni, nj)

                        ptI = MPoint(maxF * LIM_SIGN_DIRS(idir_t2), - maxF * LIM_SIGN_DIRS(idir_t2))

                        ptA = MPoint(ptI%freqI() - (deltaI_S2_2 * LIM_SIGN_DIRS(idir_t2)), ptI%freqJ())
                        ptB = MPoint(ptI%freqI(), ptI%freqJ() + (deltaI_S2_2 * LIM_SIGN_DIRS(idir_t2)))

                        call tz%define(ptI, ptA, ptB)
                        call tz%compute()
#ifndef _BSA_USE_CACHED_POD_DATA
                        n_bfm_pts_pre_ = n_bfm_pts_pre_ + tz%zoneTotNPts()
#endif

                     endif ! closing triangle, if not warn_zone_over_limits

#ifndef _BSA_USE_CACHED_POD_DATA
                     !$omp critical
                     call logger_debug%logZonePremeshingTotTime(&
                        zone_title, timer%time(), n_bfm_pts_pre_, .true.)
                     !$omp end critical
#endif
                  endif ! pre mesh mode

               enddo ! idir
               !$omp end parallel do

               print '(1x, a, a, i0, a/)', &
                  INFOMSG, 'Done with   ', msh_NZones, '  pre meshing zones.'

            end block


         endif ! ipre_mesh_type




         !===============================
         ! EXTERNAL PADDING
         !
         if (.not. warn_zone_over_limits .and. (settings%i_full_coverage_ == 1)) then

#ifdef _OPENMP
            if (allocated(zone_title)) deallocate(zone_title)
#endif

            ! BUG: check if it is ok setting this interest modes' pointer.
            call rz%setInterestModeIndexPtr(id_im_last)
            
            pol = MPolicy_PAD_ZONE_EXTERN
            call rz%setPolicy(pol)

            df_I = df_I_ref * pol%delta_fI_fct_
            df_J = df_J_ref * pol%delta_fJ_fct_

            allocate(maxext_sym_(4))  ! BUG: default allocate with max (4) elements

            if (settings%i_bisp_sym_ == BSA_SPATIAL_SYM_HALF) then
               basePts(1) = MPoint(0._bsa_real_t, maxF)
               basePts(2) = MPoint( maxF, maxF)
               basePts(3) = MPoint( maxF,-maxF)

               maxext_sym_(1) = max_ext
               maxext_sym_(2) = max_ext
               maxext_sym_(3) = 0._bsa_real_t

               n_dirs_ = N_DIRS_HALF
            else
               basePts(1) = MPoint(-maxF, maxF)
               basePts(2) = MPoint( maxF, maxF)
               basePts(3) = MPoint( maxF,-maxF)
               basePts(4) = MPoint(-maxF,-maxF)

               maxext_sym_(:) = max_ext
            endif

            !$omp parallel do  &
            !$omp   default(firstprivate),   &
            !$omp   private(zone_title),     &
            !$omp   shared(DIRS_LABELS, df_I, df_J, basePts, ROTATIONS  &
            !$omp          , struct_data, wd, settings, logger_debug    &
            !$omp          , LIM_SIGN_DIRS, LEFT_RZ_SIGNS, max_ext, maxext_sym_     &  
            !$omp          , NFREQS, NNODES, NNODESL, NLIBS, NLIBSL, NMODES, DIRS   &
            !$omp          , NMODES_EFF, MODES, NPSDEL, NTCOMPS, NDIRS, TCOMPS      &
            !$omp          , MSHR_SVD_INFO, MSHR_SVD_LWORK, MSHR_SVD_WORK           &
            !$omp          , msh_NZones, msh_bfmpts_pre_, msh_max_zone_NPts, m3mf_msh_ptr_), &
            !$omp   num_threads(n_dirs_)
            do idir = 1, n_dirs_

               call timer%init()
               zone_title = 'External rect Padding '//DIRS_LABELS(idir)
               call rz%zoneName(zone_title)

               call rz%setRotation(ROTATIONS(idir + 1))

               call rz%defineFromDeltas(basePts(idir), 'i', &
                  df_I, df_J, &
                  LIM_SIGN_DIRS(idir) * maxext_sym_(idir), &
                  LEFT_RZ_SIGNS(idir) * max_ext, force=.true.)

               call rz%compute()

#ifndef _BSA_USE_CACHED_POD_DATA
               !$omp critical
               call logger_debug%logZonePremeshingTotTime(&
                  zone_title, timer%time(), rz%zoneTotNPts(), .true.)
               !$omp end critical
#endif
            enddo ! ndirs
            !$omp end parallel do

            deallocate(maxext_sym_)

         endif ! (.not. warn_zone_over_limits .and. (settings%i_full_coverage_))

         print '(1x, a, a, i0, a/)', &
            INFOMSG, 'Done with   ', msh_NZones, '  pre meshing zones.'

      endblock



      998 continue


      ! NOTE: Dump modal info if needed, before rewinding..
      write(unit_dump_bfm_) settings%i_dump_modal_
      if (settings%i_dump_modal_ == 1) then

         print '(1x, a, a /)', &
            WARNMSG, 'Including modal info in dump file. Check if this can be avoided.'
         
         ! write kept modes, might serve after as well.
         ! NOTE: in fact, nm_eff_ is the VERY FIRST thing which is dumped!
         !       But, at the very beginning, we don't yet how many modes will be kept, 
         !       this is why it is done now here.
         write(unit_dump_bfm_) struct_data%modal_%modes_

         write(unit_dump_bfm_) &
            struct_data%modal_%Mm_(struct_data%modal_%modes_)
         write(unit_dump_bfm_) &
            struct_data%modal_%Cm_(struct_data%modal_%modes_, struct_data%modal_%modes_)
         write(unit_dump_bfm_) &
            struct_data%modal_%Km_(struct_data%modal_%modes_)

#ifdef _BSA_DEBUG
         print '(1x, a, a)', INFOMSG, 'Modal info dumped -- ok.'
#endif
      endif

      ! NOTE: Ok, now that premesh has finished, before going to actual meshing, 
      !       rewind dump file and rewrite actual needed head information.
      rewind(unit_dump_bfm_)
      write (unit_dump_bfm_) struct_data%modal_%nm_eff_
      write (unit_dump_bfm_) dimM_bisp_
      write (unit_dump_bfm_) msh_NZones
      write (unit_dump_bfm_) msh_max_zone_NPts

#ifdef _BSA_DEBUG
      write(unit_debug_, *) ' @BsaMesherImpl::PreMesh() : Init BSA-Mesher pre meshing phase -- ok.'
#endif
   end subroutine PreMesh















#ifdef _BSA_USE_CACHED_POD_DATA
# define __bfm_undump__
# define __bfm_undump_interp__
#else
# define __bfm_undump__  ,bfm_undump
# define __bfm_undump_interp__ bfm_undump,
#endif
   subroutine Mesh()
      !! Post meshing phase.
      !! Once data has been dumped from PreMeshing phase, 
      !! retrieve and process it.
      !! BFM data is interpolated based on interpolation method.
      !! Supported methods:
      !!    - HTPC : Head-Tail-Previous-Current
      use BsaLib_MZone, only: MZone_t, MZone_ID, UndumpZone
      integer(int32) :: izone_id, izone, ival2
      class(MZone_t), pointer     :: z => null()
      type(MRectZone_t), target   :: rz
      type(MTriangZone_t), target :: tz

#ifndef _BSA_USE_CACHED_POD_DATA
      real(bsa_real_t), allocatable :: bfm_undump(:, :)
# ifndef _OPENMP
      character(len = 128)          :: emsg
# endif
#endif

      class(*), pointer :: brm_export_data_ => null()
      type(BrmExportBaseData_t), allocatable, target :: brm_export_base_data_
      logical :: do_export_brm_base_
      

      ! skip them, we already have them stored in module variables
      ! However, there since they might serve outside this scope
      ! (i.e. if undumping file from another program)
      rewind(unit_dump_bfm_)        ! make sure to rewind unit before starting reading it
      read(unit_dump_bfm_) izone    ! n modes effective
      read(unit_dump_bfm_) izone    ! dimM_bisp
      read(unit_dump_bfm_) ival2    ! n. of meshing zones
      read(unit_dump_bfm_) ival2    ! msh_max_zone_PTS

#ifdef _BSA_DEBUG
      if (.not. izone == dimM_bisp_) &
         call bsa_Abort('First undumped value does not match BISP dimension.')
      if (.not. ival2 == msh_max_zone_NPts) &
         call bsa_Abort('Second undumped value has a wrong value.')
#endif


#ifndef _BSA_USE_CACHED_POD_DATA
# ifndef _OPENMP
      ! allocate BFM tmp variable to hold data for at most the zone with max n. of points.
      ! NOTE: in case of OMP parallelisation, this allocation will be resized internally if 
      !       needed. Otherwise, for serial case, at most  msh_max_zone_NPts  memory needed (once).
      allocate(bfm_undump(dimM_bisp_, msh_max_zone_NPts), stat=izone_id, errmsg=emsg)
      if (izone_id /= 0) call allocKOMsg('bfm_undump', izone_id, emsg)
# endif
#endif

      
      print '(1x, a)', '-----------------------------------------------------------'
      print '(1x, a)', '--------------------    POST - MESH    --------------------'
      print '(1x, a)', '-----------------------------------------------------------'

      do_export_brm_base_ = do_export_brm_ .and. i_brmexport_mode_ == BSA_EXPORT_BRM_MODE_BASE
      if (do_export_brm_base_) then
         allocate(brm_export_base_data_)
         brm_export_base_data_%nm_     =  struct_data%modal_%nm_eff_
         brm_export_base_data_%ncomb_  =  dimM_bisp_
         brm_export_base_data_%ispsym_ =  settings%i_bisp_sym_
         brm_export_base_data_%nzones_ =  msh_NZones
         brm_export_base_data_%modes_  => struct_data%modal_%modes_

         brm_export_data_ => brm_export_base_data_
      endif

      ! NOTE: Undump main Rect BKG Peak zone separately
      !
      read(unit_dump_bfm_) izone_id
      print '(1x, a, a, i6, a, i0 )', &
         INFOMSG, 'Interpolating zone n. ', 1, ', with ID=  ', izone_id
      if (do_export_brm_) brm_export_base_data_%idZone_ = izone_id
      call UndumpZone( rz   __bfm_undump__)
      call rz%interpolate(__bfm_undump_interp__   brm_export_data_)
      if (do_export_brm_base_) brm_export_base_data_%i_doNotPrintGenHeader_ = 1   ! now we can disable gen header print

#ifndef _BSA_USE_CACHED_POD_DATA
# ifdef _OPENMP
      deallocate(bfm_undump)  !<-- better to copy a null pointer than a whole bunch of memory.
      if (associated(brm_export_data_)) brm_export_data_ => null()
# endif
#endif


      ! NOTE: no need to check for EOF. We know how many zones we have dumped.
      !
      !$omp parallel do &
      !$omp   firstprivate(brm_export_base_data_, brm_export_data_), &
      !$omp   private(z, rz, tz, izone_id), &
#ifndef _BSA_USE_CACHED_POD_DATA
      !$omp   private(bfm_undump), &
#endif
      !$omp   shared(struct_data, wd, settings, logger_debug      &
      !$omp          , NFREQS, NNODES, NNODESL, NLIBS, NLIBSL        &
      !$omp          , NMODES, NMODES_EFF, MODES                     &
      !$omp          , NPSDEL, NTCOMPS, NDIRS, TCOMPS, DIRS          &
      !$omp          , MSHR_SVD_INFO, MSHR_SVD_LWORK, MSHR_SVD_WORK  &
      !$omp          , MZone_ID, msh_NZones, m3mr_msh_ptr_        & 
      !$omp          , msh_ZoneLimsInterestModes, do_validate_deltas_      &
      !$omp          , msh_bfmpts_post_, msh_brmpts_post_, unit_dump_bfm_  &
      !$omp          , bkg_peakw_ &
      !$omp          , dimM_bisp_, getBFM_msh, getBRM_msh, write_brm_fptr_),       &
      !$omp   num_threads(8)
      do izone = 2, msh_NZones

         !$omp critical
         read(unit_dump_bfm_) izone_id   ! fetch zone type ID

         print '(1x, a, a, i6, a, i0 )', &
            INFOMSG, 'Interpolating zone n. ', izone, ', with ID=  ', izone_id
            
         if (izone_id == MZone_ID%RECTANGLE) then
            call UndumpZone( rz   __bfm_undump__)
            z => rz
         elseif (izone_id == MZone_ID%TRIANGLE) then
            call UndumpZone( tz   __bfm_undump__)
            z => tz
         endif
         !$omp end critical

         ! brm_export_base_data_%i_doNotPrintZonHeader_ = 0
         if (do_export_brm_base_) then
            brm_export_base_data_%idZone_ = izone_id
            if (.not. associated(brm_export_data_)) brm_export_data_ => brm_export_base_data_
         endif
         call z%interpolate(__bfm_undump_interp__   brm_export_data_)
      enddo ! nZones
      !$omp end parallel do

#ifndef _BSA_USE_CACHED_POD_DATA
      if (allocated(bfm_undump)) deallocate(bfm_undump)
#endif
      99 return
   end subroutine Mesh
#undef __bfm_undump__
#undef __bfm_undump_interp__











   function getEquivalentLooperIterator_char(dim, pattern) result(LoopIter)
      integer, intent(in) :: dim
      character(len = *), intent(in) :: pattern
      character(len = 1)             :: LoopIter(dim)

      integer(int32) :: lpat, npat, ipat, id, i

      lpat = len(pattern)
      if (lpat > dim) call bsa_Abort('String length is greater than required Iterator length.')

      ! compute how many (integer) times pattern length stays inside dim
      npat = dim / lpat

      id = 1
      do ipat = 1, npat
         do i = 1, lpat
            LoopIter(id) = pattern(i:i)
            id = id + 1
         enddo
      enddo
      ! trailing
      do i = 1, dim - (lpat * npat)
         LoopIter(id) = pattern(i:i)
         id = id + 1
      enddo
   end function getEquivalentLooperIterator_char


   function getEquivalentLooperIterator_real(dim, vals) result(LoopIter)
      integer, intent(in) :: dim
      real(bsa_real_t), intent(in) :: vals(:)
      real(bsa_real_t)             :: LoopIter(dim)

      integer(int32) :: nvals, nint, i, j, id

      nvals = size(vals)
      if (nvals > dim) call bsa_Abort('Num of values greater than required Iterator length.')

      nint = dim / nvals

      id = 1
      do j = 1, nint
         do i = 1, nvals
            LoopIter(id) = vals(i)
            id = id + 1
         enddo
      enddo
      ! trailing
      do i = 1, dim - (nvals*nint)
         LoopIter(id) = vals(i)
         id = id + 1
      enddo
   end function getEquivalentLooperIterator_real






   pure elemental function getMaxSpaceExtension_() result(max_ext)
      real(bsa_real_t) :: max_ext
      
      max_ext = maxval(struct_data%modal_%nat_freqs_)
      max_ext = max_ext * settings%max_area_extension_
   end function getMaxSpaceExtension_






   subroutine prefetchZoneLimits_(bpw_ext_2, limits, policies, NLims, inter_modes)
      real(bsa_real_t), intent(in)  :: bpw_ext_2
      real(bsa_real_t), allocatable, intent(out) :: limits(:)
      type(MPolicy_t), allocatable, intent(out)  :: policies(:)
      integer(int32), intent(out)                :: NLims
      integer(int32), allocatable, intent(out)   :: inter_modes(:)

      integer(int32)   :: skip, imodesout
      real(bsa_real_t) :: peak_ext_lims_(2, NM_EFF__)


      ! get actual peak extensions, for each mode (BACK and FORTH limits)
      peak_ext_lims_ = getActualPeakZoneExtensionLimits_()

      ! search for modes that FALL (entirely) IN BKG PEAK AREA
      skip = 1
      do while (peak_ext_lims_(2, skip) <= bpw_ext_2)
         skip = skip + 1
         if (skip .gt. NM_EFF__) exit
      enddo
      ! skip = findloc(peak_ext_lims_(2, 1:NM_EFF__) > bpw_ext_2, .true.)
      ! if (skip == 0) then
      !    skip = NM_EFF__
      ! else
      !    skip = skip - 1
      ! endif
      skip = skip - 1     !<-- how many modes actually FULLY included

      ! NOTE: at worst, there will be 2*NM limits (clean case)
      imodesout = NM_EFF__ - skip

      if (imodesout == 0) then
         print '(1x, a, a)', WARNMSG, 'All resonant peak fall within BKG peak !'
         NLims = 0
         return
      endif

      print '(1x, a, i0, a, i0, a)', &
         INFOMSG, skip, '  res peak (out of ', NM_EFF__, ') fall(s) in BKG peak area.'


      ! warn if too many modes fall in BKG PEAK ZONE
      if (skip >= ceiling(real(NM_EFF__) / N_RES_PEAK_IN_BKG_ZONE_DIV_FCT_)) &
         print '(1x, 2a, i0, a/)', &
            WARNMSG, 'More than  1/', N_RES_PEAK_IN_BKG_ZONE_DIV_FCT_, &
            '  of resonant peaks fall entirely within BKG peak area.'


      block
         ! local instances, to move using move_alloc()
         real(bsa_real_t), allocatable :: limits_(:)
         type(MPolicy_t), allocatable  :: policies_(:)
         integer(int32)                :: NLims_
         integer(int32), allocatable   :: inter_modes_(:)

         integer(int32) :: itmp, iim, jim, im, istat, nmode
         character(len = 256) :: emsg

         ! peak zone's BACK and FORTH Frontiers
         real(bsa_real_t) :: pzBF, pzFF


         ! allocate results
         itmp = imodesout * 2 + 1
         allocate(limits_(itmp), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('limits_', istat, emsg)
         limits_ = 0._bsa_real_t

         allocate(policies_(itmp), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('policies_', istat, emsg)
         policies_(:) = MPolicy_DEF

         ! BUG: do not hard code dimension !!??
         !      Instead try to relate to NM and itmp
         !      Might throw run time error if we try to access
         !      out of bound!!
         allocate(inter_modes_(200), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('inter_modes_', istat, emsg)
         inter_modes_ = 0


         ! START
         
         skip = skip + 1  !<-- start from next (completely/partially out) resonant peak zone.
         if (bpw_ext_2 < peak_ext_lims_(1, skip)) then ! two separate peak zones -> PRE_PEAK in between

            ! we set first limits by DEFAULT.
            ! NOTE: they might be overridden from following res peak limit zones.
            limits_(1:2) = peak_ext_lims_(:, skip)
            NLims_       = 2
            policies_(1) = MPolicy_PRE_PEAK_1
            policies_(2) = MPolicy_PEAK

            ! Setting interest modes
            if (skip == 1) then
               inter_modes_(1) = CODE_PRE_PEAK_OK
            else
               inter_modes_(1) = CODE_PRE_PEAK_KO
            endif

            iim = 2
            inter_modes_(iim) = 1 ! this tells us how many interest modes to be read

            ! this are instead the real interest modes.
            jim = 3
            inter_modes_(jim) = struct_data%modal_%modes_(skip)

         else ! BKG extension falls in the MIDDLE of the PEAK ZONE (or also maybe after it!)

            limits_(1)   = peak_ext_lims_(2, skip)
            policies_(1) = MPolicy_PEAK
            NLims_       = 1

            ! setting interest modes, only PEAK zone
            inter_modes_(1) = 1
            inter_modes_(2) = struct_data%modal_%modes_(skip)
            iim = 1
            jim = 2
         endif

         ! now, continue starting from next mode's peak zone.
         skip = skip + 1
         do im = skip, NM_EFF__

            nmode = struct_data%modal_%modes_(im)

            pzBF = peak_ext_lims_(1, im) ! peak zone's BACK  FRONTIER
            pzFF = peak_ext_lims_(2, im) ! peak zone's FORTH FRONTIER

            if (pzBF < limits_(NLims_)) then ! there is OVERLAP!
                                             ! this mode's BACK frontier, falls BEHIND
                                             ! previous mode's FORTH frontier

               ! BUG: check this branch
               if (NLims_ == 2 .and. pzBF <= limits_(NLims_ - 1)) then  
                  
                  ! it is actually a FULL COVERAGE meaning that this mode's BACK frontier
                  ! entirely covers previously defined zone (i.e. resonance peak very very close)

                  itmp          = NLims_ - 1
                  limits_(itmp) = pzBF

                  policies_(itmp) = MPolicy_PRE_PEAK_2
                  
                  limits_(NLims_)   = pzFF
                  policies_(NLims_) = MPolicy_PEAK

               else  ! PARTIAL overlap
                     ! i.e. update only FORTH frontier, keeping BACK unchanged.

                  limits_(NLims_)   = pzFF
                  policies_(NLims_) = MPolicy_PEAK
               endif


               ! NOTE: since now, for PRE_PEAK zones we use "references"
               !       to interest modes of ADJACENT PEAK ZONES, 
               !       no need to update its list, but only for current peak zone.

               ! add an interest mode (NOTE: do not update its index)
               inter_modes_(iim) = inter_modes_(iim) + 1

               ! append current interest mode
               jim = jim + 1
               inter_modes_(jim) = nmode

            else ! NO COVERAGE at all, normal flow

               NLims_            = NLims_ + 1
               limits_(NLims_)   = pzBF
               policies_(NLims_) = MPolicy_PRE_PEAK_2

               NLims_            = NLims_ + 1
               limits_(NLims_)   = pzFF
               policies_(NLims_) = MPolicy_PEAK


               ! NOTE: in such case, a NEW INTEREST MODES COUNTING
               !       is initialised, since mode "im" and "im-1"
               !       are well distant and separated, so not to have 
               !       any interaction
               iim = jim + 1
               inter_modes_(iim) = CODE_PRE_PEAK_OK

               iim = iim + 1 ! this is for the new peak zone

               ! initialise counting to 1, might be incremented if overlapping to the right
               inter_modes_(iim) = 1

               jim = iim + 1
               inter_modes_(jim) = nmode
            endif
         enddo ! modes


         ! NOTE: appending interest mode (only last) for 
         !       that little padding zone added after last peak zone, for better shaping 
         !       (and for not losing important info)
         ! BUG: can be removed ??
         iim = jim + 1
         inter_modes_(iim) = 1
         jim = iim + 1
         inter_modes_(jim) = struct_data%modal_%modes_(NM_EFF__) ! only last mode is of interest


         ! move memory before leaving block
         NLims        = NLims_

         ! NOTE: this is needed because of covering after last peak.
         NLims_       = NLims_ + 1
         limits_      = limits_(1 : NLims_)
         call move_alloc(limits_, limits)


         inter_modes_ = inter_modes_(1 : jim) ! TODO: maybe do it here as well +1 ??
         call move_alloc(inter_modes_, inter_modes)


         allocate(policies(NLims_), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('policies', istat, emsg)
         policies = policies_(1 : NLims_)
      end block
   end subroutine ! prefetch zone limits




   function getActualPeakZoneExtensionLimits_() result(peak_exts_lims)
!DIR$ ATTRIBUTES FORCEINLINE :: getActualPeakZoneExtensionLimits_
      real(bsa_real_t) :: peak_exts_lims(2, NM_EFF__)
      real(bsa_real_t) :: cst, modf, ext
      integer(int32)   :: im, nmode
      integer(int32), parameter :: I_PEAK_EXT_DIV_ = 1

      cst = real(settings%gen_peak_area_extension_ / I_PEAK_EXT_DIV_, kind=bsa_real_t)

      allocate(peak_exts_(NM__))
      peak_exts_ = -1.
      do im = 1, NM_EFF__  ! TODO: implement do concurrent

         nmode = struct_data%modal_%modes_(im)
         
         modf = struct_data%modal_%nat_freqs_(nmode)
         ext  = struct_data%modal_%xsi_(nmode) * modf
         ext  = ext * cst

         peak_exts_lims(1, im) = modf - ext
         peak_exts_lims(2, im) = modf + ext

         peak_exts_(nmode) = ext + ext
      enddo
   end function


end submodule