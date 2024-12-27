!! This file is part of BsaLib.
!! Copyright (C) 2024  Michele Esposito Marzino 
!!
!! BsaLib is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BsaLib is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BsaLib.  If not, see <https://www.gnu.org/licenses/>.
submodule(BsaLib_Functions) BsaLib_FunctionsImpl

   use BsaLib_Utility
   use BsaLib_IO
   use BsaLib_CONSTANTS, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG, NOTEMSG
   use BsaLib_Data,      only: bsa_Abort  &
                              , do_trunc_POD_, POD_trunc_lim_     &
                              , do_export_POD_info_               &
                              , do_export_POD_trunc_              &
                              , nPODmodes_set_, nmodes_POD_
   implicit none (type, external)


#ifdef BSA_SINGLE_FLOATING_PRECISION
# ifdef BSA_USE_SVD_METHOD
#  define POD__ sgesvd
# else
#  define POD__ ssyev
# endif
#else
# ifdef BSA_USE_SVD_METHOD
#  define POD__ dgesvd
# else
#  define POD__ dsyev
# endif
#endif



   interface
#ifdef BSA_USE_SVD_METHOD
      SUBROUTINE POD__( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
                           VT, LDVT, WORK, LWORK, INFO )
         import bsa_real_t
         !    .. Scalar Arguments ..
         CHARACTER          JOBU, JOBVT
         INTEGER            LDA, LDU, LDVT, M, N
         integer            info, lwork
         !     .. Array Arguments ..
         real(bsa_real_t)   A( LDA, * ), S( * ), U( LDU, * ),&
                               VT( LDVT, * ), WORK( * )
      end subroutine
#else
      SUBROUTINE POD__( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
         import bsa_real_t
         !     .. Scalar Arguments ..
         CHARACTER          JOBZ, UPLO
         INTEGER            INFO, LDA, LWORK, N
         !     ..
         !     .. Array Arguments ..
         real(bsa_real_t)   A( LDA, * ), W( * ), WORK( * )
      END SUBROUTINE
#endif
   end interface



contains


   module subroutine setBsaFunctionLocalVars()

      NFREQS  = settings%nfreqs_

      ! nodal
      NNODES  = struct_data%nn_
      NNODESL = struct_data%nn_load_
      NLIBS   = struct_data%nlibs_    ! tot n. of LIBs per node
      NLIBSL  = struct_data%nlibs_load_    ! actual n. of loaded LIBs

      ! modal
      NMODES     = struct_data%modal_%nm_
      NMODES_EFF = struct_data%modal_%nm_eff_
      MODES      = struct_data%modal_%modes_

      ! wind
      NTCOMPS = wd%i_ntc_
      TCOMPS  = wd%tc_
      NDIRS   = wd%i_ndirs_
      DIRS    = wd%dirs_
      NPSDEL  = NNODESL * NTCOMPS * NDIRS
   end subroutine





   module subroutine prefetchSVDWorkDim_()
      real(bsa_real_t) :: tmpmat(NNODESL, NNODESL)
      ! real(bsa_real_t), allocatable :: tmpmat(:, :)

      real(bsa_real_t), dimension(NNODESL) :: tmpv

      real(bsa_real_t), dimension(1) :: optWork
      real(bsa_real_t), dimension(1) :: tmp1arr

      integer(int32) :: istat
      character(len = 256) :: emsg

      if (.not. allocated(MSHR_SVD_INFO)) then
         allocate(MSHR_SVD_INFO, stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('MSHR_SVD_INFO', istat, emsg)
      endif

      MSHR_SVD_INFO = 0
#ifdef BSA_USE_SVD_METHOD
   call POD__(&
           'O' &        ! min(M,N) columns of U are returned in array U
         , 'N' &        ! no rows of V are computed
         , NNODESL &    ! n. of rows M
         , NNODESL &    ! n. of cols N
         , tmpmat  &    ! A matrix
         , NNODESL &
         , tmpv    &
         , tmp1arr &    ! U array
         , 1       &
         , tmp1arr &
         , 1       &
         , optWork &
         , MSHR_SVD_LWORK &
         , MSHR_SVD_INFO  &
      )
#else
      call POD__('V', 'L', &
         NNODESL, tmpmat, NNODESL, tmp1arr, optWork, MSHR_SVD_LWORK, MSHR_SVD_INFO)
#endif

      if (MSHR_SVD_INFO == 0) then

         MSHR_SVD_LWORK = int(optWork(1))

#ifdef BSA_DEBUG
         print '(1x, 2a, i0 /)', &
            DBGMSG, 'WORK query ok. Optimal work dimension = ', MSHR_SVD_LWORK
#endif

         if (allocated(MSHR_SVD_WORK)) then
            if (size(MSHR_SVD_WORK) /= MSHR_SVD_LWORK) then
               deallocate(MSHR_SVD_WORK)
               allocate(MSHR_SVD_WORK(MSHR_SVD_LWORK), stat=istat, errmsg=emsg)
               if (istat /= 0) call allocKOMsg('MSHR_SVD_WORK', istat, emsg)
            endif
         else
            allocate(MSHR_SVD_WORK(MSHR_SVD_LWORK), stat=istat, errmsg=emsg)
            if (istat /= 0) call allocKOMsg('MSHR_SVD_WORK', istat, emsg)
         endif
         return ! correct execution flow
      endif


      print '(1x, 2a, i0)', &
         ERRMSG, 'WORK query for SVD decomposition returned code   ', MSHR_SVD_INFO
      print '(1x, 2a)', &
         MSGCONT, 'Please, check again.'
      call bsa_Abort()
   end subroutine




   module subroutine cleanSVDWorkInfo_()
      integer(int32) :: istat
      character(len = 256) :: emsg

      ! NOTE: reset to -1 so that next call is going to query again.
      MSHR_SVD_LWORK = - 1

      if (allocated(MSHR_SVD_INFO)) then
         deallocate(MSHR_SVD_INFO, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('MSHR_SVD_INFO', istat, emsg)
      endif

      if (allocated(MSHR_SVD_WORK)) then
         deallocate(MSHR_SVD_WORK, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('MSHR_SVD_WORK', istat, emsg)
      endif

#ifdef BSA_DEBUG
      print '(1x, 2a)', &
         INFOMSG, 'SVD related data cleaned -- ok.'
#endif
   end subroutine




   integer(int32) pure function getNPODModesByThreshold_(eigvals, rlim) result(nPODmodes)
#if (_WIN32 & __INTEL_COMPILER)
!DIR$ ATTRIBUTES FORCEINLINE :: getNPODModesByThreshold_
#endif
      real(bsa_real_t), intent(in), contiguous :: eigvals(:)
      real(bsa_real_t), intent(in) :: rlim
      real(bsa_real_t) :: limval

#ifdef BSA_USE_SVD_METHOD
# define __IDX 1
# define __OP  +
#else
# define __IDX NNODESL
# define __OP  -
#endif

      nPODmodes = 1
      if (any(eigvals /= eigvals(__IDX))) then
         limval    = rlim * eigvals(__IDX)
         nPODmodes = 2
         do while (eigvals(nPODmodes) >= limval)
            nPODmodes = nPODmodes __OP 1
         enddo
         nPODmodes = nPODmodes - 1
      endif
#undef __IDX
#undef __OP
   end function getNPODModesByThreshold_







!!========================================================================================
!!========================================================================================
!!
!!  mesher (implicit scalar)
!!
!!========================================================================================
!!========================================================================================


   module subroutine getFM_full_tnm_scalar_msh_(bfm, fi, fj)
      real(bsa_real_t), intent(inout), contiguous :: bfm(:, :)
      real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)

      real(bsa_real_t) :: fiPfj(1), abs_fi, abs_fj, abs_fiPfj

      ! indexes
      integer(int32) :: itc, tc, tcP3, iposM
      integer(int32) :: iposNK, iposNJ, iposNI
      integer(int32) :: ink, inj, ini
      integer(int32) :: nk, nj, ni
      integer(int32) :: posKi, posJi, posIi
      ! integer(int32) :: posKe, posJe, posIe
      integer(int32) :: imk, imj, imi
      integer(int32) :: ilk

      ! modal matrix slices
      real(bsa_real_t) :: phiJ(1, NLIBSL), phiI(NLIBSL, 1) !, phiK(1, 1, NLIBSL)
      real(bsa_real_t) :: phiK_(NMODES_EFF), phiK

      ! wind forces coeffs
      real(bsa_real_t), dimension(1, NLIBSL)    :: ajU, aj
      real(bsa_real_t), dimension(NLIBSL, 1)    :: aiU, ai, akU, ak

      ! basic PSDs
      real(bsa_real_t), dimension(1, NNODESL) :: S_IJK_fi, S_IJK_fj, S_IJK_fiPfj
      real(bsa_real_t) :: S_K_fi, S_K_fj, S_K_fiPfj
      real(bsa_real_t) :: S_J_fi, S_J_fj, S_J_fiPfj
      real(bsa_real_t) :: S_I_fi, S_I_fj, S_I_fiPfj

      ! nodal spactial correlations
      real(bsa_real_t) :: corrIJ, corrIK, corrJK

      ! crossed PSDs
      real(bsa_real_t) :: S_JK_fi, S_JK_fj
      real(bsa_real_t) :: term1, term2, term3, BF_IJK_ijk(NLIBSL, NLIBSL)
      real(bsa_real_t), dimension(NLIBSL, NLIBSL) :: tmp1, tmp2, tmp3

      ! frequencies values
      fiPfj(1)  = fi(1) + fj(1)
      abs_fi    = abs(fi(1))
      abs_fj    = abs(fj(1))
      abs_fiPfj = abs(fiPfj(1))


      ! NOTE: preinitalise to 0, to avoid uninitialised precision errors
      !       Like using memset() in C.
      bfm = 0._bsa_real_t

      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3

         ! prefetch wind turbulence PSDs for all loaded nodes
         S_IJK_fi = wd%evalPSD(1, fi, NNODESL, struct_data%n_load_, 1, tc)

         S_IJK_fj = wd%evalPSD(1, fj, NNODESL, struct_data%n_load_, 1, tc)

         S_IJK_fiPfj = wd%evalPSD(1, fiPfj, NNODESL, struct_data%n_load_, 1, tc)


         iposNK = 1
         do ink = 1, NNODESL

            nk = struct_data%n_load_(ink)

            ! BUG: must be this because of anoher bug ahead..
            posKi = (nk - 1) * NLIBS
            ! posKi = (nk - 1) * NLIBSL + 1
            ! posKe = nk * NLIBSL


            S_K_fi    = S_IJK_fi   (1, iposNK)
            S_K_fj    = S_IJK_fj   (1, iposNK)
            S_K_fiPfj = S_IJK_fiPfj(1, iposNK)


            ! NOTE: use node index and NOT values
            !       since we only store actually loaded ones.
            !       This applies also to actually loaded LIBs.
            !       Since we use it all -> : syntax
            ! akU(1, 1, :) = wd%wfc_(tc, :, ink)
            ! ak(1, 1, :)  = wd%wfc_(tcP3, :, ink)


            akU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ink)
            ak (:, 1) = wd%wfc_(struct_data%libs_load_, tcP3, ink)


            iposNJ = 1
            do inj = 1, NNODESL

               nj = struct_data%n_load_(inj)

               posJi = (nj - 1) * NLIBS
               ! posJi = (nj - 1) * NLIBSL + 1
               ! posJe = nj * NLIBSL

               S_J_fi    = S_IJK_fi   (1, iposNJ)
               S_J_fj    = S_IJK_fj   (1, iposNJ)
               S_J_fiPfj = S_IJK_fiPfj(1, iposNJ)

               ! BUG: should we put what??  DIR??  TC??
               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               ! NOTE: we can precompute for perf
               S_JK_fi = corrJK**(abs_fi) * sqrt(S_J_fi * S_K_fi)
               S_JK_fj = corrJK**(abs_fj) * sqrt(S_J_fj * S_K_fj)


               ajU(1, :) = wd%wfc_(struct_data%libs_load_, tc,   inj)
               aj (1, :) = wd%wfc_(struct_data%libs_load_, tcP3, inj)



               iposNI = 1
               do ini = 1, NNODESL

                  ni = struct_data%n_load_(ini)

                  posIi = (ni - 1) * NLIBS
                  ! posIi = (ni - 1) * NLIBSL + 1
                  ! posIe = ni * NLIBSL

                  S_I_fi    = S_IJK_fi   (1, iposNI)
                  S_I_fj    = S_IJK_fj   (1, iposNI)
                  S_I_fiPfj = S_IJK_fiPfj(1, iposNI)


                  aiU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ini)
                  ai (:, 1) = wd%wfc_(struct_data%libs_load_, tcP3, ini)


                  corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)
                  corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)


                  ! term1
                  term1 = corrIJ**(abs_fi) * sqrt(S_I_fi * S_J_fi)
                  term1 = corrIK**(abs_fj) * sqrt(S_I_fj * S_K_fj) * term1

                  tmp1 = matmul(ai, ajU) * term1


                  ! term2
                  term2 = corrIJ**(abs_fiPfj) * sqrt(S_I_fiPfj * S_J_fiPfj)
                  term2 = S_JK_fj * term2

                  tmp2 = matmul(aiU, aj) * term2


                  ! term3
                  term3 = corrIK**(abs_fiPfj) * sqrt(S_I_fiPfj * S_K_fiPfj)
                  term3 = S_JK_fi * term3

                  tmp3 = matmul(aiU, ajU) * term3


                  ! BUG: apparently, cannot make a 3D product..
                  do ilk = 1, NLIBSL


                     phiK_ = struct_data%modal_%phi_(posKi + struct_data%libs_load_(ilk), MODES)


                     ! BUG: this formulation DOES NOT account for
                     !      interaction between turbulent components
                     !      (i.e. uv, uw, vw)
                     BF_IJK_ijk = 2 * &
                        ( &
                           (tmp1 * akU(ilk, 1)) + &
                           (tmp2 * akU(ilk, 1)) + &
                           (tmp3 * ak (ilk, 1))   &
                        )


                     iposM = 1
                     do imk = 1, NMODES_EFF

                        phiK = phiK_(imk)

                        do imj = 1, NMODES_EFF

                           phiJ(1, :) = struct_data%modal_%phi_(posJi + struct_data%libs_load_, MODES(imj))

                           do imi = 1, NMODES_EFF

                              phiI(:, 1) = struct_data%modal_%phi_(posIi + struct_data%libs_load_, MODES(imi))

                              ! bfm(iposM) = bfm(iposM) + &
                              !    sum(BF_IJK_ijk(:, :, :) * (matmul(phiI, phiJ) * phiK(:, :, :)))

                              bfm(iposM, 1) = bfm(iposM, 1) + &
                                 sum(BF_IJK_ijk(:, :) * (matmul(phiI, phiJ) * phiK))

                              iposM = iposM + 1
                           enddo ! modes I

                        enddo ! modes J

                     enddo ! modes K


                  enddo ! lib K


                  iposNI = iposNI + 1
               enddo ! nodes I

               iposNJ = iposNJ + 1
            enddo ! nodes J

            iposNK = iposNK + 1
         enddo ! nodes K

      enddo ! itc

   end subroutine getFM_full_tnm_scalar_msh_








   module subroutine getFM_full_tm_scalar_msh_POD_(bfm, fi, fj)
#ifdef _OPENMP
      !$ use omp_lib, only: omp_get_thread_num
# define __export_POD_trunc_id__  omp_get_thread_num()+1
#else
# define __export_POD_trunc_id__  1
#endif
use BsaLib_Data, only: iun_POD_trunc_, do_export_POD_info_
      real(bsa_real_t), intent(inout), contiguous :: bfm(:, :)
      real(bsa_real_t), intent(in),    contiguous :: fi(:), fj(:)

#ifdef BSA_USE_POD_DATA_CACHING
# define __EGVL_w2 D_S_uvw_w2_ptr
# define __EGVT_w2 S_uvw_w2_ptr

      integer(int32)   :: nfi_, nfj_, ifi, ifj
      real(bsa_real_t) :: fi_
      integer(int32), allocatable   :: npodw2(:)

      real(bsa_real_t), allocatable, target  :: D_S_uvw_w2(:, :)
      real(bsa_real_t),              pointer :: D_S_uvw_w2_ptr(:) => null()

      real(bsa_real_t), allocatable, target  :: S_uvw_w2(:, :, :)
      real(bsa_real_t),              pointer :: S_uvw_w2_ptr(:, :) => null()
#else
# define __EGVL_w2 D_S_uvw_w2
# define __EGVT_w2 S_uvw_w2

      real(bsa_real_t)              :: D_S_uvw_w2(NNODESL)
      real(bsa_real_t), allocatable :: S_uvw_w2(:, :)
#endif

      real(bsa_real_t)              :: D_S_uvw_w1(NNODESL)   ! singular values vectors (DECREASING ordering!)
      real(bsa_real_t), allocatable :: S_uvw_w1(:, :)

      real(bsa_real_t)              :: D_S_uvw_w1w2(NNODESL)
      real(bsa_real_t), allocatable :: S_uvw_w1w2(:, :)

      ! for SVD related routines
      integer :: lwork, info
      real(bsa_real_t), allocatable :: work_arr(:)

      real(bsa_real_t) :: tmpv(1, NNODESL)   ! tmp vec for interfacing, BUG: might be avoided?

      real(bsa_real_t), dimension(NNODESL, 1)    :: eigvp, eigvq
      real(bsa_real_t), dimension(NMODES_EFF, 1) :: tmpm1, tmpm2, tmpm3
      real(bsa_real_t) :: tmpDp, tmpTq, tmpo, tmpn

      ! wind turbulent comps indexes
      integer(int32) :: itc, tc, tcP3

      ! n. of kept modes from wind fields decomposition
      integer(int32) :: nmw1, nmw2, nmw1w2
      integer(int32) :: p, q
      integer(int32) :: m, n, o, posm, posf

      real(bsa_real_t) :: fiPfj_(1)


      bfm = 0._bsa_real_t

#ifdef BSA_USE_POD_DATA_CACHING
      nfi_ = size(fi)
      nfj_ = size(fj)
      if (do_trunc_POD_) allocate(npodw2(nfj_))

      allocate(D_S_uvw_w2(NNODESL,          nfj_))
      allocate(S_uvw_w2  (NNODESL, NNODESL, nfj_))

#else
      allocate(S_uvw_w2  (NNODESL, NNODESL))
#endif
      allocate(S_uvw_w1  (NNODESL, NNODESL))
      allocate(S_uvw_w1w2(NNODESL, NNODESL))


      ! make a (private) copy of shared WORK array (avoid dependencies!)
      ! NOTE: at this point, MSHR_SVD_WORK should be allocated.
      work_arr = MSHR_SVD_WORK
      lwork    = MSHR_SVD_LWORK
      info     = 0


      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3

         posf = 1


#ifdef BSA_USE_POD_DATA_CACHING
         do ifj = 1, nfj_

            S_uvw_w2(:, 1:1, ifj) = &
               reshape(wd%evalPSD(1, fj(ifj:ifj), NNODESL, struct_data%n_load_, 1, tc), &
                  [NNODESL, 1])
            S_uvw_w2(:, :, ifj) = &
               wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w2(:, 1, ifj), fj(ifj), 1)

# ifdef BSA_USE_SVD_METHOD
            call POD__(&
                 'O' &                    ! min(M,N) columns of U are overwritten on array A (saves memory)
               , 'N' &                    ! no rows of V are computed
               , NNODESL    &             ! n. of rows M
               , NNODESL    &             ! n. of cols N
               , S_uvw_w2(:, :, ifj) &    ! A matrix (overwritten with left-singular vectors)
               , NNODESL             &
               , D_S_uvw_w2(:, ifj)  &    ! singular values
               , tmpv       &    ! U
               , 1          &
               , tmpv       &    ! VT
               , 1          &
               , work_arr, lwork, info)
# else
            call POD__('V', 'L', &
                        NNODESL, S_uvw_w2(:, :, ifj), NNODESL, D_S_uvw_w2(:, ifj), &
                           work_arr, lwork, info)
# endif
            if (info /= 0) then
               print '(1x, 2a, i0)', &
                  ERRMSG, 'Error applying SVD to S_uvw_w2. Exit code  ', info
               call bsa_Abort()
            endif

            if (do_trunc_POD_) npodw2(ifj) = getNPODModesByThreshold_(D_S_uvw_w2(:, ifj), POD_trunc_lim_)
         enddo ! nfj


#else  ! BSA_USE_POD_DATA_CACHING  not defined

         !
         ! NODAL WIND TURBULENCEs PSDs (for given tc)
         !
         S_uvw_w1(:, 1:1) = &
            reshape(wd%evalPSD(1, fi,   NNODESL, struct_data%n_load_, 1, tc), [NNODESL, 1])
         S_uvw_w2(:, 1:1) = &
            reshape(wd%evalPSD(1, fj,   NNODESL, struct_data%n_load_, 1, tc), [NNODESL, 1])

         !
         ! applying spatial nodal coherence
         !
         S_uvw_w1 = wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w1(:, 1), fi(1), 1)
         S_uvw_w2 = wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w2(:, 1), fj(1), 1)


# ifdef BSA_USE_SVD_METHOD
         call POD__(&
              'O' &           ! min(M,N) columns of U are overwritten on array A (saves memory)
            , 'N' &           ! no rows of V are computed
            , NNODESL    &    ! n. of rows M
            , NNODESL    &    ! n. of cols N
            , S_uvw_w1   &    ! A matrix (overwritten with left-singular vectors)
            , NNODESL    &
            , D_S_uvw_w1 &    ! singular values
            , tmpv       &    ! U
            , 1          &
            , tmpv       &    ! VT
            , 1          &
            , work_arr, lwork, info)
# else
         call POD__('V', 'L', &
                     NNODESL, S_uvw_w1, NNODESL, D_S_uvw_w1, &
                        work_arr, lwork, info)
# endif
         if (info /= 0) then
            print '(1x, 2a, i0)', &
               ERRMSG, 'Error applying SVD to S_uvw_w1. Exit code  ', info
            call bsa_Abort()
         endif


# ifdef BSA_USE_SVD_METHOD
         call POD__(&
              'O' &           ! min(M,N) columns of U are overwritten on array A (saves memory)
            , 'N' &           ! no rows of V are computed
            , NNODESL    &    ! n. of rows M
            , NNODESL    &    ! n. of cols N
            , S_uvw_w2   &    ! A matrix (overwritten with left-singular vectors)
            , NNODESL    &
            , D_S_uvw_w2 &    ! singular values
            , tmpv       &    ! U
            , 1          &
            , tmpv       &    ! VT
            , 1          &
            , work_arr, lwork, info)
# else
         call POD__('V', 'L', &
            NNODESL, S_uvw_w2, NNODESL, D_S_uvw_w2, &
               work_arr, lwork, info)
# endif
         if (info /= 0) then
            print '(1x, 2a, i0)', &
               ERRMSG, 'Error applying SVD to S_uvw_w2. Exit code  ', info
            call bsa_Abort()
         endif

#endif  ! BSA_USE_POD_DATA_CACHING



#ifdef BSA_USE_POD_DATA_CACHING

   if (.not. do_trunc_POD_) then
      if (nPODmodes_set_) then
         nmw1 = nmodes_POD_
         nmw2 = nmodes_POD_
      else
         nmw1 = NNODESL
         nmw2 = NNODESL
      endif
   endif


   do ifi = 1, nfi_

      S_uvw_w1(:, 1:1) = &
               reshape(wd%evalPSD(1, fi(ifi:ifi), NNODESL, struct_data%n_load_, 1, tc), &
                  [NNODESL, 1])

      S_uvw_w1 = &
         wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w1(:, 1), fi(ifi), 1)


# ifdef BSA_USE_SVD_METHOD
      call POD__(&
           'O' &           ! min(M,N) columns of U are overwritten on array A (saves memory)
         , 'N' &           ! no rows of V are computed
         , NNODESL    &    ! n. of rows M
         , NNODESL    &    ! n. of cols N
         , S_uvw_w1   &    ! A matrix (overwritten with left-singular vectors)
         , NNODESL    &
         , D_S_uvw_w1 &    ! singular values
         , tmpv       &    ! U
         , 1          &
         , tmpv       &    ! VT
         , 1          &
         , work_arr, lwork, info)
# else
      call POD__('V', 'L', &
         NNODESL, S_uvw_w1, NNODESL, D_S_uvw_w1, &
            work_arr, lwork, info)
# endif
      if (info /= 0) then
         print '(1x, 2a, i0)', &
            ERRMSG, 'Error applying SVD to S_uvw_w1. Exit code  ', info
         call bsa_Abort()
      endif

      if (do_trunc_POD_) nmw1 = getNPODModesByThreshold_(D_S_uvw_w1, POD_trunc_lim_)
      fi_ = fi(ifi)

      do ifj = 1, nfj_

         D_S_uvw_w2_ptr => D_S_uvw_w2(:, ifj)
         S_uvw_w2_ptr   => S_uvw_w2(:, :, ifj)
         if (do_trunc_POD_) nmw2 = npodw2(ifj)

         fiPfj_(1) = fi_ + fj(ifj)
#else  ! BSA_USE_POD_DATA_CACHING  not defined
         fiPfj_(1) = fi(1) + fj(1)
#endif ! BSA_USE_POD_DATA_CACHING

         S_uvw_w1w2(:, 1:1) = &
            reshape(wd%evalPSD(1, fiPfj_, NNODESL, struct_data%n_load_, 1, tc), [NNODESL, 1])
         S_uvw_w1w2 = wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w1w2(:, 1), fiPfj_(1), 1)

#ifdef BSA_USE_SVD_METHOD
         call POD__(&
              'O' &             ! min(M,N) columns of U are overwritten on array A (saves memory)
            , 'N' &             ! no rows of V are computed
            , NNODESL      &    ! n. of rows M
            , NNODESL      &    ! n. of cols N
            , S_uvw_w1w2   &    ! A matrix (overwritten with left-singular vectors)
            , NNODESL      &
            , D_S_uvw_w1w2 &    ! singular values
            , tmpv         &    ! U
            , 1            &
            , tmpv         &    ! VT
            , 1            &
            , work_arr, lwork, info)
#else
         call POD__('V', 'L', &
            NNODESL, S_uvw_w1w2, NNODESL, D_S_uvw_w1w2, &
               work_arr, lwork, info)
#endif ! BSA_USE_SVD_METHOD
         if (info /= 0) then
            print '(1x, 2a, i0)', &
               ERRMSG, 'Error applying SVD to S_uvw_w1w2. Exit code  ', info
            call bsa_Abort()
         endif


         if (do_trunc_POD_) then
#ifndef BSA_USE_POD_DATA_CACHING
            nmw1   = getNPODModesByThreshold_(D_S_uvw_w1, POD_trunc_lim_)
            nmw2   = getNPODModesByThreshold_(D_S_uvw_w2, POD_trunc_lim_)
#endif
            nmw1w2 = getNPODModesByThreshold_(D_S_uvw_w1w2, POD_trunc_lim_)
         else
            if (nPODmodes_set_) then
#ifndef BSA_USE_POD_DATA_CACHING
               nmw1   = nmodes_POD_
               nmw2   = nmodes_POD_
#endif
               nmw1w2 = nmodes_POD_
            else
#ifndef BSA_USE_POD_DATA_CACHING
               nmw1   = NNODESL
               nmw2   = NNODESL
#endif
               nmw1w2 = NNODESL
            endif
         endif


         if (do_export_POD_info_) then
            if (itc == 1 .and. do_export_POD_trunc_(__export_POD_trunc_id__)) then
               !$omp critical
               write(iun_POD_trunc_) int(__export_POD_trunc_id__, kind=int32)
               write(iun_POD_trunc_) size(fj)
               write(iun_POD_trunc_) real(fj, kind=real64)
               write(iun_POD_trunc_) int(size(D_S_uvw_w2), kind=int32)
               write(iun_POD_trunc_) D_S_uvw_w2
               !$omp end critical
            endif
         endif
         ! NOTE: we could place a barrier ?


         ! 5-2-2, 2-2-5 (6-3-3, 3-3-6) (7-4-4, 4-4-7)
#ifdef BSA_USE_SVD_METHOD
         do p = 1, nmw1
#else
         do p = NNODESL, NNODESL-(nmw1-1), -1
#endif

            eigvp(:, 1) = S_uvw_w1(:, p)

            ! V_lin_w1
            ! tmpm1 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvp)
            do q = 1, NMODES_EFF
               tmpm1(q, 1) = sum(wd%phi_times_A_ndegw_(q, :, tc) * eigvp(:, 1))
            enddo

            tmpDp = D_S_uvw_w1(p) ! D_p_w1


            ! 5-2-2 (6-3-3, 7-4-4)
#ifdef BSA_USE_SVD_METHOD
            do q = 1, nmw2
#else
            do q = NNODESL, NNODESL-(nmw2-1), -1
#endif

               eigvq(:, 1) = __EGVT_w2(:, q)

               ! VZ_quad_w1w2
               tmpm2 = matmul(wd%phi_times_A_ndegw_(:, :, tcP3), eigvp * eigvq)

               ! Z_lin_w2
               tmpm3 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvq)

               ! T_q_w2
               tmpTq = __EGVL_w2(q)


               posm = 1
               do o = 1, NMODES_EFF

                  tmpo = tmpm3(o, 1)

                  do n = 1, NMODES_EFF

                     tmpn = tmpm1(n, 1)

                     do m = 1, NMODES_EFF

                        bfm(posm, posf) = bfm(posm, posf) + &
                           (2 * tmpm2(m, 1) * &
                              tmpn * &
                              tmpo * tmpDp * tmpTq)

                        posm = posm + 1
                     enddo ! m modes
                  enddo ! n modes
               enddo ! o modes

            enddo ! q = 1, nmw2



            ! 2-2-5 (3-3-6, 4-4-7)
#ifdef BSA_USE_SVD_METHOD
            do q = 1, nmw1w2
#else
            do q = NNODESL, NNODESL-(nmw1w2-1), -1
#endif

               eigvq(:, 1) = S_uvw_w1w2(:, q)

               ! Z_lin_w1w2
               tmpm2 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvq)

               ! VZ_quad_w1w1w2
               tmpm3 = matmul(wd%phi_times_A_ndegw_(:, :, tcP3), eigvp * eigvq)

               ! T_q_w1w2
               tmpTq = D_S_uvw_w1w2(q)

               posm = 1
               do o = 1, NMODES_EFF

                  tmpo = tmpm3(o, 1)

                  do n = 1, NMODES_EFF

                     tmpn = tmpm1(n, 1)

                     do m = 1, NMODES_EFF

                        bfm(posm, posf) = bfm(posm, posf) + &
                           (2 * tmpm2(m, 1) * &
                              tmpn * &
                              tmpo * tmpDp * tmpTq)

                        posm = posm + 1
                     enddo ! m modes
                  enddo ! n modes
               enddo ! o modes

            enddo ! q = 1, nmw1w2
         enddo ! p = 1, nmw1


         ! 2-5-2 (3-6-3, 4-7-4)
#ifdef BSA_USE_SVD_METHOD
         do p = 1, nmw1w2
#else
         do p = NNODESL, NNODESL-(nmw1w2-1), -1
#endif

            eigvp(:, 1) = S_uvw_w1w2(:, p)

            ! V_lin_w1w2
            tmpm1 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvp)

            tmpDp = D_S_uvw_w1w2(p)

#ifdef BSA_USE_SVD_METHOD
            do q = 1, nmw2
#else
            do q = NNODESL, NNODESL-(nmw2-1), -1
#endif

               eigvq(:, 1) = __EGVT_w2(:, q)

               ! Z_lin_w2
               tmpm2 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvq)

               ! VZ_quad_w1w2w2
               tmpm3 = matmul(wd%phi_times_A_ndegw_(:, :, tcP3), eigvp * eigvq)

               tmpTq = __EGVL_w2(q)

               posm = 1
               do o = 1, NMODES_EFF

                  tmpo = tmpm2(o, 1)

                  do n = 1, NMODES_EFF

                     tmpn = tmpm3(n, 1)

                     do m = 1, NMODES_EFF

                        bfm(posm, posf) = bfm(posm, posf) + &
                           (2 * tmpm1(m, 1) * &
                              tmpn * &
                              tmpo * tmpDp * tmpTq)

                        posm = posm + 1
                     enddo ! m modes
                  enddo ! n modes
               enddo ! o modes

            enddo ! q = 1, nmw2
         enddo ! p = 1, nmw1w2


#ifdef BSA_USE_POD_DATA_CACHING
         posf = posf + 1
      enddo ! nfj_
   enddo ! nfi_
#endif

      enddo ! itc = 1, NTCOMPS
      99 return

#ifdef BSA_USE_POD_DATA_CACHING
      if (do_trunc_POD_) then
         if (allocated(npodw2)) deallocate(npodw2)
      endif
      if (allocated(D_S_uvw_w2)) deallocate(D_S_uvw_w2)
#endif
      if (allocated(S_uvw_w1))   deallocate(S_uvw_w1)
      if (allocated(S_uvw_w2))   deallocate(S_uvw_w2)
      if (allocated(S_uvw_w1w2)) deallocate(S_uvw_w1w2)
      if (allocated(work_arr))   deallocate(work_arr)
   end subroutine getFM_full_tm_scalar_msh_POD_











   module subroutine getRM_full_scalar_msh_(brm, fi, fj, bfm)
      real(bsa_real_t), intent(inout), contiguous :: brm(:, :)
      real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)
      real(bsa_real_t), intent(in), contiguous :: bfm(:, :)

      integer(int32)   :: i, j
      real(bsa_real_t) :: wi, wj, wiPwj
      integer(int32)   :: posm, posf, imk, imj, imi

      real(bsa_real_t), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(bsa_real_t), dimension(NMODES_EFF) :: H1r, H1i
      real(bsa_real_t), dimension(NMODES_EFF) :: H2r, H2i
      real(bsa_real_t), dimension(NMODES_EFF) :: H12r, H12i

      real(bsa_real_t) :: H12k_r, H12k_i, H2j_r, H2j_i


      posf = 1
      do i = 1, size(fi)

         wi = fi(i) * CST_PIt2

         do j = 1, size(fj)

            wj    = fj(j) * CST_PIt2
            wiPwj = wi + wj


            ! pre evaluate TFs (per mode)

            ! H1
            rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
            do imi = 1, NMODES_EFF
               Cdiag(imi) = struct_data%modal_%Cm_(MODES(imi), MODES(imi))
            enddo
            ipart = Cdiag * wi
            htmp  = rpart*rpart + ipart*ipart
            H1r   =   rpart / htmp
            H1i   = - ipart / htmp

            ! H2
            rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
            ipart = Cdiag * wj
            htmp  = rpart*rpart + ipart*ipart
            H2r   =   rpart / htmp
            H2i   = - ipart / htmp


            ! H12
            rpart = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
            ipart = Cdiag * wiPwj
            htmp  = rpart*rpart + ipart*ipart
            H12r  =   rpart / htmp
            H12i  = - ipart / htmp

            posm = 1
            do imk = 1, NMODES_EFF

               H12k_r = H12r(imk)
               H12k_i = H12i(imk)

               do imj = 1, NMODES_EFF

                  H2j_r = H2r(imj)
                  H2j_i = H2i(imj)

                  do imi = 1, NMODES_EFF

                     brm(posm, posf) = bfm(posm, posf) * &
                        (&
                           H1r(imi) * H2j_r * H12k_r + &
                           H1r(imi) * H2j_i * H12k_i + &
                           H1i(imi) * H2j_r * H12k_i - &
                           H1i(imi) * H2j_i * H12k_r   &
                        )

                     posm = posm + 1
                  enddo ! imi
               enddo ! imj
            enddo ! imk

            posf = posf + 1
         enddo
      enddo
   end subroutine getRM_full_scalar_msh_







   module subroutine getFM_diag_tnm_scalar_msh_(bfm, fi, fj)
      real(bsa_real_t), intent(inout), contiguous :: bfm(:, :)
      real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)

      real(bsa_real_t) :: fiPfj(1)

      integer(int32) :: itc, tc, tcP3, imode
      integer(int32) :: posi, inode, node, ilibk

      real(bsa_real_t), dimension(1, NNODESL) :: Suvw_fi, Suvw_fj, Suvw_fiPfj
      real(bsa_real_t), dimension(NNODESL) :: Suvw_IJ, Suvw_IJI, Suvw_IJJ

      real(bsa_real_t) :: akU, ak, phik(NMODES_EFF)
      real(bsa_real_t), dimension(NLIBSL, 1) :: aiU, ai
      real(bsa_real_t), dimension(1, NLIBSL) :: ajU, aj
      real(bsa_real_t), dimension(NLIBSL, NMODES_EFF) :: phi_

      real(bsa_real_t) :: BF_ijk_I(NLIBSL, NLIBSL)
      real(bsa_real_t), dimension(NLIBSL, NLIBSL) :: tmp1, tmp2, tmp3

      bfm = 0._bsa_real_t

      fiPfj(1) = fi(1) + fj(1)

      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3      ! quadratic term coeff

         Suvw_fi = wd%evalPSD(1, fi, NNODESL, struct_data%n_load_, 1, tc)

         Suvw_fj = wd%evalPSD(1, fj, NNODESL, struct_data%n_load_, 1, tc)

         Suvw_fiPfj = wd%evalPSD(1, fiPfj, NNODESL, struct_data%n_load_, 1, tc)


         ! precompute for perf
         Suvw_IJ  = Suvw_fi(1, :) * Suvw_fj(1, :)
         Suvw_IJI = Suvw_IJ(:) * Suvw_fi(1, :)
         Suvw_IJJ = Suvw_IJ(:) * Suvw_fj(1, :)


         do inode = 1, NNODESL

            node = int(struct_data%n_load_(inode), kind=int32)

            posi = (node - 1) * NLIBS
            phi_ = struct_data%modal_%phi_(posi + struct_data%libs_load_, MODES)


            ajU(1, :) = wd%wfc_(struct_data%libs_load_, tc,   inode)
            aj (1, :) = wd%wfc_(struct_data%libs_load_, tcP3, inode)

            aiU(:, 1) = ajU(1, :)
            ai (:, 1) = aj (1, :)


            ! NOTE: this are tmp values !!!!!
            !       Done for performance.
            tmp1 = Suvw_IJ (inode) * matmul(ai , ajU)
            tmp2 = Suvw_IJJ(inode) * matmul(aiU, aj )
            tmp3 = Suvw_IJI(inode) * matmul(aiU, ajU)


            do ilibk = 1, NLIBSL

               phik = phi_(ilibk, :)

               akU = aiU(ilibk, 1)
               ak  = ai (ilibk, 1)


               BF_ijk_I = 2 * (&
                  tmp1 * akU + &
                  tmp2 * akU + &
                  tmp3 * ak    &
               )


               do imode = 1, NMODES_EFF

                  bfm(imode, 1) = bfm(imode, 1) + &
                     sum( BF_ijk_I * &
                           (matmul(phi_(:, imode:imode), transpose(phi_(:, imode:imode))) &
                              * phik(imode) &
                           ) &
                     )
               enddo ! modes

            enddo ! libs loaded (k)
         enddo ! nodes loaded
      enddo ! n turb comps

   end subroutine getFM_diag_tnm_scalar_msh_







   module subroutine getRM_diag_scalar_msh_(brm, fi, fj, bfm)
      real(bsa_real_t), intent(inout), contiguous :: brm(:, :)
      real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)
      real(bsa_real_t), intent(in), contiguous :: bfm(:, :)

      integer(int32)     :: i, j, posf
      real(bsa_real_t)   :: wi, wj, wiPwj
      integer(bsa_int_t) :: imi

      real(bsa_real_t), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(bsa_real_t), dimension(NMODES_EFF) :: H1r, H1i
      real(bsa_real_t), dimension(NMODES_EFF) :: H2r, H2i
      real(bsa_real_t), dimension(NMODES_EFF) :: H12r, H12i

      posf = 1
      do i = 1, size(fi)

         wi = fi(i) * CST_PIt2

         do j = 1, size(fj)

            wj    = fj(j) * CST_PIt2
            wiPwj = wi + wj

            ! pre evaluate TFs (per mode)

            ! H1
            rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
            do imi = 1, NMODES_EFF
               Cdiag(imi) = struct_data%modal_%Cm_(MODES(imi), MODES(imi))
            enddo
            ipart = Cdiag * wi
            htmp  = rpart*rpart + ipart*ipart
            H1r   =   rpart / htmp
            H1i   = - ipart / htmp

            ! H2
            rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
            ipart = Cdiag * wj
            htmp  = rpart*rpart + ipart*ipart
            H2r   =   rpart / htmp
            H2i   = - ipart / htmp


            ! H12
            rpart  = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
            ipart  = Cdiag * wiPwj
            htmp   = rpart*rpart + ipart*ipart
            H12r   =   rpart / htmp
            H12i   = - ipart / htmp

            brm(:, posf) = bfm(:, posf) * (&
               H1r * H2r * H12r + &
               H1r * H2i * H12i + &
               H1i * H2r * H12i - &
               H1i * H2i * H12r   &
            )

            posf = posf + 1
         enddo
      enddo

   end subroutine getRM_diag_scalar_msh_











!!========================================================================================
!!========================================================================================
!!
!!  classic - VECTORISED
!!
!!========================================================================================
!!========================================================================================




   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   module subroutine getFM_full_tnm_vect_cls_(f, Suvw, psd, bisp)
      real(bsa_real_t), intent(in) :: f(NFREQS)
      real(bsa_real_t), intent(in) :: Suvw(NFREQS, NPSDEL)
      real(bsa_real_t), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)

      integer(int32) :: iin, ien, itmp, ifrj, i_cycle
#ifdef _OPENMP
      integer(int32) :: itmp_
#endif
      integer(int32) :: i_n_pad, i_pad_len

      ! turb components related
      integer(int32) :: itc, tc, tc_posN, tc_pk, tc_pj, tc_pi

      ! nodes indexed values
      integer(int32) :: pos_nk, pos_nj, pos_ni
      integer(int32) :: ink, inj, ini
      integer(int32) :: ni, nj, nk

      ! modes indexed values
      real(bsa_real_t), dimension(NMODES_EFF, 2) :: phik_, phij_, phii_
      real(bsa_real_t) :: phij_Ub_, phij_u_, phik_Ub_, phik_u_
      integer(int32) :: posm_
      integer(int32) :: imk, imj, imi

      real(bsa_real_t) :: f_abs(NFREQS)

      ! local nodal correlations
      real(bsa_real_t) :: corrJK, corrIK, corrIJ

      ! wfc extractions
      integer(bsa_int_t)   :: tcP3

      ! PSDs local
      real(bsa_real_t), allocatable :: S_uvw_i(:), S_uvw_j(:), S_uvw_k(:), PSDF_jk_JK_w(:)
      real(bsa_real_t), allocatable :: S_uvw_JK(:), S_uvw_IK(:), S_uvw_IJ(:)
      real(bsa_real_t), allocatable :: S_uvw_IK_w1w2(:), S_uvw_IJ_w1w2(:)

      ! BF local
      real(bsa_real_t), allocatable :: BF_ijk_IJK_w_w2(:), tmp1(:), tmp2(:), tmp3(:)

      character(len = 256) :: emsg
      !========================================================================

      f_abs = abs(f)

      ! getting padded length and relative init/end indices (non zero zone)
      itmp      = NFREQS - 1     ! do not consider 0 (point of symmetry)
      i_n_pad   = itmp / 2       ! spread it on the two sides (left / right)
      ien       = i_n_pad + NFREQS
      iin       = i_n_pad + 1
      i_pad_len = itmp + NFREQS


#ifdef BSA_DEBUG
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : i pad length = ', i_pad_len
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : init index   = ', iin
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : end  index   = ', ien
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : pad range    = ', ien - iin + 1
#endif


      ! these are needed regardlessly of if PSDs or BISPs
      allocate(psd(NFREQS, dimM_psd_), stat=itc, errmsg=emsg)
      if (itc /= 0) call allocKOMsg('psd', itc, emsg)
      psd = 0._bsa_real_t

      allocate(S_uvw_k(NFREQS), stat=itc, errmsg=emsg)
      if (itc /= 0) call allocKOMsg('S_uvw_k', itc, emsg)

      allocate(S_uvw_j(NFREQS), stat=itc, errmsg=emsg)
      if (itc /= 0) call allocKOMsg('S_uvw_j', itc, emsg)

      allocate(S_uvw_JK(NFREQS), stat=itc, errmsg=emsg)
      if (itc /= 0) call allocKOMsg('S_uvw_JK', itc, emsg)

      allocate(PSDF_jk_JK_w(NFREQS), stat=itc, errmsg=emsg)
      if (itc /= 0) call allocKOMsg('PSDF_jk_JK_w', itc, emsg)

      if (settings%i_compute_bisp_ == 1) then

         allocate(bisp(NFREQS, NFREQS, dimM_bisp_), stat=itc, errmsg=emsg)
         if (itc /= 0) call allocKOMsg('bisp', itc, emsg)
         bisp = 0._bsa_real_t

         allocate(bf_ijk_IJK_w_w2(NFREQS), stat=itc, errmsg=emsg)
         if (itc /= 0) call allocKOMsg('bf_ijk_IJK_w_w2', itc, emsg)

         allocate(S_uvw_i(NFREQS), stat=itc, errmsg=emsg)
         if (itc /= 0) call allocKOMsg('S_uvw_i', itc, emsg)

         allocate(S_uvw_IK(NFREQS), stat=itc, errmsg=emsg)
         if (itc /= 0) call allocKOMsg('S_uvw_IK', itc, emsg)

         allocate(S_uvw_IJ(NFREQS), stat=itc, errmsg=emsg)
         if (itc /= 0) call allocKOMsg('S_uvw_IJ', itc, emsg)

         allocate(S_uvw_IK_w1w2(i_pad_len), stat=itc, errmsg=emsg)
         if (itc /= 0) call allocKOMsg('S_uvw_IK_w1w2', itc, emsg)

         allocate(S_uvw_IJ_w1w2(i_pad_len), stat=itc, errmsg=emsg)
         if (itc /= 0) call allocKOMsg('S_uvw_IJ_w1w2', itc, emsg)
      endif ! i bisp allocation


      i_cycle = 0

      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, NTCOMPS

         tc      = wd%tc_(itc) ! get actual turbulent component
         tcP3    = tc + 3   ! quadratic term coeff
         tc_posN = (itc - 1) * NNODESL

         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + ink

            phik_(:, 1) = wd%phi_times_A_ndegw_(:, ink, tc  )
            phik_(:, 2) = wd%phi_times_A_ndegw_(:, ink, tcP3)

            S_uvw_k = Suvw(:, tc_pk)

            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + inj

               phij_(:, 1) = wd%phi_times_A_ndegw_(:, inj, tc  )
               phij_(:, 2) = wd%phi_times_A_ndegw_(:, inj, tcP3)

               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               S_uvw_j  = Suvw(:, tc_pj)
               S_uvw_JK = corrJK**(f_abs) * sqrt(S_uvw_k * S_uvw_j)


               !! BISPs
               if (settings%i_compute_bisp_ == 1) then

                  do ini = 1, NNODESL

                     ni     = struct_data%n_load_(ini)
                     pos_ni = (ni - 1) * NLIBS
                     tc_pi  = tc_posN + ini

                     phii_(:, 1) = wd%phi_times_A_ndegw_(:, ini, tc  )
                     phii_(:, 2) = wd%phi_times_A_ndegw_(:, ini, tcP3)

                     corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)
                     corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)

                     S_uvw_i  = Suvw(:, tc_pi)

                     S_uvw_IK = corrIK**(f_abs) * sqrt(S_uvw_i * S_uvw_k)
                     S_uvw_IK_w1w2(iin : ien) = S_uvw_IK

                     S_uvw_IJ = corrIJ**(f_abs) * sqrt(S_uvw_i * S_uvw_j)
                     S_uvw_IJ_w1w2(iin : ien) = S_uvw_IJ


                     ! loop on frequencies (second dimension, j)
#ifdef _OPENMP
                     itmp_ = NFREQS - 1

                     !$omp parallel do                &
                     !$omp   default(firstprivate),   &
                     !$omp   shared(bisp, itmp_, NMODES_EFF, NFREQS  &
                     !$omp      , S_uvw_IJ, S_uvw_IK, S_uvw_JK       &
                     !$omp      , S_uvw_IJ_w1w2, S_uvw_IK_w1w2       &
                     !$omp      , phik_, phij_, phii_),              &
                     !$omp   num_threads(8)
#else
                     itmp = NFREQS
#endif
                     do ifrj = 1, NFREQS

#ifdef _OPENMP
                        itmp = itmp_ + ifrj
#endif

                        ! make use of implicit allocation
                        ! NOTE: in OpenMP, better to copy a nullptr
                        tmp1 = S_uvw_IJ                   * S_uvw_IK(ifrj)
                        tmp2 = S_uvw_IJ_w1w2(ifrj : itmp) * S_uvw_JK(ifrj)
                        tmp3 = S_uvw_JK                   * S_uvw_IK_w1w2(ifrj : itmp)

                        posm_ = 1
                        do imk = 1, NMODES_EFF

                           phik_Ub_ = phik_(imk, 1)
                           phik_u_  = phik_(imk, 2)

                           do imj = 1, NMODES_EFF

                              phij_Ub_ = phij_(imj, 1)
                              phij_u_  = phij_(imj, 2)


                              ! BUG: this loop can be optimised??
                              !       Could we use SIMD lanes?
                              !       Maybe not, loop acts on non-contiguous memory...
                              do imi = 1, NMODES_EFF

                                 bisp(:, ifrj, posm_) = bisp(:, ifrj, posm_) + &
                                    2 * ( &
                                          (phii_(imi, 2) * phij_Ub_ * phik_Ub_ * tmp1)  &
                                       +  (phii_(imi, 1) * phij_u_  * phik_Ub_ * tmp2)  &
                                       +  (phii_(imi, 1) * phij_Ub_ * phik_u_  * tmp3)  &
                                    )

                                 posm_ = posm_ + 1
                              enddo ! i mode
                           enddo ! j mode
                        enddo ! k mode

#ifndef _OPENMP
                        itmp = itmp + 1
#endif
                     enddo ! n freqs j
#ifdef _OPENMP
                     !$omp end parallel do
#endif

                  enddo ! i node

               endif ! bisp computation


               posm_ = 1
               do imk = 1, NMODES_EFF

                  phik_Ub_ = phik_(imk, 1)

                  do imj = 1, NMODES_EFF

                     psd(:, posm_) = psd(:, posm_) + &
                        (phik_Ub_ * phij_(imj, 1) * S_uvw_JK)

                     posm_ = posm_ + 1
                  enddo ! j mode
               enddo ! k mode

            enddo ! j node

            i_cycle = i_cycle + 1
            print '(1x, a, " done ", f10.4, " %")', &
               INFOMSG, real(i_cycle, bsa_real_t)/real(NNODESL*NTCOMPS, bsa_real_t) * 100._bsa_real_t
        enddo ! k node
      enddo ! itc


      ! deallocation
      if (allocated(S_uvw_i))          deallocate(S_uvw_i)
      if (allocated(S_uvw_j))          deallocate(S_uvw_j)
      if (allocated(S_uvw_k))          deallocate(S_uvw_k)
      if (allocated(PSDF_jk_JK_w))     deallocate(PSDF_jk_JK_w)
      if (allocated(S_uvw_JK))         deallocate(S_uvw_JK)
      if (allocated(S_uvw_IK))         deallocate(S_uvw_IK)
      if (allocated(S_uvw_IJ))         deallocate(S_uvw_IJ)
      if (allocated(S_uvw_IK_w1w2))    deallocate(S_uvw_IK_w1w2)
      if (allocated(S_uvw_IJ_w1w2))    deallocate(S_uvw_IJ_w1w2)
      if (allocated(BF_ijk_IJK_w_w2))  deallocate(BF_ijk_IJK_w_w2)

#ifdef BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : computing modal forces spectra -- ok.'
#endif
   end subroutine getFM_full_tnm_vect_cls_










   module subroutine getRM_full_vect_cls_(f, psd, bisp)
      real(bsa_real_t), intent(in)                 :: f(NFREQS)
      real(bsa_real_t), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)

      integer(int32) :: ifrj
      integer(int32) :: posm_psd = 1, posm_bisp = 1

      ! modal indexed
      integer(int32) :: imi, imj, imk, mi

      real(bsa_real_t) :: omegas(NFREQS, 1)
      real(bsa_real_t), allocatable :: r_part(:, :), i_part(:, :), h_tmp(:, :), h_tmp2(:, :)
      real(bsa_real_t), allocatable :: Hr_w(:, :), Hi_w(:, :)
      real(bsa_real_t), allocatable :: Hr_w1w2(:, :, :), Hi_w1w2(:, :, :)
      !=====================================================================================

      ! BUG: check logic if correct
      if (settings%i_compute_bisp_ == 1) then

         allocate(r_part(NFREQS, NFREQS))
         allocate(i_part(NFREQS, NFREQS))
         allocate(h_tmp(NFREQS, NFREQS))
         allocate(h_tmp2(NFREQS, NFREQS))
         allocate(Hr_w1w2(NFREQS, NFREQS, NMODES_EFF))
         allocate(Hi_w1w2(NFREQS, NFREQS, NMODES_EFF))

      elseif (settings%i_compute_psd_ == 1) then

         allocate(r_part(NFREQS, 1))
         allocate(i_part(NFREQS, 1))
         allocate(h_tmp(NFREQS, 1))

      else ! none of them ??
         return
      endif

      ! these allocations are neede regardlessly
      allocate(Hr_w(NFREQS, NMODES_EFF))
      allocate(Hi_w(NFREQS, NMODES_EFF))

      omegas(:, 1) = f * CST_PIt2
      do imi = 1, NMODES_EFF

         mi = MODES(imi)

         ! H1 / H2
         r_part(:, 1) = - (omegas(:, 1)*omegas(:, 1) * struct_data%modal_%Mm_(mi)) + struct_data%modal_%Km_(mi)
         i_part(:, 1) = omegas(:, 1) * struct_data%modal_%Cm_(mi, mi)
         h_tmp(:, 1)  = r_part(:, 1)*r_part(:, 1) + i_part(:, 1)*i_part(:, 1)

         Hr_w(:, imi) =   r_part(:, 1) / h_tmp(:, 1)
         Hi_w(:, imi) = - i_part(:, 1) / h_tmp(:, 1)

         if (settings%i_compute_bisp_ == 1) then
            ! H12
#ifdef _OPENMP
            !$omp parallel do shared(NFREQS, omegas, h_tmp) num_threads(8)
#endif
            do ifrj = 1, NFREQS
               h_tmp(:, ifrj) = omegas(:, 1) + omegas(ifrj, 1)
            enddo
#ifdef _OPENMP
            !$omp end parallel do
#endif
            r_part = - (h_tmp*h_tmp * struct_data%modal_%Mm_(mi)) + struct_data%modal_%Km_(mi)
            i_part = h_tmp * struct_data%modal_%Cm_(mi, mi)
            h_tmp  = r_part*r_part + i_part*i_part

            Hr_w1w2(:, :, imi) =   r_part / h_tmp
            Hi_w1w2(:, :, imi) = - i_part / h_tmp
         endif
      enddo

      ! free mem as soon as possible
      if (allocated(r_part)) deallocate(r_part)
      if (allocated(i_part)) deallocate(i_part)


      ! BUG: optimise memory accesses !!!
      do imk = 1, NMODES_EFF

         if (settings%i_compute_bisp_ == 1) then
            h_tmp = Hr_w1w2(:, :, imk)
            h_tmp2= Hi_w1w2(:, :, imk)
         endif

         do imj = 1, NMODES_EFF

            psd(:, posm_psd) = &
               psd(:, posm_psd) * &
               (  Hr_w(:, imk) * Hr_w(:, imj) + &
                  Hi_w(:, imk) * Hi_w(:, imj) )

            ! BISPs
            if (settings%i_compute_bisp_ == 1) then

               do imi = 1, NMODES_EFF

#ifdef _OPENMP
                  !$omp parallel do &
                  !$omp   shared(NFREQS, posm_bisp, bisp    &
                  !$omp      , Hr_w, Hi_w, h_tmp, h_tmp2 ), &
                  !$omp   num_threads(8)
#endif
                  do ifrj = 1, NFREQS

                     bisp(:, ifrj, posm_bisp) = &
                        bisp(:, ifrj, posm_bisp) * &
                        (  Hr_w(:, imi) * Hr_w(ifrj, imj) * h_tmp(:, ifrj) + &
                           Hr_w(:, imi) * Hi_w(ifrj, imj) * h_tmp2(:, ifrj) + &
                           Hi_w(:, imi) * Hr_w(ifrj, imj) * h_tmp2(:, ifrj) - &
                           Hi_w(:, imi) * Hi_w(ifrj, imj) * h_tmp(:, ifrj) )
                  enddo
#ifdef _OPENMP
                  !$omp end parallel do
#endif

                  posm_bisp = posm_bisp + 1
               enddo ! i mode
            endif ! i bisp

            posm_psd = posm_psd + 1
         enddo ! j mode
      enddo ! k mode

      if (allocated(Hr_w)) deallocate(Hr_w)
      if (allocated(Hi_w)) deallocate(Hi_w)
      if (allocated(h_tmp)) deallocate(h_tmp)

      if (settings%i_compute_bisp_ == 1) then
         if (allocated(h_tmp2)) deallocate(h_tmp2)
         if (allocated(Hr_w1w2)) deallocate(Hr_w1w2)
         if (allocated(Hi_w1w2)) deallocate(Hi_w1w2)
      endif


#ifdef BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getRM_full_vect_cls_() : computing modal responses spectra -- ok.'
#endif
   end subroutine getRM_full_vect_cls_








   module subroutine getFM_diag_tnlm_vect_cls_(f, Suvw, psd, bisp)
      real(bsa_real_t), intent(in) :: f(NFREQS)
      real(bsa_real_t), intent(in) :: Suvw(NFREQS, NPSDEL)
      real(bsa_real_t), intent(inout), allocatable :: psd(:, :), bisp(:, :, :)

      integer(int32) :: iin = 0, ien = 0, itmp = 0, ifrj = 0
      integer(int32) :: i_n_pad = 0, i_pad_len = 0

      ! turb components related
      integer(int32) :: itc   = 0, tc_pos = 0, tc_posN = 0
      integer(int32) :: tc, tcP3

      ! node related indexes
      integer(int32) :: in, n, posN

      ! libs related vars
      integer(int32) :: ilk, ilj, ili, lk, lj, li

      ! modal matrix related
      real(bsa_real_t), dimension(NMODES_EFF) :: phik, phij, phii
      integer(int32) :: posk, posj, posi

      ! wind forces coeffs
      real(bsa_real_t) :: ai, aiU, aj, ajU, ak, akU

      ! modes related
      integer(int32) :: im

      !
      ! real(bsa_real_t) :: S_N_curr(settings%nfreqs_)
      real(bsa_real_t), allocatable :: Suvw_N_T(:, :), BF_ijk_III_w1w2(:, :)
      real(bsa_real_t), allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
      real(bsa_real_t), allocatable :: Suvw_N_w1w2(:, :), Suvw_N_pad(:)
      real(bsa_real_t), allocatable :: PSDF_jk_JJ_w(:)


#ifdef BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_diag_tnlm_vect_cls_() : computing modal forces spectra...'
#endif

      ! getting padded length and relative init/end indices (non zero zone)
      itmp      = NFREQS - 1    ! do not consider 0 (point of symmetry)
      i_n_pad   = itmp / 2                ! spread it on the two sides (left / right)
      ien       = i_n_pad + NFREQS
      iin       = i_n_pad + 1
      i_pad_len = itmp + NFREQS


      if (settings%i_compute_psd_ == 1) then

         allocate(PSDF_jk_JJ_w(NFREQS))
         allocate(psd(NFREQS, dimM_psd_))
         psd = 0._bsa_real_t
      endif

      if (settings%i_compute_bisp_ == 1) then

         allocate(Suvw_N_T(1, NFREQS))
         allocate(Suvw_N_pad(i_pad_len))
         allocate(Suvw_N_w1w2(NFREQS, NFREQS))

         allocate(tmp1(NFREQS, NFREQS))
         allocate(tmp2(NFREQS, NFREQS))
         allocate(tmp3(NFREQS, NFREQS))

         allocate(BF_ijk_III_w1w2(NFREQS, NFREQS))
         allocate(bisp(NFREQS, NFREQS, dimM_bisp_))
         bisp = 0._bsa_real_t
      endif


      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, NTCOMPS

         tc      = wd%tc_(itc)
         tcP3    = tc + 3
         tc_posN = (itc - 1) * NNODESL

         do in = 1, NNODESL

            n      = struct_data%n_load_(in)
            posN   = (n - 1) * NLIBS
            tc_pos = tc_posN + in


            ! BISP
            if (settings%i_compute_bisp_ == 1) then

               Suvw_N_pad(iin : ien) = Suvw(:, tc_pos)
               Suvw_N_T(1, :)        = Suvw(:, tc_pos)  ! storing transpose once for all

               itmp = 0

               do ifrj = 1, NFREQS

                  Suvw_N_w1w2(:, ifrj) = Suvw_N_pad(ifrj : NFREQS+itmp)

                  tmp1(:, ifrj) = Suvw_N_w1w2(:, ifrj) * Suvw_N_T(1, ifrj)
                  tmp2(:, ifrj) = Suvw_N_w1w2(:, ifrj) * Suvw(:, tc_pos)
                  ! tmp3(:, ifrj) = Suvw(:, tc_pos) * Suvw_N_T(1, ifrj)

                  itmp = itmp + 1
               enddo ! omega j

               ! original version
               ! NOTE: had stack overflow not using /heap-arrays0 compile option
               tmp3 = matmul(Suvw(:, tc_pos:tc_pos), Suvw_N_T)


               do ilk = 1, NLIBSL

                  lk   = struct_data%libs_load_(ilk)
                  posk = posN + lk
                  phik = struct_data%modal_%phi_(posk, MODES)

                  akU = wd%wfc_(lk, tc,   in)
                  ak  = wd%wfc_(lk, tcP3, in)


                  do ilj = 1, NLIBSL

                     lj   = struct_data%libs_load_(ilj)
                     posj = posN + lj
                     phij = struct_data%modal_%phi_(posj, MODES)

                     ajU = wd%wfc_(lj, tc,   in)
                     aj  = wd%wfc_(lj, tcP3, in)

                     do ili = 1, NLIBSL

                        li   = struct_data%libs_load_(ili)
                        posi = posN + li
                        phii = struct_data%modal_%phi_(posi, MODES)

                        aiU = wd%wfc_(li, tc,   in)
                        ai  = wd%wfc_(li, tcP3, in)


                        BF_ijk_III_w1w2 = 2 * (&
                           ai  * ajU * akU * tmp3 + &
                           aiU * aj  * akU * tmp1 + &
                           aiU * ajU * ak  * tmp2 &
                        &)


                        ! NOTE: no need to retrieve actual mode
                        !       since we are supposed to store only
                        !       kept mode for all related variables.
#ifdef BSA_DEBUG
                        do im = 1, NMODES_EFF
                           bisp(:, :, im) = bisp(:, :, im) + &
                              phik(im) * phij(im) * phii(im) * BF_ijk_III_w1w2
                        enddo
#else
                        do concurrent (im = 1:NMODES_EFF)
                           bisp(:, :, im) = bisp(:, :, im) + &
                              phik(im) * phij(im) * phii(im) * BF_ijk_III_w1w2
                        enddo
#endif

                     enddo ! i lib
                  enddo ! j lib
               enddo ! k lib

            endif ! do bisp computation


            ! PSDs
            do ilk = 1, NLIBSL

               lk   = struct_data%libs_load_(ilk)
               posk = posN + lk
               phik = struct_data%modal_%phi_(posk, MODES)

               akU = wd%wfc_(lk, tc, in)

               do ilj = 1, NLIBSL

                  lj   = struct_data%libs_load_(ilj)
                  posj = posN + lj
                  phij = struct_data%modal_%phi_(posj, MODES)

                  ajU = wd%wfc_(lj, tc, in)


                  PSDF_jk_JJ_w = ajU * akU * Suvw(:, tc_pos)

                  do im = 1, NMODES_EFF

                     psd(:, im) = psd(:, im) + &
                        phij(im) * phik(im) * PSDF_jk_JJ_w
                  enddo
               enddo ! j lib
            enddo ! k lib


#ifdef BSA_DEBUG
            print '(1x, 2a, f10.4, " %")', &
               INFOMSG, '  done  ', real(in, bsa_real_t) / NNODESL * 100
#endif
         enddo ! nodes
      enddo ! turb comps

      if (allocated(PSDF_jk_JJ_w))    deallocate(PSDF_jk_JJ_w)
      if (allocated(Suvw_N_T))        deallocate(Suvw_N_T)
      if (allocated(Suvw_N_pad))      deallocate(Suvw_N_pad)
      if (allocated(Suvw_N_w1w2))     deallocate(Suvw_N_w1w2)
      if (allocated(tmp1))            deallocate(tmp1)
      if (allocated(tmp2))            deallocate(tmp2)
      if (allocated(tmp3))            deallocate(tmp3)
      if (allocated(BF_ijk_III_w1w2)) deallocate(BF_ijk_III_w1w2)

#ifdef BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_diag_tnlm_vect_cls_() : computing modal forces spectra -- ok.'
#endif
   end subroutine






   module subroutine getRM_diag_vect_cls_(f, psd, bisp)
      real(bsa_real_t), intent(in)                 :: f(NFREQS)
      real(bsa_real_t), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)

      integer(int32) :: ifrj
      integer(int32) :: pos = 1

      ! modal indexed
      integer(int32) :: imi, mi

      real(bsa_real_t) :: omegas(NFREQS, 1)
      real(bsa_real_t) :: Mgi, Kgi, Cgi
      real(bsa_real_t), allocatable :: r_part(:, :), i_part(:, :), h_tmp(:, :)
      real(bsa_real_t), allocatable :: Hr_w(:), Hi_w(:)
      real(bsa_real_t), allocatable :: Hr_w1w2(:, :), Hi_w1w2(:, :)


#ifdef BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getRM_diag_vect_cls_() : computing modal responses spectra...'
#endif


      ! BUG: check logic if correct
      if (settings%i_compute_bisp_ == 1) then

         allocate(r_part(NFREQS, NFREQS))
         allocate(i_part(NFREQS, NFREQS))
         allocate(h_tmp(NFREQS, NFREQS))
         allocate(Hr_w1w2(NFREQS, NFREQS))
         allocate(Hi_w1w2(NFREQS, NFREQS))

      elseif (settings%i_compute_psd_ == 1) then

         allocate(r_part(NFREQS, 1))
         allocate(i_part(NFREQS, 1))
         allocate(h_tmp(NFREQS, 1))

      else ! none of them. NOTE: if so, we don't ever get here..
         return
      endif


      ! TRANSFER FUNCTION COMPUTATION

      ! pulsations
      omegas(:, 1) = f * CST_PIt2


      allocate(Hr_w(NFREQS))
      allocate(Hi_w(NFREQS))


      do imi = 1, NMODES_EFF

         mi = MODES(imi)

         ! NOTE: assumes we are storing only kept modes
         !       related modal info
         Mgi = struct_data%modal_%Mm_(mi)
         Kgi = struct_data%modal_%Km_(mi)
         Cgi = struct_data%modal_%Cm_(mi, mi)

         ! H1 / H2
         r_part(:, 1) = - (omegas(:, 1)*omegas(:, 1) * Mgi) + Kgi
         i_part(:, 1) = omegas(:, 1) * Cgi
         h_tmp(:, 1)  = r_part(:, 1)*r_part(:, 1) + i_part(:, 1)*i_part(:, 1)

         Hr_w =   r_part(:, 1) / h_tmp(:, 1)
         Hi_w = - i_part(:, 1) / h_tmp(:, 1)


         ! PSD
         psd(:, pos) = &
            psd(:, pos) * &
            (  Hr_w * Hr_w + &
               Hi_w * Hi_w )


         if (settings%i_compute_bisp_ == 1) then

            ! w1+w2^T
            do ifrj = 1, NFREQS
               h_tmp(:, ifrj) = omegas(:, 1) + omegas(ifrj, 1)
            enddo

            r_part = - (h_tmp*h_tmp * Mgi) + Kgi
            i_part = h_tmp * Cgi
            h_tmp  = r_part*r_part + i_part*i_part

            Hr_w1w2 =   r_part / h_tmp
            Hi_w1w2 = - i_part / h_tmp

            ! BISP
            do ifrj = 1, NFREQS

               bisp(:, ifrj, pos) = &
                  bisp(:, ifrj, pos) * &
                  (  Hr_w * Hr_w(ifrj) * Hr_w1w2(:, ifrj) + &
                     Hr_w * Hi_w(ifrj) * Hi_w1w2(:, ifrj) + &
                     Hi_w * Hr_w(ifrj) * Hi_w1w2(:, ifrj) - &
                     Hi_w * Hi_w(ifrj) * Hr_w1w2(:, ifrj) )
            enddo
         endif

         pos = pos + 1
      enddo ! i mode


      ! deallocate
      if (allocated(r_part)) deallocate(r_part)
      if (allocated(i_part)) deallocate(i_part)
      if (allocated(Hr_w))   deallocate(Hr_w)
      if (allocated(Hi_w))   deallocate(Hi_w)
      if (allocated(h_tmp))  deallocate(h_tmp)

      if (settings%i_compute_bisp_ == 1) then
         if (allocated(Hr_w1w2)) deallocate(Hr_w1w2)
         if (allocated(Hi_w1w2)) deallocate(Hi_w1w2)
      endif


#ifdef BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getRM_diag_vect_cls_() : computing modal responses spectra -- ok.'
#endif
   end subroutine getRM_diag_vect_cls_








!!========================================================================================
!!========================================================================================
!!
!!  classic - SCALAR
!!
!!========================================================================================
!!========================================================================================


   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   pure module subroutine getFM_full_tnm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
      integer(bsa_int_t), intent(in)  :: ii, ij
      real(bsa_real_t), intent(in)    :: fi, fj
      real(bsa_real_t), intent(in)    :: Suvw(NFREQS, NPSDEL)
      real(bsa_real_t), intent(in)    :: Suvw_pad(NPSDEL)
      real(bsa_real_t), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)

      ! turb components related
      integer(int32) :: itc, tc_posN, tc_pk, tc_pj

      ! freqs related
      real(bsa_real_t) :: fiabs, fjabs, fiPfj, fiPfjabs

      ! nodes indexed values
      integer(int32) :: pos_nk, pos_nj
      integer(int32) :: ink, inj
      integer(int32) :: nj, nk

      ! modes indexed values
      real(bsa_real_t), dimension(NMODES_EFF, 2) :: phik_, phij_
      real(bsa_real_t) :: phij_Ub_, phij_u_, phik_Ub_, phik_u_
      integer   :: posm
      integer   :: imk, imj, imi

      ! local nodal correlations
      real(bsa_real_t) :: corrJK

      ! wfc extractions
      integer(int32) :: tc, tcP3

      ! PSDs local
      real(bsa_real_t) :: S_uvw_j_i,  S_uvw_j_j,  S_uvw_j_ij
      real(bsa_real_t) :: S_uvw_k_i,  S_uvw_k_j,  S_uvw_k_ij
      real(bsa_real_t) :: S_uvw_JK_i, S_uvw_JK_j
      !========================================================================

      psd  = 0._bsa_real_t
      bisp = 0._bsa_real_t

      fiabs = abs(fi)
      fjabs = abs(fj)

      fiPfj    = fi + fj
      fiPfjabs = abs(fiPfj)


      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, wd%i_ntc_

         tc_posN = (itc - 1) * NNODESL
         tc      = wd%tc_(itc) ! get actual turbulent component
         tcP3    = tc + 3      ! quadratic term coeff


         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + ink

            phik_(:, 1) = wd%phi_times_A_ndegw_(:, ink, tc  )
            phik_(:, 2) = wd%phi_times_A_ndegw_(:, ink, tcP3)

            S_uvw_k_i  = Suvw(ii, tc_pk)
            S_uvw_k_j  = Suvw(ij, tc_pk)
            S_uvw_k_ij = Suvw_pad(tc_pk)


            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + inj

               phij_(:, 1) = wd%phi_times_A_ndegw_(:, inj, tc  )
               phij_(:, 2) = wd%phi_times_A_ndegw_(:, inj, tcP3)

               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               S_uvw_j_i  = Suvw(ii, tc_pj)
               S_uvw_j_j  = Suvw(ij, tc_pj)
               S_uvw_j_ij = Suvw_pad(tc_pj)


               S_uvw_JK_i = corrJK**(fiabs) * sqrt(S_uvw_j_i * S_uvw_k_i)
               S_uvw_JK_j = corrJK**(fjabs) * sqrt(S_uvw_j_j * S_uvw_k_j)


               !! BISPs

               if (settings%i_compute_bisp_ == 1) then

                  block
                     integer(int32) :: ini, ni, pos_ni
                     integer(int32) :: tc_pi

                     real(bsa_real_t) :: corrIK, corrIJ
                     real(bsa_real_t) :: S_uvw_i_i,  S_uvw_i_j,  S_uvw_i_ij
                     real(bsa_real_t) :: S_uvw_IK_j, S_uvw_IK_ij
                     real(bsa_real_t) :: S_uvw_IJ_i, S_uvw_IJ_ij
                     real(bsa_real_t) :: tmp1, tmp2, tmp3

                     real(bsa_real_t), dimension(NMODES_EFF, 2) :: phii_


                     do ini = 1, NNODESL

                        ni     = struct_data%n_load_(ini)
                        pos_ni = (ni - 1) * NLIBS
                        tc_pi  = tc_posN + ini

                        phii_(:, 1) = wd%phi_times_A_ndegw_(:, ini, tc  )
                        phii_(:, 2) = wd%phi_times_A_ndegw_(:, ini, tcP3)

                        corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)
                        corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)

                        S_uvw_i_i  = Suvw(ii, tc_pi)
                        S_uvw_i_j  = Suvw(ij, tc_pi)
                        S_uvw_i_ij = Suvw_pad(tc_pi)


                        S_uvw_IK_j  = corrIK**(fjabs)    * sqrt(S_uvw_i_j  * S_uvw_k_j )
                        S_uvw_IK_ij = corrIK**(fiPfjabs) * sqrt(S_uvw_i_ij * S_uvw_k_ij)

                        S_uvw_IJ_i  = corrIJ**(fiabs)    * sqrt(S_uvw_i_i  * S_uvw_j_i )
                        S_uvw_IJ_ij = corrIJ**(fiPfjabs) * sqrt(S_uvw_i_ij * S_uvw_j_ij)


                        tmp1 = S_uvw_IJ_i  * S_uvw_IK_j
                        tmp2 = S_uvw_IJ_ij * S_uvw_JK_j
                        tmp3 = S_uvw_JK_i  * S_uvw_IK_ij


                        posm = 1
                        do imk = 1, NMODES_EFF

                           phik_Ub_ = phik_(imk, 1)
                           phik_u_  = phik_(imk, 2)

                           do imj = 1, NMODES_EFF

                              phij_Ub_ = phij_(imj, 1)
                              phij_u_  = phij_(imj, 2)

                              do imi = 1, NMODES_EFF

                                 bisp(posm) = bisp(posm) +  &
                                    2 * ( &
                                          (phii_(imi, 2) * phij_Ub_ * phik_Ub_ * tmp1)  &
                                       +  (phii_(imi, 1) * phij_u_  * phik_Ub_ * tmp2)  &
                                       +  (phii_(imi, 1) * phij_Ub_ * phik_u_  * tmp3)  &
                                    )

                                 posm = posm + 1
                              enddo ! i mode
                           enddo ! j mode
                        enddo ! k mode

                     enddo ! i node

                  end block

               endif ! bisp computation


               !! PSDs
               if (ij == 1) then

                  posm = 1
                  do imk = 1, NMODES_EFF

                     phik_Ub_ = phik_(imk, 1)

                     do imj = 1, NMODES_EFF

                        psd(posm) = psd(posm) + &
                           (phik_Ub_ * phij_(imj, 1) * S_uvw_JK_i)

                        posm = posm + 1
                     enddo ! j mode
                  enddo ! k mode

               endif ! PSD computation


            enddo ! j node

         enddo ! k node

      enddo ! itc

   end subroutine getFM_full_tnm_scalar_cls_










   module subroutine getRM_full_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
      integer(bsa_int_t), intent(in)  :: ii, ij
      real(bsa_real_t), intent(in)    :: fi, fj
      real(bsa_real_t), intent(in)    :: psdin(dimM_psd_), bispin(dimM_bisp_)
      real(bsa_real_t), intent(out)   :: psdout(dimM_psd_), bispout(dimM_bisp_)

      real(bsa_real_t) :: wi
      integer(int32) :: posm, imk, imj

      real(bsa_real_t), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(bsa_real_t), dimension(NMODES_EFF) :: H1r, H1i

      real(bsa_real_t) :: H2j_r, H2j_i


      wi = fi * CST_PIt2

      ! pre evaluate TFs (per mode)

      ! H1
      rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      do imk = 1, NMODES_EFF
         Cdiag(imk) = struct_data%modal_%Cm_(MODES(imk), MODES(imk))
      enddo
      ipart = Cdiag * wi
      htmp  = rpart*rpart + ipart*ipart
      H1r   =   rpart / htmp
      H1i   = - ipart / htmp


      if (ij == 1) then ! also PSDr

         posm = 1
         do imk = 1, NMODES_EFF

            H2j_r = H1r(imk)
            H2j_i = H1i(imk)

            do imj = 1, NMODES_EFF

               psdout(posm) = psdin(posm) * (&
                  H1r(imj) * H2j_r + &
                  H1i(imj) * H2j_i   &
               )

               posm = posm + 1
            enddo ! j modes
         enddo ! k modes

         if (settings%i_compute_bisp_ == 0) return
      endif


      ! BISPr

      block
         real(bsa_real_t), dimension(NMODES_EFF) :: H2r, H2i
         real(bsa_real_t), dimension(NMODES_EFF) :: H12r, H12i
         real(bsa_real_t) :: H12k_r, H12k_i
         real(bsa_real_t) :: wj, wiPwj, tmp1, tmp2, tmp3, tmp4
         integer   :: imi

         wj    = fj * CST_PIt2
         wiPwj = wi + wj

         ! H2
         rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart = Cdiag * wj
         htmp  = rpart*rpart + ipart*ipart
         H2r   =   rpart / htmp
         H2i   = - ipart / htmp


         ! H12
         rpart  = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart  = Cdiag * wiPwj
         htmp   = rpart*rpart + ipart*ipart
         H12r   =   rpart / htmp
         H12i   = - ipart / htmp


         posm = 1
         do imk = 1, NMODES_EFF

            H12k_r = H12r(imk)
            H12k_i = H12i(imk)

            do imj = 1, NMODES_EFF

               H2j_r = H2r(imj)
               H2j_i = H2i(imj)

               tmp1 = H2j_r * H12k_r  ! real
               tmp2 = H2j_i * H12k_i  ! real
               tmp3 = H2j_r * H12k_i  ! imag
               tmp4 = H2j_i * H12k_r  ! imag

               do imi = 1, NMODES_EFF

                  bispout(posm) = bispin(posm) * (&
                     H1r(imi) * tmp1 + &
                     H1r(imi) * tmp2 + &
                     H1i(imi) * tmp3 - &
                     H1i(imi) * tmp4   &
                  )

                  posm = posm + 1
               enddo ! i modes
            enddo ! j modes
         enddo ! k modes

      end block

   end subroutine getRM_full_scalar_cls_














   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   pure module subroutine getFM_diag_tnlm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
      integer(bsa_int_t), intent(in)  :: ii, ij
      real(bsa_real_t), intent(in)    :: fi, fj
      real(bsa_real_t), intent(in)    :: Suvw(NFREQS, NPSDEL)
      real(bsa_real_t), intent(in)    :: Suvw_pad(NPSDEL)
      real(bsa_real_t), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)

      integer(int32)   :: itc, tc, tcP3, iposN, inode, n, imode, posNi
      real(bsa_real_t) :: Suvw_i

      real(bsa_real_t), dimension(NLIBSL, NMODES_EFF) :: phiN_
      real(bsa_real_t), dimension(NLIBSL, NLIBSL) :: PSDF_jk_JJ_w, tmp3, phiIJ_
      real(bsa_real_t), dimension(NLIBSL) :: aNU, aN


      psd  = 0._bsa_real_t
      bisp = 0._bsa_real_t


      iposN = 1
      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3


         do inode = 1, NNODESL

            n     = struct_data%n_load_(inode)
            posNi = (n - 1) * NLIBS

            phiN_ = struct_data%modal_%phi_(posNi + struct_data%libs_load_, MODES)

            aNU   = wd%wfc_(struct_data%libs_load_, tc,   inode)
            an    = wd%wfc_(struct_data%libs_load_, tcP3, inode)

            Suvw_i = Suvw(ii, iposN)

            ! PSD comp
            if (ij == 1) then

               do imode = 1, NLIBSL
                  tmp3(:, imode) = aNU * aNU(imode)
               enddo
               PSDF_jk_JJ_w = tmp3 * Suvw_i

               do imode = 1, NMODES_EFF

                  phiIJ_ = matmul(phiN_(:, imode:imode), transpose(phiN_(:, imode:imode)))

                  psd(imode) = psd(imode) + &
                     sum(phiIJ_ * PSDF_jk_JJ_w)
               enddo ! n modes
            endif ! ij == 1


            if (settings%i_compute_bisp_ == 1) then

               block
                  real(bsa_real_t) :: Suvw_j, Suvw_ij
                  integer   :: ilk

                  real(bsa_real_t) :: BF_ijk_III_wiwj(NLIBSL, NLIBSL)
                  real(bsa_real_t), dimension(NLIBSL, NLIBSL) :: tmp2, tmp1

                  Suvw_j  = Suvw(ij, iposN)
                  Suvw_ij = Suvw_pad(iposN)

                  tmp1 = matmul(reshape(aN , [NLIBSL, 1]), reshape(aNU, [1, NLIBSL]))
                  tmp2 = matmul(reshape(aNU, [NLIBSL, 1]), reshape(aN , [1, NLIBSL]))
                  tmp3 = matmul(reshape(aNU, [NLIBSL, 1]), reshape(aNU, [1, NLIBSL]))

                  tmp1 = tmp1 * (Suvw_i  * Suvw_j)
                  tmp2 = tmp2 * (Suvw_ij * Suvw_j)
                  tmp3 = tmp3 * (Suvw_ij * Suvw_i)

                  do ilk = 1, NLIBSL

                     BF_ijk_III_wiwj = 2 * (&
                        tmp1 * aNU(ilk) + &
                        tmp2 * aNU(ilk) + &
                        tmp3 * aN (ilk)   &
                     )

                     do imode = 1, NMODES_EFF

                        phiIJ_ = matmul(phiN_(:, imode:imode), transpose(phiN_(:, imode:imode)))

                        bisp(imode) = bisp(imode) + &
                           sum(phiIJ_ * BF_ijk_III_wiwj * phiN_(ilk, imode))
                     enddo ! n modes

                  enddo ! k libs

               end block

            endif ! compute BISP

            iposN = iposN + 1
         enddo ! nodes loaded

      enddo ! turb comps

   end subroutine getFM_diag_tnlm_scalar_cls_






   module subroutine getRM_diag_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
      integer(bsa_int_t), intent(in)  :: ii, ij  ! freqs indexes
      real(bsa_real_t), intent(in)    :: fi, fj
      real(bsa_real_t), intent(in)    :: psdin(dimM_psd_), bispin(dimM_bisp_)
      real(bsa_real_t), intent(out)   :: psdout(dimM_psd_), bispout(dimM_bisp_)

      real(bsa_real_t)   :: wi
      integer(bsa_int_t) :: imi

      real(bsa_real_t), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(bsa_real_t), dimension(NMODES_EFF) :: H1r, H1i


      wi = fi * CST_PIt2


      ! pre evaluate TFs (per mode)

      ! H1
      rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      do imi = 1, NMODES_EFF
         Cdiag(imi) = struct_data%modal_%Cm_(MODES(imi), MODES(imi))
      enddo
      ipart = Cdiag * wi
      htmp  = rpart*rpart + ipart*ipart
      H1r   =   rpart / htmp
      H1i   = - ipart / htmp


      if (ij == 1) psdout = psdin * (H1r * H1r + H1i * H1i)


      block
         real(bsa_real_t) :: wj, wiPwj
         real(bsa_real_t), dimension(NMODES_EFF) :: H2r, H2i
         real(bsa_real_t), dimension(NMODES_EFF) :: H12r, H12i
         real(bsa_real_t) :: H12k_r, H12k_i, H2j_r, H2j_i


         wj    = fj * CST_PIt2
         wiPwj = wi + wj

         ! H2
         rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart = Cdiag * wj
         htmp  = rpart*rpart + ipart*ipart
         H2r   =   rpart / htmp
         H2i   = - ipart / htmp


         ! H12
         rpart = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart = Cdiag * wiPwj
         htmp  = rpart*rpart + ipart*ipart
         H12r  =   rpart / htmp
         H12i  = - ipart / htmp


         bispout = bispin * (&
            H1r * H2r * H12r + &
            H1r * H2i * H12i + &
            H1i * H2r * H12i - &
            H1i * H2i * H12r   &
         )

      end block

   end subroutine getRM_diag_scalar_cls_












   pure module subroutine getBR_SFm_val_(nm, Suvw, fnat, im, m, psd)
      !! BUG: very unoptimised..
      !!      is basically a small copy of "getFM_full_tnlm_scalar_cls_"
      integer(bsa_int_t), intent(in)  :: im, m, nm
      real(bsa_real_t), intent(in)    :: Suvw(nm, NPSDEL), fnat
      real(bsa_real_t), intent(inout) :: psd

      ! turb components related
      integer(int32) :: itc, tc_posN, tc_pk, tc_pj

      ! freqs related
      real(bsa_real_t) :: fiabs

      ! nodes indexed values
      integer(int32) :: i_pos_nk, i_pos_nj
      integer(int32) :: pos_nk, pos_nj
      integer(int32) :: ink, inj
      integer(int32) :: nj, nk

      ! local nodal correlations
      real(bsa_real_t) :: corrJK

      integer(int32) :: tc
      real(bsa_real_t), dimension(NLIBSL, 1) :: akU, phij_
      real(bsa_real_t), dimension(1, NLIBSL) :: ajU, phik_

      ! PSDs local
      real(bsa_real_t) :: S_uvw_j_i
      real(bsa_real_t) :: S_uvw_k_i
      real(bsa_real_t) :: S_uvw_JK_i

      real(bsa_real_t), dimension(NLIBSL, NLIBSL) :: PSD_jk_JK_w
      !========================================================================

      psd   = 0._bsa_real_t
      fiabs = abs(fnat)


      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, wd%i_ntc_

         tc_posN = (itc - 1) * NNODESL
         tc      = wd%tc_(itc) ! get actual turbulent component

         i_pos_nk = 1
         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + i_pos_nk

            phik_  = transpose(struct_data%modal_%phi_(pos_nk + struct_data%libs_load_, m:m))

            akU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ink)

            S_uvw_k_i = Suvw(im, tc_pk)



            i_pos_nj = 1
            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + i_pos_nj

               phij_  = struct_data%modal_%phi_(pos_nj + struct_data%libs_load_, m:m)

               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               S_uvw_j_i  = Suvw(im, tc_pj)
               S_uvw_JK_i = corrJK**(fiabs) * sqrt(S_uvw_j_i * S_uvw_k_i)


               ajU(1, :) = wd%wfc_(struct_data%libs_load_, tc,   inj)


               ! PSD f
               PSD_jk_JK_w = matmul(akU, ajU) * S_uvw_JK_i

               psd = psd + sum(matmul(phij_, phik_) * PSD_jk_JK_w)

               i_pos_nj = i_pos_nj + 1
            enddo ! j node

            i_pos_nk = i_pos_nk + 1
         enddo ! k node

      enddo ! itc
   end subroutine



end submodule BsaLib_FunctionsImpl
