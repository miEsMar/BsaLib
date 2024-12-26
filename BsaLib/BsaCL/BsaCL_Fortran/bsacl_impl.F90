submodule(BsaCL) BsaCL_impl

#include "_msgtypes.h"

   use, intrinsic :: iso_c_binding
   use BsaCL_c
   implicit none


#ifdef BSACL_ENABLE_EVALFCT_PTR
   procedure(bsacl_EvalFunc_C_to_Fortran__), pointer :: evalfct_ptr_ => null()
#endif


contains


   module integer(IK) function bsacl_Init(n_threads)
      integer(IK), intent(in), value :: n_threads
      bsacl_Init = int( bsaclInit__(int(n_threads, c_int)), IK )
   end function


   module integer(IK) function bsacl_InitDeviceMemory()
      bsacl_InitDeviceMemory = int( bsaclInitDeviceMemory__(), IK )
   end function


   module integer(IK) function bsacl_Run(i_thread, res)
      integer(IK), intent(in), value :: i_thread
      real(RK), pointer, intent(in)  :: res(:)
      bsacl_Run = int( bsaclRun__(int(i_thread, c_int), int(size(res), c_size_t), c_loc(res)), IK )
   end function


   module function bsacl_SetKernelID(kid) result(ierr)
      integer(IK), value, intent(in) :: kid
      integer(IK) :: ierr

      ierr = int(bsaclSetKernelID__(int(kid, c_int)))
   end function



   module subroutine bsacl_AcquirePSDId(psdid)
      integer(IK), intent(in) :: psdid
      call bsaclAcquirePSDId__(int(psdid, c_int))
   end subroutine


   module subroutine bsacl_AcquireStructModMat(modmat, natf)
      real(RK), intent(in), target :: modmat(:, :), natf(:)
      integer(c_int) :: ndofs_c_, nmodes_c_

      ndofs_c_  = size(modmat, dim=1, kind=c_int)
      nmodes_c_ = size(modmat, dim=2, kind=c_int)
      if (size(natf, kind=c_int) /= nmodes_c_) then
         print '(// 1x, 2a)', INFO_MSG, "N. of modes does not match between modal matrix and nat freqs vector."
         call bsaclAbort__(-200)
      endif
      call bsaclAcquireStructModMat__(c_loc(modmat), c_loc(natf), ndofs_c_, nmodes_c_)
   end subroutine



   module subroutine bsacl_AcquireModalMatrices(Mg, Cg, Kg)
      real(RK), intent(in), dimension(:),    target :: Mg, Kg
      real(RK), intent(in), dimension(:, :), target :: Cg

      call bsaclAcquireModalMatrices__(c_loc(Mg), c_loc(Cg), c_loc(Kg))
   end subroutine





   module subroutine bsacl_AcquireLoadedNodesList(nodes_load)
      integer(IK), intent(in), target :: nodes_load(:)
      integer(c_int) :: nnodes_l_
      nnodes_l_ = size(nodes_load, kind=c_int)
      call bsaclAcquireLoadedNodesList__(c_loc(nodes_load), nnodes_l_)
   end subroutine


   module subroutine bsacl_AcquireTotalNOfNodes(nn)
      integer(IK), intent(in) :: nn
      call bsaclAcquireTotalNOfNodes__(int(nn, kind=c_int))
   end subroutine


   module subroutine bsacl_AcquireUsedModesList(modes)
      integer(IK), intent(in), target :: modes(:)
      integer(c_int) :: nmodes_eff_
      nmodes_eff_ = size(modes, kind = c_int)
      call bsaclAcquireUsedModesList__(c_loc(modes), nmodes_eff_)
   end subroutine


   module subroutine bsacl_AcquireWindCoeffs(wfc)
      real(RK), intent(in), target :: wfc(:, :, :)
      integer(c_int) :: nlibs_, ndegw_, nnodes_l_

      nlibs_    = size(wfc, 1, c_int)
      ndegw_    = size(wfc, 2, c_int)
      nnodes_l_ = size(wfc, 3, c_int)
      call bsaclAcquireWindCoeffs__(c_loc(wfc), nnodes_l_, nlibs_, ndegw_)
   end subroutine


   module subroutine bsacl_AcquireTurbComponentsList(tc)
      integer(IK), intent(in), target :: tc(:)
      integer(c_int) :: ntc_
      ntc_ = size(tc, kind=c_int)
      call bsaclAcquireTurbComponentsList__(c_loc(tc), ntc_)
   end subroutine


   module subroutine bsacl_AcquirePhiTimesCMat(phi_T_c)
      real(RK), intent(in), target :: phi_T_c(:, :, :)
      integer(c_int) :: nm_eff_, nnodes_l_, ndegw_

      nm_eff_   = size(phi_T_c, 1, c_int)
      nnodes_l_ = size(phi_T_c, 2, c_int)
      ndegw_    = size(phi_T_c, 3, c_int)
      call bsaclAcquirePhiTimesCMat__(c_loc(phi_T_c), nm_eff_, nnodes_l_, ndegw_)
! #ifdef BSA_DEBUG
!       block
!          integer(IK) :: m, n, d
!          integer(IK) :: id = 0
!          do d = 1, ndegw_
!             do n = 1, nnodes_l_
!                write(4532, '(2x, i5, *(2x, f15.5))') id, phi_T_c(:, n, d)
!                id = id + nm_eff_
!             enddo
!             write(4532, *) ''
!          enddo
!       end block
! #endif
   end subroutine


   module subroutine bsacl_AcquireNodalCorrelation(nod_corr)
      real(RK), intent(in), target :: nod_corr(:, :)
      integer(c_int) :: nnod_corr
      nnod_corr = size(nod_corr, dim=1)
      call bsaclAcquireNodalCorrelation__(c_loc(nod_corr), nnod_corr)
! #ifdef BSA_DEBUG
!       print '(/1x, 2a, i/)', INFO_MSG, 'n. nod_corr  = ', nnod_corr
! #endif
   end subroutine



   module subroutine bsacl_AcquireWindNodalVelocities(nod_vel)
      real(RK), intent(in), target :: nod_vel(:)

      call bsaclAcquireWindNodalVelocities__(c_loc(nod_vel))
   end subroutine


   module subroutine bsacl_AcquireWindNodalWindZones(nod_wz)
      integer(IK), intent(in), target :: nod_wz(:)

      call bsaclAcquireWindNodalWindZones__(c_loc(nod_wz))
   end subroutine


   module subroutine bsacl_AcquireWindTurbScales(wt_scl, nwz)
      real(RK), intent(in), target :: wt_scl(:, :, :)
      integer(IK), intent(in)      :: nwz

      call bsaclAcquireWindTurbScales__(c_loc(wt_scl), int(nwz, c_int))
   end subroutine


   module subroutine bsacl_AcquireWindTurbStd(wt_std, nwz)
      real(RK), intent(in), target :: wt_std(:, :)
      integer(IK), intent(in)      :: nwz

      call bsaclAcquireWindTurbStd__(c_loc(wt_std), int(nwz, c_int))
   end subroutine





#ifdef BSACL_ENABLE_EVALFCT_PTR
   ! This is what gets actually called in the C core
   module subroutine evalFunc_C_to_F_wrapper_(itc, nf, f, innl, res) bind(c)
      integer(c_int), value         :: itc
      integer(c_int), value         :: nf, innl
      real(c_double), intent(in), target :: f
      real(c_double), intent(inout) :: res(*)

      real(RK), allocatable   :: res_(:, :)

      ! locals for Fortran procedure
      integer(IK) :: nf_, itc_, i_
      real(RK), pointer :: f_(:)

      nf_  = int(nf,  kind=IK)
      itc_ = int(itc, kind=IK)
      call c_f_pointer(c_loc(f), f_, [nf_])

      ! Then, calling internal Fortran function pointer to actual eval function
      res_  = evalfct_ptr_(f_, nf_, itc_)
      do i_ = 1, innl
         res((i_-1)*nf_ + 1 : (i_*nf_)) = real(res_(:, i_), kind=c_double)
      enddo
      if (allocated(res_)) deallocate(res_)
   end subroutine


   module subroutine bsacl_AcquireEvaluationFunc(evalfct)
      procedure(bsacl_EvalFunc_C_to_Fortran__), pointer, intent(in) :: evalfct

      evalfct_ptr_ => evalfct
      call bsaclAcquireEvaluationFunc__(c_funloc(evalFunc_C_to_F_wrapper_))
   end subroutine


   module subroutine bsacl_AcquireEvalFuncByStrings(strings)
      character(len = *), target, intent(in) :: strings

      call bsaclAcquireEvalFuncByStrings__(c_loc(strings))
   end subroutine


   module subroutine bsacl_AcquireEvalFuncByFile(filename)
      character(len = *), target, intent(in) :: filename

      call bsaclAcquireEvalFuncByFile__(c_loc(filename))
   end subroutine
#endif



   module subroutine bsacl_AcquireComputationFreqs(i_thread, nfi, fi, nfj, fj)
      integer(IK), intent(in)      :: i_thread, nfi, nfj
      real(RK), intent(in), target :: fi(..), fj(..)
      integer(c_int) :: nfi_, nfj_

      nfi_ = int(nfi, kind=c_int)
      nfj_ = int(nfj, kind=c_int)
      call bsaclAcquireComputationFreqs__(int(i_thread, c_int), nfi_, c_loc(fi), nfj_, c_loc(fj))
   end subroutine


   module subroutine bsacl_AcquireBaseWindTurbPSD(S_uvw)
      real(RK), intent(in), target :: S_uvw(:, :)
      call bsaclAcquireBaseWindTurbPSD__(c_loc(S_uvw))
   end subroutine



   module subroutine bsacl_SetDeviceType(itype)
      integer(IK), intent(in) :: itype
      call bsaclSetDeviceType__(int(itype, kind=c_int))
   end subroutine



   module logical function bsacl_VerifyMaxAllocCondition(idim) result(bool)
      integer(int64), intent(in) :: idim
      integer(c_int64_t) :: idim_
      integer(c_int), target :: i_can_ = -1_c_int

      idim_ = int(idim, kind=c_int64_t)
      call bsaclVerifyMaxAllocCondition__(idim_, c_loc(i_can_))
      if (i_can_ == 0_c_int) then
         bool = .false.
      else
         bool = .true.
      endif
   end function



   module subroutine bsacl_Abort(ierr)
      integer(IK), intent(in) :: ierr
      call bsaclAbort__(int(ierr, kind=c_int))
   end subroutine


   module subroutine bsacl_Finalise()
      call bsaclFinalise__()
   end subroutine

end submodule
