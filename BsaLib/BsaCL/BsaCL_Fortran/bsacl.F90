module BsaCL

#include "_devtypes.h"
#include "_errtypes.h"

   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   implicit none
   public
#ifdef BSACL_ENABLE_EVALFCT_PTR
   private :: evalFunc_C_to_F_wrapper_
#endif

#ifdef INT32
   integer, parameter, public :: IK = int32
#elif INT64
   integer, parameter, public :: IK = int64
#else
   integer, parameter, public :: IK = int32
#endif


#ifdef BSA_SINGLE_FLOATING_PRECISION
   integer, parameter, public :: RK = real32
#else
   integer, parameter, public :: RK = real64
#endif


#ifdef BSACL_ENABLE_EVALFCT_PTR
   abstract interface
      function bsacl_EvalFunc_C_to_Fortran__(f_, nf_, itc_) result(res)
         import IK, RK
         integer(IK), intent(in) :: nf_, itc_
         real(RK), intent(in)    :: f_(nf_)
         real(RK), allocatable, target :: res(:, :)
      end function
   end interface
#endif


   integer(kind = IK), parameter :: BSACL_SUCCESS                      = 0_IK
   integer(kind = IK), parameter :: BSACL_DEVICE_TYPE_CPU              = BSACL_DEVICE_TYPE_CPU_
   integer(kind = IK), parameter :: BSACL_DEVICE_TYPE_GPU              = BSACL_DEVICE_TYPE_GPU_
   integer(kind = IK), parameter :: BSACL_DEVICE_TYPE_ACC              = BSACL_DEVICE_TYPE_ACC_
   integer(kind = IK), parameter :: BSACL_DEVICE_TYPE_DEF              = BSACL_DEVICE_TYPE_DEF_
   integer(kind = IK), parameter :: BSACL_VALUE_MISMATCH_ERROR         = BSACL_VALUE_MISMATCH_ERROR_
   integer(kind = IK), parameter :: BSACL_PLATFORM_FIND_ERROR          = BSACL_PLATFORM_FIND_ERROR_
   integer(kind = IK), parameter :: BSACL_DEF_PLATFORM_ACQUIRE_ERROR   = BSACL_DEF_PLATFORM_ACQUIRE_ERROR_
   integer(kind = IK), parameter :: BSACL_CQUEUES_CREATION_ERROR       = BSACL_CQUEUES_CREATION_ERROR_
   integer(kind = IK), parameter :: BSACL_INVALID_DEVICE_TYPE          = BSACL_INVALID_DEVICE_TYPE_


   interface
      module integer(IK) function bsacl_Init(n_threads)
         integer(IK), intent(in), value :: n_threads
      end function


      module integer(IK) function bsacl_InitDeviceMemory()
      end function


      module integer(IK) function bsacl_Run(i_thread, res)
         integer(IK), intent(in), value :: i_thread
         real(RK), pointer, intent(in)  :: res(:)
      end function


      module function bsacl_SetKernelID(kid) result(ierr)
         integer(IK), value, intent(in) :: kid
         integer(IK) :: ierr
      end function


      module subroutine bsacl_AcquirePSDId(psdid)
         integer(IK), intent(in) :: psdid
      end subroutine


      module subroutine bsacl_AcquireStructModMat(modmat, natf)
         real(RK), intent(in), target :: modmat(:, :), natf(:)
      end subroutine


      module subroutine bsacl_AcquireModalMatrices(Mg, Cg, Kg)
         real(RK), intent(in), dimension(:),    target :: Mg, Kg
         real(RK), intent(in), dimension(:, :), target :: Cg
      end subroutine


      module subroutine bsacl_AcquireLoadedNodesList(nodes_load)
         integer(IK), intent(in), target :: nodes_load(:)
      end subroutine


      module subroutine bsacl_AcquireTotalNOfNodes(nn)
         integer(IK), intent(in) :: nn
      end subroutine


      module subroutine bsacl_AcquireUsedModesList(modes)
         integer(IK), intent(in), target :: modes(:)
      end subroutine


      module subroutine bsacl_AcquireWindCoeffs(wfc)
         real(RK), intent(in), target :: wfc(:, :, :)
      end subroutine


      module subroutine bsacl_AcquireTurbComponentsList(tc)
         integer(IK), intent(in), target :: tc(:)
      end subroutine


      module subroutine bsacl_AcquirePhiTimesCMat(phi_T_c)
         real(RK), intent(in), target :: phi_T_c(:, :, :)
      end subroutine


      module subroutine bsacl_AcquireNodalCorrelation(nod_corr)
         real(RK), intent(in), target :: nod_corr(:, :)
      end subroutine



      module subroutine bsacl_AcquireWindNodalVelocities(nod_vel)
         real(RK), intent(in), target :: nod_vel(:)
      end subroutine


      module subroutine bsacl_AcquireWindNodalWindZones(nod_wz)
         integer(IK), intent(in), target :: nod_wz(:)
      end subroutine


      module subroutine bsacl_AcquireWindTurbScales(wt_scl, nwz)
         real(RK), intent(in), target :: wt_scl(:, :, :)
         integer(IK), intent(in)      :: nwz
      end subroutine


      module subroutine bsacl_AcquireWindTurbStd(wt_std, nwz)
         real(RK), intent(in), target :: wt_std(:, :)
         integer(IK), intent(in)      :: nwz
      end subroutine



#ifdef BSACL_ENABLE_EVALFCT_PTR
      module subroutine evalFunc_C_to_F_wrapper_(itc, nf, f, innl, res) bind(c)
         import c_int, c_double
         integer(c_int), value         :: itc
         integer(c_int), value         :: nf, innl
         real(c_double), intent(in), target :: f
         real(c_double), intent(inout) :: res(*)
      end subroutine


      module subroutine bsacl_AcquireEvaluationFunc(evalfct)
         procedure(bsacl_EvalFunc_C_to_Fortran__), pointer, intent(in) :: evalfct
      end subroutine


      module subroutine bsacl_AcquireEvalFuncByStrings(strings)
         character(len = *), target, intent(in) :: strings
      end subroutine


      module subroutine bsacl_AcquireEvalFuncByFile(filename)
         character(len = *), target, intent(in) :: filename
      end subroutine
#endif


      module subroutine bsacl_AcquireComputationFreqs(i_thread, nfi, fi, nfj, fj)
         integer(IK), intent(in)      :: i_thread, nfi, nfj
         real(RK), intent(in), target :: fi(..), fj(..)
      end subroutine


      module subroutine bsacl_AcquireBaseWindTurbPSD(S_uvw)
         real(RK), intent(in), target :: S_uvw(:, :)
      end subroutine


      module subroutine bsacl_SetDeviceType(itype)
         integer(IK), intent(in) :: itype
      end subroutine


      module logical function bsacl_VerifyMaxAllocCondition(idim) result(bool)
         integer(int64), intent(in) :: idim
      end function


      module subroutine bsacl_Abort(ierr)
         integer(IK), intent(in) :: ierr
      end subroutine


      module subroutine bsacl_Finalise()
      end subroutine
   end interface

end module
