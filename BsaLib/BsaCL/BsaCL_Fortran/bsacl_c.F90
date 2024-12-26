module BsaCL_c
   use, intrinsic :: iso_c_binding
   implicit none
   public
   interface
      integer function bsaclInit__(n_threads) bind(c, name="bsaclInit")
         import c_int
         integer(c_int), value :: n_threads
      end function


      integer function bsaclInitDeviceMemory__() bind(c, name="bsaclInitDeviceMemory")
      end function


      integer function bsaclRun__(i_thread, dim, res) bind(c, name="bsaclRun")
         import c_int, c_ptr, c_size_t
         integer(c_int),    value :: i_thread
         integer(c_size_t), value :: dim
         type(c_ptr),       value :: res
      end function


      function bsaclSetKernelID__(kid) result(ierr) bind(c, name="bsaclSetKernelID")
         import c_int
         integer(c_int), value :: kid
         integer(c_int) :: ierr
      end function


      subroutine bsaclAcquirePSDId__(psdid) bind(c, name="bsaclAcquirePSDId")
         import c_int
         integer(c_int), value :: psdid
      end subroutine


      subroutine bsaclAcquireStructModMat__(modmat, natf, ndofs, nmodes) bind(c, name="bsaclAcquireStructModMat")
         import c_ptr, c_int
         type(c_ptr), value    :: modmat, natf
         integer(c_int), value :: ndofs, nmodes
      end subroutine


      subroutine bsaclAcquireModalMatrices__(Mg, Cg, Kg) bind(c, name="bsaclAcquireModalMatrices")
         import c_ptr
         type(c_ptr), value :: Mg, Cg, Kg
      end subroutine


      subroutine bsaclAcquireLoadedNodesList__(nodes_load, nnodes_l) &
            bind(c, name="bsaclAcquireLoadedNodesList")
         import c_ptr, c_int
         type(c_ptr), value    :: nodes_load
         integer(c_int), value :: nnodes_l
      end subroutine


      subroutine bsaclAcquireTotalNOfNodes__(nn) bind(c, name="bsaclAcquireTotalNOfNodes")
         import c_int
         integer(c_int), value :: nn
      end subroutine


      subroutine bsaclAcquireUsedModesList__(modes, nmodes_eff) bind(c, name="bsaclAcquireUsedModesList")
         import c_ptr, c_int
         type(c_ptr), value    :: modes
         integer(c_int), value :: nmodes_eff
      end subroutine


      subroutine bsaclAcquireWindCoeffs__(wfc, nnodes_l, nlibs, ndegw) bind(c, name="bsaclAcquireWindCoeffs")
         import c_ptr, c_int
         type(c_ptr), value    :: wfc
         integer(c_int), value :: nnodes_l, nlibs, ndegw
      end subroutine


      subroutine bsaclAcquireTurbComponentsList__(tc, ntc) bind(c, name="bsaclAcquireTurbComponentsList")
         import c_ptr, c_int
         type(c_ptr), value    :: tc
         integer(c_int), value :: ntc
      end subroutine


      subroutine bsaclAcquirePhiTimesCMat__(phi_T_c, nmodes, nnodes_l, ndegw) bind(c, name="bsaclAcquirePhiTimesCMat")
         import c_ptr, c_int
         type(c_ptr), value    :: phi_T_c
         integer(c_int), value :: nnodes_l, nmodes, ndegw
      end subroutine


      subroutine bsaclAcquireNodalCorrelation__(nod_corr, nnod_corr) bind(c, name="bsaclAcquireNodalCorrelation")
         import c_ptr, c_int
         type(c_ptr), value    :: nod_corr
         integer(c_int), value :: nnod_corr
      end subroutine



      subroutine bsaclAcquireWindNodalVelocities__(nod_vel) bind(c, name="bsaclAcquireWindNodalVelocities")
         import c_ptr
         type(c_ptr), value :: nod_vel
      end subroutine


      subroutine bsaclAcquireWindNodalWindZones__(nod_wz) bind(c, name="bsaclAcquireWindNodalWindZones")
         import c_ptr
         type(c_ptr), value :: nod_wz
      end subroutine


      subroutine bsaclAcquireWindTurbScales__(wt_scl, nwz) bind(c, name="bsaclAcquireWindTurbScales")
         import c_ptr, c_int
         type(c_ptr), value    :: wt_scl
         integer(c_int), value :: nwz
      end subroutine


      subroutine bsaclAcquireWindTurbStd__(wt_std, nwz) bind(c, name="bsaclAcquireWindTurbStd")
         import c_ptr, c_int
         type(c_ptr), value    :: wt_std
         integer(c_int), value :: nwz
      end subroutine



#ifdef BSACL_ENABLE_EVALFCT_PTR
      subroutine bsaclAcquireEvaluationFunc__(fct) bind(c, name="bsaclAcquireEvalFunc")
         import c_funptr
         type(c_funptr), value :: fct
      end subroutine


      subroutine bsaclAcquireEvalFuncByStrings__(strings) bind(c, name="bsaclAcquireEvalFuncByStrings")
         import c_ptr
         type(c_ptr), value :: strings
      end subroutine


      subroutine bsaclAcquireEvalFuncByFile__(filename) bind(c, name="bsaclAcquireEvalFuncByFile")
         import c_ptr
         type(c_ptr), value :: filename
      end subroutine
#endif


      subroutine bsaclAcquireComputationFreqs__(i_thread, nfi, fi, nfj, fj) bind(c, name="bsaclAcquireComputationFreqs")
         import c_ptr, c_int
         integer(c_int), value :: i_thread, nfi, nfj
         type(c_ptr), value    :: fi, fj
      end subroutine


      subroutine bsaclAcquireBaseWindTurbPSD__(S_uvw) bind(c, name="bsaclAcquireBaseWindTurbPSD")
         import c_ptr
         type(c_ptr), value :: S_uvw
      end subroutine


      subroutine bsaclSetDeviceType__(itype) bind(c, name="bsaclSetDeviceType")
         import c_int
         integer(c_int), value :: itype
      end subroutine


      subroutine bsaclVerifyMaxAllocCondition__(idim, i_can_alloc) bind(c, name="bsaclVerifyMaxAllocCondition")
         import c_ptr, c_int64_t
         integer(c_int64_t), value :: idim
         type(c_ptr), value        :: i_can_alloc
      end subroutine


      subroutine bsaclAbort__(ierr) bind(c, name="bsaclAbort")
         import c_int
         integer(c_int), value :: ierr
      end subroutine


      subroutine bsaclFinalise__() bind(c, name="bsaclFinalise")
      end subroutine
   end interface

end module
