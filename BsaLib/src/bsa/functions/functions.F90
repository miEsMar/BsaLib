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
module BsaLib_Functions

   use BsaLib_CONSTANTS
   use BsaLib_Data, only: wd, struct_data, settings, dimM_bisp_, dimM_psd_
   implicit none (type, external)
   public
   private :: wd, struct_data, settings, dimM_bisp_, dimM_psd_

   ! make a local internal copy
   integer(bsa_int_t) :: NFREQS, NNODES, NNODESL, NLIBS, NLIBSL
   integer(bsa_int_t) :: NMODES, NMODES_EFF
   integer(bsa_int_t) :: NPSDEL, NTCOMPS, NDIRS = 1
   integer(bsa_int_t), allocatable :: MODES(:)
   integer(bsa_int_t), allocatable :: TCOMPS(:), DIRS(:)

   integer              :: MSHR_SVD_LWORK = - 1
   integer, allocatable :: MSHR_SVD_INFO
   real(bsa_real_t), allocatable :: MSHR_SVD_WORK(:)


   interface

      module subroutine setBsaFunctionLocalVars()
      end subroutine


      module subroutine prefetchSVDWorkDim_()
      end subroutine

      module subroutine cleanSVDWorkInfo_()
      end subroutine



      module subroutine getFM_full_tnm_scalar_msh_(bfm, fi, fj)
         real(bsa_real_t), intent(inout), contiguous :: bfm(:, :)
         real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)
      end subroutine


      module subroutine getFM_full_tm_scalar_msh_POD_(bfm, fi, fj)
         real(bsa_real_t), intent(inout), contiguous :: bfm(:, :)
         real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)
      end subroutine


      module subroutine getRM_full_scalar_msh_(brm, fi, fj, bfm)
         real(bsa_real_t), intent(inout), contiguous :: brm(:, :)
         real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)
         real(bsa_real_t), intent(in), contiguous :: bfm(:, :)
      end subroutine



      module subroutine getFM_diag_tnm_scalar_msh_(bfm, fi, fj)
         real(bsa_real_t), intent(inout), contiguous :: bfm(:, :)
         real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)
      end subroutine


      module subroutine getRM_diag_scalar_msh_(brm, fi, fj, bfm)
         real(bsa_real_t), intent(inout), contiguous :: brm(:, :)
         real(bsa_real_t), intent(in), contiguous :: fi(:), fj(:)
         real(bsa_real_t), intent(in), contiguous :: bfm(:, :)
      end subroutine









      module subroutine getFM_full_tnm_vect_cls_(f, Suvw, psd, bisp)
         real(bsa_real_t), intent(in) :: f(NFREQS)
         real(bsa_real_t), intent(in) :: Suvw(NFREQS, NPSDEL)
         real(bsa_real_t), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine



      module subroutine getRM_full_vect_cls_(f, psd, bisp)
         real(bsa_real_t), intent(in)                 :: f(NFREQS)
         real(bsa_real_t), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine



      module subroutine getFM_diag_tnlm_vect_cls_(f, Suvw, psd, bisp)
         real(bsa_real_t), intent(in) :: f(NFREQS)
         real(bsa_real_t), intent(in) :: Suvw(NFREQS, NPSDEL)
         real(bsa_real_t), intent(inout), allocatable :: psd(:, :), bisp(:, :, :)
      end subroutine   



      module subroutine getRM_diag_vect_cls_(f, psd, bisp)
         real(bsa_real_t), intent(in)                 :: f(NFREQS)
         real(bsa_real_t), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine







      pure module subroutine getFM_full_tnm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         integer(bsa_int_t), intent(in)  :: ii, ij
         real(bsa_real_t), intent(in)    :: fi, fj
         real(bsa_real_t), intent(in)    :: Suvw(NFREQS, NPSDEL)
         real(bsa_real_t), intent(in)    :: Suvw_pad(NPSDEL)
         real(bsa_real_t), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine



      module subroutine getRM_full_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         integer(bsa_int_t), intent(in) :: ii, ij
         real(bsa_real_t), intent(in)   :: fi, fj
         real(bsa_real_t), intent(in)   :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real(bsa_real_t), intent(out)  :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine


      !> BUG: this routine is adapted to the case where we use
      !>      convention on PULSATION.
      !>      Please, adapt it to the case of convention over FREQUENCIES.
      pure module subroutine getFM_diag_tnlm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         integer(bsa_int_t), intent(in)  :: ii, ij
         real(bsa_real_t), intent(in)    :: fi, fj
         real(bsa_real_t), intent(in)    :: Suvw(NFREQS, NPSDEL)
         real(bsa_real_t), intent(in)    :: Suvw_pad(NPSDEL)
         real(bsa_real_t), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine



      module subroutine getRM_diag_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         integer(bsa_int_t), intent(in) :: ii, ij  ! freqs indexes
         real(bsa_real_t), intent(in)   :: fi, fj
         real(bsa_real_t), intent(in)   :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real(bsa_real_t), intent(out)  :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine




      pure module subroutine getBR_SFm_val_(nm, Suvw, fnat, im, m, psd)
         integer(bsa_int_t), intent(in)  :: im, m, nm
         real(bsa_real_t), intent(in)    :: Suvw(nm, NPSDEL), fnat
         real(bsa_real_t), intent(inout) :: psd
      end subroutine

   end interface

end module BsaLib_Functions
