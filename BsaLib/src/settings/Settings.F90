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
module BsaLib_Settings

   use BsaLib_CONSTANTS, only: bsa_int_t, bsa_real_t, int32
   implicit none (type, external)
   private

   !> Minimum rounding precision to guarantee.
   integer(int32), public :: i_min_round_prec_ = 10_int32

   type, public :: settings_t

      ! === GENERAL settings

      !> Subanalysis type
      !> ==1, classic
      !> ==2, mesher
      !> ==3, BOTH (comparison) 
      integer(bsa_int_t) :: i_suban_type_ = 2

      !> NOTE: only for "classic" suban type.
      !> Manage which version to use (specially for dev testing).
      !> ==1 uses old adapted spectra
      !> ==2, uses new spectra, sistematic
      integer(bsa_int_t) :: i_vers_ = 2

      !> Turbulence PSDs scaling convention
      !> ==1, pulsations  (-infty, +infty)
      !> ==2, frequencies (0, +infty)
      integer(bsa_int_t) :: i_def_scaling_ = 1

      !> If ==0, do not compute PSDs
      integer(bsa_int_t) :: i_compute_psd_ = 1

      !> If ==0, do not compute BISPs
      integer(bsa_int_t) :: i_compute_bisp_ = 1

      !> Whether computing FULL matrices or not.
      !> True if ==1.
      integer(bsa_int_t) :: i_only_diag_ = 0

      !> Activate for testing some new features.
      integer(bsa_int_t) :: i_test_mode_ = 0



      ! === CLASSIC settings


      !> If suban=="classic", number of sistematic freqs.
      integer(bsa_int_t) :: nfreqs_ = 0

      !> If suban=="classic", constant delta frequency.
      real(bsa_real_t) :: df_ = 0._bsa_real_t

      !> If ==0, using VECTORISED functions version.
      !> Otherwise (==1), using SCALAR versions.
      integer(bsa_int_t) :: i_scalar_vers_ = 0

      !> Bisp symmetry case.
      !> 0 = full
      !> 2 = half
      !> 4 = fourth
      integer(bsa_int_t) :: i_bisp_sym_ = 0

      !> 3D bisp matrix symmetry exploitation.
      !> 0 = no
      !> 1 = yes
      !> NOTE: if i_bisp_sym_==4, automatically 0
      integer(bsa_int_t) :: i_spctr_sym_ = 0




      ! === MESHER settings


      !> If ==1, using SVD to S_uvw matrices
      integer(bsa_int_t) :: i_use_svd_ = 1

      !> How many points (per side) for meshing 
      !> main central BKG peak zone.
      integer(bsa_int_t) :: bkg_base_rfmnt_ = 20

      !> 
      real(bsa_real_t) :: max_area_ext_ = 2._bsa_real_t

      !> How much to extend BKG peak area influence.
      real(bsa_real_t) :: bkg_area_ext_ = 2._bsa_real_t

      !> How much to extend general peak area influence.
      real(bsa_real_t) :: peak_area_ext_ = 2._bsa_real_t

      !> If true, we get up to 2*max_freq.
      integer(bsa_int_t) :: i_full_coverage_ = 1

      !> Controls wheter to include modal info when
      !> writing to dump file.
      !> Unactive by default!
      !> Warn if gets activated.
      integer(bsa_int_t) :: i_dump_modal_ = 0


   contains

      procedure, public, pass :: SetSubanType
      procedure, public, pass :: SetVersion
      procedure, public, pass :: SetScalingType
      procedure, public, pass :: ActivateSpectraComputation
      procedure, public, pass :: SetExtension
      procedure, public, pass :: TestMode
      procedure, public, pass :: setClsSettings
      procedure, public, pass :: SetMshrSetts
   end type settings_t



   interface

      !> Sets sub analysis type
      module subroutine SetSubanType(this, isuban)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), intent(in)   :: isuban
      end subroutine

      !> Set version
      module subroutine SetVersion(this, ivers)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), intent(in)   :: ivers
      end subroutine

      !> Sets PSDs scaling convention.
      !>   ==1, default, pulsation (-infty, +infty)
      !>   ==2, frequencies [WARNING] (0, +infty)
      module subroutine SetScalingType(this, idefsc)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), intent(in)   :: idefsc
      end subroutine


      !> Allows to control spectra comptation.
      !> Pass 0 to deactivate.
      module subroutine ActivateSpectraComputation(this, ipsd, ibisp)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), value :: ipsd, ibisp
      end subroutine


      !> Allows to specify whether full 2d/3d matrices are to be computed.
      !> If 0, only main diagonal elements are computed (uncorrelated case).
      module subroutine SetExtension(this, ionlydiag)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), intent(in)   :: ionlydiag
      end subroutine


      !> Controls whether testing mode is active.
      module subroutine TestMode(this, itest)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), intent(in)   :: itest
      end subroutine



      !> Sets main Classic suban settings.
      module subroutine setClsSettings(this, nfreqs, df)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), intent(in)   :: nfreqs
         real(bsa_real_t), intent(in) :: df
      end subroutine


      !> Sets main Mesher suban settings.
      module subroutine SetMshrSetts(this, isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
         class(settings_t), intent(inout) :: this
         integer(bsa_int_t), value :: isvd, bkgrfmt, ifcov, idumpmod
         real(bsa_real_t),   value :: bkgaext, genpaext, maxaext
      end subroutine

   end interface


end module
