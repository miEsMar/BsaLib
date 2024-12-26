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
module BsaLib_Structure

   use BsaLib_CONSTANTS, only: bsa_int_t, bsa_real_t, int32
   implicit none (type, external)
   private


   type, public :: StructureModalData_t

      !> n. of (in memory) kept structure vibration modes
      integer(bsa_int_t) :: nm_ = 0

      !> n. of "usable" structure vibration modes
      !> To account for non 1-normalised modes, if any.
      integer(bsa_int_t) :: nm_eff_ = 0

      !> List of "usable" structural vibration modes.
      !> They might be less than the actual computed, 
      !> since some of them might refer to torsional
      !> modes, etc..  (NM_EFF)
      integer(bsa_int_t), allocatable :: modes_(:)

      !> structural natural frequencies (NM)
      real(bsa_real_t), dimension(:), pointer    :: nat_freqs_;

      !> generalised modal matrix (NDOFs, NM)
      real(bsa_real_t), dimension(:, :), pointer :: phi_;

      !> modal damping ratios (NM)
      real(bsa_real_t), dimension(:), pointer :: xsi_;

      !> generalised mass matrix (NM)
      real(bsa_real_t), dimension(:), pointer :: Mm_;

      !> generalised damping matrix (Rayleigh) (NM)
      real(bsa_real_t), dimension(:, :), pointer :: Cm_;

      !> generalised stiffness matrix (NM, NM)
      real(bsa_real_t), dimension(:), pointer :: Km_;
   end type StructureModalData_t



   type, public :: StructureData_t

      !> n. of all nodes
      integer(bsa_int_t) :: nn_ = 0

      !> n. of all libs (per node)
      integer(bsa_int_t) :: nlibs_ = 0

      !> n. of total DOFs
      integer(bsa_int_t) :: ndofs_ = 0

      !> n. of actually loaded nodes
      integer(bsa_int_t) :: nn_load_ = 0

      !> n. of actually loaded DOFs (per loaded node)
      integer(bsa_int_t) :: nlibs_load_ = 0

      !> list of actual loaded nodes
      integer(bsa_int_t), pointer :: n_load_(:) => null();

      !> list of actual loaded DOFs (per node)
      integer(bsa_int_t), pointer :: libs_load_(:) => null();

      ! NOTE: pointer since we do not want to copy.
      !> Nodal coordinates.
      real(bsa_real_t), pointer :: coords_(:, :) => null()

      !> modal structure info
      type(StructureModalData_t) :: modal_


      ! NOTE: following allocatables since are local instances !

      !> structure time scales
      real(bsa_real_t), dimension(:), allocatable :: str_time_scales_

      !> background peak widths
      real(bsa_real_t), dimension(:, :), allocatable :: bkg_peak_width_

      !> resonant peak widths
      real(bsa_real_t), dimension(:), allocatable :: res_peak_width_

   contains

      procedure, public, pass :: SetNodalCoords
      procedure, public, pass :: SetNOfNodalDOFs
      procedure, public, pass :: SetTotalNOfNodes
      procedure, public, pass :: SetLoadedNodalDOFs
      procedure, public, pass :: SetLoadedNodes
      procedure, public, pass :: SetModalInfo
      procedure, public, pass :: SetKeptModes
      procedure, public, pass :: SetKeptModesDefault
      procedure, public, pass :: SetModalMatrices
      procedure, public, pass :: SetTotDamping
      procedure, public, pass :: ComputeResPeakWidths
      procedure, public, pass :: computeBKGPeakWidths
      procedure, public, pass :: clean
   end type StructureData_t





   interface

      module subroutine SetNodalCoords(this, nn, coords)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in)        :: nn
         real(bsa_real_t), target, contiguous  :: coords(:, :)
      end subroutine

      module subroutine SetNOfNodalDOFs(this, nlibs)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in) :: nlibs
      end subroutine


      module subroutine SetTotalNOfNodes(this, nn)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in) :: nn
      end subroutine



      module subroutine SetLoadedNodalDOFs(this, nlib, lib)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in) :: nlib
         integer(bsa_int_t), target, intent(in) :: lib(:)
      end subroutine



      module subroutine SetLoadedNodes(this, nnl, nl)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in) :: nnl
         integer(bsa_int_t), target, intent(in) :: nl(:)
      end subroutine






      module subroutine SetModalInfo(this, ndofs, nm, Phi, natf)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in) :: ndofs, nm
         real(bsa_real_t), intent(in), target :: Phi(ndofs, nm), natf(nm)
      end subroutine


      module subroutine SetKeptModes(this, modes)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in) :: modes(:)
      end subroutine


      !> NOTE: assumes that modes are not allocated yet.
      module subroutine SetKeptModesDefault(this)
         class(StructureData_t), intent(inout) :: this
      end subroutine



      module subroutine SetModalMatrices(this, nm, Mg, Kg, Cg)
         class(StructureData_t), intent(inout) :: this
         integer(bsa_int_t), intent(in) :: nm
         real(bsa_real_t), intent(in), target, dimension(nm) :: Mg, Kg
         real(bsa_real_t), intent(in), target :: Cg(nm, nm)
      end subroutine



      module subroutine SetTotDamping(this, xsi)
         class(StructureData_t), intent(inout) :: this
         real(bsa_real_t), target, intent(in) :: xsi(this%modal_%nm_)
      end subroutine






      module subroutine ComputeResPeakWidths(this)
         class(StructureData_t), intent(inout) :: this
      end subroutine



      module subroutine computeBKGPeakWidths(this, wind_scales)
         class(StructureData_t), intent(inout) :: this
         real(bsa_real_t), intent(in) :: wind_scales(:, :)
      end subroutine




      module subroutine clean(this)
         class(StructureData_t) :: this
      end subroutine

   end interface


end module BsaLib_Structure
