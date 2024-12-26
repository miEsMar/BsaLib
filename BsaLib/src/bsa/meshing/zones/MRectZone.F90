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
module BsaLib_MRectZone

   use BsaLib_CONSTANTS,    only: bsa_real_t
   use BsaLib_MPoint,       only: MPoint_t
   use BsaLib_M2DPolygZone, only: M2DPolygZone_t
   implicit none (type, external)
   private

   type, public, extends(M2DPolygZone_t) :: MRectZone_t

      !> Init point (from where everything else is made reference)
      type(MPoint_t) :: Ipt_

      !> End point
      type(MPoint_t) :: Ept_

      !> Rect base along I-dir (x axis)
      real(bsa_real_t) :: base_I_ = 0._bsa_real_t

      !> Rect base along J-dir (y axis)
      real(bsa_real_t) :: base_J_ = 0._bsa_real_t

      !> Delta freq [Hz] along I-dir (x axis)
      real(bsa_real_t) :: deltaf_I_ = 0._bsa_real_t

      !> Delta freq [Hz] along J-dir (y axis)
      real(bsa_real_t) :: deltaf_J_ = 0._bsa_real_t


      logical, private :: refmts_set_ = .false.
      logical, private :: deltas_set_ = .false.

   contains
      private
      procedure, pass, public :: baseI => baseI_rct
      procedure, pass, public :: baseJ => baseJ_rct
      procedure, pass, public :: deduceDeltas => deduceDeltas_rect

      procedure, pass, public :: getAPoint
      procedure, pass, public :: getBPoint
      procedure, pass, public :: setDeltas
      procedure, pass, public :: deduceRefinements
      generic, public :: defineFromDeltas => defineFromDeltas_refinements, defineFromDeltas_maxvalues
      procedure, pass :: defineFromDeltas_refinements, defineFromDeltas_maxvalues
      generic, public :: defineFromEndPtCoordAndBase => defineFromEndPtCoordAndBase_norm, defineFromEndPtCoordAndBase_forceDeltas
      procedure, pass :: defineFromEndPtCoordAndBase_norm, defineFromEndPtCoordAndBase_forceDeltas
      procedure, pass, public  :: define
      procedure, pass :: setIEpts
      procedure, pass :: getOtherBase
      ! procedure, pass :: getIJfsteps
      ! procedure, pass :: validateDeltas
      ! procedure, pass :: reconstructZoneBaseMesh
      procedure, pass, public  :: compute => compute_rz_
      procedure, pass, private :: compute_rz_
      procedure, pass, public  :: getNthQuadVtx
      procedure, pass, public  :: dump => dumpRZ
      procedure, pass, public  :: undump => undumpRZ
      procedure, pass, public  :: interpolate => interpolateRZ
   end type MRectZone_t



   interface MRectZone
      module procedure MRectZone_t_custom_constructor
   end interface
   public :: MRectZone



   ! main interface to module procedures
   interface
      
      module function MRectZone_t_custom_constructor(rot, name) result(this)
         real(bsa_real_t), intent(in), optional :: rot
         character(len=*), intent(in), optional :: name
         ! NOTE: compiler uses default initialisation here (built-in)
         type(MRectZone_t) :: this
      end function


      !> Gets rect base along I-dir
      elemental module function baseI_rct(this) result(res)
         class(MRectZone_t), intent(in) :: this
         real(bsa_real_t) :: res
      end function


      !> Gets rect base along J-dir
      elemental module function baseJ_rct(this) result(res)
         class(MRectZone_t), intent(in) :: this
         real(bsa_real_t) :: res
      end function



      !> Get A point.
      !> A point is the point defined, starting from I point, 
      !> along the J-dir (Y-axis) parallel side.
      pure module function getAPoint(this) result(pt)
         class(MRectZone_t), intent(in) :: this
         type(MPoint_t) :: pt
      end function


      !> Get B point.
      !> B point is the point defined, starting from I point, 
      !> along the I-dir (X-axis) parallel side.
      pure module function getBPoint(this) result(pt)
         class(MRectZone_t), intent(in) :: this
         type(MPoint_t) :: pt
      end function


      !> If refinements are set, deltas are deduced
      module subroutine deduceDeltas_rect(this)
         class(MRectZone_t), intent(inout) :: this
      end subroutine


      !> Set frequency deltas
      module subroutine setDeltas(this, dfi, dfj, adapt)
         class(MRectZone_t), intent(inout) :: this
         real(bsa_real_t), intent(in)             :: dfi, dfj
         logical, intent(in), optional     :: adapt
      end subroutine


      !> If deltas are set, deduces refinements
      module subroutine deduceRefinements(this, adapt)
         class(MRectZone_t), intent(inout) :: this
         logical, intent(in) :: adapt
      end subroutine



      !> Defines zone bases (sides) from given Deltas, specifying refinement.
      !> NOTE: Init and End zone points MUST be known.
      module subroutine defineFromDeltas_refinements(this, pt, loc, dfi, dfj, ni, nj)
         class(MRectZone_t), intent(inout) :: this

         !> Point specification at location 'loc'
         class(MPoint_t), intent(in) :: pt

         !> Location of the specified point.
         !> Can be:
         !>   'i' -> init point
         !>   'c' -> centre point
         !>   'e' -> end point
         character(len=1) :: loc

         !> Delta values
         real(bsa_real_t), intent(in) :: dfi, dfj

         !> Refinements
         integer, value :: ni, nj
      end subroutine





      !> Defines zone bases (sides) from given Deltas, specifying max deltas values.
      !> NOTE: Init and End zone points MUST be known.
      !> If force is present and true, forces deltas to reach max values specified
      !> by 'valI' and 'valJ'. In this case, deltas are then slightly adjusted to fit.
      module subroutine defineFromDeltas_maxvalues(this, pt, loc, dfi, dfj, maxF_i, maxF_j, force, exceed)
         class(MRectZone_t), intent(inout) :: this

         !> Point specification, which location is specified by 'loc'
         class(MPoint_t), intent(in)  :: pt

         !> Point location:
         !>   'i' -> Init location
         !>   'c' -> centre (NOT YET IMPLEMENTED)
         !>   'e' -> end    (NOT YET IMPLEMENTED)
         character(len=1), intent(in) :: loc

         !> Delta values
         real(bsa_real_t), value :: dfi, dfj

         !> Max deltas values
         real(bsa_real_t) :: maxF_i, maxF_j

         !> adjusts deltas to max values specified.
         logical, intent(in), optional :: force


         logical, intent(in), optional :: exceed
      end subroutine





      module subroutine define(this, pt, loc, base_i, base_j)
         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in) :: pt
         character(len=1), optional, intent(in) :: loc
         real(bsa_real_t), intent(in), optional  :: base_i, base_j
      end subroutine






      module subroutine defineFromEndPtCoordAndBase_norm(&
         this, Pi, coord_val, coord_ty_ch, baseval, base_dir, called)

         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in) :: Pi
         real(bsa_real_t), intent(in) :: coord_val
         character(len = 1), intent(in) :: coord_ty_ch
         real(bsa_real_t), intent(in) :: baseval
         character(len = 1), intent(in) :: base_dir
         logical, intent(in)            :: called
      end subroutine



      module subroutine defineFromEndPtCoordAndBase_forceDeltas(&
         this, Pi, coord_val, coord_ty_ch, baseval, base_dir, dfi, dfj)

         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in) :: Pi
         real(bsa_real_t), intent(in) :: coord_val
         character(len = 1), intent(in) :: coord_ty_ch
         real(bsa_real_t), intent(in) :: baseval
         character(len = 1), intent(in) :: base_dir
         real(bsa_real_t), intent(in)          :: dfi, dfj
      end subroutine







      !> Sets, based on passed point location
      !> Init and End zone points.
      module subroutine setIEpts(this, pt, loc)
         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in)  :: pt
         character(len=1), intent(in) :: loc
      end subroutine



      ! !> Avoid setting a delta smaller than given limit
      ! elemental module impure subroutine validateDeltas(this, lval)
      !    class(MRectZone_t), intent(inout) :: this
      !    real(bsa_real_t), intent(in) :: lval
      ! end subroutine


      !> Automatically computes the second remaining (unknown) rect base
      !> based on the point's coordinates and the base that we already defined.
      !> BUG: cannot see when a base is negative (i.e. when END pt is "behind" INIT one).
      module subroutine getOtherBase(this, pt, base_dir, known_coord, coord_val)
         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in)  :: pt
         character(len=1), intent(in) :: base_dir, known_coord
         real(bsa_real_t), intent(in) :: coord_val
      end subroutine



      ! !> Gets actualised frequency deltas along two main
      ! !> sides directions (I, J), actualised based on *this
      ! !> zone rotation w.r.t. GRS.
      ! module subroutine getIJfsteps(this, dfIi, dfIj, dfJi, dfJj)
      !    class(MRectZone_t), intent(in) :: this
      !    real(bsa_real_t), intent(out)  :: dfIi, dfIj, dfJi, dfJj
      ! end subroutine



      ! !> Returns the whole zone reconstructed mesh
      ! !> Maybe used for visualising.
      ! module function reconstructZoneBaseMesh(this) result(msh)
      !    class(MRectZone_t), intent(in) :: this
      !    !> BUG: might be 2-rank array instead of 3!
      !    real(bsa_real_t) :: msh(2, this%nj_, this%ni_)
      ! end function


      !> Actual zone comutation (pre phase).
      module subroutine compute_rz_(this)
         class(MRectZone_t), intent(inout) :: this
      end subroutine



      !> Gets vertex point pt coordinates in Nth quadrant
      !> w.r.t. Center point, in zone's LRS.
      pure module function getNthQuadVtx(this, iquad) result(pt)
         class(MRectZone_t), intent(in) :: this
         integer, intent(in) :: iquad
         type(MPoint_t) :: pt
      end function





      !> Dumps a RECTANGULAR zone.
      !>
      !> NOTE: Each specific zone dumping method is called 
      !>       from the STATIC MZone_t procedure DumpZone().
      module subroutine dumpRZ(this)
         class(MRectZone_t), intent(in) :: this
      end subroutine



      !> Undumps a RECTANGULAR zone.
      !>
      !> NOTE: Each specific zone dumping method is called 
      !>       from the STATIC MZone_t procedure UndumpZone().
      !> NOTE: This is the equivalent as reconstruct() of MATLAB.
      module subroutine undumpRZ(this)
         class(MRectZone_t), intent(inout) :: this
      end subroutine



      !> Implementation of rect zone interpolation methods
      module subroutine interpolateRZ( this &
#ifndef BSA_USE_POD_DATA_CACHING
         & , bfm   &
#endif
         & , pdata )
         class(MRectZone_t), intent(inout) :: this
#ifndef BSA_USE_POD_DATA_CACHING
         real(bsa_real_t), intent(in) :: bfm(:, :)
#endif
         class(*), pointer, intent(in) :: pdata
      end subroutine


   end interface


end module BsaLib_MRectZone
