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
module BsaLib_MTriangZone
   
   use BsaLib_CONSTANTS,    only: bsa_real_t, bsa_int_t
   use BsaLib_MPoint,       only: MPoint_t
   use BsaLib_M2DPolygZone, only: M2DPolygZone_t
   implicit none (type, external)
   private
   

   ! TODO: triang and rect actually share almost everything..
   type, public, extends(M2DPolygZone_t) :: MTriangZone_t

      !> A pt.
      !> It is the upper point when rot is 0.
      type(MPoint_t) :: Apt_

      !> B point.
      !> It is the point to the right when the rotation is 0.
      type(MPoint_t) :: Bpt_

      !> Corner point.
      !> NOTE: we need it in case it ends being
      !>       either init/end point for meshing
      type(MPoint_t) :: Cpt_


      !> PAB abgle
      real(bsa_real_t) :: PABang_

   contains

      procedure, pass :: baseI => baseI_triang
      procedure, pass :: baseJ => baseJ_triang

      procedure, pass :: deduceDeltas => deduceDeltas_triang

      !> Returns total N of zone's meshing points.
      !> NOTE: for the moment, ONLY triangles rectangles supported
      !> NOTE: override inherited procedure.
      procedure, pass :: zoneTotNPts => zoneTotNPts_triang

      !> Sets a values for PAB angle.
      procedure, pass, private :: setPABangle

      !> Useful if we want specific deltas in I and J directions
      procedure, pass, private :: adaptToDeltas

      !> Deduce zone rotation from internal data
      procedure, pass, private :: deduceRotation

      !> Gets unary deltas along I and J dirs, based on zone rotation
      procedure, pass :: getRotatedUnaryDF

      !> Completely defines a triang zone, from points C, A, B.
      !> (NOTE: Either from refinements, or deltas)
      generic :: define => defineFromPts_norm, defineFromPts_refmtsORdeltas
      procedure, pass, private :: defineFromPts_norm, defineFromPts_refmtsORdeltas

      !> Computes a finally defined triang zone.
      procedure, pass :: compute

      !> Dumps a triang zone
      procedure, pass :: dump => dumpTZ

      !> Undumps a triang zone
      procedure, pass :: undump => undumpTZ

      procedure, pass :: interpolate => interpolateTZ
   end type


   !> MTriangZone constructors interafce
   interface MTriangZone
      module procedure MTriangZone_constr_def
   end interface
   public :: MTriangZone




   interface

      !> Default MTriangZone basic constructor
      module function MTriangZone_constr_def(name) result(this)
         character(len = *), intent(in), optional :: name
         type(MTriangZone_t) :: this
      end function



      !> Gets rect base along I-dir
      elemental module function baseI_triang(this) result(res)
         class(MTriangZone_t), intent(in) :: this
         real(bsa_real_t) :: res
      end function


      !> Gets rect base along J-dir
      elemental module function baseJ_triang(this) result(res)
         class(MTriangZone_t), intent(in) :: this
         real(bsa_real_t) :: res
      end function



      module subroutine deduceDeltas_triang(this)
         class(MTriangZone_t), intent(inout) :: this
      end subroutine



      !> Returns total N of zone's meshing points.
      !> NOTE: for the moment, ONLY triangles rectangles supported
      pure module function zoneTotNPts_triang(this) result(npt)
         class(MTriangZone_t), intent(in) :: this
         integer(bsa_int_t) :: npt
      end function



      !> Computes PAB angle, using "cosine" rule.
      !> NOTE: this only works if Cpt is the point at which
      !>       the wider angle is located!
      module subroutine setPABangle(this)
         class(MTriangZone_t), intent(inout) :: this
      end subroutine


      !> Used if we want some specific deltas
      !> along I and J directions
      module subroutine adaptToDeltas(this, dfi, dfj)
         class(MTriangZone_t), intent(inout) :: this
         real(bsa_real_t), value :: dfi, dfj
      end subroutine



      !> BUG: this routine might hide a bug.
      !>
      !> Deduce triangle zone rotation, knowing the 
      !> coordinates of the 3 triang points.
      !>
      !> Still, triangle rotation is given by the rotation
      !> needed to be done to reach the CA side, starting from
      !> Y positive axis of the GRS.
      module subroutine deduceRotation(this)
         class(MTriangZone_t), intent(inout) :: this
      end subroutine




      !> Gets actua points' delta increments based on this zone's
      !> rotation w.r.t GRS, when moving along the mesh.
      !> If NOT inverted:
      !>   - df_var : total delta along CA side (J-dir in LRS)
      !>   - df_cst : total delta along CB side (I-dir in LRS)
      module subroutine getRotatedUnaryDF(this, df_I_var, df_J_var, df_I_cst, df_J_cst)
         class(MTriangZone_t), intent(inout) :: this
         real(bsa_real_t), intent(out) :: df_I_var, df_J_var
         ! logical, intent(in), optional :: invert
         real(bsa_real_t), intent(out), optional :: df_I_cst, df_J_cst
      end subroutine




      !> Defines a triang zone from 3 triang points. 
      !> No need to specify which one is A or B point, this is 
      !> automatically deduced.
      !> Automatically choses refinements definition.
      module subroutine defineFromPts_norm(this, Cp, P1, P2)
         class(MTriangZone_t), intent(inout) :: this
         class(MPoint_t), intent(in)         :: Cp, P1, P2
      end subroutine


      !> Defines a triang zone from 3 triang points. 
      !> No need to specify which one is A or B point, this is 
      !> automatically deduced.
      !> Also, possible to specify either a value.
      !> If "is_refinement", then is either refinement along I or J direction.
      !> NOTE: this is the defult if no optional is passed.
      !>
      !> If "is_refinement" is false, but still a value is passed, 
      !> it is interpreted as a delta value, still along either I or J direction.
      module subroutine defineFromPts_refmtsORdeltas(this, Cp, P1, P2, is_refinement, val_types, val1, val2)
         class(MTriangZone_t), intent(inout) :: this
         class(MPoint_t), intent(in) :: Cp, P1, P2
         logical, intent(in)         :: is_refinement
         character(len = *), intent(in), optional :: val_types
         real(bsa_real_t), intent(in), optional   :: val1, val2
      end subroutine



      !> Computes a completely defined triang zone
      module subroutine compute(this)
         class(MTriangZone_t), intent(inout) :: this
      end subroutine



      !> Dumps a triang zone data for later reconstruction
      module subroutine dumpTZ(this)
         class(MTriangZone_t), intent(in) :: this
      end subroutine



      !> Undumps a triang zone data for later reconstruction
      module subroutine undumpTZ(this)
         class(MTriangZone_t), intent(inout) :: this
      end subroutine



      !> Implementation of triang zone interpolation methods
      module subroutine interpolateTZ( this &
#ifndef BSA_USE_POD_DATA_CACHING
         & , bfm   &
#endif
         & , pdata )
         class(MTriangZone_t), intent(inout) :: this
#ifndef BSA_USE_POD_DATA_CACHING
         real(bsa_real_t), intent(in)  :: bfm(:, :)
#endif
         class(*), pointer, intent(in) :: pdata
      end subroutine


   end interface

end module BsaLib_MTriangZone
