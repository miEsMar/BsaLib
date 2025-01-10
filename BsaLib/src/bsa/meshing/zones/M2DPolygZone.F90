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
module BsaLib_M2DPolygZone

   use BsaLib_CONSTANTS, only: int32, bsa_real_t, bsa_int_t &
                              , CST_PIGREC, CST_PIt2, CST_PId2, CST_PIt3d2
   use BsaLib_MPoint,    only: MPoint_t
   use BsaLib_MZone,     only: MZone_t
   implicit none (type, external)
   private

   type, public, abstract, extends(MZone_t) :: M2DPolygZone_t

      !> Refinement (n meshing pts) along I-dir
      integer(int32) :: ni_ = 0

      !> Refinement (n meshing pts) along J-dir
      integer(int32) :: nj_ = 0

      !> Rotation angle w.r.t to XY plane axes (CLOCKWISE)
      real(bsa_real_t) :: rot_ = 0._bsa_real_t

   contains

      ! By default rect case. Others -> specialise
      procedure, pass :: zoneTotNPts => zoneTotNPts_pol2D_

      procedure, pass :: refinements
      procedure, pass :: setRotation
      procedure, pass :: isGRSAligned
      procedure, pass :: setRefinements

      procedure, pass :: getUnaryDeltas

      procedure(VoidSub),     pass, deferred :: deduceDeltas
      procedure(RealVoidFct), pass, deferred :: baseI
      procedure(RealVoidFct), pass, deferred :: baseJ
   end type



   abstract interface
      elemental function RealVoidFct(this) result(res)
         import :: M2DPolygZone_t, bsa_real_t
         class(M2DPolygZone_t), intent(in) :: this
         real(bsa_real_t) :: res
      end function

      subroutine VoidSub(this)
         import :: M2DPolygZone_t
         class(M2DPolygZone_t), intent(inout) :: this
      end subroutine
   end interface




contains


   !> Gets total number of zone's meshing points
   pure function zoneTotNPts_pol2D_(this) result(np)
      class(M2DPolygZone_t), intent(in) :: this
      integer(int32) :: np

      np = this%ni_ * this%nj_
   end function



   !> Get this zone's refinements
   pure function refinements(this) result(refmts)
      class(M2DPolygZone_t), intent(in) :: this
      integer(int32) :: refmts(2)

      refmts(1) = this%ni_
      refmts(2) = this%nj_
   end function



   subroutine setRotation(this, rot, deg)
      class(M2DPolygZone_t), intent(inout) :: this
      real(bsa_real_t), value       :: rot
      logical, intent(in), optional :: deg
      logical :: is_deg = .false.
      integer(int32) :: n2pirot = 0

      if (present(deg)) then
         if (deg) is_deg = .true.
      endif

      if (is_deg) rot = rot / 180 * CST_PIGREC

      ! we want rot € [0, 2*pi)
      ! 2*pi -> 0. (this why ==) 
      if (rot >= CST_PIt2) then
         n2pirot = floor(rot / CST_PIt2)
         rot     = rot - n2pirot * CST_PIt2
      endif

      this%rot_ = rot
   end subroutine





   !> Verifies if zone is not diagonal w.r.t.
   !> to the orientation of the GRS (Global Reference System)
   !> BUG: might need usage of rounding precision.
   elemental function isGRSAligned(this) result(bool)
      class(M2DPolygZone_t), intent(in) :: this
      logical :: bool

      bool = .false.

      if (this%rot_ == 0._bsa_real_t .or. &
          this%rot_ == CST_PId2      .or. &
          this%rot_ == CST_PIGREC    .or. &
          this%rot_ == CST_PIt3d2) bool = .true.
   end function




   !> Sets zone refinements.
   subroutine setRefinements(this, ni, nj, force)
      class(M2DPolygZone_t), intent(inout) :: this
      integer(bsa_int_t), intent(in) :: ni, nj
      logical, intent(in), optional  :: force
      logical :: do_force = .false.

      if (present(force)) then
         if (force) do_force = .true.
      endif

      if (do_force) then ! usually, when reconstructing

         this%ni_ = ni
         this%nj_ = nj

      else ! check for oddness.

         this%ni_ = ni
         if (mod(ni, 2) == 0) this%ni_ = this%ni_ + 1

         this%nj_ = nj
         if (mod(nj, 2) == 0) this%nj_ = this%nj_ + 1

      endif
      call this%deduceDeltas()
   end subroutine





   subroutine getUnaryDeltas(this, dfIi, dfIj, dfJi, dfJj)
      class(M2DPolygZone_t), intent(inout) :: this
      real(bsa_real_t), intent(out) :: dfIi, dfIj
      real(bsa_real_t), intent(out) :: dfJi, dfJj

      real(bsa_real_t) :: ang, c, s
      real(bsa_real_t), intrinsic :: sin, cos

      if (this%rot_ < CST_PId2) then ! FIRST quadrant

         c = cos(this%rot_)
         s = sin(this%rot_)

         dfJi =  s
         dfJj =  c

         dfIi =  c
         dfIj = -s

      elseif (this%rot_ < CST_PIGREC) then ! SECOND quadrant

         ang = this%rot_ - CST_PId2
         c   = cos(ang)
         s   = sin(ang)

         dfJi = -s   !  c
         dfJj =  c   ! -s

         dfIi = -c   ! -s
         dfIj = -s   ! -c

      elseif (this%rot_ < CST_PIt3d2) then ! THIRD quadrant

         ang = this%rot_ - CST_PIGREC
         c   = cos(ang)
         s   = sin(ang)

         dfJi = -s
         dfJj = -c

         dfIi = -c
         dfIj =  s

      elseif (this%rot_ < CST_PIt2) then ! FOURTH quadrant

         ang = this%rot_ - CST_PIt3d2
         c   = cos(ang)
         s   = sin(ang)

         dfJi =  s   ! -c
         dfJj = -c   !  s

         dfIi =  c   !  s
         dfIj =  s   !  c
      endif
   end subroutine


end module BsaLib_M2DPolygZone
