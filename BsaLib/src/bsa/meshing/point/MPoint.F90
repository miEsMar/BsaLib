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
module BsaLib_MPoint

   use BsaLib_CONSTANTS
   implicit none (type, external)
   private
   public :: getPointsDistance, MPoint


   type, public :: MPoint_t
   
      real(bsa_real_t), private :: fi_ = 0._bsa_real_t
      real(bsa_real_t), private :: fj_ = 0._bsa_real_t
   
   contains

      procedure, pass :: setFreqs
      procedure, pass :: move
      procedure, pass :: freqI
      procedure, pass :: freqJ
      procedure, pass :: getDistanceI => getDistanceIfromCoord, getDistanceIfromPt
      procedure, pass :: getDistanceJ => getDistanceJfromCoord, getDistanceJfromPt
      procedure, pass :: getNewPointFromDistAndRot
      procedure, pass :: scale => scaleInt, scaleReal
   end type MPoint_t

   ! class constructors
   interface MPoint
      module procedure MPoint_as_compiler
      module procedure MPoint_from_ints
      module procedure MPoint_from_MPoint
   end interface


   interface operator(==)
      module procedure PointsAreEqual
   end interface
   public :: operator(==)

   interface operator(/=)
      module procedure PointsAreNotEqual
   end interface
   public :: operator(/=)


   interface operator(*)
      module procedure ScaleByInt
      module procedure ScaleByReal
   end interface
   public :: operator(*)


   interface assignment(=)
      module procedure assignFromPt
      module procedure assignFromInt
      module procedure assignFromReal
   end interface
   public :: assignment(=)



contains



   pure function MPoint_as_compiler(fi, fj) result(this)
      real(bsa_real_t), intent(in) :: fi, fj
      type(MPoint_t) :: this

      ! BUG: maybe here needed to use rounding precision
      this%fi_ = fi
      this%fj_ = fj
   end function


   pure function MPoint_from_ints(fi, fj) result(this)
      integer, intent(in) :: fi, fj
      type(MPoint_t) :: this

      this = MPoint_as_compiler(real(fi, bsa_real_t), real(fj, bsa_real_t))
   end function 


   pure function MPoint_from_MPoint(p) result(this)
      class(MPoint_t), intent(in) :: p
      type(MPoint_t) :: this

      ! make a copy of instance member variables.
      this%fi_ = p%fi_
      this%fj_ = p%fj_
   end function







   elemental function freqI(this) result(val)
      class(MPoint_t), intent(in) :: this
      real(bsa_real_t) :: val

      val = this%fi_
   end function


   elemental function freqJ(this) result(val)
      class(MPoint_t), intent(in) :: this
      real(bsa_real_t) :: val

      val = this%fj_
   end function



   subroutine setFreqs(this, fi, fj)
      class(MPoint_t), intent(inout) :: this
      real(bsa_real_t), intent(in) :: fi, fj

      this%fi_ = fi
      this%fj_ = fj
   end subroutine




   elemental function getPointsDistance(p1, p2) result(dist)
      class(MPoint_t), intent(in) :: p1, p2
      real(bsa_real_t) :: dist

      real(bsa_real_t) :: dx, dy

      dx = abs(p1%fi_ - p2%fi_)
      dy = abs(p1%fj_ - p2%fj_)

      dist = sqrt(dx*dx + dy*dy)
   end function getPointsDistance





   !> Moves a point by specified x-y deltas.
   subroutine move(this, di, dj)
      class(MPoint_t), intent(inout) :: this
      real(bsa_real_t), intent(in) :: di, dj

      this%fi_ = this%fi_ + di
      this%fj_ = this%fj_ + dj
   end subroutine move



   elemental function getDistanceIfromCoord(this, i_coord) result(dist)
      class(MPoint_t), intent(in) :: this
      real(bsa_real_t), intent(in) :: i_coord
      real(bsa_real_t) :: dist

      dist = abs(this%fi_ - i_coord)
   end function getDistanceIfromCoord

   elemental function getDistanceIfromPt(this, p) result(dist)
      class(MPoint_t), intent(in) :: this
      class(MPoint_t), intent(in) :: p
      real(bsa_real_t) :: dist

      dist = abs(this%fi_ - p%fi_)
   end function getDistanceIfromPt


   elemental function getDistanceJfromCoord(this, j_coord) result(dist)
      class(MPoint_t), intent(in) :: this
      real(bsa_real_t), intent(in) :: j_coord
      real(bsa_real_t) :: dist

      dist = abs(this%fj_ - j_coord)
   end function getDistanceJfromCoord

   elemental function getDistanceJfromPt(this, p) result(dist)
      class(MPoint_t), intent(in) :: this
      class(MPoint_t), intent(in) :: p
      real(bsa_real_t) :: dist

      dist = abs(this%fj_ - p%fj_)
   end function getDistanceJfromPt



   pure recursive function getNewPointFromDistAndRot(this, dist, rot) result(P)
      !! Returns a new Point located by distance and rotation 
      !! from current Point.
      class(MPoint_t), intent(in) :: this
      real(bsa_real_t), intent(in) :: dist, rot
      type(MPoint_t) :: P

      real(bsa_real_t) :: ang, dI, dJ

      if (rot < CST_PId2) then

         ang = rot
         dI  = sin(ang)
         dJ  = cos(ang)

      elseif (rot < CST_PIGREC) then

         ang = real(rot, 8) - CST_PId2
         dI  = cos(ang)
         dJ  = - sin(ang)

      elseif (rot < CST_PIt3d2) then

         ang = rot - CST_PIGREC
         dI  = - sin(ang)
         dJ  = - cos(ang)

      elseif (rot <= CST_PIt2) then

         ang = rot - CST_PIt3d2
         dI  = - cos(ang)
         dJ  = sin(ang)

      else

         ! NOTE: call this function recursively,
         !       after having subtracted 2*PI.
         ang = rot - CST_PIt2
         P   = this%getNewPointFromDistAndRot(dist, ang)
      endif

      ! scale unitary deltas (given by rotation)
      dI = dI * dist
      dJ = dJ * dist

      ! apply deltas to current point coords
      dI = dI + this%fi_
      dJ = dJ + this%fj_

      ! get new point location
      P = MPoint_t(dI, dJ)
   end function getNewPointFromDistAndRot




   subroutine scaleInt(this, i)
      class(MPoint_t), intent(inout) :: this
      integer, intent(in) :: i

      this%fi_ = this%fi_ * i
      this%fj_ = this%fj_ * i
   end subroutine scaleInt

   subroutine scaleReal(this, i)
      class(MPoint_t), intent(inout) :: this
      real(bsa_real_t), intent(in) :: i

      this%fi_ = this%fi_ * i
      this%fj_ = this%fj_ * i
   end subroutine scaleReal
   





   elemental function PointsAreEqual(p1, p2) result(eq)
      class(MPoint_t), intent(in) :: p1, p2
      logical :: eq

      real(bsa_real_t) :: dfi, dfj

      dfi = abs(p1%fi_ - p2%fi_)
      dfj = abs(p1%fj_ - p2%fj_)
      eq  = (dfi <= MACHINE_PRECISION .and. dfj <= MACHINE_PRECISION)
   end function PointsAreEqual


   elemental function PointsAreNotEqual(p1, p2) result(neq)
      class(MPoint_t), intent(in) :: p1, p2
      logical :: neq

      neq = .not. PointsAreEqual(p1, p2)
   end function PointsAreNotEqual




   
   pure function ScaleByInt(p, i) result(res)
      class(MPoint_t), intent(in) :: p
      integer, intent(in) :: i
      type(MPoint_t) :: res

      res = MPoint_t(p%fi_ * i, p%fj_ * i)
   end function ScaleByInt

   pure function ScaleByReal(p, i) result(res)
      class(MPoint_t), intent(in) :: p
      real(bsa_real_t), intent(in) :: i
      type(MPoint_t) :: res

      res = MPoint_t(p%fi_ * i, p%fj_ * i)
   end function ScaleByReal






   pure subroutine assignFromPt(lhs, rhs)
      type(MPoint_t), intent(out)  :: lhs
      class(MPoint_t), intent(in)  :: rhs

      lhs%fi_ = rhs%fi_
      lhs%fj_ = rhs%fj_
   end subroutine


   pure subroutine assignFromInt(lhs, rhs)
      type(MPoint_t), intent(out)  :: lhs
      integer, intent(in)          :: rhs
      real(bsa_real_t) :: rval

      rval = real(rhs, bsa_real_t)
      lhs  = MPoint_t(rval, rval)
   end subroutine


   pure subroutine assignFromReal(lhs, rhs)
      type(MPoint_t), intent(out)  :: lhs
      real(bsa_real_t), intent(in)        :: rhs

      lhs = MPoint_t(rhs, rhs)
end subroutine


end module BsaLib_MPoint
