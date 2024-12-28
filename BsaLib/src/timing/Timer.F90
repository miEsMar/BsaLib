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
module BsaLib_Timing

   use, intrinsic :: iso_fortran_env, only: int64, real64
   implicit none (type, external)
   private
   public :: timingInit

   type, public :: timer_t
      private
      integer(int64) :: t_init_     = 0._int64
      integer(int64) :: t_last_     = 0._int64
   contains
      procedure, public, pass(this)   :: init  => timerInit
      procedure, public, pass(this)   :: time  => timerTime
      procedure, public, nopass       :: total => timerGetTotal
      procedure, public, pass(this)   :: reset => timerReset
   end type timer_t


   integer(int64) :: irate_  = 0_int64
   integer(int64) :: t_start = 0_int64

contains


   subroutine timingInit()
      call system_clock(count_rate=irate_)
      call system_clock(count=t_start)
   end subroutine


   real(real64) function get_seconds(end, ini) result(secs)
      integer(int64), intent(in) :: end, ini

      secs = real(end - ini, kind=real64) / real(irate_, kind=real64)
   end function


   subroutine timerInit(this)
      class(timer_t) :: this

      call system_clock(count=this%t_init_)
      this%t_last_ = this%t_init_
   end subroutine timerInit



   real(real64) function timerTime(this) result(dt)
      class(timer_t) :: this
      integer(int64) :: t_now

      call system_clock(count=t_now)
      dt = get_seconds(t_now, this%t_last_)

      this%t_last_ = t_now
   end function timerTime




   real(real64) function timerGetTotal() result(tot)
      integer(int64) :: t_now

      call system_clock(count=t_now)
      tot = get_seconds(t_now, t_start)
   end function timerGetTotal




   subroutine timerReset(this)
      class(timer_t) :: this

      this%t_init_     = 0._int64
      this%t_last_     = 0._int64
   end subroutine timerReset



end module
