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

   use, intrinsic :: iso_fortran_env, only: real64
   implicit none (type, external)
   private

   type, public :: timer_t
      private
      real(real64) :: t_init_     = 0._real64
      real(real64) :: t_last_     = 0._real64
      real(real64) :: t_tot_      = 0._real64
      real(real64) :: t_tot_prev_ = 0._real64
   contains
      procedure, public, pass(this) :: init  => InitTimer
      procedure, public, pass(this) :: time  => ClockTimer
      procedure, public, pass(this) :: total => GetTimerTotal
      procedure, public, pass(this) :: reset => ResetTimer
   end type timer_t


contains


   subroutine InitTimer(this)
      class(timer_t) :: this

      call cpu_time(this%t_init_)
      this%t_last_ = this%t_init_
   end subroutine InitTimer




   function ClockTimer(this) result(dt)
      class(timer_t) :: this
      real(real64) :: dt_tmp
      real(real64) :: dt

      call cpu_time(dt_tmp)
      dt = dt_tmp - this%t_last_

      ! update total
      this%t_tot_ = this%t_tot_ + dt

      ! update last cpu_time call
      this%t_last_ = dt_tmp
   end function ClockTimer





   pure elemental function GetTimerTotal(this) result(tot)
      class(timer_t), intent(in) :: this
      real(real64) :: tot

      tot = this%t_tot_
   end function GetTimerTotal




   subroutine ResetTimer(this)
      class(timer_t) :: this

      this%t_init_     = 0._real64
      this%t_last_     = 0._real64
      this%t_tot_      = 0._real64
      this%t_tot_prev_ = 0._real64
   end subroutine ResetTimer



end module
