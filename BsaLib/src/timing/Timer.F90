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
   public

   integer(int64), private :: irate_  = 0_int64
   integer(int64), private :: t_start = 0_int64

contains


   subroutine timing_init()
      call system_clock(count_rate=irate_)
      call system_clock(count=t_start)
   end subroutine



   integer(int64) function timing_clock() result(t)
      call system_clock(count=t)
   end function



   real(real64) function timing_getElapsedSeconds(end, ini) result(secs)
      integer(int64), intent(in) :: end, ini

      secs = real(end - ini, kind=real64) / real(irate_, kind=real64)
   end function


   real(real64) function timing_total() result(tot)
      integer(int64) :: t_now

      call system_clock(count=t_now)
      tot = timing_getElapsedSeconds(t_now, t_start)
   end function


end module
