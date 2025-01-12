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
module BsaLib_Utility

   use BsaLib_CONSTANTS, only: INFOMSG, ERRMSG, MSGCONT, WARNMSG, NOTEMSG, int32
   implicit none (type, external)
   public

contains


   function util_createDirIfNotExist(dirname) result(ierr)
      !!
      !! Creates a directory if it does not exists.
      !!
      character(len = *), intent(in)  :: dirname
      character(len = 128) :: cmd
      integer(int32) :: ierr

      ierr = 0
#ifdef __INTEL_COMPILER
      block
         logical :: lflag

         inquire(directory=dirname, exist=lflag)
         if (lflag) then
            print '(1x, 2a)', &
               INFOMSG, 'Directory  "'//dirname//'"  already exists.'
            return
         endif
      end block
#endif

#ifdef _WIN32
      cmd = dirname
      do ierr = 1, len_trim(dirname)
         if (cmd(ierr:ierr) == '/') then
            cmd(ierr:ierr) = '\'
         endif
      enddo
      cmd = 'if not exist "' // trim(cmd) // '"  mkdir  '// trim(cmd)
#else
      cmd = 'mkdir -p ' // dirname
#endif
      call execute_command_line(cmd(1 : len_trim(cmd)), .true., ierr)
      if (ierr == 0) goto 10

      print '(1x, 3a)', &
         ERRMSG,  'Cannot create directory  ', trim(dirname)

      10 return
   end function


   pure elemental function util_getCorrVectIndex(ni, nj, tot) result(id)
      !! Returns the equivalent index when storing spatial nodal
      !! correlation as a vector (avoiding storing duplicates, symmetric)
      !! NOTE: assumes that ni is the leading node, nj >= ni.
      !!       Otherwise, swaps them (makes use of symmetry).
      integer(int32), intent(in) :: ni, nj
      integer(int32), intent(in) :: tot
      integer(int32) :: id

      if (nj >= ni) then
         id = (ni - 1) * tot + nj - int((ni*ni - ni) / 2., kind=int32)
      else
         id = (nj - 1) * tot + ni - int((nj*nj - nj) / 2., kind=int32)
      endif
   end function

end module BsaLib_Utility
