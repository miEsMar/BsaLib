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
module BsaLib_IO

   use BsaLib_CONSTANTS
   implicit none (type, external)
   public
   private :: setDefFileNameFromUnitNum_


   !**************************************************************************
   !  I/O  UNITs (mutables)
   !**************************************************************************
   ! dumpfile
   integer(int32) :: unit_dump_bfm_ = 999_int32
   integer(int32) :: un_export_bisp_cls_ = 1203_int32
   integer(int32) :: un_export_bisp_msh_ = 1204_int32

   ! debug
   integer(int32)      :: unit_debug_ = 99999_int32
   character(len = 64) :: BSA_DEBUG_FNAME = 'bsadebug.bsa'


   !**************************************************************************
   !  EXPORTING STATE VARIABLES
   !**************************************************************************
   logical                         :: export_in_cwd_ = .true.
   character(len = :), allocatable :: exp_dir_     ! either export, or default (out) if none set.

   character(len = :), private, allocatable :: export_file_access_
   character(len = :), private, allocatable :: export_file_action_
   character(len = :), private, allocatable :: export_file_async_
   character(len = :), private, allocatable :: export_file_form_
   character(len = :), private, allocatable :: export_file_position_
   character(len = :), private, allocatable :: export_file_status_


   interface
      module subroutine allocKOMsg(name_, istat, emsg)
         character(len = *), intent(in) :: name_, emsg
         integer, intent(in) :: istat
      end subroutine

      module subroutine deallocKOMsg(name_, istat, emsg)
         character(len = *), intent(in) :: name_, emsg
         integer, intent(in) :: istat
      end subroutine

      module function io_readBsaInputParams() result(ret)
         integer(bsa_int_t) :: ret
      end function

      module subroutine io_printUserData()
      end subroutine
   end interface



contains


   subroutine io_setExportAppendMode(imode)
      integer(int32), intent(in) :: imode

      if (imode == BSA_EXPORT_MODE_REPLACE) then
         export_file_status_   = IO_STATUS_REPLACE  ! overrides if exists
      elseif (imode == BSA_EXPORT_MODE_APPEND) then
         export_file_position_ = IO_POSITION_APPEND
      endif
   end subroutine


   subroutine io_setExportDirectory(dirname)
      !! Sets export directory to a different path than outdir.
      character(len = *), intent(in) :: dirname

      exp_dir_ = io_appendFilesep(dirname)
   end subroutine


   subroutine io_setExportInCurrDir()
      export_in_cwd_ = .true.
   end subroutine




   subroutine io_setExportFileFormat(iform)
      integer(int32), intent(in) :: iform

      if (iform == BSA_EXPORT_FORMAT_FORMATTED) then
         export_file_form_ = IO_FORM_FORMATTED
      elseif (iform == BSA_EXPORT_FORMAT_UNFORMATTED) then
         export_file_form_ = IO_FORM_UNFORMATTED
      endif
   end subroutine




   subroutine io_exportMomentToFile(fname, vec, form)
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in)   :: vec(:)
      character(len = *), intent(in), optional :: form
      integer(int32) :: iun, i, dim

      if (present(form)) export_file_form_ = form

      iun = io_openExportFileByName(fname)
      if (iun == 0) return
      dim = size(vec)
      if (export_file_form_ == IO_FORM_FORMATTED) then
         write(iun, *) dim
         do i = 1, dim
            write(iun, *) vec(i)
         enddo
      else
         write(iun) dim
         do i = 1, dim
            write(iun) vec(i)
         enddo
      endif
      close(iun)

      if (present(form)) call io_setExportDefaultSpecifiers()
   end subroutine




   subroutine io_setExportSpecifiers()
      if (.not. allocated(export_file_access_))   export_file_access_   = IO_ACCESS_SEQUEN
      if (.not. allocated(export_file_action_))   export_file_action_   = IO_ACTION_WRITE
      if (.not. allocated(export_file_async_))    export_file_async_    = IO_ASYNC_NO
      if (.not. allocated(export_file_form_))     export_file_form_     = IO_FORM_FORMATTED
      if (.not. allocated(export_file_position_)) export_file_position_ = IO_POSITION_ASIS
      if (.not. allocated(export_file_status_))   export_file_status_   = IO_STATUS_UNKNOWN
   end subroutine




   subroutine io_setExportDefaultSpecifiers()
      export_file_access_   = IO_ACCESS_SEQUEN
      export_file_action_   = IO_ACTION_WRITE
      export_file_async_    = IO_ASYNC_NO
      export_file_form_     = IO_FORM_FORMATTED
      export_file_position_ = IO_POSITION_ASIS
      export_file_status_   = IO_STATUS_UNKNOWN
   end subroutine




   function io_openExportFileByName(file) result(iun)
      !! Opens a file, returning its intenal integer unit descriptor.
      character(len = *), intent(in) :: file
      integer(int32) :: iun
      integer(int32) :: ierr_

      if (allocated(exp_dir_)) then
         open(newunit=iun, file=exp_dir_//file  &
            , iostat=ierr_                      &
            , access=export_file_access_        &
            , action=export_file_action_        &
            , asynchronous=export_file_async_   &
            , form=export_file_form_            &
            , position=export_file_position_    &
            , status=export_file_status_)
      else
         open(newunit=iun, file=file          &
            , iostat=ierr_                    &
            , access=export_file_access_      &
            , action=export_file_action_      &
            , asynchronous=export_file_async_ &
            , form=export_file_form_          &
            , position=export_file_position_  &
            , status=export_file_status_)
      endif


      if (ierr_ == 0) return

      iun = 0
      print '(/ 1x, 3a, """.")', &
         ERRMSG,  '@IO::io_openExportFileByName() : error trying opening export file  "', file
      print '(1x, 2a, i0)', &
         MSGCONT, 'Exiting with status ', ierr_
   end function



   subroutine io_getVerifiedFile(iun, fname, openfile)
      integer(int32), intent(inout)     :: iun
      character(len = *), intent(inout) :: fname
      logical, intent(in), optional     :: openfile
      logical :: is_opn

      ! BUG: maybe throw an error
      if (iun == 0) return

      ! VERIFY UNIT
      ! IF is open, it mean it is already in use -> use it!
      ! BUG: get its actual name to avoid unwanted outcomes!
      inquire(unit=iun, opened=is_opn)
      if (is_opn) then
         inquire(unit=iun, name=fname)
         return
      endif

      ! Unit available. Verify filename (and unit by consequence)
      inquire(file=fname, opened=is_opn)
      if (is_opn) then
         do while (is_opn)
            fname = setDefFileNameFromUnitNum_(iun)
            inquire(file=fname, opened=is_opn)
            if (.not. is_opn) inquire(unit=iun, opened=is_opn)
            iun = iun + 1
         enddo
         iun = iun - 1
      endif

      ! BUG: customise ??
      if (present(openfile)) then
         if (openfile) then
            open(unit=iun                    &
               , file=fname                  &
               , form=IO_FORM_UNFORMATTED    &
               , access=IO_ACCESS_STREAM     &
               , status=IO_STATUS_REPLACE    &
               , action=IO_ACTION_READWRITE)
         endif
      endif
   end subroutine




   function setDefFileNameFromUnitNum_(iun) result(fname)
      integer(int32), intent(in) :: iun
      character(len = 64) :: fname
      character(len = 64) :: tmpfname

      write(unit=tmpfname, fmt='(a, i0, a)') BSA_OUT_FILENAME_PREFIX_DEFAULT_, iun, '.bsa'
      fname = tmpfname(1 : len_trim(tmpfname))
   end function




   function io_appendFilesep(path) result(res)
      character(len = *), intent(in)  :: path
      character(len = :), allocatable :: res
      character(len = 1), parameter   :: filesep = &
#ifdef _WIN32
         & '\'
#else
         & '/'
#endif
      integer(int32) :: ilen

      ilen = len_trim(path)
      if (path(ilen:ilen) == filesep) then
         res = path
      else
         res = path//filesep
      endif
   end function


end module
