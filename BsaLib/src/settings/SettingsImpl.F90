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
submodule(BsaLib_Settings) BsaLib_SettingsImpl

   use BsaLib_IO,        only:  unit_debug_
   use BsaLib_Data,      only:  bsa_Abort
   use BsaLib_CONSTANTS, only:  INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG, BSA_SPATIAL_SYM_NONE
   implicit none (type, external)

contains



   module subroutine SetSubanType(this, isuban)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)   :: isuban

      if (isuban < 0 .or. isuban > 3) call bsa_Abort('Invalid "sub-an" value.')
      this%i_suban_type_ = isuban
   end subroutine


   module subroutine SetVersion(this, ivers)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)   :: ivers

      if (ivers < 0 .or. ivers > 2) call bsa_Abort('Invalid "ivers" value.')
      this%i_vers_ = ivers
   end subroutine


   module subroutine SetScalingType(this, idefsc)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)   :: idefsc

      if (idefsc < 0 .or. idefsc > 2) call bsa_Abort('Invalid "idefsc" value.')
      this%i_def_scaling_ = idefsc
   end subroutine


   module subroutine ActivateSpectraComputation(this, ipsd, ibisp)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), value :: ipsd, ibisp

      if (ipsd < 0 .or. ipsd > 1) then
         print '(1x, 2a)', WARNMSG, "Invalid  ipsd  value. Set to DEFAULT (1)"
         ipsd = 1
      endif
      this%i_compute_psd_ = ipsd

      if (ibisp < 0 .or. ibisp > 1) then
         print '(1x, 2a)', WARNMSG, "Invalid  ibisp  value. Set to DEFAULT (1)"
         ibisp = 1
      endif
      this%i_compute_bisp_ = ibisp
   end subroutine



   module subroutine SetExtension(this, ionlydiag)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)   :: ionlydiag

      if (ionlydiag < 0 .or. ionlydiag > 1) call bsa_Abort('Invalid "ionlydiag" value.')
      this%i_only_diag_ = ionlydiag
   end subroutine



   module subroutine TestMode(this, itest)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)   :: itest

      if (itest < 0 .or. itest > 1) call bsa_Abort('Invalid "itest" value.')
      this%i_test_mode_ = itest
   end subroutine



   module subroutine setClsSettings(this, nfreqs, df)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)   :: nfreqs
      real(bsa_real_t), intent(in) :: df

      if (nfreqs <= 0) call bsa_Abort('Invalid "nfreqs" value.')
      this%nfreqs_ = nfreqs

      if (df <= 0._bsa_real_t) call bsa_Abort('Invalid "df" value.')
      this%df_ = df
   end subroutine



   module subroutine SetMshrSetts(this, isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
      class(settings_t), intent(inout) :: this
      integer(bsa_int_t), value :: isvd, bkgrfmt, ifcov, idumpmod
      real(bsa_real_t),   value :: bkgaext, genpaext, maxaext

      this%i_use_svd_       = isvd
      this%bkg_base_rfmnt_  = bkgrfmt
      this%bkg_area_ext_    = bkgaext
      this%peak_area_ext_   = genpaext
      this%max_area_ext_    = maxaext
      this%i_full_coverage_ = ifcov
      this%i_dump_modal_    = idumpmod
   end subroutine


end submodule
