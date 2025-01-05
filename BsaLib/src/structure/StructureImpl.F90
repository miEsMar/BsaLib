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
submodule(BsaLib_Structure) BsaLib_StructDataImpl

   use BsaLib_IO,        only: unit_debug_, allocKOMsg, deallocKOMsg
   use BsaLib_Data,      only: bsa_Abort
   use BsaLib_CONSTANTS, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG
   implicit none (type, external)


contains


   module subroutine SetNodalCoords(this, nn, coords)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)        :: nn
      real(bsa_real_t), target, contiguous  :: coords(:, :)

      if (this%nn_ == 0) then
         this%nn_ = nn
      else
         if (this%nn_ /= nn) &
            call bsa_Abort('Nodal info does not match in setting nodal coordinates. Check again.')
      endif

      this%coords_ => coords
   end subroutine






   module subroutine SetNOfNodalDOFs(this, nlibs)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)        :: nlibs

      this%nlibs_ = nlibs
   end subroutine SetNOfNodalDOFs






   module subroutine SetTotalNOfNodes(this, nn)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in)        :: nn

      this%nn_ = nn
   end subroutine SetTotalNOfNodes




   module subroutine SetLoadedNodalDOFs(this, nlib, lib)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in) :: nlib
      integer(bsa_int_t), target, intent(in) :: lib(:)


      this%nlibs_load_ = nlib
      this%libs_load_  => lib
   end subroutine SetLoadedNodalDOFs





   module subroutine SetLoadedNodes(this, nnl, nl)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in) :: nnl
      integer(bsa_int_t), target, intent(in) :: nl(:)

      this%nn_load_ = nnl
      this%n_load_  => nl
   end subroutine SetLoadedNodes





   module subroutine SetModalInfo(this, ndofs, nm, Phi, natf)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in) :: ndofs, nm
      real(bsa_real_t), intent(in), target :: Phi(ndofs, nm), natf(nm)

      this%ndofs_     = ndofs
      this%modal_%nm_ = nm
      this%modal_%phi_       => Phi
      this%modal_%nat_freqs_ => natf
   end subroutine SetModalInfo





   module subroutine SetKeptModes(this, modes)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in) :: modes(:)
      integer(int32) :: istat, nmodes
      character(len = 256) :: emsg

      nmodes = size(modes)
      if (this%modal_%nm_eff_ == 0) then
         this%modal_%nm_eff_ = nmodes
      else
         if (nmodes /= this%modal_%nm_eff_) call bsa_Abort('Sizes do not match.')
      endif


      if (.not. allocated(this%modal_%modes_)) then
         allocate(this%modal_%modes_(this%modal_%nm_eff_), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('this % modal_ % modes_', istat, emsg)
      endif

      this%modal_%modes_ = modes
   end subroutine SetKeptModes




   module subroutine SetKeptModesDefault(this)
      class(StructureData_t), intent(inout) :: this
      character(len = 256) :: emsg
      integer(int32) :: istat

      if (this%modal_%nm_ == 0) call bsa_Abort('Trying to allocate modes when NM==0 yet.')

      if (.not. allocated(this%modal_%modes_)) then
         allocate(this%modal_%modes_(this%modal_%nm_), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('this % modal_%modes_', istat, emsg)
      endif

      do concurrent (istat = 1 : this%modal_%nm_)
         this%modal_%modes_(istat) = istat
      enddo
   end subroutine SetKeptModesDefault





   module subroutine SetModalMatrices(this, nm, Mg, Kg, Cg)
      class(StructureData_t), intent(inout) :: this
      integer(bsa_int_t), intent(in) :: nm
      real(bsa_real_t), intent(in), target :: Mg(nm), Kg(nm)
      real(bsa_real_t), intent(in), target :: Cg(nm, nm)

      if (.not. this%modal_%nm_ == 0) then
         if (.not. this%modal_%nm_ == nm) &
            call bsa_Abort('You passed a value of "nm" which differs from previously set.')
      else ! == 0
         this%modal_%nm_  = nm
      endif

      this%modal_%Mm_ => Mg
      this%modal_%Km_ => Kg
      this%modal_%Cm_ => Cg

   end subroutine SetModalMatrices





   module subroutine SetTotDamping(this, xsi)
      class(StructureData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: xsi(this%modal_%nm_)

      this%modal_%xsi_ => xsi
   end subroutine





   module subroutine ComputeResPeakWidths(this)
      class(StructureData_t), intent(inout) :: this
      integer(int32) :: istat
      character(len = 256) :: emsg

      if (allocated(this%res_peak_width_)) then

         if (.not. all(this%res_peak_width_ == 0._bsa_real_t)) return
      else

         allocate(this%res_peak_width_(this%modal_%nm_), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('this % res_peak_width_', istat, emsg)
      endif

      this%res_peak_width_ = this%modal_%xsi_ * this%modal_%nat_freqs_

#ifdef BSA_DEBUG
      write(unit_debug_, *) DBGMSG, 'Resonant peak widths = '
      write(unit_debug_, '(*(g12.5, 2x))') this%res_peak_width_
#endif
   end subroutine ComputeResPeakWidths






   module subroutine computeBKGPeakWidths(this, wind_scales)
      class(StructureData_t), intent(inout) :: this
      real(bsa_real_t), intent(in) :: wind_scales(:, :)
      integer   :: j, i
      integer(int32) :: istat
      character(len = 256) :: emsg


      if (.not. allocated(this%bkg_peak_width_)) then
         allocate(this%bkg_peak_width_(3, 3), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('this % bkg_peak_width_', istat, emsg)
      endif

      this%bkg_peak_width_ = 0._bsa_real_t

#ifdef BSA_DEBUG
      do j = 1, 3
         do i = 1, 3
            if (wind_scales(i, j) == 0._bsa_real_t) cycle
            this%bkg_peak_width_(i, j) = 1._bsa_real_t / wind_scales(i, j)
         enddo
      enddo
#else
      do concurrent (j = 1:3, i = 1:3)
         if (.not. wind_scales(i, j) == 0._bsa_real_t) then
            this%bkg_peak_width_(i, j) = 1._bsa_real_t / wind_scales(i, j)
         endif
      enddo
#endif

#ifdef BSA_DEBUG
      write(unit_debug_, *) &
         DBGMSG, '@StructImpl::computeBKGPeakWidths() : bkg peak widths = '
      write(unit_debug_, '(3(g12.5, 2x))') this%bkg_peak_width_
#endif
   end subroutine computeBKGPeakWidths








   module subroutine clean(this)
      class(StructureData_t) :: this
      integer(int32) :: istat
      character(len = 256) :: emsg


      if (associated(this%n_load_))    nullify(this%n_load_)
      if (associated(this%libs_load_)) nullify(this%libs_load_)
      if (associated(this%coords_))    nullify(this%coords_)


      if (allocated(this%modal_%modes_)) then
         deallocate(this%modal_%modes_, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('this%modal_%modes_', istat, emsg)
      endif

      if (associated(this%modal_%Cm_))        nullify(this%modal_%Cm_)
      if (associated(this%modal_%Km_))        nullify(this%modal_%Km_)
      if (associated(this%modal_%Mm_))        nullify(this%modal_%Mm_)
      if (associated(this%modal_%nat_freqs_)) nullify(this%modal_%nat_freqs_)
      if (associated(this%modal_%phi_))       nullify(this%modal_%phi_)
      if (associated(this%modal_%xsi_))       nullify(this%modal_%xsi_)


      if (allocated(this%str_time_scales_)) then
         deallocate(this%str_time_scales_, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('this%str_time_scales_', istat, emsg)
      endif

      if (allocated(this%bkg_peak_width_)) then
         deallocate(this%bkg_peak_width_, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('this%bkg_peak_width_', istat, emsg)
      endif

      if (allocated(this%res_peak_width_)) then
         deallocate(this%res_peak_width_, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('this%res_peak_width_', istat, emsg)
      endif
   end subroutine clean


end submodule
