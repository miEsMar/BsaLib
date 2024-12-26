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
submodule(BsaLib_WindData) BsaLib_WindDataImpl

   use BsaLib_Data,      only: bsa_Abort
   implicit none (type, external)

contains



   module subroutine SetWindvertProf(this, iwprof)
      class(WindData_t), intent(inout) :: this
      integer(bsa_int_t), value :: iwprof

      if (iwprof == 5) iwprof = 1
      if (iwprof < 1 .or. iwprof > 3) call bsa_Abort('Invalid "iwprof" value.')

      this%i_wind_prof_ = iwprof

#ifdef BSA_DEBUG
      write(unit_debug_, '(1x, 3a)') &
         DBGMSG, &
         'Wind vertical profile set to  ', trim(CST_WIND_V_PROFILES(this%i_wind_prof_))
#endif
   end subroutine SetWindvertProf




   module subroutine SetMainvertDir(this, ivert)
      class(WindData_t), intent(inout) :: this
      integer(bsa_int_t), value :: ivert

      if (ivert < 1 .or. ivert > 3) call bsa_Abort('Invalid "ivert" value.')
      this%i_vert_ = ivert

#ifdef BSA_DEBUG
      write(unit_debug_, '(1x, 3a)') &
         DBGMSG, &
         'Wind vertical axis set to  ', CST_WIND_VERT_DIRS(this%i_vert_)
#endif
   end subroutine SetMainvertDir




   module subroutine SetWindZoneLimits(this, lim, ilim)
      class(WindData_t), intent(inout) :: this
#if  ((defined(__INTEL_COMPILER_BUILD_DATE)) && (__INTEL_COMPILER_BUILD_DATE >= 20221019))
      real(bsa_real_t), intent(in), target     :: lim(..)
      integer(bsa_int_t), intent(in), optional :: ilim(..)
#else
      real(bsa_real_t), intent(in), target     :: lim(:)
      integer(bsa_int_t), intent(in), optional :: ilim(:)   ! limits' index passed
#endif

#if  ((defined(__INTEL_COMPILER_BUILD_DATE)) && (__INTEL_COMPILER_BUILD_DATE >= 20221019))
      select rank (lim)
      rank (0)

         if (.not. present(ilim)) call bsa_Abort('Please provide "ilim" value.')

         select rank (ilim)
         rank (0)
            if (ilim < 1 .or. ilim > this%nz_ + 1) call bsa_Abort('Invalid "ilim" value.')
            this%limits_wz_(ilim) = lim
         rank default
            call bsa_Abort('"ilim" must be an integer scalar value.')
         endselect

      rank (1)
#endif

         if (present(ilim)) then

#if  ((defined(__INTEL_COMPILER_BUILD_DATE)) && (__INTEL_COMPILER_BUILD_DATE >= 20221019))
            select rank (ilim)
            rank (1)
#endif
               if (.not. (size(lim(:)) == size(ilim(:)))) &
                  call bsa_Abort('sizes of "lim" and "ilim" do not match.')

               this%limits_wz_(ilim(:)) = lim(:)

#if  ((defined(__INTEL_COMPILER_BUILD_DATE)) && (__INTEL_COMPILER_BUILD_DATE >= 20221019))
            rank default
               call bsa_Abort('expecting "ilim" to be a 1-rank array.')
            endselect
#endif

         else

            if (.not. this%nz_ == 0) then
               if (.not. (size(lim(:)) == this%nz_ + 1)) & 
                  call bsa_Abort('size of "lim" does not match number of wind zones.')
            else
               this%nz_ = size(lim(:)) - 1
            endif
            this%limits_wz_ => lim(:)
         endif ! present ilim

#if  ((defined(__INTEL_COMPILER_BUILD_DATE)) && (__INTEL_COMPILER_BUILD_DATE >= 20221019))
      rank default
         call bsa_Abort('Expeting "lim" either to be a scalar or a 1-rank array.')
      endselect
#endif
   end subroutine





   module subroutine SetAirDensity(aird)
      real(bsa_real_t), value :: aird

      if (aird < 0._bsa_real_t) call bsa_Abort('Air density has a negative value.')
      air_dens_ = aird
   end subroutine SetAirDensity




   module subroutine SetGlobalW2G(mat)
      real(bsa_real_t), intent(in) :: mat(3, 3)
      integer(int32)       :: istat
      character(len = 256) :: emsg

      if (.not. allocated(rot_W2G_)) then
         allocate(rot_W2G_(3, 3), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('rot_W2G_', istat, emsg)
      endif
      rot_W2G_ = mat
   end subroutine SetGlobalW2G




   module subroutine SetWZMeanWindVel(this, UBref)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: UBref(this%nz_)

      this%u_mean_ref_wz_ => UBref
   end subroutine SetWZMeanWindVel




   module subroutine SetWZRefAlt(this, Zref)
      class(WindData_t), intent(inout)     :: this
      real(bsa_real_t), target, intent(in) :: Zref(this%nz_)

      this%Zref_wz_ => Zref
   end subroutine SetWZRefAlt




   module subroutine SetTurbWindScales(this, L)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: L(3, 3, this%nz_)

      this%turb_scales_wz_ => L
   end subroutine SetTurbWindScales




   module subroutine SetTurbWindSDT(this, wtstd)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: wtstd(3, this%nz_)

      this%sigmaUVW_wz_ => wtstd
   end subroutine SetTurbWindSDT




   module subroutine SetWindCorrCoeffs(this, corrcoeff)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: corrcoeff(3, 3, this%nz_)

      this%corrCoeffs_wz_ => corrcoeff
   end subroutine SetWindCorrCoeffs




   module subroutine SetWindCorrExpnts(this, correxpn)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: correxpn(3, 3, this%nz_)

      this%corrExp_wz_ => correxpn
   end subroutine SetWindCorrExpnts





   module subroutine SetIncidenceAngles(this, incang)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: incang(this%nz_)

      this%incAng_wz_ => incang
   end subroutine SetIncidenceAngles





   module subroutine SetLocalRotMatW2G(this, rotW2G_L)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: rotW2G_L(3, 3, this%nz_)

      this%rot_LW2G_wz_ => rotW2G_L
   end subroutine SetLocalRotMatW2G




   module subroutine setWindDirections(this, dirs, ndirs)
      class(WindData_t) :: this
      integer(bsa_int_t), intent(in) :: dirs(:)
      integer(bsa_int_t), value, optional :: ndirs
      integer(int32) :: itmp

      itmp = size(dirs)

      if (present(ndirs)) then
         if (itmp /= ndirs) then
            print '(1x, 2a/)', &
               ERRMSG, 'Size mismatch in setting spatial directions.'
            call bsa_Abort()
         endif
      endif

      this%i_ndirs_ = itmp
      this%dirs_    = dirs(1 : itmp)
   end subroutine



   module subroutine setTurbComps(this, tc, ntc)
      class(WindData_t) :: this
      integer(bsa_int_t), intent(in) :: tc(:)
      integer(bsa_int_t), value, optional :: ntc
      integer(int32) :: itmp

      itmp = size(tc)

      if (present(ntc)) then
         if (itmp /= ntc) then
            print '(1x, 2a/)', &
               ERRMSG, 'Size mismatch in setting turbulent components.'
            call bsa_Abort()
         endif
      endif

      this%i_ntc_ = itmp
      this%tc_    = tc(1 : itmp)
   end subroutine





   module subroutine SetTurbCompsAndDirsDefault(this)
      class(WindData_t), intent(inout) :: this
      integer(int32) :: itmp, istat
      character(len = 256) :: emsg

      this%i_ntc_ = 1
      if (.not. allocated(this%tc_)) then
         allocate(this%tc_(1), stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('this % tc_', istat, emsg)
      else
         itmp = size(this%tc_)
         if (itmp /= 1) then
            deallocate(this%tc_)
            allocate(this%tc_(1), stat=istat, errmsg=emsg)
            if (istat /= 0) call allocKOMsg('this % tc_', istat, emsg)
         endif
      endif
      this%tc_(1) = 1

      this%i_ndirs_ = 1
      if (.not. allocated(this%dirs_)) then
         allocate(this%dirs_(1))
      else
         itmp = size(this%dirs_)
         if (itmp /= 1) then
            deallocate(this%dirs_)
            allocate(this%dirs_(1))
         endif
      endif
      this%dirs_(1) = 1

#ifdef BSA_DEBUG
      write(unit_debug_, '(1x, 2a)') &
         INFOMSG, &
         'Default direction and turbulent components set -- ok.'
#endif
   end subroutine SetTurbCompsAndDirsDefault





   module subroutine SetNodalVel(this, Unod)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: Unod(:)

      this%u_node_ => Unod
   end subroutine SetNodalVel




   module subroutine SetNodalWindZones(this, NodWZ)
      class(WindData_t), intent(inout) :: this
      integer(bsa_int_t), target, intent(in) :: NodWZ(:)

      this%wz_node_ => NodWZ
   end subroutine SetNodalWindZones





   module subroutine SetNodalWindAltitudes(this, WnodAlt)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: WnodAlt(:)

      this%wAlt_node_ => WnodAlt
   end subroutine SetNodalWindAltitudes






   module subroutine SetSpatialNodalCorr(this, nodCorr)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: nodCorr(:, :)

      this%nod_corr_ => nodCorr
   end subroutine SetSpatialNodalCorr



   module subroutine getFull2DNodCorrMat(this, nn, nodcorr2d)
      use BsaLib_Utility, only: util_getCorrVectIndex
      class(WindData_t), intent(in)  :: this
      integer(bsa_int_t), value      :: nn
      real(bsa_real_t), allocatable, intent(inout) :: nodcorr2d(:, :)
      integer(int32) :: i_ = 0, j_, id_

      if (.not. associated(this%nod_corr_)) return

      if (.not. allocated(nodcorr2d)) allocate(nodcorr2d(nn, nn), stat=i_)
      if (i_ /= 0) return
      do j_ = 1, nn
         do i_ = 1, nn
            id_ = util_getCorrVectIndex(i_, j_, nn)
            nodcorr2d(i_, j_) = this%nod_corr_(id_, 1)
         enddo
      enddo
   end subroutine




   module subroutine SetWindFCoeffs(this, wfc)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: wfc(:, :, :)

      this%wfc_ => wfc
   end subroutine SetWindFCoeffs





   module subroutine SetPhitimesC(this, phiTc)
      class(WindData_t), intent(inout) :: this
      real(bsa_real_t), target, intent(in) :: phiTc(:, :, :)


      this%phi_times_A_ndegw_ => phiTc
   end subroutine SetPhitimesC




   module subroutine clean(this)
      class(WindData_t) :: this

      if (allocated(rot_W2G_)) deallocate(rot_W2G_)

      if (associated(this%u_node_))            nullify(this%u_node_)
      if (associated(this%wz_node_))           nullify(this%wz_node_)
      if (associated(this%nod_corr_))          nullify(this%nod_corr_)
      if (associated(this%wfc_))               nullify(this%wfc_)
      if (associated(this%phi_times_A_ndegw_)) nullify(this%phi_times_A_ndegw_)

      if (associated(this%sigmaUVW_wz_))    nullify(this%sigmaUVW_wz_)
      if (associated(this%turb_scales_wz_)) nullify(this%turb_scales_wz_)
      if (associated(this%corrCoeffs_wz_))  nullify(this%corrCoeffs_wz_)
      if (associated(this%corrExp_wz_))     nullify(this%corrExp_wz_)
      if (associated(this%u_mean_ref_wz_))  nullify(this%u_mean_ref_wz_)
      if (associated(this%Zref_wz_))        nullify(this%Zref_wz_)
      if (associated(this%rot_LW2G_wz_))    nullify(this%rot_LW2G_wz_)
      if (associated(this%limits_wz_))      nullify(this%limits_wz_)
      if (associated(this%incAng_wz_))      nullify(this%incAng_wz_)

      if (allocated(this%psd_)) deallocate(this%psd_)
   end subroutine clean


end submodule
