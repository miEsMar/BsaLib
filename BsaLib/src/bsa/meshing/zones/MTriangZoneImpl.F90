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
submodule(BsaLib_MTriangZone) BsaLib_MTriangZoneImpl

! #ifndef _BSA_M3MF_ONLY_PREMESH
! # define _BSA_M3MF_ONLY_PREMESH 0
! #else
! # if (_BSA_M3MF_ONLY_PREMESH != 0 && _BSA_M3MF_ONLY_PREMESH != 1)
! #  undef _BSA_M3MF_ONLY_PREMESH
! #  define _BSA_M3MF_ONLY_PREMESH 0
! # endif
! #endif

   use BsaLib_CONSTANTS
   use BsaLib_MPoint,    only: MPoint_t, MPoint, getPointsDistance, assignment(=), operator(==)
   use BsaLib_MZone,     only: MZone_ID, DUmpZone
   use BsaLib_Data,      only: bsa_Abort, msh_max_zone_NPts
   use BsaLib_IO,        only: io_units_bfmdump
   implicit none (type, external)


contains



   !> Default MTriangZone basic constructor
   module function MTriangZone_constr_def(name) result(this)
      character(len = *), intent(in), optional :: name
      type(MTriangZone_t) :: this

   end function




   !> Gets rect base along I-dir
   elemental module function baseI_triang(this) result(res)
      class(MTriangZone_t), intent(in) :: this
      real(bsa_real_t) :: res

      res = getPointsDistance(this%Cpt_, this%Bpt_)
   end function


   !> Gets rect base along J-dir
   elemental module function baseJ_triang(this) result(res)
      class(MTriangZone_t), intent(in) :: this
      real(bsa_real_t) :: res

      res = getPointsDistance(this%Cpt_, this%Apt_)
   end function




   module subroutine deduceDeltas_triang(this)
      class(MTriangZone_t), intent(inout) :: this

      ! NOTE: do nothing
      !       DONE JUST BECAUSE WE NEED TO IMPLEMENT THIS
      !       METHOD SINCE IT IS AN ABSTRACT COMING FROM 
      !       PARENT CLASS.
   end subroutine





   !> Returns total N of zone's meshing points.
   !> NOTE: for the moment, ONLY triangles rectangles supported
   !> NOTE: overrides default inherited from parent class.
   pure module function zoneTotNPts_triang(this) result(npt)
      class(MTriangZone_t), intent(in) :: this
      integer :: npt

      npt = this%ni_ * this%ni_

      ! BUG: might cast to real and back to int ??
      npt = (npt + this%ni_ ) / 2
   end function





   !> Computes PAB angle, using "cosine" rule.
   !> NOTE: this only works if Cpt is the point at which
   !>       the wider angle is located!
   module subroutine setPABangle(this)
      class(MTriangZone_t), intent(inout) :: this

      real(bsa_real_t) :: CA, AB, CB, cang

      CA = getPointsDistance(this%Cpt_, this%Apt_)
      AB = getPointsDistance(this%Apt_, this%Bpt_)
      CB = getPointsDistance(this%Cpt_, this%Bpt_)

      cang = (CA*CA + AB*AB - CB*CB) / (2._bsa_real_t * CA * AB)
      this%PABang_ = acos(cang)
   end subroutine




   !> Used if we want some specific deltas
   !> along I and J directions
   module subroutine adaptToDeltas(this, dfi, dfj)
      class(MTriangZone_t), intent(inout) :: this
      real(bsa_real_t), value :: dfi, dfj

      real(bsa_real_t) :: CA, CB
      integer   :: ni, nj

      CA = this%baseJ()
      CB = this%baseI()

      ! n. of segments
      ni = floor(CB / dfi)
      nj = floor(CA / dfj)

      ! BUG: is this necessary??
      ! readapt deltas to fit in FIXED sides (bases)
      dfi = CB / ni
      dfj = CA / nj

      ! now get actual n. of points (refinements)
      this%ni_ = ni + 1
      this%nj_ = nj + 1
   end subroutine





   !> BUG: this routine might hide a bug.
   !>
   !> Deduce triangle zone rotation, knowing the 
   !> coordinates of the 3 triang points.
   !>
   !> Still, triangle rotation is given by the rotation
   !> needed to be done to reach the CA side, starting from
   !> Y positive axis of the GRS.
   module subroutine deduceRotation(this)
      class(MTriangZone_t), intent(inout) :: this

      real(bsa_real_t) :: dx, dy, ang

      ! BUG: not always 0 when supposed to.
      !      machine precision + finite floating point repr instabilities
      ! NOTE: do not call distance methods 
      !       since WE WANT TO KEEP SIGNS. 
      dx = this%Apt_%freqI() - this%Cpt_%freqI()
      dy = this%Apt_%freqJ() - this%Cpt_%freqJ()

      ang = atan(abs(dx) / abs(dy))

      if (.not. (dx >= 0._bsa_real_t .and. dy >= 0._bsa_real_t)) then ! NOT 1st quadrant

         if (dy < 0._bsa_real_t) then

            if (dx >= 0._bsa_real_t) then ! 2nd quadrant

               ang = CST_PIGREC - ang
            else ! 3rd quadrant

               ang = CST_PIGREC + ang
            endif

         else ! 4th quadrant

            ! BUG: precision instability.
            !      Tru to mitigate it, vary bad!!
            if (ang > 1e-14_bsa_real_t) then ! fixed tolerance acceptance

               ang = CST_PIt2 - ang
            else ! assume it 0

               ang = 0._bsa_real_t
            endif
         endif
      endif

      this%rot_ = ang
   end subroutine deduceRotation





   !> Gets actua points' delta increments based on this zone's
   !> rotation w.r.t GRS, when moving along the mesh.
   !> If NOT inverted:
   !>   - df_var : total delta along CA side (J-dir in LRS)
   !>   - df_cst : total delta along CB side (I-dir in LRS)
   module subroutine getRotatedUnaryDF(this, df_I_var, df_J_var, df_I_cst, df_J_cst)
      class(MTriangZone_t), intent(inout) :: this
      real(bsa_real_t), intent(out)       :: df_I_var, df_J_var
      real(bsa_real_t), intent(out), optional :: df_I_cst, df_J_cst

      logical :: do_invert = .false.
      real(bsa_real_t) :: ang


      if (present(df_I_cst) .and. present(df_J_cst)) do_invert = .true.


      if (this%rot_ <= CST_PId2) then ! 1st quadrant

         df_I_var = sin(this%rot_)
         df_J_var = cos(this%rot_)

      elseif (this%rot_ <= CST_PIGREC) then ! 2nd quadrant

         ang = this%rot_ - CST_PId2
         df_I_var =   cos(ang)
         df_J_var = - sin(ang)

      elseif (this%rot_ <= CST_PIt3d2) then ! 3rd quadrant

         ang = this%rot_ - CST_PIGREC
         df_I_var = - sin(ang)
         df_J_var = - cos(ang)

      elseif (this%rot_ <= CST_PIt2) then ! 4th quadrant

         ang = this%rot_ - CST_PIt3d2
         df_I_var = - cos(ang)
         df_J_var =   sin(ang)

      endif

      ! BUG: needed to avoid machine rounding precision
      if (abs(df_I_var) < MACHINE_PRECISION) df_I_var = 0._bsa_real_t
      if (abs(df_J_var) < MACHINE_PRECISION) df_J_var = 0._bsa_real_t


      if (do_invert) then ! BUG: ?????????

         df_I_cst =   df_J_var
         df_J_cst = - df_I_var

         ! BUG: needed to avoid machine rounding precision
         if (abs(df_I_cst) < MACHINE_PRECISION) df_I_cst = 0._bsa_real_t
         if (abs(df_J_cst) < MACHINE_PRECISION) df_J_cst = 0._bsa_real_t
      endif

   end subroutine getRotatedUnaryDF









   !> Defines a triang zone from 3 triang points. 
   !> No need to specify which one is A or B point, this is 
   !> automatically deduced.
   !> Automatically choses refinements definition.
   module subroutine defineFromPts_norm(this, Cp, P1, P2)
      class(MTriangZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)         :: Cp, P1, P2

      if (this%ni_ == 0 .or. this%nj_ == 0) & 
         call bsa_Abort('Missing refinement information. Aborting.')


      block
         integer(int32) :: ni, nj
         character(len = *), parameter :: msg_segm = &
            ERRMSG//'Triangle collapsed into a segment. Aborting.'

         ni = this%ni_
         nj = this%nj_


         ! Try to deduce which one is A and B points.
         ! NOTE: this relies on the fact that we know
         !       a priori which is the corner point.
         if (P1%freqJ() > P2%freqJ() .or. P1%freqJ() < P2%freqJ()) then ! NE, E, SE || SW, W, NW

            this%Apt_ = P1
            this%Bpt_ = P2

         else ! same FJ

            if (P1%freqJ() > Cp%freqJ()) then ! w2 > 0

               if (P1%freqI() < P2%freqI()) then

                  this%Apt_ = P1
                  this%Bpt_ = P2

               elseif (P2%freqI() < P1%freqI()) then

                  this%Apt_ = P2
                  this%Bpt_ = P1
               else
                  call bsa_Abort(msg_segm)
               endif

            elseif (P1%freqJ() < Cp%freqJ()) then ! w2 < 0

               if (P1%freqI() > P2%freqI()) then

                  this%Apt_ = P1
                  this%Bpt_ = P2

               elseif (P2%freqI() > P1%freqI()) then

                  this%Apt_ = P2
                  this%Bpt_ = P1
               else
                  call bsa_Abort(msg_segm)
               endif

            else
               call bsa_Abort('Three points are aligned on the same line. Aborting.')
            endif
         endif

      end block


      ! NOTE: make a copy
      this%Cpt_ = MPoint(Cp)

      call this%setPABangle()
      call this%deduceRotation()

   end subroutine defineFromPts_norm









   !> Defines a triang zone from 3 triang points. 
   !> No need to specify which one is A or B point, this is 
   !> automatically deduced.
   !> Also, possible to specify either a value.
   !> If "is_refinement", then is either refinement along I or J direction.
   !> NOTE: this is the defult if no optional is passed.
   !>
   !> If "is_refinement" is false, but still a value is passed, 
   !> it is interpreted as a delta value, still along either I or J direction.
   module subroutine defineFromPts_refmtsORdeltas(&
         this, Cp, P1, P2, is_refinement, val_types, val1, val2)
      class(MTriangZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)    :: Cp, P1, P2
      logical, intent(in)            :: is_refinement
      character(len = *), intent(in), optional :: val_types
      real(bsa_real_t), intent(in), optional   :: val1, val2

      integer(int32)   :: ni, nj
      real(bsa_real_t) :: dfi, dfj


      if (is_refinement) then

         if (present(val_types) .and. present(val1) .and. present(val2)) then ! is present, acquire

            block
               integer :: slen, i

               slen = len(val_types)
               if (slen == 0 .or. slen > 2) &
                  call bsa_Abort('"val_types" is either empty or too long.')

               do i = 1, slen

                  if (val_types(i:i) == 'i') then

                     if (this%ni_ == 0) then ! we need to set it
                        if (i == 1) then
                           if (.not.(present(val1))) &
                              call bsa_Abort('Missing base information along I-dir')
                           ni = val1
                        else ! == 2
                           if (.not.(present(val2))) &
                              call bsa_Abort('Missing base information along I-dir')
                           ni = val2
                        endif
                        this%ni_ = ni
                     endif

                  elseif (val_types(i:i) == 'j') then

                     if (this%nj_ == 0) then ! we need to set it
                        if (i == 1) then
                           if (.not.(present(val1))) &
                              call bsa_Abort('Missing base information along J-dir')
                           nj = val1
                        else ! == 2
                           if (.not.(present(val2))) &
                              call bsa_Abort('Missing base information along J-dir')
                           nj = val2
                        endif
                        this%nj_ = nj
                     endif

                  endif
               enddo

            end block

         endif ! optional vals present


      else ! use deltas


         if (.not. (present(val_types) .and. present(val1) .and. present(val2))) &
            call bsa_Abort('Please provide both deltas.')

         block
            integer :: slen, i

            slen = len(val_types)
            if (slen == 0 .or. slen > 2) &
               call bsa_Abort('"val_types" is either empty or too long.')

            do i = 1, slen

               if (val_types(i:i) == 'i') then

                  if (i == 1) then
                     if (.not.(present(val1))) &
                        call bsa_Abort('Missing delta information along I-dir')
                     dfi = val1
                  else ! == 2
                     if (.not.(present(val2))) &
                        call bsa_Abort('Missing delta information along I-dir')
                     dfi = val2
                  endif

               elseif (val_types(i:i) == 'j') then

                  if (i == 1) then
                     if (.not.(present(val1))) &
                        call bsa_Abort('Missing delta information along J-dir')
                     dfj = val1
                  else ! == 2
                     if (.not.(present(val2))) &
                        call bsa_Abort('Missing delta information along J-dir')
                     dfj = val2
                  endif

               endif
            enddo
         end block

      endif ! use refinements / deltas

      ! Common definition.
      call this%defineFromPts_norm(Cp, P1, P2)

      if (.not. is_refinement) call this%adaptToDeltas(dfi, dfj)
   end subroutine defineFromPts_refmtsORdeltas








   !> Computes a completely defined triang zone
   module subroutine compute(this)
      class(MTriangZone_t), intent(inout) :: this

      real(bsa_real_t) :: df_CA, df_CB


      ! get (full) deltas in straight directions
      df_CA = getPointsDistance(this%Cpt_, this%Apt_) ! base J
      df_CB = getPointsDistance(this%Cpt_, this%Bpt_) ! base I


      ! BUG: Allow for different triangle sides (NOT isosceles)

      ! BUG: it might be that they're not equal 
      !      when they should. So, no direct comparison
      !      instead use MACHINE PRECISION
      if (abs(df_CA - df_CB) <= MACHINE_PRECISION) then
         call compute_TZ_iso_(this)
      else
         call bsa_Abort('Triangle must be "isosceles". Aborting.')
      endif
   end subroutine compute





   module subroutine dumpTZ(this)
      !! Dumps a triang zone data for later reconstruction
      class(MTriangZone_t), intent(in) :: this

      write(io_units_bfmdump(1)) MZone_ID%TRIANGLE

      ! 3 pts
      write(io_units_bfmdump(1)) this%Cpt_%freqI(), this%Cpt_%freqJ()
      write(io_units_bfmdump(1)) this%Apt_%freqI(), this%Apt_%freqJ()
      write(io_units_bfmdump(1)) this%Bpt_%freqI(), this%Bpt_%freqJ()

      ! NOTE: useless, since rot might reconstructed from points
      write(io_units_bfmdump(1)) this%rot_
      write(io_units_bfmdump(1)) this%ni_, this%nj_
   end subroutine dumpTZ




   module subroutine undumpTZ(this)
      !! Undumps a triang zone data for later reconstruction
      class(MTriangZone_t), intent(inout) :: this
      real(bsa_real_t) :: rval1, rval2

      ! 3 pts
      read(io_units_bfmdump(1)) rval1, rval2
      call this%Cpt_%setFreqs(rval1, rval2)

      read(io_units_bfmdump(1)) rval1, rval2
      call this%Apt_%setFreqs(rval1, rval2)

      read(io_units_bfmdump(1)) rval1, rval2
      call this%Bpt_%setFreqs(rval1, rval2)


      read(io_units_bfmdump(1)) this%rot_
      read(io_units_bfmdump(1)) this%ni_, this%nj_
   end subroutine undumpTZ



   elemental function getTriangZoneEquivNPts(ni, nj) result(npt)
      !! Returns Triang zone n of points if it had ni-nj refinements
      integer(int32), intent(in) :: ni, nj
      integer(int32) :: npt

      ! NOTE: for the moment ni==nj
      npt = ni * nj

      ! BUG: do we really need to divide by real??
      npt = int((npt + ni) / 2._real32, kind=int32)
   end function





#if (defined(BSA_USE_POD_DATA_CACHING)) || (defined(_OPENMP))
# define __new_interp_proc__
#endif

   module subroutine interpolateTZ( this &
#ifndef BSA_USE_POD_DATA_CACHING
# define __bfm_undump__ bfm, 
      & , bfm &
#else
# define __bfm_undump__
#endif
      & , pdata )
      !! Implementation of triang zone interpolation methods
      class(MTriangZone_t), intent(inout) :: this
#ifndef BSA_USE_POD_DATA_CACHING
      real(bsa_real_t), intent(in)  :: bfm(:, :)
#endif
      class(*), pointer, intent(in) :: pdata

      ! NOTE: for the moment only supporting HTPC method
      call interpolateTZ_HTPC_v3(this, __bfm_undump__  pdata)
   end subroutine



#define __triang_zone__
#include "_shared_poly2d.fi"

end submodule
