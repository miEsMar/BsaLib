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
submodule(BsaLib_MRectZone) BsaLib_MRectZoneImpl

! #ifndef _BSA_M3MF_ONLY_PREMESH
! # define _BSA_M3MF_ONLY_PREMESH 0
! #else
! # if (_BSA_M3MF_ONLY_PREMESH != 0 && _BSA_M3MF_ONLY_PREMESH != 1)
! #  undef _BSA_M3MF_ONLY_PREMESH
! #  define _BSA_M3MF_ONLY_PREMESH 0
! # endif
! #endif

   use BsaLib_CONSTANTS
   use BsaLib_Data,     only: bsa_Abort, msh_max_zone_NPts
   use BsaLib_IO,       only: unit_debug_, io_units_bfmdump
   use BsaLib_MPoint,   only: MPoint, MPoint_t, getPointsDistance
   use BsaLib_MZone,    only: MZone_ID, DumpZone
   implicit none (type, external)


contains


   !> Gets rect base along I-dir
   elemental module function baseI_rct(this) result(res)
      class(MRectZone_t), intent(in) :: this
      real(bsa_real_t) :: res

      res = this%base_I_
   end function


   !> Gets rect base along J-dir
   elemental module function baseJ_rct(this) result(res)
      class(MRectZone_t), intent(in) :: this
      real(bsa_real_t) :: res

      res = this%base_J_
   end function






   !> Get A point.
   !> A point is the point defined, starting from I point, 
   !> along the J-dir (Y-axis) parallel side.
   pure module function getAPoint(this) result(pt)
      class(MRectZone_t), intent(in) :: this
      type(MPoint_t) :: pt

      pt = this%Ipt_%getNewPointFromDistAndRot(this%base_J_, this%rot_)
   end function



   !> Get B point.
   !> B point is the point defined, starting from I point, 
   !> along the I-dir (X-axis) parallel side.
   pure module function getBPoint(this) result(pt)
      class(MRectZone_t), intent(in) :: this
      type(MPoint_t)   :: pt
      real(bsa_real_t) :: ang

      ! in 4-th quadrant, decrement by 3/2*pi -> first quadrant
      if (this%rot_ > CST_PIt3d2) then
         ang = this%rot_ - CST_PIt3d2
      else ! otherwise, increment by 1/2*pi
         ang = this%rot_ + CST_PId2
      endif
      pt = this%Ipt_%getNewPointFromDistAndRot(this%base_I_, ang)
   end function



   module subroutine deduceDeltas_rect(this)
      class(MRectZone_t), intent(inout) :: this
      integer(int32) :: nsegi, nsegj

      !
      this%refmts_set_ = .true.

      ! n. of segments along I-dir
      nsegi = this%ni_ - 1

      ! n. of segments along J-dir
      nsegj = this%nj_ - 1

      this%deltaf_I_ = this%base_I_ / nsegi
      this%deltaf_J_ = this%base_J_ / nsegj
   end subroutine






   module subroutine setDeltas(this, dfi, dfj, adapt)
      class(MRectZone_t), intent(inout) :: this
      real(bsa_real_t), intent(in)      :: dfi, dfj
      logical, intent(in), optional     :: adapt
      logical :: do_adapt = .false.

      if (present(adapt)) then
         if (adapt) do_adapt = .true.
      endif

      this%deltaf_I_   = dfi
      this%deltaf_J_   = dfj
      this%deltas_set_ = .true.
      call this%deduceRefinements(do_adapt)
   end subroutine



   module subroutine deduceRefinements(this, adapt)
      class(MRectZone_t), intent(inout) :: this
      logical, intent(in) :: adapt
      integer(int32) :: ni, nj

      ! NOTE: as such, they mean N. OF SEGMENTS
      ni = floor(this%base_I_ / this%deltaf_I_)
      nj = floor(this%base_J_ / this%deltaf_J_)

      ! NOTE: we want even N. of segments -> odd N. of points
      if (.not. mod(ni, 2) == 0) ni = ni + 1
      if (.not. mod(nj, 2) == 0) nj = nj + 1


      ! BUG: if we imposed deltas greater than bases -> at least 1 segment!
      if (ni == 0) ni = 1
      if (nj == 0) nj = 1


      if (adapt) then ! we readapt deltas, to fit in base

         this%deltaf_I_ = this%base_I_ / ni
         this%deltaf_J_ = this%base_J_ / nj

      else ! recompute bases, such to preserve desired deltas.

         this%base_I_ = this%deltaf_I_ * ni
         this%base_J_ = this%deltaf_J_ * nj
      endif

      ! Now, get actual number of points along sides (odd)
      this%ni_ = ni + 1
      this%nj_ = nj + 1
   end subroutine







   !> Defines zone bases (sides) from given Deltas, specifying refinement.
   !> NOTE: Init and End zone points MUST be known.
   module subroutine defineFromDeltas_refinements(this, pt, loc, dfi, dfj, ni, nj)
      class(MRectZone_t), intent(inout) :: this

      !> Point specification at location 'loc'
      class(MPoint_t), intent(in) :: pt

      !> Location of the specified point.
      !> Can be:
      !>   'i' -> init point
      !>   'c' -> centre point
      !>   'e' -> end point
      character(len=1) :: loc

      !> Delta values
      real(bsa_real_t), intent(in) :: dfi, dfj

      !> Refinements
      integer, value :: ni, nj

      if (.not. this%isGRSAligned()) call bsa_Abort('Rect zone is not GRS aligned.')


      ! NOTE: force them odd
      if (mod(ni, 2) == 0) ni = ni + 1
      if (mod(nj, 2) == 0) nj = nj + 1

      this%ni_ = ni
      this%nj_ = nj

      ! here we refer to number of segments (even)
      ni = ni - 1
      nj = nj - 1

      ! deduce bases
      this%base_I_ = ni * dfi
      this%base_J_ = nj * dfj

      this%refmts_set_ = .true.
      this%deltaf_I_   = dfi
      this%deltaf_J_   = dfj

      call this%define(pt, loc=loc)
   end subroutine






   !> Defines zone bases (sides) from given Deltas, specifying max deltas values.
   !> NOTE: Init and End zone points MUST be known.
   !> If force is present and true, forces deltas to reach max values specified
   !> by 'valI' and 'valJ'. In this case, deltas are then slightly adjusted to fit.
   module subroutine defineFromDeltas_maxvalues(this, pt, loc, dfi, dfj, maxF_i, maxF_j, force, exceed)
      class(MRectZone_t), intent(inout) :: this

      !> Point specification, which location is specified by 'loc'
      class(MPoint_t), intent(in) :: pt

      !> Point location:
      !>   'i' -> Init location
      !>   'c' -> centre (NOT YET IMPLEMENTED)
      !>   'e' -> end    (NOT YET IMPLEMENTED)
      character(len=1), intent(in) :: loc

      !> Delta values
      real(bsa_real_t), value :: dfi, dfj

      !> Max deltas values
      real(bsa_real_t) :: maxF_i, maxF_j

      !> adjusts deltas to max values specified.
      logical, intent(in), optional :: force

      logical, intent(in), optional :: exceed


      if (.not. this%isGRSAligned()) &
         call bsa_Abort('Zone is not aligned with GRS.')

      if (.not. loc == 'i') &
         call bsa_Abort(&
            'Cannot define deltas from max values if given point location is not "i" (Init).')


      block
         character(len = *), parameter :: invalid_max_vals = &
            ERRMSG//'Invalid max values (out of allowed bounds).'
         logical :: do_force
         real(bsa_real_t) :: fi, fj, bi, bj
         integer(int32)   :: ni, nj

         do_force = .false.
         if (present(force)) then
            if (force) do_force = .true.
         endif

         ! initialise
         bi = 0._bsa_real_t
         bj = 0._bsa_real_t
         fi = 0._bsa_real_t
         fj = 0._bsa_real_t
         ni = 0
         nj = 0

         if (do_force) then

            fi = pt%freqI()
            fj = pt%freqJ()

            if (this%rot_ == 0._bsa_real_t) then
               if (maxF_i <= fi .or. maxF_j <= fj) call bsa_Abort(invalid_max_vals)

               bi = maxF_i - fi
               bj = maxF_j - fj

            elseif (this%rot_ == CST_PId2) then ! 1/2 * pi
               if (maxF_i <= fi .or. maxF_j >= fj) call bsa_Abort(invalid_max_vals)

               bj = maxF_i - fi
               bi = fj - maxF_j

            elseif (this%rot_ == CST_PIGREC) then
               if (maxF_i >= fi .or. maxF_j >= fj) call bsa_Abort(invalid_max_vals)

               bi = fi - maxF_i
               bj = fj - maxF_j

            elseif (this%rot_ == CST_PIt3d2) then ! 3/2 * pi
               if (maxF_i >= fi .or. maxF_j <= fj) call bsa_Abort(invalid_max_vals)

               bj = fi - maxF_i
               bi = maxF_j - fj

            end if


            ! invert deltas if needed
            if (.not. (this%rot_ == 0 .or. this%rot_ == CST_PIGREC)) then
               block
                  real(bsa_real_t) :: tmp

                  tmp = dfi
                  dfi = dfj
                  dfj = tmp
               end block
            endif

            ! get approx n. of segments
            ni = ceiling(bi / dfi)
            nj = ceiling(bj / dfj)

            ! get refactored deltas
            dfi = bi / ni
            dfj = bj / nj

            ! actual n. of points (odd)
            ni = ni + 1
            nj = nj + 1



         else ! DO NOT force (keep those deltas)


            block
               real(bsa_real_t) :: di, dj
               logical :: do_exceed = .false.


               ! defaults
               if (present(exceed)) then
                  if (exceed) do_exceed = .true.
               endif
               di = dfi
               dj = dfj
               ni = 1
               nj = 1


               if (this%rot_ == 0._bsa_real_t) then

                  fi = pt%freqI() + di
                  fj = pt%freqJ() + dj

                  if (maxF_i < fi .or. maxF_j < fj) call bsa_Abort(invalid_max_vals)

                  do while (fi <= maxF_i)
                     ni = ni + 1
                     fi = fi + di
                  enddo
                  do while (fj <= maxF_j)
                     nj = nj + 1
                     fj = fj + dj
                  enddo


               elseif (this%rot_ == CST_PId2) then ! 1/2 * pi
                  fi = pt%freqI() + dj
                  fj = pt%freqJ() - di

                  if (maxF_i < fi .or. maxF_j > fj) call bsa_Abort(invalid_max_vals)

                  do while (fi <= maxF_i)
                     nj = nj + 1
                     fi = fi + dj
                  enddo
                  do while (fj >= maxF_j)
                     ni = ni + 1
                     fj = fj - di
                  enddo


               elseif (this%rot_ == CST_PIGREC) then
                  di = - di
                  dj = - dj
                  fi = pt%freqI() + di
                  fj = pt%freqJ() + dj

                  if (maxF_i > fi .or. maxF_j > fj) call bsa_Abort(invalid_max_vals)

                  do while (fi >= maxF_i)
                     ni = ni + 1
                     fi = fi + di
                  enddo
                  do while (fj >= maxF_j)
                     nj = nj + 1
                     fj = fj + dj
                  enddo

                  di = abs(di)
                  dj = abs(dj)


               ! TODO: check this
               elseif (this%rot_ == CST_PIt3d2) then ! 3/2 * pi
                  fi = pt%freqI() - dj
                  fj = pt%freqJ() + di

                  if (maxF_i > fi .or. maxF_j < fj) call bsa_Abort(invalid_max_vals)

                  do while (fi >= maxF_i)
                     nj = nj + 1
                     fi = fi + dj
                  enddo
                  do while (fj <= maxF_j)
                     ni = ni + 1
                     fj = fj - di
                  enddo

               end if ! rot


               if (ni == 1 .or. nj == 1) call bsa_Abort('At least one max value is too small.')


               if (this%rot_ == 0._bsa_real_t) then
                  bi = (ni - 1) * di
                  bj = (nj - 1) * dj

                  if (do_exceed) then
                     if (pt%freqI() + bi < maxF_i) then
                        bi = ni + di
                        ni = ni + 1
                     endif
                     if (pt%freqJ() + bj < maxF_j) then
                        bj = nj * dj
                        nj = nj + 1
                     endif
                  endif


               elseif (this%rot_ == CST_PId2) then
                  bi = (ni - 1) * dj
                  bj = (nj - 1) * di

                  if (do_exceed) then
                     if (pt%freqI() + bj < maxF_i) then
                        bj = nj * di
                        nj = nj + 1
                     endif
                     if (pt%freqJ() - bi > maxF_j) then
                        bi = ni * dj
                        ni = ni + 1
                     endif
                  endif


               elseif (this%rot_ == CST_PIGREC) then
                  bi = (ni - 1) * di
                  bj = (nj - 1) * dj

                  if (do_exceed) then
                     if (pt%freqI() - bi > maxF_i) then
                        bi = ni * di
                        ni = ni + 1
                     endif
                     if (pt%freqJ() - bj > maxF_j) then
                        bj = nj * dj
                        nj = nj + 1
                     endif
                  endif


               elseif (this%rot_ == CST_PIt3d2) then
                  bi = (ni - 1) * dj
                  bj = (nj - 1) * di

                  if (do_exceed) then
                     if (pt%freqI() - bj > maxF_i) then
                        bj = nj * di
                        nj = nj + 1
                     endif
                     if (pt%freqJ() + bi < maxF_j) then
                        bi = ni * dj
                        ni = ni + 1
                     endif
                  endif

               end if ! rot
            end block

         endif ! force

         ! backup actual n. of points
         this%ni_ = ni
         this%nj_ = nj


         ! BUG: this is a code copy
         this%refmts_set_ = .true.
         this%deltaf_I_   = dfi
         this%deltaf_J_   = dfj

         call this%define(pt, loc, bi, bj)
      end block
   end subroutine






   module subroutine defineFromEndPtCoordAndBase_norm(&
      this, Pi, coord_val, coord_ty_ch, baseval, base_dir, called)

      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)       :: Pi
      real(bsa_real_t), intent(in)      :: coord_val
      character(len = 1), intent(in)    :: coord_ty_ch
      real(bsa_real_t), intent(in)      :: baseval
      character(len = 1), intent(in)    :: base_dir
      logical, intent(in)               :: called

      real(bsa_real_t) :: ang
      type(MPoint_t)   :: pt

      if (.not. (coord_ty_ch == 'i' .or. coord_ty_ch == 'j')) &
         call bsa_Abort('Unvalid coordinate type identifier. Must be one of "i"/"j".')

      if (.not. (base_dir == 'i' .or. base_dir == 'j')) &
         call bsa_Abort('Unvalid base direction identifier. Must be one of "i"/"j".')

      ! backup init point
      this%Ipt_ = MPoint(Pi)


      ! get A or B point
      if (base_dir == 'i') then ! B point

         this%base_I_ = baseval

         ! if 4th quadrant, go back to 1st one
         if (this%rot_ >= CST_PIt3d2) then
            ang = this%rot_ - CST_PIt3d2
         else
            ang = this%rot_ + CST_PId2
         endif

      else ! A pt

         this%base_J_ = baseval
         ang = this%rot_
      endif


      ! NOTE: here pt is either A or B, depending on base passed
      pt = this%Ipt_%getNewPointFromDistAndRot(baseval, ang)

      ! deduce missing end point coordinate
      ! From it, compute missing rect base for complete definition
      call this%getOtherBase(pt, base_dir, coord_ty_ch, coord_val)

      if (called) return

      ! once we have both bases, redefine deltas
      ! NOTE: this assumes refinements have been already set
      call this%deduceDeltas()
   end subroutine defineFromEndPtCoordAndBase_norm




   module subroutine defineFromEndPtCoordAndBase_forceDeltas(&
      this, Pi, coord_val, coord_ty_ch, baseval, base_dir, dfi, dfj)

      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)       :: Pi
      real(bsa_real_t), intent(in)      :: coord_val
      character(len = 1), intent(in)    :: coord_ty_ch
      real(bsa_real_t), intent(in)      :: baseval
      character(len = 1), intent(in)    :: base_dir
      real(bsa_real_t), intent(in)      :: dfi, dfj


      call this%defineFromEndPtCoordAndBase_norm(&
         Pi, coord_val, coord_ty_ch, baseval, base_dir, .true.)

      ! forcing deltas
      call this%setDeltas(dfi, dfj, .true.)
   end subroutine defineFromEndPtCoordAndBase_forceDeltas










   module subroutine define(this, pt, loc, base_i, base_j)
      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)       :: pt
      character(len=1), optional, intent(in) :: loc
      real(bsa_real_t), intent(in), optional :: base_i, base_j

      character(len=1) :: location = 'i'

      if (present(loc)) location = loc

      if (location == 'c') then

         block
            real(bsa_real_t) :: bid2, bjd2

            if (present(base_i)) then
               this%base_I_ = base_i
               bid2         = base_i / 2
            else
               bid2         = this%base_I_ / 2
            endif

            if (present(base_j)) then
               this%base_J_ = base_j
               bjd2         = base_j / 2
            else
               bjd2         = this%base_J_ / 2
            endif

            this%Ipt_ = MPoint(pt%freqI() - bid2, pt%freqJ() - bjd2)
            this%Ept_ = MPoint(pt%freqI() + bid2, pt%freqJ() + bjd2)
         end block

      else ! loc == 'i' or 'e'

         ! NOTE: don't forget to update bases if passed!
         if (present(base_i)) this%base_I_ = base_i
         if (present(base_j)) this%base_J_ = base_j

         call this%setIEpts(pt, loc)
      endif
   end subroutine define






   module subroutine setIEpts(this, pt, loc)
      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)       :: pt
      character(len=1), intent(in)      :: loc

      real(bsa_real_t) :: c, s, ang, di, dj

      if (.not. (loc == 'i' .or. loc == 'e')) &
         call bsa_Abort('Invalid point location.')


      if (this%rot_ < CST_PId2) then ! FIRST quadrant

         c = cos(this%rot_)
         s = sin(this%rot_)

         di = this%base_I_ * c + this%base_J_ * s
         dj = this%base_J_ * c - this%base_I_ * s

      elseif (this%rot_ < CST_PIGREC) then ! SECOND quadrant

         ang = this%rot_ - CST_PId2
         c   = cos(ang)
         s   = sin(ang)

         di = this%base_J_ * c - this%base_I_ * s
         dj = - (this%base_J_ * s + this%base_I_ * c)

      elseif (this%rot_ < CST_PIt3d2) then ! THIRD quadrant

         ang = this%rot_ - CST_PIGREC
         c   = cos(ang)
         s   = sin(ang)

         di = - (this%base_J_ * s + this%base_I_ * c)
         dj = - this%base_J_ * c  + this%base_I_ * s


      elseif (this%rot_ < CST_PIt2) then ! FOURTH quadrant

         ang = this%rot_ - CST_PIt3d2
         c   = cos(ang)
         s   = sin(ang)

         di = - this%base_J_ * c + this%base_I_ * s
         dj =   this%base_J_ * s + this%base_I_ * c
      endif



      select case (loc)

      case ('i')
         this%Ipt_ = MPoint(pt)

         this%Ept_ = MPoint(pt%freqI() + di, pt%freqJ() + dj)

      case ('e')
         this%Ept_ = MPoint(pt)

         this%Ipt_ = MPoint(pt%freqI() - di, pt%freqJ() - dj)
      end select
   end subroutine





   ! !> Avoid setting a delta smaller than given limit
   ! elemental module impure subroutine validateDeltas(this, lval)
   !    class(MRectZone_t), intent(inout) :: this
   !    real(bsa_real_t), intent(in) :: lval

   !    real(bsa_real_t) :: dfi, dfj
   !    logical :: coarsen = .false.

   !    dfi = this%deltaf_I_
   !    dfj = this%deltaf_J_
   !    if (dfi < lval .or. dfi > lval) then
   !       dfi     = lval
   !       coarsen = .true.
   !    endif
   !    if (dfj < lval .or. dfj > lval) then
   !       dfj     = lval
   !       coarsen = .true.
   !    endif
   !    if (coarsen) then
   !       print '( 1x, 2a, g10.5, " [Hz])." )', &
   !          WARNMSG, 'Detected at least one zone deltas too  small/big. Setting to optimal value  (', &
   !          lval
   !       call this%setDeltas(dfi, dfj, .true.)
   !    endif
   ! end subroutine




   ! !> Returns euqivalent rotation in the FIRST quadrant.
   ! elemental function getFirstQuadEquivRot(rot) result(rot1)
   !    real(bsa_real_t), intent(in) :: rot
   !    real(bsa_real_t) :: rot1

   !    if (rot < CST_PId2) then ! FIRST quad, keep it
   !       rot1 = rot
   !    elseif (rot < CST_PIGREC) then ! SECOND quad
   !       rot1 = rot - CST_PId2
   !    elseif (rot < CST_3d2) then ! THIRD quad
   !       rot1 = rot - CST_PIGREC
   !    elseif (rot < CST_PIt2) then ! FOURTH quad
   !       rot1 = rot - CST_3d2
   !    endif
   ! end function





   !> Automatically computes the second remaining (unknown) rect base
   !> based on the point's coordinates and the base that we already defined.
   !> BUG: cannot see when a base is negative (i.e. when END pt is "behind" INIT one).
   module subroutine getOtherBase(this, pt, base_dir, known_coord, coord_val)
      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)  :: pt
      character(len=1), intent(in) :: base_dir, known_coord
      real(bsa_real_t), intent(in) :: coord_val

      real(bsa_real_t) :: kd, rot, cd, fi, fj
      type(MPoint_t) :: Pe


      ! NOTE: preinitialise to avoid errors !!!!!!!!!!!!
      cd = 0._bsa_real_t


      if (base_dir == 'i') then ! we passed B point, we know side along I-dir

         if (known_coord == 'i') then ! we search DJ

            ! BUG: not always 0 when it should be
            kd = abs(coord_val - pt%freqI())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then

#ifdef BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(1) kd < machine precision. Assuming kd == 0.'
#endif
            else

               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd / tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd / tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = kd * tan(rot)
               endif
            endif

            fj = pt%freqJ() + cd
            Pe = MPoint(coord_val, fj)


         elseif (known_coord == 'j') then ! we search DI

            kd = abs(coord_val - pt%freqJ())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then

#ifdef BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(2) kd < machine precision. Assuming kd == 0.'
#endif
            else

               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd * tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = kd / tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = - kd / tan(rot)
               endif
            endif

            fi = pt%freqI() + cd
            Pe = MPoint(fi, coord_val)

         endif ! known_coord

         this%base_J_ = getPointsDistance(pt, Pe)



      elseif (base_dir == 'j') then ! we pass A pt, we know base along J-dir


         if (known_coord == 'i') then

            kd = abs(coord_val - pt%freqI())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then

#ifdef BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(3) kd < machine precision. Assuming kd == 0.'
#endif
            else

               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd * tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = - kd / tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = kd / tan(rot)
               endif
            endif

            fj = pt%freqJ() + cd
            Pe = MPoint(coord_val, fj)


         elseif (known_coord == 'j') then


            kd = abs(coord_val - pt%freqJ())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then

#ifdef BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(4) kd < machine precision. Assuming kd == 0.'
#endif
            else

               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd / tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd / tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = kd * tan(rot)
               endif
            endif

            fi = pt%freqI() + cd
            Pe = MPoint(fi, coord_val)

         endif ! known_coord


         this%base_I_ = getPointsDistance(pt, Pe)


      endif ! base_dir

      ! Now once defined save End point
      this%Ept_ = Pe
   end subroutine getOtherBase





   ! !> Gets actualised frequency deltas along two main
   ! !> sides directions (I, J), actualised based on *this
   ! !> zone rotation w.r.t. GRS.
   ! module subroutine getIJfsteps(this, dfIi, dfIj, dfJi, dfJj)
   !    class(MRectZone_t), intent(in) :: this
   !    real(bsa_real_t), intent(out)  :: dfIi, dfIj, dfJi, dfJj

   !    real(bsa_real_t) :: c, s, ang

   !    if (this%rot_ < CST_PId2) then ! FIRST quadrant

   !       c = cos(this%rot_)
   !       s = sin(this%rot_)

   !       dfIi =   this%deltaf_I_ * c
   !       dfIj = - this%deltaf_I_ * s

   !       dfJi = this%deltaf_J_ * s
   !       dfJj = this%deltaf_J_ * c

   !    elseif (this%rot_ < CST_PIGREC) then ! SECOND quadrant

   !       ang = this%rot_ - CST_PId2
   !       c   = cos(ang)
   !       s   = sin(ang)

   !       dfIi = - this%deltaf_I_ * c   ! - this%deltaf_I_ * s
   !       dfIj = - this%deltaf_I_ * s   ! - this%deltaf_I_ * c

   !       dfJi = - this%deltaf_J_ * s   !   this%deltaf_J_ * c
   !       dfJj =   this%deltaf_J_ * c   ! - this%deltaf_J_ * s

   !    elseif (this%rot_ < CST_PIt3d2) then ! THIRD quadrant

   !       ang = this%rot_ - CST_PIGREC
   !       c   = cos(ang)
   !       s   = sin(ang)

   !       dfIi = - this%deltaf_I_ * c
   !       dfIj =   this%deltaf_I_ * s

   !       dfJi = - this%deltaf_J_ * s
   !       dfJj = - this%deltaf_J_ * c

   !    elseif (this%rot_ < CST_PIt2) then ! FOURTH quadrant

   !       ang = this%rot_ - CST_PIt3d2
   !       c   = cos(ang)
   !       s   = sin(ang)

   !       dfIi =   this%deltaf_I_ * c   ! this%deltaf_I_ * s
   !       dfIj =   this%deltaf_I_ * s   ! this%deltaf_I_ * c

   !       dfJi =   this%deltaf_J_ * s   ! - this%deltaf_J_ * c
   !       dfJj = - this%deltaf_J_ * c   !   this%deltaf_J_ * s

   !    endif
   ! end subroutine







   ! module function reconstructZoneBaseMesh(this) result(msh)
   !    class(MRectZone_t), intent(in) :: this
   !    !> BUG: might be 2-rank array instead of 3!
   !    real(bsa_real_t) :: msh(2, this%nj_, this%ni_)

   !    real(bsa_real_t) :: dfIi, dfIj, dfJi, dfJj
   !    real(bsa_real_t) :: base_fi, base_fj, fi, fj
   !    integer(int32)   :: i, j

   !    call this%getIJfsteps(dfIi, dfIj, dfJi, dfJj)

   !    fi = this%Ipt_%freqI()
   !    fj = this%Ipt_%freqJ()
   !    base_fi = fi
   !    base_fj = fj

   !    msh(:, 1, 1) = [fj, fi]

   !    ! internal lines
   !    do j = 2, this%nj_
   !       fi = fi + dfJi
   !       fj = fj + dfJj
   !       msh(:, j, 1) = [fj, fi]
   !    enddo ! pj_head


   !    ! internal columns
   !    do i = 2, this%ni_

   !       base_fi = base_fi + dfIi
   !       base_fj = base_fj + dfIj

   !       fi = base_fi
   !       fj = base_fj

   !       msh(:, 1, i) = [fj, fi]

   !       do j = 2, this%nj_
   !          fi = fi + dfJi
   !          fj = fj + dfJj
   !          msh(:, j, i) = [fj, fi]
   !       enddo ! pj_head
   !    enddo ! pi_head
   ! end function 







   !> Gets vertex point pt coordinates in Nth quadrant
   !> w.r.t. Center point, in zone's LRS.
   pure module function getNthQuadVtx(this, iquad) result(pt)
      class(MRectZone_t), intent(in) :: this
      integer(int32), intent(in)     :: iquad
      type(MPoint_t)   :: pt
      real(bsa_real_t) :: c, s, rot, di, dj

      if (iquad < 1 .or. iquad > 4) return

      if (this%rot_ < CST_PId2) then

         c = cos(this%rot_)
         s = sin(this%rot_)

         select case (iquad)
         case (1)
            di = this%base_I_ * c + this%base_J_ * s
            dj = this%base_J_ * c - this%base_I_ * s

         case (2)
            di =   this%base_I_ * c
            dj = - this%base_I_ * s

         case (3)
            pt = this%Ipt_
            return

         case (4)
            di = this%base_J_ * s
            dj = this%base_J_ * c

         end select


      elseif (this%rot_ < CST_PIGREC) then

         rot = this%rot_ - CST_PId2
         c   = cos(rot)
         s   = sin(rot)

         select case (iquad)
         case (1)
            di =   this%base_J_ * c
            dj = - this%base_J_ * s

         case (2)
            di =    this%base_J_ * c - this%base_I_ * s
            dj = - (this%base_J_ * s + this%base_I_ * c)

         case (3)
            di = - this%base_I_ * s
            dj = - this%base_I_ * c

         case (4)
            pt = this%Ipt_
            return

         end select

      elseif (this%rot_ < CST_PIt3d2) then

         rot = this%rot_ - CST_PIGREC
         c   = cos(rot)
         s   = sin(rot)

         select case (iquad)
         case (1)
            pt = this%Ipt_
            return

         case (2)
            di = - this%base_J_ * s
            dj = - this%base_J_ * c

         case (3)
            di = - (this%base_J_ * s + this%base_I_ * c)
            dj = - this%base_J_ * c  + this%base_I_ * s

         case (4)
            di = - this%base_I_ * c
            dj =   this%base_I_ * s
         end select

      elseif (this%rot_ < CST_PIt2) then

         rot = this%rot_ - CST_PIt3d2
         c   = cos(rot)
         s   = sin(rot)

         select case (iquad)
         case (1)
            di = this%base_I_ * s
            dj = this%base_I_ * c

         case (2)
            pt = this%Ipt_
            return

         case (3)
            di = - this%base_J_ * c
            dj =   this%base_J_ * s

         case (4)
            di = - this%base_J_ * c + this%base_I_ * s
            dj =   this%base_J_ * s + this%base_I_ * c
         end select

      endif

      pt = MPoint(this%Ipt_%freqI() + di, this%Ipt_%freqJ() + dj)
   end function






   module subroutine dumpRZ(this)
      !! Dumps a RECTANGULAR zone.
      !!
      !! NOTE: Each specific zone dumping method is called 
      !!       from the STATIC MSaver_t procedure dump().
      !!       This makes no need for specific saver to have
      !!       their own dump() implementation, since this
      !!       current imlementation is what we are looking for.
      class(MRectZone_t), intent(in) :: this

      write(io_units_bfmdump(1)) MZone_ID%RECTANGLE

      ! init pt
      write(io_units_bfmdump(1)) this%Ipt_%freqI(), this%Ipt_%freqJ()

      write(io_units_bfmdump(1)) this%rot_
      write(io_units_bfmdump(1)) this%base_I_, this%base_J_

      ! NOTE: maybe useless ?
      write(io_units_bfmdump(1)) this%ni_, this%nj_
   end subroutine dumpRZ





   module subroutine undumpRZ(this)
      !! Undumps a RECTANGULAR zone.
      !!
      !! NOTE: Each specific zone dumping method is called 
      !!       from the STATIC procedure UndumpZone() in MZone Module.
      class(MRectZone_t), intent(inout) :: this

      real(bsa_real_t)   :: rval1, rval2
      integer(bsa_int_t) :: ival1, ival2

      ! init point
      read(io_units_bfmdump(1)) rval1, rval2
      call this%Ipt_%setfreqs(rval1, rval2)

      ! rotation
      read(io_units_bfmdump(1)) rval1
      this%rot_ = rval1

      ! sides
      read(io_units_bfmdump(1)) rval1, rval2
      this%base_I_ = rval1
      this%base_J_ = rval2

      ! refinements
      read(io_units_bfmdump(1)) ival1, ival2
      call this%setRefinements(ival1, ival2, .true.)
   end subroutine undumpRZ




#if (defined(BSA_USE_POD_DATA_CACHING)) || (defined(_OPENMP))
# define __new_interp_proc__
#endif

   module subroutine interpolateRZ( this &
#ifndef BSA_USE_POD_DATA_CACHING
# define __bfm_undump__ bfm,
      & , bfm &
#else
# define __bfm_undump__
#endif
      & , pdata )
      !! Implementation of rect zone interpolation wrapper routine
      class(MRectZone_t), intent(inout) :: this
#ifndef BSA_USE_POD_DATA_CACHING
      real(bsa_real_t), intent(in)  :: bfm(:, :)
#endif
      class(*), pointer, intent(in) :: pdata

      ! NOTE: for the moment only supporting HTPC method
      call interpolateRZ_HTPC_v3(this, __bfm_undump__  pdata)
   end subroutine



#include "_shared_poly2d.fi"

end submodule
