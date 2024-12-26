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
module BsaLib_MZone

   use BsaLib_CONSTANTS, only: bsa_int_t, bsa_real_t, int32
   use BsaLib_Data,      only: bsa_Abort, test_no_bfm_mlr_
   use BsaLib_IO,        only: unit_dump_bfm_
   use BsaLib_MPolicy,   only: MPolicy_t, MPolicy_NULL, assignment(=), operator(==)
   implicit none (type, external)
   private
   public :: DefaultInitBaseZone, DumpZone, UndumpZone

   type, public :: MZoneEnum_t
      integer(int32) :: NULL      = 0
      integer(int32) :: RECTANGLE = 1
      integer(int32) :: TRIANGLE  = 2
      integer(int32) :: LINE      = 3
   end type MZoneEnum_t
   type(MZoneEnum_t), public, parameter :: MZone_ID = MZoneEnum_t()
   ! integer, public, parameter :: MZone_RECTANGLE = 1
   ! integer, public, parameter :: MZone_TRIANGLE  = 2
   ! integer, public, parameter :: MZone_LINEAR    = 3


   type, abstract, public :: MZone_t

      character(len = 64) :: name_   = ''
      type(MPolicy_t)     :: policy_

      !> Pointer to index of zone's interest modes
      integer(bsa_int_t), public :: id_im_

   contains
      procedure, pass :: zoneName
      procedure, pass :: setPolicy
      procedure, pass :: policy
      procedure, pass :: setInterestModeIndexPtr
      procedure, pass :: disableZonePolicyBfmMLR
      procedure(intf_MZoneIntFct_),    pass, deferred :: zoneTotNPts
      procedure(intf_MZoneGenIn_),     pass, deferred :: dump
      procedure(intf_MZoneGenInOut_),  pass, deferred :: undump
      procedure(intf_MZoneInterp_),    pass, deferred :: interpolate
   end type MZone_t



   abstract interface
      pure function intf_MZoneIntFct_(this) result(res)
         import :: MZone_t
         class(MZone_t), intent(in) :: this
         integer :: res
      end function

      subroutine intf_MZoneGenIn_(this)
         import :: MZone_t
         class(MZone_t), intent(in) :: this
      end subroutine

      subroutine intf_MZoneGenInOut_(this)
         import :: MZone_t
         class(MZone_t), intent(inout) :: this
      end subroutine

      subroutine intf_MZoneInterp_(this &
#ifndef BSA_USE_POD_DATA_CACHING
         & , bfm &
#endif
         &, pdata)
         import :: MZone_t, bsa_real_t
         class(MZone_t), intent(inout) :: this
#ifndef BSA_USE_POD_DATA_CACHING
         real(bsa_real_t), intent(in)  :: bfm(:, :)
#endif
         class(*), pointer, intent(in) :: pdata
      end subroutine
   end interface



contains



   subroutine DefaultInitBaseZone(this)
      class(MZone_t), intent(inout) :: this

      this%policy_ = MPolicy_NULL
   end subroutine


   subroutine zoneName(this, name_in)
#if (_WIN32 & __INTEL_COMPILER)
!DIR$ ATTRIBUTES FORCEINLINE :: zoneName
#endif
      class(MZone_t), intent(inout) :: this
      character(len=*), intent(in)  :: name_in

      this%name_ = name_in(1:len_trim(name_in))
   end subroutine




   subroutine setInterestModeIndexPtr(this, id)
      class(MZone_t), intent(inout)  :: this
      integer(bsa_int_t), intent(in) :: id

      this%id_im_ = id
   end subroutine




   subroutine setPolicy(this, var_in)
      class(MZone_t), intent(inout)  :: this
      class(*), intent(in) :: var_in
      select type (var_in)
         class is (MPolicy_t)
            this%policy_ = var_in
         type is (integer)
            this%policy_ = var_in
         class default
            call bsa_Abort('Unsupported type. Must be either "integer" or "MPolicy_t".')
      end select
   end subroutine

   function policy(this) result(pol_out)
      class(MZone_t), intent(inout)  :: this
      type(MPolicy_t) :: pol_out

      pol_out = this%policy_
   end function





#ifdef BSA_USE_POD_DATA_CACHING
# define __bfm_dump__
# define __decl__
#else
# define __bfm_dump__ ,data
# define __decl__ real(bsa_real_t), intent(in) :: data(:, :)
#endif
   subroutine DumpZone(z  __bfm_dump__ )
      class(MZone_t), intent(in)   :: z
      __decl__
#undef __bfm_dump__
#undef __decl__
      ! integer             :: tot

      ! dump specific zone data
      ! NOTE: keep this first since 
      !       we want to read as first parameter,
      !       the actual zone type identifier, so that
      !       we can directly specialise undumping in Mesh()
      !       routine.
      call z%dump()

      ! write common zone data
      write(unit_dump_bfm_) z%name_

      ! policy
      write(unit_dump_bfm_) z%policy_%getID()

      ! zone interest modes index ptr
      write(unit_dump_bfm_) z%id_im_


      ! Dump BFM data.
      !
      ! ! write how many bytes in total to be read
      ! ! afterwards. Then, dimBISP is automatically
      ! ! deferred knowing num of zone's meshing points
      ! tot = size(data)
      ! write(unit_dump_bfm_) tot
#ifndef BSA_USE_POD_DATA_CACHING
      write(unit_dump_bfm_) data ! NOTE: dimBISP first, then nj * ni
#endif
   end subroutine DumpZone



#ifdef BSA_USE_POD_DATA_CACHING
# define __bfm_dump__
# define __decl__
#else
# define __bfm_dump__  ,bfm_undump
# define __decl__ real(bsa_real_t), allocatable, intent(inout) :: bfm_undump(:, :)
#endif
   subroutine UndumpZone(z  __bfm_dump__ )
#ifndef BSA_USE_POD_DATA_CACHING
# ifdef _OPENMP
      use BsaLib_Data, only: dimM_bisp_
# endif
#endif
      class(MZone_t), intent(inout) :: z
      __decl__
      character(len = 64) :: name_hdr
      integer(int32)      :: zNp

#undef __bfm_dump__
#undef __decl__

      call z%undump()  ! read zone's specific data first

      ! read zone common data
      read(unit_dump_bfm_) name_hdr
      call z%zoneName(name_hdr(1:len_trim(name_hdr)))

      ! policy (ID)
      read(unit_dump_bfm_) zNp
      call z%setPolicy(zNp)
      if (test_no_bfm_mlr_) call z%disableZonePolicyBfmMLR()

      ! zone interest modes index ptr
      read(unit_dump_bfm_) z%id_im_


#ifndef BSA_USE_POD_DATA_CACHING
      ! once zone is undumped, get its N. of points
      ! NOTE: needed in order to correctly index into bfm_undump !
      zNp = z%zoneTotNPts()

# ifdef _OPENMP
      if (.not. allocated(bfm_undump)) then
         allocate(bfm_undump(dimM_bisp_, zNp))
      else
         ! reallocate if more space needed
         if (zNp > size(bfm_undump, 2)) then
            deallocate(bfm_undump)
            allocate(bfm_undump(dimM_bisp_, zNp))
         endif
      endif
# endif

      ! read actual BFM dumped data
      ! NOTE: in second dimension, nj leading over ni
      !       laydown.
      read(unit_dump_bfm_) bfm_undump(:, 1 : zNp)
#endif
   end subroutine UndumpZone




   elemental pure subroutine disableZonePolicyBfmMLR(z)
      class(MZone_t), intent(inout) :: z

      z%policy_%n_interp_bfm_lvs_ = 0
      z%policy_%bfm_pol_%i_fct_   = 1
      z%policy_%bfm_pol_%j_fct_   = 1
   end subroutine


end module
