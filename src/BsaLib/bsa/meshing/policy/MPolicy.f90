!! This file is part of BSA Library.
!! Copyright (C) 2023  Michele Esposito Marzino 
!!
!! BSA Library is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BSA Library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BSA Library.  If not, see <https://www.gnu.org/licenses/>.
module BsaLib_MPolicy
   
   use BsaLib_CONSTANTS, only: int32
   implicit none (type, external)
   private

   ! enum, bind(C)
   !    enumerator :: MPolicy_NULL  = 0
   !    enumerator :: MPolicy_DEF   = 1
   !    enumerator :: MPolicy_CONST = 2
   !    enumerator :: MPolicy_PRE_PEAK_1 = 3
   !    enumerator :: MPolicy_PRE_PEAK_2 = 4
   !    enumerator :: MPolicy_PEAK  = 5
   !    enumerator :: MPolicy_CREST = 6
   !    enumerator :: MPolicy_BASIN = 7
   !    enumerator :: MPolicy_PAD_ZONE_INTERN = 8
   !    enumerator :: MPolicy_PAD_ZONE_EXTERN = 9
   ! end enum
   integer(int32), public, parameter :: MPolicy_NULL  = 0
   integer(int32), public, parameter :: MPolicy_DEF   = 1
   integer(int32), public, parameter :: MPolicy_CONST = 2
   integer(int32), public, parameter :: MPolicy_PRE_PEAK_1 = 3
   integer(int32), public, parameter :: MPolicy_PRE_PEAK_2 = 4
   integer(int32), public, parameter :: MPolicy_PEAK  = 5
   integer(int32), public, parameter :: MPolicy_CREST = 6
   integer(int32), public, parameter :: MPolicy_BASIN = 7
   integer(int32), public, parameter :: MPolicy_PAD_ZONE_INTERN = 8
   integer(int32), public, parameter :: MPolicy_PAD_ZONE_EXTERN = 9

   
   type, public :: MPolicy_t

      integer(int32) :: delta_fI_fct_     = 0
      integer(int32) :: delta_fJ_fct_     = 0
      integer(int32) :: interp_bfm_I_fct_ = 0
      integer(int32) :: interp_bfm_J_fct_ = 0
      integer(int32) :: interp_I_fct_     = 0
      integer(int32) :: interp_J_fct_     = 0
      integer(int32) :: n_interp_bfm_lvs_ = 0
      
      integer(int32), private :: id_pol_ = 0
   contains

      procedure :: getID

      ! generic :: assignment(=) => MPolicy_fromPol_sub, MPolicy_fromID_sub
      ! procedure, private, pass :: MPolicy_fromID_sub, MPolicy_fromPol_sub
      ! generic :: operator(==) => MPolicy_isID_order1, MPolicy_isID_order2
      ! procedure, pass(pol), private :: MPolicy_isID_order1, MPolicy_isID_order2
   end type MPolicy_t


   interface MPolicy_t
      module procedure MPolicy_constructor_integer
      module procedure MPolicy_fromID
   end interface


   ! NOTE: COMMENT THESE INTERFACES OUT IF USING 
   !       TYPE-BOUND PROCEDURES.
   interface assignment(=)
      module procedure MPolicy_fromPol_sub
      module procedure MPolicy_fromID_sub
   end interface assignment(=)
   public :: assignment(=)

   interface operator(==)
      module procedure MPolicy_isID_order1
      module procedure MPolicy_isID_order2
   end interface operator(==)
   public :: operator(==)


contains


   elemental pure function MPolicy_constructor_integer(&
      dfi, dfj, interp_bfm_i, interp_bfm_j, interpi, interpj, nlevs, id) result(pol)
      integer, intent(in) :: dfi, dfj, interp_bfm_i, interp_bfm_j, interpi, interpj, nlevs, id
      type(MPolicy_t) :: pol

      pol%delta_fI_fct_     = int(dfi,          kind=int32)
      pol%delta_fJ_fct_     = int(dfj,          kind=int32)
      pol%interp_bfm_I_fct_ = int(interp_bfm_i, kind=int32)
      pol%interp_bfm_J_fct_ = int(interp_bfm_j, kind=int32)
      pol%interp_I_fct_     = int(interpi,      kind=int32)
      pol%interp_J_fct_     = int(interpj,      kind=int32)
      pol%n_interp_bfm_lvs_ = int(nlevs,        kind=int32)
      pol%id_pol_           = int(id,           kind=int32)
   end function



   elemental pure function MPolicy_fromID(mpol) result(pol)
      integer, intent(in) :: mpol
      type(MPolicy_t)     :: pol

      select case (mpol)
         case (MPolicy_NULL)
            pol = MPolicy_t(0, 0,   0, 0,   0, 0,    0, MPolicy_NULL)
         
         case (MPolicy_DEF)
            pol = MPolicy_t(2, 2,   2, 2,   2, 2,    1, MPolicy_DEF)
         
         case (MPolicy_CONST)
            pol = MPolicy_t(1, 1,   1, 1,   1, 1,    1, MPolicy_CONST)
         
         case (MPolicy_PRE_PEAK_1)
            pol = MPolicy_t(1, 4,   4, 4,   4, 2,    1, MPolicy_PRE_PEAK_1)
         
         case (MPolicy_PRE_PEAK_2)
            pol = MPolicy_t(1, 8,   2, 2,   4, 2,    1, MPolicy_PRE_PEAK_2)
         
         case (MPolicy_PEAK)
            pol = MPolicy_t(1, 1,   4, 4,   4, 4,    1, MPolicy_PEAK)
         
         case (MPolicy_CREST)
            pol = MPolicy_t(1, 4,   4, 4,   4, 2,    1, MPolicy_CREST)
         
         case (MPolicy_BASIN)
            pol = MPolicy_t(1, 4,   2, 2,   2, 2,    1, MPolicy_BASIN)
         
         case (MPolicy_PAD_ZONE_INTERN)
            pol = MPolicy_t(4, 4,   1, 1,   2, 2,    1, MPolicy_PAD_ZONE_INTERN)
         
         case (MPolicy_PAD_ZONE_EXTERN)
            pol = MPolicy_t(8, 8,   1, 1,   2, 2,    1, MPolicy_PAD_ZONE_EXTERN)
         
         case default
            pol = MPolicy_t(2, 2,   2, 2,   2, 2,    1, MPolicy_DEF)
      end select
   end function MPolicy_fromID



   elemental pure subroutine MPolicy_fromPol_sub(lhs, rhs_pol)
      type(MPolicy_t), intent(out)  :: lhs
      class(MPolicy_t), intent(in)  :: rhs_pol

      lhs%delta_fI_fct_     = rhs_pol%delta_fI_fct_
      lhs%delta_fJ_fct_     = rhs_pol%delta_fJ_fct_
      lhs%interp_bfm_I_fct_ = rhs_pol%interp_bfm_I_fct_
      lhs%interp_bfm_J_fct_ = rhs_pol%interp_bfm_J_fct_
      lhs%interp_I_fct_     = rhs_pol%interp_I_fct_
      lhs%interp_J_fct_     = rhs_pol%interp_J_fct_
      lhs%n_interp_bfm_lvs_ = rhs_pol%n_interp_bfm_lvs_
      lhs%id_pol_           = rhs_pol%id_pol_
   end subroutine MPolicy_fromPol_sub



   elemental pure subroutine MPolicy_fromID_sub(lhs, rhs_id)
      type(MPolicy_t), intent(out) :: lhs
      integer, intent(in) :: rhs_id

      lhs = MPolicy_t(rhs_id)
   end subroutine MPolicy_fromID_sub



   elemental pure function getID(this) result(id)
      class(MPolicy_t), intent(in) :: this
      integer(int32) :: id

      id = this%id_pol_
   end function getID





   elemental pure function MPolicy_isID_order1(pol, id) result(isID)
      class(MPolicy_t), intent(in) :: pol
      integer, intent(in)          :: id
      logical :: isID

      isID = pol%id_pol_ == id
   end function MPolicy_isID_order1


   elemental pure function MPolicy_isID_order2(id, pol) result(isID)
      integer, intent(in)          :: id
      class(MPolicy_t), intent(in) :: pol
      logical :: isID

      isID = pol%id_pol_ == id
   end function MPolicy_isID_order2


end module