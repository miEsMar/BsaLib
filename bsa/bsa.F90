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
module data

   use BsaLib
   implicit none (type, external)
   public
   integer(int32), parameter :: IUN_EXTDATA = 22223
   logical :: l_formmode = .false.
   logical :: ext_data_read_ = .false.

   integer(int32) :: i_nnodes, i_nlibs, i_nnodesl, i_nlibsl
   integer(int32), target, allocatable :: nodesl(:), libsl(:)
   real(real64),   target, allocatable :: nod_cords(:, :)

   integer(int32) :: i_varu, i_su, i_vert, i_degw
   integer(int32) :: i_nzones
   real(real64)   :: r_aird
   real(real64)   :: r_rotW2G(3, 3)
   real(real64), target, allocatable :: r_Zref_z(:), r_UBref_z(:), r_alph_z(:), r_lims_z(:)
   real(real64), target, allocatable :: r_L_z(:, :, :), r_std_z(:, :), r_corrC_z(:, :, :)
   real(real64), target, allocatable :: r_corrEx_z(:, :, :), r_rotW2G_z(:, :, :), r_incang_z(:)

   integer(int32), target, allocatable :: i_wzNod(:)
   real(real64),   target, allocatable :: r_wAltNod(:), r_UBnod(:), r_corrNod(:, :)

   real(real64), target, allocatable :: r_wfc(:, :, :)

   integer(int32) :: i_nm, i_ndofs
   real(real64), target, allocatable :: r_natf(:), r_modm(:, :)
   real(real64), target, allocatable :: r_Mg(:), r_Kg(:), r_Cg(:, :)
   real(real64), target, allocatable :: r_xsist(:), r_xsiad(:)

   integer(int32) :: i_exprt_mode_ = BSA_EXPORT_MODE_REPLACE
   integer(int32) :: i_exprt_form_ = BSA_EXPORT_FORMAT_FORMATTED
   logical :: export_results_to_files_ = .true.

   logical :: is_visual_       = .false.
   logical :: is_visual_nodal_ = .false.

end module data



!=========================================================================================
!=========================================================================================
! MAIN PROGRAM
!=========================================================================================
!=========================================================================================
program bsa

   use data

   implicit none (type, external)

   integer(bsa_int_t) :: ierr

   ! local, used to retrieve single run BSA results
   real(bsa_real_t), allocatable :: bkg_(:),  res_(:)
   real(bsa_real_t), target, allocatable :: m2mf_(:), m2mr_(:), m2o2mr_(:)   ! NOTE: implicitly classic
   real(bsa_real_t), target, allocatable :: m3mf_cls_(:), m3mr_cls_(:)
   real(bsa_real_t), target, allocatable :: m3mf_msh_(:), m3mr_msh_(:)

   real(bsa_real_t), allocatable, dimension(:) :: m2_r_diag, m3_r_diag, sk_r_diag, m2o2_r_diag
   real(bsa_real_t), allocatable, dimension(:) :: m2_r_full, m3_r_full, sk_r_full, m2o2_r_full

   real(bsa_real_t), allocatable :: peak_pos_r_diag_g(:), peak_pos_r_diag_ng(:)
   real(bsa_real_t), allocatable :: peak_pos_r_full_g(:), peak_pos_r_full_ng(:)
   real(bsa_real_t), allocatable :: extr_pos_r_diag_g(:), extr_pos_r_diag_ng(:)
   real(bsa_real_t), allocatable :: extr_pos_r_full_g(:), extr_pos_r_full_ng(:)

   real(bsa_real_t), allocatable :: peak_neg_r_diag_ng(:)
   real(bsa_real_t), allocatable :: peak_neg_r_full_ng(:)
   real(bsa_real_t), allocatable :: extr_neg_r_diag_ng(:)
   real(bsa_real_t), allocatable :: extr_neg_r_full_ng(:)

   character(len = *), parameter   :: exp_prfx = 'export_'
   character(len = *), parameter   :: cls_sffx = 'cls', msh_sffx = 'msh'
   character(len = *), parameter   :: exp_fext = '.txt', udscr = '_'
   character(len = :), allocatable :: fname, cmb_sffx


!=========================================================================================
! MAIN BODY
!=========================================================================================
   call bsa_printBSAHeader()
   call printTool()
   call parseArgs()
   call readDataFiles()
   call setup()

   call bsa_setTotDamping(r_xsist)

   ierr = 0_bsa_int_t

   ! BUG: allow bsa_Run to accept already allocated entities
   !      (check for size match).
   call bsa_Run(m2mf_, m2mr_, m2o2mr_, m3mf_msh_, m3mr_msh_, m3mf_cls_, m3mr_cls_, ierr=ierr)
#ifdef BSALIB_SAFE_RETURN
   if (0_bsa_int_t /= ierr) goto 99
#endif

   ! POST-PROCESSING
   if (export_results_to_files_ .and. .not.is_visual_) then

      if (bsa_isFullComp()) then
         cmb_sffx = 'full' // BSA_FILE_NAME_CL_SUFFIX
      else
         cmb_sffx = 'diag' // BSA_FILE_NAME_CL_SUFFIX
      endif

      block
         integer(bsa_int_t), allocatable :: modes_(:)
         integer(bsa_int_t) :: nmodes_

         modes_  = bsa_getUsedModeShapes()
         nmodes_ = size(modes_)

         if (allocated(m2mf_)) then
            call bsa_computeBRdecomp(m2mf_, bkg_, res_)
            fname = exp_prfx // 'm2_BR_decomp.txt'
            call bsa_exportBRdecomp(fname, bkg_, res_, r_xsist(modes_))

            fname = exp_prfx // 'm2_mf_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m2mf_)

            if (allocated(m3mf_cls_)) &
               call bsa_exportSkewness(exp_prfx // 'sk_mf_' // cls_sffx // udscr // cmb_sffx // exp_fext, nmodes_, m2mf_, m3mf_cls_)
            if (allocated(m3mf_msh_)) &
               call bsa_exportSkewness(exp_prfx // 'sk_mf_' // msh_sffx // udscr // cmb_sffx // exp_fext, nmodes_, m2mf_, m3mf_msh_)
         endif


         if (allocated(m2mr_)) then

            fname = exp_prfx // 'm2_mr_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m2mr_)

            if (allocated(m3mr_msh_)) then
               fname = exp_prfx // 'sk_mr_' // msh_sffx // udscr // cmb_sffx // exp_fext
               call bsa_exportSkewness(fname, nmodes_, m2mr_, m3mr_msh_)
               call modalRecombination(r_modm(:, modes_), m2mr_, m3mr_msh_, m2o2mr_)
            endif

            if (allocated(m3mr_cls_)) then
               fname = exp_prfx // 'sk_mr_' // cls_sffx // udscr // cmb_sffx // exp_fext
               call bsa_exportSkewness(fname, nmodes_, m2mr_, m3mr_cls_)
               if (.not. allocated(m3mr_msh_)) &
                  call modalRecombination(r_modm(:, modes_), m2mr_, m3mr_cls_, m2o2mr_)
            endif


            ! TODO: introduce proper T variable
            call bsa_computePeakFactors(m2_r_diag, m2o2_r_diag, 600._bsa_real_t, &
               peak_pos_r_diag_g, sk_r_diag, peak_pos_r_diag_ng, peak_neg_r_diag_ng)
            call bsa_computePeakFactors(m2_r_full, m2o2_r_full, 600._bsa_real_t, &
               peak_pos_r_full_g, sk_r_full, peak_pos_r_full_ng, peak_neg_r_full_ng)


            if (allocated(sk_r_full)) then
               fname = exp_prfx // 'sk_r_full' // exp_fext
               call bsa_exportSkewness(fname, sk_r_full)
            endif


            ! TODO: maybe check if file already exists
            call bsa_saveCoordinatesToFile('coordinates.txt')

            if (bsa_isFullComp()) then
               if (allocated(peak_pos_r_diag_g))  then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'D' // '_g' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_pos_r_diag_g)
                  extr_pos_r_diag_g = peak_pos_r_diag_g  * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'D' // '_g' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_pos_r_diag_g)
               endif

               if (allocated(peak_pos_r_diag_ng)) then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_pos_r_diag_ng)
                  extr_pos_r_diag_ng = peak_pos_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_pos_r_diag_ng)
               endif
               if (allocated(peak_neg_r_diag_ng)) then
                  fname = exp_prfx // 'peak_neg_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_neg_r_diag_ng)
                  extr_neg_r_diag_ng = peak_neg_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_neg_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_neg_r_diag_ng)
               endif

               if (allocated(peak_pos_r_full_g))  then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'F' // '_g' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_pos_r_full_g)
                  extr_pos_r_full_g  = peak_pos_r_full_g  * sqrt(m2_r_full)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'F' // '_g' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_pos_r_full_g)
               endif

               if (allocated(peak_pos_r_full_ng)) then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_pos_r_full_ng)
                  extr_pos_r_full_ng = peak_pos_r_full_ng * sqrt(m2_r_full)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_pos_r_full_ng)
               endif
               if (allocated(peak_neg_r_full_ng)) then
                  fname = exp_prfx // 'peak_neg_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_neg_r_full_ng)
                  extr_neg_r_full_ng = peak_neg_r_full_ng * sqrt(m2_r_full)
                  fname = exp_prfx // 'extr_neg_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_neg_r_full_ng)
               endif

            else

               if (allocated(peak_pos_r_diag_g))  then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // '_g' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_pos_r_diag_g)
                  extr_pos_r_diag_g  = peak_pos_r_diag_g  * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // '_g' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_pos_r_diag_g)
               endif


               if (allocated(peak_pos_r_diag_ng)) then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_pos_r_diag_ng)
                  extr_pos_r_diag_ng = peak_pos_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_pos_r_diag_ng)
               endif
               if (allocated(peak_neg_r_diag_ng)) then
                  fname = exp_prfx // 'peak_neg_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, peak_neg_r_diag_ng)
                  extr_neg_r_diag_ng = peak_neg_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_neg_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportExtremeValuesToFile(fname, extr_neg_r_diag_ng)
               endif
            endif

         endif ! allocated m2mr
         if (allocated(modes_)) deallocate(modes_)

         if (allocated(m3mf_cls_)) then
            fname = exp_prfx // 'm3_mf_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mf_cls_)
         endif
         if (allocated(m3mr_cls_)) then
            fname = exp_prfx // 'm3_mr_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mr_cls_)
         endif
         if (allocated(m3mf_msh_)) then
            fname = exp_prfx // 'm3_mf_' // msh_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mf_msh_)
         endif
         if (allocated(m3mr_msh_)) then
            fname = exp_prfx // 'm3_mr_' // msh_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mr_msh_)
         endif
      end block
   endif

   if (allocated(m2mf_))      deallocate(m2mf_)
   if (allocated(m2mr_))      deallocate(m2mr_)
   if (allocated(m3mf_cls_))  deallocate(m3mf_cls_)
   if (allocated(m3mf_msh_))  deallocate(m3mf_msh_)
   if (allocated(m3mr_cls_))  deallocate(m3mr_cls_)
   if (allocated(m3mr_msh_))  deallocate(m3mr_msh_)

   if (allocated(peak_pos_r_diag_g))   deallocate(peak_pos_r_diag_g)
   if (allocated(peak_pos_r_diag_ng))  deallocate(peak_pos_r_diag_ng)
   if (allocated(peak_pos_r_full_g))   deallocate(peak_pos_r_full_g)
   if (allocated(peak_pos_r_full_ng))  deallocate(peak_pos_r_full_ng)
   if (allocated(extr_pos_r_diag_g))   deallocate(extr_pos_r_diag_g)
   if (allocated(extr_pos_r_diag_ng))  deallocate(extr_pos_r_diag_ng)
   if (allocated(extr_pos_r_full_g))   deallocate(extr_pos_r_full_g)
   if (allocated(extr_pos_r_full_ng))  deallocate(extr_pos_r_full_ng)


   99 call releaseMemory(ierr)

!=========================================================================================
!=========================================================================================
! END MAIN PROGRAM
!=========================================================================================
!=========================================================================================



contains ! utility procedures




   subroutine printTool()
      print *, '                          BSA  -  main program'
      print *, '   Bispectral Stochastic Analysis of MDOFs systems under stationary actions.'
      print '(/)'
   end subroutine


   subroutine usage()
      print *
      print *, ' Syntax:'
      print *, '    bsa.exe  [[options] file]'
      print *
      print *, ' NOTE: every command option can be preceded by "-", "--" or "/"'
      print *
      print *, ' Options:'
      print *, '   readmode  <mode>'
      print *, '        valid  <mode>  options:'
      print *, '           "formatted", "unformatted" (default)'
      print *
      print *, '   append-exports'
      print *, '        Appends to instead of overriding existing export files.'
      print *
      print *, '   export-binary'
      print *, '        Exports results using unformatted (binary) files.'
      print *
      print *, '   no-export'
      print *, '        Do not export results to files.'
      print *
      print *, '   out-file  <outfilename>'
      print *, '        If specified, uses  <outfilename>  as output file name.'
      print *
      print *, '   visual  [-m,-n]  idx1[-ix2-idx3]'
      print *, '        Enables visual mode (exporting Bispetrum to file).'
      print *, '        Flags:'
      print *, '           -m  modal'
      print *, '           -n  nodal'
      print *, '        If no flag is given, "-m" is assumed.'
      print *, '        idx{1,2,3} is an integer.'

      call releaseMemory(1)
   end subroutine




   function removeInputArgPlaceholder_(arg) result(arg_)
#if (_WIN32 & __INTEL_COMPILER)
!DIR$ ATTRIBUTES FORCEINLINE :: removeInputArgPlaceholder_
#endif
      character(len = *), intent(in)  :: arg
      character(len = :), allocatable :: arg_
      character(len = *), parameter   :: MINUS = '-'

      if (arg(1:1) == '/') then

         arg_ = arg(2 : len_trim(arg))

      elseif (arg(1:1) == MINUS) then

         if (arg(2:2) == MINUS) then
            arg_ = arg(3 : len_trim(arg))
         else
            arg_ = arg(2 : len_trim(arg))
         endif

      else ! take it as such

         arg_ = arg(1 : len_trim(arg))
      endif
   end function





   subroutine parseArgs()
      !! Parse input arguments, if any.
      integer :: argc

      argc = command_argument_count()
      if (argc == 0) return

      block
         integer :: currargc, istat
         character(len = 64) :: arg
         character(len = :), allocatable :: arg_

         currargc = 1

         ! main parsing loop
         do
            call get_command_argument(currargc, arg)
            arg_ = removeInputArgPlaceholder_(arg)

            select case (arg_)

               case ('readmode')
                  currargc = currargc + 1
                  call get_command_argument(currargc, arg, status=istat)
                  if (istat /= 0) call usage()

                  select case (arg(1:len_trim(arg)))
                     case ('formatted')
                        l_formmode = .true.
                     case ('unformatted')
                     case default
                        call usage()
                  end select

               case ('append-exports')
                  i_exprt_mode_ = BSA_EXPORT_MODE_APPEND

               case ('export-binary')
                  i_exprt_form_ = BSA_EXPORT_FORMAT_UNFORMATTED

               case ('no-export')
                  export_results_to_files_ = .false.


               case ('out-file')
                  currargc = currargc + 1
                  call get_command_argument(currargc, arg, status=istat)
                  if (istat /= 0) call usage()
                  call bsa_setOutFileName(arg(1 : len_trim(arg)))

               case ('visual')
                  currargc = currargc + 1
                  call getVisualCLIInfo_(currargc)

            end select

            if (currargc == argc) exit ! parsing finished
            currargc = currargc + 1    ! go to next input arg
         enddo
      end block

   end subroutine parseArgs




   subroutine getVisualCLIInfo_(iarg)
      integer, intent(inout) :: iarg
      character(len = 64)    :: arg
      integer :: i, istat, list(3)

      ! first flag must be either -m or -n (modal/nodal)
      call get_command_argument(iarg, arg, status=istat)
      if (istat /= 0) call releaseMemory(-12)
      i = 1
      do while ((arg(i:i)=='-' .or. arg(i:i)=='/'))
         i = i + 1
      enddo
      if (i > 1) then
         iarg = iarg + 1
         if (arg(i:) == 'n') is_visual_nodal_ = .true.
         ! flag was passed
         call get_command_argument(iarg, arg, status=istat)
         if (istat /= 0) call releaseMemory(-13)
      endif

      ! extract list
      list = 0
      call extractListFromString_(arg, list)
      if (list(1) == 0) goto 98
      if (is_visual_nodal_) then
         if (list(2) == 0) goto 98
      else
         if (list(2) == 0 .or. list(3) == 0) list = list(1)
      endif
      call bsa_setVisualIndexes(list, .not. is_visual_nodal_)
      goto 10

      98 continue
      print *, ERRMSG, "Error while parsing visual indexes."
      call releaseMemory(-34)
      return

      10 is_visual_ = .true.
      return
   end subroutine



   subroutine extractListFromString_(string, list)
      character(len = *), intent(in) :: string
      integer, intent(out) :: list(3)
      integer :: idx, p, i, l

      l = len_trim(string)
      i = 1
      p = 1
      do
         if (p > l .or. i == 4) exit
         idx = scan(string(p : l), '-')
         if (idx == 0) then
            idx = l+1
         else
            idx = idx + p - 1
         endif
         read(unit=string(p : idx-1), fmt=*) list(i)
         p = idx + 1
         i = i + 1
      enddo
   end subroutine



   subroutine setup()
      !! This routine collects all the 
      !! calls to the BsaLibrary interface routines,
      !! before launch of the main computing process.


      ! NOTE: everything before bsa_Init() does not affect global (internal) vars.

      call bsa_setExportInCurrDir()  ! NOTE: this to not break plotter functioning..
      call bsa_closeUnitsAtEnd()
      call bsa_setExportFileFormat(i_exprt_form_)
      call bsa_setExportAppendMode(i_exprt_mode_)

      if (is_visual_) then
         call bsa_enableVisualMode()
      else
         call bsa_setBRMExportDefaultMode(BSA_EXPORT_BRM_MODE_NONE)
      endif

      call bsa_doValidateZoneDeltas(.false.)
      call bsa_setBfmMLR(.false.)

      call bsa_forceBsaClsExecution(.true.)

      call bsa_setMaxBkgPeakRestriction(.true.)

      ! call bsa_setPODTruncationThreshold(100.d0)
      call bsa_setPODNOfModesKept(0)

      call bsa_Init()  ! This initialises all necessary instances.

      call bsa_readInputParamsFromBSAFile()
      call bsa_exportPODInfo()

      ! NODAL
      call bsa_setTotalNOfNodes(i_nnodes)
      call bsa_setNumOfNodalDOFs(i_nlibs)
      call bsa_setLoadedNodes(nodesl)
      call bsa_setLoadedNodalDOFs(libsl)
      call bsa_setNodalCoords(i_nnodes, nod_cords)

      ! WIND
      call bsa_setWindVertProf(i_varu)
      call bsa_setPSDType(i_su)
      call bsa_setWindAltDir(i_vert)
      call bsa_setWindZoneLimits(r_lims_z)
      call bsa_setAirDensity(real(r_aird, kind=bsa_real_t))
      call bsa_setGlobalRotMatW2G(real(r_rotW2G, kind=bsa_real_t))

      call bsa_setWZMeanWindVel(r_UBref_z)
      call bsa_setWZRefAlt(r_Zref_z)
      call bsa_setTurbWindScales(r_L_z)
      call bsa_setTurbWindSDT(r_std_z)
      call bsa_setWindCorrCoeffs(r_corrC_z)
      call bsa_setWindCorrExpnts(r_corrEx_z)
      call bsa_setWZRotMatW2G(r_rotW2G_z)
      call bsa_setIncidenceAngles(r_incang_z)
      call bsa_setNodalVel(r_UBnod)
      call bsa_setNodalWindZones(i_wzNod)
      call bsa_setNodalWindAltitudes(r_wAltNod)
      call bsa_setSpatialNodalCorr(r_corrNod)

      call bsa_setWindFCoeffs(r_wfc)

      ! MODAL
      call bsa_setModalInfo(i_ndofs, i_nm, r_modm, r_natf)
      call bsa_setModalMatrices(i_nm, r_Mg, r_Kg, r_Cg)

   end subroutine setup






   subroutine errAllocVarMsg_(varname, istat, emsg)
      character(len = *), intent(in) :: varname, emsg
      integer, intent(in) :: istat

      print '(/ 1x, 4a, i0)', &
         ERRMSG, 'Cannot allocate  "', varname, '".  Error code  ', istat
      print '(1x, 3a /)', &
         MSGCONT, 'Erorr message:  ', emsg(1 : len_trim(emsg))

      call releaseMemory(2)
   end subroutine


   subroutine errDeallocVarMsg_(varname, istat, emsg)
      character(len = *), intent(in) :: varname, emsg
      integer, intent(in) :: istat

      print '(/ 1x, 4a, i0)', &
         ERRMSG, 'Cannot de-allocate  "', varname, '".  Error code  ', istat
      print '(1x, 3a /)', &
         MSGCONT, 'Erorr message:  ', emsg(1 : len_trim(emsg))
      error stop
   end subroutine






   subroutine readDataFiles()

      call getExtData()

      if (.not. ext_data_read_) then
         print '(1x, 2a)', &
            ERRMSG, 'Error reading input data from files.'
         call releaseMemory(5)
      endif

#ifdef _BSA_DEBUG
      print '(1x, 2a)', &
         INFOMSG, ' BSA   data read correctly.'
#endif
   end subroutine



   subroutine openExtInputFile()
      integer :: istat
      character(len = :), allocatable :: form_, access_

      if (l_formmode) then
         form_   = 'formatted'
         access_ = 'sequential'
      else
         form_   = 'unformatted'
         access_ = 'stream'
      endif

      open(unit=IUN_EXTDATA    &
         , file=EXT_DATA_FNAME &
         , iostat=istat   &
         , form=form_     &
         , access=access_ &
         , action=IO_ACTION_READ)

      if (istat == 0) then
#ifdef _DEBUG
         print '(1x, 2a)', DBGMSG, 'Input file correctly opened.'
#endif
         rewind(IUN_EXTDATA)
         return
      endif

      print *
      print '(/1x, 2a, i0)', &
         ERRMSG, 'Error opening input file. Error code  ', istat
      error stop
   end subroutine



   subroutine getExtData()
      integer :: istat, itmp
      character(len = 132) :: emsg

      call openExtInputFile()

      if (l_formmode) then
         read(IUN_EXTDATA, *) i_nnodes
         read(IUN_EXTDATA, *) i_nlibs
         read(IUN_EXTDATA, *) i_nnodesl
         read(IUN_EXTDATA, *) i_nlibsl
      else
         read(IUN_EXTDATA) i_nnodes
         read(IUN_EXTDATA) i_nlibs
         read(IUN_EXTDATA) i_nnodesl
         read(IUN_EXTDATA) i_nlibsl
      endif

      allocate(nodesl(i_nnodesl), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('nodesl', istat, emsg)

      allocate(libsl(i_nlibsl), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('libsl', istat, emsg)

      allocate(nod_cords(3, i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('nod_cords', istat, emsg)

      if (l_formmode) then
         read(IUN_EXTDATA, *) nodesl
         read(IUN_EXTDATA, *) libsl
         read(IUN_EXTDATA, *) nod_cords
      else
         read(IUN_EXTDATA) nodesl
         read(IUN_EXTDATA) libsl
         read(IUN_EXTDATA) nod_cords
      endif


      if (l_formmode) then
         read(IUN_EXTDATA, *) i_varu
         read(IUN_EXTDATA, *) i_su
         read(IUN_EXTDATA, *) i_vert
         read(IUN_EXTDATA, *) i_degw
         read(IUN_EXTDATA, *) r_aird
         read(IUN_EXTDATA, *) r_rotW2G
         read(IUN_EXTDATA, *) i_nzones
      else
         read(IUN_EXTDATA) i_varu
         read(IUN_EXTDATA) i_su
         read(IUN_EXTDATA) i_vert
         read(IUN_EXTDATA) i_degw
         read(IUN_EXTDATA) r_aird
         read(IUN_EXTDATA) r_rotW2G
         read(IUN_EXTDATA) i_nzones
      endif

      allocate(r_Zref_z(i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Zref_z', istat, emsg)
      allocate(r_UBref_z(i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_UBref_z', istat, emsg)
      ! allocate(r_alph_z(i_nzones), stat=istat, errmsg=emsg)
      ! if (istat /= 0) call errAllocVarMsg_('r_alph_z', istat, emsg)
      allocate(r_L_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_L_z', istat, emsg)
      allocate(r_std_z(3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_std_z', istat, emsg)
      allocate(r_corrC_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_corrC_z', istat, emsg)
      allocate(r_corrEx_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_corrEx_z', istat, emsg)
      allocate(r_lims_z(i_nzones+1), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_lims_z', istat, emsg)
      allocate(r_rotW2G_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_rotW2G_z', istat, emsg)
      allocate(r_incang_z(i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_incang_z', istat, emsg)

      allocate(i_wzNod(i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('i_wzNod', istat, emsg)
      allocate(r_wAltNod(i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_wAltNod', istat, emsg)
      allocate(r_UBnod(i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_UBnod', istat, emsg)
      itmp = i_nnodes * i_nnodes
      itmp = (itmp + i_nnodes) / 2
      allocate(r_corrNod(itmp, 3), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_corrNod', istat, emsg)

      allocate(r_wfc(i_nlibsl, i_degw+3, i_nnodesl), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_wfc', istat, emsg)

      if (l_formmode) then
         read(IUN_EXTDATA, *) r_Zref_z
         read(IUN_EXTDATA, *) r_UBref_z
         ! read(IUN_EXTDATA, *) r_alph_z
         read(IUN_EXTDATA, *) r_L_z
         read(IUN_EXTDATA, *) r_std_z
         read(IUN_EXTDATA, *) r_corrC_z
         read(IUN_EXTDATA, *) r_corrEx_z
         read(IUN_EXTDATA, *) r_lims_z
         read(IUN_EXTDATA, *) r_rotW2G_z
         read(IUN_EXTDATA, *) r_incang_z
         read(IUN_EXTDATA, *) i_wzNod
         read(IUN_EXTDATA, *) r_wAltNod
         read(IUN_EXTDATA, *) r_UBnod
         read(IUN_EXTDATA, *) r_corrNod

         read(IUN_EXTDATA, *) r_wfc
      else
         read(IUN_EXTDATA) r_Zref_z
         read(IUN_EXTDATA) r_UBref_z
         ! read(IUN_EXTDATA) r_alph_z
         read(IUN_EXTDATA) r_L_z
         read(IUN_EXTDATA) r_std_z
         read(IUN_EXTDATA) r_corrC_z
         read(IUN_EXTDATA) r_corrEx_z
         read(IUN_EXTDATA) r_lims_z
         read(IUN_EXTDATA) r_rotW2G_z
         read(IUN_EXTDATA) r_incang_z
         read(IUN_EXTDATA) i_wzNod
         read(IUN_EXTDATA) r_wAltNod
         read(IUN_EXTDATA) r_UBnod
         read(IUN_EXTDATA) r_corrNod

         read(IUN_EXTDATA) r_wfc
      endif


      if (l_formmode) then
         read(IUN_EXTDATA, *) i_ndofs
         read(IUN_EXTDATA, *) i_nm
      else
         read(IUN_EXTDATA) i_ndofs
         read(IUN_EXTDATA) i_nm
      endif
      allocate(r_natf(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_natf', istat, emsg)
      allocate(r_modm(i_ndofs, i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_modm', istat, emsg)
      allocate(r_Mg(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Mg', istat, emsg)
      allocate(r_Kg(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Kg', istat, emsg)
      allocate(r_Cg(i_nm, i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Cg', istat, emsg)
      allocate(r_xsist(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_xsist', istat, emsg)
      ! allocate(r_xsiad(i_nm), stat=istat, errmsg=emsg)
      ! if (istat /= 0) call errAllocVarMsg_('r_xsiad', istat, emsg)

      if (l_formmode) then
         read(IUN_EXTDATA, *) r_natf
         read(IUN_EXTDATA, *) r_modm
         read(IUN_EXTDATA, *) r_Mg
         read(IUN_EXTDATA, *) r_Kg
         read(IUN_EXTDATA, *) r_Cg
         read(IUN_EXTDATA, *) r_xsist
         ! read(IUN_EXTDATA, *) r_xsiad
      else
         read(IUN_EXTDATA) r_natf
         read(IUN_EXTDATA) r_modm
         read(IUN_EXTDATA) r_Mg
         read(IUN_EXTDATA) r_Kg
         read(IUN_EXTDATA) r_Cg
         read(IUN_EXTDATA) r_xsist
         ! read(IUN_EXTDATA) r_xsiad
      endif


      ext_data_read_ = .true.
      close(IUN_EXTDATA)
#ifdef _BSA_DEBUG
      print '(1x, 2a)', &
         INFOMSG, 'Ext data read correctly.'
#endif
   end subroutine




   subroutine modalRecombination(r_Phi, m2mr, m3mr, m2o2mr)
      real(bsa_real_t), intent(in) :: r_Phi(:, :), m2mr(:), m3mr(:)
      real(bsa_real_t), intent(in), allocatable :: m2o2mr(:)
      integer :: ndofs, nmodes, nm3
      integer :: idof, imode
      integer :: imodeM1, posm2, posm3, itmp1
      logical :: is_diag
      real(bsa_real_t) :: phi1, r_tmp


      ndofs  = size(r_Phi, 1)
      nmodes = size(r_Phi, 2)
      nm3    = size(m3mr)

      is_diag = nm3 == nmodes

      allocate(m2_r_diag(ndofs))
      m2_r_diag = 0
      allocate(m3_r_diag(ndofs))
      m3_r_diag = 0
      allocate(m2o2_r_diag(ndofs))
      m2o2_r_diag = 0

      ! UNILATERAL (do it in any case)
      do imode = 1, nmodes

         imodeM1 = imode - 1

         if (is_diag) then
            posm2 = imode
            posm3 = imode
         else
            itmp1 = imodeM1 * nmodes
            posm2 = itmp1 + imode
            posm3 = itmp1*nmodes + posm2
         endif

         ! TODO: use do concurrent
         do idof = 1, ndofs

            phi1  = r_Phi(idof, imode)
            r_tmp = phi1 * phi1


            m2_r_diag(idof) = m2_r_diag(idof) + m2mr(posm2) * r_tmp
            if (allocated(m2o2mr)) m2o2_r_diag(idof) = m2o2_r_diag(idof) + m2o2mr(posm2) * r_tmp

            r_tmp = r_tmp * phi1
            m3_r_diag(idof) = m3_r_diag(idof) + m3mr(posm3) * r_tmp
         enddo ! idof
      enddo ! imode

      sk_r_diag = m3_r_diag / (m2_r_diag)**(CST_3d2)


      if (is_diag) return


      ! CCC (Complete Cubic Combination)
      allocate(m2_r_full(ndofs))
      m2_r_full = 0
      allocate(m3_r_full(ndofs))
      m3_r_full = 0
      allocate(m2o2_r_full(ndofs))
      m2o2_r_full = 0

      block
         integer :: jmode, kmode, jmodeM1, kmodeM1
         integer :: itmp12, itmp13, itmp23
         real(bsa_real_t), dimension(ndofs) :: phi2_, phi3_

         do kmode = 1, nmodes

            kmodeM1 = kmode - 1
            phi3_   = r_Phi(:, kmode)

            itmp12 = kmodeM1 * nmodes
            itmp13 = itmp12  * nmodes

            do jmode = 1, nmodes

               jmodeM1 = jmode - 1
               phi2_   = r_Phi(:, jmode)
               itmp23  = jmodeM1 * nmodes
               posm2   = itmp12 + jmode

               ! 2nd order
               ! TODO: use do concurrent
               do idof = 1, ndofs

                  m2_r_full(idof) = m2_r_full(idof) + &
                     m2mr(posm2) * phi2_(idof) * phi3_(idof)

                  if (allocated(m2o2mr)) &
                     m2o2_r_full(idof) = m2o2_r_full(idof) + &
                        m2o2mr(posm2) * phi2_(idof) * phi3_(idof)
               enddo ! idof


               do imode = 1, nmodes

                  posm3 = itmp13 + itmp23 + imode

                  ! TODO: use do concurrent
                  do idof = 1, ndofs

                     phi1 = r_Phi(idof, imode)

                     m3_r_full(idof) = m3_r_full(idof) + &
                        m3mr(posm3) * phi1 * phi2_(idof) * phi3_(idof)
                  enddo ! idof
               enddo ! imode
            enddo ! jmode
         enddo ! kmode

         sk_r_full = m3_r_full / (m2_r_full)**(CST_3d2)

      end block
   end subroutine





   subroutine releaseMemory(iexit)
      integer, intent(in) :: iexit
      integer :: istat
      character(len = 132) :: emsg
      logical :: lflag

      inquire(unit=IUN_EXTDATA, opened=lflag)
      if (lflag) close(IUN_EXTDATA)

      istat = 0

      if (allocated(nodesl)) deallocate(nodesl, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('nodesl', istat, emsg)

      if (allocated(libsl)) deallocate(libsl, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('libsl', istat, emsg)

      if (allocated(nod_cords)) deallocate(nod_cords, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('nod_cords', istat, emsg)


      if (allocated(r_Zref_z)) deallocate(r_Zref_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Zref_z', istat, emsg)
      if (allocated(r_UBref_z)) deallocate(r_UBref_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_UBref_z', istat, emsg)
      ! if (allocated(r_alph_z)) deallocate(r_alph_z, stat=istat, errmsg=emsg)
      ! if (istat /= 0) call errDeallocVarMsg_('r_alph_z', istat, emsg)
      if (allocated(r_L_z)) deallocate(r_L_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_L_z', istat, emsg)
      if (allocated(r_std_z)) deallocate(r_std_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_std_z', istat, emsg)
      if (allocated(r_corrC_z)) deallocate(r_corrC_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrC_z', istat, emsg)
      if (allocated(r_corrEx_z)) deallocate(r_corrEx_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrEx_z', istat, emsg)
      if (allocated(r_lims_z)) deallocate(r_lims_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_lims_z', istat, emsg)
      if (allocated(r_rotW2G_z)) deallocate(r_rotW2G_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_rotW2G_z', istat, emsg)
      if (allocated(r_incang_z)) deallocate(r_incang_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_incang_z', istat, emsg)

      if (allocated(i_wzNod)) deallocate(i_wzNod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('i_wzNod', istat, emsg)
      if (allocated(r_wAltNod)) deallocate(r_wAltNod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_wAltNod', istat, emsg)
      if (allocated(r_UBnod)) deallocate(r_UBnod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_UBnod', istat, emsg)
      if (allocated(r_corrNod)) deallocate(r_corrNod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrNod', istat, emsg)

      if (allocated(r_wfc)) deallocate(r_wfc, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrNod', istat, emsg)  


      if (allocated(r_natf)) deallocate(r_natf, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_natf', istat, emsg)
      if (allocated(r_modm)) deallocate(r_modm, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_modm', istat, emsg)
      if (allocated(r_Mg)) deallocate(r_Mg, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Mg', istat, emsg)
      if (allocated(r_Kg)) deallocate(r_Kg, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Kg', istat, emsg)
      if (allocated(r_Cg)) deallocate(r_Cg, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Cg', istat, emsg)
      if (allocated(r_xsist)) deallocate(r_xsist, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_xsist', istat, emsg)

      call bsa_Finalise()

      if (iexit == 0) then
         print '(/ 1x, 2a)', INFOMSG, 'BSA terminated correctly.'
      else
         print '(/ 1x, 2a, i0)', &
            ERRMSG, 'BSA terminated with error. Exit status code  ', iexit
         error stop
      endif
   end subroutine releaseMemory


end program
