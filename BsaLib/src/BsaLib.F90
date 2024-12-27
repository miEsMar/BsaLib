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
!! along with BsaLib. If not, see <https://www.gnu.org/licenses/>.
module BsaLib
!! author: Michele Esposito Marzino
!! version: 0.1.0

   use, intrinsic :: iso_fortran_env
   implicit none (type, external)
   public


#include "CONSTANTS/_CONSTANTS.F90"


   interface bsa_exportBRdecomp
      module procedure bsa_exportBR_nocompute_
   end interface


   interface bsa_exportSkewness
      module procedure bsa_exportSkewness_compute_
      module procedure bsa_exportSkewness_nocompute_
   end interface




! **************************************************************************************
!    INTERAFCE FOR  PUBLIC  PROCEDURES
! **************************************************************************************
   interface



   ! **************************************
   !    GENERAL
   ! **************************************

      module subroutine bsa_printBSAHeader()
         !! Prints  `B.S.A.`  header to `stdout`.
      end subroutine


      module subroutine bsa_enableGPU()
         !! Enables GPU code.
      end subroutine


      module subroutine bsa_readInputParamsFromBSAFile()
      end subroutine


      module subroutine bsa_setSpatialSymmetry(isym)
         !# <span style="white-space: pre-line">
         ! Sets spatial symmetry value. 
         ! Valid options are:
         ! <ul>
         !    <li> <code>BSA_SPATIAL_SYM_NONE (0)</code>:  <i>NONE</i>, full space computed  (DEFAULT)</li>
         !    <li> <code>BSA_SPATIAL_SYM_HALF (2)</code>:  <i>HALF</i> space computed.</li>
         !    <li> <code>BSA_SPATIAL_SYM_FOUR (4)</code>:  <i>1-FOURTH</i> space computed.</li>
         ! </ul>
         ! </span>

         integer(bsa_int_t), value :: isym
      end subroutine



      module subroutine bsa_setSpectraSymmetries(ispctrsym)
         !# <span style="white-space: pre-line">
         ! Controls whether Spectra and Bispectra (tensorial) symmetries are exploited, or not.
         ! In many cases, specially when the loading process is a real process, 
         ! spectal tensors of Spectra and Bispectra are symmetric with respect to their main diagonal, that is: 
         !  $$ B_{ijk}(\omega_1,\omega_2) \equiv 
         !       B_{kji}(\omega_1,\omega_2) \quad \forall \{k, j, i\} \in \mathrm{P}(3, 3), $$
         !  where \(\mathrm{P}(3, 3)\) represents the set of all permutations of the 3 indexes, taken in groups of 3.
         ! <ol>
         !    <li> <code>NO  (0)</code>:  no tensor-elements symmetry is used <i>NONE</i>, full space computed  (DEFAULT)</li>
         !    <li> <code>YES (1)</code>:  tensor-elements symetry is used.
         !       This means that approximately \(\frac{N}{2}\) elements of the 2D trensor of Spectra,
         !       and approximately \(\frac{N}{6}\) elements of the 3D tensor of Bispectra are effectively computed,
         !       where \(N\) is the characteristic dimension of the square tensors (both 2D and 3D).
         !       </li>
         ! </ol>
         !
         ! @note
         ! For the case of tensor of Bispectra, if spatial symmetry is set to 
         ! <code>BSA_SPATIAL_SYM_FOUR [[bsalib(module):bsa_setspatialsymmetry(interface)]]</code>, 
         ! tensor-elements symmetry is automatically disabled.
         ! <span>

         integer(bsa_int_t), value :: ispctrsym
      end subroutine




      module subroutine bsa_setPremeshScheme(itype)
         !# <span style="white-space: pre-line">
         ! Set Pre-mesh Scheme. 
         ! Valid options:
         ! <ul>
         !   <li> <code>BSA_PREMESH_TYPE_DIAG_CREST_NO (0)</code>:  No zones to cover Diagonal crests in 2-4 quadrants
         !        (DEFAULT)</li>
         !   <li> <code>BSA_PREMESH_TYPE_DIAG_CREST_YES (1)</code>:  Diagonal crests in 2-4 quadrants are explicitly meshed</li>
         ! </ul>
         !
         ! @warning
         ! <code>itype = 1</code>  is <b>highly NOT recommended</b>  as per the current implementation.
         ! @endwarning
         ! </span>

         integer(bsa_int_t), value :: itype
      end subroutine


      module subroutine bsa_setPremeshMode(imode)
         !# <span style="white-space: pre-line">
         ! Sets Pre-mesh mode. 
         ! Valid options:
         ! <ul>
         !   <li><code>BSA_PREMESH_MODE_BASE (0)</code></li>
         !   <li><code>BSA_PREMESH_MODE_ZONE_REFINED (1)</code>  (DEFAULT)</li>
         ! </ul>
         ! </span>

         integer(bsa_int_t), value :: imode
      end subroutine



      module subroutine bsa_enableVisualMode()
         !# <span style="white-space: pre-line">
         !  Enables <i>visual</i> mode.
         !  Visual mode allows the user to write bispectra of modal/nodal responses.
         !  When in visual mode, Pre-mesh phase is skipped, and data read from file (<code>dumpfile</code>).
         !  Only the Post-mesh phase is done (mostly <i>interpolation</i>), and data written to files.
         !  </span>
      end subroutine



      module subroutine bsa_generateBSAInputFiles(run)
         !# <span style="white-space: pre-line">
         !  Enables generation of <code>BSA</code> (built-in executable) compatible input files.
         !
         !  This is particularly useful if needing to execute <code>BsaLib</code> in environments 
         !  (e.g. <i>clusters</i>) where the hosting program/library does not have the (license)
         !  rights to reside.
         !  For this exact reason, <code>BsaLib</code> is shipped with its built-in executable, 
         !  that can be easily compiled and run in restricted environments.
         !
         !  This routine can be then called to let <code>BsaLib</code> generate the input files 
         !  read by <code>BSA</code>.
         !
         !  </span>

         logical, value :: run
         !# If <code>.false.</code>, exits after files generation. 
         !  Otherwise, keeps running <code>BsaLib</code> normally.
      end subroutine



      module subroutine bsa_setVisualIndexes(indexes, modal)
         !# <span style="white-space: pre-line">
         ! 
         !  </span>

         integer(bsa_int_t), intent(in) :: indexes(3)
            !# <span style="white-space: pre-line">
            ! Array of indexes.
            ! Specify indexing combination of the bispectra to be exported.
            ! Resulting index \(i\) value is computed as:
            ! <ul type="1">
            !   <li>modal: $$i = id_1 + id_2*\mathtt{M} + id_3*\mathtt{M}^2$$ 
            !       where \(id_i \in [1, \dots, \mathtt{M}]\), \(\mathtt{M}=\) n. of kept modes.
            !   </li>
            !   <li>nodal: $$i = (id_1 - 1)*\mathtt{NNDOFs} + id_2$$
            !       where \(id_1\) refers to the <i>node</i> index, 
            !       \(id_2\) to the nodal <i>degree-of-freedom</i> index,
            !       \(\mathtt{NNDOFs}\) being the total number of nodal degrees-of-freedom.
            !   </li>
            ! </ul>
            !
            ! @note
            ! In the modal case, if <code>ONLY_DIAG</code> is <code>.true.</code> 
            ! (i.e. only the elements on the main tensor diagonal are computed), 
            ! the index is computed as simply $$ i = id_1 $$
            ! @endnote
            !  </span>
         logical, value :: modal
            !# <span style="white-space: pre-line">
            ! If <code>.true.</code>, indicates to export <i>modal</i> response bispectra.
            ! Otherwise, <i>nodal</i> bispectra is exported.
            !
            ! @note
            ! In a future implementation, other kinds of response bispectra could be allowed to 
            ! be exported (internal efforts, support reactions, etc.)
            ! @endnote
            ! </span>
      end subroutine



      module subroutine bsa_enableOnlyPremesh()
         !! If called, exits after Pre-mesh phase. Skips Post-mesh.
      end subroutine




      module subroutine bsa_doValidateModalData(bool)
         !# <span style="white-space: pre-line">
         !  Enables (or disables, default option) validation of modal (structural) data.
         ! </span>

         logical, intent(in) :: bool
            !# <span style="white-space: pre-line">
            !  If <code>.false.</code> (DEFAULT), do not perform any modal data validation.
            !  If <code>.true.</code>, checks if some modes in the structural modal matrix \(\boldsymbol{\Phi}\)
            !  are not 1-normalised, and removes them from the <i>effectively kept</i> modes.
            !  This is to avoid having statistical moments of modal responses that are not computed 
            !  with the same normalisation across all vibration modes (e.g. vertical and torsional modes).
            !
            ! @bug
            ! As per the current implementation, all non 1-normalised modes are discarded.
            ! However, this might not be the desired behaviour, if for example, modes are 
            ! normalised with a criteria different than the 1-normalisation approach, 
            ! or if only torsional modes (which are usually normalised to rotations and not displacements) 
            ! are to be considered.
            ! A fix is certainly needed to remove this limitation.
            ! @endbug
            ! </span>
      end subroutine



      module subroutine bsa_doValidateZoneDeltas(bool)
         !# <span style="white-space: pre-line">
         ! If <code>.true.</code> is given, enables zone's deltas validation.
         ! </span>
         logical, intent(in) :: bool
      end subroutine



      module subroutine bsa_setPolicyIDValidationValues(id, i_bfm, j_bfm, i_brm, j_brm)
         !# <span style="white-space: pre-line">
         ! Allows to set custom delta validation values for a given built-in Policy ID.
         !
         ! @note
         ! This API call has effect if delta validation is set to ON via the 
         ! <code>[[bsalib(module):bsa_dovalidatezonedeltas(interface)]]</code> API call.
         ! endnote
         ! </span>

         integer(int32), value :: id
            !! Built-in Policy ID
         integer(int32), value :: i_bfm
            !! BFM i-direction validation factor
         integer(int32), value :: j_bfm
            !! BFM j-direction validation factor
         integer(int32), value :: i_brm
            !! BRM i-direction validation factor
         integer(int32), value :: j_brm
            !! BRM j-direction validation factor
      end subroutine



      module subroutine bsa_setBfmMLR(bool)
         !# <span style="white-space: pre-line">
         ! Set BFM MLR (Multi-Level-Refinement) <code>ON/OFF</code>.
         ! If <code>.true. (ON)</code>, exact BFM points are added in the Post-Meshing phase, before actual interpolation steps 
         ! If <code>.false. (OFF)</code>, no other points than those coming from the Pre-Meshing phase are used.
         ! </span>

         logical, intent(in) :: bool
      end subroutine




      module subroutine bsa_Init()
         !! Initialises `BsaLib` runtime internals
      end subroutine



      module subroutine bsa_forceBsaClsExecution(bool)
         !# <span style="white-space: pre-line">
         ! Controls the forced execution of a \(\mathtt{Classic}\) approach.
         ! </span>
         logical, intent(in) :: bool
      end subroutine


      module subroutine bsa_setMaxBkgPeakRestriction(bool)
         !# <span style="white-space: pre-line">
         ! The width of the background (i.e. <i>quasi-static</i>) peak is equal to:
         ! $$ W_{\mathrm{bkg}} = \frac{\overline{U}}{L} \quad \mathrm{[Hz]}$$
         ! where \(\overline{U}\) is the mean wind speed, \(L\) the turbulence length scale.
         ! In a 3D-spatial turbulence, there are in total 9 turbulence scales (lengths), 
         ! organised in a \(3 \times 3\) matrix of turbulence scales:
         ! $$ \mathbf{L} = 
         !    \begin{Bmatrix}
         !       \mathbf{l}_1 && \mathbf{l}_2 && \mathbf{l}_3
         !    \end{Bmatrix} = 
         !    \begin{Bmatrix}
         !       L_{ux} && L_{vx} && L_{wx}\\L_{uy} && L_{vy} && L_{wy}\\L_{uz} && L_{vz} && L_{wz}
         !    \end{Bmatrix} $$
         ! By default, $$ L = \mathrm{maxval}(\mathbf{L}). $$
         ! If a restriction is applied (i.e. <code>.true.</code> is passed as argument), 
         ! then the maximum value of the matrix is reduced to the maximum value along the first column:
         !  $$ L = \mathrm{maxval}(\mathbf{l}_1) $$
         ! </span>

         logical, intent(in) :: bool
      end subroutine



      module subroutine bsa_setPODTruncationThreshold(rval)
         !# <span style="white-space: pre-line">
         ! Sets the threshold limit \(r_{lim}\) of total energy to be kept in the truncation of POD modes, 
         ! when decomposing the Cross-Spectral-Density-Matrix of base turbulence \(\mathbf{S}(\omega)\).
         !
         ! @warning
         ! If invoking this API call, a call to <code>[[bsalib(module):bsa_setpodnofmodeskept(interface)]]</code> will be automatically 
         ! invalidated, since this API call has higher precedence (more accurate, specifies energy content).
         ! </span>

         real(bsa_real_t), value :: rval
            !! Energy threshold limit, in the range \(\mathtt{rval} \in [0, \dots, 1] \).
            !! If `0` is passed, all energy is kept (i.e. no truncation, as \(\mathtt{rval} = 1\)).
            !! If \(\mathtt{rval}\) is passed in the range \([0, \dots, 100]\%\), it is automatically rescaled
            !! to the absolute range \([0, \dots, 1]\).
      end subroutine



      module subroutine bsa_setPODNOfModesKept(nmodes)
         !# <span style="white-space: pre-line">
         ! Contrary to <code>[[bsalib(module):bsa_setpodtruncationthreshold(interface)]]</code>, this API call sets the desired number of 
         ! POD modes to be kept in the truncation procedure.
         !
         ! @warning
         ! To make this API call effective, please make sure to not call
         ! <code>[[bsalib(module):bsa_setpodtruncationthreshold(interface)]]</code>
         ! as well, since that call would automatically invalidate this one (higher precedence).
         ! </span>

         integer(bsa_int_t), value :: nmodes
            !! number of POD modes (\(\le \mathtt{NNL}\)), with \(\mathtt{NNL}\) 
            !! number of effectively loaded nodes.
      end subroutine



      module subroutine bsa_exportPODInfo()
         !# <span style="white-space: pre-line">
         ! Enables runtime exporting of POD truncation information.
         ! This might serve to extract information about how much energy 
         ! each POD mode captures, at the various frequencies.
         !
         ! </span>
      end subroutine





      module subroutine bsa_Run(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_msh, m3mr_msh, m3mf_cls, m3mr_cls)
         !# <span style="white-space: pre-line">
         ! <b>Main</b> `BsaLib` procedure. 
         ! This should be the <b>last</b> API call to make, after all settings has been done.
         ! At return, post-processing API calls can be done to export/elaborate results.
         !
         ! @note
         ! For any argument that is not desired (provided that internal settings are compatible), 
         ! a <code>null()</code> argument can be passed.
         ! </span>

         real(bsa_real_t), target, allocatable, dimension(:) :: m2mf_cls
            !! \(2^{\mathrm{nd}}\) order moments, modal forces, \(\mathtt{Classic}\) approach
         real(bsa_real_t), target, allocatable, dimension(:) :: m2mr_cls
            !! \(2^{\mathrm{nd}}\) order moments, modal responses, \(\mathtt{Classic}\) approach
         real(bsa_real_t), target, allocatable, dimension(:) :: m2o2mr_cls
            !! \(2^{\mathrm{nd}}\) order spectral moments, modal responses, \(\mathtt{Classic}\) approach
         real(bsa_real_t), target, allocatable, dimension(:) :: m3mf_msh
            !! \(3^{\mathrm{rd}}\) order moments, modal forces, \(\mathtt{Mesher}\) approach
         real(bsa_real_t), target, allocatable, dimension(:) :: m3mr_msh
            !! \(3^{\mathrm{rd}}\) order moments, modal responses, \(\mathtt{Mesher}\) approach
         real(bsa_real_t), target, allocatable, dimension(:) :: m3mf_cls
            !! \(3^{\mathrm{rd}}\) order moments, modal forces, \(\mathtt{Classic}\) approach
         real(bsa_real_t), target, allocatable, dimension(:) :: m3mr_cls
            !! \(3^{\mathrm{rd}}\) order moments, modal responses, \(\mathtt{Classic}\) approach
      end subroutine



      module subroutine bsa_Finalise()
         !! Once `BsaLib` has finished, call this API procedure to ensure the cleanup of internal memory.
      end subroutine



      logical pure module function bsa_isCleaned()
         !! Query if `BsaLib` has cleaned (<code>.true.</code>) its internal memory, or not (<code>.false.</code>)
      end function





   ! **************************************
   !    SETTINGS 
   ! **************************************

      elemental module function bsa_isFullComp() result(bool)
         !# <span style="white-space: pre-line">
         ! Query if a FULL computation is issued.
         ! By FULL computation it is meant that all elements of the 2D and 3D tensors of 
         ! Spectra \(\mathbf{S}\) and Bispectra \(\mathbf{B}\) are computed, i.e. both 
         ! diagonal and outer-diagonal elements.
         ! If this API call returns <code>.false.</code>, it means that only diagonal elements of 
         ! such tensors are computed (<code>ONLY_DIAG = .true.</code>)
         !
         ! @warning
         ! Issuing a non-FULL computation is indeed good in terms of 
         ! <ul>
         !   <li>speed: much less computations required</li>
         !   <li>memory: also, much less memory allocated (no need to allocate the whole tensors, 
         !       if only the diagonal elements are used)</li>
         ! </ul>
         ! On the other hand, non-FULL computations are <b>highly</b> discouraged 
         ! in terms of results accuracy.
         ! </span>
         logical :: bool
      end function


      module subroutine bsa_setAnalysisType(isuban)
         !# <span style="white-space: pre-line">
         ! Sets analysis type.
         ! Valid options:
         ! <ol>
         !   <li>\(\mathtt{Classic}\)</li>
         !   <li>\(\mathtt{Mesher}\)</li>
         !   <li>\(\mathtt{Both}\)</li>
         ! </ol>
         ! </span>
         integer(bsa_int_t), value :: isuban
      end subroutine


      module subroutine bsa_setClassicMode(i_mode)
         !# <span style="white-space: pre-line">
         ! Sets \(\mathtt{Classic}\) computation mode.
         ! Valid options:
         ! <ol>
         !  <li> <code> BSA_CLASSIC_MODE_VECTOR </code> 
         !     The internal <i>vectorised</i> implementation is used (DEFAULT).
         !     This is indeed the preferred option in terms of speed.
         !     There is however a limitation of this approach:
         !     since it requires a considerable amount of allocated memory, 
         !     and considered the limit of memory that is requirable to the runtime 
         !     (before going to unoptimised mechanisms, such swap partition)
         !     if a given limit (\(\approx 8\) Gb) is exceeded, the `BsaLib` runtime 
         !     automatically switches to a <i>scalar</i> implementation.
         !   </li>
         !  <li> <code> BSA_CLASSIC_MODE_SCALAR </code> 
         !     The internal <i>scalar</i> implementation is used.
         !     While this implementation is slower than its vectorised counterpart, 
         !     it is certainly the one to be used to limit the memory footprint.
         !     Also, for very big cases, the `BsaLib` runtime might automatically switch 
         !     to this implementation is too much memory allocation is required.
         !   </li>
         ! </ol>
         ! </span>

         integer(bsa_int_t), value :: i_mode
      end subroutine


      module subroutine bsa_setVersion(ivers)
         !! @warning
         !! Deprecated.

         integer(bsa_int_t), value :: ivers
      end subroutine



      module subroutine bsa_setScalingConv(iconv)
         !# <span style="white-space: pre-line">
         ! Set Power-Spectral-Density (PSD) integration convention.
         ! There are two fundamental integration conventions:
         ! <ol>
         !  <li>
         !    The <i>frequency</i> convention (<code>BSA_PSD_CONVENTION_FREQ</code>), 
         !    where the integral of the PSD from \(0\) to \(+\infty\) gives the variance
         !    $$ m_2 = \sigma^2 = \int_0^{+\infty} S(f) df. $$
         !     </li>
         !  <li>
         !    The <i>pulsation</i> convention (<code>BSA_PSD_CONVENTION_PULS</code>), 
         !    where the integral of the PSD from \(-\infty\) to \(+\infty\) gives the variance
         !    $$ m_2 = \sigma^2 = \int_{-\infty}^{+\infty} S(\omega) d\omega. $$
         !     </li>
         ! </ol>
         !
         ! @note
         ! By DEFAULT, <code>BsaLib</code> uses the convention on <i>pulsations</i> (<code>BSA_PSD_CONVENTION_FREQ</code>).
         ! </span>

         integer(bsa_int_t), value :: iconv
      end subroutine



      module subroutine bsa_setSpectraComputation(ipsd, ibisp)
         !! Gives control on which kind of spectral features are computed.

         integer(bsa_int_t), value :: ipsd
            !! If \(1\), activates computation of Spectra (\(2^{\mathrm{nd}}\) order statisitcs)
         integer(bsa_int_t), value :: ibisp
            !! If \(1\), activates computation of Bispectra (\(3^{\mathrm{rd}}\) order statisitcs)
      end subroutine



      module subroutine bsa_setSpectraExtension(ionlydiag)
         !# <span style="white-space: pre-line">
         ! Sets extension of statistical information.
         ! There are two possible scenarios:
         ! <ul>
         !   <li><code>FULL</code>
         !     All elements of the 2D and 3D tensors of Spectra and Bispectra are computed.
         !     This is indeed the most accurate analysis.
         !      </li>
         !   <li><code>ONLY_DIAG</code>
         !     Only elements of the main diagonal are computed.
         !     While this approach is faster, and consumes much less memory, it is not 
         !     recommended for results accuracy issues.
         !      </li>
         ! </ul>
         !
         ! @note
         ! By DEFAULT, <code>BsaLib</code> performs a <code>FULL</code> analysis.
         ! </span>

         integer(bsa_int_t), value :: ionlydiag
      end subroutine



      module subroutine bsa_setTestMode(itest)
         !# <span style="white-space: pre-line">
         ! Enables testing mode.
         !
         ! @warning
         ! If testing mode is enabled, some important features and checks are turned OFF.
         ! Use with care.
         ! <span>

         integer(bsa_int_t), value :: itest
      end subroutine




      module subroutine bsa_setupClassic(nfreqs, df)
         !# <span style="white-space: pre-line">
         ! Main setup data for the \(\mathtt{Classic}\) approach.
         !
         ! @warning
         ! If testing mode is OFF, these values might be internally modified in order to 
         ! meet some minimum accuracy requirements.
         ! </span>

         integer(bsa_int_t), value :: nfreqs
            !# <span style="white-space: pre-line">
            !  number of discretisation points (i.e. n. of frequencies) 
            !
            ! @note 
            ! Refers to the interval \([0, f_{max}]\). 
            ! If the <i>pulsations</i> convention is used (<code>BSA_PSD_CONVENTION_PULS</code>),
            ! then this number will be automatically actualised to the interval \([-f_{max},f_{max}]\).
            ! <span>
         real(bsa_real_t), value :: df
            !# <span style="white-space: pre-line">
            ! Spacing between discretisation points (i.e. delta \(\Delta f\))
            ! </span>
      end subroutine



      module subroutine bsa_setupMesher(isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
         !# <span style="white-space: pre-line">
         ! Main setup data for the \(\mathtt{Mesher}\) approach.
         ! </span>

         integer(bsa_int_t), value :: isvd
            !! If `1`, enables use of POD techniques to decompose 
            !! Cross-Spectral-Density-Matrices of base wind turbulence.
         integer(bsa_int_t), value :: bkgrfmt
            !! Defines the base n. of refinement points \(N_{bkg}\) used to mesh the 
            !! background (quasi-static) zone, placed at the origin \((0, 0)\).
            !!
            !! @note
            !! Several other zones' discretisation depend on this value.
         real(bsa_real_t), value :: bkgaext
            !! Defines the factor by which the background zone is extended.
            !! This is done to allow the user extending this zone, avoiding cutting
            !! the zone's extensions where gradients are still important.
            !! Acts as a safety factor.
         real(bsa_real_t), value :: genpaext
            !! Defines the factor by which the any peak zone is extended.
            !! Same reasons as for the background zone.
         real(bsa_real_t), value :: maxaext
            !! Defines the factor by which the total covered area \(f_{max}\) is extended.
            !! In this case, this is done to avoid loosing information coming from the 
            !! secondary peaks, which are usually placed at extensions up to \(2\cdot f_i\),
            !! where \(f_i\) is the \(i-\)th modal frequency.
         integer(bsa_int_t), value :: ifcov
            !! @warning
            !! Deprecated
         integer(bsa_int_t), value :: idumpmod
            !! If `1`, includes modal data in the `dumpfile`.
            !!
            !! @note
            !! While it is not optimal in terms of disk usage,
            !! it is certainly recommended in order to keep needed information all in one place.
      end subroutine





   ! **************************************
   !    WIND
   ! **************************************

      module subroutine bsa_setWindDirections(dirs, ndirs)
         integer(bsa_int_t), intent(in) :: dirs(:)
         integer(bsa_int_t), value, optional :: ndirs
      end subroutine


      module subroutine bsa_setWindTurbComps(tc, ntc)
         !# <span style="white-space: pre-line">
         ! Specifies which turbulent components, among \(\{u, v, w\}\), should be considered
         ! in the definition of wind loads, and for which PSDs will be computed.
         ! </span>

         integer(bsa_int_t), intent(in) :: tc(:)
            !! Array of components. Each array element value should 
            !! be equal to one of \(\{1:u,\ 2:v,\ 3:w\}\).
         integer(bsa_int_t), value, optional :: ntc
            !! total n. of turbulent components to consider.
      end subroutine


      module subroutine bsa_setWindVertProf(iwprof)
         !# <span style="white-space: pre-line">
         ! Sets the wind vertical profile.
         ! Valid options:
         ! <ol>
         !   <li> <code>BSA_WIND_VERT_PROFILE_POWER</code>: uses a power law
         !         $$ \overline{U}(z) = \left({\frac{z}{z_{ref}}}\right)^{\alpha} $$
         !       </li>
         !   <li> <code>BSA_WIND_VERT_PROFILE_LOG</code>: uses a logarithmic law
         !         $$ \overline{U}(z) = {\frac{1}{k}}\sqrt{\frac{\tau_0}{\rho}}\ln{\frac{z}{z_0}} $$
         !       </li>
         ! </ol>
         ! </span>
         integer(bsa_int_t), value :: iwprof
      end subroutine



      module subroutine bsa_setPSDType(ipsd)
         !# <span style="white-space: pre-line">
         ! Sets base wind turbulence PSD (Power Spectral Density).
         ! Valid options:
         ! <ul>
         !   <li> <code>BSA_WIND_PSD_VONKARMAN</code>: VonKarman spectum:
         !       </li>
         !   <li> <code>BSA_WIND_PSD_KAIMAL</code>: Kaimal spectrum
         !       </li>
         !   <li> <code>BSA_WIND_PSD_DAVENPORT</code>: Davenport spectrum
         !       </li>
         ! </ul>
         ! </span>
         integer(bsa_int_t), value :: ipsd
      end subroutine



      module subroutine bsa_setWindAltDir(ivert)
         !# <span style="white-space: pre-line">
         ! Specifies which of the 3 Euclidean axes of the Wind Reference coordinate System
         ! is to be accounted as vertical axis.
         ! 
         ! @note
         ! This is an important detail, which affect how altitudes and so wind speeds 
         ! are determined.
         ! </span>

         integer(bsa_int_t), value :: ivert
            !! Axis index, with value among one of \(\{1:x,2:y,3:z\}\).
      end subroutine



      module subroutine bsa_setWindZoneLimits(lim, ilim)
         !# <span style="white-space: pre-line">
         ! Sets limits defining wind zones.
         ! In real cases, if the structural system spans long distances horizontally (bridges), or vertically (skyscrapers), 
         ! it might happen that the wind characteristics (wind speed, turbulence scales, etc.) change, so that 
         ! applying a unique wind all along would be physically incorrect.
         ! To do so, more than one <q>wind zone</q> can be defined, to provide differentiation of wind characteristics along the structure.
         !
         ! @warning
         ! Limits are usually considered to be defined along the X-axis of the 
         ! Global Reference coordinate System (GRS), assuming horizontal structural 
         ! systems, such as long-span bridges.
         ! In a later implementation, freedom to choose on which axis defining the wind 
         ! zones MUST be given, for flexibility and correctness reasons.
         ! </span>

#if  ((defined(__INTEL_COMPILER_BUILD_DATE)) && (__INTEL_COMPILER_BUILD_DATE >= 20221019))
         real(bsa_real_t), intent(in) :: lim(..)
         integer(bsa_int_t), intent(in), optional :: ilim(..)
#else
         real(bsa_real_t), intent(in), target     :: lim(:)
         integer(bsa_int_t), intent(in), optional :: ilim(:)   ! limits' index passed
#endif
      end subroutine



      module subroutine bsa_setAirDensity(aird)
         !! Specifies a custom value for air density \(\rho_{air}\). 
         !! Defaults to \(1.225 \ \mathrm{kg/m^3}\).

         real(bsa_real_t), value :: aird
      end subroutine



      module subroutine bsa_setGlobalRotMatW2G(rotW2G)
         !# <span style="white-space: pre-line">
         ! Sets global rotation matrix from WRS to GRS.
         ! Different are the rotation matrices from the local WRS (local to a wind zone), and the GRS.
         ! For those, see <code>[[bsalib(module):bsa_setwzrotmatw2g(interface)]]</code>.
         ! </span>

         real(bsa_real_t), intent(in) :: rotW2G(3, 3)
      end subroutine



      module subroutine bsa_setWZMeanWindVel(mat)
         !! Defines mean wind speeds \(\overline{U}\), for each wind zone.

         real(bsa_real_t), target, intent(in) :: mat(:)
      end subroutine



      module subroutine bsa_setWZRefAlt(Zref)
         !! Defines reference altitudes \(z_{ref}\), for each wind zone.

         real(bsa_real_t), target, intent(in) :: Zref(:)
      end subroutine



      module subroutine bsa_setTurbWindScales(L)
         !# <span style="white-space: pre-line">
         ! Defines wind turbulence scales \(\mathbf{L}\)
         ! $$ \begin{Bmatrix}
         !       L_{ux} && L_{vx} && L_{wx}\\L_{uy} && L_{vy} && L_{wy}\\L_{uz} && L_{vz} && L_{wz}
         !    \end{Bmatrix} $$
         ! for each wind zone.
         ! </span>

         real(bsa_real_t), target, intent(in) :: L(3, 3, *)
      end subroutine



      module subroutine bsa_setTurbWindSDT(sigma)
         !# <span style="white-space: pre-line">
         ! Defines wind turbulence standard deviation \(\boldsymbol{\sigma}\), for each wind zone.
         ! This is directly linked to the <i>turbulence intensity</i> \(I_{\{u,v,w\}}\), defined as:
         ! $$ I_{\{u,v,w\}} = \frac{\sigma_{\{u,v,w\}}}{\overline{U}} $$
         ! </span>

         real(bsa_real_t), target, intent(in) :: sigma(3, *)
      end subroutine



      module subroutine bsa_setWindCorrCoeffs(ccoeffs)
         !# <span style="white-space: pre-line">
         ! Defines the coefficients \(C_{\mu\xi}\)
         ! $$ \mathbf{C} = \begin{Bmatrix}
         !       C_{ux} && C_{vx} && C_{wx}\\C_{uy} && C_{vy} && C_{wy}\\C_{uz} && C_{vz} && C_{wz}
         !    \end{Bmatrix} $$
         ! used in the decreasing-exponential formulation of wind spatial coherence.
         ! </span>

         real(bsa_real_t), target, intent(in) :: ccoeffs(3, 3, *)
      end subroutine



      module subroutine bsa_setWindCorrExpnts(cexpn)
         !# <span style="white-space: pre-line">
         ! Defines the exponent coefficients \(p_{\mu\xi}\)
         ! $$ \mathbf{P} = \begin{Bmatrix}
         !       p_{ux} && p_{vx} && p_{wx}\\p_{uy} && p_{vy} && p_{wy}\\p_{uz} && p_{vz} && p_{wz}
         !    \end{Bmatrix} $$
         ! used in the decreasing-exponential formulation of wind spatial coherence.
         ! </span>

         real(bsa_real_t), target, intent(in) :: cexpn(3, 3, *)
      end subroutine



      module subroutine bsa_setIncidenceAngles(incang)
         !# <span style="white-space: pre-line">
         ! Defines wind mean incidence angle \(\overline{i}\), for each wind zone.
         ! </span>

         real(bsa_real_t), target, intent(in) :: incang(:)
      end subroutine



      module subroutine bsa_setWZRotMatW2G(rotW2G_L)
         !# <span style="white-space: pre-line">
         ! Defines the rotation matrix from the wind zone's local wind reference system (WRSl), to the global one (GRS).
         ! This to allow different principal wind flow directions for each wond zone.
         ! </span>

         real(bsa_real_t), target, intent(in) :: rotW2G_L(3, 3, *)
      end subroutine



      module subroutine bsa_setNodalVel(Unod)
         !# <span style="white-space: pre-line">
         ! Provides mean wind speeds \(\overline{U}_{i}\) at all structural nodes \(i\).
         !
         ! @note
         ! If this API call is made, then all the wind characteristics info (e.g. turbulence scales, turbulence intensities, etc.)
         ! are not needed, since they serve computing the final nodal wind speeds, which 
         ! this API call provides already.
         ! @endnote
         !
         ! @warning
         ! Currently, wind speeds at all <b>nodes</b> must be given, even if 
         ! not all of them are actually under the wind lading action.
         ! Info on which nodes are effectively loaded is needed, so that internally, 
         ! the distinction can be done.
         !
         ! @warning
         ! <b>Will be deprecated</b>
         ! </span>

         real(bsa_real_t), target, intent(in) :: Unod(:)
      end subroutine



      module subroutine bsa_setNodalWindZones(NodWZ)
         !# <span style="white-space: pre-line">
         ! Provides info on each node's wind zone. That is, in which wind zone each node is located.
         !
         ! @warning
         ! <b>Will be deprecated</b>
         ! </span>

         integer(bsa_int_t), target, intent(in) :: NodWZ(:)
      end subroutine


      module subroutine bsa_setNodalWindAltitudes(WnodAlt)
         !# <span style="white-space: pre-line">
         ! Sets altitude for each node.
         !
         ! @warning
         ! In a later implementation, this should be internally deferred based on 
         ! the selected vertical wind axis, and nodal coordinates.
         ! @endwarning
         !
         ! @warning
         ! <b>Will be deprecated</b>
         ! </span>

         real(bsa_real_t), target, intent(in) :: WnodAlt(:)
      end subroutine




      module subroutine bsa_setSpatialNodalCorr(nodCorr)
         !! @warning
         !! <b>Will be deprecated</b>

         real(bsa_real_t), target, intent(in) :: nodCorr(:, :)
      end subroutine



      module subroutine bsa_setWindFCoeffs(wfc)
         !! @warning
         !! <b>Will be deprecated</b>

         real(bsa_real_t), target, intent(in) :: wfc(:, :, :)
            !! Dimensions should be [nlibs_l, ndegw+3, nnodes_l]
      end subroutine


      module subroutine bsa_setPhitimesC(phiTc)
         !! @warning
         !! <b>Will be deprecated</b>

         real(bsa_real_t), target, intent(in) :: phiTc(:, :, :)
      end subroutine






   ! **************************************
   !    STRUCTURE
   ! **************************************

      module subroutine bsa_setNodalCoords(nn, coords)
         !! Provides nodal spatial coordinates.

         integer(bsa_int_t), value :: nn
            !! total number of nodes.
         real(bsa_real_t), target, contiguous :: coords(:, :)
            !! Dimensions `[3, nn]`, where `3` refers to the 3 spatial 
            !! directions \(\{x,y,z\}\).
      end subroutine


      module subroutine bsa_setTotalNumOfDOFs(ndofs)
         !# <span style="white-space: pre-line">
         ! Provides the total number of structural degrees-of-freedom \(\mathtt{NDOFs}\).
         ! $$ \mathtt{NDOFs} = \mathtt{NN} \cdot \mathtt{NNDOFs} $$
         ! where \(\mathtt{NN}\) is the total number of structural nodes, 
         ! \(\mathtt{NNDOFs}\) the number of degrees-of-freedom per node.
         ! </span>
         integer(bsa_int_t), value :: ndofs
      end subroutine


      module subroutine bsa_setNumOfNodalDOFs(nndofs)
         !# <span style="white-space: pre-line">
         ! Sets number of total degrees-of-freedom per node \(\mathtt{NNDOFs}\).
         ! </span>

         integer(bsa_int_t), value :: nndofs
      end subroutine



      module subroutine bsa_setTotalNOfNodes(nn)
         !# <span style="white-space: pre-line">
         ! Sets number of total structural nodes \(\mathtt{NN}\).
         ! </span>

         integer(bsa_int_t), value :: nn
      end subroutine



      module subroutine bsa_setLoadedNodalDOFs(libs_l, nlibs_l)
         !# <span style="white-space: pre-line">
         ! Sets number of nodal DOFs effectively loaded \(\mathtt{NNDOFsL}\), 
         ! where the ending \(\mathtt{L}\) signifies <q>loaded</q>.
         !
         ! @note
         ! This API call is optional. If not called, `BsaLib` automatically sets the number of 
         ! loaded nodal DOFs to the maximum, \(\mathtt{NNDOFs}\) (see <code>[[bsalib(module):bsa_setnumofnodaldofs(interface)]]</code>).
         ! </span>

         integer(bsa_int_t), intent(in), target, allocatable :: libs_l(:)
             !! Array of loaded nodal DOFs.
             !! @note
             !! Each array element value must be included in the range 
             !! \{[1, \mathtt{NNDOFs}]\}, where \(\mathtt{NNDOFs}\) is the total 
             !! n. of nodal degrees-of-freedom (see <code>[[bsalib(module):bsa_setnumofnodaldofs(interface)]]</code>).
         integer(bsa_int_t), value, optional :: nlibs_l
             !! N. of passed loaded nodal DOFs.
      end subroutine



      module subroutine bsa_setLoadedNodes(nodes_l, nn_l)
         !# <span style="white-space: pre-line">
         ! Sets number of structural nodes effectively loaded \(\mathtt{NNL}\), 
         ! where the ending \(\mathtt{L}\) signifies <q>loaded</q>.
         !
         ! @note
         ! This API call is optional. If not called, `BsaLib` automatically sets the number of 
         ! loaded nodes to match all nodes, \(\mathtt{NN}\) (see <code>[[bsalib(module):bsa_settotalnofnodes(interface)]]</code>).
         ! </span>

         integer(bsa_int_t), intent(in), target, allocatable :: nodes_l(:)
         integer(bsa_int_t), value, optional :: nn_l
      end subroutine



      module subroutine bsa_setModalInfo(ndofs, nm, Phi, natf)
         !# <span style="white-space: pre-line">
         ! Provides info about structural modal vibrations.
         ! Namely, the structural modal matrix \(\boldsymbol{\Phi}\) (eigenmatrix) <code>Phi</code>
         ! and the relative natural frequencies \(\mathbf{f}\) (eigenvalues) <code>natf</code>.
         !
         ! @note
         ! As a layer of internal verification and correctness checks, this API call forces to specify
         ! dimensions of the modal matrix (in order).
         ! Usually, they are the total number of structural degrees-of-freedom (see
         ! <code>[[bsalib(module):bsa_settotalnumofdofs(interface)]]</code>)
         ! and the number of vibration modes.
         ! If for any reason, a value mismatch is caught, <code>BsaLib</code> wil throw an internal error.
         ! </span>

         integer(bsa_int_t), value :: ndofs, nm
         real(bsa_real_t), intent(in), target :: Phi(ndofs, nm), natf(nm)
      end subroutine



      module subroutine bsa_setKeptModalShapes(modes)
         !# <span style="white-space: pre-line">
         ! This API call allows to specify a list of vibration modes to be actually accounted for in the computation.
         ! It might be tought to directly ask the user to provide <b>only</b> the modes to be used, 
         ! when using the <code>[[bsalib(module):bsa_setmodalinfo(interface)]]</code> API call.
         ! However, this might limit API flexibility, so that the user is given the chance to 
         ! give a list of indexes, and restrict the computation to those indexes only, 
         ! without requiring to bother to much on how the modal vibration information is given.
         ! </span>

         integer(bsa_int_t), intent(in) :: modes(:)
      end subroutine



      module subroutine bsa_setModalMatrices(nm, Mgen, Kgen, Cgen)
         !# <span style="white-space: pre-line">
         ! Provides structural modal matrices of mass, stiffness and damping, respectively.
         !
         ! @note
         ! At the current stage, <code>BsaLib</code> depends on some basic functionalities of 
         ! common Finite-Element (FE) software, such as determining the modal matrix of a structure.
         ! This choice has been made in order not to bind <code>BsaLib</code> to any specific 
         ! implementation, and let the final user the choice of the method used to determine these matrices.
         ! @endnote
         !
         ! @warning
         ! As per the current implementation, only symmetric modal mass and stiffness matrices are 
         ! considered, meaning that in fact, the modal matrix \(\boldsymbol{\Phi}\) given in 
         ! <code>[[bsalib(module):bsa_setmodalinfo(interface)]]</code> is obtained by solving the 
         ! following Eigenvalue problem:
         ! $$ \mathbf{K}\boldsymbol{\Phi} = \mathbf{M}\boldsymbol{\Phi}\boldsymbol{\Lambda} $$
         ! where \(\mathbf{K}\) and \(\mathbf{M}\) are the structural stiffness and mass matrices, 
         ! \(\boldsymbol{\Lambda}\) a diagonal matrix containing the coefficients of the characteristic 
         ! polynomial equation.
         ! In such cases, the resulting modal mass and stiffness matrices
         ! $$ \mathbf{K}^\star = \boldsymbol{\Phi}^T\mathbf{K}\boldsymbol{\Phi}; \quad \quad 
         !    \mathbf{M}^\star = \boldsymbol{\Phi}^T\mathbf{M}\boldsymbol{\Phi} $$
         ! will be diagonalised by \(\boldsymbol{\Phi}\).
         ! Therefore, only the diagonal elements are required by <code>BsaLib</code>.
         ! </span>

         integer(bsa_int_t), value :: nm
            !! n. of vibration modes
            !! @note Should match with the one given in <code>[[bsalib(module):bsa_setmodalinfo(interface)]]</code>.
         real(bsa_real_t), intent(in), target, dimension(nm) :: Mgen, Kgen
            !! 
         real(bsa_real_t), intent(in), target :: Cgen(nm, nm)
      end subroutine



      module subroutine bsa_setTotDamping(xsi)
         !# <span style="white-space: pre-line">
         ! Provides modal damping ratios \(\xi_i \ \ \forall i=1,\dots,\mathtt{NM}\).
         !
         ! @note
         ! This should include any possible source of damping, from structural damping 
         ! to aerodynamic damping in case of wind loading, etc.
         ! </span>

         real(bsa_real_t), target, intent(in) :: xsi(:)
      end subroutine



      pure module function bsa_getUsedModeShapes() result(modes)
         !# <span style="white-space: pre-line">
         ! Returns a copy of the array of effectively used vibration modes.
         ! </span>
         ! This is directly linked to <code>[[bsalib(module):bsa_dovalidatemodaldata(interface)]]</code>, 
         ! or <code>[[bsalib(module):bsa_setkeptmodalshapes(interface)]]</code> (even though this last 
         ! API call contains the same info, provided by the user).
         integer(bsa_int_t), allocatable :: modes(:)
      end function





   ! **************************************
   !    COMPUTE 
   ! **************************************

      module subroutine bsa_computeBRdecomp(m2mf, bkg, res)
         !# Returns the Background-Resonant decomposition, given the second order statistical moments 
         ! \(\mathbf{m}_2\) of the modal loads.
         !
         ! @note
         ! This is an API call to be used in a post-processing phase.

         real(bsa_real_t), intent(in)  :: m2mf(:)
         real(bsa_real_t), allocatable, intent(out) :: bkg(:), res(:)
      end subroutine


      module subroutine bsa_computePeakFactors(&
            m2, m2o2, obs_time, peak_g, sk, peak_ng_pos, peak_ng_neg)
         !# Computes Gaussian and Non-Gaussian Peak factors.
         ! Formulations are taken from CITE!.
         !
         ! @note
         ! This is an API call to be used in a post-processing phase.

         real(bsa_real_t), intent(in)  :: m2(:)
         real(bsa_real_t), intent(in)  :: m2o2(:)
         real(bsa_real_t), intent(in)  :: obs_time
         real(bsa_real_t), allocatable, intent(inout) :: peak_g(:)
         real(bsa_real_t), allocatable, intent(in)    :: sk(:)
         real(bsa_real_t), allocatable, intent(inout) :: peak_ng_pos(:)
         real(bsa_real_t), allocatable, intent(inout), optional :: peak_ng_neg(:)
      end subroutine






   ! **************************************
   !    I/O  -  EXPORTING
   ! **************************************

      module subroutine bsa_setOutputDirectory(dirname)
         !! Sets output directory.
         character(len=*), intent(in) :: dirname
      end subroutine



      module subroutine bsa_setOutFileName(fname)
         !! @warning Deprecated.
         character(len=*), intent(in) :: fname
      end subroutine



      module subroutine bsa_setOutUnit(iunit)
         !! Sets unit to be used for the output file. 
         !! Useful if needed to use an already opened unit from hosting program unit.

         integer(bsa_int_t), value :: iunit
      end subroutine



      module subroutine bsa_closeUnitsAtEnd()
         !! If called, tells `BsaLib` to close internal unit handles at exit.
         !! This is not a wanted behaviour if, for example, output unit is taken from 
         !! hosting program unit (see <code>[[bsalib(module):bsa_setoutunit(interface)]]</code>), 
         !! and might be used even after `BsaLib` returns from main (<code>[[bsalib(module):bsa_run(interface)]]</code>).
      end subroutine



      module subroutine bsa_setExportFileFormat(iform)
         !! Defines the desired format (e.g. `FORMATTED/UNFORMATTED`) to be used when exporting 
         !! statistical moments to files.

         integer(bsa_int_t), value :: iform
            !! Format flag, see `BSA_EXPORT_FORMAT_*`.
      end subroutine



      module subroutine bsa_setExportAppendMode(imode)
         !! @warning Deprecated.

         integer(bsa_int_t), value :: imode
      end subroutine



      module subroutine bsa_setExportDirectory(dirname)
         !! Defines a directory to place all the exporting of computed statistical moments, 
         !! when invoking any of the related API calls 
         !! (see for instance <code>[[bsalib(module):bsa_exportmomenttofile(interface)]]</code>).

         character(len = *), intent(in) :: dirname
      end subroutine



      module subroutine bsa_setExportInCurrDir()
         !! Differently from <code>[[bsalib(module):bsa_setexportdirectory(interface)]]</code>, 
         !! this API call telss `BsaLib` to place all exports in the current working directory.
      end subroutine



      module subroutine bsa_exportBR_nocompute_(fname, bkg, res, xsi)
         !! Exports the Background-Resonant decomposition, computing it internally 
         !! (see <code>[[bsalib(module):bsa_computebrdecomp(interface)]]</code>).

         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in)   :: bkg(:), res(:), xsi(:)
      end subroutine



      module subroutine bsa_exportMomentToFile(fname, vec)
         !! Exports given statistical moments to file. `fname` is the file name.
         !!
         !! @note
         !! Uses a `FORMATTED` output.

         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in)   :: vec(:)
      end subroutine



      module subroutine bsa_exportSkewness_nocompute_(fname, sk)
         !! Exports skewness to file, no internal computation.
         !!
         !! @note
         !! Uses a `FORMATTED` output.

         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in)   :: sk(:)
      end subroutine


      module subroutine bsa_exportSkewness_compute_(fname, dim, m2, m3)
         !! Exports skewness to file, with internal computation prior to write.
         !!
         !! @note
         !! Uses a `FORMATTED` output.

         character(len = *), intent(in)  :: fname
         integer(bsa_int_t), value       :: dim
            !# @note `dim` represents the principal dimension of the square tensors of 
            ! 2nd and 3rd statistical moments.
            ! So that:
            ! <ul>
            !   <li>if a non-FULL computation (see <code>[[bsalib(module):bsa_isfullcomp(interface)]]</code>) 
            !         is issued, then `dim` is equal to the array extensions. </li>
            !   <li> otherwise, `dim` refers to the characteristic square tensors dimension 
            !         (i.e. as if they were computed in a non-FULL case).
            !         Then, the correct indexes are computed internally. </li>
            ! </ul>
         real(bsa_real_t), intent(in) :: m2(:)
         real(bsa_real_t), intent(in) :: m3(:)
      end subroutine



      module subroutine bsa_exportPSDToFile(fname, psd, f)
         !! Exports series of PSDs (defined at frequencies `f`) to file.
         !!
         !! @note
         !! Uses a `FORMATTED` output.

         character(len = *), intent(in) :: fname
         real(bsa_real_t),   intent(in) :: psd(:, :)
         real(bsa_real_t),   intent(in), optional :: f(:)
      end subroutine



      module subroutine bsa_exportBispToFile(fname, bisp)
         !! Exports series of Bispectra to file.
         !!
         !! @note
         !! Uses a `FORMATTED` output.

         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in)   :: bisp(:, :, :)
      end subroutine



      module subroutine bsa_saveCoordinatesToFile(fname, coords)
         !! Exports nodal coordinates to file
         !!
         !! @note
         !! Uses a `FORMATTED` output.

         character(len = *), intent(in)  :: fname
         real(bsa_real_t), intent(in), target, optional :: coords(:, :)
      end subroutine


      module subroutine bsa_exportModalData()
         !! Exports modal data to file names `modal.txt`.
         !!
         !! @note
         !! Uses a `FORMATTED` output.
         !!
         !! @warning Will be deprecated.

      end subroutine


      module subroutine bsa_exportExtremeValuesToFile(fname, rvar)
         !! Exports extreme (or peak) values to file.
         !!
         !! @note
         !! Uses a `FORMATTED` output.

         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in) :: rvar(:)
      end subroutine



      module subroutine bsa_setBRMExportDefaultMode(imode)
         !! Specifies the desired mode of exporting Bispectra (assumes bispectra computation is turned ON).
         !!
         !! @note
         !! Real data is exported as `real(real32)` from the `iso_fortran_env` intrinsic module.
         !!
         !! @warning
         !! Contrary to the <q>visual mode</q> (see <code>[[bsalib(module):bsa_enablevisualmode(interface)]]</code>), 
         !! this API call enables exporting when actually doing the real computation.
         !! That is, in the case of the \(\mathtt{Mesher}\) approach, the Pre-mesh phase is not skipped.
         !! <b>All</b> bispectra are exported, so expect to have big files generated.
         !! @endwarning

         integer(bsa_int_t), value :: imode
      end subroutine


      module subroutine bsa_setBispExportCallback(fptr)
         !! Allows to specify a custom function callback for exporting bispectra.
         !! Applies to both cases of 
         !! <code>[[bsalib(module):bsa_enablevisualmode(interface)]]</code> and 
         !! <code>[[bsalib(module):bsa_setbrmexportdefaultmode(interface)]]</code>

         procedure(exportInterf_vect_), pointer, intent(in) :: fptr
            !! Function callback.
      end subroutine


   end interface




! **************************************************************************************
!    INTERFACE FOR PRIVATE PROCEDURES
! **************************************************************************************
   interface
      module subroutine mainClassic_(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_cls, m3mr_cls)
         real(bsa_real_t), allocatable, intent(inout), target :: &
            m2mf_cls(:), m2mr_cls(:), m2o2mr_cls(:), m3mf_cls(:), m3mr_cls(:)
      end subroutine


      module subroutine mainMesher_(m3mf_msh, m3mr_msh)
         real(bsa_real_t), target, allocatable :: m3mf_msh(:), m3mr_msh(:)
      end subroutine
   end interface


   private :: mainClassic_, mainMesher_


   ! BUG: these might be moved, and called via an internal function pointer.
   private :: bsa_exportSkewness_compute_, bsa_exportSkewness_nocompute_
   private :: bsa_exportBR_nocompute_


end module BsaLib
