/**
 * This file is part of BsaLib.
 * Copyright (C) 2024  Michele Esposito Marzino 
 *
 * BsaLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BsaLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BsaLib.  If not, see <https://www.gnu.org/licenses/>.
 * */
#pragma once

#include <stdint.h>
#include "_base.h"


#ifdef __cplusplus
extern "C" {
#endif

int  bsaclInit(unsigned n_threads);
int  bsaclInitDeviceMemory(void);
int  bsaclRun(unsigned i_thread, const size_t dim, __real * __RESTRICT_PTR __EXT_PTR_CONST h_res);
void bsaclAcquirePSDId(const uint32_t psdid);
void bsaclAcquireStructModMat(__real * __RESTRICT_PTR __EXT_PTR_CONST modmat, __real * __RESTRICT_PTR __EXT_PTR_CONST natf, const uint32_t ndofs, const uint32_t nmodes);
void bsaclAcquireModalMatrices(__real * __RESTRICT_PTR __EXT_PTR_CONST Mg, __real * __RESTRICT_PTR __EXT_PTR_CONST Cg, __real * __RESTRICT_PTR __EXT_PTR_CONST Kg);
void bsaclAcquireLoadedNodesList(int * __RESTRICT_PTR __EXT_PTR_CONST nodes_load, const uint32_t nnodes_l);
void bsaclAcquireUsedModesList(int * __RESTRICT_PTR __EXT_PTR_CONST modes, const uint32_t nmodes_eff);
void bsaclAcquireWindCoeffs(__real * __RESTRICT_PTR __EXT_PTR_CONST wfc, const uint32_t nnodes_l, const uint32_t nlibs, const uint32_t ndegw);
void bsaclAcquireTurbComponentsList(int * __RESTRICT_PTR __EXT_PTR_CONST tc, const uint32_t ntc);
void bsaclAcquirePhiTimesCMat(__real * __RESTRICT_PTR __EXT_PTR_CONST phi_T_c, const uint32_t nmodes_eff, const uint32_t nnodes_l, const uint32_t ndegw);
void bsaclAcquireNodalCorrelation(__real * __RESTRICT_PTR __EXT_PTR_CONST nod_corr, const uint32_t nnod_corr);
void bsaclAcquireWindNodalVelocities(__real * __RESTRICT_PTR __EXT_PTR_CONST nod_vel);
void bsaclAcquireWindNodalWindZones(int * __RESTRICT_PTR __EXT_PTR_CONST nod_wz);
void bsaclAcquireWindTurbScales(__real * __RESTRICT_PTR __EXT_PTR_CONST wt_scl, const uint32_t nwz);
void bsaclAcquireWindTurbStd(__real * __RESTRICT_PTR __EXT_PTR_CONST wt_std, const uint32_t nwz);
#ifdef BSACL_ENABLE_EVALFCT_PTR
void bsaclAcquireEvalFunc(evalFct_t fct);
void bsaclAcquireEvalFuncByStrings(char ** __RESTRICT_PTR __EXT_PTR_CONST strings);
void bsaclAcquireEvalFuncByFile(char * __RESTRICT_PTR __EXT_PTR_CONST filename);
#endif
void bsaclAcquireComputationFreqs(const unsigned i_thread,
   const uint32_t nfi, __real * __RESTRICT_PTR __EXT_PTR_CONST fi, const uint32_t nfj, __real * __RESTRICT_PTR __EXT_PTR_CONST fj);
void bsaclAcquireBaseWindTurbPSD(__real * __RESTRICT_PTR __EXT_PTR_CONST S_uvw);
void bsaclSetDeviceType(const uint32_t itype);
void bsaclVerifyMaxAllocCondition(size_t idim, uint32_t * __RESTRICT_PTR __EXT_PTR_CONST ican);
void bsaclAbort(const int ierr);
void bsaclFinalise(void);

#ifdef __cplusplus
}
#endif
