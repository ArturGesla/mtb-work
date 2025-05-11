/*
 * MATLAB Compiler: 8.4 (R2022a)
 * Date: Sun Feb 16 20:19:37 2025
 * Arguments: "-B""macro_default""-W""lib:magiclibrary""evalRhsAndJac.m"
 */

#ifndef magiclibrary_h
#define magiclibrary_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_magiclibrary_C_API 
#define LIB_magiclibrary_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_magiclibrary_C_API 
bool MW_CALL_CONV magiclibraryInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_magiclibrary_C_API 
bool MW_CALL_CONV magiclibraryInitialize(void);

extern LIB_magiclibrary_C_API 
void MW_CALL_CONV magiclibraryTerminate(void);

extern LIB_magiclibrary_C_API 
void MW_CALL_CONV magiclibraryPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_magiclibrary_C_API 
bool MW_CALL_CONV mlxEvalRhsAndJac(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

/* C INTERFACE -- MLF WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_magiclibrary_C_API bool MW_CALL_CONV mlfEvalRhsAndJac();

#ifdef __cplusplus
}
#endif
/* C INTERFACE -- MLF WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#endif
