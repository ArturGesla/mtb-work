/*
 * MATLAB Compiler: 8.4 (R2022a)
 * Date: Sun Feb 16 20:19:37 2025
 * Arguments: "-B""macro_default""-W""lib:magiclibrary""evalRhsAndJac.m"
 */

#define EXPORTING_magiclibrary 1
#include "magiclibrary.h"

static HMCRINSTANCE _mcr_inst = NULL; /* don't use nullptr; this may be either C or C++ */

#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultPrintHandler(const char *s)
{
    return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern C block */
#endif

#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultErrorHandler(const char *s)
{
    int written = 0;
    size_t len = 0;
    len = strlen(s);
    written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
    if (len > 0 && s[ len-1 ] != '\n')
        written += mclWrite(2 /* stderr */, "\n", sizeof(char));
    return written;
}

#ifdef __cplusplus
} /* End extern C block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_magiclibrary_C_API
#define LIB_magiclibrary_C_API /* No special import/export declaration */
#endif

LIB_magiclibrary_C_API 
bool MW_CALL_CONV magiclibraryInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
    if (_mcr_inst)
        return true;
    if (!mclmcrInitialize())
        return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream((void *)(magiclibraryInitializeWithHandlers));
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(&_mcr_inst,
                                                             error_handler, 
                                                             print_handler,
                                                             ctfStream);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
    return true;
}

LIB_magiclibrary_C_API 
bool MW_CALL_CONV magiclibraryInitialize(void)
{
    return magiclibraryInitializeWithHandlers(mclDefaultErrorHandler, 
                                            mclDefaultPrintHandler);
}

LIB_magiclibrary_C_API 
void MW_CALL_CONV magiclibraryTerminate(void)
{
    if (_mcr_inst)
        mclTerminateInstance(&_mcr_inst);
}

LIB_magiclibrary_C_API 
void MW_CALL_CONV magiclibraryPrintStackTrace(void) 
{
    char** stackTrace;
    int stackDepth = mclGetStackTrace(&stackTrace);
    int i;
    for(i=0; i<stackDepth; i++)
    {
        mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
        mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
    }
    mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_magiclibrary_C_API 
bool MW_CALL_CONV mlxEvalRhsAndJac(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    return mclFeval(_mcr_inst, "evalRhsAndJac", nlhs, plhs, nrhs, prhs);
}

LIB_magiclibrary_C_API 
bool MW_CALL_CONV mlfEvalRhsAndJac()
{
    return mclMlfFeval(_mcr_inst, "evalRhsAndJac", 0, 0, 0);
}

