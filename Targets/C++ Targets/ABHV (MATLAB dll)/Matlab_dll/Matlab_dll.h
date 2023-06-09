//
// MATLAB Compiler: 8.0 (R2020a)
// Date: Mon Sep 12 18:00:26 2022
// Arguments:
// "-B""macro_default""-W""cpplib:Matlab_dll,legacy,version=1.0""-T""link:lib""-
// d""C:\Users\rahee\OneDrive\Desktop\Soved-Problems\Complete_Problem\iTarget
// Library\Matlab_dll\for_testing""-v""C:\Users\rahee\OneDrive\Desktop\Soved-Pro
// blems\Complete_Problem\Matlab_dll.m"
//

#ifndef Matlab_dll_h
#define Matlab_dll_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_Matlab_dll_C_API 
#define LIB_Matlab_dll_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_Matlab_dll_C_API 
bool MW_CALL_CONV Matlab_dllInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_Matlab_dll_C_API 
bool MW_CALL_CONV Matlab_dllInitialize(void);

extern LIB_Matlab_dll_C_API 
void MW_CALL_CONV Matlab_dllTerminate(void);

extern LIB_Matlab_dll_C_API 
void MW_CALL_CONV Matlab_dllPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_Matlab_dll_C_API 
bool MW_CALL_CONV mlxMatlab_dll(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_Matlab_dll
#define PUBLIC_Matlab_dll_CPP_API __declspec(dllexport)
#else
#define PUBLIC_Matlab_dll_CPP_API __declspec(dllimport)
#endif

#define LIB_Matlab_dll_CPP_API PUBLIC_Matlab_dll_CPP_API

#else

#if !defined(LIB_Matlab_dll_CPP_API)
#if defined(LIB_Matlab_dll_C_API)
#define LIB_Matlab_dll_CPP_API LIB_Matlab_dll_C_API
#else
#define LIB_Matlab_dll_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_Matlab_dll_CPP_API void MW_CALL_CONV Matlab_dll(int nargout, mwArray& status, mwArray& message, mwArray& solution_out, mwArray& start_states, mwArray& targets, const mwArray& FID, const mwArray& config_file, const mwArray& solution_in, const mwArray& current_time);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
