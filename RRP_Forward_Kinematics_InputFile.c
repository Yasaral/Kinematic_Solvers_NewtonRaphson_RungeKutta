#define S_FUNCTION_NAME  RRP_Forward_Kinematics
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>
#include <stdio.h>

static void mdlInitializeSizes(SimStruct *S)
{
    /* See sfuntmpl_doc.c for more details on the macros below */

    ssSetNumSFcnParams(S, 0);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 2)) return;
	
   /*
    * Configure the input port 1 -- Independent Coordinates
    */
	ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortWidth(S, 0, 3);
	ssSetInputPortComplexSignal(S, 0, COMPLEX_NO);
    ssSetInputPortRequiredContiguous(S, 0, true);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
	ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
	
   /*
    * Configure the input port 2 -- Dependent Coordinates (Previous Step)
    */
	ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortWidth(S, 1, 6);
	ssSetInputPortComplexSignal(S, 1, COMPLEX_NO);
    ssSetInputPortRequiredContiguous(S, 1, true);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
	ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
	
   /*
    * Set the number of output ports.
    */
    if (!ssSetNumOutputPorts(S, 1)) return;

    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortWidth(S, 0, 6);
    ssSetOutputPortComplexSignal(S, 0, COMPLEX_NO);
	ssSetOutputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
	
   /*
    * Register reserved identifiers to avoid name conflict
    */
    if (ssGetSimMode(S) == SS_SIMMODE_RTWGEN) {
      /*
       * Register reserved identifier for OutputFcnSpec
       */
      ssRegMdlInfo(S, "RRP_Forward_Kinematics_wrapper", MDL_INFO_ID_RESERVED, 0, 0,
                   (void*) ssGetPath(S));
    }

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

	ssSetOptions(S,
				 SS_OPTION_USE_TLC_WITH_ACCELERATOR |
				 SS_OPTION_CAN_BE_CALLED_CONDITIONALLY |
                 SS_OPTION_EXCEPTION_FREE_CODE |
                 SS_OPTION_WORKS_WITH_CODE_REUSE |      
                 SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
	ssSetOffsetTime(S, 0, FIXED_IN_MINOR_STEP_OFFSET);

#if defined(ssSetModelReferenceSampleTimeDefaultInheritance)

	ssSetModelReferenceSampleTimeDefaultInheritance(S);

#endif

}



#undef MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)

  static void mdlInitializeConditions(SimStruct *S)
  {
  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#undef MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {
  }
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
	const real_T *posi  = (const real_T*) ssGetInputPortSignal(S,0);
	const real_T *posd_old  = (const real_T*) ssGetInputPortSignal(S,1);
	real_T       *posd = ssGetOutputPortSignal(S,0);
    
	RRP_Forward_Kinematics_wrapper( posi, posd_old, posd );
}



#undef MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif /* MDL_UPDATE */



#undef MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  static void mdlDerivatives(SimStruct *S)
  {
  }
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
}


/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
