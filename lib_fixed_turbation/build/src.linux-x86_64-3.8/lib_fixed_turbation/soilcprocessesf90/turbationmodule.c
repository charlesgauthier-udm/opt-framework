/* File: turbationmodule.c
 * This file is auto-generated with f2py (version:1.24.4).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Fri Jan 19 13:40:23 2024
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif /* PY_SSIZE_T_CLEAN */

/* Unconditionally included */
#include <Python.h>
#include <numpy/npy_os.h>

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *turbation_error;
static PyObject *turbation_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
/* See fortranobject.h for definitions. The macros here are provided for BC. */
#define rank f2py_rank
#define shape f2py_shape
#define fshape f2py_shape
#define len f2py_len
#define flen f2py_flen
#define slen f2py_slen
#define size f2py_size

#define CHECKSCALAR(check,tcheck,name,show,var)\
    if (!(check)) {\
        char errstring[256];\
        sprintf(errstring, "%s: "show, "("tcheck") failed for "name, var);\
        PyErr_SetString(turbation_error,errstring);\
        /*goto capi_fail;*/\
    } else 
#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int
double_from_pyobj(double* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;
    if (PyFloat_Check(obj)) {
        *v = PyFloat_AsDouble(obj);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    tmp = PyNumber_Float(obj);
    if (tmp) {
        *v = PyFloat_AsDouble(tmp);
        Py_DECREF(tmp);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj)) {
        PyErr_Clear();
        tmp = PyObject_GetAttrString(obj,"real");
    }
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
        /*pass*/;
    }
    else if (PySequence_Check(obj)) {
        PyErr_Clear();
        tmp = PySequence_GetItem(obj, 0);
    }

    if (tmp) {
        if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = turbation_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}

static int
int_from_pyobj(int* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;

    if (PyLong_Check(obj)) {
        *v = Npy__PyLong_AsInt(obj);
        return !(*v == -1 && PyErr_Occurred());
    }

    tmp = PyNumber_Long(obj);
    if (tmp) {
        *v = Npy__PyLong_AsInt(tmp);
        Py_DECREF(tmp);
        return !(*v == -1 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj)) {
        PyErr_Clear();
        tmp = PyObject_GetAttrString(obj,"real");
    }
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
        /*pass*/;
    }
    else if (PySequence_Check(obj)) {
        PyErr_Clear();
        tmp = PySequence_GetItem(obj, 0);
    }

    if (tmp) {
        if (int_from_pyobj(v, tmp, errmess)) {
            Py_DECREF(tmp);
            return 1;
        }
        Py_DECREF(tmp);
    }

    {
        PyObject* err = PyErr_Occurred();
        if (err == NULL) {
            err = turbation_error;
        }
        PyErr_SetString(err, errmess);
    }
    return 0;
}

static int
float_from_pyobj(float* v, PyObject *obj, const char *errmess)
{
    double d=0.0;
    if (double_from_pyobj(&d,obj,errmess)) {
        *v = (float)d;
        return 1;
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************* turbation *********************************/
static char doc_f2py_rout_turbation_soilcprocesses_turbation[] = "\
litter_out,soilc_out = turbation(il1,il2,iccp1,cryodiffus,biodiffus,kterm,zbotw,isand,actlyr,spinfast,ipeatland,litter,soilc,[ilg,ignd,iccp2])\n\nWrapper for ``turbation``.\
\n\nParameters\n----------\n"
"il1 : input int\n"
"il2 : input int\n"
"iccp1 : input int\n"
"cryodiffus : input float\n"
"biodiffus : input float\n"
"kterm : input float\n"
"zbotw : input rank-2 array('f') with bounds (f2py_zbotw_d0,f2py_zbotw_d1)\n"
"isand : input rank-2 array('i') with bounds (f2py_isand_d0,f2py_isand_d1)\n"
"actlyr : input rank-1 array('f') with bounds (f2py_actlyr_d0)\n"
"spinfast : input int\n"
"ipeatland : input rank-1 array('i') with bounds (f2py_ipeatland_d0)\n"
"litter : input rank-3 array('f') with bounds (ilg,iccp2,ignd)\n"
"soilc : input rank-3 array('f') with bounds (ilg,iccp2,ignd)\n"
"\nOther Parameters\n----------------\n"
"ilg : input int, optional\n    Default: shape(litter, 0)\n"
"ignd : input int, optional\n    Default: shape(litter, 2)\n"
"iccp2 : input int, optional\n    Default: shape(litter, 1)\n"
"\nReturns\n-------\n"
"litter_out : rank-3 array('f') with bounds (ilg,iccp2,ignd)\n"
"soilc_out : rank-3 array('f') with bounds (ilg,iccp2,ignd)";
/* #declfortranroutine# */
static PyObject *f2py_rout_turbation_soilcprocesses_turbation(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,int*,int*,int*,int*,int*,float*,float*,float*,float*,int*,float*,int*,int*,float*,float*,float*,float*,int*,int*,int*,int*,int*,int*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    int il1 = 0;
    PyObject *il1_capi = Py_None;
    int il2 = 0;
    PyObject *il2_capi = Py_None;
    int ilg = 0;
    PyObject *ilg_capi = Py_None;
    int ignd = 0;
    PyObject *ignd_capi = Py_None;
    int iccp2 = 0;
    PyObject *iccp2_capi = Py_None;
    int iccp1 = 0;
    PyObject *iccp1_capi = Py_None;
    float cryodiffus = 0;
    PyObject *cryodiffus_capi = Py_None;
    float biodiffus = 0;
    PyObject *biodiffus_capi = Py_None;
    float kterm = 0;
    PyObject *kterm_capi = Py_None;
    float *zbotw = NULL;
    npy_intp zbotw_Dims[2] = {-1, -1};
    const int zbotw_Rank = 2;
    PyArrayObject *capi_zbotw_as_array = NULL;
    int capi_zbotw_intent = 0;
    PyObject *zbotw_capi = Py_None;
    int *isand = NULL;
    npy_intp isand_Dims[2] = {-1, -1};
    const int isand_Rank = 2;
    PyArrayObject *capi_isand_as_array = NULL;
    int capi_isand_intent = 0;
    PyObject *isand_capi = Py_None;
    float *actlyr = NULL;
    npy_intp actlyr_Dims[1] = {-1};
    const int actlyr_Rank = 1;
    PyArrayObject *capi_actlyr_as_array = NULL;
    int capi_actlyr_intent = 0;
    PyObject *actlyr_capi = Py_None;
    int spinfast = 0;
    PyObject *spinfast_capi = Py_None;
    int *ipeatland = NULL;
    npy_intp ipeatland_Dims[1] = {-1};
    const int ipeatland_Rank = 1;
    PyArrayObject *capi_ipeatland_as_array = NULL;
    int capi_ipeatland_intent = 0;
    PyObject *ipeatland_capi = Py_None;
    float *litter = NULL;
    npy_intp litter_Dims[3] = {-1, -1, -1};
    const int litter_Rank = 3;
    PyArrayObject *capi_litter_as_array = NULL;
    int capi_litter_intent = 0;
    PyObject *litter_capi = Py_None;
    float *soilc = NULL;
    npy_intp soilc_Dims[3] = {-1, -1, -1};
    const int soilc_Rank = 3;
    PyArrayObject *capi_soilc_as_array = NULL;
    int capi_soilc_intent = 0;
    PyObject *soilc_capi = Py_None;
    float *litter_out = NULL;
    npy_intp litter_out_Dims[3] = {-1, -1, -1};
    const int litter_out_Rank = 3;
    PyArrayObject *capi_litter_out_as_array = NULL;
    int capi_litter_out_intent = 0;
    float *soilc_out = NULL;
    npy_intp soilc_out_Dims[3] = {-1, -1, -1};
    const int soilc_out_Rank = 3;
    PyArrayObject *capi_soilc_out_as_array = NULL;
    int capi_soilc_out_intent = 0;
    int f2py_zbotw_d0 = 0;
    int f2py_zbotw_d1 = 0;
    int f2py_isand_d0 = 0;
    int f2py_isand_d1 = 0;
    int f2py_actlyr_d0 = 0;
    int f2py_ipeatland_d0 = 0;
    static char *capi_kwlist[] = {"il1","il2","iccp1","cryodiffus","biodiffus","kterm","zbotw","isand","actlyr","spinfast","ipeatland","litter","soilc","ilg","ignd","iccp2",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOOOOOOOOOOO|OOO:turbation.soilcprocesses.turbation",\
        capi_kwlist,&il1_capi,&il2_capi,&iccp1_capi,&cryodiffus_capi,&biodiffus_capi,&kterm_capi,&zbotw_capi,&isand_capi,&actlyr_capi,&spinfast_capi,&ipeatland_capi,&litter_capi,&soilc_capi,&ilg_capi,&ignd_capi,&iccp2_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable il1 */
        f2py_success = int_from_pyobj(&il1,il1_capi,"turbation.soilcprocesses.turbation() 1st argument (il1) can't be converted to int");
    if (f2py_success) {
    /* Processing variable il2 */
        f2py_success = int_from_pyobj(&il2,il2_capi,"turbation.soilcprocesses.turbation() 2nd argument (il2) can't be converted to int");
    if (f2py_success) {
    /* Processing variable iccp1 */
        f2py_success = int_from_pyobj(&iccp1,iccp1_capi,"turbation.soilcprocesses.turbation() 3rd argument (iccp1) can't be converted to int");
    if (f2py_success) {
    /* Processing variable cryodiffus */
        f2py_success = float_from_pyobj(&cryodiffus,cryodiffus_capi,"turbation.soilcprocesses.turbation() 4th argument (cryodiffus) can't be converted to float");
    if (f2py_success) {
    /* Processing variable biodiffus */
        f2py_success = float_from_pyobj(&biodiffus,biodiffus_capi,"turbation.soilcprocesses.turbation() 5th argument (biodiffus) can't be converted to float");
    if (f2py_success) {
    /* Processing variable kterm */
        f2py_success = float_from_pyobj(&kterm,kterm_capi,"turbation.soilcprocesses.turbation() 6th argument (kterm) can't be converted to float");
    if (f2py_success) {
    /* Processing variable spinfast */
        f2py_success = int_from_pyobj(&spinfast,spinfast_capi,"turbation.soilcprocesses.turbation() 10th argument (spinfast) can't be converted to int");
    if (f2py_success) {
    /* Processing variable ipeatland */
    ;
    capi_ipeatland_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the 11st argument `ipeatland`";
    capi_ipeatland_as_array = ndarray_from_pyobj(  NPY_INT,1,ipeatland_Dims,ipeatland_Rank,  capi_ipeatland_intent,ipeatland_capi,capi_errmess);
    if (capi_ipeatland_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        ipeatland = (int *)(PyArray_DATA(capi_ipeatland_as_array));

    /* Processing variable zbotw */
    ;
    capi_zbotw_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the 7th argument `zbotw`";
    capi_zbotw_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,zbotw_Dims,zbotw_Rank,  capi_zbotw_intent,zbotw_capi,capi_errmess);
    if (capi_zbotw_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        zbotw = (float *)(PyArray_DATA(capi_zbotw_as_array));

    /* Processing variable isand */
    ;
    capi_isand_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the 8th argument `isand`";
    capi_isand_as_array = ndarray_from_pyobj(  NPY_INT,1,isand_Dims,isand_Rank,  capi_isand_intent,isand_capi,capi_errmess);
    if (capi_isand_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        isand = (int *)(PyArray_DATA(capi_isand_as_array));

    /* Processing variable actlyr */
    ;
    capi_actlyr_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the 9th argument `actlyr`";
    capi_actlyr_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,actlyr_Dims,actlyr_Rank,  capi_actlyr_intent,actlyr_capi,capi_errmess);
    if (capi_actlyr_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        actlyr = (float *)(PyArray_DATA(capi_actlyr_as_array));

    /* Processing variable litter */
    ;
    capi_litter_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the 12nd argument `litter`";
    capi_litter_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,litter_Dims,litter_Rank,  capi_litter_intent,litter_capi,capi_errmess);
    if (capi_litter_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        litter = (float *)(PyArray_DATA(capi_litter_as_array));

    /* Processing variable ilg */
    if (ilg_capi == Py_None) ilg = shape(litter, 0); else
        f2py_success = int_from_pyobj(&ilg,ilg_capi,"turbation.soilcprocesses.turbation() 1st keyword (ilg) can't be converted to int");
    if (f2py_success) {
    CHECKSCALAR(shape(litter, 0) == ilg,"shape(litter, 0) == ilg","1st keyword ilg","turbation:ilg=%d",ilg) {
    /* Processing variable ignd */
    if (ignd_capi == Py_None) ignd = shape(litter, 2); else
        f2py_success = int_from_pyobj(&ignd,ignd_capi,"turbation.soilcprocesses.turbation() 2nd keyword (ignd) can't be converted to int");
    if (f2py_success) {
    CHECKSCALAR(shape(litter, 2) == ignd,"shape(litter, 2) == ignd","2nd keyword ignd","turbation:ignd=%d",ignd) {
    /* Processing variable iccp2 */
    if (iccp2_capi == Py_None) iccp2 = shape(litter, 1); else
        f2py_success = int_from_pyobj(&iccp2,iccp2_capi,"turbation.soilcprocesses.turbation() 3rd keyword (iccp2) can't be converted to int");
    if (f2py_success) {
    CHECKSCALAR(shape(litter, 1) == iccp2,"shape(litter, 1) == iccp2","3rd keyword iccp2","turbation:iccp2=%d",iccp2) {
    /* Processing variable soilc */
    soilc_Dims[0]=ilg,soilc_Dims[1]=iccp2,soilc_Dims[2]=ignd;
    capi_soilc_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the 13rd argument `soilc`";
    capi_soilc_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,soilc_Dims,soilc_Rank,  capi_soilc_intent,soilc_capi,capi_errmess);
    if (capi_soilc_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        soilc = (float *)(PyArray_DATA(capi_soilc_as_array));

    /* Processing variable litter_out */
    litter_out_Dims[0]=ilg,litter_out_Dims[1]=iccp2,litter_out_Dims[2]=ignd;
    capi_litter_out_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the hidden `litter_out`";
    capi_litter_out_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,litter_out_Dims,litter_out_Rank,  capi_litter_out_intent,Py_None,capi_errmess);
    if (capi_litter_out_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        litter_out = (float *)(PyArray_DATA(capi_litter_out_as_array));

    /* Processing variable soilc_out */
    soilc_out_Dims[0]=ilg,soilc_out_Dims[1]=iccp2,soilc_out_Dims[2]=ignd;
    capi_soilc_out_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.turbation: failed to create array from the hidden `soilc_out`";
    capi_soilc_out_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,soilc_out_Dims,soilc_out_Rank,  capi_soilc_out_intent,Py_None,capi_errmess);
    if (capi_soilc_out_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        soilc_out = (float *)(PyArray_DATA(capi_soilc_out_as_array));

    /* Processing variable f2py_zbotw_d0 */
    f2py_zbotw_d0 = shape(zbotw, 0);
    /* Processing variable f2py_zbotw_d1 */
    f2py_zbotw_d1 = shape(zbotw, 1);
    /* Processing variable f2py_isand_d0 */
    f2py_isand_d0 = shape(isand, 0);
    /* Processing variable f2py_isand_d1 */
    f2py_isand_d1 = shape(isand, 1);
    /* Processing variable f2py_actlyr_d0 */
    f2py_actlyr_d0 = shape(actlyr, 0);
    /* Processing variable f2py_ipeatland_d0 */
    f2py_ipeatland_d0 = shape(ipeatland, 0);
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
    (*f2py_func)(&il1,&il2,&ilg,&ignd,&iccp2,&iccp1,&cryodiffus,&biodiffus,&kterm,zbotw,isand,actlyr,&spinfast,ipeatland,litter,soilc,litter_out,soilc_out,&f2py_zbotw_d0,&f2py_zbotw_d1,&f2py_isand_d0,&f2py_isand_d1,&f2py_actlyr_d0,&f2py_ipeatland_d0);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("NN",capi_litter_out_as_array,capi_soilc_out_as_array);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    /* End of cleaning variable f2py_ipeatland_d0 */
    /* End of cleaning variable f2py_actlyr_d0 */
    /* End of cleaning variable f2py_isand_d1 */
    /* End of cleaning variable f2py_isand_d0 */
    /* End of cleaning variable f2py_zbotw_d1 */
    /* End of cleaning variable f2py_zbotw_d0 */
    }  /* if (capi_soilc_out_as_array == NULL) ... else of soilc_out */
    /* End of cleaning variable soilc_out */
    }  /* if (capi_litter_out_as_array == NULL) ... else of litter_out */
    /* End of cleaning variable litter_out */
    if((PyObject *)capi_soilc_as_array!=soilc_capi) {
        Py_XDECREF(capi_soilc_as_array); }
    }  /* if (capi_soilc_as_array == NULL) ... else of soilc */
    /* End of cleaning variable soilc */
    } /*CHECKSCALAR(shape(litter, 1) == iccp2)*/
    } /*if (f2py_success) of iccp2*/
    /* End of cleaning variable iccp2 */
    } /*CHECKSCALAR(shape(litter, 2) == ignd)*/
    } /*if (f2py_success) of ignd*/
    /* End of cleaning variable ignd */
    } /*CHECKSCALAR(shape(litter, 0) == ilg)*/
    } /*if (f2py_success) of ilg*/
    /* End of cleaning variable ilg */
    if((PyObject *)capi_litter_as_array!=litter_capi) {
        Py_XDECREF(capi_litter_as_array); }
    }  /* if (capi_litter_as_array == NULL) ... else of litter */
    /* End of cleaning variable litter */
    if((PyObject *)capi_actlyr_as_array!=actlyr_capi) {
        Py_XDECREF(capi_actlyr_as_array); }
    }  /* if (capi_actlyr_as_array == NULL) ... else of actlyr */
    /* End of cleaning variable actlyr */
    if((PyObject *)capi_isand_as_array!=isand_capi) {
        Py_XDECREF(capi_isand_as_array); }
    }  /* if (capi_isand_as_array == NULL) ... else of isand */
    /* End of cleaning variable isand */
    if((PyObject *)capi_zbotw_as_array!=zbotw_capi) {
        Py_XDECREF(capi_zbotw_as_array); }
    }  /* if (capi_zbotw_as_array == NULL) ... else of zbotw */
    /* End of cleaning variable zbotw */
    if((PyObject *)capi_ipeatland_as_array!=ipeatland_capi) {
        Py_XDECREF(capi_ipeatland_as_array); }
    }  /* if (capi_ipeatland_as_array == NULL) ... else of ipeatland */
    /* End of cleaning variable ipeatland */
    } /*if (f2py_success) of spinfast*/
    /* End of cleaning variable spinfast */
    } /*if (f2py_success) of kterm*/
    /* End of cleaning variable kterm */
    } /*if (f2py_success) of biodiffus*/
    /* End of cleaning variable biodiffus */
    } /*if (f2py_success) of cryodiffus*/
    /* End of cleaning variable cryodiffus */
    } /*if (f2py_success) of iccp1*/
    /* End of cleaning variable iccp1 */
    } /*if (f2py_success) of il2*/
    /* End of cleaning variable il2 */
    } /*if (f2py_success) of il1*/
    /* End of cleaning variable il1 */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/****************************** end of turbation ******************************/

/********************************** tridiag **********************************/
static char doc_f2py_rout_turbation_soilcprocesses_tridiag[] = "\
u = tridiag(a,b,c,r)\n\nWrapper for ``tridiag``.\
\n\nParameters\n----------\n"
"a : input rank-1 array('f') with bounds (f2py_a_d0)\n"
"b : input rank-1 array('f') with bounds (f2py_b_d0)\n"
"c : input rank-1 array('f') with bounds (f2py_c_d0)\n"
"r : input rank-1 array('f') with bounds (f2py_r_d0)\n"
"\nReturns\n-------\n"
"u : rank-1 array('f') with bounds (f2py_u_d0)";
/* #declfortranroutine# */
static PyObject *f2py_rout_turbation_soilcprocesses_tridiag(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(float*,float*,float*,float*,float*,int*,int*,int*,int*,int*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    float *a = NULL;
    npy_intp a_Dims[1] = {-1};
    const int a_Rank = 1;
    PyArrayObject *capi_a_as_array = NULL;
    int capi_a_intent = 0;
    PyObject *a_capi = Py_None;
    float *b = NULL;
    npy_intp b_Dims[1] = {-1};
    const int b_Rank = 1;
    PyArrayObject *capi_b_as_array = NULL;
    int capi_b_intent = 0;
    PyObject *b_capi = Py_None;
    float *c = NULL;
    npy_intp c_Dims[1] = {-1};
    const int c_Rank = 1;
    PyArrayObject *capi_c_as_array = NULL;
    int capi_c_intent = 0;
    PyObject *c_capi = Py_None;
    float *r = NULL;
    npy_intp r_Dims[1] = {-1};
    const int r_Rank = 1;
    PyArrayObject *capi_r_as_array = NULL;
    int capi_r_intent = 0;
    PyObject *r_capi = Py_None;
    float *u = NULL;
    npy_intp u_Dims[1] = {-1};
    const int u_Rank = 1;
    PyArrayObject *capi_u_as_array = NULL;
    int capi_u_intent = 0;
    int f2py_a_d0 = 0;
    int f2py_b_d0 = 0;
    int f2py_c_d0 = 0;
    int f2py_r_d0 = 0;
    int f2py_u_d0 = 0;
    static char *capi_kwlist[] = {"a","b","c","r",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOO|:turbation.soilcprocesses.tridiag",\
        capi_kwlist,&a_capi,&b_capi,&c_capi,&r_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable a */
    ;
    capi_a_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.tridiag: failed to create array from the 1st argument `a`";
    capi_a_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,a_Dims,a_Rank,  capi_a_intent,a_capi,capi_errmess);
    if (capi_a_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        a = (float *)(PyArray_DATA(capi_a_as_array));

    /* Processing variable b */
    ;
    capi_b_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.tridiag: failed to create array from the 2nd argument `b`";
    capi_b_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,b_Dims,b_Rank,  capi_b_intent,b_capi,capi_errmess);
    if (capi_b_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        b = (float *)(PyArray_DATA(capi_b_as_array));

    /* Processing variable c */
    ;
    capi_c_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.tridiag: failed to create array from the 3rd argument `c`";
    capi_c_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,c_Dims,c_Rank,  capi_c_intent,c_capi,capi_errmess);
    if (capi_c_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        c = (float *)(PyArray_DATA(capi_c_as_array));

    /* Processing variable r */
    ;
    capi_r_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.tridiag: failed to create array from the 4th argument `r`";
    capi_r_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,r_Dims,r_Rank,  capi_r_intent,r_capi,capi_errmess);
    if (capi_r_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        r = (float *)(PyArray_DATA(capi_r_as_array));

    /* Processing variable u */
    ;
    capi_u_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
    const char * capi_errmess = "turbation.turbation.soilcprocesses.tridiag: failed to create array from the hidden `u`";
    capi_u_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,u_Dims,u_Rank,  capi_u_intent,Py_None,capi_errmess);
    if (capi_u_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = turbation_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        u = (float *)(PyArray_DATA(capi_u_as_array));

    /* Processing variable f2py_a_d0 */
    f2py_a_d0 = shape(a, 0);
    /* Processing variable f2py_b_d0 */
    f2py_b_d0 = shape(b, 0);
    /* Processing variable f2py_c_d0 */
    f2py_c_d0 = shape(c, 0);
    /* Processing variable f2py_r_d0 */
    f2py_r_d0 = shape(r, 0);
    /* Processing variable f2py_u_d0 */
    f2py_u_d0 = shape(u, 0);
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
    (*f2py_func)(a,b,c,r,u,&f2py_a_d0,&f2py_b_d0,&f2py_c_d0,&f2py_r_d0,&f2py_u_d0);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("N",capi_u_as_array);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    /* End of cleaning variable f2py_u_d0 */
    /* End of cleaning variable f2py_r_d0 */
    /* End of cleaning variable f2py_c_d0 */
    /* End of cleaning variable f2py_b_d0 */
    /* End of cleaning variable f2py_a_d0 */
    }  /* if (capi_u_as_array == NULL) ... else of u */
    /* End of cleaning variable u */
    if((PyObject *)capi_r_as_array!=r_capi) {
        Py_XDECREF(capi_r_as_array); }
    }  /* if (capi_r_as_array == NULL) ... else of r */
    /* End of cleaning variable r */
    if((PyObject *)capi_c_as_array!=c_capi) {
        Py_XDECREF(capi_c_as_array); }
    }  /* if (capi_c_as_array == NULL) ... else of c */
    /* End of cleaning variable c */
    if((PyObject *)capi_b_as_array!=b_capi) {
        Py_XDECREF(capi_b_as_array); }
    }  /* if (capi_b_as_array == NULL) ... else of b */
    /* End of cleaning variable b */
    if((PyObject *)capi_a_as_array!=a_capi) {
        Py_XDECREF(capi_a_as_array); }
    }  /* if (capi_a_as_array == NULL) ... else of a */
    /* End of cleaning variable a */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of tridiag *******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/

static FortranDataDef f2py_soilcprocesses_def[] = {
  {"turbation",-1,{{-1}},0,0,NULL,(void *)f2py_rout_turbation_soilcprocesses_turbation,doc_f2py_rout_turbation_soilcprocesses_turbation},
  {"tridiag",-1,{{-1}},0,0,NULL,(void *)f2py_rout_turbation_soilcprocesses_tridiag,doc_f2py_rout_turbation_soilcprocesses_tridiag},
  {NULL}
};

static void f2py_setup_soilcprocesses(char *turbation,char *tridiag) {
  int i_f2py=0;
  f2py_soilcprocesses_def[i_f2py++].data = turbation;
  f2py_soilcprocesses_def[i_f2py++].data = tridiag;
}
extern void F_FUNC(f2pyinitsoilcprocesses,F2PYINITSOILCPROCESSES)(void (*)(char *,char *));
static void f2py_init_soilcprocesses(void) {
  F_FUNC(f2pyinitsoilcprocesses,F2PYINITSOILCPROCESSES)(f2py_setup_soilcprocesses);
}

/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {

/*eof routine_defs*/
    {NULL}
};

static PyMethodDef f2py_module_methods[] = {

    {NULL,NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "turbation",
    NULL,
    -1,
    f2py_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_turbation(void) {
    int i;
    PyObject *m,*d, *s, *tmp;
    m = turbation_module = PyModule_Create(&moduledef);
    Py_SET_TYPE(&PyFortran_Type, &PyType_Type);
    import_array();
    if (PyErr_Occurred())
        {PyErr_SetString(PyExc_ImportError, "can't initialize module turbation (failed to import numpy)"); return m;}
    d = PyModule_GetDict(m);
    s = PyUnicode_FromString("1.24.4");
    PyDict_SetItemString(d, "__version__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString(
        "This module 'turbation' is auto-generated with f2py (version:1.24.4).\nFunctions:\n"
"Fortran 90/95 modules:\n""  soilcprocesses --- turbation(),tridiag()"".");
    PyDict_SetItemString(d, "__doc__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString("1.24.4");
    PyDict_SetItemString(d, "__f2py_numpy_version__", s);
    Py_DECREF(s);
    turbation_error = PyErr_NewException ("turbation.error", NULL, NULL);
    /*
     * Store the error object inside the dict, so that it could get deallocated.
     * (in practice, this is a module, so it likely will not and cannot.)
     */
    PyDict_SetItemString(d, "_turbation_error", turbation_error);
    Py_DECREF(turbation_error);
    for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
        tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
        PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
        Py_DECREF(tmp);
    }


/*eof initf2pywraphooks*/
  PyDict_SetItemString(d, "soilcprocesses", PyFortranObject_New(f2py_soilcprocesses_def,f2py_init_soilcprocesses));
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
    if (! PyErr_Occurred())
        on_exit(f2py_report_on_exit,(void*)"turbation");
#endif
    return m;
}
#ifdef __cplusplus
}
#endif
