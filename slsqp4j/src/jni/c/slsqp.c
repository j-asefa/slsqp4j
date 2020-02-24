/*
 Copyright (c) 1988 Dieter Kraft

Copyright (c) 1994 Association for Computing Machinery

Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2019 SciPy Developers.
All rights reserved.

Copyright (c) 2020, Skew Ltd.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <jni.h>

extern void slsqp_(
    int* m, // standard int
    int* meq, // standard int
    int* la, // la = len(c). check len(c) >= la
    int* n, // n = len(x). check len(x) >= n
    double *x, // array of length n    ---  value is returned to caller
    double *xl, // array of length n
    double *xu, // array of length n
    double* f, // standard double
    double *c, // array of length la
    double *g, // array of length n + 1
    double *a, // matrix of dims (la, n + 1)
    double *acc, // standard double --- value is returned to caller
    int *iter, // standard int -- value is returned to caller
    int *mode, // standard int -- value is returned to caller
    double *w, // array of length l_w
    int* l_w, // standard int. check len(w) >= l_w
    int *jw, // array of length l_jw
    int* l_jw, // standard int. check len(jw) >= l_jw
    double *alpha, // standard double  -- value is returned to caller
    double *f0, // standard double -- value is returned to caller
    double *gs, // standard double -- value is returned to caller
    double *h1, // standard double -- value is returned to caller
    double *h2, // standard double -- value is returned to caller
    double *h3, // standard double -- value is returned to caller
    double *h4, // standard double -- value is returned to caller
    double *t, // standard double -- value is returned to caller
    double *t0, // standard double -- value is returned to caller
    double *tol, // standard double -- value is returned to caller
    int *iexact, // standard int -- value is returned to caller
    int *incons, // standard int -- value is returned to caller
    int *ireset, // standard int -- value is returned to caller
    int *itermx, // standard int -- value is returned to caller
    int *line, // standard int -- value is returned to caller
    int *n1, // standard int -- value is returned to caller
    int *n2, // standard int -- value is returned to caller
    int *n3 // standard int -- value is returned to caller
);

JNIEXPORT jint JNICALL Java_com_skew_slsqp4j_util_NativeUtils_slsqp(
      JNIEnv *env,         /* interface pointer */
      jclass cl,         /* "this" pointer */
      jint m, // standard int
      jint meq, // standard int
      jint la, // la = len(c). check len(c) >= la
      jint n, // n = len(x). check len(x) >= n
      jdoubleArray x, // array of length n    ---  value is returned to caller
      jdoubleArray xl, // array of length n
      jdoubleArray xu, // array of length n
      jdouble f, // standard double
      jdoubleArray c, // array of length la
      jdoubleArray g, // array of length n + 1
      jobjectArray a, // 2D array of dims (n + 1, la) -- this array is passed in *column-major* order
      jdoubleArray acc, // standard double --- value is returned to caller
      jintArray iter, // standard int -- value is returned to caller
      jintArray mode, // standard int -- value is returned to caller
      jdoubleArray w, // array of length l_w
      jint l_w, // standard int. check len(w) >= l_w
      jintArray jw, // array of length l_jw
      jint l_jw, // standard int. check len(jw) >= l_jw
      jdoubleArray alpha, // standard double  -- value is returned to caller
      jdoubleArray f0, // standard double -- value is returned to caller
      jdoubleArray gs, // standard double -- value is returned to caller
      jdoubleArray h1, // standard double -- value is returned to caller
      jdoubleArray h2, // standard double -- value is returned to caller
      jdoubleArray h3, // standard double -- value is returned to caller
      jdoubleArray h4, // standard double -- value is returned to caller
      jdoubleArray t, // standard double -- value is returned to caller
      jdoubleArray t0, // standard double -- value is returned to caller
      jdoubleArray tol, // standard double -- value is returned to caller
      jintArray iexact, // standard int -- value is returned to caller
      jintArray incons, // standard int -- value is returned to caller
      jintArray ireset, // standard int -- value is returned to caller
      jintArray itermx, // standard int -- value is returned to caller
      jintArray line, // standard int -- value is returned to caller
      jintArray n1, // standard int -- value is returned to caller
      jintArray n2, // standard int -- value is returned to caller
      jintArray n3 // standard int -- value is returned to caller
)
{
    jsize x_array_length = (*env)->GetArrayLength(env, x);
    if (x_array_length < n)
    {
        exit(-1);
    }
    jdouble *x_array = (*env)->GetDoubleArrayElements(env, x, 0);
    jdouble *xl_array = (*env)->GetDoubleArrayElements(env, xl, 0);
    jdouble *xu_array = (*env)->GetDoubleArrayElements(env, xu, 0);

    jsize c_array_length = (*env)->GetArrayLength(env, c);
    if (c_array_length < la)
    {
        exit(-1);
    }
    jdouble *c_array = (*env)->GetDoubleArrayElements(env, c, 0);
    jdouble *g_array = (*env)->GetDoubleArrayElements(env, g, 0);

    jdouble *acc_array = (*env)->GetDoubleArrayElements(env, acc, 0);

    jint *iter_array = (*env)->GetIntArrayElements(env, iter, 0);
    jint *mode_array = (*env)->GetIntArrayElements(env, mode, 0);

    jsize w_array_length = (*env)->GetArrayLength(env, w);
    if (w_array_length < l_w)
    {
        exit(-1);
    }
    jdouble* w_array = (*env)->GetDoubleArrayElements(env, w, 0);

    jsize jw_array_length = (*env)->GetArrayLength(env, jw);
    if (jw_array_length < l_jw)
    {
        exit(-1);
    }
    jint *jw_array = (*env)->GetIntArrayElements(env, jw, 0);

    jdouble *alpha_array = (*env)->GetDoubleArrayElements(env, alpha, 0);
    jdouble *f0_array = (*env)->GetDoubleArrayElements(env, f0, 0);
    jdouble *gs_array = (*env)->GetDoubleArrayElements(env, gs, 0);
    jdouble *h1_array = (*env)->GetDoubleArrayElements(env, h1, 0);
    jdouble *h2_array = (*env)->GetDoubleArrayElements(env, h2, 0);
    jdouble *h3_array = (*env)->GetDoubleArrayElements(env, h3, 0);
    jdouble *h4_array = (*env)->GetDoubleArrayElements(env, h4, 0);
    jdouble *t_array = (*env)->GetDoubleArrayElements(env, t, 0);
    jdouble *t0_array = (*env)->GetDoubleArrayElements(env, t0, 0);
    jdouble *tol_array = (*env)->GetDoubleArrayElements(env, tol, 0);

    jint *iexact_array = (*env)->GetIntArrayElements(env, iexact, 0);
    jint *incons_array = (*env)->GetIntArrayElements(env, incons, 0);
    jint *ireset_array = (*env)->GetIntArrayElements(env, ireset, 0);
    jint *itermx_array = (*env)->GetIntArrayElements(env, itermx, 0);
    jint *line_array = (*env)->GetIntArrayElements(env, line, 0);
    jint *n1_array = (*env)->GetIntArrayElements(env, n1, 0);
    jint *n2_array = (*env)->GetIntArrayElements(env, n2, 0);
    jint *n3_array = (*env)->GetIntArrayElements(env, n3, 0);

    double *local2Darray = (double*) malloc(sizeof(double) * la * (n + 1));
    int i, j;
    for (i = 0; i < n + 1; i++)
    {
        jdoubleArray oneDim = (jdoubleArray)(*env)->GetObjectArrayElement(env, a, i);
        for (j = 0; j < la; j++) {
            jdouble *element = (*env)->GetDoubleArrayElements(env, oneDim, 0);
            local2Darray[(i * la) + j] = element[j];
            (*env)->ReleaseDoubleArrayElements(env, oneDim, element, 0);
        }
        (*env)->DeleteLocalRef(env, oneDim);
    }

    // Call the Fortran routine.
    slsqp_(
        &m,
        &meq,
        &la,
        &n,
        x_array,
        xl_array,
        xu_array,
        &f,
        c_array,
        g_array,
        local2Darray,
        acc_array,
        iter_array,
        mode_array,
        w_array,
        &l_w,
        jw_array,
        &l_jw,
        alpha_array,
        f0_array,
        gs_array,
        h1_array,
        h2_array,
        h3_array,
        h4_array,
        t_array,
        t0_array,
        tol_array,
        iexact_array,
        incons_array,
        ireset_array,
        itermx_array,
        line_array,
        n1_array,
        n2_array,
        n3_array
    );

    for (i = 0; i < n + 1; i++) {
        jdoubleArray doubleArray = (*env)->NewDoubleArray(env, la);
        for (j = 0; j < la; j++)
        {
            double val = local2Darray[(i * la) + j];
            (*env)->SetDoubleArrayRegion(env, doubleArray, j, 1, &val);
        }

        (*env)->SetObjectArrayElement(env, a, i, doubleArray);
        (*env)->DeleteLocalRef(env, doubleArray);
    }

    free(local2Darray);

    (*env)->SetDoubleArrayRegion(env, x, 0, n, x_array);
    (*env)->SetDoubleArrayRegion(env, xl, 0, n, xl_array);
    (*env)->SetDoubleArrayRegion(env, xu, 0, n, xu_array);
    (*env)->SetDoubleArrayRegion(env, c, 0, la, c_array);
    (*env)->SetDoubleArrayRegion(env, g, 0, n + 1, g_array);
    (*env)->SetDoubleArrayRegion(env, acc, 0, 1, acc_array);
    (*env)->SetIntArrayRegion(env, iter, 0, 1, iter_array);
    (*env)->SetIntArrayRegion(env, mode, 0, 1, mode_array);
    (*env)->SetDoubleArrayRegion(env, w, 0, l_w, w_array);
    (*env)->SetIntArrayRegion(env, jw, 0, l_jw, jw_array);
    (*env)->SetDoubleArrayRegion(env, alpha, 0, 1, alpha_array);
    (*env)->SetDoubleArrayRegion(env, f0, 0, 1, f0_array);
    (*env)->SetDoubleArrayRegion(env, gs, 0, 1, gs_array);
    (*env)->SetDoubleArrayRegion(env, h1, 0, 1, h1_array);
    (*env)->SetDoubleArrayRegion(env, h2, 0, 1, h2_array);
    (*env)->SetDoubleArrayRegion(env, h3, 0, 1, h3_array);
    (*env)->SetDoubleArrayRegion(env, h4, 0, 1, h4_array);
    (*env)->SetDoubleArrayRegion(env, t, 0, 1, t_array);
    (*env)->SetDoubleArrayRegion(env, t0, 0, 1, t0_array);
    (*env)->SetDoubleArrayRegion(env, tol, 0, 1, tol_array);
    (*env)->SetIntArrayRegion(env, iexact, 0, 1, iexact_array);
    (*env)->SetIntArrayRegion(env, incons, 0, 1, incons_array);
    (*env)->SetIntArrayRegion(env, ireset, 0, 1, ireset_array);
    (*env)->SetIntArrayRegion(env, itermx, 0, 1, itermx_array);
    (*env)->SetIntArrayRegion(env, line, 0, 1, line_array);
    (*env)->SetIntArrayRegion(env, n1, 0, 1, n1_array);
    (*env)->SetIntArrayRegion(env, n2, 0, 1, n2_array);
    (*env)->SetIntArrayRegion(env, n3, 0, 1, n3_array);


    (*env)->ReleaseDoubleArrayElements(env, x, x_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, xl, xl_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, xu, xu_array, JNI_ABORT);

    (*env)->ReleaseDoubleArrayElements(env, c, c_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, g, g_array, JNI_ABORT);

    (*env)->ReleaseDoubleArrayElements(env, acc, acc_array, JNI_ABORT);

    (*env)->ReleaseIntArrayElements(env, iter, iter_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, mode, mode_array, JNI_ABORT);

    (*env)->ReleaseDoubleArrayElements(env, w, w_array, JNI_ABORT);

    (*env)->ReleaseIntArrayElements(env, jw, jw_array, JNI_ABORT);

    (*env)->ReleaseDoubleArrayElements(env, alpha, alpha_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, f0, f0_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, gs, gs_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, h1, h1_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, h2, h2_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, h3, h3_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, h4, h4_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, t, t_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, t0, t0_array, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, tol, tol_array, JNI_ABORT);

    (*env)->ReleaseIntArrayElements(env, iexact, iexact_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, incons, incons_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, ireset, ireset_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, itermx, itermx_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, line, line_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, n1, n1_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, n2, n2_array, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, n3, n3_array, JNI_ABORT);

    return 0;
}
