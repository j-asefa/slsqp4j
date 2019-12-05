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
    double **a, // matrix of dims (la, n + 1)
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

JNIEXPORT jint JNICALL Java_com_example_slsqp_Slsqp_slsqp(
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
      jobjectArray a, // 2D array of dims (la, n + 1)
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


    double **local2Darray = (double**) malloc(sizeof(double*) * la);
    int i;
    for (i = 0; i < la; i++)
    {
        local2Darray[i] = (double*) malloc(sizeof(double) * (n + 1));
    }

    int j;
    for (i = 0; i < la; i++) {
         jdoubleArray oneDim = (jdoubleArray)(*env)->GetObjectArrayElement(env, a, i);
         jdouble *element = (*env)->GetDoubleArrayElements(env, oneDim, 0);
        for (j = 0; j < n + 1; j++) {
            local2Darray[i][j] = element[j];
        }
        (*env)->ReleaseDoubleArrayElements(env, oneDim, element, 0);
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

    for (i = 0; i < la; i++)
    {
        free(local2Darray[i]);
    }

    free(local2Darray);

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
