#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>

extern void slsqp_(
    int m, // standard int
    int meq, // standard int
    int la, // la = len(c). check len(c) >= la
    int n, // n = len(x). check len(x) >= n
    double *x, // array of length n    ---  value is returned to caller
    double *xl, // array of length n
    double *xu, // array of length n
    double f, // standard double
    double *c, // array of length la
    double *g, // array of length n + 1
    double **a, // matrix of dims (la, n + 1)
    double *acc, // standard double --- value is returned to caller
    int *iter, // standard int -- value is returned to caller
    int *mode, // standard int -- value is returned to caller
    double *w, // array of length l_w
    int l_w, // standard int. check len(w) >= l_w
    int *jw, // array of length l_jw
    int l_jw, // standard int. check len(jw) >= l_jw
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

int main(int argc, char *argv[])
{
    int la = 1; // la = len(c). check len(c) >= la
    int n = 1; // n = len(x). check len(x) >= n

    double **local2Darray = (double**) malloc(sizeof(double*) * la);
    int i,j;
    for (i = 0; i < la; i++)
    {
        local2Darray[i] = (double*) malloc(sizeof(double) * (n + 1));
        for (j = 0; j < n + 1; j++)
        {
            local2Darray[i][j] = j;
        }
    }

    int m = 1; // standard int
    int meq = 1; // standard int

    double x = 1; // array of length n    ---  value is returned to caller
    double xl = 1; // array of length n
    double xu = 1; // array of length n
    double f = 1; // standard double
    double c = 1; // array of length la
    double g[2] = {0, 1}; // array of length n + 1
    double **a = local2Darray; // matrix of dims (la, n + 1)
    double acc = 1; // standard double --- value is returned to caller
    int iter = 1; // standard int -- value is returned to caller
    int mode = 1; // standard int -- value is returned to caller
    double w = 1; // array of length l_w
    int l_w = 1; // standard int. check len(w) >= l_w
    int jw = 1; // array of length l_jw
    int l_jw = 1; // standard int. check len(jw) >= l_jw
    double alpha = 1; // standard double  -- value is returned to caller
    double f0 = 1; // standard double -- value is returned to caller
    double gs = 1; // standard double -- value is returned to caller
    double h1 = 1; // standard double -- value is returned to caller
    double h2 = 1; // standard double -- value is returned to caller
    double h3 = 1; // standard double -- value is returned to caller
    double h4 = 1; // standard double -- value is returned to caller
    double t = 1; // standard double -- value is returned to caller
    double t0 = 1; // standard double -- value is returned to caller
    double tol = 1; // standard double -- value is returned to caller
    int iexact = 1; // standard int -- value is returned to caller
    int incons = 1; // standard int -- value is returned to caller
    int ireset = 1; // standard int -- value is returned to caller
    int itermx = 1; // standard int -- value is returned to caller
    int line = 1; // standard int -- value is returned to caller
    int n1 = 1; // standard int -- value is returned to caller
    int n2 = 1; // standard int -- value is returned to caller
    int n3 = 1; // standard int -- value is returned to caller

    // Call the Fortran routine.
    slsqp_(
        m,
        meq,
        la,
        n,
        &x,
        &xl,
        &xu,
        f,
        &c,
        g,
        local2Darray,
        &acc,
        &iter,
        &mode,
        &w,
        l_w,
        &jw,
        l_jw,
        &alpha,
        &f0,
        &gs,
        &h1,
        &h2,
        &h3,
        &h4,
        &t,
        &t0,
        &tol,
        &iexact,
        &incons,
        &ireset,
        &itermx,
        &line,
        &n1,
        &n2,
        &n3
    );

    for (i = 0; i < la; i++)
    {
        free(local2Darray[i]);
    }

    free(local2Darray);

    return 0;
}
