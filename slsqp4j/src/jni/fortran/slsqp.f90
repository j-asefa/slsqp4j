!      Copyright (c) 1988 Dieter Kraft
!
!      Copyright (c) 1994 Association for Computing Machinery
!
!      Copyright (c) 2001, 2002 Enthought, Inc.
!      All rights reserved.
!
!      Copyright (c) 2003-2019 SciPy Developers.
!      All rights reserved.
!
!      Copyright (c) 2020, Skew Ltd.
!      All rights reserved.
!
!      Redistribution and use in source and binary forms, with or without
!      modification, are permitted provided that the following conditions are met:
!
!      1. Redistributions of source code must retain the above copyright notice, this
!      list of conditions and the following disclaimer.
!
!      2. Redistributions in binary form must reproduce the above copyright notice,
!      this list of conditions and the following disclaimer in the documentation
!      and/or other materials provided with the distribution.
!
!      3. Neither the name of the copyright holder nor the names of its
!      contributors may be used to endorse or promote products derived from
!      this software without specific prior written permission.
!
!      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!      AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!      IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!      DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!      FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!      DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!      SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!      CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!      OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!      --------------------------------------------------------------------------------
!
!      ALGORITHM 733, COLLECTED ALGORITHMS FROM ACM.
!      TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 262-281.
!      https://doi.org/10.1145/192115.192124
!
!
!      https://web.archive.org/web/20170106155705/http://permalink.gmane.org/gmane.comp.python.scientific.devel/6725
!      ------
!      From: Deborah Cotton <cotton@hq.acm.org>
!      Date: Fri, 14 Sep 2007 12:35:55 -0500
!      Subject: RE: Algorithm License requested
!      To: Alan Isaac
!
!      Prof. Issac,
!
!      In that case, then because the author consents to [the ACM] releasing
!      the code currently archived at http://www.netlib.org/toms/733 under the
!      BSD license, the ACM hereby releases this code under the BSD license.
!
!      Regards,
!
!      Deborah Cotton, Copyright & Permissions
!      ACM Publications
!      2 Penn Plaza, Suite 701**
!      New York, NY 10121-0701
!      permissions@acm.org
!      212.869.7440 ext. 652
!      Fax. 212.869.0481
!      ------
!

!***********************************************************************
!                              optimizer                               *
!***********************************************************************

SUBROUTINE slsqp (m, meq, la, n, x, xl, xu, f, c, g, a, &
        acc, iter, mode, w, l_w, jw, l_jw, &
        alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, &
        iexact, incons, ireset, itermx, line, &
        n1, n2, n3)

    !   SLSQP       S EQUENTIAL  L EAST  SQ UARES  P ROGRAMMING
    !            TO SOLVE GENERAL NONLINEAR OPTIMIZATION PROBLEMS

    !***********************************************************************
    !*                                                                     *
    !*                                                                     *
    !*            A NONLINEAR PROGRAMMING METHOD WITH                      *
    !*            QUADRATIC  PROGRAMMING  SUBPROBLEMS                      *
    !*                                                                     *
    !*                                                                     *
    !*  THIS SUBROUTINE SOLVES THE GENERAL NONLINEAR PROGRAMMING PROBLEM   *
    !*                                                                     *
    !*            MINIMIZE    F(X)                                         *
    !*                                                                     *
    !*            SUBJECT TO  C (X) .EQ. 0  ,  J = 1,...,MEQ               *
    !*                         J                                           *
    !*                                                                     *
    !*                        C (X) .GE. 0  ,  J = MEQ+1,...,M             *
    !*                         J                                           *
    !*                                                                     *
    !*                        XL .LE. X .LE. XU , I = 1,...,N.             *
    !*                          I      I       I                           *
    !*                                                                     *
    !*  THE ALGORITHM IMPLEMENTS THE METHOD OF HAN AND POWELL              *
    !*  WITH BFGS-UPDATE OF THE B-MATRIX AND L1-TEST FUNCTION              *
    !*  WITHIN THE STEPLENGTH ALGORITHM.                                   *
    !*                                                                     *
    !*    PARAMETER DESCRIPTION:                                           *
    !*    ( * MEANS THIS PARAMETER WILL BE CHANGED DURING CALCULATION )    *
    !*                                                                     *
    !*    M              IS THE TOTAL NUMBER OF CONSTRAINTS, M .GE. 0      *
    !*    MEQ            IS THE NUMBER OF EQUALITY CONSTRAINTS, MEQ .GE. 0 *
    !*    LA             SEE A, LA .GE. MAX(M,1)                           *
    !*    N              IS THE NUMBER OF VARIBLES, N .GE. 1               *
    !*  * X()            X() STORES THE CURRENT ITERATE OF THE N VECTOR X  *
    !*                   ON ENTRY X() MUST BE INITIALIZED. ON EXIT X()     *
    !*                   STORES THE SOLUTION VECTOR X IF MODE = 0.         *
    !*    XL()           XL() STORES AN N VECTOR OF LOWER BOUNDS XL TO X.  *
    !*                   ELEMENTS MAY BE NAN TO INDICATE NO LOWER BOUND.   *
    !*    XU()           XU() STORES AN N VECTOR OF UPPER BOUNDS XU TO X.  *
    !*                   ELEMENTS MAY BE NAN TO INDICATE NO UPPER BOUND.   *
    !*    F              IS THE VALUE OF THE OBJECTIVE FUNCTION.           *
    !*    C()            C() STORES THE M VECTOR C OF CONSTRAINTS,         *
    !*                   EQUALITY CONSTRAINTS (IF ANY) FIRST.              *
    !*                   DIMENSION OF C MUST BE GREATER OR EQUAL LA,       *
    !*                   which must be GREATER OR EQUAL MAX(1,M).          *
    !*    G()            G() STORES THE N VECTOR G OF PARTIALS OF THE      *
    !*                   OBJECTIVE FUNCTION; DIMENSION OF G MUST BE        *
    !*                   GREATER OR EQUAL N+1.                             *
    !*    A(),LA,M,N     THE LA BY N + 1 ARRAY A() STORES                  *
    !*                   THE M BY N MATRIX A OF CONSTRAINT NORMALS.        *
    !*                   A() HAS FIRST DIMENSIONING PARAMETER LA,          *
    !*                   WHICH MUST BE GREATER OR EQUAL MAX(1,M).          *
    !*    F,C,G,A        MUST ALL BE SET BY THE USER BEFORE EACH CALL.     *
    !*  * ACC            ABS(ACC) CONTROLS THE FINAL ACCURACY.             *
    !*                   IF ACC .LT. ZERO AN EXACT LINESEARCH IS PERFORMED,*
    !*                   OTHERWISE AN ARMIJO-TYPE LINESEARCH IS USED.      *
    !*  * ITER           PRESCRIBES THE MAXIMUM NUMBER OF ITERATIONS.      *
    !*                   ON EXIT ITER INDICATES THE NUMBER OF ITERATIONS.  *
    !*  * MODE           MODE CONTROLS CALCULATION:                        *
    !*                   REVERSE COMMUNICATION IS USED IN THE SENSE THAT   *
    !*                   THE PROGRAM IS INITIALIZED BY MODE = 0; THEN IT IS*
    !*                   TO BE CALLED REPEATEDLY BY THE USER UNTIL A RETURN*
    !*                   WITH MODE .NE. IABS(1) TAKES PLACE.               *
    !*                   IF MODE = -1 GRADIENTS HAVE TO BE CALCULATED,     *
    !*                   WHILE WITH MODE = 1 FUNCTIONS HAVE TO BE CALCULATED
    !*                   MODE MUST NOT BE CHANGED BETWEEN SUBSEQUENT CALLS *
    !*                   OF SQP.                                           *
    !*                   EVALUATION MODES:                                 *
    !*        MODE = -1: GRADIENT EVALUATION, (G&A)                        *
    !*                0: ON ENTRY: INITIALIZATION, (F,G,C&A)               *
    !*                   ON EXIT : REQUIRED ACCURACY FOR SOLUTION OBTAINED *
    !*                1: FUNCTION EVALUATION, (F&C)                        *
    !*                                                                     *
    !*                   FAILURE MODES:                                    *
    !*                2: NUMBER OF EQUALITY CONSTRAINTS LARGER THAN N      *
    !*                3: MORE THAN 3*N ITERATIONS IN LSQ SUBPROBLEM        *
    !*                4: INEQUALITY CONSTRAINTS INCOMPATIBLE               *
    !*                5: SINGULAR MATRIX E IN LSQ SUBPROBLEM               *
    !*                6: SINGULAR MATRIX C IN LSQ SUBPROBLEM               *
    !*                7: RANK-DEFICIENT EQUALITY CONSTRAINT SUBPROBLEM HFTI*
    !*                8: POSITIVE DIRECTIONAL DERIVATIVE FOR LINESEARCH    *
    !*                9: MORE THAN ITER ITERATIONS IN SQP                  *
    !*             >=10: WORKING SPACE W OR JW TOO SMALL,                  *
    !*                   W SHOULD BE ENLARGED TO L_W=MODE/1000             *
    !*                   JW SHOULD BE ENLARGED TO L_JW=MODE-1000*L_W       *
    !*  * W(), L_W       W() IS A ONE DIMENSIONAL WORKING SPACE,           *
    !*                   THE LENGTH L_W OF WHICH SHOULD BE AT LEAST        *
    !*                   (3*N1+M)*(N1+1)                        for LSQ    *
    !*                  +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ         for LSI    *
    !*                  +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI   *
    !*                  + N1*N/2 + 2*M + 3*N + 3*N1 + 1         for SLSQPB *
    !*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          *
    !*        NOTICE:    FOR PROPER DIMENSIONING OF W IT IS RECOMMENDED TO *
    !*                   COPY THE FOLLOWING STATEMENTS INTO THE HEAD OF    *
    !*                   THE CALLING PROGRAM (AND REMOVE THE COMMENT C)    *
    !#######################################################################
    !     INTEGER LEN_W, LEN_JW, M, N, N1, MEQ, MINEQ
    !     PARAMETER (M=... , MEQ=... , N=...  )
    !     PARAMETER (N1= N+1, MINEQ= M-MEQ+N1+N1)
    !     PARAMETER (LEN_W=
    !    $           (3*N1+M)*(N1+1)
    !    $          +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
    !    $          +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1
    !    $          +(N+1)*N/2 + 2*M + 3*N + 3*N1 + 1,
    !    $           LEN_JW=MINEQ)
    !     DOUBLE PRECISION W(LEN_W)
    !     INTEGER          JW(LEN_JW)
    !#######################################################################
    !*                   THE FIRST M+N+N*N1/2 ELEMENTS OF W MUST NOT BE    *
    !*                   CHANGED BETWEEN SUBSEQUENT CALLS OF SLSQP.        *
    !*                   ON RETURN W(1) ... W(M) CONTAIN THE MULTIPLIERS   *
    !*                   ASSOCIATED WITH THE GENERAL CONSTRAINTS, WHILE    *
    !*                   W(M+1) ... W(M+N(N+1)/2) STORE THE CHOLESKY FACTOR*
    !*                   L*D*L(T) OF THE APPROXIMATE HESSIAN OF THE        *
    !*                   LAGRANGIAN COLUMNWISE DENSE AS LOWER TRIANGULAR   *
    !*                   UNIT MATRIX L WITH D IN ITS 'DIAGONAL' and        *
    !*                   W(M+N(N+1)/2+N+2 ... W(M+N(N+1)/2+N+2+M+2N)       *
    !*                   CONTAIN THE MULTIPLIERS ASSOCIATED WITH ALL       *
    !*                   ALL CONSTRAINTS OF THE QUADRATIC PROGRAM FINDING  *
    !*                   THE SEARCH DIRECTION TO THE SOLUTION X*           *
    !*  * JW(), L_JW     JW() IS A ONE DIMENSIONAL INTEGER WORKING SPACE   *
    !*                   THE LENGTH L_JW OF WHICH SHOULD BE AT LEAST       *
    !*                   MINEQ                                             *
    !*                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          *
    !*                                                                     *
    !*  THE USER HAS TO PROVIDE THE FOLLOWING SUBROUTINES:                 *
    !*     LDL(N,A,Z,SIG,W) :   UPDATE OF THE LDL'-FACTORIZATION.          *
    !*     LINMIN(A,B,F,TOL) :  LINESEARCH ALGORITHM IF EXACT = 1          *
    !*     LSQ(M,MEQ,LA,N,NC,C,D,A,B,XL,XU,X,LAMBDA,W,....) :              *
    !*                                                                     *
    !*        SOLUTION OF THE QUADRATIC PROGRAM                            *
    !*                QPSOL IS RECOMMENDED:                                *
    !*     PE GILL, W MURRAY, MA SAUNDERS, MH WRIGHT:                      *
    !*     USER'S GUIDE FOR SOL/QPSOL:                                     *
    !*     A FORTRAN PACKAGE FOR QUADRATIC PROGRAMMING,                    *
    !*     TECHNICAL REPORT SOL 83-7, JULY 1983                            *
    !*     DEPARTMENT OF OPERATIONS RESEARCH, STANFORD UNIVERSITY          *
    !*     STANFORD, CA 94305                                              *
    !*     QPSOL IS THE MOST ROBUST AND EFFICIENT QP-SOLVER                *
    !*     AS IT ALLOWS WARM STARTS WITH PROPER WORKING SETS               *
    !*                                                                     *
    !*     IF IT IS NOT AVAILABLE USE LSEI, A CONSTRAINT LINEAR LEAST      *
    !*     SQUARES SOLVER IMPLEMENTED USING THE SOFTWARE HFTI, LDP, NNLS   *
    !*     FROM C.L. LAWSON, R.J.HANSON: SOLVING LEAST SQUARES PROBLEMS,   *
    !*     PRENTICE HALL, ENGLEWOOD CLIFFS, 1974.                          *
    !*     LSEI COMES WITH THIS PACKAGE, together with all necessary SR's. *
    !*                                                                     *
    !*     TOGETHER WITH A COUPLE OF SUBROUTINES FROM BLAS LEVEL 1         *
    !*                                                                     *
    !*     SQP IS HEAD SUBROUTINE FOR BODY SUBROUTINE SQPBDY               *
    !*     IN WHICH THE ALGORITHM HAS BEEN IMPLEMENTED.                    *
    !*                                                                     *
    !*  IMPLEMENTED BY: DIETER KRAFT, DFVLR OBERPFAFFENHOFEN               *
    !*  as described in Dieter Kraft: A Software Package for               *
    !*                                Sequential Quadratic Programming     *
    !*                                DFVLR-FB 88-28, 1988                 *
    !*  which should be referenced if the user publishes results of SLSQP  *
    !*                                                                     *
    !*  DATE:           APRIL - OCTOBER, 1981.                             *
    !*  STATUS:         DECEMBER, 31-ST, 1984.                             *
    !*  STATUS:         MARCH   , 21-ST, 1987, REVISED TO FORTRAN 77       *
    !*  STATUS:         MARCH   , 20-th, 1989, REVISED TO MS-FORTRAN       *
    !*  STATUS:         APRIL   , 14-th, 1989, HESSE   in-line coded       *
    !*  STATUS:         FEBRUARY, 28-th, 1991, FORTRAN/2 Version 1.04      *
    !*                                         accepts Statement Functions *
    !*  STATUS:         MARCH   ,  1-st, 1991, tested with SALFORD         *
    !*                                         FTN77/386 COMPILER VERS 2.40*
    !*                                         in protected mode           *
    !*                                                                     *
    !***********************************************************************
    !*                                                                     *
    !*  Copyright 1991: Dieter Kraft, FHM                                  *
    !*                                                                     *
    !***********************************************************************

    INTEGER          il, im, ir, is, iter, iu, iv, iw, ix, l_w, l_jw, &
            jw(l_jw), la, m, meq, mineq, mode, n

    DOUBLE PRECISION acc, a(la, n + 1), c(la), f, g(n + 1), &
            x(n), xl(n), xu(n), w(l_w)

    INTEGER          iexact, incons, ireset, itermx, line, n1, n2, n3

    DOUBLE PRECISION alpha, f0, gs, h1, h2, h3, h4, t, t0, tol

    !     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ
    !                    +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ          for LSI
    !                    +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1        for LSEI
    !                    + N1*N/2 + 2*M + 3*N +3*N1 + 1           for SLSQPB
    !                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1

    !   CHECK LENGTH OF WORKING ARRAYS

    n1 = n + 1
    mineq = m - meq + n1 + n1
    il = (3 * n1 + m) * (n1 + 1) + &
            (n1 - meq + 1) * (mineq + 2) + 2 * mineq + &
            (n1 + mineq) * (n1 - meq) + 2 * meq + &
            n1 * n / 2 + 2 * m + 3 * n + 4 * n1 + 1
    im = MAX(mineq, n1 - meq)

    IF (l_w .LT. il .OR. l_jw .LT. im) THEN
        mode = 1000 * MAX(10, il)
        mode = mode + MAX(10, im)
        RETURN
    ENDIF

    !   PREPARE DATA FOR CALLING SQPBDY  -  INITIAL ADDRESSES IN W

    im = 1
    il = im + MAX(1, m)
    il = im + la
    ix = il + n1 * n / 2 + 1
    ir = ix + n
    is = ir + n + n + MAX(1, m)
    is = ir + n + n + la
    iu = is + n1
    iv = iu + n1
    iw = iv + n1

    CALL slsqpb  (m, meq, la, n, x, xl, xu, f, c, g, a, acc, iter, &
            mode, w(ir), w(il), w(ix), w(im), w(is), w(iu), w(iv), w(iw), jw, &
            alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, &
            iexact, incons, ireset, itermx, line, &
            n1, n2, n3)

END

SUBROUTINE slsqpb (m, meq, la, n, x, xl, xu, f, c, g, a, acc, &
        iter, mode, r, l, x0, mu, s, u, v, w, iw, &
        alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, &
        iexact, incons, ireset, itermx, line, &
        n1, n2, n3)

    !   NONLINEAR PROGRAMMING BY SOLVING SEQUENTIALLY QUADRATIC PROGRAMS

    !        -  L1 - LINE SEARCH,  POSITIVE DEFINITE  BFGS UPDATE  -

    !                      BODY SUBROUTINE FOR SLSQP

    INTEGER          iw(*), i, iexact, incons, ireset, iter, itermx, &
            k, j, la, line, m, meq, mode, n, n1, n2, n3
    LOGICAL          badlin

    DOUBLE PRECISION a(la, n + 1), c(la), g(n + 1), l((n + 1) * (n + 2) / 2), &
            mu(la), r(m + n + n + 2), s(n + 1), u(n + 1), v(n + 1), w(*), &
            x(n), xl(n), xu(n), x0(n), &
            ddot_sl, dnrm2_, linmin, &
            acc, alfmin, alpha, f, f0, gs, h1, h2, h3, h4, &
            hun, one, t, t0, ten, tol, two, ZERO

    !     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ
    !                     +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
    !                     +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI
    !                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1

    DATA             ZERO /0.0d0/, one /1.0d0/, alfmin /1.0d-1/, &
            hun /1.0d+2/, ten /1.0d+1/, two /2.0d0/

    !     The badlin flag keeps track whether the SQP problem on the current
    !     iteration was inconsistent or not.
    badlin = .false.

    IF (mode) 260, 100, 220

    100 itermx = iter
    IF (acc.GE.ZERO) THEN
        iexact = 0
    ELSE
        iexact = 1
    ENDIF
    acc = ABS(acc)
    tol = ten * acc
    iter = 0
    ireset = 0
    n1 = n + 1
    n2 = n1 * n / 2
    n3 = n2 + 1
    s(1) = ZERO
    mu(1) = ZERO
    CALL dcopy_(n, s(1), 0, s, 1)
    CALL dcopy_(m, mu(1), 0, mu, 1)

    !   RESET BFGS MATRIX

    110 ireset = ireset + 1
    IF (ireset.GT.5) GO TO 255
    l(1) = ZERO
    CALL dcopy_(n2, l(1), 0, l, 1)
    j = 1
    DO 120 i = 1, n
        l(j) = one
        j = j + n1 - i
    120 CONTINUE

    !   MAIN ITERATION : SEARCH DIRECTION, STEPLENGTH, LDL'-UPDATE

    130 iter = iter + 1
    mode = 9
    IF (iter.GT.itermx) GO TO 330

    !   SEARCH DIRECTION AS SOLUTION OF QP - SUBPROBLEM

    CALL dcopy_(n, xl, 1, u, 1)
    CALL dcopy_(n, xu, 1, v, 1)
    CALL daxpy_sl(n, -one, x, 1, u, 1)
    CALL daxpy_sl(n, -one, x, 1, v, 1)
    h4 = one
    CALL lsq (m, meq, n, n3, la, l, g, a, c, u, v, s, r, w, iw, mode)

    !   AUGMENTED PROBLEM FOR INCONSISTENT LINEARIZATION
    !
    !   If it turns out that the original SQP problem is inconsistent,
    !   disallow termination with convergence on this iteration,
    !   even if the augmented problem was solved.

    badlin = .false.
    IF (mode.EQ.6) THEN
        IF (n.EQ.meq) THEN
            mode = 4
        ENDIF
    ENDIF
    IF (mode.EQ.4) THEN
        badlin = .true.
        DO 140 j = 1, m
            IF (j.LE.meq) THEN
                a(j, n1) = -c(j)
            ELSE
                a(j, n1) = MAX(-c(j), ZERO)
            ENDIF
        140     CONTINUE
        s(1) = ZERO
        CALL dcopy_(n, s(1), 0, s, 1)
        h3 = ZERO
        g(n1) = ZERO
        l(n3) = hun
        s(n1) = one
        u(n1) = ZERO
        v(n1) = one
        incons = 0
        150     CALL lsq (m, meq, n1, n3, la, l, g, a, c, u, v, s, r, &
                w, iw, mode)
        h4 = one - s(n1)
        IF (mode.EQ.4) THEN
            l(n3) = ten * l(n3)
            incons = incons + 1
            IF (incons.GT.5) GO TO 330
            GOTO 150
        ELSE IF (mode.NE.1) THEN
            GOTO 330
        ENDIF
    ELSE IF (mode.NE.1) THEN
        GOTO 330
    ENDIF

    !   UPDATE MULTIPLIERS FOR L1-TEST

    DO 160 i = 1, n
        v(i) = g(i) - ddot_sl(m, a(1, i), 1, r, 1)
    160 CONTINUE
    f0 = f
    CALL dcopy_(n, x, 1, x0, 1)
    gs = ddot_sl(n, g, 1, s, 1)
    h1 = ABS(gs)
    h2 = ZERO
    DO 170 j = 1, m
        IF (j.LE.meq) THEN
            h3 = c(j)
        ELSE
            h3 = ZERO
        ENDIF
        h2 = h2 + MAX(-c(j), h3)
        h3 = ABS(r(j))
        mu(j) = MAX(h3, (mu(j) + h3) / two)
        h1 = h1 + h3 * ABS(c(j))
    170 CONTINUE

    !   CHECK CONVERGENCE

    mode = 0
    IF (h1.LT.acc .AND. h2.LT.acc .AND. .NOT. badlin&
            .AND. f .EQ. f) GO TO 330
    h1 = ZERO
    DO 180 j = 1, m
        IF (j.LE.meq) THEN
            h3 = c(j)
        ELSE
            h3 = ZERO
        ENDIF
        h1 = h1 + mu(j) * MAX(-c(j), h3)
    180 CONTINUE
    t0 = f + h1
    h3 = gs - h1 * h4
    mode = 8
    IF (h3.GE.ZERO) GO TO 110

    !   LINE SEARCH WITH AN L1-TESTFUNCTION

    line = 0
    alpha = one
    IF (iexact.EQ.1) GOTO 210

    !   INEXACT LINESEARCH

    190     line = line + 1
    h3 = alpha * h3
    CALL dscal_sl(n, alpha, s, 1)
    CALL dcopy_(n, x0, 1, x, 1)
    CALL daxpy_sl(n, one, s, 1, x, 1)
    mode = 1
    GO TO 330
    200         IF (h1.LE.h3 / ten .OR. line.GT.10) GO TO 240
    alpha = MAX(h3 / (two * (h3 - h1)), alfmin)
    GO TO 190

    !   EXACT LINESEARCH

    210 IF (line.NE.3) THEN
        alpha = linmin(line, alfmin, one, t, tol)
        CALL dcopy_(n, x0, 1, x, 1)
        CALL daxpy_sl(n, alpha, s, 1, x, 1)
        mode = 1
        GOTO 330
    ENDIF
    CALL dscal_sl(n, alpha, s, 1)
    GOTO 240

    !   CALL FUNCTIONS AT CURRENT X

    220     t = f
    DO 230 j = 1, m
        IF (j.LE.meq) THEN
            h1 = c(j)
        ELSE
            h1 = ZERO
        ENDIF
        t = t + mu(j) * MAX(-c(j), h1)
    230     CONTINUE
    h1 = t - t0
    GOTO (200, 210) iexact + 1

    !   CHECK CONVERGENCE

    240 h3 = ZERO
    DO 250 j = 1, m
        IF (j.LE.meq) THEN
            h1 = c(j)
        ELSE
            h1 = ZERO
        ENDIF
        h3 = h3 + MAX(-c(j), h1)
    250 CONTINUE
    IF ((ABS(f - f0).LT.acc .OR. dnrm2_(n, s, 1).LT.acc) .AND. h3.LT.acc&
            .AND. .NOT. badlin .AND. f .EQ. f)&
            THEN
        mode = 0
    ELSE
        mode = -1
    ENDIF
    GO TO 330

    !   CHECK relaxed CONVERGENCE in case of positive directional derivative

    255 CONTINUE
    h3 = ZERO
    DO 256 j = 1, m
        IF (j.LE.meq) THEN
            h1 = c(j)
        ELSE
            h1 = ZERO
        ENDIF
        h3 = h3 + MAX(-c(j), h1)
    256 CONTINUE
    IF ((ABS(f - f0).LT.tol .OR. dnrm2_(n, s, 1).LT.tol) .AND. h3.LT.tol&
            .AND. .NOT. badlin .AND. f .EQ. f)&
            THEN
        mode = 0
    ELSE
        mode = 8
    ENDIF
    GO TO 330

    !   CALL JACOBIAN AT CURRENT X

    !   UPDATE CHOLESKY-FACTORS OF HESSIAN MATRIX BY MODIFIED BFGS FORMULA

    260 DO 270 i = 1, n
        u(i) = g(i) - ddot_sl(m, a(1, i), 1, r, 1) - v(i)
    270 CONTINUE

    !   L'*S

    k = 0
    DO 290 i = 1, n
        h1 = ZERO
        k = k + 1
        DO 280 j = i + 1, n
            k = k + 1
            h1 = h1 + l(k) * s(j)
        280    CONTINUE
        v(i) = s(i) + h1
    290 CONTINUE

    !   D*L'*S

    k = 1
    DO 300 i = 1, n
        v(i) = l(k) * v(i)
        k = k + n1 - i
    300 CONTINUE

    !   L*D*L'*S

    DO 320 i = n, 1, -1
        h1 = ZERO
        k = i
        DO 310 j = 1, i - 1
            h1 = h1 + l(k) * v(j)
            k = k + n - j
        310    CONTINUE
        v(i) = v(i) + h1
    320 CONTINUE

    h1 = ddot_sl(n, s, 1, u, 1)
    h2 = ddot_sl(n, s, 1, v, 1)
    h3 = 0.2d0 * h2
    IF (h1.LT.h3) THEN
        h4 = (h2 - h3) / (h2 - h1)
        h1 = h3
        CALL dscal_sl(n, h4, u, 1)
        CALL daxpy_sl(n, one - h4, v, 1, u, 1)
    ENDIF
    IF (h1.EQ.0 .or. h2.EQ.0) THEN
        !         Singular update: reset hessian.
        GO TO 110
    end if
    CALL ldl(n, l, u, +one / h1, v)
    CALL ldl(n, l, v, -one / h2, u)

    !   END OF MAIN ITERATION

    GO TO 130

    !   END OF SLSQPB

330 END


SUBROUTINE lsq(m, meq, n, nl, la, l, g, a, b, xl, xu, x, y, w, jw, mode)

    !   MINIMIZE with respect to X

    !             ||E*X - F||
    !                                      1/2  T
    !   WITH UPPER TRIANGULAR MATRIX E = +D   *L ,

    !                                      -1/2  -1
    !                     AND VECTOR F = -D    *L  *G,

    !  WHERE THE UNIT LOWER TRIDIANGULAR MATRIX L IS STORED COLUMNWISE
    !  DENSE IN THE N*(N+1)/2 ARRAY L WITH VECTOR D STORED IN ITS
    ! 'DIAGONAL' THUS SUBSTITUTING THE ONE-ELEMENTS OF L

    !   SUBJECT TO

    !             A(J)*X - B(J) = 0 ,         J=1,...,MEQ,
    !             A(J)*X - B(J) >=0,          J=MEQ+1,...,M,
    !             XL(I) <= X(I) <= XU(I),     I=1,...,N,
    !     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS L, G, A, B, XL, XU.
    !     WITH DIMENSIONS: L(N*(N+1)/2), G(N), A(LA,N), B(M), XL(N), XU(N)
    !     THE WORKING ARRAY W MUST HAVE AT LEAST THE FOLLOWING DIMENSION:
    !     DIM(W) =        (3*N+M)*(N+1)                        for LSQ
    !                    +(N-MEQ+1)*(MINEQ+2) + 2*MINEQ        for LSI
    !                    +(N+MINEQ)*(N-MEQ) + 2*MEQ + N        for LSEI
    !                      with MINEQ = M - MEQ + 2*N
    !     ON RETURN, NO ARRAY WILL BE CHANGED BY THE SUBROUTINE.
    !     X     STORES THE N-DIMENSIONAL SOLUTION VECTOR
    !     Y     STORES THE VECTOR OF LAGRANGE MULTIPLIERS OF DIMENSION
    !           M+N+N (CONSTRAINTS+LOWER+UPPER BOUNDS)
    !     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
    !               5: MATRIX E IS NOT OF FULL RANK
    !               6: MATRIX C IS NOT OF FULL RANK
    !               7: RANK DEFECT IN HFTI

    !     coded            Dieter Kraft, april 1987
    !     revised                        march 1989

    DOUBLE PRECISION l, g, a, b, w, xl, xu, x, y, &
            diag, ZERO, one, ddot_sl, xnorm

    INTEGER          jw(*), i, ic, id, ie, IF, ig, ih, il, im, ip, iu, iw, &
            i1, i2, i3, i4, la, m, meq, mineq, mode, m1, n, nl, n1, n2, n3, &
            nancnt, j

    DIMENSION        a(la, n), b(la), g(n), l(nl), &
            w(*), x(n), xl(n), xu(n), y(m + n + n)

    DATA             ZERO/0.0d0/, one/1.0d0/

    n1 = n + 1
    mineq = m - meq
    m1 = mineq + n + n

    !  determine whether to solve problem
    !  with inconsistent linerarization (n2=1)
    !  or not (n2=0)

    n2 = n1 * n / 2 + 1
    IF (n2.EQ.nl) THEN
        n2 = 0
    ELSE
        n2 = 1
    ENDIF
    n3 = n - n2

    !  RECOVER MATRIX E AND VECTOR F FROM L AND G

    i2 = 1
    i3 = 1
    i4 = 1
    ie = 1
    IF = n * n + 1
    DO 10 i = 1, n3
        i1 = n1 - i
        diag = SQRT (l(i2))
        w(i3) = ZERO
        CALL dcopy_ (i1, w(i3), 0, w(i3), 1)
        CALL dcopy_ (i1 - n2, l(i2), 1, w(i3), n)
        CALL dscal_sl (i1 - n2, diag, w(i3), n)
        w(i3) = diag
        w(IF - 1 + i) = (g(i) - ddot_sl (i - 1, w(i4), 1, w(IF), 1)) / diag
        i2 = i2 + i1 - n2
        i3 = i3 + n1
        i4 = i4 + n
    10 CONTINUE
    IF (n2.EQ.1) THEN
        w(i3) = l(nl)
        w(i4) = ZERO
        CALL dcopy_ (n3, w(i4), 0, w(i4), 1)
        w(IF - 1 + n) = ZERO
    ENDIF
    CALL dscal_sl (n, - one, w(IF), 1)

    ic = IF + n
    id = ic + meq * n

    IF (meq .GT. 0) THEN

        !  RECOVER MATRIX C FROM UPPER PART OF A

        DO 20 i = 1, meq
            CALL dcopy_ (n, a(i, 1), la, w(ic - 1 + i), meq)
        20     CONTINUE

        !  RECOVER VECTOR D FROM UPPER PART OF B

        CALL dcopy_ (meq, b(1), 1, w(id), 1)
        CALL dscal_sl (meq, - one, w(id), 1)

    ENDIF

    ig = id + meq

    !  RECOVER MATRIX G FROM LOWER PART OF A
    !  The matrix G(mineq+2*n,m1) is stored at w(ig)
    !  Not all rows will be filled if some of the upper/lower
    !  bounds are unbounded.

    IF (mineq .GT. 0) THEN

        DO 30 i = 1, mineq
            CALL dcopy_ (n, a(meq + i, 1), la, w(ig - 1 + i), m1)
        30     CONTINUE

    ENDIF

    ih = ig + m1 * n
    iw = ih + mineq + 2 * n

    IF (mineq .GT. 0) THEN

        !  RECOVER H FROM LOWER PART OF B
        !  The vector H(mineq+2*n) is stored at w(ih)

        CALL dcopy_ (mineq, b(meq + 1), 1, w(ih), 1)
        CALL dscal_sl (mineq, - one, w(ih), 1)

    ENDIF

    !  AUGMENT MATRIX G BY +I AND -I, AND,
    !  AUGMENT VECTOR H BY XL AND XU
    !  NaN value indicates no bound

    ip = ig + mineq
    il = ih + mineq
    nancnt = 0

    DO 40 i = 1, n
        if (xl(i).eq.xl(i)) then
            w(il) = xl(i)
            do 41 j = 1, n
                w(ip + m1 * (j - 1)) = 0
            41         continue
            w(ip + m1 * (i - 1)) = 1
            ip = ip + 1
            il = il + 1
        else
            nancnt = nancnt + 1
        end if
    40 CONTINUE

    DO 50 i = 1, n
        if (xu(i).eq.xu(i)) then
            w(il) = -xu(i)
            do 51 j = 1, n
                w(ip + m1 * (j - 1)) = 0
            51         continue
            w(ip + m1 * (i - 1)) = -1
            ip = ip + 1
            il = il + 1
        else
            nancnt = nancnt + 1
        end if
    50   CONTINUE

    CALL lsei (w(ic), w(id), w(ie), w(IF), w(ig), w(ih), MAX(1, meq), &
            meq, n, n, m1, m1 - nancnt, n, x, xnorm, w(iw), jw, mode)

    IF (mode .EQ. 1) THEN

        !   restore Lagrange multipliers (only for user-defined variables)

        CALL dcopy_ (m, w(iw), 1, y(1), 1)

        !   set rest of the multipliers to nan (they are not used)

        IF (n3 .GT. 0) THEN
            y(m + 1) = 0
            y(m + 1) = 0 / y(m + 1)
            do 60 i = m + 2, m + n3 + n3
                y(i) = y(m + 1)
            60          continue
        ENDIF

    ENDIF
    call bound(n, x, xl, xu)

    !   END OF SUBROUTINE LSQ

END


SUBROUTINE lsei(c, d, e, f, g, h, lc, mc, LE, me, lg, mg, n, x, xnrm, w, jw, mode)

    !     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF
    !     EQUALITY & INEQUALITY CONSTRAINED LEAST SQUARES PROBLEM LSEI :

    !                MIN ||E*X - F||
    !                 X

    !                S.T.  C*X  = D,
    !                      G*X >= H.

    !     USING QR DECOMPOSITION & ORTHOGONAL BASIS OF NULLSPACE OF C
    !     CHAPTER 23.6 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS.

    !     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM
    !     ARE NECESSARY
    !     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N)
    !     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  )
    !     DIM(C) :   FORMAL (LC,N),    ACTUAL (MC,N)
    !     DIM(D) :   FORMAL (LC  ),    ACTUAL (MC  )
    !     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N)
    !     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  )
    !     DIM(X) :   FORMAL (N   ),    ACTUAL (N   )
    !     DIM(W) :   2*MC+ME+(ME+MG)*(N-MC)  for LSEI
    !              +(N-MC+1)*(MG+2)+2*MG     for LSI
    !     DIM(JW):   MAX(MG,L)
    !     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS C, D, E, F, G, AND H.
    !     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE.
    !     X     STORES THE SOLUTION VECTOR
    !     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM
    !     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST
    !           MC+MG ELEMENTS
    !     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
    !               5: MATRIX E IS NOT OF FULL RANK
    !               6: MATRIX C IS NOT OF FULL RANK
    !               7: RANK DEFECT IN HFTI

    !     18.5.1981, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN
    !     20.3.1987, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN

    INTEGER          jw(*), i, ie, IF, ig, iw, j, k, krank, l, lc, LE, lg, &
            mc, mc1, me, mg, mode, n
    DOUBLE PRECISION c(lc, n), e(LE, n), g(lg, n), d(lc), f(LE), h(lg), x(n), &
            w(*), t, ddot_sl, xnrm, dnrm2_, epmach, ZERO
    DATA             epmach/2.22d-16/, ZERO/0.0d+00/

    mode = 2
    IF(mc.GT.n)                      GOTO 75
    l = n - mc
    mc1 = mc + 1
    iw = (l + 1) * (mg + 2) + 2 * mg + mc
    ie = iw + mc + 1
    IF = ie + me * l
    ig = IF + me

    !  TRIANGULARIZE C AND APPLY FACTORS TO E AND G

    DO 10 i = 1, mc
        j = MIN(i + 1, lc)
        CALL h12(1, i, i + 1, n, c(i, 1), lc, w(iw + i), c(j, 1), lc, 1, mc - i)
        CALL h12(2, i, i + 1, n, c(i, 1), lc, w(iw + i), e, LE, 1, me)
    10     CALL h12(2, i, i + 1, n, c(i, 1), lc, w(iw + i), g, lg, 1, mg)

    !  SOLVE C*X=D AND MODIFY F

    mode = 6
    DO 15 i = 1, mc
        IF(ABS(c(i, i)).LT.epmach)    GOTO 75
        x(i) = (d(i) - ddot_sl(i - 1, c(i, 1), lc, x, 1)) / c(i, i)
    15 CONTINUE
    mode = 1
    w(mc1) = ZERO
    CALL dcopy_ (mg - mc, w(mc1), 0, w(mc1), 1)

    IF(mc.EQ.n)                      GOTO 50

    DO 20 i = 1, me
    20     w(IF - 1 + i) = f(i) - ddot_sl(mc, e(i, 1), LE, x, 1)

    !  STORE TRANSFORMED E & G

    DO 25 i = 1, me
    25     CALL dcopy_(l, e(i, mc1), LE, w(ie - 1 + i), me)
    DO 30 i = 1, mg
    30     CALL dcopy_(l, g(i, mc1), lg, w(ig - 1 + i), mg)

    IF(mg.GT.0)                      GOTO 40

    !  SOLVE LS WITHOUT INEQUALITY CONSTRAINTS

    mode = 7
    k = MAX(LE, n)
    t = SQRT(epmach)
    CALL hfti (w(ie), me, me, l, w(IF), k, 1, t, krank, xnrm, w, w(l + 1), jw)
    CALL dcopy_(l, w(IF), 1, x(mc1), 1)
    IF(krank.NE.l)                   GOTO 75
    mode = 1
    GOTO 50
    !  MODIFY H AND SOLVE INEQUALITY CONSTRAINED LS PROBLEM

    40 DO 45 i = 1, mg
    45     h(i) = h(i) - ddot_sl(mc, g(i, 1), lg, x, 1)
    CALL lsi&
            (w(ie), w(IF), w(ig), h, me, me, mg, mg, l, x(mc1), xnrm, w(mc1), jw, mode)
    IF(mc.EQ.0)                      GOTO 75
    t = dnrm2_(mc, x, 1)
    xnrm = SQRT(xnrm * xnrm + t * t)
    IF(mode.NE.1)                    GOTO 75

    !  SOLUTION OF ORIGINAL PROBLEM AND LAGRANGE MULTIPLIERS

    50 DO 55 i = 1, me
    55     f(i) = ddot_sl(n, e(i, 1), LE, x, 1) - f(i)
    DO 60 i = 1, mc
    60     d(i) = ddot_sl(me, e(1, i), 1, f, 1) - ddot_sl(mg, g(1, i), 1, w(mc1), 1)

    DO 65 i = mc, 1, -1
    65     CALL h12(2, i, i + 1, n, c(i, 1), lc, w(iw + i), x, 1, 1, 1)

    DO 70 i = mc, 1, -1
        j = MIN(i + 1, lc)
        w(i) = (d(i) - ddot_sl(mc - i, c(j, i), 1, w(j), 1)) / c(i, i)
    70 CONTINUE

    !  END OF SUBROUTINE LSEI

75                                  END


SUBROUTINE lsi(e, f, g, h, LE, me, lg, mg, n, x, xnorm, w, jw, mode)

    !     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF
    !     INEQUALITY CONSTRAINED LINEAR LEAST SQUARES PROBLEM:

    !                    MIN ||E*X-F||
    !                     X

    !                    S.T.  G*X >= H

    !     THE ALGORITHM IS BASED ON QR DECOMPOSITION AS DESCRIBED IN
    !     CHAPTER 23.5 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS

    !     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM
    !     ARE NECESSARY
    !     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N)
    !     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  )
    !     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N)
    !     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  )
    !     DIM(X) :   N
    !     DIM(W) :   (N+1)*(MG+2) + 2*MG
    !     DIM(JW):   LG
    !     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS E, F, G, AND H.
    !     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE.
    !     X     STORES THE SOLUTION VECTOR
    !     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM
    !     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST
    !           MG ELEMENTS
    !     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
    !               5: MATRIX E IS NOT OF FULL RANK

    !     03.01.1980, DIETER KRAFT: CODED
    !     20.03.1987, DIETER KRAFT: REVISED TO FORTRAN 77

    INTEGER          i, j, LE, lg, me, mg, mode, n, jw(lg)
    DOUBLE PRECISION e(LE, n), f(LE), g(lg, n), h(lg), x(n), w(*), &
            ddot_sl, xnorm, dnrm2_, epmach, t, one
    DATA             epmach/2.22d-16/, one/1.0d+00/

    !  QR-FACTORS OF E AND APPLICATION TO F

    DO 10 i = 1, n
        j = MIN(i + 1, n)
        CALL h12(1, i, i + 1, me, e(1, i), 1, t, e(1, j), 1, LE, n - i)
    10 CALL h12(2, i, i + 1, me, e(1, i), 1, t, f, 1, 1, 1)

    !  TRANSFORM G AND H TO GET LEAST DISTANCE PROBLEM

    mode = 5
    DO 30 i = 1, mg
        DO 20 j = 1, n
            IF (.NOT.(ABS(e(j, j)).GE.epmach)) GOTO 50
        20         g(i, j) = (g(i, j) - ddot_sl(j - 1, g(i, 1), lg, e(1, j), 1)) / e(j, j)
    30     h(i) = h(i) - ddot_sl(n, g(i, 1), lg, f, 1)

    !  SOLVE LEAST DISTANCE PROBLEM

    CALL ldp(g, lg, mg, n, h, x, xnorm, w, jw, mode)
    IF (mode.NE.1)                     GOTO 50

    !  SOLUTION OF ORIGINAL PROBLEM

    CALL daxpy_sl(n, one, f, 1, x, 1)
    DO 40 i = n, 1, -1
        j = MIN(i + 1, n)
    40     x(i) = (x(i) - ddot_sl(n - i, e(i, j), LE, x(j), 1)) / e(i, i)
    j = MIN(n + 1, me)
    t = dnrm2_(me - n, f(j), 1)
    xnorm = SQRT(xnorm * xnorm + t * t)

    !  END OF SUBROUTINE LSI

50                                    END

SUBROUTINE ldp(g, mg, m, n, h, x, xnorm, w, INDEX, mode)

    !                     T
    !     MINIMIZE   1/2 X X    SUBJECT TO   G * X >= H.

    !       C.L. LAWSON, R.J. HANSON: 'SOLVING LEAST SQUARES PROBLEMS'
    !       PRENTICE HALL, ENGLEWOOD CLIFFS, NEW JERSEY, 1974.

    !     PARAMETER DESCRIPTION:

    !     G(),MG,M,N   ON ENTRY G() STORES THE M BY N MATRIX OF
    !                  LINEAR INEQUALITY CONSTRAINTS. G() HAS FIRST
    !                  DIMENSIONING PARAMETER MG
    !     H()          ON ENTRY H() STORES THE M VECTOR H REPRESENTING
    !                  THE RIGHT SIDE OF THE INEQUALITY SYSTEM

    !     REMARK: G(),H() WILL NOT BE CHANGED DURING CALCULATIONS BY LDP

    !     X()          ON ENTRY X() NEED NOT BE INITIALIZED.
    !                  ON EXIT X() STORES THE SOLUTION VECTOR X IF MODE=1.
    !     XNORM        ON EXIT XNORM STORES THE EUCLIDIAN NORM OF THE
    !                  SOLUTION VECTOR IF COMPUTATION IS SUCCESSFUL
    !     W()          W IS A ONE DIMENSIONAL WORKING SPACE, THE LENGTH
    !                  OF WHICH SHOULD BE AT LEAST (M+2)*(N+1) + 2*M
    !                  ON EXIT W() STORES THE LAGRANGE MULTIPLIERS
    !                  ASSOCIATED WITH THE CONSTRAINTS
    !                  AT THE SOLUTION OF PROBLEM LDP
    !     INDEX()      INDEX() IS A ONE DIMENSIONAL INTEGER WORKING SPACE
    !                  OF LENGTH AT LEAST M
    !     MODE         MODE IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING
    !                  MEANINGS:
    !          MODE=1: SUCCESSFUL COMPUTATION
    !               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N.LE.0)
    !               3: ITERATION COUNT EXCEEDED BY NNLS
    !               4: INEQUALITY CONSTRAINTS INCOMPATIBLE

    DOUBLE PRECISION g, h, x, xnorm, w, u, v, &
            ZERO, one, fac, rnorm, dnrm2_, ddot_sl, diff
    INTEGER          INDEX, i, IF, iw, iwdual, iy, iz, j, m, mg, mode, n, n1
    DIMENSION        g(mg, n), h(m), x(n), w(*), INDEX(m)
    diff(u, v) = u - v
    DATA             ZERO, one/0.0d0, 1.0d0/

    mode = 2
    IF(n.LE.0)                    GOTO 50

    !  STATE DUAL PROBLEM

    mode = 1
    x(1) = ZERO
    CALL dcopy_(n, x(1), 0, x, 1)
    xnorm = ZERO
    IF(m.EQ.0)                    GOTO 50
    iw = 0
    DO 20 j = 1, m
        DO 10 i = 1, n
            iw = iw + 1
        10         w(iw) = g(j, i)
        iw = iw + 1
    20     w(iw) = h(j)
    IF = iw + 1
    DO 30 i = 1, n
        iw = iw + 1
    30     w(iw) = ZERO
    w(iw + 1) = one
    n1 = n + 1
    iz = iw + 2
    iy = iz + n1
    iwdual = iy + m

    !  SOLVE DUAL PROBLEM

    CALL nnls (w, n1, n1, m, w(IF), w(iy), rnorm, w(iwdual), w(iz), INDEX, mode)

    IF(mode.NE.1)                 GOTO 50
    mode = 4
    IF(rnorm.LE.ZERO)             GOTO 50

    !  COMPUTE SOLUTION OF PRIMAL PROBLEM

    fac = one - ddot_sl(m, h, 1, w(iy), 1)
    IF(.NOT.(diff(one + fac, one).GT.ZERO)) GOTO 50
    mode = 1
    fac = one / fac
    DO 40 j = 1, n
    40     x(j) = fac * ddot_sl(m, g(1, j), 1, w(iy), 1)
    xnorm = dnrm2_(n, x, 1)

    !  COMPUTE LAGRANGE MULTIPLIERS FOR PRIMAL PROBLEM

    w(1) = ZERO
    CALL dcopy_(m, w(1), 0, w, 1)
    CALL daxpy_sl(m, fac, w(iy), 1, w, 1)

    !  END OF SUBROUTINE LDP

50                               END


SUBROUTINE nnls (a, mda, m, n, b, x, rnorm, w, z, INDEX, mode)

    !     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY:
    !     'SOLVING LEAST SQUARES PROBLEMS'. PRENTICE-HALL.1974

    !      **********   NONNEGATIVE LEAST SQUARES   **********

    !     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B, COMPUTE AN
    !     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM

    !                  A*X = B  SUBJECT TO  X >= 0

    !     A(),MDA,M,N
    !            MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE ARRAY,A().
    !            ON ENTRY A()  CONTAINS THE M BY N MATRIX,A.
    !            ON EXIT A() CONTAINS THE PRODUCT Q*A,
    !            WHERE Q IS AN M BY M ORTHOGONAL MATRIX GENERATED
    !            IMPLICITLY BY THIS SUBROUTINE.
    !            EITHER M>=N OR M<N IS PERMISSIBLE.
    !            THERE IS NO RESTRICTION ON THE RANK OF A.
    !     B()    ON ENTRY B() CONTAINS THE M-VECTOR, B.
    !            ON EXIT B() CONTAINS Q*B.
    !     X()    ON ENTRY X() NEED NOT BE INITIALIZED.
    !            ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR.
    !     RNORM  ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE
    !            RESIDUAL VECTOR.
    !     W()    AN N-ARRAY OF WORKING SPACE.
    !            ON EXIT W() WILL CONTAIN THE DUAL SOLUTION VECTOR.
    !            W WILL SATISFY W(I)=0 FOR ALL I IN SET P
    !            AND W(I)<=0 FOR ALL I IN SET Z
    !     Z()    AN M-ARRAY OF WORKING SPACE.
    !     INDEX()AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
    !            ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS
    !            P AND Z AS FOLLOWS:
    !            INDEX(1)    THRU INDEX(NSETP) = SET P.
    !            INDEX(IZ1)  THRU INDEX (IZ2)  = SET Z.
    !            IZ1=NSETP + 1 = NPP1, IZ2=N.
    !     MODE   THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANING:
    !            1    THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
    !            2    THE DIMENSIONS OF THE PROBLEM ARE WRONG,
    !                 EITHER M <= 0 OR N <= 0.
    !            3    ITERATION COUNT EXCEEDED, MORE THAN 3*N ITERATIONS.

    INTEGER          i, ii, ip, iter, itmax, iz, izmax, iz1, iz2, j, jj, jz, &
            k, l, m, mda, mode, n, npp1, nsetp, INDEX(n)

    DOUBLE PRECISION a(mda, n), b(m), x(n), w(n), z(m), asave, diff, &
            factor, ddot_sl, ZERO, one, wmax, alpha, &
            c, s, t, u, v, up, rnorm, unorm, dnrm2_

    diff(u, v) = u - v

    DATA             ZERO, one, factor/0.0d0, 1.0d0, 1.0d-2/

    !     revised          Dieter Kraft, March 1983

    mode = 2
    IF(m.LE.0.OR.n.LE.0)            GOTO 290
    mode = 1
    iter = 0
    itmax = 3 * n

    ! STEP ONE (INITIALIZE)

    DO 100 i = 1, n
    100    INDEX(i) = i
    iz1 = 1
    iz2 = n
    nsetp = 0
    npp1 = 1
    x(1) = ZERO
    CALL dcopy_(n, x(1), 0, x, 1)

    ! STEP TWO (COMPUTE DUAL VARIABLES)
    ! .....ENTRY LOOP A

    110 IF(iz1.GT.iz2.OR.nsetp.GE.m)    GOTO 280
    DO 120 iz = iz1, iz2
        j = INDEX(iz)
    120    w(j) = ddot_sl(m - nsetp, a(npp1, j), 1, b(npp1), 1)

    ! STEP THREE (TEST DUAL VARIABLES)

    130 wmax = ZERO
    DO 140 iz = iz1, iz2
        j = INDEX(iz)
        IF(w(j).LE.wmax)             GOTO 140
        wmax = w(j)
        izmax = iz
    140 CONTINUE

    ! .....EXIT LOOP A

    IF(wmax.LE.ZERO)                GOTO 280
    iz = izmax
    j = INDEX(iz)

    ! STEP FOUR (TEST INDEX J FOR LINEAR DEPENDENCY)

    asave = a(npp1, j)
    CALL h12(1, npp1, npp1 + 1, m, a(1, j), 1, up, z, 1, 1, 0)
    unorm = dnrm2_(nsetp, a(1, j), 1)
    t = factor * ABS(a(npp1, j))
    IF(diff(unorm + t, unorm).LE.ZERO) GOTO 150
    CALL dcopy_(m, b, 1, z, 1)
    CALL h12(2, npp1, npp1 + 1, m, a(1, j), 1, up, z, 1, 1, 1)
    IF(z(npp1) / a(npp1, j).GT.ZERO)   GOTO 160
    150 a(npp1, j) = asave
    w(j) = ZERO
    GOTO 130
    ! STEP FIVE (ADD COLUMN)

    160 CALL dcopy_(m, z, 1, b, 1)
    INDEX(iz) = INDEX(iz1)
    INDEX(iz1) = j
    iz1 = iz1 + 1
    nsetp = npp1
    npp1 = npp1 + 1
    DO 170 jz = iz1, iz2
        jj = INDEX(jz)
    170    CALL h12(2, nsetp, npp1, m, a(1, j), 1, up, a(1, jj), 1, mda, 1)
    k = MIN(npp1, mda)
    w(j) = ZERO
    CALL dcopy_(m - nsetp, w(j), 0, a(k, j), 1)

    ! STEP SIX (SOLVE LEAST SQUARES SUB-PROBLEM)
    ! .....ENTRY LOOP B

    180 DO 200 ip = nsetp, 1, -1
        IF(ip.EQ.nsetp)              GOTO 190
        CALL daxpy_sl(ip, -z(ip + 1), a(1, jj), 1, z, 1)
        190    jj = INDEX(ip)
    200    z(ip) = z(ip) / a(ip, jj)
    iter = iter + 1
    IF(iter.LE.itmax)               GOTO 220
    210 mode = 3
    GOTO 280
    ! STEP SEVEN TO TEN (STEP LENGTH ALGORITHM)

    220 alpha = one
    jj = 0
    DO 230 ip = 1, nsetp
        IF(z(ip).GT.ZERO)            GOTO 230
        l = INDEX(ip)
        t = -x(l) / (z(ip) - x(l))
        IF(alpha.LT.t)               GOTO 230
        alpha = t
        jj = ip
    230 CONTINUE
    DO 240 ip = 1, nsetp
        l = INDEX(ip)
    240    x(l) = (one - alpha) * x(l) + alpha * z(ip)

    ! .....EXIT LOOP B

    IF(jj.EQ.0)                     GOTO 110

    ! STEP ELEVEN (DELETE COLUMN)

    i = INDEX(jj)
    250 x(i) = ZERO
    jj = jj + 1
    DO 260 j = jj, nsetp
        ii = INDEX(j)
        INDEX(j - 1) = ii
        CALL dsrotg(a(j - 1, ii), a(j, ii), c, s)
        t = a(j - 1, ii)
        CALL dsrot(n, a(j - 1, 1), mda, a(j, 1), mda, c, s)
        a(j - 1, ii) = t
        a(j, ii) = ZERO
    260    CALL dsrot(1, b(j - 1), 1, b(j), 1, c, s)
    npp1 = nsetp
    nsetp = nsetp - 1
    iz1 = iz1 - 1
    INDEX(iz1) = i
    IF(nsetp.LE.0)                  GOTO 210
    DO 270 jj = 1, nsetp
        i = INDEX(jj)
        IF(x(i).LE.ZERO)             GOTO 250
    270 CONTINUE
    CALL dcopy_(m, b, 1, z, 1)
    GOTO 180
    ! STEP TWELVE (SOLUTION)

    280 k = MIN(npp1, m)
    rnorm = dnrm2_(m - nsetp, b(k), 1)
    IF(npp1.GT.m) THEN
        w(1) = ZERO
        CALL dcopy_(n, w(1), 0, w, 1)
    ENDIF

    ! END OF SUBROUTINE NNLS

290                                 END

SUBROUTINE hfti(a, mda, m, n, b, mdb, nb, tau, krank, rnorm, h, g, ip)

    !     RANK-DEFICIENT LEAST SQUARES ALGORITHM AS DESCRIBED IN:
    !     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
    !     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974

    !     A(*,*),MDA,M,N   THE ARRAY A INITIALLY CONTAINS THE M x N MATRIX A
    !                      OF THE LEAST SQUARES PROBLEM AX = B.
    !                      THE FIRST DIMENSIONING PARAMETER MDA MUST SATISFY
    !                      MDA >= M. EITHER M >= N OR M < N IS PERMITTED.
    !                      THERE IS NO RESTRICTION ON THE RANK OF A.
    !                      THE MATRIX A WILL BE MODIFIED BY THE SUBROUTINE.
    !     B(*,*),MDB,NB    IF NB = 0 THE SUBROUTINE WILL MAKE NO REFERENCE
    !                      TO THE ARRAY B. IF NB > 0 THE ARRAY B() MUST
    !                      INITIALLY CONTAIN THE M x NB MATRIX B  OF THE
    !                      THE LEAST SQUARES PROBLEM AX = B AND ON RETURN
    !                      THE ARRAY B() WILL CONTAIN THE N x NB SOLUTION X.
    !                      IF NB>1 THE ARRAY B() MUST BE DOUBLE SUBSCRIPTED
    !                      WITH FIRST DIMENSIONING PARAMETER MDB>=MAX(M,N),
    !                      IF NB=1 THE ARRAY B() MAY BE EITHER SINGLE OR
    !                      DOUBLE SUBSCRIPTED.
    !     TAU              ABSOLUTE TOLERANCE PARAMETER FOR PSEUDORANK
    !                      DETERMINATION, PROVIDED BY THE USER.
    !     KRANK            PSEUDORANK OF A, SET BY THE SUBROUTINE.
    !     RNORM            ON EXIT, RNORM(J) WILL CONTAIN THE EUCLIDIAN
    !                      NORM OF THE RESIDUAL VECTOR FOR THE PROBLEM
    !                      DEFINED BY THE J-TH COLUMN VECTOR OF THE ARRAY B.
    !     H(), G()         ARRAYS OF WORKING SPACE OF LENGTH >= N.
    !     IP()             INTEGER ARRAY OF WORKING SPACE OF LENGTH >= N
    !                      RECORDING PERMUTATION INDICES OF COLUMN VECTORS

    INTEGER          i, j, jb, k, kp1, krank, l, ldiag, lmax, m, &
            mda, mdb, n, nb, ip(n)
    DOUBLE PRECISION a(mda, n), b(mdb, nb), h(n), g(n), rnorm(nb), factor, &
            tau, ZERO, hmax, diff, tmp, ddot_sl, dnrm2_, u, v
    diff(u, v) = u - v
    DATA             ZERO/0.0d0/, factor/1.0d-3/

    k = 0
    ldiag = MIN(m, n)
    IF(ldiag.LE.0)                  GOTO 270

    !   COMPUTE LMAX

    DO 80 j = 1, ldiag
        IF(j.EQ.1)                  GOTO 20
        lmax = j
        DO 10 l = j, n
            h(l) = h(l) - a(j - 1, l)**2
        10         IF(h(l).GT.h(lmax)) lmax = l
        IF(diff(hmax + factor * h(lmax), hmax).GT.ZERO)&
                GOTO 50
        20     lmax = j
        DO 40 l = j, n
            h(l) = ZERO
            DO 30 i = j, m
            30             h(l) = h(l) + a(i, l)**2
        40         IF(h(l).GT.h(lmax)) lmax = l
        hmax = h(lmax)

        !   COLUMN INTERCHANGES IF NEEDED

        50     ip(j) = lmax
        IF(ip(j).EQ.j)              GOTO 70
        DO 60 i = 1, m
            tmp = a(i, j)
            a(i, j) = a(i, lmax)
        60         a(i, lmax) = tmp
        h(lmax) = h(j)

        !   J-TH TRANSFORMATION AND APPLICATION TO A AND B

        70     i = MIN(j + 1, n)
        CALL h12(1, j, j + 1, m, a(1, j), 1, h(j), a(1, i), 1, mda, n - j)
    80     CALL h12(2, j, j + 1, m, a(1, j), 1, h(j), b, 1, mdb, nb)

    !   DETERMINE PSEUDORANK

    DO 90 j = 1, ldiag
    90     IF(ABS(a(j, j)).LE.tau)      GOTO 100
    k = ldiag
    GOTO 110
    100 k = j - 1
    110 kp1 = k + 1

    !   NORM OF RESIDUALS

    DO 130 jb = 1, nb
    130     rnorm(jb) = dnrm2_(m - k, b(kp1, jb), 1)
    IF(k.GT.0)                      GOTO 160
    DO 150 jb = 1, nb
        DO 150 i = 1, n
        150         b(i, jb) = ZERO
    GOTO 270
    160 IF(k.EQ.n)                      GOTO 180

    !   HOUSEHOLDER DECOMPOSITION OF FIRST K ROWS

    DO 170 i = k, 1, -1
    170     CALL h12(1, i, kp1, n, a(i, 1), mda, g(i), a, mda, 1, i - 1)
    180 DO 250 jb = 1, nb

        !   SOLVE K*K TRIANGULAR SYSTEM

        DO 210 i = k, 1, -1
            j = MIN(i + 1, n)
        210         b(i, jb) = (b(i, jb) - ddot_sl(k - i, a(i, j), mda, b(j, jb), 1)) / a(i, i)

        !   COMPLETE SOLUTION VECTOR

        IF(k.EQ.n)                  GOTO 240
        DO 220 j = kp1, n
        220         b(j, jb) = ZERO
        DO 230 i = 1, k
        230         CALL h12(2, i, kp1, n, a(i, 1), mda, g(i), b(1, jb), 1, mdb, 1)

        !   REORDER SOLUTION ACCORDING TO PREVIOUS COLUMN INTERCHANGES

        240     DO 250 j = ldiag, 1, -1
            IF(ip(j).EQ.j)          GOTO 250
            l = ip(j)
            tmp = b(l, jb)
            b(l, jb) = b(j, jb)
            b(j, jb) = tmp
        250 CONTINUE
        270 krank = k
END

SUBROUTINE h12 (mode, lpivot, l1, m, u, iue, up, c, ice, icv, ncv)

    !     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12
    !     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974

    !     CONSTRUCTION AND/OR APPLICATION OF A SINGLE
    !     HOUSEHOLDER TRANSFORMATION  Q = I + U*(U**T)/B

    !     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 .
    !     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
    !     L1,M   IF L1 <= M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
    !            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.
    !            IF L1 > M THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
    !     U(),IUE,UP
    !            ON ENTRY TO H1 U() STORES THE PIVOT VECTOR.
    !            IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS.
    !            ON EXIT FROM H1 U() AND UP STORE QUANTITIES DEFINING
    !            THE VECTOR U OF THE HOUSEHOLDER TRANSFORMATION.
    !            ON ENTRY TO H2 U() AND UP
    !            SHOULD STORE QUANTITIES PREVIOUSLY COMPUTED BY H1.
    !            THESE WILL NOT BE MODIFIED BY H2.
    !     C()    ON ENTRY TO H1 OR H2 C() STORES A MATRIX WHICH WILL BE
    !            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER
    !            TRANSFORMATION IS TO BE APPLIED.
    !            ON EXIT C() STORES THE SET OF TRANSFORMED VECTORS.
    !     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
    !     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
    !     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED.
    !            IF NCV <= 0 NO OPERATIONS WILL BE DONE ON C().

    INTEGER          incr, ice, icv, iue, lpivot, l1, mode, ncv
    INTEGER          i, i2, i3, i4, j, m
    DOUBLE PRECISION u, up, c, cl, clinv, b, sm, one, ZERO
    DIMENSION        u(iue, *), c(*)
    DATA             one/1.0d+00/, ZERO/0.0d+00/

    IF (0.GE.lpivot.OR.lpivot.GE.l1.OR.l1.GT.m) GOTO 80
    cl = ABS(u(1, lpivot))
    IF (mode.EQ.2)                              GOTO 30

    !     ****** CONSTRUCT THE TRANSFORMATION ******

    DO 10 j = l1, m
        sm = ABS(u(1, j))
    10     cl = MAX(sm, cl)
    IF (cl.LE.ZERO)                             GOTO 80
    clinv = one / cl
    sm = (u(1, lpivot) * clinv)**2
    DO 20 j = l1, m
    20     sm = sm + (u(1, j) * clinv)**2
    cl = cl * SQRT(sm)
    IF (u(1, lpivot).GT.ZERO) cl = -cl
    up = u(1, lpivot) - cl
    u(1, lpivot) = cl
    GOTO 40
    !     ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C ******

    30 IF (cl.LE.ZERO)                             GOTO 80
    40 IF (ncv.LE.0)                               GOTO 80
    b = up * u(1, lpivot)
    IF (b.GE.ZERO)                              GOTO 80
    b = one / b
    i2 = 1 - icv + ice * (lpivot - 1)
    incr = ice * (l1 - lpivot)
    DO 70 j = 1, ncv
        i2 = i2 + icv
        i3 = i2 + incr
        i4 = i3
        sm = c(i2) * up
        DO 50 i = l1, m
            sm = sm + c(i3) * u(1, i)
        50         i3 = i3 + ice
        IF (sm.EQ.ZERO)                         GOTO 70
        sm = sm * b
        c(i2) = c(i2) + sm * up
        DO 60 i = l1, m
            c(i4) = c(i4) + sm * u(1, i)
        60         i4 = i4 + ice
    70     CONTINUE
80                                             END

SUBROUTINE ldl (n, a, z, sigma, w)
    !   LDL     LDL' - RANK-ONE - UPDATE

    !   PURPOSE:
    !           UPDATES THE LDL' FACTORS OF MATRIX A BY RANK-ONE MATRIX
    !           SIGMA*Z*Z'

    !   INPUT ARGUMENTS: (* MEANS PARAMETERS ARE CHANGED DURING EXECUTION)
    !     N     : ORDER OF THE COEFFICIENT MATRIX A
    !   * A     : POSITIVE DEFINITE MATRIX OF DIMENSION N;
    !             ONLY THE LOWER TRIANGLE IS USED AND IS STORED COLUMN BY
    !             COLUMN AS ONE DIMENSIONAL ARRAY OF DIMENSION N*(N+1)/2.
    !   * Z     : VECTOR OF DIMENSION N OF UPDATING ELEMENTS
    !     SIGMA : SCALAR FACTOR BY WHICH THE MODIFYING DYADE Z*Z' IS
    !             MULTIPLIED

    !   OUTPUT ARGUMENTS:
    !     A     : UPDATED LDL' FACTORS

    !   WORKING ARRAY:
    !     W     : VECTOR OP DIMENSION N (USED ONLY IF SIGMA .LT. ZERO)

    !   METHOD:
    !     THAT OF FLETCHER AND POWELL AS DESCRIBED IN :
    !     FLETCHER,R.,(1974) ON THE MODIFICATION OF LDL' FACTORIZATION.
    !     POWELL,M.J.D.      MATH.COMPUTATION 28, 1067-1078.

    !   IMPLEMENTED BY:
    !     KRAFT,D., DFVLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME
    !               D-8031  OBERPFAFFENHOFEN

    !   STATUS: 15. JANUARY 1980

    !   SUBROUTINES REQUIRED: NONE

    INTEGER          i, ij, j, n
    DOUBLE PRECISION a(*), t, v, w(*), z(*), u, tp, one, beta, four, &
            ZERO, alpha, delta, gamma, sigma, epmach
    DATA ZERO, one, four, epmach /0.0d0, 1.0d0, 4.0d0, 2.22d-16/

    IF(sigma.EQ.ZERO)               GOTO 280
    ij = 1
    t = one / sigma
    IF(sigma.GT.ZERO)               GOTO 220
    ! PREPARE NEGATIVE UPDATE
    DO 150 i = 1, n
    150     w(i) = z(i)
    DO 170 i = 1, n
        v = w(i)
        t = t + v * v / a(ij)
        DO 160 j = i + 1, n
            ij = ij + 1
        160         w(j) = w(j) - v * a(ij)
    170     ij = ij + 1
    IF(t.GE.ZERO) t = epmach / sigma
    DO 210 i = 1, n
        j = n + 1 - i
        ij = ij - i
        u = w(j)
        w(j) = t
    210     t = t - u * u / a(ij)
    220 CONTINUE
    ! HERE UPDATING BEGINS
    DO 270 i = 1, n
        v = z(i)
        delta = v / a(ij)
        IF(sigma.LT.ZERO) tp = w(i)
        IF(sigma.GT.ZERO) tp = t + delta * v
        alpha = tp / t
        a(ij) = alpha * a(ij)
        IF(i.EQ.n)                  GOTO 280
        beta = delta / tp
        IF(alpha.GT.four)           GOTO 240
        DO 230 j = i + 1, n
            ij = ij + 1
            z(j) = z(j) - v * a(ij)
        230         a(ij) = a(ij) + beta * z(j)
        GOTO 260
        240     gamma = t / tp
        DO 250 j = i + 1, n
            ij = ij + 1
            u = a(ij)
            a(ij) = gamma * u + beta * z(j)
        250         z(j) = z(j) - v * u
        260     ij = ij + 1
    270     t = tp
    280 RETURN
    ! END OF LDL
END

DOUBLE PRECISION FUNCTION linmin (mode, ax, bx, f, tol)
    !   LINMIN  LINESEARCH WITHOUT DERIVATIVES

    !   PURPOSE:

    !  TO FIND THE ARGUMENT LINMIN WHERE THE FUNCTION F TAKES IT'S MINIMUM
    !  ON THE INTERVAL AX, BX.
    !  COMBINATION OF GOLDEN SECTION AND SUCCESSIVE QUADRATIC INTERPOLATION.

    !   INPUT ARGUMENTS: (* MEANS PARAMETERS ARE CHANGED DURING EXECUTION)

    ! *MODE   SEE OUTPUT ARGUMENTS
    !  AX     LEFT ENDPOINT OF INITIAL INTERVAL
    !  BX     RIGHT ENDPOINT OF INITIAL INTERVAL
    !  F      FUNCTION VALUE AT LINMIN WHICH IS TO BE BROUGHT IN BY
    !         REVERSE COMMUNICATION CONTROLLED BY MODE
    !  TOL    DESIRED LENGTH OF INTERVAL OF UNCERTAINTY OF FINAL RESULT

    !   OUTPUT ARGUMENTS:

    !  LINMIN ABSCISSA APPROXIMATING THE POINT WHERE F ATTAINS A MINIMUM
    !  MODE   CONTROLS REVERSE COMMUNICATION
    !         MUST BE SET TO 0 INITIALLY, RETURNS WITH INTERMEDIATE
    !         VALUES 1 AND 2 WHICH MUST NOT BE CHANGED BY THE USER,
    !         ENDS WITH CONVERGENCE WITH VALUE 3.

    !   WORKING ARRAY:

    !  NONE

    !   METHOD:

    !  THIS FUNCTION SUBPROGRAM IS A SLIGHTLY MODIFIED VERSION OF THE
    !  ALGOL 60 PROCEDURE LOCALMIN GIVEN IN
    !  R.P. BRENT: ALGORITHMS FOR MINIMIZATION WITHOUT DERIVATIVES,
    !              PRENTICE-HALL (1973).

    !   IMPLEMENTED BY:

    !     KRAFT, D., DFVLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME
    !                D-8031  OBERPFAFFENHOFEN

    !   STATUS: 31. AUGUST  1984

    !   SUBROUTINES REQUIRED: NONE

    INTEGER          mode
    DOUBLE PRECISION f, tol, a, b, c, d, e, p, q, r, u, v, w, x, m, &
            fu, fv, fw, fx, eps, tol1, tol2, ZERO, ax, bx
    DATA             c /0.381966011d0/, eps /1.5d-8/, ZERO /0.0d0/

    !  EPS = SQUARE - ROOT OF MACHINE PRECISION
    !  C = GOLDEN SECTION RATIO = (3-SQRT(5))/2

    GOTO (10, 55), mode

    !  INITIALIZATION

    a = ax
    b = bx
    e = ZERO
    v = a + c * (b - a)
    w = v
    x = w
    linmin = x
    mode = 1
    GOTO 100

    !  MAIN LOOP STARTS HERE

    10 fx = f
    fv = fx
    fw = fv
    20 m = 0.5d0 * (a + b)
    tol1 = eps * ABS(x) + tol
    tol2 = tol1 + tol1

    !  TEST CONVERGENCE

    IF (ABS(x - m) .LE. tol2 - 0.5d0 * (b - a)) GOTO 90
    r = ZERO
    q = r
    p = q
    IF (ABS(e) .LE. tol1) GOTO 30

    !  FIT PARABOLA

    r = (x - w) * (fx - fv)
    q = (x - v) * (fx - fw)
    p = (x - v) * q - (x - w) * r
    q = q - r
    q = q + q
    IF (q .GT. ZERO) p = -p
    IF (q .LT. ZERO) q = -q
    r = e
    e = d

    !  IS PARABOLA ACCEPTABLE

    30 IF (ABS(p) .GE. 0.5d0 * ABS(q * r) .OR.&
            p .LE. q * (a - x) .OR. p .GE. q * (b - x)) GOTO 40

    !  PARABOLIC INTERPOLATION STEP

    d = p / q

    !  F MUST NOT BE EVALUATED TOO CLOSE TO A OR B

    IF (u - a .LT. tol2) d = SIGN(tol1, m - x)
    IF (b - u .LT. tol2) d = SIGN(tol1, m - x)
    GOTO 50

    !  GOLDEN SECTION STEP

    40 IF (x .GE. m) e = a - x
    IF (x .LT. m) e = b - x
    d = c * e

    !  F MUST NOT BE EVALUATED TOO CLOSE TO X

    50 IF (ABS(d) .LT. tol1) d = SIGN(tol1, d)
    u = x + d
    linmin = u
    mode = 2
    GOTO 100
    55 fu = f

    !  UPDATE A, B, V, W, AND X

    IF (fu .GT. fx) GOTO 60
    IF (u .GE. x) a = x
    IF (u .LT. x) b = x
    v = w
    fv = fw
    w = x
    fw = fx
    x = u
    fx = fu
    GOTO 85
    60 IF (u .LT. x) a = u
    IF (u .GE. x) b = u
    IF (fu .LE. fw .OR. w .EQ. x) GOTO 70
    IF (fu .LE. fv .OR. v .EQ. x .OR. v .EQ. w) GOTO 80
    GOTO 85
    70 v = w
    fv = fw
    w = u
    fw = fu
    GOTO 85
    80 v = u
    fv = fu
    85 GOTO 20

    !  END OF MAIN LOOP

    90 linmin = x
    mode = 3
    100 RETURN

    !  END OF LINMIN

END

!## Following a selection from BLAS Level 1

SUBROUTINE daxpy_sl(n, da, dx, incx, dy, incy)

    !     CONSTANT TIMES A VECTOR PLUS A VECTOR.
    !     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION dx(*), dy(*), da
    INTEGER i, incx, incy, ix, iy, m, mp1, n

    IF(n.LE.0)RETURN
    IF(da.EQ.0.0d0)RETURN
    IF(incx.EQ.1.AND.incy.EQ.1)GO TO 20

    !        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
    !        NOT EQUAL TO 1

    ix = 1
    iy = 1
    IF(incx.LT.0)ix = (-n + 1) * incx + 1
    IF(incy.LT.0)iy = (-n + 1) * incy + 1
    DO 10 i = 1, n
        dy(iy) = dy(iy) + da * dx(ix)
        ix = ix + incx
        iy = iy + incy
    10 CONTINUE
    RETURN

    !        CODE FOR BOTH INCREMENTS EQUAL TO 1

    !        CLEAN-UP LOOP

    20 m = MOD(n, 4)
    IF(m .EQ. 0) GO TO 40
    DO 30 i = 1, m
        dy(i) = dy(i) + da * dx(i)
    30 CONTINUE
    IF(n .LT. 4) RETURN
    40 mp1 = m + 1
    DO 50 i = mp1, n, 4
        dy(i) = dy(i) + da * dx(i)
        dy(i + 1) = dy(i + 1) + da * dx(i + 1)
        dy(i + 2) = dy(i + 2) + da * dx(i + 2)
        dy(i + 3) = dy(i + 3) + da * dx(i + 3)
    50 CONTINUE
    RETURN
END

SUBROUTINE  dcopy_(n, dx, incx, dy, incy)

    !     COPIES A VECTOR, X, TO A VECTOR, Y.
    !     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION dx(*), dy(*)
    INTEGER i, incx, incy, ix, iy, m, mp1, n

    IF(n.LE.0)RETURN
    IF(incx.EQ.1.AND.incy.EQ.1)GO TO 20

    !        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
    !        NOT EQUAL TO 1

    ix = 1
    iy = 1
    IF(incx.LT.0)ix = (-n + 1) * incx + 1
    IF(incy.LT.0)iy = (-n + 1) * incy + 1
    DO 10 i = 1, n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
    10 CONTINUE
    RETURN

    !        CODE FOR BOTH INCREMENTS EQUAL TO 1

    !        CLEAN-UP LOOP

    20 m = MOD(n, 7)
    IF(m .EQ. 0) GO TO 40
    DO 30 i = 1, m
        dy(i) = dx(i)
    30 CONTINUE
    IF(n .LT. 7) RETURN
    40 mp1 = m + 1
    DO 50 i = mp1, n, 7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
    50 CONTINUE
    RETURN
END

DOUBLE PRECISION FUNCTION ddot_sl(n, dx, incx, dy, incy)

    !     FORMS THE DOT PRODUCT OF TWO VECTORS.
    !     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION dx(*), dy(*), dtemp
    INTEGER i, incx, incy, ix, iy, m, mp1, n

    ddot_sl = 0.0d0
    dtemp = 0.0d0
    IF(n.LE.0)RETURN
    IF(incx.EQ.1.AND.incy.EQ.1)GO TO 20

    !        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
    !          NOT EQUAL TO 1

    ix = 1
    iy = 1
    IF(incx.LT.0)ix = (-n + 1) * incx + 1
    IF(incy.LT.0)iy = (-n + 1) * incy + 1
    DO 10 i = 1, n
        dtemp = dtemp + dx(ix) * dy(iy)
        ix = ix + incx
        iy = iy + incy
    10 CONTINUE
    ddot_sl = dtemp
    RETURN

    !        CODE FOR BOTH INCREMENTS EQUAL TO 1

    !        CLEAN-UP LOOP

    20 m = MOD(n, 5)
    IF(m .EQ. 0) GO TO 40
    DO 30 i = 1, m
        dtemp = dtemp + dx(i) * dy(i)
    30 CONTINUE
    IF(n .LT. 5) GO TO 60
    40 mp1 = m + 1
    DO 50 i = mp1, n, 5
        dtemp = dtemp + dx(i) * dy(i) + dx(i + 1) * dy(i + 1) + &
                dx(i + 2) * dy(i + 2) + dx(i + 3) * dy(i + 3) + dx(i + 4) * dy(i + 4)
    50 CONTINUE
    60 ddot_sl = dtemp
    RETURN
END

DOUBLE PRECISION FUNCTION dnrm1(n, x, i, j)
    INTEGER n, i, j, k
    DOUBLE PRECISION snormx, sum, x(n), ZERO, one, scale, temp
    DATA ZERO/0.0d0/, one/1.0d0/

    !      DNRM1 - COMPUTES THE I-NORM OF A VECTOR
    !              BETWEEN THE ITH AND THE JTH ELEMENTS

    !      INPUT -
    !      N       LENGTH OF VECTOR
    !      X       VECTOR OF LENGTH N
    !      I       INITIAL ELEMENT OF VECTOR TO BE USED
    !      J       FINAL ELEMENT TO USE

    !      OUTPUT -
    !      DNRM1   NORM

    snormx = ZERO
    DO 10 k = i, j
    10      snormx = MAX(snormx, ABS(x(k)))
    dnrm1 = snormx
    IF (snormx.EQ.ZERO) RETURN
    scale = snormx
    IF (snormx.GE.one) scale = SQRT(snormx)
    sum = ZERO
    DO 20 k = i, j
        temp = ZERO
        IF (ABS(x(k)) + scale .NE. scale) temp = x(k) / snormx
        IF (one + temp.NE.one) sum = sum + temp * temp
    20      CONTINUE
    sum = SQRT(sum)
    dnrm1 = snormx * sum
    RETURN
END

DOUBLE PRECISION FUNCTION dnrm2_ (n, dx, incx)
    INTEGER          n, i, j, nn, next, incx
    DOUBLE PRECISION dx(*), cutlo, cuthi, hitest, sum, xmax, ZERO, one
    DATA             ZERO, one /0.0d0, 1.0d0/

    !     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
    !     INCREMENT INCX .
    !     IF    N .LE. 0 RETURN WITH RESULT = 0.
    !     IF N .GE. 1 THEN INCX MUST BE .GE. 1

    !           C.L.LAWSON, 1978 JAN 08

    !     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
    !     HOPEFULLY APPLICABLE TO ALL MACHINES.
    !         CUTLO = MAXIMUM OF  SQRT(U/EPS)   OVER ALL KNOWN MACHINES.
    !         CUTHI = MINIMUM OF  SQRT(V)       OVER ALL KNOWN MACHINES.
    !     WHERE
    !         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
    !         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
    !         V   = LARGEST  NO.            (OVERFLOW  LIMIT)

    !     BRIEF OUTLINE OF ALGORITHM..

    !     PHASE 1    SCANS ZERO COMPONENTS.
    !     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
    !     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
    !     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
    !     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.

    !     VALUES FOR CUTLO AND CUTHI..
    !     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
    !     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
    !     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
    !                   UNIVAC AND DEC AT 2**(-103)
    !                   THUS CUTLO = 2**(-51) = 4.44089E-16
    !     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
    !                   THUS CUTHI = 2**(63.5) = 1.30438E19
    !     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
    !                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
    !     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
    !     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
    !     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
    DATA cutlo, cuthi / 8.232d-11, 1.304d19 /

    IF(n .GT. 0) GO TO 10
    dnrm2_ = ZERO
    GO TO 300

    10 assign 30 to next
    sum = ZERO
    nn = n * incx
    !                       BEGIN MAIN LOOP
    i = 1
    20    GO TO next, (30, 50, 70, 110)
    30 IF(ABS(dx(i)) .GT. cutlo) GO TO 85
    assign 50 to next
    xmax = ZERO

    !                        PHASE 1.  SUM IS ZERO

    50 IF(dx(i) .EQ. ZERO) GO TO 200
    IF(ABS(dx(i)) .GT. cutlo) GO TO 85

    !                        PREPARE FOR PHASE 2.

    assign 70 to next
    GO TO 105

    !                        PREPARE FOR PHASE 4.

    100 i = j
    assign 110 to next
    sum = (sum / dx(i)) / dx(i)
    105 xmax = ABS(dx(i))
    GO TO 115

    !                   PHASE 2.  SUM IS SMALL.
    !                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.

    70 IF(ABS(dx(i)) .GT. cutlo) GO TO 75

    !                   COMMON CODE FOR PHASES 2 AND 4.
    !                   IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.

    110 IF(ABS(dx(i)) .LE. xmax) GO TO 115
    sum = one + sum * (xmax / dx(i))**2
    xmax = ABS(dx(i))
    GO TO 200

    115 sum = sum + (dx(i) / xmax)**2
    GO TO 200

    !                  PREPARE FOR PHASE 3.

    75 sum = (sum * xmax) * xmax

    !     FOR REAL OR D.P. SET HITEST = CUTHI/N
    !     FOR COMPLEX      SET HITEST = CUTHI/(2*N)

    85 hitest = cuthi / float(n)

    !                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.

    DO 95 j = i, nn, incx
        IF(ABS(dx(j)) .GE. hitest) GO TO 100
    95    sum = sum + dx(j)**2
    dnrm2_ = SQRT(sum)
    GO TO 300

    200 CONTINUE
    i = i + incx
    IF (i .LE. nn) GO TO 20

    !              END OF MAIN LOOP.

    !              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.

    dnrm2_ = xmax * SQRT(sum)
    300 CONTINUE
    RETURN
END

SUBROUTINE  dsrot (n, dx, incx, dy, incy, c, s)

    !     APPLIES A PLANE ROTATION.
    !     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION dx(*), dy(*), dtemp, c, s
    INTEGER i, incx, incy, ix, iy, n

    IF(n.LE.0)RETURN
    IF(incx.EQ.1.AND.incy.EQ.1)GO TO 20

    !       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
    !         TO 1

    ix = 1
    iy = 1
    IF(incx.LT.0)ix = (-n + 1) * incx + 1
    IF(incy.LT.0)iy = (-n + 1) * incy + 1
    DO 10 i = 1, n
        dtemp = c * dx(ix) + s * dy(iy)
        dy(iy) = c * dy(iy) - s * dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
    10 CONTINUE
    RETURN

    !       CODE FOR BOTH INCREMENTS EQUAL TO 1

    20 DO 30 i = 1, n
        dtemp = c * dx(i) + s * dy(i)
        dy(i) = c * dy(i) - s * dx(i)
        dx(i) = dtemp
    30 CONTINUE
    RETURN
END

SUBROUTINE dsrotg(da, db, c, s)

    !     CONSTRUCT GIVENS PLANE ROTATION.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    !                    MODIFIED 9/27/86.

    DOUBLE PRECISION da, db, c, s, roe, scale, r, z, one, ZERO
    DATA one, ZERO /1.0d+00, 0.0d+00/

    roe = db
    IF(ABS(da) .GT. ABS(db)) roe = da
    scale = ABS(da) + ABS(db)
    IF(scale .NE. ZERO) GO TO 10
    c = one
    s = ZERO
    r = ZERO
    GO TO 20
    10 r = scale * SQRT((da / scale)**2 + (db / scale)**2)
    r = SIGN(one, roe) * r
    c = da / r
    s = db / r
    20 z = s
    IF(ABS(c) .GT. ZERO .AND. ABS(c) .LE. s) z = one / c
    da = r
    db = z
    RETURN
END

SUBROUTINE  dscal_sl(n, da, dx, incx)

    !     SCALES A VECTOR BY A CONSTANT.
    !     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
    !     JACK DONGARRA, LINPACK, 3/11/78.

    DOUBLE PRECISION da, dx(*)
    INTEGER i, incx, m, mp1, n, nincx

    IF(n.LE.0)RETURN
    IF(incx.EQ.1)GO TO 20


    !        CODE FOR INCREMENT NOT EQUAL TO 1

    nincx = n * incx
    DO 10 i = 1, nincx, incx
        dx(i) = da * dx(i)
    10 CONTINUE
    RETURN

    !        CODE FOR INCREMENT EQUAL TO 1

    !        CLEAN-UP LOOP

    20 m = MOD(n, 5)
    IF(m .EQ. 0) GO TO 40
    DO 30 i = 1, m
        dx(i) = da * dx(i)
    30 CONTINUE
    IF(n .LT. 5) RETURN
    40 mp1 = m + 1
    DO 50 i = mp1, n, 5
        dx(i) = da * dx(i)
        dx(i + 1) = da * dx(i + 1)
        dx(i + 2) = da * dx(i + 2)
        dx(i + 3) = da * dx(i + 3)
        dx(i + 4) = da * dx(i + 4)
    50 CONTINUE
    RETURN
END

subroutine bound(n, x, xl, xu)
    integer n, i
    double precision x(n), xl(n), xu(n)
    do i = 1, n
        !        Note that xl(i) and xu(i) may be NaN to indicate no bound
        if(xl(i).eq.xl(i).and.x(i) < xl(i))then
            x(i) = xl(i)
        else if(xu(i).eq.xu(i).and.x(i) > xu(i))then
            x(i) = xu(i)
        end if
    end do
end subroutine bound