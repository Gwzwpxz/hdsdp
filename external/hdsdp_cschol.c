#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "external/hdsdp_cs.h"
#else
#include "hdsdp_utils.h"
#include "hdsdp_cs.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

/* CSparse routines for Sparse Chokesly. Referenced from Tim Davis Suite Sparse */

static int dcs_wclear (int mark, int lemax, int *w, int n) {
    
    int k ;
    if (mark < 2 || (mark + lemax < 0))
    {
        for (k = 0 ; k < n ; k++) if (w [k] != 0) w [k] = 1 ;
        mark = 2 ;
    }
    return (mark) ;     /* at this point, w [0..n-1] < mark holds */
}

/* keep off-diagonal entries; drop diagonal entries */
static int dcs_diag (int i, int j, double aij, void *other) {
    
    return (i != j) ;
}

/* free a numeric factorization */
dcsn *dcs_nfree (dcsn *N) {
    
    if (!N) return (NULL) ;     /* do nothing if N already NULL */
    dcs_spfree (N->L) ;
    dcs_spfree (N->U) ;
    dcs_free (N->pinv) ;
    dcs_free (N->B) ;
    return ((dcsn *) dcs_free (N)) ;  /* free the csn struct and return NULL */
}

/* free a symbolic factorization */
dcss *dcs_sfree (dcss *S) {
    
    if (!S) return (NULL) ;     /* do nothing if S already NULL */
    dcs_free (S->pinv) ;
    dcs_free (S->q) ;
    dcs_free (S->parent) ;
    dcs_free (S->cp) ;
    dcs_free (S->leftmost) ;
    return ((dcss *) dcs_free (S)) ;  /* free the css struct and return NULL */
}

dcsn *dcs_ndone (dcsn *N, dcs *C, void *w, void *x, int ok) {
    
    dcs_spfree (C) ;                     /* free temporary matrix */
    dcs_free (w) ;                       /* free workspace */
    dcs_free (x) ;
    return (ok ? N : dcs_nfree (N)) ;    /* return result if OK, else free it */
}

dcs *dcs_add (const dcs *A, const dcs *B, double alpha, double beta) {
    
    int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
    double *x, *Bx, *Cx ;
    dcs *C ;
    if (!IS_CSC (A) || !IS_CSC (B)) return (NULL) ;         /* check inputs */
    if (A->m != B->m || A->n != B->n) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
    w = dcs_calloc (m, sizeof (int)) ;                       /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? dcs_malloc (m, sizeof (double)) : NULL ;    /* get workspace */
    C = dcs_spalloc (m, n, anz + bnz, values, 0) ;           /* allocate result*/
    if (!C || !w || (values && !x)) return (dcs_done (C, w, x, 0)) ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)
    {
        Cp [j] = nz ;                   /* column j of C starts here */
        nz = dcs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        nz = dcs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    dcs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (dcs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/* p = amd(A+A') if symmetric is true, or amd(A'A) otherwise */
int *dcs_amd (int order, const dcs *A)  /* order 0:natural, 1:Chol, 2:LU, 3:QR */ {
    
    dcs *C, *A2, *AT ;
    int *Cp, *Ci, *last, *W, *len, *nv, *next, *P, *head, *elen, *degree, *w,
        *hhead, *ATp, *ATi, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
        k2, k3, jlast, ln, dense, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
        ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, n, m, t ;
    int h ;
    /* --- Construct matrix C ----------------------------------------------- */
    if (!IS_CSC (A) || order <= 0 || order > 3) return (NULL) ; /* check */
    AT = dcs_transpose (A, 0) ;              /* compute A' */
    if (!AT) return (NULL) ;
    m = A->m ; n = A->n ;
    dense = CS_MAX (16, 10 * sqrt ((double) n)) ;   /* find dense threshold */
    dense = CS_MIN (n-2, dense) ;
    if (order == 1 && n == m)
    {
        C = dcs_add (A, AT, 0, 0) ;          /* C = A+A' */
    }
    else if (order == 2)
    {
        ATp = AT->p ;                       /* drop dense columns from AT */
        ATi = AT->i ;
        for (p2 = 0, j = 0 ; j < m ; j++)
        {
            p = ATp [j] ;                   /* column j of AT starts here */
            ATp [j] = p2 ;                  /* new column j starts here */
            if (ATp [j+1] - p > dense) continue ;   /* skip dense col j */
            for ( ; p < ATp [j+1] ; p++) ATi [p2++] = ATi [p] ;
        }
        ATp [m] = p2 ;                      /* finalize AT */
        A2 = dcs_transpose (AT, 0) ;         /* A2 = AT' */
        C = A2 ? dcs_multiply (AT, A2) : NULL ;  /* C=A'*A with no dense rows */
        dcs_spfree (A2) ;
    }
    else
    {
        C = dcs_multiply (AT, A) ;           /* C=A'*A */
    }
    dcs_spfree (AT) ;
    if (!C) return (NULL) ;
    dcs_fkeep (C, &dcs_diag, NULL) ;          /* drop diagonal entries */
    Cp = C->p ;
    cnz = Cp [n] ;
    P = dcs_malloc (n+1, sizeof (int)) ;     /* allocate result */
    W = dcs_malloc (8*(n+1), sizeof (int)) ; /* get workspace */
    t = cnz + cnz/5 + 2*n ;                 /* add elbow room to C */
    if (!P || !W || !dcs_sprealloc (C, t)) return (dcs_idone (P, C, W, 0)) ;
    len  = W           ; nv     = W +   (n+1) ; next   = W + 2*(n+1) ;
    head = W + 3*(n+1) ; elen   = W + 4*(n+1) ; degree = W + 5*(n+1) ;
    w    = W + 6*(n+1) ; hhead  = W + 7*(n+1) ;
    last = P ;                              /* use P as workspace for last */
    /* --- Initialize quotient graph ---------------------------------------- */
    for (k = 0 ; k < n ; k++) len [k] = Cp [k+1] - Cp [k] ;
    len [n] = 0 ;
    nzmax = C->nzmax ;
    Ci = C->i ;
    for (i = 0 ; i <= n ; i++)
    {
        head [i] = -1 ;                     /* degree list i is empty */
        last [i] = -1 ;
        next [i] = -1 ;
        hhead [i] = -1 ;                    /* hash list i is empty */
        nv [i] = 1 ;                        /* node i is just one node */
        w [i] = 1 ;                         /* node i is alive */
        elen [i] = 0 ;                      /* Ek of node i is empty */
        degree [i] = len [i] ;              /* degree of node i */
    }
    mark = dcs_wclear (0, 0, w, n) ;         /* clear w */
    elen [n] = -2 ;                         /* n is a dead element */
    Cp [n] = -1 ;                           /* n is a root of assembly tree */
    w [n] = 0 ;                             /* n is a dead element */
    /* --- Initialize degree lists ------------------------------------------ */
    for (i = 0 ; i < n ; i++)
    {
        d = degree [i] ;
        if (d == 0)                         /* node i is empty */
        {
            elen [i] = -2 ;                 /* element i is dead */
            nel++ ;
            Cp [i] = -1 ;                   /* i is a root of assembly tree */
            w [i] = 0 ;
        }
        else if (d > dense)                 /* node i is dense */
        {
            nv [i] = 0 ;                    /* absorb i into element n */
            elen [i] = -1 ;                 /* node i is dead */
            nel++ ;
            Cp [i] = CS_FLIP (n) ;
            nv [n]++ ;
        }
        else
        {
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;           /* put node i in degree list d */
            head [d] = i ;
        }
    }
    while (nel < n)                         /* while (selecting pivots) do */
    {
        /* --- Select node of minimum approximate degree -------------------- */
        for (k = -1 ; mindeg < n && (k = head [mindeg]) == -1 ; mindeg++) ;
        if (next [k] != -1) last [next [k]] = -1 ;
        head [mindeg] = next [k] ;          /* remove k from degree list */
        elenk = elen [k] ;                  /* elenk = |Ek| */
        nvk = nv [k] ;                      /* # of nodes k represents */
        nel += nvk ;                        /* nv[k] nodes of A eliminated */
        /* --- Garbage collection ------------------------------------------- */
        if (elenk > 0 && cnz + mindeg >= nzmax)
        {
            for (j = 0 ; j < n ; j++)
            {
                if ((p = Cp [j]) >= 0)      /* j is a live node or element */
                {
                    Cp [j] = Ci [p] ;       /* save first entry of object */
                    Ci [p] = CS_FLIP (j) ;  /* first entry is now CS_FLIP(j) */
                }
            }
            for (q = 0, p = 0 ; p < cnz ; ) /* scan all of memory */
            {
                if ((j = CS_FLIP (Ci [p++])) >= 0)  /* found object j */
                {
                    Ci [q] = Cp [j] ;       /* restore first entry of object */
                    Cp [j] = q++ ;          /* new pointer to object j */
                    for (k3 = 0 ; k3 < len [j]-1 ; k3++) Ci [q++] = Ci [p++] ;
                }
            }
            cnz = q ;                       /* Ci [cnz...nzmax-1] now free */
        }
        /* --- Construct new element ---------------------------------------- */
        dk = 0 ;
        nv [k] = -nvk ;                     /* flag k as in Lk */
        p = Cp [k] ;
        pk1 = (elenk == 0) ? p : cnz ;      /* do in place if elen[k] == 0 */
        pk2 = pk1 ;
        for (k1 = 1 ; k1 <= elenk + 1 ; k1++)
        {
            if (k1 > elenk)
            {
                e = k ;                     /* search the nodes in k */
                pj = p ;                    /* list of nodes starts at Ci[pj]*/
                ln = len [k] - elenk ;      /* length of list of nodes in k */
            }
            else
            {
                e = Ci [p++] ;              /* search the nodes in e */
                pj = Cp [e] ;
                ln = len [e] ;              /* length of list of nodes in e */
            }
            for (k2 = 1 ; k2 <= ln ; k2++)
            {
                i = Ci [pj++] ;
                if ((nvi = nv [i]) <= 0) continue ; /* node i dead, or seen */
                dk += nvi ;                 /* degree[Lk] += size of node i */
                nv [i] = -nvi ;             /* negate nv[i] to denote i in Lk*/
                Ci [pk2++] = i ;            /* place i in Lk */
                if (next [i] != -1) last [next [i]] = last [i] ;
                if (last [i] != -1)         /* remove i from degree list */
                {
                    next [last [i]] = next [i] ;
                }
                else
                {
                    head [degree [i]] = next [i] ;
                }
            }
            if (e != k)
            {
                Cp [e] = CS_FLIP (k) ;      /* absorb e into k */
                w [e] = 0 ;                 /* e is now a dead element */
            }
        }
        if (elenk != 0) cnz = pk2 ;         /* Ci [cnz...nzmax] is free */
        degree [k] = dk ;                   /* external degree of k - |Lk\i| */
        Cp [k] = pk1 ;                      /* element k is in Ci[pk1..pk2-1] */
        len [k] = pk2 - pk1 ;
        elen [k] = -2 ;                     /* k is now an element */
        /* --- Find set differences ----------------------------------------- */
        mark = dcs_wclear (mark, lemax, w, n) ;  /* clear w if necessary */
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan 1: find |Le\Lk| */
        {
            i = Ci [pk] ;
            if ((eln = elen [i]) <= 0) continue ;/* skip if elen[i] empty */
            nvi = -nv [i] ;                      /* nv [i] was negated */
            wnvi = mark - nvi ;
            for (p = Cp [i] ; p <= Cp [i] + eln - 1 ; p++)  /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] >= mark)
                {
                    w [e] -= nvi ;          /* decrement |Le\Lk| */
                }
                else if (w [e] != 0)        /* ensure e is a live element */
                {
                    w [e] = degree [e] + wnvi ; /* 1st time e seen in scan 1 */
                }
            }
        }
        /* --- Degree update ------------------------------------------------ */
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan2: degree update */
        {
            i = Ci [pk] ;                   /* consider node i in Lk */
            p1 = Cp [i] ;
            p2 = p1 + elen [i] - 1 ;
            pn = p1 ;
            for (h = 0, d = 0, p = p1 ; p <= p2 ; p++)    /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] != 0)             /* e is an unabsorbed element */
                {
                    dext = w [e] - mark ;   /* dext = |Le\Lk| */
                    if (dext > 0)
                    {
                        d += dext ;         /* sum up the set differences */
                        Ci [pn++] = e ;     /* keep e in Ei */
                        h += e ;            /* compute the hash of node i */
                    }
                    else
                    {
                        Cp [e] = CS_FLIP (k) ;  /* aggressive absorb. e->k */
                        w [e] = 0 ;             /* e is a dead element */
                    }
                }
            }
            elen [i] = pn - p1 + 1 ;        /* elen[i] = |Ei| */
            p3 = pn ;
            p4 = p1 + len [i] ;
            for (p = p2 + 1 ; p < p4 ; p++) /* prune edges in Ai */
            {
                j = Ci [p] ;
                if ((nvj = nv [j]) <= 0) continue ; /* node j dead or in Lk */
                d += nvj ;                  /* degree(i) += |j| */
                Ci [pn++] = j ;             /* place j in node list of i */
                h += j ;                    /* compute hash for node i */
            }
            if (d == 0)                     /* check for mass elimination */
            {
                Cp [i] = CS_FLIP (k) ;      /* absorb i into k */
                nvi = -nv [i] ;
                dk -= nvi ;                 /* |Lk| -= |i| */
                nvk += nvi ;                /* |k| += nv[i] */
                nel += nvi ;
                nv [i] = 0 ;
                elen [i] = -1 ;             /* node i is dead */
            }
            else
            {
                degree [i] = CS_MIN (degree [i], d) ;   /* update degree(i) */
                Ci [pn] = Ci [p3] ;         /* move first node to end */
                Ci [p3] = Ci [p1] ;         /* move 1st el. to end of Ei */
                Ci [p1] = k ;               /* add k as 1st element in of Ei */
                len [i] = pn - p1 + 1 ;     /* new len of adj. list of node i */
                h = ((h<0) ? (-h):h) % n ;  /* finalize hash of i */
                next [i] = hhead [h] ;      /* place i in hash bucket */
                hhead [h] = i ;
                last [i] = h ;              /* save hash of i in last[i] */
            }
        }                                   /* scan2 is done */
        degree [k] = dk ;                   /* finalize |Lk| */
        lemax = CS_MAX (lemax, dk) ;
        mark = dcs_wclear (mark+lemax, lemax, w, n) ;    /* clear w */
        /* --- Supernode detection ------------------------------------------ */
        for (pk = pk1 ; pk < pk2 ; pk++)
        {
            i = Ci [pk] ;
            if (nv [i] >= 0) continue ;         /* skip if i is dead */
            h = last [i] ;                      /* scan hash bucket of node i */
            i = hhead [h] ;
            hhead [h] = -1 ;                    /* hash bucket will be empty */
            for ( ; i != -1 && next [i] != -1 ; i = next [i], mark++)
            {
                ln = len [i] ;
                eln = elen [i] ;
                for (p = Cp [i]+1 ; p <= Cp [i] + ln-1 ; p++) w [Ci [p]] = mark;
                jlast = i ;
                for (j = next [i] ; j != -1 ; ) /* compare i with all j */
                {
                    ok = (len [j] == ln) && (elen [j] == eln) ;
                    for (p = Cp [j] + 1 ; ok && p <= Cp [j] + ln - 1 ; p++)
                    {
                        if (w [Ci [p]] != mark) ok = 0 ;    /* compare i and j*/
                    }
                    if (ok)                     /* i and j are identical */
                    {
                        Cp [j] = CS_FLIP (i) ;  /* absorb j into i */
                        nv [i] += nv [j] ;
                        nv [j] = 0 ;
                        elen [j] = -1 ;         /* node j is dead */
                        j = next [j] ;          /* delete j from hash bucket */
                        next [jlast] = j ;
                    }
                    else
                    {
                        jlast = j ;             /* j and i are different */
                        j = next [j] ;
                    }
                }
            }
        }
        /* --- Finalize new element------------------------------------------ */
        for (p = pk1, pk = pk1 ; pk < pk2 ; pk++)   /* finalize Lk */
        {
            i = Ci [pk] ;
            if ((nvi = -nv [i]) <= 0) continue ;/* skip if i is dead */
            nv [i] = nvi ;                      /* restore nv[i] */
            d = degree [i] + dk - nvi ;         /* compute external degree(i) */
            d = CS_MIN (d, n - nel - nvi) ;
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;               /* put i back in degree list */
            last [i] = -1 ;
            head [d] = i ;
            mindeg = CS_MIN (mindeg, d) ;       /* find new minimum degree */
            degree [i] = d ;
            Ci [p++] = i ;                      /* place i in Lk */
        }
        nv [k] = nvk ;                      /* # nodes absorbed into k */
        if ((len [k] = p-pk1) == 0)         /* length of adj list of element k*/
        {
            Cp [k] = -1 ;                   /* k is a root of the tree */
            w [k] = 0 ;                     /* k is now a dead element */
        }
        if (elenk != 0) cnz = p ;           /* free unused space in Lk */
    }
    /* --- Postordering ----------------------------------------------------- */
    for (i = 0 ; i < n ; i++) Cp [i] = CS_FLIP (Cp [i]) ;/* fix assembly tree */
    for (j = 0 ; j <= n ; j++) head [j] = -1 ;
    for (j = n ; j >= 0 ; j--)              /* place unordered nodes in lists */
    {
        if (nv [j] > 0) continue ;          /* skip if j is an element */
        next [j] = head [Cp [j]] ;          /* place j in list of its parent */
        head [Cp [j]] = j ;
    }
    for (e = n ; e >= 0 ; e--)              /* place elements in lists */
    {
        if (nv [e] <= 0) continue ;         /* skip unless e is an element */
        if (Cp [e] != -1)
        {
            next [e] = head [Cp [e]] ;      /* place e in list of its parent */
            head [Cp [e]] = e ;
        }
    }
    for (k = 0, i = 0 ; i <= n ; i++)       /* postorder the assembly tree */
    {
        if (Cp [i] == -1) k = dcs_tdfs (i, k, head, next, P, w) ;
    }
    return (dcs_idone (P, C, W, 1)) ;
}

#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
static void init_ata (dcs *AT, const int *post, int *w, int **head, int **next) {
    
    int i, k, p, m = AT->n, n = AT->m, *ATp = AT->p, *ATi = AT->i ;
    *head = w+4*n, *next = w+5*n+1 ;
    for (k = 0 ; k < n ; k++) w [post [k]] = k ;    /* invert post */
    for (i = 0 ; i < m ; i++)
    {
        for (k = n, p = ATp[i] ; p < ATp[i+1] ; p++) k = CS_MIN (k, w [ATi[p]]);
        (*next) [i] = (*head) [k] ;     /* place row i in linked list k */
        (*head) [k] = i ;
    }
}
int *dcs_counts (const dcs *A, const int *parent, const int *post, int ata) {
    
    int i, j, k, n, m, J, s, p, q, jleaf, *ATp, *ATi, *maxfirst, *prevleaf,
        *ancestor, *head = NULL, *next = NULL, *colcount, *w, *first, *delta ;
    dcs *AT ;
    if (!IS_CSC (A) || !parent || !post) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ;
    s = 4*n + (ata ? (n+m+1) : 0) ;
    delta = colcount = dcs_malloc (n, sizeof (int)) ;    /* allocate result */
    w = dcs_malloc (s, sizeof (int)) ;                   /* get workspace */
    AT = dcs_transpose (A, 0) ;                          /* AT = A' */
    if (!AT || !colcount || !w) return (dcs_idone (colcount, AT, w, 0)) ;
    ancestor = w ; maxfirst = w+n ; prevleaf = w+2*n ; first = w+3*n ;
    for (k = 0 ; k < s ; k++) w [k] = -1 ;      /* clear workspace w [0..s-1] */
    for (k = 0 ; k < n ; k++)                   /* find first [j] */
    {
        j = post [k] ;
        delta [j] = (first [j] == -1) ? 1 : 0 ;  /* delta[j]=1 if j is a leaf */
        for ( ; j != -1 && first [j] == -1 ; j = parent [j]) first [j] = k ;
    }
    ATp = AT->p ; ATi = AT->i ;
    if (ata) init_ata (AT, post, w, &head, &next) ;
    for (i = 0 ; i < n ; i++) ancestor [i] = i ; /* each node in its own set */
    for (k = 0 ; k < n ; k++)
    {
        j = post [k] ;          /* j is the kth node in postordered etree */
        if (parent [j] != -1) delta [parent [j]]-- ;    /* j is not a root */
        for (J = HEAD (k,j) ; J != -1 ; J = NEXT (J))   /* J=j for LL'=A case */
        {
            for (p = ATp [J] ; p < ATp [J+1] ; p++)
            {
                i = ATi [p] ;
                q = dcs_leaf (i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
                if (jleaf >= 1) delta [j]++ ;   /* A(i,j) is in skeleton */
                if (jleaf == 2) delta [q]-- ;   /* account for overlap in q */
            }
        }
        if (parent [j] != -1) ancestor [j] = parent [j] ;
    }
    for (j = 0 ; j < n ; j++)           /* sum up delta's of each child */
    {
        if (parent [j] != -1) colcount [parent [j]] += colcount [j] ;
    }
    return (dcs_idone (colcount, AT, w, 1)) ;    /* success: free workspace */
}


int dcs_ereach (const dcs *A, int k, const int *parent, int *s, int *w) {
    
    int i, p, n, len, top, *Ap, *Ai ;
    if (!IS_CSC (A) || !parent || !s || !w) return (-1) ;   /* check inputs */
    top = n = A->n ; Ap = A->p ; Ai = A->i ;
    CS_MARK (w, k) ;                /* mark node k as visited */
    for (p = Ap [k] ; p < Ap [k+1] ; p++)
    {
        i = Ai [p] ;                /* A(i,k) is nonzero */
        if (i > k) continue ;       /* only use upper triangular part of A */
        for (len = 0 ; !CS_MARKED (w,i) ; i = parent [i]) /* traverse up etree*/
        {
            s [len++] = i ;         /* L(k,i) is nonzero */
            CS_MARK (w, i) ;        /* mark i as visited */
        }
        while (len > 0) s [--top] = s [--len] ; /* push path onto stack */
    }
    for (p = top ; p < n ; p++) CS_MARK (w, s [p]) ;    /* unmark all nodes */
    CS_MARK (w, k) ;                /* unmark node k */
    return (top) ;                  /* s [top..n-1] contains pattern of L(k,:)*/
}

int *dcs_etree (const dcs *A, int ata) {
    
    int i, k, p, m, n, inext, *Ap, *Ai, *w, *parent, *ancestor, *prev ;
    if (!IS_CSC (A)) return (NULL) ;        /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ;
    parent = dcs_malloc (n, sizeof (int)) ;              /* allocate result */
    w = dcs_malloc (n + (ata ? m : 0), sizeof (int)) ;   /* get workspace */
    if (!w || !parent) return (dcs_idone (parent, NULL, w, 0)) ;
    ancestor = w ; prev = w + n ;
    if (ata) for (i = 0 ; i < m ; i++) prev [i] = -1 ;
    for (k = 0 ; k < n ; k++)
    {
        parent [k] = -1 ;                   /* node k has no parent yet */
        ancestor [k] = -1 ;                 /* nor does k have an ancestor */
        for (p = Ap [k] ; p < Ap [k+1] ; p++)
        {
            i = ata ? (prev [Ai [p]]) : (Ai [p]) ;
            for ( ; i != -1 && i < k ; i = inext)   /* traverse from i to k */
            {
                inext = ancestor [i] ;              /* inext = ancestor of i */
                ancestor [i] = k ;                  /* path compression */
                if (inext == -1) parent [i] = k ;   /* no anc., parent is k */
            }
            if (ata) prev [Ai [p]] = k ;
        }
    }
    return (dcs_idone (parent, NULL, w, 1)) ;
}

int dcs_fkeep (dcs *A, int (*fkeep) (int, int, double, void *), void *other) {
    
    int j, p, nz = 0, n, *Ap, *Ai ;
    double *Ax ;
    if (!IS_CSC (A) || !fkeep) return (-1) ;    /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        Ap [j] = nz ;                       /* record new location of col j */
        for ( ; p < Ap [j+1] ; p++)
        {
            if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other))
            {
                if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
                Ai [nz++] = Ai [p] ;
            }
        }
    }
    Ap [n] = nz ;                           /* finalize A */
    dcs_sprealloc (A, 0) ;                   /* remove extra space from A */
    return (nz) ;
}

int dcs_ipvec (const int *p, const double *b, double *x, int n) {
    
    int k ;
    if (!x || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) x [p ? p [k] : k] = b [k] ;
    return (1) ;
}

int dcs_leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf) {
    
    int q, s, sparent, jprev ;
    if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1) ;
    *jleaf = 0 ;
    if (i <= j || first [j] <= maxfirst [i]) return (-1) ;  /* j not a leaf */
    maxfirst [i] = first [j] ;      /* update max first[j] seen so far */
    jprev = prevleaf [i] ;          /* jprev = previous leaf of ith subtree */
    prevleaf [i] = j ;
    *jleaf = (jprev == -1) ? 1: 2 ; /* j is first or subsequent leaf */
    if (*jleaf == 1) return (i) ;   /* if 1st leaf, q = root of ith subtree */
    for (q = jprev ; q != ancestor [q] ; q = ancestor [q]) ;
    for (s = jprev ; s != q ; s = sparent)
    {
        sparent = ancestor [s] ;    /* path compression */
        ancestor [s] = q ;
    }
    return (q) ;                    /* q = least common ancester (jprev,j) */
}

int dcs_lsolve (const dcs *L, double *x) {
    
    int p, j, n, *Lp, *Li ;
    double *Lx ;
    if (!IS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = 0 ; j < n ; j++)
    {
        x [j] /= Lx [Lp [j]] ;
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [Li [p]] -= Lx [p] * x [j] ;
        }
    }
    return (1) ;
}

int dcs_ltsolve (const dcs *L, double *x) {
    
    int p, j, n, *Lp, *Li ;
    double *Lx ;
    if (!IS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [j] -= Lx [p] * x [Li [p]] ;
        }
        x [j] /= Lx [Lp [j]] ;
    }
    return (1) ;
}

dcs *dcs_multiply (const dcs *A, const dcs *B) {
    
    int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    double *x, *Bx, *Cx ;
    dcs *C ;
    if (!IS_CSC (A) || !IS_CSC (B)) return (NULL) ;      /* check inputs */
    if (A->n != B->m) return (NULL) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
    w = dcs_calloc (m, sizeof (int)) ;                    /* get workspace */
    values = (A->x != NULL) && (Bx != NULL) ;
    x = values ? dcs_malloc (m, sizeof (double)) : NULL ; /* get workspace */
    C = dcs_spalloc (m, n, anz + bnz, values, 0) ;        /* allocate result */
    if (!C || !w || (values && !x)) return (dcs_done (C, w, x, 0)) ;
    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !dcs_sprealloc (C, 2*(C->nzmax)+m))
        {
            return (dcs_done (C, w, x, 0)) ;             /* out of memory */
        }
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = dcs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    dcs_sprealloc (C, 0) ;               /* remove extra space from C */
    return (dcs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

int *dcs_pinv (int const *p, int n) {
    
    int k, *pinv ;
    if (!p) return (NULL) ;                     /* p = NULL denotes identity */
    pinv = dcs_malloc (n, sizeof (int)) ;        /* allocate result */
    if (!pinv) return (NULL) ;                  /* out of memory */
    for (k = 0 ; k < n ; k++) pinv [p [k]] = k ;/* invert the permutation */
    return (pinv) ;                             /* return result */
}

int dcs_pvec (const int *p, const double *b, double *x, int n) {
    
    int k ;
    if (!x || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) x [k] = b [p ? p [k] : k] ;
    return (1) ;
}

int *dcs_post (const int *parent, int n) {
    int j, k = 0, *post, *w, *head, *next, *stack ;
    if (!parent) return (NULL) ;                        /* check inputs */
    post = dcs_malloc (n, sizeof (int)) ;                /* allocate result */
    w = dcs_malloc (3*n, sizeof (int)) ;                 /* get workspace */
    if (!w || !post) return (dcs_idone (post, NULL, w, 0)) ;
    head = w ; next = w + n ; stack = w + 2*n ;
    for (j = 0 ; j < n ; j++) head [j] = -1 ;           /* empty linked lists */
    for (j = n-1 ; j >= 0 ; j--)            /* traverse nodes in reverse order*/
    {
        if (parent [j] == -1) continue ;    /* j is a root */
        next [j] = head [parent [j]] ;      /* add j to list of its parent */
        head [parent [j]] = j ;
    }
    for (j = 0 ; j < n ; j++)
    {
        if (parent [j] != -1) continue ;    /* skip j if it is not a root */
        k = dcs_tdfs (j, k, head, next, post, stack) ;
    }
    return (dcs_idone (post, NULL, w, 1)) ;  /* success; free w, return post */
}

int dcs_scatter (const dcs *A, int j, double beta, int *w, double *x, int mark, dcs *C, int nz) {
    
    int i, p, *Ap, *Ai, *Ci ;
    double *Ax ;
    if (!IS_CSC (A) || !w || !IS_CSC (C)) return (-1) ;     /* check inputs */
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}

dcs *dcs_symperm (const dcs *A, const int *pinv, int values) {
    
    int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w ;
    double *Cx, *Ax ;
    dcs *C ;
    if (!IS_CSC (A)) return (NULL) ;                    /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = dcs_spalloc (n, n, Ap [n], values && (Ax != NULL), 0) ; /* alloc result*/
    w = dcs_calloc (n, sizeof (int)) ;                   /* get workspace */
    if (!C || !w) return (dcs_done (C, w, NULL, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)           /* count entries in each column of C */
    {
        j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            if (i > j) continue ;       /* skip lower triangular part of A */
            i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
            w [CS_MAX (i2, j2)]++ ;     /* column count of C */
        }
    }
    dcs_cumsum (Cp, w, n) ;              /* compute column pointers of C */
    for (j = 0 ; j < n ; j++)
    {
        j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            if (i > j) continue ;       /* skip lower triangular part of A*/
            i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
            Ci [q = w [CS_MAX (i2, j2)]++] = CS_MIN (i2, j2) ;
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return (dcs_done (C, w, NULL, 1)) ;  /* success; free workspace, return C */
}

int dcs_tdfs (int j, int k, int *head, const int *next, int *post, int *stack) {
    
    int i, p, top = 0 ;
    if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
    stack [0] = j ;                 /* place j on the stack */
    while (top >= 0)                /* while (stack is not empty) */
    {
        p = stack [top] ;           /* p = top of stack */
        i = head [p] ;              /* i = youngest child of p */
        if (i == -1)
        {
            top-- ;                 /* p has no unordered children left */
            post [k++] = p ;        /* node p is the kth postordered node */
        }
        else
        {
            head [p] = next [i] ;   /* remove i from children of p */
            stack [++top] = i ;     /* start dfs on child node i */
        }
    }
    return (k) ;
}

dcs *dcs_transpose (const dcs *A, int values) {
    
    int p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
    double *Cx, *Ax ;
    dcs *C ;
    if (!IS_CSC (A)) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = dcs_spalloc (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
    w = dcs_calloc (m, sizeof (int)) ;                      /* get workspace */
    if (!C || !w) return (dcs_done (C, w, NULL, 0)) ;       /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    dcs_cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return (dcs_done (C, w, NULL, 1)) ;  /* success; free w and return C */
}

/* Sparse Cholesky routines */
dcss *dcs_schol (int order, const dcs *A) {
    
    int n, *c, *post, *P ;
    dcs *C ;
    dcss *S ;
    if (!IS_CSC (A)) return (NULL) ;        /* check inputs */
    n = A->n ;
    S = dcs_calloc (1, sizeof (dcss)) ;       /* allocate result S */
    if (!S) return (NULL) ;                 /* out of memory */
    P = dcs_amd (order, A) ;                 /* P = amd(A+A'), or natural */
    S->pinv = dcs_pinv (P, n) ;              /* find inverse permutation */
    dcs_free (P) ;
    if (order && !S->pinv) return (dcs_sfree (S)) ;
    C = dcs_symperm (A, S->pinv, 0) ;        /* C = spones(triu(A(P,P))) */
    S->parent = dcs_etree (C, 0) ;           /* find etree of C */
    post = dcs_post (S->parent, n) ;         /* postorder the etree */
    c = dcs_counts (C, S->parent, post, 0) ; /* find column counts of chol(C) */
    dcs_free (post) ;
    dcs_spfree (C) ;
    S->cp = dcs_malloc (n+1, sizeof (int)) ; /* allocate result S->cp */
    S->unz = S->lnz = dcs_cumsum (S->cp, c, n) ; /* find column pointers for L */
    dcs_free (c) ;
    return ((S->lnz >= 0) ? S : dcs_sfree (S)) ;
}

dcsn *dcs_chol (const dcs *A, const dcss *S) {
    
    double d, lki, *Lx, *x, *Cx ;
    int top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci ;
    dcs *L, *C, *E ;
    dcsn *N ;
    if (!IS_CSC (A) || !S || !S->cp || !S->parent) return (NULL) ;
    n = A->n ;
    N = dcs_calloc (1, sizeof (dcsn)) ;       /* allocate result */
    c = dcs_malloc (2*n, sizeof (int)) ;     /* get int workspace */
    x = dcs_malloc (n, sizeof (double)) ;    /* get double workspace */
    cp = S->cp ; pinv = S->pinv ; parent = S->parent ;
    C = pinv ? dcs_symperm (A, pinv, 1) : ((dcs *) A) ;
    E = pinv ? C : NULL ;           /* E is alias for A, or a copy E=A(p,p) */
    if (!N || !c || !x || !C) return (dcs_ndone (N, E, c, x, 0)) ;
    s = c + n ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    N->L = L = dcs_spalloc (n, n, cp [n], 1, 0) ;    /* allocate result */
    if (!L) return (dcs_ndone (N, E, c, x, 0)) ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (k = 0 ; k < n ; k++) Lp [k] = c [k] = cp [k] ;
    for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
    {
        /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        top = dcs_ereach (C, k, parent, s, c) ;      /* find pattern of L(k,:) */
        x [k] = 0 ;                                 /* x (0:k) is now zero */
        for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
        {
            if (Ci [p] <= k) x [Ci [p]] = Cx [p] ;
        }
        d = x [k] ;                     /* d = C(k,k) */
        x [k] = 0 ;                     /* clear x for k+1st iteration */
        /* --- Triangular solve --------------------------------------------- */
        for ( ; top < n ; top++)    /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        {
            i = s [top] ;               /* s [top..n-1] is pattern of L(k,:) */
            lki = x [i] / Lx [Lp [i]] ; /* L(k,i) = x (i) / L(i,i) */
            x [i] = 0 ;                 /* clear x for k+1st iteration */
            for (p = Lp [i] + 1 ; p < c [i] ; p++)
            {
                x [Li [p]] -= Lx [p] * lki ;
            }
            d -= lki * lki ;            /* d = d - L(k,i)*L(k,i) */
            p = c [i]++ ;
            Li [p] = k ;                /* store L(k,i) in column i */
            Lx [p] = lki ;
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        if (d <= 0) return (dcs_ndone (N, E, c, x, 0)) ; /* not pos def */
        p = c [k]++ ;
        Li [p] = k ;                /* store L(k,k) = sqrt (d) in column k */
        Lx [p] = sqrt (d) ;
    }
    Lp [n] = cp [n] ;               /* finalize L */
    return (dcs_ndone (N, E, c, x, 1)) ; /* success: free E,s,x; return N */
}

int dcs_cholsol (int order, const dcs *A, double *b) {
    
    double *x ;
    dcss *S ;
    dcsn *N ;
    int n, ok ;
    if (!IS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;
    S = dcs_schol (order, A) ;               /* ordering and symbolic analysis */
    N = dcs_chol (A, S) ;                    /* numeric Cholesky factorization */
    x = dcs_malloc (n, sizeof (double)) ;    /* get workspace */
    ok = (S && N && x) ;
    if (ok)
    {
        dcs_ipvec (S->pinv, b, x, n) ;   /* x = P*b */
        dcs_lsolve (N->L, x) ;           /* x = L\x */
        dcs_ltsolve (N->L, x) ;          /* x = L'\x */
        dcs_pvec (S->pinv, x, b, n) ;    /* b = P'*x */
    }
    dcs_free (x) ;
    dcs_sfree (S) ;
    dcs_nfree (N) ;
    return (ok) ;
}
