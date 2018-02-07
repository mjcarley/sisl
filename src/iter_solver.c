/* iter_solver.c
 * 
 * Copyright (C) 2006 Michael Carley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/**
 * @defgroup solver Iterative solvers
 * @brief  Various iterative solvers for linear systems.
 * 
 * Functions implementing (some of) the iterative solvers in Barrett,
 * R. et al, `Templates for the solution of linear systems: Building
 * blocks for iterative methods', SIAM, 1994.
 * 
 * @{
 */

#include <stdio.h>
#include <stdlib.h>

#include <glib.h>
#include <wmpi.h>

#include "sisl.h"
#include "sisl_private.h"
#include "iter_solver.h"

/** 
 * Allocate a workspace for iterative solution of linear systems
 * 
 * @return pointer to newly allocated workspace
 */

sisl_solver_workspace_t *sisl_solver_workspace_new(void)

{
  sisl_solver_workspace_t *w ;

  w = (sisl_solver_workspace_t *)g_malloc(sizeof(sisl_solver_workspace_t)) ;
  
  w->v = g_ptr_array_new() ;

  return w ;
}

gint sisl_solver_workspace_free(sisl_solver_workspace_t *w)

{
  gint i ;

  for ( i = 0 ; i < SISL_WORKSPACE_SIZE(w) ; i ++ )
    sisl_vector_free((sisl_vector_t *)g_ptr_array_index(w->v,i)) ;
  g_ptr_array_free(w->v, TRUE) ;
  
  g_free(w) ;

  return 0 ;
}

gint sisl_solver_workspace_size(sisl_solver_workspace_t *w, gint n)

{
  gint i, j ;

  SISL_CHECK_ARGUMENT_NULL(w) ;

  if ( n < (j = SISL_WORKSPACE_SIZE(w)) ) return SISL_SUCCESS ;
  
  g_ptr_array_set_size(w->v, n) ;
  for ( i = j ; i < SISL_WORKSPACE_SIZE(w) ; i ++ ) 
    g_ptr_array_index(w->v,i) = sisl_vector_new(SISL_REAL) ;

  return SISL_SUCCESS ;
}

gint sisl_solve_bicgstab(sisl_matrix_t *A, sisl_vector_t *x, 
			 sisl_vector_t *b,
			 gdouble tol, gint niter,
			 sisl_solver_workspace_t *w,
			 sisl_solver_performance_t *perf)

     /*
       Stabilized biconjugate gradient solution of Ax=b,
       from Barrett, R. et al, `Templates for the solution of
       linear systems: Building blocks for iterative methods', 
       SIAM, 1994, figure 2.10, but without preconditioning.
     */

{
  sisl_vector_t *v, *r, *rt, *p, *ph, *s, *sh ;
  sisl_norm_t error_norm ;
  gdouble rho, rho_old, alpha, omega, beta, bnorm, rnorm_2, rnorm_inf ;
  gint i, j ;

  error_norm = SISL_NORM_2 ;

  /*unpack the workspace*/
  sisl_solver_workspace_size(w, 7) ;
  v  = SISL_WORKSPACE_VECTOR(w,0) ;
  r  = SISL_WORKSPACE_VECTOR(w,1) ;
  rt = SISL_WORKSPACE_VECTOR(w,2) ;
  p  = SISL_WORKSPACE_VECTOR(w,3) ;
  ph = SISL_WORKSPACE_VECTOR(w,4) ;
  s  = SISL_WORKSPACE_VECTOR(w,5) ;
  sh = SISL_WORKSPACE_VECTOR(w,6) ;

  v->c = r->c = rt->c = p->c = ph->c = s->c = sh->c = A->c ;

  sisl_vector_clear(v) ; 
  sisl_vector_clear(r) ; sisl_vector_clear(rt) ;
  sisl_vector_clear(p) ; sisl_vector_clear(ph) ; 
  sisl_vector_clear(s) ; sisl_vector_clear(sh) ;

  /*size the vectors*/
  i = sisl_matrix_row_number(A) ;
  j = sisl_matrix_column_number(A) ;
  if ( i != j ) 
    g_error("%s: %dx%d matrix A must be square", __FUNCTION__, i, j) ;
  if ( sisl_vector_length(x) != j) 
    g_error("%s: length %d vector x not compatible with %dx%d matrix A",
	    __FUNCTION__, sisl_vector_length(x), i, j) ;
  if ( sisl_vector_length(b) != i) 
    g_error("%s: length %d vector b not compatible with %dx%d matrix A",
	    __FUNCTION__, sisl_vector_length(b), i, j) ;
  
  sisl_vector_set_length(v, i) ; 
  sisl_vector_set_length(r, i) ; sisl_vector_set_length(rt, i) ;
  sisl_vector_set_length(p, i) ; sisl_vector_set_length(ph, i) ;
  sisl_vector_set_length(s, i) ; sisl_vector_set_length(sh, i) ;

  /* initial residual r = b - Ax */
  sisl_matrix_vector_mul(A, x, r) ;

  sisl_vector_sub(r, b) ; sisl_vector_scale(r, -1.0) ;
  sisl_vector_copy(rt, r) ;

  bnorm = sisl_vector_norm(b, SISL_NORM_2) ; rho = omega = alpha = 0.0 ;
  if ( wmpi_rank() == 0 )
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "bnorm: %1.2e", bnorm) ;
  perf->niter = 0 ;

  rnorm_2 = rnorm_inf = G_MAXDOUBLE ;
  for ( i = 0 ; i < niter ; i ++ ) {
    rho_old = rho ; rho = sisl_vector_dot(rt, r) ;

    if ( rho == 0.0 ) { 
      if ( wmpi_rank() == 0 )
	g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, 
	      "%s failure: rho=0.0", __FUNCTION__) ;
      return -1 ; 
    }

    if ( i == 0 ) sisl_vector_copy(p, r) ;
    else {
      beta = (rho/rho_old)*(alpha/omega) ;
      /* p = r + beta*(p-omega*v) */
      sisl_vector_scale(v, -omega) ; sisl_vector_add(p, v) ;
      sisl_vector_scale(p, beta) ; sisl_vector_add(p, r) ;
    }

    /* preconditioner point */
    sisl_vector_copy(ph, p) ;
    
    sisl_matrix_vector_mul(A, ph, v) ;

    alpha = sisl_vector_dot(rt, v) ;
    alpha = rho/alpha ;

    sisl_vector_copy(s, v) ; sisl_vector_scale(s, -alpha) ;
    sisl_vector_add(s, r) ;

    omega = sisl_vector_norm(s, error_norm) ;
    if ( omega < tol ) {
      sisl_vector_scale(ph, alpha) ; sisl_vector_add(x, ph) ;
      if ( wmpi_rank() == 0 )
	g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "norm: %e", omega) ;
      break ;
    }

    /* preconditioner point */
    sisl_vector_copy(sh, s) ;

    /* use r to store t for a few lines */
    sisl_matrix_vector_mul(A, sh, r) ;
  
    omega = sisl_vector_norm(r, error_norm) ;
    rnorm_2 = sisl_vector_dot(r, s) ; 
    omega = rnorm_2/omega/omega ;

    sisl_vector_scale(sh, omega) ; sisl_vector_add(x, sh) ;
    sisl_vector_scale(ph, alpha) ; sisl_vector_add(x, ph) ;

    /* r = s - omega t */
    sisl_vector_scale(r, -omega) ; sisl_vector_add(r, s) ;

    rnorm_2 = sisl_vector_norm(r, error_norm) ;
    rnorm_inf = sisl_vector_norm(r, SISL_NORM_INFINITY) ;
    if ( wmpi_rank() == 0 )
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, 
	    "iter: %4u; |r|_2/|b|_2=%1.2e; |r|_inf=%1.2e", 
	    i, rnorm_2/bnorm, rnorm_inf) ;
    if ( omega == 0.0 ) break ;
    if ( rnorm_2 < tol*bnorm ) break ;
  }

  perf->r_1 = sisl_vector_norm(r, SISL_NORM_1) ;
  perf->r_2 = sisl_vector_norm(r, SISL_NORM_2) ;
  perf->r_inf = sisl_vector_norm(r, SISL_NORM_INFINITY) ;
  perf->niter = i ;
  perf->bnorm = bnorm ;

  if ( (rnorm_2 < tol*bnorm) || (rnorm_inf < tol) ) return 0 ;

  return -1 ;
}

#if 0
gint sisl_solve_bicg(sisl_matrix_t *A, sisl_vector_t *x, 
		     sisl_vector_t *b,
		     gdouble tol, gint niter,
		     sisl_solver_workspace_t *w,
		     sisl_solver_performance_t *perf)

     /*
       Biconjugate gradient solution of Ax=b, from Barrett, R. et al,
       `Templates for the solution of linear systems: Building blocks
       for iterative methods', SIAM, 1994, figure 2.7, but without
       preconditioning.  
     */

{
  sisl_vector_t *r, *rt, *z, *zt, *p, *pt, *q, *qt ;
  sisl_norm_t error_norm ;
  gdouble rho, rho_old, alpha, beta, bnorm, rnorm_2, rnorm_inf ;
  gint i, j ;

  error_norm = SISL_NORM_2 ;

  /*unpack the workspace*/
  r = w->v1 ; rt = w->v2 ; z = w->v3 ; zt = w->v4 ; 
  p = w->v5 ; pt = w->v6 ; q = w->v7 ; qt = w->v8 ;

  r->c = rt->c = z->c = zt->c = p->c = pt->c = q->c = qt->c = A->c ;

  sisl_vector_clear(r) ; sisl_vector_clear(rt) ; 
  sisl_vector_clear(z) ; sisl_vector_clear(zt) ; 
  sisl_vector_clear(p) ; sisl_vector_clear(pt) ; 
  sisl_vector_clear(q) ; sisl_vector_clear(qt) ; 

  /*size the vectors*/
  i = sisl_matrix_row_number(A) ;
  j = sisl_matrix_column_number(A) ;
  g_assert(sisl_vector_length(x) == i) ; 
  g_assert(sisl_vector_length(b) == i) ;
  
  sisl_vector_set_length(r, i) ; sisl_vector_set_length(rt, i) ; 
  sisl_vector_set_length(z, i) ; sisl_vector_set_length(zt, i) ; 
  sisl_vector_set_length(p, i) ; sisl_vector_set_length(pt, i) ; 
  sisl_vector_set_length(q, i) ; sisl_vector_set_length(qt, i) ; 

  /* initial residual r = b - Ax */
  sisl_matrix_vector_mul(A, x, r) ;

  sisl_vector_sub(r, b) ; sisl_vector_scale(r, -1.0) ;
  sisl_vector_copy(rt, r) ;

  sisl_vector_norm(b, SISL_NORM_2, &bnorm) ; rho = alpha = 0.0 ;
  if ( wmpi_rank() == 0 )
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "bnorm: %1.16e", bnorm) ;
  perf->niter = 0 ;

  rnorm_2 = rnorm_inf = G_MAXDOUBLE ;
  for ( i = 0 ; i < niter ; i ++ ) {
    rho_old = rho ;
    /*preconditioner point*/
    sisl_vector_copy(z, r) ; 
    sisl_vector_copy(zt, rt) ; 

    rho = sisl_vector_dot(z, rt) ;

    if ( rho == 0.0 ) { 
      if ( wmpi_rank() == 0 )
	g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "%s failure: rho=0.0",
	      __FUNCTION__) ;
      return -1 ; 
    }

    if ( i == 0 ) {
      sisl_vector_copy(p, z) ; sisl_vector_copy(pt, zt) ; 
    } else {
      beta = rho/rho_old ;
      sisl_vector_scale(p, beta) ; sisl_vector_add(p, z) ;
      sisl_vector_scale(pt, beta) ; sisl_vector_add(pt, zt) ;
    }

    sisl_matrix_vector_mul(A, p, q) ;
    sisl_matrix_transpose_vector_mul(A, pt, qt) ;

    alpha = sisl_vector_dot(pt, q) ;
    alpha = rho/alpha ;

    sisl_vector_add_weighted(x, p, alpha) ;
    sisl_vector_add_weighted(r, q, -alpha) ;
    sisl_vector_add_weighted(rt, qt, -alpha) ;

    sisl_vector_norm(r, error_norm, &rnorm_2) ;
    sisl_vector_norm(r, SISL_NORM_INFINITY, &rnorm_inf) ;
    if ( wmpi_rank() == 0 )
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, 
	    "iter: %4u; |r|_2/|b|_2=%1.2e; |r|_inf=%1.2e", 
	    i, rnorm_2/bnorm, rnorm_inf) ;
    if ( rnorm_2 < tol*bnorm ) break ;
  }

  perf->r_1   = sisl_vector_norm(r, SISL_NORM_1) ;
  perf->r_2   = sisl_vector_norm(r, SISL_NORM_2) ;
  perf->r_inf = sisl_vector_norm(r, SISL_NORM_INFINITY) ;
  perf->niter = i ;
  perf->bnorm = bnorm ;

  if ( (rnorm_2 < tol*bnorm) || (rnorm_inf < tol) ) return 0 ;

  return -1 ;
}
#endif

gint sisl_solve_cgs(sisl_matrix_t *A, sisl_vector_t *x, 
		    sisl_vector_t *b,
		    gdouble tol, gint niter,
		    sisl_solver_workspace_t *w,
		    sisl_solver_performance_t *perf)

     /*
       Conjugate gradient squared solution of Ax=b, from Barrett,
       R. et al, `Templates for the solution of linear systems:
       Building blocks for iterative methods', SIAM, 1994, figure 2.9,
       but without preconditioning.  
     */

{
  sisl_vector_t *u, *uh, *v, *vh, *r, *rt, *p, *ph, *q, *qh ;
  sisl_norm_t error_norm ;
  gdouble rho, rho_old, alpha, beta, bnorm, rnorm_2, rnorm_inf ;
  gint i ;

  error_norm = SISL_NORM_2 ;

  /*unpack the workspace*/
  sisl_solver_workspace_size(w, 10) ;
  v  = SISL_WORKSPACE_VECTOR(w,0) ;
  r  = SISL_WORKSPACE_VECTOR(w,1) ;
  rt = SISL_WORKSPACE_VECTOR(w,2) ;
  p  = SISL_WORKSPACE_VECTOR(w,3) ;
  ph = SISL_WORKSPACE_VECTOR(w,4) ;
  u  = SISL_WORKSPACE_VECTOR(w,5) ;
  uh = SISL_WORKSPACE_VECTOR(w,6) ;
  q  = SISL_WORKSPACE_VECTOR(w,7) ;
  qh = SISL_WORKSPACE_VECTOR(w,8) ;
  vh = SISL_WORKSPACE_VECTOR(w,9) ;

  v->c = r->c = rt->c = p->c = ph->c = u->c = uh->c = q->c = qh->c = vh->c
    = A->c ;

  sisl_vector_clear(v) ; sisl_vector_clear(vh) ;
  sisl_vector_clear(r) ; sisl_vector_clear(rt) ;
  sisl_vector_clear(p) ; sisl_vector_clear(ph) ; 
  sisl_vector_clear(u) ; sisl_vector_clear(uh) ;
  sisl_vector_clear(q) ; sisl_vector_clear(qh) ;

  /*size the vectors*/
  i = sisl_matrix_row_number(A) ;
  g_assert(sisl_vector_length(x) == i) ; 
  g_assert(sisl_vector_length(b) == i) ;
  
  sisl_vector_set_length(v, i) ; sisl_vector_set_length(vh, i) ; 
  sisl_vector_set_length(r, i) ; sisl_vector_set_length(rt, i) ;
  sisl_vector_set_length(p, i) ; sisl_vector_set_length(ph, i) ;
  sisl_vector_set_length(u, i) ; sisl_vector_set_length(uh, i) ;

  /* initial residual r = b - Ax */
  sisl_matrix_vector_mul(A, x, r) ;

  sisl_vector_sub(r, b) ; sisl_vector_scale(r, -1.0) ;
  sisl_vector_copy(rt, r) ;

  bnorm = sisl_vector_norm(b, SISL_NORM_2) ; rho = alpha = 0.0 ;
  g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "bnorm: %1.16e", bnorm) ;
  perf->niter = 0 ;

  rnorm_2 = rnorm_inf = G_MAXDOUBLE ;
  for ( i = 0 ; i < niter ; i ++ ) {
    rho_old = rho ; rho = sisl_vector_dot(rt, r) ;

    if ( rho == 0.0 ) { 
      if ( wmpi_rank() == 0 )
	g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "failure: rho=0.0") ;
      return -1 ; 
    }

    if ( i == 0 ) {
      sisl_vector_copy(u, r) ;
      sisl_vector_copy(p, u) ;
    } else {
      beta = rho/rho_old ;
      /* p = res + beta*(p-omega) */
      sisl_vector_copy(u, q) ; sisl_vector_scale(u, beta) ;
      sisl_vector_add(u, r) ;

      sisl_vector_scale(p, beta) ; sisl_vector_add(p, q) ;
      sisl_vector_scale(p, beta) ; sisl_vector_add(p, u) ;
    }

    /* preconditioner point */
    sisl_vector_copy(ph, p) ;
    
    sisl_matrix_vector_mul(A, ph, vh) ;

    alpha = sisl_vector_dot(rt, vh) ;
    alpha = rho/alpha ;

    sisl_vector_copy(q, vh) ; sisl_vector_scale(q, -alpha) ;
    sisl_vector_add(q, u) ;

    /* preconditioner point */
    sisl_vector_copy(uh,u) ; sisl_vector_add(uh,q) ;

    sisl_vector_copy(qh, uh) ; sisl_vector_scale(qh, alpha) ;
    sisl_vector_add(x, qh) ;

    sisl_matrix_vector_mul(A, uh, qh) ;
    sisl_vector_scale(qh, -alpha) ;
    sisl_vector_add(r, qh) ;

    rnorm_2 = sisl_vector_norm(r, error_norm) ;
    rnorm_inf = sisl_vector_norm(r, SISL_NORM_INFINITY) ;
    if ( wmpi_rank() == 0 )    
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, 
	    "iter: %4u; |r|_2/|b|_2=%1.2e; |r|_inf=%1.2e", 
	    i, rnorm_2/bnorm, rnorm_inf) ;
    if ( rnorm_2 < tol*bnorm ) break ;
  }

  perf->r_1 = sisl_vector_norm(r, SISL_NORM_1) ;
  perf->r_2 = sisl_vector_norm(r, SISL_NORM_2) ;
  perf->r_inf = sisl_vector_norm(r, SISL_NORM_INFINITY) ;
  perf->niter = i ;
  perf->bnorm = bnorm ;

  if ( (rnorm_2 < tol*bnorm) || (rnorm_inf < tol) ) return 0 ;

  return -1 ;
}

#if 0
gint sisl_solve_cg(sisl_matrix_t *A, sisl_vector_t *x, sisl_vector_t *b,
		   gdouble tol, gint niter,
		   sisl_solver_workspace_t *w,
		   sisl_solver_performance_t *perf)

     /*
       Conjugate gradient solution of Ax=b, from Barrett, R. et al,
       `Templates for the solution of linear systems: Building blocks
       for iterative methods', SIAM, 1994, figure 2.5, but without
       preconditioning.  
     */

{
  sisl_vector_t *r, *p, *q, *z ;
  sisl_norm_t error_norm ;
  gdouble rho, rho_old, alpha, beta, bnorm, rnorm_2, rnorm_inf ;
  gint i, j ;

  error_norm = SISL_NORM_2 ;

  /*unpack the workspace*/
  r = w->v1 ; p = w->v2 ; q = w->v3 ; z = w->v4 ;

  sisl_vector_clear(r) ; sisl_vector_clear(p) ; 
  sisl_vector_clear(q) ; sisl_vector_clear(z) ; 

  /*size the vectors*/
  sisl_mat_size(A, &i, &j) ; g_assert(i == j) ;
  g_assert(sisl_vector_length(x) == i) ; 
  g_assert(sisl_vector_length(b) == i) ;
  
  sisl_vector_set_length(r, i) ; sisl_vector_set_length(p, i) ; 
  sisl_vector_set_length(q, i) ; sisl_vector_set_length(z, i) ;

  /* initial residual r = b - Ax */
  sisl_mat_vector_multiply(A, x, r) ;

  sisl_vector_sub(r, b) ; sisl_vector_scale(r, -1.0) ;

  sisl_vector_norm(b, SISL_NORM_2, &bnorm) ;
  g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "bnorm: %1.16e", bnorm) ;
  perf->niter = 0 ;

  rnorm_2 = rnorm_inf = G_MAXDOUBLE ;
  for ( i = 0 ; i < niter ; i ++ ) {
    /*preconditioner point*/
    sisl_vector_copy(z, r) ;

    rho_old = rho ; 
    sisl_vector_inner_product(r, z, &rho) ;

    if ( i == 0 ) {
      sisl_vector_copy(p, z) ;
    } else {
      beta = rho/rho_old ;
      /* p = res + beta*(p-omega) */
      sisl_vector_scale(p, beta) ;
      sisl_vector_add(p, z) ;
    }

    sisl_mat_vector_multiply(A, p, q) ;

    sisl_vector_inner_product(p, q, &alpha) ;
    alpha = rho/alpha ;

    sisl_vector_add_weighted(x, p, alpha) ;
    sisl_vector_add_weighted(r, q, -alpha) ;

    sisl_vector_norm(r, error_norm, &rnorm_2) ;
    sisl_vector_norm(r, SISL_NORM_INFINITY, &rnorm_inf) ;
    if ( wmpi_rank() == 0 )
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, 
	    "iter: %4u; |r|_2/|b|_2=%1.2e; |r|_inf=%1.2e", 
	    i, rnorm_2/bnorm, rnorm_inf) ;
    if ( rnorm_2 < tol*bnorm ) break ;
  }

  sisl_vector_norm(r, SISL_NORM_1, &(perf->r_1)) ;
  sisl_vector_norm(r, SISL_NORM_2, &(perf->r_2)) ;
  sisl_vector_norm(r, SISL_NORM_INFINITY, &(perf->r_inf)) ;
  perf->niter = i ;
  perf->bnorm = bnorm ;

  if ( (rnorm_2 < tol*bnorm) || (rnorm_inf < tol) ) return 0 ;

  return -1 ;
}
#endif

/** 
 * Iterative solution of a real or complex linear system using one of
 * the methods in Barrett, R. et al, `Templates for the solution of
 * linear systems: Building blocks for iterative methods', SIAM, 1994. 
 * Currently implemented methods are: ::SISL_SOLVER_CG,
 * ::SISL_SOLVER_BICG, ::SISL_SOLVER_CGS, ::SISL_SOLVER_BICGSTAB. 
 * 
 * @param solver iterative solver to use (::sisl_solver_t);
 * @param A left hand side matrix;
 * @param x solution;
 * @param b right hand side;
 * @param tol tolerance for solution;
 * @param niter maximum number of iterations;
 * @param w workspace allocated with ::sisl_solver_workspace_new;
 * @param perf solution performance data (convergence tolerance, etc.);
 * 
 * @return 0 on success
 */

gint sisl_solve(sisl_solver_t solver,
		sisl_matrix_t *A, sisl_vector_t *x, 
		sisl_vector_t *b,
		gdouble tol, gint niter,
		sisl_solver_workspace_t *w,
		sisl_solver_performance_t *perf)

{
  gint r ;

  r = 0 ;
  switch (solver) {
  default: g_assert_not_reached() ; break ;
/*   case SISL_SOLVER_JACOBI: g_assert_not_reached() ; break ; */
/*   case SISL_SOLVER_GAUSS_SEIDEL: g_assert_not_reached() ; break ; */
/*   case SISL_SOLVER_SOR: g_assert_not_reached() ; break ; */
/*   case SISL_SOLVER_SSOR: g_assert_not_reached() ; break ; */
/*   case SISL_SOLVER_CG:  */
/*     if ( sisl_is_real(A) ) */
/*       r = sisl_solve_cg(A, x, b, tol, niter, w, perf) ; */
/*     else */
/*       g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, */
/* 	    "%s: complex CG not supported", __FUNCTION__) ; */
/*     break ; */
/*   case SISL_SOLVER_MINRES: g_assert_not_reached() ; break ; */
/*   case SISL_SOLVER_SYMMLQ: g_assert_not_reached() ; break ; */
/*   case SISL_SOLVER_GMRES: g_assert_not_reached() ; break ; */
/*   case SISL_SOLVER_BICG:  */
/*     if ( sisl_is_real(A) ) */
/*       r = sisl_solve_bicg(A, x, b, tol, niter, w, perf) ; */
/*     else */
/*       g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, */
/* 	    "%s: complex BICG not supported", __FUNCTION__) ; */
/*     break ; */
/*   case SISL_SOLVER_QMR: g_assert_not_reached() ; break ; */
  case SISL_SOLVER_CGS:
    if ( sisl_is_real(A) )
      r = sisl_solve_cgs(A, x, b, tol, niter, w, perf) ;
    else
      r = sisl_solve_cgs_c(A, x, b, tol, niter, w, perf) ;
    break ;
  case SISL_SOLVER_BICGSTAB: 
    if ( sisl_is_real(A) )
      r = sisl_solve_bicgstab(A, x, b, tol, niter, w, perf) ;
    else
      r = sisl_solve_bicgstab_c(A, x, b, tol, niter, w, perf) ;
    break ;
  }

  return r ;
}

/**
 * @}
 * 
 */
