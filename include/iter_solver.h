#ifndef _SISL_ITER_SOLVER_H_INCLUDED_
#define _SISL_ITER_SOLVER_H_INCLUDED_

#include <glib.h>

#include "vector.h"
#include "matrix.h"

typedef enum {
  SISL_ITER_JACOBI,
  SISL_ITER_GAUSS_SEIDEL,
  SISL_ITER_SOR,
  SISL_ITER_SSOR,
  SISL_ITER_CG,
  SISL_ITER_MINRES,
  SISL_ITER_SYMMLQ,
  SISL_ITER_GMRES,
  SISL_ITER_BICG,
  SISL_ITER_QMR,
  SISL_ITER_CGS,
  SISL_ITER_BICGSTAB
} sisl_solver_t ;

typedef struct {
  sisl_vector_t *v1, *v2, *v3, *v4, *v5, *v6, *v7, *v8, *v9, *v10 ;
} sisl_solver_workspace_t ;

typedef struct {
  gdouble r_1, r_2, r_inf, bnorm ;
  guint niter ;
} sisl_solver_performance_t ;

sisl_solver_workspace_t *sisl_solver_workspace_new(guint n,
						   data_complex_t rc,
						   sisl_dist_t d) ;
gint sisl_solve(sisl_solver_t solver,
		sisl_matrix_t *A, sisl_vector_t *x, sisl_vector_t *b,
		gdouble tol, guint niter,
		sisl_solver_workspace_t *w,
		sisl_solver_performance_t *perf) ;

#endif /*_SISL_ITER_SOLVER_H_INCLUDED_*/
