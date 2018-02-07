/* sisl.h
 * 
 * Copyright (C) 2009, 2017 Michael Carley
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

#ifndef SISL_H_INCLUDED
#define SISL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_complex.h>

#include <sislconfig.h>

  /**
   * @defgroup status Status codes
   * @{
   */

  /**
   * Return codes for SISL functions.
   * 
   */

  typedef enum {
    SISL_SUCCESS = 0,       /**< success */
    SISL_ERROR = 1,         /**< error*/
    SISL_NULL_ARGUMENT = 2  /**< an argument was NULL and should not have been*/
  } sisl_status_t ;

  /**
   * @}
   * 
   */

  /**
   * Matrix and vector norms
   * @ingroup vector
   */

  typedef enum {
    SISL_NORM_1,        /**< the 1 norm*/
    SISL_NORM_2,        /**< the 2 norm*/
    SISL_NORM_INFINITY  /**< the \f$\infty\f$ norm*/
  } sisl_norm_t ;

  /**
   * Matrix types, currently sparse matrices are not supported.
   * @ingroup matrix
   */

  typedef enum {
    SISL_MATRIX_DENSE,        /**< dense matrix*/
    SISL_MATRIX_SPARSE,	      /**< sparse matrix */
    SISL_MATRIX_USER_DEFINED  /**< user-defined matrix (e.g. for matrix-free
				 methods)*/
  } sisl_matrix_density_t ;

  /**
   * Self-explanatory.
   * @ingroup vector
   */

  typedef enum {
    SISL_REAL,    /**< real matrix or vector*/
    SISL_COMPLEX  /**< complex matrix or vector*/
  } sisl_complex_t ;

  /**
   * Distribution of a matrix
   * @ingroup matrix
   */

  typedef enum {
    SISL_SINGLE = 0,        /**< single processor holds the whole matrix*/
    SISL_DISTRIBUTED = 1    /**< matrix is distributed over more than one 
			       processor*/
  } sisl_distribution_t ;

/**
 * @ingroup solver
 *
 * Iterative solvers
 * 
 */
typedef enum {
  SISL_SOLVER_JACOBI,		
  SISL_SOLVER_GAUSS_SEIDEL,	
  SISL_SOLVER_SOR,
  SISL_SOLVER_SSOR,
  SISL_SOLVER_CG,			/**< Conjugate gradient */
  SISL_SOLVER_MINRES,		
  SISL_SOLVER_SYMMLQ,
  SISL_SOLVER_GMRES,
  SISL_SOLVER_BICG,		/**< Biconjugate gradient */
  SISL_SOLVER_QMR,
  SISL_SOLVER_CGS,		/**< Conjugate gradient squared */
  SISL_SOLVER_BICGSTAB		/**< Stabilized biconjugate gradient */
} sisl_solver_t ;

#define SISL_ROW_NUMBER_OFFSET 0
#define SISL_COL_NUMBER_OFFSET 1
#define SISL_ROW_START_OFFSET  2
#define SISL_ROW_END_OFFSET    3

  /**
   * The basic SISL matrix type. This should only be accessed through
   * the provided macros and functions.
   * @ingroup matrix
   * @hideinitializer
   */

  typedef struct _sisl_matrix_t sisl_matrix_t ;

  /**
   * The basic SISL vector type. This should only be accessed through
   * the provided macros and functions.
   * @ingroup vector
   * @hideinitializer
   */

  typedef struct _sisl_vector_t sisl_vector_t ;

  typedef struct _sisl_solver_workspace_t sisl_solver_workspace_t ;
  typedef struct _sisl_solver_performance_t sisl_solver_performance_t ;
  typedef struct _sisl_matrix_user_defined_t sisl_matrix_user_defined_t ;

  typedef gint (* sisl_matrix_vector_multiply_func)
    (sisl_matrix_t *, sisl_vector_t *, sisl_vector_t *) ;

  struct _sisl_matrix_t {
    sisl_matrix_density_t d ;
    sisl_complex_t c ;
    sisl_distribution_t dist ;
    gpointer m ;
  } ;

  struct _sisl_matrix_user_defined_t {
    gint size[4] ;
    sisl_matrix_vector_multiply_func multiply ;
    gpointer data ;
  } ;

  struct _sisl_vector_t {
    sisl_complex_t c ;
    GArray *x ;
  } ;

  struct _sisl_solver_workspace_t {
    GPtrArray *v ;
  } ;

  struct _sisl_solver_performance_t {
    gdouble r_1, r_2, r_inf, bnorm ;
    gint niter ;
  } ;

#define sisl_complex_string(c)			\
  (((c) == SISL_REAL) ? "real" : "complex")
#define sisl_vector_data(_v)				\
  ((gdouble *)(&(g_array_index(((_v)->x),gdouble,0))))
#define sisl_vector_data_length(v) (((v)->x->len))

  /**
   * @addtogroup vector
   * @{
   */


  /**
   * Check if a SISL type is real or complex, returning TRUE if real,
   * FALSE otherwise.
   * @hideinitializer
   */

#define sisl_is_real(v) ((v)->c == SISL_REAL)

  /**
   * Find the length of a ::sisl_vector_t.
   *
   * @hideinitializer
   */

#define sisl_vector_length(v)		\
  (((sisl_is_real(v)) ? ((v)->x->len) : ((v)->x->len/2)))

  /**
   * @}
   * 
   */


/*   /\** */
/*    * Set the length of a ::sisl_vector_t. */
/*    * @addtogroup vector */
/*    * @hideinitializer */
/*    *\/ */

/* #define sisl_vector_set_length(v,n)		\ */
/*   (((sisl_is_real(v)) ? (g_array_set_size((v)->x,(n))) :	\ */
/*     (g_array_set_size((v)->x,(2*n))))) */

/**
 * @addtogroup matrix
 * @{
 * 
 */


  /**
   * The number of rows in a ::sisl_matrix_t.
   *
   * @hideinitializer
   */

#define sisl_matrix_row_number(A)		\
  (((gint *)((A)->m))[SISL_ROW_NUMBER_OFFSET])

  /**
   * The number of columns in a ::sisl_matrix_t
   * 
   * @hideinitializer
   */

#define sisl_matrix_column_number(A)		\
  (((gint *)((A)->m))[SISL_COL_NUMBER_OFFSET])

  /**
   * The first locally held row of a ::sisl_matrix_t
   * 
   * @hideinitializer
   */

#define sisl_matrix_local_row_start(A) \
  (((gint *)((A)->m))[SISL_ROW_START_OFFSET])

  /**
   * The first locally held row of the next block of a ::sisl_matrix_t.
   * 
   * @hideinitializer
   */
#define sisl_matrix_local_row_end(A) \
  (((gint *)((A)->m))[SISL_ROW_END_OFFSET])

  /**
   * Check if a ::sisl_matrix_t stores row \a i locally.
   * 
   * @hideinitializer
   */
#define sisl_matrix_has_local_row(A,i)		\
  (((i)>=sisl_matrix_local_row_start((A)) &&	\
    ((i)<sisl_matrix_local_row_end((A)))))

#define sisl_matrix_user_defined_multiply(A)		\
  ((sisl_matrix_user_defined_t *)((A)->m))->multiply
#define sisl_matrix_user_defined_data(A)		\
  ((sisl_matrix_user_defined_t *)((A)->m))->data

  /**
   * The density of a ::sisl_matrix_t, a ::sisl_matrix_density_t
   * 
   * @hideinitializer
   */
#define sisl_matrix_density(m) (m->d)

  /**
   * The distribution of a ::sisl_matrix_t, a ::sisl_distribution_t.
   * 
   * @hideinitializer
   */

#define sisl_matrix_distribution(m) ((m)->dist)

  /**
   * @}
   * 
   */


  sisl_vector_t *sisl_vector_new(sisl_complex_t c) ;
  gint sisl_vector_free(sisl_vector_t *v) ;
  gdouble sisl_vector_get(sisl_vector_t *v, gint i) ;
  gsl_complex sisl_vector_get_complex(sisl_vector_t *v, gint i) ;
  gint sisl_vector_set(sisl_vector_t *v, gint i, gdouble x) ;
  gint sisl_vector_set_complex(sisl_vector_t *v, gint i, gsl_complex x) ;
  gint sisl_vector_element_add(sisl_vector_t *v, gint i, gdouble x) ;
  gint sisl_vector_element_add_complex(sisl_vector_t *v, gint i, 
				       gsl_complex x) ;
  gint sisl_vector_element_sub(sisl_vector_t *v, gint i, gdouble x) ;
  gint sisl_vector_element_sub_complex(sisl_vector_t *v, gint i, 
				       gsl_complex x) ;
  gint sisl_vector_set_length(sisl_vector_t *v, gint n) ;
  gint sisl_vector_multiprocessor_sum(sisl_vector_t *v) ;

  gdouble sisl_vector_norm(sisl_vector_t *v, sisl_norm_t n) ;
  gdouble sisl_vector_dot(sisl_vector_t *v, sisl_vector_t *w) ;
  gsl_complex sisl_vector_dot_complex(sisl_vector_t *v, sisl_vector_t *w) ;
  gint sisl_vector_clear(sisl_vector_t *v) ;
  gint sisl_vector_copy(sisl_vector_t *v, sisl_vector_t *w) ;
  gint sisl_vector_add(sisl_vector_t *v, sisl_vector_t *w) ;
  gint sisl_vector_sub(sisl_vector_t *v, sisl_vector_t *w) ;
  gint sisl_vector_scale(sisl_vector_t *v, gdouble x) ;
  gint sisl_vector_scale_complex(sisl_vector_t *v, gsl_complex x) ;
  gint sisl_vector_set_all(sisl_vector_t *v, gdouble x) ;
  gint sisl_vector_set_all_complex(sisl_vector_t *v, gsl_complex x) ;

  sisl_matrix_t *sisl_matrix_new(sisl_complex_t c, sisl_matrix_density_t d) ;
  gint sisl_matrix_free(sisl_matrix_t *m) ;
  gint sisl_matrix_set_block_size(sisl_matrix_t *m, gint rows, gint columns) ;
  gint sisl_matrix_vector_mul(sisl_matrix_t *m, sisl_vector_t *v,
			      sisl_vector_t *w) ;
  gint sisl_matrix_set(sisl_matrix_t *m, gint i, gint j, gdouble x) ;
  gint sisl_matrix_set_complex(sisl_matrix_t *m, gint i, gint j, 
			       gsl_complex x) ;
  gdouble sisl_matrix_get(sisl_matrix_t *m, gint i, gint j) ;
  gsl_complex sisl_matrix_get_complex(sisl_matrix_t *m, gint i, gint j) ;

  sisl_solver_workspace_t *sisl_solver_workspace_new(void) ;
  gint sisl_solver_workspace_size(sisl_solver_workspace_t *w, gint n) ;
  gint sisl_solver_workspace_free(sisl_solver_workspace_t *w) ;
  gint sisl_solve(sisl_solver_t solver,
		  sisl_matrix_t *A, sisl_vector_t *x, 
		  sisl_vector_t *b,
		  gdouble tol, gint niter,
		  sisl_solver_workspace_t *w,
		  sisl_solver_performance_t *perf) ;

gint sisl_logging_init(FILE *f, gchar *p, 
		       GLogLevelFlags log_level,
		       gpointer exit_func) ;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*SISL_H_INCLUDED*/
