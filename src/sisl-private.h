/* sisl_private.h
 * 
 * Copyright (C) 2009, 2012 Michael Carley
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

#ifndef SISL_PRIVATE_H_INCLUDED
#define SISL_PRIVATE_H_INCLUDED

typedef struct _sisl_matrix_dense_t sisl_matrix_dense_t ;

struct _sisl_matrix_dense_t {
  gint size[4] ;
  GArray *x ;
} ;

typedef struct _sisl_matrix_sparse_t sisl_matrix_sparse_t ;

struct _sisl_matrix_sparse_t {
  gint size[4] ;
  GArray *ij, *x ;
} ;

typedef struct _sisl_matrix_diagonal_t sisl_matrix_diagonal_t ;

struct _sisl_matrix_diagonal_t {
  gint size[4] ;
  GArray *x ;
} ;

/*declarations of BLAS routines*/
extern gint c_0, c_1, c_2 ;
extern gdouble zero[] ;
extern gdouble one[] ;
extern gdouble minus_one[] ;
extern gdouble d_0, d_1, d_m1 ;

extern void dgemv_(gchar *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;
extern void zgemv_(gchar *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;

extern void dgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;
extern void zgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;

extern void zgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;
extern void dgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;

extern void dgemm_(gchar *transa, gchar *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;
extern void zgemm_(gchar *transa, gchar *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;

extern gdouble dasum_ (gint *n, gdouble *x, gint *incx) ;
extern gdouble dzasum_(gint *n, gdouble *x, gint *incx) ;
extern gdouble dnrm2_ (gint *n, gdouble *x, gint *incx) ;
extern gdouble dznrm2_(gint *n, gdouble *x, gint *incx) ;
extern gint    idamax_(gint *n, gdouble *x, gint *incx) ;
extern gint    izamax_(gint *n, gdouble *x, gint *incx) ;
extern gdouble ddot_  (gint *n, gdouble *x, gint *incx, 
		       gdouble *y, gint *incy) ;
extern gsl_complex zdotu_ (gint *n, gdouble *x, gint *incx, 
			   gdouble *y, gint *incy) ;
extern void    dcopy_(gint *n, 
		      gdouble *x, const gint *incx,
		      gdouble *y, const gint *incy) ;
extern void    daxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;
extern void    zaxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;

#define SISL_MATRIX_REAL_STRING(rc)		\
  (((rc) == SISL_REAL ? "real" : "complex"))

#define SISL_MATRIX_DISTRIBUTION_STRING(d)	\
  (((d) == SISL_SINGLE ? "single" : "distributed"))

#define SISL_MATRIX_DENSE_DATA(A) ((sisl_matrix_dense_t *)(A->m))
#define SISL_MATRIX_SPARSE_DATA(A) ((sisl_matrix_sparse_t *)(A->m))
#define SISL_MATRIX_DIAGONAL_DATA(A) ((sisl_matrix_diagonal_t *)(A->m))

#define SISL_WORKSPACE_SIZE(w) (w->v->len)
#define SISL_WORKSPACE_VECTOR(w,i) (g_ptr_array_index(w->v,i))

#define SISL_CHECK_ARGUMENT_NULL(a)		\
  if ( a == NULL )				\
    g_error("%s: argument %s may not be NULL", __FUNCTION__, #a) 

#define SISL_CHECK_ARGUMENT_REAL(a)			\
  if ( !sisl_is_real(a) )				\
    g_error("%s: %s must be real", __FUNCTION__, #a) ;
#define SISL_CHECK_ARGUMENT_COMPLEX(a)				\
  if ( sisl_is_real(a) )					\
    g_error("%s: %s must be complex", __FUNCTION__, #a) ;

gint sisl_sparse_mul(sisl_matrix_sparse_t *m, gdouble *v, gdouble *w) ;
gint sisl_sparse_mul_complex(sisl_matrix_sparse_t *m, gdouble *v, gdouble *w) ;
gint sisl_dense_sparse_mul_w(sisl_matrix_dense_t *A,
			     sisl_matrix_sparse_t *B,
			     gdouble a, gdouble c,
			     sisl_matrix_dense_t *C) ;
gint sisl_dense_sparse_mul_w_complex(sisl_matrix_dense_t *A,
				     sisl_matrix_sparse_t *B,
				     gsl_complex a, gsl_complex c,
				     sisl_matrix_dense_t *C) ;

gint sisl_diagonal_mul(sisl_matrix_diagonal_t *m, gdouble *v,
		       gdouble *w) ;
gint sisl_diagonal_mul_complex(sisl_matrix_diagonal_t *m, gdouble *v,
			 gdouble *w) ;
gint sisl_diagonal_invert(sisl_matrix_diagonal_t *m) ;
gint sisl_diagonal_invert_complex(sisl_matrix_diagonal_t *m) ;

#endif /*SISL_PRIVATE_H_INCLUDED*/
