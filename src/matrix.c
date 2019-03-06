/* matrix.c
 * 
 * Copyright (C) 2009, 2019 Michael Carley
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

#include <string.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <sisl.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <wmpi.h>

#include "sisl-private.h"

gdouble zero[] = {0.0, 0.0} ;
gdouble one[] = {1.0, 0.0} ;
gdouble minus_one[] = {-1.0, 0.0} ;
gint c_0 = 0 ;
gint c_1 = 1 ;
gint c_2 = 2 ;
gdouble d_0  =  0.0 ;
gdouble d_1  =  1.0 ;
gdouble d_m1 = -1.0 ;

static gint array_search(gint *i, gint len, gint s)

{
  gint j0, j, j1 ;

  if ( len == 0 ) return 0 ;
  if ( i[len-1] < s ) return len ;
  if ( i[0] > s ) return 0 ;

  j = j0 = 0 ; j1 = len-1 ;

  while ( j0 < j1 ) {
    j = j0 + (j1 - j0)/2 ;
    if ( i[j] == s ) return j ;
    if ( i[j] < s ) j0 = j+1 ;
    else j1 = j-1 ;
  }
  /*index not found but we need to find where to insert it*/
  j = MIN(j0,MIN(j,j1)) ;
  while ( i[j] < s ) j ++ ;
  return j ;
}

/**
 * @defgroup matrix SISL matrices
 *
 * The basic matrix type in SISL is the ::sisl_matrix_t. This can be
 * real or complex and have numerical values or be of a user-defined
 * type, e.g. for matrix-free computations. It should only be accessed
 * using the functions provided, with the possible exception of
 * user-defined matrices (but if you know enough to be able to define
 * your own matrix methods, you should know enough to handle them like
 * an adult). Built-in matrix types will reallocate memory
 * automatically when the user resizes them or adds elements to a
 * sparse matrix. The matrix operations are actually wrappers for the
 * corresponding BLAS types, so they should be fairly efficient. The
 * underlying data type is a glib GArray, which allows for efficient
 * allocation and resizing.
 *
 * @{
 */

/** 
 * Allocate a new ::sisl_matrix_t, by default with distribution
 * ::SISL_SINGLE.
 * 
 * @param c SISL_REAL or SISL_COMPLEX;
 * @param d SISL_MATRIX_DENSE or SISL_MATRIX_USER_DEFINED.
 * 
 * @return a newly allocated ::sisl_matrix_t. You will need to set its
 * size and, possibly allocate its data block
 * (::sisl_matrix_set_block_size) before you can use it.
 */

sisl_matrix_t *sisl_matrix_new(sisl_complex_t c, sisl_matrix_density_t d)

{
  sisl_matrix_t *m ;

  m = g_malloc(sizeof(sisl_matrix_t)) ;
  m->c = c ; m->d = d ;

  switch ( d ) {
  default: 
    g_error("%s: unrecognized matrix density type %d", __FUNCTION__, d) ;
    break ;
  case SISL_MATRIX_USER_DEFINED:
    m->m = g_malloc(sizeof(sisl_matrix_user_defined_t)) ;
    sisl_matrix_row_number(m) = 0 ;
    sisl_matrix_column_number(m) = 0 ;
    sisl_matrix_local_row_start(m) = 0 ;
    sisl_matrix_local_row_end(m) = 0 ;
    sisl_matrix_user_defined_multiply(m) = NULL ;
    sisl_matrix_user_defined_data(m) = NULL ;
    sisl_matrix_distribution(m) = SISL_SINGLE ;
    break ;
  case SISL_MATRIX_DENSE:
    m->m = g_malloc(sizeof(sisl_matrix_dense_t)) ;
    ((sisl_matrix_dense_t *)(m->m))->x = 
      g_array_new(FALSE, FALSE, sizeof(gdouble)) ;
    sisl_matrix_row_number(m) = 0 ;
    sisl_matrix_column_number(m) = 0 ;
    sisl_matrix_local_row_start(m) = 0 ;
    sisl_matrix_local_row_end(m) = 0 ;
    sisl_matrix_distribution(m) = SISL_SINGLE ;
    break ;
  case SISL_MATRIX_SPARSE:
    m->m = g_malloc(sizeof(sisl_matrix_sparse_t)) ;
    ((sisl_matrix_sparse_t *)(m->m))->x = 
      g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
    ((sisl_matrix_sparse_t *)(m->m))->ij = 
      g_array_new(TRUE, TRUE, sizeof(gint)) ;
    sisl_matrix_row_number(m) = 0 ;
    sisl_matrix_column_number(m) = 0 ;
    sisl_matrix_local_row_start(m) = 0 ;
    sisl_matrix_local_row_end(m) = 0 ;
    sisl_matrix_distribution(m) = SISL_SINGLE ;
    break ;
  case SISL_MATRIX_DIAGONAL:
    m->m = g_malloc(sizeof(sisl_matrix_diagonal_t)) ;
    ((sisl_matrix_diagonal_t *)(m->m))->x = 
      g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
    sisl_matrix_row_number(m) = 0 ;
    sisl_matrix_column_number(m) = 0 ;
    sisl_matrix_local_row_start(m) = 0 ;
    sisl_matrix_local_row_end(m) = 0 ;
    sisl_matrix_distribution(m) = SISL_SINGLE ;
    break ;
  }

  return m ;
}

/** 
 * Free a ::sisl_matrix_t and the underlying memory.
 * 
 * @param m a ::sisl_matrix_t.
 * 
 * @return SISL_SUCCESS on success.
 */

gint sisl_matrix_free(sisl_matrix_t *m)

{
  SISL_CHECK_ARGUMENT_NULL(m) ;

  switch ( sisl_matrix_density(m) ) {
  default: 
    g_error("%s: unrecognized matrix density type %d", __FUNCTION__, 
	    sisl_matrix_density(m)) ;
    break ;
  case SISL_MATRIX_USER_DEFINED:
    if ( sisl_matrix_user_defined_data(m) != NULL ) 
      g_warning("%s: non-NULL user data in user-defined matrix", 
		__FUNCTION__) ;
    g_free(m) ;
    break ;
  case SISL_MATRIX_DENSE:
    g_array_free(SISL_MATRIX_DENSE_DATA(m)->x, TRUE) ;
    g_free(m->m) ;
    g_free(m) ;
    break ;
  }

  return SISL_SUCCESS ;
}

/** 
 * Allocate the data block for a ::sisl_matrix_t. On a single
 * processor, this will be the same as the size of the matrix itself;
 * on a distributed system, this will be the size of the locally held
 * part of the matrix.
 * 
 * @param m a ::sisl_matrix_t;
 * @param rows the number of rows to allocate;
 * @param columns the number of columns to allocate.
 * 
 * @return SISL_SUCCESS on success.
 */

gint sisl_matrix_set_block_size(sisl_matrix_t *m, gint rows, gint columns)

{
  SISL_CHECK_ARGUMENT_NULL(m) ;

  if ( (rows < 0) || (columns < 0) ) 
    g_error("%s: invalid matrix size (%dx%d) requested",
	    __FUNCTION__, rows, columns) ;

  switch ( sisl_matrix_density(m) ) {
  default:
    g_error("%s: unrecognized matrix density type %d in matrix m", 
	    __FUNCTION__, sisl_matrix_density(m)) ;
    break ;
  case SISL_MATRIX_DENSE:
    if ( sisl_is_real(m) )
      g_array_set_size(SISL_MATRIX_DENSE_DATA(m)->x, rows*columns) ;
    else
      g_array_set_size(SISL_MATRIX_DENSE_DATA(m)->x, 2*rows*columns) ;
    break ;
  case SISL_MATRIX_DIAGONAL:
    if ( sisl_is_real(m) )
      g_array_set_size(SISL_MATRIX_DIAGONAL_DATA(m)->x, rows) ;
    else
      g_array_set_size(SISL_MATRIX_DIAGONAL_DATA(m)->x, 2*rows) ;
    break ;
  }

  return SISL_SUCCESS ;
}

/** 
 * Set an element of a ::sisl_matrix_t, \f$m_{ij}=x\f$. If the indices
 * \f$(i,j)\f$ do not lie in the data block of the matrix, the
 * function returns silently without setting any element.
 * 
 * @param m a ::sisl_matrix_t;
 * @param i row index;
 * @param j column index;
 * @param x value for \f$m_{ij}\f$.
 * 
 * @return SISL_SUCCESS on success, an error code if the indices are
 * not valid for the matrix as a whole.
 */

gint sisl_matrix_set(sisl_matrix_t *m, gint i, gint j, gdouble x)

{
  gint k, idx, *ij ;

  SISL_CHECK_ARGUMENT_NULL(m) ;
  SISL_CHECK_ARGUMENT_REAL(m) ;

  if ( (i < 0) || (i > sisl_matrix_row_number(m)) ||
       (j < 0) || (j > sisl_matrix_column_number(m)) ) 
    g_error("%s: element (%d,%d) is out of range for %dx%d matrix",
	    __FUNCTION__, i, j, 
	    sisl_matrix_row_number(m),
	    sisl_matrix_column_number(m)) ;

  if ( i < sisl_matrix_local_row_start(m) ||
       i >= sisl_matrix_local_row_end(m) ) return SISL_SUCCESS ;

  switch ( sisl_matrix_density(m) ) {
  default:
    g_error("%s: unrecognized matrix density type %d in matrix m", 
	    __FUNCTION__, sisl_matrix_density(m)) ;
    break ;
  case SISL_MATRIX_DENSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return SISL_SUCCESS ;
    g_array_index(SISL_MATRIX_DENSE_DATA(m)->x,
		  gdouble,
		  ((i-sisl_matrix_local_row_start(m))*
		   sisl_matrix_column_number(m)+j)) = x ;
    break ;
  case SISL_MATRIX_SPARSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return SISL_SUCCESS ;
    idx = i*sisl_matrix_column_number(m)+j ;
    ij = (gint *)(SISL_MATRIX_SPARSE_DATA(m)->ij->data) ;
    k = array_search(ij, SISL_MATRIX_SPARSE_DATA(m)->ij->len, idx) ;
    if ( k == SISL_MATRIX_SPARSE_DATA(m)->ij->len ) {
      g_array_append_val(SISL_MATRIX_SPARSE_DATA(m)->ij,idx) ;
      g_array_append_val(SISL_MATRIX_SPARSE_DATA(m)->x,x) ;
      return SISL_SUCCESS ;
    }
    if ( g_array_index(SISL_MATRIX_SPARSE_DATA(m)->ij,gint,k) != idx ) {
      g_array_insert_val(SISL_MATRIX_SPARSE_DATA(m)->ij,k,idx) ;
      g_array_insert_val(SISL_MATRIX_SPARSE_DATA(m)->x,k,x) ;
      return SISL_SUCCESS ;
    }
    g_array_index(SISL_MATRIX_SPARSE_DATA(m)->x,gdouble,k) = x ;
    return SISL_SUCCESS ;
    break ;
  }

  return SISL_SUCCESS ;
}

/** 
 * Set an element of a complex ::sisl_matrix_t, \f$m_{ij}=z\f$. If the
 * indices \f$(i,j)\f$ do not lie in the data block of the matrix, the
 * function returns silently without setting any element.
 * 
 * @param m a ::sisl_matrix_t;
 * @param i row index;
 * @param j column index;
 * @param z a ::gsl_complex number.
 * 
 * @return SISL_SUCCESS on success, an error code if the indices are
 * not valid for the matrix as a whole.
 */

gint sisl_matrix_set_complex(sisl_matrix_t *m, gint i, gint j, 
			     gsl_complex z)
  
{
  gint idx, k, *ij ;

  SISL_CHECK_ARGUMENT_NULL(m) ;
  SISL_CHECK_ARGUMENT_COMPLEX(m) ;

  if ( (i < 0) || (i > sisl_matrix_row_number(m)) ||
       (j < 0) || (j > sisl_matrix_column_number(m)) ) 
    g_error("%s: element (%d,%d) is out of range for %dx%d matrix",
	    __FUNCTION__, i, j, 
	    sisl_matrix_row_number(m),
	    sisl_matrix_column_number(m)) ;

  if ( i < sisl_matrix_local_row_start(m) ||
       i >= sisl_matrix_local_row_end(m) ) return SISL_SUCCESS ;

  switch ( sisl_matrix_density(m) ) {
  default:
    g_error("%s: unrecognized matrix density type %d in matrix m", 
	    __FUNCTION__, sisl_matrix_density(m)) ;
    break ;
  case SISL_MATRIX_DENSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return SISL_SUCCESS ;
    g_array_index(SISL_MATRIX_DENSE_DATA(m)->x,
		  gdouble,
		  (2*((i-sisl_matrix_local_row_start(m))*
		      sisl_matrix_column_number(m)+j)+0)) = GSL_REAL(z) ;
    g_array_index(SISL_MATRIX_DENSE_DATA(m)->x,
		  gdouble,
		  (2*((i-sisl_matrix_local_row_start(m))*
		      sisl_matrix_column_number(m)+j)+1)) = GSL_IMAG(z) ;
    break ;
  case SISL_MATRIX_SPARSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return SISL_SUCCESS ;
    idx = i*sisl_matrix_column_number(m)+j ;
    ij = (gint *)(SISL_MATRIX_SPARSE_DATA(m)->ij->data) ;
    k = array_search(ij, SISL_MATRIX_SPARSE_DATA(m)->ij->len, idx) ;
    if ( k == SISL_MATRIX_SPARSE_DATA(m)->ij->len ) {
      g_array_append_val(SISL_MATRIX_SPARSE_DATA(m)->ij,idx) ;
      g_array_append_vals(SISL_MATRIX_SPARSE_DATA(m)->x,&(GSL_REAL(z)),2) ;
      return SISL_SUCCESS ;
    }
    if ( g_array_index(SISL_MATRIX_SPARSE_DATA(m)->ij,gint,k) != idx ) {
      g_array_insert_val(SISL_MATRIX_SPARSE_DATA(m)->ij,k,idx) ;
      g_array_insert_vals(SISL_MATRIX_SPARSE_DATA(m)->x,2*k,&(GSL_REAL(z)),2) ;
      return SISL_SUCCESS ;
    }
    g_array_index(SISL_MATRIX_SPARSE_DATA(m)->x,gdouble,2*k+0) = 
      GSL_REAL(z) ;
    g_array_index(SISL_MATRIX_SPARSE_DATA(m)->x,gdouble,2*k+1) = 
      GSL_IMAG(z) ;
    return SISL_SUCCESS ;
    break ;
  case SISL_MATRIX_DIAGONAL:
    if ( i != j )
      g_error("%s: row (%d) and column (%d) indices must be equal "
	      "for diagonal matrix", __FUNCTION__, i, j) ;
    g_array_index(SISL_MATRIX_DIAGONAL_DATA(m)->x,gdouble,2*i+0) = 
      GSL_REAL(z) ;
    g_array_index(SISL_MATRIX_DIAGONAL_DATA(m)->x,gdouble,2*i+1) = 
      GSL_IMAG(z) ;    
    break ;    
  }

  return SISL_SUCCESS ;
}

/** 
 * Multiply a vector by a matrix, \f$w=mv\f$. The result vector is
 * automatically sized to hold the result.
 * 
 * @param m a ::sisl_matrix_t;
 * @param v a ::sisl_vector_t;
 * @param w a ::sisl_vector_t, to contain the result of the multiplication.
 * 
 * @return SISL_SUCCESS on success
 */

gint sisl_matrix_vector_mul(sisl_matrix_t *m, sisl_vector_t *v,
			    sisl_vector_t *w)

{
  gint nr, nc ;

  SISL_CHECK_ARGUMENT_NULL(m) ;
  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_NULL(w) ;

  g_debug("%s: %dx%d matrix with distribution %s",
	  __FUNCTION__, 
	  sisl_matrix_row_number(m),
	  sisl_matrix_column_number(m),
	  SISL_MATRIX_DISTRIBUTION_STRING(sisl_matrix_distribution(m))) ;

  if ( sisl_matrix_column_number(m) != sisl_vector_length(v) )
    g_error("%s: cannot multiply length %d vector v by %dx%d matrix m",
	    __FUNCTION__, sisl_vector_length(v), 
	    sisl_matrix_row_number(m), sisl_matrix_column_number(m)) ;

  if ( (sisl_matrix_local_row_start(m) < 0) ||
       (sisl_matrix_local_row_start(m) >= sisl_matrix_local_row_end(m)) ||
       (sisl_matrix_local_row_end(m) > sisl_matrix_row_number(m)) )
    g_error("%s: local row indices %d--%d are invalid for %dx%d matrix", 
	    __FUNCTION__,
	    sisl_matrix_local_row_start(m), sisl_matrix_local_row_end(m),
	    sisl_matrix_row_number(m), sisl_matrix_column_number(m)) ;

  if ( m->c != v->c ) 
    g_error("%s: matrix m and vector v must be both real or both complex",
	    __FUNCTION__) ;

  w->c = v->c ; sisl_vector_clear(w) ;
  sisl_vector_set_length(w, sisl_matrix_row_number(m)) ;
  nr = sisl_matrix_local_row_end(m) - sisl_matrix_local_row_start(m) ;
  nc = sisl_matrix_column_number(m) ;

  if ( sisl_is_real(m) ) sisl_vector_set_all(w, 0.0) ;
  else sisl_vector_set_all_complex(w, GSL_COMPLEX_ZERO) ;

  switch ( sisl_matrix_density(m) ) {
  default:
    g_error("%s: unrecognized matrix density type %d in matrix m", 
	    __FUNCTION__, sisl_matrix_density(m)) ;
    break ;
  case SISL_MATRIX_USER_DEFINED:
    ((sisl_matrix_user_defined_t *)(m->m))->multiply(m, v, w) ;
    break ;
  case SISL_MATRIX_DENSE:
    if ( sisl_is_real(m) ) 
      dgemv_("T", &nc, &nr, &d_1, 
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(m)->x->data),
	     &nc, sisl_vector_data(v), &c_1, &d_0, 
	     &((sisl_vector_data(w))[sisl_matrix_local_row_start(m)]),
	     &c_1) ;
    else 
      zgemv_("T", &nc, &nr, one,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(m)->x->data),
	     &nc, sisl_vector_data(v), &c_1, zero,
	     &((sisl_vector_data(w))[2*sisl_matrix_local_row_start(m)]),
	     &c_1) ;
    break ;
  case SISL_MATRIX_SPARSE:
    if ( sisl_matrix_distribution(m) != SISL_SINGLE ) 
      g_error("%s: parallel sparse matrix multiplication is "
	      "not implemented yet", __FUNCTION__) ;
    if ( sisl_is_real(m) ) 
      sisl_sparse_mul(SISL_MATRIX_SPARSE_DATA(m), 
		      sisl_vector_data(v), sisl_vector_data(w)) ;
    else
      sisl_sparse_mul_complex(SISL_MATRIX_SPARSE_DATA(m), 
			      sisl_vector_data(v), sisl_vector_data(w)) ;
    break ;
  case SISL_MATRIX_DIAGONAL:
    if ( sisl_matrix_distribution(m) != SISL_SINGLE ) 
      g_error("%s: parallel diagonal matrix multiplication is "
	      "not implemented yet", __FUNCTION__) ;
    if ( sisl_is_real(m) ) 
      sisl_diagonal_mul(SISL_MATRIX_DIAGONAL_DATA(m), 
			sisl_vector_data(v), sisl_vector_data(w)) ;
    else
      sisl_diagonal_mul_complex(SISL_MATRIX_DIAGONAL_DATA(m), 
				sisl_vector_data(v), sisl_vector_data(w)) ;
    break ;
  }

  if ( sisl_matrix_distribution(m) == SISL_SINGLE ) return SISL_SUCCESS ;

  sisl_vector_multiprocessor_sum(w) ;

  return SISL_SUCCESS ;
}

/** 
 * Get an element of a a ::sisl_matrix_t.
 * 
 * @param m a ::sisl_matrix_t;
 * @param i row index;
 * @param j column index.
 * 
 * @return element (\a i, \a j) of \a m.
 */

gdouble sisl_matrix_get(sisl_matrix_t *m, gint i, gint j)

{
  gint k, idx ;

  SISL_CHECK_ARGUMENT_NULL(m) ;

  if ( (i < 0) || (i > sisl_matrix_row_number(m)) ||
       (j < 0) || (j > sisl_matrix_column_number(m)) ) 
    g_error("%s: element (%d,%d) is out of range for %dx%d matrix",
	    __FUNCTION__, i, j, 
	    sisl_matrix_row_number(m),
	    sisl_matrix_column_number(m)) ;

  if ( i < sisl_matrix_local_row_start(m) ||
       i >= sisl_matrix_local_row_end(m) ) 
    g_error("%s: element (%d,%d) is not in data block of matrix, rows %d--%d",
	    __FUNCTION__, i, j, 
	    sisl_matrix_local_row_start(m), sisl_matrix_local_row_end(m)) ;

  switch ( sisl_matrix_density(m) ) {
  default:
    g_error("%s: unrecognized matrix density type %d in matrix m", 
	    __FUNCTION__, sisl_matrix_density(m)) ;
    break ;
  case SISL_MATRIX_DENSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return 0.0 ;
    return g_array_index(SISL_MATRIX_DENSE_DATA(m)->x,
			 gdouble,
			 ((i-sisl_matrix_local_row_start(m))*
			  sisl_matrix_column_number(m)+j)) ;
    break ;
  case SISL_MATRIX_SPARSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return 0.0 ;
    idx = i*sisl_matrix_column_number(m)+j ;
    k = array_search((gint *)(SISL_MATRIX_SPARSE_DATA(m)->ij->data),
		     SISL_MATRIX_SPARSE_DATA(m)->ij->len, idx) ;
    if ( g_array_index(SISL_MATRIX_SPARSE_DATA(m)->ij,gint,k) != idx )
      return 0.0 ;
    return g_array_index(SISL_MATRIX_SPARSE_DATA(m)->x, gdouble, k) ;
    break ;
  }

  return 0.0 ;
}

/** 
 * Get an element of a complex ::sisl_matrix_t.
 * 
 * @param m a ::sisl_matrix_t;
 * @param i row index;
 * @param j column index.
 * 
 * @return the complex element (\a i, \a j) of \a m.
 */

gsl_complex sisl_matrix_get_complex(sisl_matrix_t *m, gint i, gint j)

{
  gsl_complex x ;
  gint k, idx ;

  SISL_CHECK_ARGUMENT_NULL(m) ;
  GSL_SET_COMPLEX(&x,0,0) ;

  if ( (i < 0) || (i > sisl_matrix_row_number(m)) ||
       (j < 0) || (j > sisl_matrix_column_number(m)) ) 
    g_error("%s: element (%d,%d) is out of range for %dx%d matrix",
	    __FUNCTION__, i, j, 
	    sisl_matrix_row_number(m),
	    sisl_matrix_column_number(m)) ;
  if ( i < sisl_matrix_local_row_start(m) ||
       i >= sisl_matrix_local_row_end(m) ) 
    g_error("%s: element (%d,%d) is not in data block of matrix, rows %d--%d",
	    __FUNCTION__, i, j, 
	    sisl_matrix_local_row_start(m), sisl_matrix_local_row_end(m)) ;

  switch ( sisl_matrix_density(m) ) {
  default:
    g_error("%s: unrecognized matrix density type %d in matrix m", 
	    __FUNCTION__, sisl_matrix_density(m)) ;
    break ;
  case SISL_MATRIX_DENSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return x ;
    GSL_REAL(x) = g_array_index(SISL_MATRIX_DENSE_DATA(m)->x,
				gdouble,
				(2*((i-sisl_matrix_local_row_start(m))*
				    sisl_matrix_column_number(m)+j)+0)) ;
    GSL_IMAG(x) = g_array_index(SISL_MATRIX_DENSE_DATA(m)->x,
				gdouble,
				(2*((i-sisl_matrix_local_row_start(m))*
				    sisl_matrix_column_number(m)+j)+1)) ;
    return x ;
    break ;
  case SISL_MATRIX_SPARSE:
    if ( (i >= sisl_matrix_local_row_end(m)) ||
	 (i < sisl_matrix_local_row_start(m)) ) return x ;
    idx = i*sisl_matrix_column_number(m)+j ;
    k = array_search((gint *)(SISL_MATRIX_SPARSE_DATA(m)->ij->data),
		     SISL_MATRIX_SPARSE_DATA(m)->ij->len, idx) ;
    if ( g_array_index(SISL_MATRIX_SPARSE_DATA(m)->ij,gint,k) != idx )
      return x ;
    GSL_REAL(x) = g_array_index(SISL_MATRIX_SPARSE_DATA(m)->x,gdouble,
				2*k+0) ;
    GSL_IMAG(x) = g_array_index(SISL_MATRIX_SPARSE_DATA(m)->x,gdouble,
				2*k+1) ;
    return x ;
    break ;
  case SISL_MATRIX_DIAGONAL:
    if ( i != j )
      g_error("%s: row (%d) and column (%d) indices must be equal "
	      "for diagonal matrix", __FUNCTION__, i, j) ;
    GSL_REAL(x) = g_array_index(SISL_MATRIX_DIAGONAL_DATA(m)->x,gdouble,2*i+0) ;
    GSL_IMAG(x) = g_array_index(SISL_MATRIX_DIAGONAL_DATA(m)->x,gdouble,2*i+1) ;
    break ;    
  }  

  return x ;
}

/** 
 * Compute the inverse of a matrix in-place, overwriting the original.
 * 
 * @param A a ::sisl_matrix_t which is replaced with its inverse on exit. 
 * 
 * @return 0 on success.
 */

gint sisl_matrix_invert(sisl_matrix_t *A)

{
  /*stopgap until I implement some proper workspaces*/
  gdouble work[4096] ;
  gint ip[4096], info, nr, nc, lwork = 4096 ;
  
  SISL_CHECK_ARGUMENT_NULL(A) ;

  if ( (nr = sisl_matrix_row_number(A)) != (nc = sisl_matrix_column_number(A)) )
    g_error("%s: matrix (%dx%d) must be square to compute matrix",
	    __FUNCTION__, nr, nc) ;

  if ( wmpi_process_number() != 1 )
    g_error("%s: only implemented for single-process calculations",
	    __FUNCTION__) ;

  g_assert(sisl_matrix_row_number(A) <= 4096 ) ;
  if ( sisl_matrix_density(A) == SISL_MATRIX_SPARSE )
    g_error("%s: cannot invert sparse matrix", __FUNCTION__) ;

  if ( sisl_matrix_density(A) == SISL_MATRIX_DENSE ) {
    if ( sisl_is_real(A) ) {
      dgetrf_(&nr, &nc, (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data),
	      &nc, ip, &info) ;
      dgetri_(&nr, (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data), &nr,
	      ip, work, &lwork, &info) ;
      return info ;
    }

    lwork /= 2 ;
    zgetrf_(&nr, &nc, (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data),
	    &nc, ip, &info) ;
    zgetri_(&nr, (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data), &nr,
	    ip, work, &lwork, &info) ;
    
    return info ;
  }

  if ( sisl_matrix_density(A) == SISL_MATRIX_DIAGONAL ) {
    if ( sisl_is_real(A) ) {
      sisl_diagonal_invert(SISL_MATRIX_DIAGONAL_DATA(A)) ;
      return 0 ;
    }

    sisl_diagonal_invert_complex(SISL_MATRIX_DIAGONAL_DATA(A)) ;
    return 0 ;
  }
  
  g_assert_not_reached() ;
  return 0 ;
}

/** 
 * Copy one matrix into another, \f$A=B\f$. Matrices must be both
 * real, or both complex. If \a A is dense, \a B can have arbitrary
 * density; if \a B is not dense, \a A must have the same density
 * (sparse, diagonal, etc.) as \a B.
 * 
 * @param A destination ::sisl_matrix_t
 * @param B source ::sisl_matrix_t
 * 
 * @return 0 on success.
 */

gint sisl_matrix_copy(sisl_matrix_t *A, sisl_matrix_t *B)

{
  GArray *dA, *dB ;
  gint nr, nc, i, k ;
  gdouble *dat ;
  sisl_matrix_sparse_t *Bs ;
  sisl_matrix_diagonal_t *Bd ;
  
  SISL_CHECK_ARGUMENT_NULL(A) ;
  SISL_CHECK_ARGUMENT_NULL(B) ;

  /* A->d = B->d ; */
  A->c = B->c ; A->dist = B->dist ;

  sisl_matrix_row_number(A) = sisl_matrix_row_number(B) ;
  sisl_matrix_column_number(A) = sisl_matrix_column_number(B) ;
  sisl_matrix_local_row_start(A) = sisl_matrix_local_row_start(B) ;
  sisl_matrix_local_row_end(A) = sisl_matrix_local_row_end(B) ;
  sisl_matrix_distribution(A) = sisl_matrix_distribution(B) ;

  if ( A->d == B->d ) {
    switch ( A->d ) {
    default: 
      g_error("%s: unrecognized matrix density type %d", __FUNCTION__, A->d) ;
    break ;
    case SISL_MATRIX_USER_DEFINED:
      g_error("%s: cannot copy user-defined matrix type", __FUNCTION__) ;
      break ;
    case SISL_MATRIX_DENSE:
      dB = ((sisl_matrix_dense_t *)(B->m))->x ;
      dA = ((sisl_matrix_dense_t *)(A->m))->x ;
      g_array_set_size(dA, dB->len) ;
      memcpy(dA->data, dB->data, (dB->len)*sizeof(gdouble)) ;
      break ;
    case SISL_MATRIX_SPARSE:
      g_assert_not_reached() ; /*untested code*/
      dB = ((sisl_matrix_sparse_t *)(B->m))->x ;
      dA = ((sisl_matrix_sparse_t *)(A->m))->x ;
      g_array_set_size(dA, dB->len) ;
      memcpy(dA->data, dB->data, (dB->len)*sizeof(gdouble)) ;
      
      dB = ((sisl_matrix_sparse_t *)(B->m))->ij ;
      dA = ((sisl_matrix_sparse_t *)(A->m))->ij ; 
      g_array_set_size(dA, dB->len) ;
      memcpy(dA->data, dB->data, (dB->len)*sizeof(gint)) ;
      break ;
    case SISL_MATRIX_DIAGONAL:
      dB = ((sisl_matrix_diagonal_t *)(B->m))->x ;
      dA = ((sisl_matrix_diagonal_t *)(A->m))->x ;
      g_array_set_size(dA, dB->len) ;
      memcpy(dA->data, dB->data, (dB->len)*sizeof(gdouble)) ;
      break ;
    }

    return 0 ;
  }

  g_assert( sisl_matrix_distribution(A) == SISL_SINGLE ) ;

  if ( sisl_matrix_density(A) != SISL_MATRIX_DENSE )
    g_error("%s: mixed density copies must have dense destination matrix",
	    __FUNCTION__) ;

  nc = sisl_matrix_column_number(A) ;
  nr = sisl_matrix_row_number(A) ;

  dA = ((sisl_matrix_dense_t *)(A->m))->x ;
  if ( sisl_is_real(A) ) {
    g_array_set_size(dA, nr*nc) ;
    memset(dA->data, 0, nr*nc*sizeof(gdouble)) ;
    dat = (gdouble *)(dA->data) ;
    switch ( B->d ) {
    default: 
      g_error("%s: unrecognized matrix density type %d", __FUNCTION__, B->d) ;
      break ;
    case SISL_MATRIX_SPARSE:
      Bs = SISL_MATRIX_SPARSE_DATA(B) ;
      for ( k = 0 ; k < Bs->ij->len ; k ++ ) {
	i = g_array_index(Bs->ij,gint,k) ;
	dat[i] += g_array_index(Bs->x,gdouble,k) ;
      }
      break ;
    case SISL_MATRIX_DIAGONAL:
      Bd = SISL_MATRIX_DIAGONAL_DATA(B) ;
      for ( k = 0 ; k < Bd->x->len ; k ++ ) {
	dat[k*nc+k] = g_array_index(Bd->x,gdouble,k) ;
      }
      break ;
    }

    return 0 ;
  }

  g_array_set_size(dA, nr*nc*2) ;
  memset(dA->data, 0, 2*nr*nc*sizeof(gdouble)) ;
  dat = (gdouble *)(dA->data) ;

  switch ( B->d ) {
  default: 
    g_error("%s: unrecognized matrix density type %d", __FUNCTION__, B->d) ;
    break ;
  case SISL_MATRIX_SPARSE:
    Bs = SISL_MATRIX_SPARSE_DATA(B) ;
    for ( k = 0 ; k < Bs->ij->len ; k ++ ) {
      i = g_array_index(Bs->ij,gint,k) ;
      dat[2*i+0] += g_array_index(Bs->x,gdouble,2*k+0) ;
      dat[2*i+1] += g_array_index(Bs->x,gdouble,2*k+1) ;
    }
    break ;
  case SISL_MATRIX_DIAGONAL:
    Bd = SISL_MATRIX_DIAGONAL_DATA(B) ;
    for ( k = 0 ; k < Bd->x->len/2 ; k ++ ) {
      dat[2*(k*nc+k)+0] = g_array_index(Bd->x,gdouble,2*k+0) ;
      dat[2*(k*nc+k)+1] = g_array_index(Bd->x,gdouble,2*k+1) ;
    }
    break ;
  }

  return 0 ;
}

/** 
 * Matrix-matrix multiplication, \f$C=AB\f$
 * 
 * @param A a ::sisl_matrix_t
 * @param B a ::sisl_matrix_t
 * @param C a ::sisl_matrix_t
 * 
 * @return 0 on success.
 */

gint sisl_matrix_matrix_mul(sisl_matrix_t *A, sisl_matrix_t *B,
			    sisl_matrix_t *C)

{
  gint nra, nca, nrb, ncb ;
  
  SISL_CHECK_ARGUMENT_NULL(A) ;
  SISL_CHECK_ARGUMENT_NULL(B) ;
  SISL_CHECK_ARGUMENT_NULL(C) ;

  nra = sisl_matrix_row_number(A) ;
  nca = sisl_matrix_column_number(A) ;
  nrb = sisl_matrix_row_number(B) ;
  ncb = sisl_matrix_column_number(B) ;

  if ( sisl_is_real(A) && !sisl_is_real(B) )
    g_error("%s: cannot multiply real matrix A and complex matrix B",
	    __FUNCTION__) ;

  if ( !sisl_is_real(A) && sisl_is_real(B) )
    g_error("%s: cannot multiply complex matrix A and real matrix B",
	    __FUNCTION__) ;
  
  if ( nca != nrb )
    g_error("%s: cannot multiply %dx%d and %dx%d matrix",
	    __FUNCTION__, nra, nca, nrb, ncb) ;
  
  if ( wmpi_process_number() != 1 )
    g_error("%s: only implemented for single-process calculations",
	    __FUNCTION__) ;

  /*for multiplication using the Fortran ordering, swap the matrices:
   C^{T} = B^{T}xA^{T}*/
  if ( sisl_matrix_density(A) == SISL_MATRIX_DENSE &&
       sisl_matrix_density(B) == SISL_MATRIX_DENSE ) {
    sisl_matrix_set_block_size(C, nra, ncb) ;
    sisl_matrix_row_number(C) = nra ; 
    sisl_matrix_column_number(C) = ncb ;
    sisl_matrix_local_row_start(C) = 0 ;
    sisl_matrix_local_row_end(C) = nra ;
    if ( sisl_is_real(A) ) 
      dgemm_("N", "N", &ncb, &nra, &nrb, &d_1,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(B)->x->data), &ncb,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data), &nca,
	     &d_0,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(C)->x->data), &ncb) ;
    else
      zgemm_("N", "N", &ncb, &nra, &nrb, one,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(B)->x->data), &ncb,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data), &nca,
	     zero,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(C)->x->data), &ncb) ;      

    return 0 ;
  }
  
  g_assert_not_reached() ;
  
  return 0 ;
}

/** 
 * Weighted matrix-matrix multiplication, \f$C=aAB + cC\f$,
 * corresponding to BLAS dgemm. Matrices may be real or complex;
 * weights are real. 
 * 
 * @param A a ::sisl_matrix_t
 * @param B a ::sisl_matrix_t
 * @param C a ::sisl_matrix_t
 * @param a real weight
 * @param c real weight
 * 
 * @return 0 on success.
 */

gint sisl_matrix_matrix_mul_w(sisl_matrix_t *A, sisl_matrix_t *B,
			      sisl_matrix_t *C, gdouble a, gdouble c)

{
  gint nra, nca, nrb, ncb ;
  gdouble alpha[2] = {0.0}, beta[2] = {0.0};
  
  SISL_CHECK_ARGUMENT_NULL(A) ;
  SISL_CHECK_ARGUMENT_NULL(B) ;
  SISL_CHECK_ARGUMENT_NULL(C) ;

  nra = sisl_matrix_row_number(A) ;
  nca = sisl_matrix_column_number(A) ;
  nrb = sisl_matrix_row_number(B) ;
  ncb = sisl_matrix_column_number(B) ;

  if ( sisl_is_real(A) && !sisl_is_real(B) )
    g_error("%s: cannot multiply real matrix A and complex matrix B",
	    __FUNCTION__) ;

  if ( !sisl_is_real(A) && sisl_is_real(B) )
    g_error("%s: cannot multiply complex matrix A and real matrix B",
	    __FUNCTION__) ;
  
  if ( nca != nrb )
    g_error("%s: cannot multiply %dx%d and %dx%d matrix",
	    __FUNCTION__, nra, nca, nrb, ncb) ;
  
  if ( wmpi_process_number() != 1 )
    g_error("%s: only implemented for single-process calculations",
	    __FUNCTION__) ;

  /*set up the Xgemm weights*/
  alpha[0] = a ; beta[0] = c ; 
  
  /*for multiplication using the Fortran ordering, swap the matrices:
   C^{T} = B^{T}xA^{T}*/
  if ( sisl_matrix_density(A) == SISL_MATRIX_DENSE &&
       sisl_matrix_density(B) == SISL_MATRIX_DENSE ) {
    sisl_matrix_set_block_size(C, nra, ncb) ;
    sisl_matrix_row_number(C) = nra ; 
    sisl_matrix_column_number(C) = ncb ;
    sisl_matrix_local_row_start(C) = 0 ;
    sisl_matrix_local_row_end(C) = nra ;
    if ( sisl_is_real(A) ) 
      dgemm_("N", "N", &ncb, &nra, &nrb, alpha,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(B)->x->data), &ncb,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data), &nca,
	     beta,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(C)->x->data), &ncb) ;
    else
      zgemm_("N", "N", &ncb, &nra, &nrb, alpha,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(B)->x->data), &ncb,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data), &nca,
	     beta,
	     (gdouble *)(SISL_MATRIX_DENSE_DATA(C)->x->data), &ncb) ;      

    return 0 ;
  }

  if ( sisl_matrix_density(A) == SISL_MATRIX_DENSE &&
       sisl_matrix_density(B) == SISL_MATRIX_SPARSE ) {
    if ( sisl_matrix_density(C) != SISL_MATRIX_DENSE )
      g_error("%s: output matrix must be dense for (dense)x(sparse) "
	      "multiplication", __FUNCTION__) ;
    sisl_matrix_set_block_size(C, nra, ncb) ;
    sisl_matrix_row_number(C) = nra ; 
    sisl_matrix_column_number(C) = ncb ;
    sisl_matrix_local_row_start(C) = 0 ;
    sisl_matrix_local_row_end(C) = nra ;
    if ( sisl_is_real(A) ) {
      sisl_dense_sparse_mul_w(SISL_MATRIX_DENSE_DATA(A),
			      SISL_MATRIX_SPARSE_DATA(B),
			      a, c,
			      SISL_MATRIX_DENSE_DATA(C)) ;
    } else {
      sisl_dense_sparse_mul_w_complex(SISL_MATRIX_DENSE_DATA(A),
				      SISL_MATRIX_SPARSE_DATA(B),
				      *((gsl_complex *)alpha), 
				      *((gsl_complex *)beta),
				      SISL_MATRIX_DENSE_DATA(C)) ;
    }

    return 0 ;
  }
  
  g_assert_not_reached() ;
  
  return 0 ;
}

/** 
 * Weighted complex matrix-matrix multiplication, \f$C=aAB + cC\f$,
 * corresponding to BLAS dgemm. \a A, \a B, and \a C must all be
 * complex.
 * 
 * @param A a ::sisl_matrix_t
 * @param B a ::sisl_matrix_t
 * @param C a ::sisl_matrix_t
 * @param a complex weight
 * @param c complex weight
 * 
 * @return 0 on success.
 */

gint sisl_matrix_matrix_mul_w_complex(sisl_matrix_t *A, sisl_matrix_t *B,
				      sisl_matrix_t *C,
				      gsl_complex a, gsl_complex c)

{
  gint nra, nca, nrb, ncb ;
  
  SISL_CHECK_ARGUMENT_NULL(A) ;
  SISL_CHECK_ARGUMENT_NULL(B) ;
  SISL_CHECK_ARGUMENT_NULL(C) ;

  nra = sisl_matrix_row_number(A) ;
  nca = sisl_matrix_column_number(A) ;
  nrb = sisl_matrix_row_number(B) ;
  ncb = sisl_matrix_column_number(B) ;

  if ( sisl_is_real(A) )
    g_error("%s: matrix A must be complex", __FUNCTION__) ;
  if ( sisl_is_real(B) )
    g_error("%s: matrix B must be complex", __FUNCTION__) ;
  if ( sisl_is_real(C) )
    g_error("%s: matrix C must be complex", __FUNCTION__) ;
  
  if ( nca != nrb )
    g_error("%s: cannot multiply %dx%d and %dx%d matrix",
	    __FUNCTION__, nra, nca, nrb, ncb) ;
  
  if ( wmpi_process_number() != 1 )
    g_error("%s: only implemented for single-process calculations",
	    __FUNCTION__) ;

  /*for multiplication using the Fortran ordering, swap the matrices:
   C^{T} = B^{T}xA^{T}*/
  if ( sisl_matrix_density(A) == SISL_MATRIX_DENSE &&
       sisl_matrix_density(B) == SISL_MATRIX_DENSE ) {
    sisl_matrix_set_block_size(C, nra, ncb) ;
    sisl_matrix_row_number(C) = nra ; 
    sisl_matrix_column_number(C) = ncb ;
    sisl_matrix_local_row_start(C) = 0 ;
    sisl_matrix_local_row_end(C) = nra ;
    zgemm_("N", "N", &ncb, &nra, &nrb, a.dat,
	   (gdouble *)(SISL_MATRIX_DENSE_DATA(B)->x->data), &ncb,
	   (gdouble *)(SISL_MATRIX_DENSE_DATA(A)->x->data), &nca,
	   c.dat,
	   (gdouble *)(SISL_MATRIX_DENSE_DATA(C)->x->data), &ncb) ;      

    return 0 ;
  }

  if ( sisl_matrix_density(A) == SISL_MATRIX_DENSE &&
       sisl_matrix_density(B) == SISL_MATRIX_SPARSE ) {
    if ( sisl_matrix_density(C) != SISL_MATRIX_DENSE )
      g_error("%s: output matrix must be dense for (dense)x(sparse) "
	      "multiplication", __FUNCTION__) ;
    sisl_matrix_set_block_size(C, nra, ncb) ;
    sisl_matrix_row_number(C) = nra ; 
    sisl_matrix_column_number(C) = ncb ;
    sisl_matrix_local_row_start(C) = 0 ;
    sisl_matrix_local_row_end(C) = nra ;

    sisl_dense_sparse_mul_w_complex(SISL_MATRIX_DENSE_DATA(A),
				    SISL_MATRIX_SPARSE_DATA(B),
				    a, c,
				    SISL_MATRIX_DENSE_DATA(C)) ;

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

/** 
 * Weighted triple matrix multiplication, \f$D=aABC + dD\f$, currently
 * only implemented for: A and C dense, B diagonal; A dense, B and C
 * diagonal; A dense, B diagonal, C sparse.
 *
 * @param A a ::sisl_matrix_t
 * @param B a ::sisl_matrix_t
 * @param C a ::sisl_matrix_t
 * @param D a ::sisl_matrix_t
 * @param a weight
 * @param d weight
 * 
 * @return 0 on success.
 */

gint sisl_matrix_triple_mul_w(sisl_matrix_t *A, sisl_matrix_t *B,
			      sisl_matrix_t *C, sisl_matrix_t *D,
			      gdouble a, gdouble d)
{
  gint nra, nca, ncc, n ;
  gdouble *dA, *dB, *dC, *dD, xt ;
  
  if ( wmpi_process_number() != 1 )
    g_error("%s: only implemented for single-process calculations",
	    __FUNCTION__) ;

  if ( !sisl_is_real(A) )
    g_error("%s: matrix A must be real", __FUNCTION__) ;

  if ( !sisl_is_real(B) )
    g_error("%s: matrix B must be real", __FUNCTION__) ;

  if ( !sisl_is_real(C) )
    g_error("%s: matrix C must be real", __FUNCTION__) ;

  if ( !sisl_is_real(D) )
    g_error("%s: matrix D must be real", __FUNCTION__) ;

  if ( sisl_matrix_density(A) != SISL_MATRIX_DENSE )
    g_error("%s: matrix A must be dense", __FUNCTION__) ;

  if ( sisl_matrix_density(A) == SISL_MATRIX_SPARSE ||
       sisl_matrix_density(B) == SISL_MATRIX_SPARSE )
    g_error("%s: not implemented for A or B sparse", __FUNCTION__) ;

  if ( sisl_matrix_density(B) != SISL_MATRIX_DIAGONAL )
    g_error("%s: matrix B must be diagonal", __FUNCTION__) ;

  if ( sisl_matrix_column_number(A) != sisl_matrix_row_number(B) )
    g_error("%s: matrices A (%dx%d) and B (%dx%d) cannot be multiplied",
	    __FUNCTION__,
	    sisl_matrix_row_number(A), sisl_matrix_column_number(A),
	    sisl_matrix_row_number(B), sisl_matrix_column_number(B)) ;
  if ( sisl_matrix_column_number(B) != sisl_matrix_row_number(C) )
    g_error("%s: matrices B (%dx%d) and C (%dx%d) cannot be multiplied",
	    __FUNCTION__,
	    sisl_matrix_row_number(B), sisl_matrix_column_number(B),
	    sisl_matrix_row_number(C), sisl_matrix_column_number(C)) ;

  if ( d == 0.0 ) {
    sisl_matrix_set_block_size(D,
			       sisl_matrix_row_number(A),
			       sisl_matrix_column_number(C)) ;
  } else {
    if ( sisl_matrix_row_number(D) != sisl_matrix_row_number(A) ||
	 sisl_matrix_column_number(D) != sisl_matrix_column_number(C) ) {
      g_error("%s: cannot perform weighted sum of (%dx%d) matrix "
	      "and (%dx%d)x(%dx%d)(%dx%d) product",
	      __FUNCTION__,
	      sisl_matrix_row_number(D), sisl_matrix_column_number(D), 
	      sisl_matrix_row_number(A), sisl_matrix_column_number(A), 
	      sisl_matrix_row_number(B), sisl_matrix_column_number(B), 
	      sisl_matrix_row_number(C), sisl_matrix_column_number(C)
	      ) ;
    }
  }
  
  dA = &(g_array_index(SISL_MATRIX_DENSE_DATA(A)->x, gdouble, 0)) ;
  dD = &(g_array_index(SISL_MATRIX_DENSE_DATA(D)->x, gdouble, 0)) ;
  nra = sisl_matrix_row_number(A) ;
  nca = sisl_matrix_column_number(A) ;
  ncc = sisl_matrix_column_number(C) ;

  for ( n = 0 ; n < nra*ncc ; n ++ ) dD[n] *= d ;
  
  if ( sisl_matrix_density(B) == SISL_MATRIX_DIAGONAL &&
       sisl_matrix_density(C) == SISL_MATRIX_DENSE ) {
    /*(dense)*(diagonal)*(dense) with square matrices*/
    gint i, j, k ;
    
    dB = &(g_array_index(SISL_MATRIX_DIAGONAL_DATA(B)->x, gdouble, 0)) ;
    dC = &(g_array_index(SISL_MATRIX_DENSE_DATA(C)->x, gdouble, 0)) ;
    
    /*this should be done in an efficient way at some point, but let's
      play safe for now*/
    for ( i = 0 ; i < nra ; i ++ ) {
      for ( j = 0 ; j < ncc ; j ++ ) {
	for ( k = 0 ; k < nca ; k ++ )
	  dD[i*ncc+j] += a*dA[i*nca+k]*dB[k]*dC[k*ncc+j] ;
      }
    }

    return 0 ;
  }

  if ( sisl_matrix_density(B) == SISL_MATRIX_DIAGONAL &&
       sisl_matrix_density(C) == SISL_MATRIX_DIAGONAL ) {
    /*(dense)*(diagonal)*(diagonal) with square matrices*/
    gint i, j ;
    
    dB = &(g_array_index(SISL_MATRIX_DIAGONAL_DATA(B)->x, gdouble, 0)) ;
    dC = &(g_array_index(SISL_MATRIX_DIAGONAL_DATA(C)->x, gdouble, 0)) ;

    for ( i = 0 ; i < nra ; i ++ ) {
      for ( j = 0 ; j < nca ; j ++ ) {
    	dD[i*nca+j] += a*dA[i*nca+j]*dB[j]*dC[j] ;
      }
    }

    return 0 ;
  }

  if ( sisl_matrix_density(B) == SISL_MATRIX_DIAGONAL &&
       sisl_matrix_density(C) == SISL_MATRIX_SPARSE ) {
    /*(dense)*(diagonal)*(sparse) with square matrices*/
    gint i, j, k ;
    sisl_matrix_sparse_t *Cs ;

    dB = &(g_array_index(SISL_MATRIX_DIAGONAL_DATA(B)->x, gdouble, 0)) ;
    Cs = SISL_MATRIX_SPARSE_DATA(C) ;

    ncc = sisl_matrix_column_number(C) ;
    nca = sisl_matrix_column_number(A) ;

    for ( n = 0 ; n < Cs->ij->len ; n ++ ) {
      k = g_array_index(Cs->ij,gint,n)/ncc ;
      j = g_array_index(Cs->ij,gint,n) % ncc ;
      xt = dB[k]*g_array_index(Cs->x,gdouble,n) ;
      for ( i = 0 ; i < nra ; i ++ ) {
    	dD[i*ncc+j] += a*dA[i*nca+k]*xt ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

/** 
 * Complex weighted triple matrix multiplication, \f$D=aABC + dD\f$,
 *
 * @param A a ::sisl_matrix_t
 * @param B a ::sisl_matrix_t
 * @param C a ::sisl_matrix_t
 * @param D a ::sisl_matrix_t
 * @param a complex weight
 * @param d complex weight
 * 
 * @return 0 on success.
 */

gint sisl_matrix_triple_mul_w_complex(sisl_matrix_t *A, sisl_matrix_t *B,
				      sisl_matrix_t *C, sisl_matrix_t *D,
				      gsl_complex a, gsl_complex d)
{
  gint nra, nca, ncc, n ;
  gsl_complex *dA, *dB, *dC, *dD, zt ;
  
  if ( wmpi_process_number() != 1 )
    g_error("%s: only implemented for single-process calculations",
	    __FUNCTION__) ;

  if ( sisl_is_real(A) )
    g_error("%s: matrix A must be complex", __FUNCTION__) ;

  if ( sisl_is_real(B) )
    g_error("%s: matrix B must be complex", __FUNCTION__) ;

  if ( sisl_is_real(C) )
    g_error("%s: matrix C must be complex", __FUNCTION__) ;

  if ( sisl_is_real(D) )
    g_error("%s: matrix D must be complex", __FUNCTION__) ;

  if ( sisl_matrix_density(A) != SISL_MATRIX_DENSE )
    g_error("%s: matrix A must be dense", __FUNCTION__) ;

  if ( sisl_matrix_density(A) == SISL_MATRIX_SPARSE ||
       sisl_matrix_density(B) == SISL_MATRIX_SPARSE )
       /* sisl_matrix_density(C) == SISL_MATRIX_SPARSE ) */
    g_error("%s: not implemented for A or B sparse", __FUNCTION__) ;
  
  if ( sisl_matrix_density(B) != SISL_MATRIX_DIAGONAL )
    g_error("%s: matrix B must be diagonal", __FUNCTION__) ;

  if ( sisl_matrix_column_number(A) != sisl_matrix_row_number(B) )
    g_error("%s: matrices A (%dx%d) and B (%dx%d) cannot be multiplied",
	    __FUNCTION__,
	    sisl_matrix_row_number(A), sisl_matrix_column_number(A),
	    sisl_matrix_row_number(B), sisl_matrix_column_number(B)) ;
  if ( sisl_matrix_column_number(B) != sisl_matrix_row_number(C) )
    g_error("%s: matrices B (%dx%d) and C (%dx%d) cannot be multiplied",
	    __FUNCTION__,
	    sisl_matrix_row_number(B), sisl_matrix_column_number(B),
	    sisl_matrix_row_number(C), sisl_matrix_column_number(C)) ;

  if ( GSL_REAL(d) == 0.0 && GSL_IMAG(d) == 0.0 ) {
    sisl_matrix_set_block_size(D,
			       sisl_matrix_row_number(A),
			       sisl_matrix_column_number(C)) ;
  } else {
    if ( sisl_matrix_row_number(D) != sisl_matrix_row_number(A) ||
	 sisl_matrix_column_number(D) != sisl_matrix_column_number(C) ) {
      g_error("%s: cannot perform weighted sum of (%dx%d) matrix "
	      "and (%dx%d)x(%dx%d)(%dx%d) product",
	      __FUNCTION__,
	      sisl_matrix_row_number(D), sisl_matrix_column_number(D), 
	      sisl_matrix_row_number(A), sisl_matrix_column_number(A), 
	      sisl_matrix_row_number(B), sisl_matrix_column_number(B), 
	      sisl_matrix_row_number(C), sisl_matrix_column_number(C)
	      ) ;
    }
  }
    
  dA = (gsl_complex *)
    (&(g_array_index(SISL_MATRIX_DENSE_DATA(A)->x, gdouble, 0))) ;
  dD = (gsl_complex *)
    (&(g_array_index(SISL_MATRIX_DENSE_DATA(D)->x, gdouble, 0))) ;
  nra = sisl_matrix_row_number(A) ;
  nca = sisl_matrix_column_number(A) ;
  ncc = sisl_matrix_column_number(C) ;

  for ( n = 0 ; n < nra*ncc ; n ++ ) dD[n] = gsl_complex_mul(dD[n], d) ;    

  if ( sisl_matrix_density(B) == SISL_MATRIX_DIAGONAL &&
       sisl_matrix_density(C) == SISL_MATRIX_DENSE ) {
    /*(dense)*(diagonal)*(dense)*/
    gint i, j, k ;
    
    dB = (gsl_complex *)
      (&(g_array_index(SISL_MATRIX_DIAGONAL_DATA(B)->x, gdouble, 0))) ;
    dC = (gsl_complex *)
      (&(g_array_index(SISL_MATRIX_DENSE_DATA(C)->x, gdouble, 0))) ;

    /*this should be done in an efficient way at some point, but let's
      play safe for now*/
    for ( i = 0 ; i < nra ; i ++ ) {
      for ( j = 0 ; j < ncc ; j ++ ) {
	for ( k = 0 ; k < nca ; k ++ ) {
	  zt = gsl_complex_mul(a, dA[i*nca+k]) ;
	  zt = gsl_complex_mul(zt, dB[k]) ;
	  zt = gsl_complex_mul(zt, dC[k*ncc+j]) ;
	  
	  dD[i*ncc+j] = gsl_complex_add(dD[i*ncc+j], zt) ;
	}
      }
    }

    return 0 ;
  }

  if ( sisl_matrix_density(B) == SISL_MATRIX_DIAGONAL &&
       sisl_matrix_density(C) == SISL_MATRIX_DIAGONAL ) {
    /*(dense)*(diagonal)*(diagonal)*/
    gint i, j ;
    
    dB = (gsl_complex *)
      (&(g_array_index(SISL_MATRIX_DIAGONAL_DATA(B)->x, gdouble, 0))) ;
    dC = (gsl_complex *)
      (&(g_array_index(SISL_MATRIX_DIAGONAL_DATA(C)->x, gdouble, 0))) ;

    for ( i = 0 ; i < nra ; i ++ ) {
      for ( j = 0 ; j < nca ; j ++ ) {
	zt = gsl_complex_mul(a, dA[i*nca+j]) ;
	zt = gsl_complex_mul(zt, dB[j]) ;
	zt = gsl_complex_mul(zt, dC[j]) ;
	
	dD[i*nca+j] = gsl_complex_add(dD[i*nca+j], zt) ;
      }
    }

    return 0 ;
  }
  
  if ( sisl_matrix_density(B) == SISL_MATRIX_DIAGONAL &&
       sisl_matrix_density(C) == SISL_MATRIX_SPARSE ) {
    /*(dense)*(diagonal)*(sparse)*/
    gint i, j, k ;
    GArray *ij ;
    
    dB = (gsl_complex *)
      (&(g_array_index(SISL_MATRIX_DIAGONAL_DATA(B)->x, gdouble, 0))) ;
    dC = (gsl_complex *)
      (&(g_array_index(SISL_MATRIX_SPARSE_DATA(C)->x, gdouble, 0))) ;
    ij = SISL_MATRIX_SPARSE_DATA(C)->ij ;

    ncc = sisl_matrix_column_number(C) ;
    nca = sisl_matrix_column_number(A) ;
    
    for ( n = 0 ; n < ij->len ; n ++ ) {
      k = g_array_index(ij,gint,n)/ncc ;
      j = g_array_index(ij,gint,n) % ncc ;
      zt = gsl_complex_mul(a, gsl_complex_mul(dB[k], dC[n])) ;
      for ( i = 0 ; i < nra ; i ++ ) {
	dD[i*ncc+j] = gsl_complex_add(dD[i*ncc+j],
				      gsl_complex_mul(dA[i*nca+k], zt)) ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

/** 
 * Generate a matrix with random entries, respecting the existing
 * distribution and complex number settings. The input matrix must be
 * dense or diagonal.
 * 
 * @param A a ::sisl_matrix_t
 * @param nr number of rows
 * @param nc number of columns
 * 
 * @return 0 on success.
 */

gint sisl_matrix_random(sisl_matrix_t *A, gint nr, gint nc)
  
{
  gint i, j, k ;
  gsl_complex z ;
  GArray *x ;
  
  if ( wmpi_process_number() != 1 )
    g_error("%s: only implemented for single-process calculations",
	    __FUNCTION__) ;

  if ( sisl_matrix_density(A) == SISL_MATRIX_DENSE ) {
    sisl_matrix_set_block_size(A, nr, nc) ;
    sisl_matrix_row_number(A) = nr ; 
    sisl_matrix_column_number(A) = nc ;
    sisl_matrix_local_row_start(A) = 0 ;
    sisl_matrix_local_row_end(A) = nr ;
    
    x = SISL_MATRIX_DENSE_DATA(A)->x ;
    for ( i = 0 ; i < x->len ; i ++ )
      g_array_index(x, gdouble, i) = g_random_double() ;

    return 0 ;
  }

  if ( sisl_matrix_density(A) == SISL_MATRIX_DIAGONAL ) {
    sisl_matrix_set_block_size(A, nr, nr) ;
    sisl_matrix_row_number(A) = nr ; 
    sisl_matrix_column_number(A) = nr ;
    sisl_matrix_local_row_start(A) = 0 ;
    sisl_matrix_local_row_end(A) = nr ;
    
    x = SISL_MATRIX_DIAGONAL_DATA(A)->x ;
    if ( sisl_is_real(A) ) g_array_set_size(x, nr) ;
    else g_array_set_size(x, 2*nr) ;

    for ( i = 0 ; i < x->len ; i ++ )
      g_array_index(x, gdouble, i) = g_random_double() ;

    return 0 ;
  }

  if ( sisl_matrix_density(A) == SISL_MATRIX_SPARSE ) {
    /*fill about 25% of the entries in the matrix*/
    sisl_matrix_row_number(A) = nr ; 
    sisl_matrix_column_number(A) = nc ;
    sisl_matrix_local_row_start(A) = 0 ;
    sisl_matrix_local_row_end(A) = nr ;

    if ( sisl_is_real(A) ) {
      for ( k = 0 ; k < nr*nc/4 ; k ++ ) {
	i = g_random_int_range(0, nr) ;
	j = g_random_int_range(0, nc) ;
	sisl_matrix_set(A, i, j, g_random_double()) ;
      }
    } else {
      for ( k = 0 ; k < nr*nc/4 ; k ++ ) {
	i = g_random_int_range(0, nr) ;
	j = g_random_int_range(0, nc) ;
	GSL_SET_COMPLEX(&z, g_random_double(), g_random_double()) ;
	sisl_matrix_set_complex(A, i, j, z) ;
      }      
    }
    return 0 ;
  }
    
    
  g_error("%s: cannot generate random matrix for matrix density %d",
	  __FUNCTION__, sisl_matrix_density(A)) ;
  return 0 ;
}

/** 
 * Set all entries of a matrix to the same value. The input matrix
 * must be diagonal (for now). If the matrix is complex, its entries
 * are set to the assigned value with the imaginary part set to zero.
 * 
 * @param m a ::sisl_matrix_t
 * @param x value to assign to entries of \a m
 * 
 * @return 0 on success.
 */

gint sisl_matrix_set_all(sisl_matrix_t *m, gdouble x)

{
  gint i ;
  gdouble *dm ;
  
  switch ( m->d ) {
  default: 
    g_error("%s: unrecognized matrix density type %d", __FUNCTION__, m->d) ;
    break ;
  case SISL_MATRIX_DIAGONAL:
    dm = &(g_array_index(SISL_MATRIX_DIAGONAL_DATA(m)->x, gdouble, 0)) ;

    if ( sisl_is_real(m) )    
      for ( i = 0 ; i < sisl_matrix_row_number(m) ; i ++ ) dm[i] = x ;
    else
      for ( i = 0 ; i < sisl_matrix_row_number(m) ; i ++ ) {
	dm[2*i+0] = x ; dm[2*i+1] = 0.0 ;
      }
    
    break ;
  }  

  return 0 ;
}

/** 
 * Set all entries of a complex matrix to the same complex value. The
 * input matrix must be diagonal (for now).
 * 
 * @param m a ::sisl_matrix_t
 * @param x value to assign to entries of \a m
 * 
 * @return 0 on success.
 */

gint sisl_matrix_set_all_complex(sisl_matrix_t *m, gsl_complex x)

{
  gint i ;
  gdouble *dm, xr, xi ;

  if ( sisl_is_real(m) )
    g_error("%s: cannot assign a complex entry to a real matrix",
	    __FUNCTION__) ;
  
  xr = GSL_REAL(x) ; xi = GSL_IMAG(x) ; 
  switch ( m->d ) {
  default: 
    g_error("%s: unrecognized matrix density type %d", __FUNCTION__, m->d) ;
    break ;
  case SISL_MATRIX_DIAGONAL:
    dm = &(g_array_index(SISL_MATRIX_DIAGONAL_DATA(m)->x, gdouble, 0)) ;
    
    for ( i = 0 ; i < sisl_matrix_row_number(m) ; i ++ ) {
      dm[2*i+0] = xr ; dm[2*i+1] = xi ;
    }
    
    break ;
  }  

  return 0 ;
}

/**
 * @}
 * 
 */

