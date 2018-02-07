/* matrix.c
 * 
 * Copyright (C) 2009 Michael Carley
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <sisl.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "sisl_private.h"

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
      g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
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
  /* gdouble alpha = 1.0, beta = 0.0 ; */

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
      sisl_sparse_mul_c(SISL_MATRIX_SPARSE_DATA(m), 
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
  }  

  return x ;
}

/**
 * @}
 * 
 */

