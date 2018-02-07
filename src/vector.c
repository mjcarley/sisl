/* vector.c
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

#ifdef HAVE_STRING_H
#include <string.h>
#endif /*HAVE_STRING_H*/

#include <sisl.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <wmpi.h>

#include "sisl_private.h"

/**
 * @defgroup vector SISL vectors
 *
 * All vectors in SISL are dense and can be real or complex. Vector
 * operations are defined for real and complex vectors (with function
 * suffix _complex) and functions will fail with an error message if
 * passed an incompatible parameter. The basic data type is
 * ::sisl_vector_t which should only be accessed using the supplied
 * functions and macros.
 *
 * @{
 */

static gdouble *buffer_double(gint n)

{
  static gdouble *b = NULL ;
  static gint len = 0 ;

  if ( n > len ) {
    len = n ;
    b = (gdouble *)g_realloc(b,len*sizeof(gdouble)) ;
    g_debug("%s: setting buffer length to %d", __FUNCTION__, len) ;
  }

  return b ;
}

/** 
 * Allocate a new ::sisl_vector_t.
 * 
 * @param c ::SISL_REAL or ::SISL_COMPLEX.
 * 
 * @return the newly allocated ::sisl_vector_t.
 */

sisl_vector_t *sisl_vector_new(sisl_complex_t c)

{
  sisl_vector_t *v ;

  v = g_malloc(sizeof(sisl_vector_t)) ;
  v->c = c ;
  v->x = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;

  return v ;
}

/** 
 * Free a ::sisl_vector_t and underlying memory.
 * 
 * @param v a ::sisl_vector_t.
 * 
 * @return ::SISL_SUCCESS on success.
 */

gint sisl_vector_free(sisl_vector_t *v)

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  g_array_free(v->x, TRUE) ;
  g_free(v) ;

  return SISL_SUCCESS ;
}

/** 
 * Get the value of a vector element.
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index.
 * 
 * @return element \a i of \a v.
 */

gdouble sisl_vector_get(sisl_vector_t *v, gint i) 

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( (v->c == SISL_REAL) && (i >= 0) && (i < v->x->len) ) 
    return g_array_index(v->x,gdouble,i) ;
  g_error("%s: cannot get real element %d of length %d %s vector",
	  __FUNCTION__, i, sisl_vector_length((v)), 
	  sisl_complex_string((v)->c)) ;
  return 0.0 ;
}

/** 
 * Get the value of a complex vector element.
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index.
 * 
 * @return element \a i of \a v.
 */

gsl_complex sisl_vector_get_complex(sisl_vector_t *v, gint i) 

{
  gsl_complex x ;
  GSL_REAL(x) = GSL_IMAG(x) = 0.0 ;

  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( (v->c == SISL_COMPLEX) && (i >= 0) && (i < v->x->len) ) {
    GSL_REAL(x) = g_array_index(v->x,gdouble,2*i) ;
    GSL_IMAG(x) = g_array_index(v->x,gdouble,2*i+1) ;
    return x ;
  }
  g_error("%s: cannot get real element %d of length %d %s vector",
	  __FUNCTION__, i, sisl_vector_length((v)), 
	  sisl_complex_string((v)->c)) ;
  return x ;
}

/** 
 * Set an element of a vector
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index;
 * @param x value for element \a i of \a v.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_set(sisl_vector_t *v, gint i, gdouble x) 

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( sisl_is_real(v) && (i >= 0) && (i < v->x->len) ) {
    g_array_index(v->x,gdouble,i) = x ;
    return SISL_SUCCESS ;
  }
  g_error("%s: cannot set real element %d of length %d %s vector",
	  __FUNCTION__, i, sisl_vector_length((v)), 
	  sisl_complex_string((v)->c)) ;
  return SISL_ERROR ;
}

/** 
 * Set an element of a complex vector.
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index;
 * @param z value for element \a i of \a v.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_set_complex(sisl_vector_t *v, gint i, gsl_complex z)

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( !sisl_is_real(v) && (i >= 0) && (i < v->x->len) ) {
    g_array_index(v->x,gdouble,2*i) = GSL_REAL(z) ;
    g_array_index(v->x,gdouble,2*i+1) = GSL_IMAG(z) ;
    return SISL_SUCCESS ;
  }
  g_error("%s: cannot set complex element %d of length %d %s vector",
	  __FUNCTION__, i, sisl_vector_length((v)), 
	  sisl_complex_string((v)->c)) ;
  return SISL_ERROR ;
}

/** 
 * Get the norm of a real or complex vector.
 * 
 * @param v a ::sisl_vector_t;
 * @param n a ::sisl_norm_t.
 * 
 * @return the prescribed norm of \a v.
 */

gdouble sisl_vector_norm(sisl_vector_t *v, sisl_norm_t n) 

{
  gint j ;

  SISL_CHECK_ARGUMENT_NULL(v) ;

  switch (n) {
  default:
    g_error("%s: unrecognized norm type %d", __FUNCTION__, n) ;
    return 0.0 ;
    break ;
  case SISL_NORM_1:
    if ( sisl_is_real(v) )
      return dasum_((gint *)&(v->x->len), (gdouble *)(v->x->data), &c_1) ;
    else {
      j = v->x->len/2 ;
      return dzasum_(&j, (gdouble *)(v->x->data), &c_1) ;
    }
    break ;
  case SISL_NORM_2:
    if ( sisl_is_real(v) )
      return dnrm2_((gint *)(&(v->x->len)), (gdouble *)(v->x->data), &c_1) ;
    else {
      j = v->x->len/2 ;
      return dznrm2_(&j, (gdouble *)(v->x->data), &c_1) ;
    }
    break ;
  case SISL_NORM_INFINITY:
    if ( sisl_is_real(v) ) {
      j = idamax_((gint *)&(v->x->len), (gdouble *)(v->x->data), &c_1) ;
      return fabs(g_array_index(v->x,gdouble,j)) ;
    } else {
      j = v->x->len/2 ;
      j = izamax_(&j, (gdouble *)(v->x->data), &c_1) ;
      return sqrt(g_array_index(v->x,gdouble,2*j)*
		  g_array_index(v->x,gdouble,2*j)+
		  g_array_index(v->x,gdouble,2*j+1)*
		  g_array_index(v->x,gdouble,2*j+1)) ;
    }
    break ;
  }
  return 0.0 ;
}

/** 
 * Find the dot (scalar) product of two vectors of the same length.
 * 
 * @param v a ::sisl_vector_t;
 * @param w another ::sisl_vector_t.
 * 
 * @return \f$\sum_{i}v_{i}w_{i}\f$.
 */

gdouble sisl_vector_dot(sisl_vector_t *v, sisl_vector_t *w)
  
{
  gdouble d ;

  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_NULL(w) ;

  if ( !sisl_is_real(v) || !sisl_is_real(w) )
    g_error("%s: vectors v and w must both be real", __FUNCTION__) ;

  if ( sisl_vector_length(v) != sisl_vector_length(w) )
    g_error("%s: vectors v and w are not the same length (%d and %d)", 
	    __FUNCTION__, sisl_vector_length(v), sisl_vector_length(w)) ;

  d = ddot_((gint *)&(v->x->len),
	    (gdouble *)(v->x->data), &c_1,
	    (gdouble *)(w->x->data), &c_1) ;
  if ( isnan(d) ) g_error("%s: result is NaN", __FUNCTION__) ;

  return d ;
}

/** 
 * Find the dot (scalar) product of two complex vectors.
 * 
 * @param v a ::sisl_vector_t;
 * @param w another ::sisl_vector_t.
 * 
 * @return \f$\sum_{i}v_{i}w_{i}\f$.
 */

gsl_complex sisl_vector_dot_complex(sisl_vector_t *v, sisl_vector_t *w)
		       
{
  gsl_complex d ;
  gint len ;

  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_NULL(w) ;

  if ( sisl_is_real(v) || sisl_is_real(w) )
    g_error("%s: vectors v and w must both be complex", __FUNCTION__) ;

  if ( sisl_vector_length(v) != sisl_vector_length(w) )
    g_error("%s: vectors v and w are not the same length (%d and %d)", 
	    __FUNCTION__, sisl_vector_length(v), sisl_vector_length(w)) ;

  GSL_SET_COMPLEX(&d,0,0) ;
  len = sisl_vector_length(v) ;
  d = zdotu_(&len, (gdouble *)(v->x->data), &c_1, 
	     (gdouble *)(w->x->data), &c_1) ;

  if ( isnan(GSL_REAL(d)) || isnan(GSL_IMAG(d)) ) 
       g_error("%s: result is NaN", __FUNCTION__) ;

  return d ;
}

/** 
 * Copy one vector into another, resizing as necessary, \f$v=w\f$.
 * 
 * @param v destination ::sisl_vector_t;
 * @param w source ::sisl_vector_t.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_copy(sisl_vector_t *v, sisl_vector_t *w) 
  /* v := w */

{
  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_NULL(w) ;

  g_array_set_size(v->x, w->x->len) ;
  dcopy_((gint *)&(v->x->len), (gdouble *)(w->x->data), &c_1,
	 (gdouble *)(v->x->data), &c_1) ;
  v->c = w->c ;

  return SISL_SUCCESS ;  
}

/** 
 * Add one vector to another, in place, \a v += \a w.
 * 
 * @param v a ::sisl_vector_t;
 * @param w another ::sisl_vector_t.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_add(sisl_vector_t *v, sisl_vector_t *w) 

/* v := v + w */

{
  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_NULL(w) ;

  if ( sisl_is_real(v) && !sisl_is_real(w) )
    g_error("%s: if w is complex, v must also be complex", __FUNCTION__) ;

  if ( sisl_vector_length(v) != sisl_vector_length(w) )
    g_error("%s: vectors v and w are not the same length (%d and %d)", 
	    __FUNCTION__, sisl_vector_length(v), sisl_vector_length(w)) ;

  if ( !sisl_is_real(v) && sisl_is_real(w) ) {
    daxpy_((gint *)&(v->x->len), &d_1, 
	   (gdouble *)(w->x->data), &c_1,
	   (gdouble *)(v->x->data), &c_2) ;

    return SISL_SUCCESS ;
  }

  daxpy_((gint *)&(v->x->len), &d_1, 
	 (gdouble *)(w->x->data), &c_1,
	 (gdouble *)(v->x->data), &c_1) ;

  return SISL_SUCCESS ;
}

/** 
 * Subtract one vector from another, \a v -= \a w.
 * 
 * @param v a ::sisl_vector_t;
 * @param w another ::sisl_vector_t.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_sub(sisl_vector_t *v, sisl_vector_t *w) 

/* v := v - w */

{
  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_NULL(w) ;

  if ( sisl_is_real(v) && !sisl_is_real(w) )
    g_error("%s: if w is complex, v must also be complex", __FUNCTION__) ;

  if ( sisl_vector_length(v) != sisl_vector_length(w) )
    g_error("%s: vectors v and w are not the same length (%d and %d)", 
	    __FUNCTION__, sisl_vector_length(v), sisl_vector_length(w)) ;

  if ( !sisl_is_real(v) && sisl_is_real(w) ) {
    daxpy_((gint *)&(v->x->len), &d_m1, 
	   (gdouble *)(w->x->data), &c_1,
	   (gdouble *)(v->x->data), &c_2) ;
    return SISL_SUCCESS ;
  }

  daxpy_((gint *)&(v->x->len), &d_m1, 
	 (gdouble *)(w->x->data), &c_1,
	 (gdouble *)(v->x->data), &c_1) ;

  return SISL_SUCCESS ;
}

/** 
 * Add a number to an element of a vector, \f$v_{i} += x\f$. 
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index;
 * @param x number to add to element \a i of \a v.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_element_add(sisl_vector_t *v, gint i, gdouble x)

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( (v->c != SISL_REAL) || (i < 0) || (i >= v->x->len) )
    g_error("%s: cannot add to real element %d of length %d %s vector",
	    __FUNCTION__, i, sisl_vector_length((v)), 
	    sisl_complex_string((v)->c)) ;
  g_array_index(v->x,gdouble,i) += x ;
  
  return SISL_SUCCESS ;
}

/** 
 * Add a number to an entry of a complex vector, \f$v_{i} += z\f$. 
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index;
 * @param z number to add to element \a i of \a v.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_element_add_complex(sisl_vector_t *v, gint i, 
				     gsl_complex z)

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( (v->c == SISL_REAL) || (i < 0) || (i >= v->x->len/2) )
    g_error("%s: cannot add to complex element %d of length %d %s vector",
	    __FUNCTION__, i, sisl_vector_length((v)), 
	    sisl_complex_string((v)->c)) ;
  g_array_index(v->x,gdouble,2*i+0) += GSL_REAL(z) ;
  g_array_index(v->x,gdouble,2*i+1) += GSL_IMAG(z) ;
  
  return SISL_SUCCESS ;
}

/** 
 * Subtract a number from an element of a vector, \f$v_{i} -= x\f$. 
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index;
 * @param x number to subtract from element \a i of \a v.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_element_sub(sisl_vector_t *v, gint i, gdouble x)

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( (v->c != SISL_REAL) || (i < 0) || (i >= v->x->len) )
    g_error("%s: cannot subtract from real element %d of "
	    "length %d %s vector",
	    __FUNCTION__, i, sisl_vector_length((v)), 
	    sisl_complex_string((v)->c)) ;
  g_array_index(v->x,gdouble,i) -= x ;
  
  return SISL_SUCCESS ;
}

/** 
 * Subtract a number from an element of a complex vector, \f$v_{i} -= z\f$. 
 * 
 * @param v a ::sisl_vector_t;
 * @param i element index;
 * @param z number to subtract from element \a i of \a v.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_element_sub_complex(sisl_vector_t *v, gint i, 
				     gsl_complex z)

{
  SISL_CHECK_ARGUMENT_NULL(v) ;

  if ( (v->c == SISL_REAL) || (i < 0) || (i >= v->x->len/2) )
    g_error("%s: cannot subtract from complex element %d of "
	    "length %d %s vector",
	    __FUNCTION__, i, sisl_vector_length((v)), 
	    sisl_complex_string((v)->c)) ;
  g_array_index(v->x,gdouble,2*i+0) -= GSL_REAL(z) ;
  g_array_index(v->x,gdouble,2*i+1) -= GSL_IMAG(z) ;
  
  return SISL_SUCCESS ;
}

/** 
 * Clear a ::sisl_vector_t, setting its length to zero.
 * 
 * @param v ::sisl_vector_t to clear.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_clear(sisl_vector_t *v)

{
  gint n ;
  SISL_CHECK_ARGUMENT_NULL(v) ;

  n = v->x->len ;
  g_array_set_size(v->x,0) ; g_array_set_size(v->x,n) ;
  
  return SISL_SUCCESS ;
}

/** 
 * Scale a ::sisl_vector_t by a real number.
 * 
 * @param v a ::sisl_vector_t;
 * @param x number to multiply entries of \a v by.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_scale(sisl_vector_t *v, gdouble x)

{
  gint i ;

  SISL_CHECK_ARGUMENT_NULL(v) ;
  if ( !sisl_is_real(v) ) 
    g_error("%s: vector v must be real", __FUNCTION__) ;

  /*there does not seem to be a compact BLAS function for this*/
  for ( i = 0 ; i < v->x->len ; i ++ ) sisl_vector_data(v)[i] *= x ;

  return SISL_SUCCESS ;
}

/** 
 * Scale the entries of a complex ::sisl_vector_t.
 * 
 * @param v a ::sisl_vector_t to be scaled;
 * @param z number to multiply entries of \a v by.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_scale_complex(sisl_vector_t *v, gsl_complex z)

{
  gint i ;

  SISL_CHECK_ARGUMENT_NULL(v) ;
  if ( sisl_is_real(v) ) 
    g_error("%s: vector v must be complex", __FUNCTION__) ;

  /*there does not seem to be a compact BLAS function for this*/
  for ( i = 0 ; i < v->x->len/2 ; i ++ ) {
    *((gsl_complex *)(&(g_array_index(v->x,gdouble,2*i)))) =
      gsl_complex_mul(*((gsl_complex *)(&(g_array_index(v->x,gdouble,2*i)))),
		      z) ;
  }

  return SISL_SUCCESS ;
}

/** 
 * Set all elements of a ::sisl_vector_t to the same value.
 * 
 * @param v a ::sisl_vector_t;
 * @param x value to set \a v to.
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_set_all(sisl_vector_t *v, gdouble x) 

{
  gint i ;

  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_REAL(v) ;

  for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) 
    (sisl_vector_data(v))[i] = x ;

  return SISL_SUCCESS ;
}

/** 
 * Set all entries of a complex vector to some value.
 * 
 * @param v a ::sisl_vector_t;
 * @param z value to set entries of \a v to. 
 * 
 * @return ::SISL_SUCCESS on success. 
 */

gint sisl_vector_set_all_complex(sisl_vector_t *v, gsl_complex z) 

{
  gint i ;

  SISL_CHECK_ARGUMENT_NULL(v) ;
  SISL_CHECK_ARGUMENT_COMPLEX(v) ;

  for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
    sisl_vector_data(v)[2*i+0] = GSL_REAL(z) ;
    sisl_vector_data(v)[2*i+1] = GSL_IMAG(z) ;
  }

  return SISL_SUCCESS ;
}

gint sisl_vector_set_length(sisl_vector_t *v, gint n)

{
  SISL_CHECK_ARGUMENT_NULL(v) ;
  if ( n < 0 ) 
    g_error("%s: vector length (%d) may not be negative", __FUNCTION__, n) ;

  if ( sisl_is_real(v) ) g_array_set_size(v->x,(n)) ;
  else g_array_set_size(v->x,(2*n)) ;

  buffer_double(sisl_vector_data_length(v)) ;

  return SISL_SUCCESS ;
}

gint sisl_vector_multiprocessor_sum(sisl_vector_t *v)

{
  gdouble *buf ;

  SISL_CHECK_ARGUMENT_NULL(v) ;

  buf = buffer_double(0) ;

  /*this is not ideal but will work for MPI 1.1; MPI 2.0 has the
    MPI_IN_PLACE option which would be a bit better*/
  wmpi_sum_all_double(buf, sisl_vector_data(v), sisl_vector_data_length(v)) ;
  g_memmove(sisl_vector_data(v), buf, 
	    sisl_vector_data_length(v)*sizeof(gdouble)) ;

  return SISL_SUCCESS ;
}

/**
 * @}
 * 
 */
