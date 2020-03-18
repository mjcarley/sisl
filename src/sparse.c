/* sparse.h
 * 
 * Copyright (C) 2009, 2018 Michael Carley
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

#include "sisl-private.h"

gint sisl_sparse_mul(sisl_matrix_sparse_t *m, gdouble *v,
		     gdouble *w)

{
  gint i, j, k ;

  for ( k = 0 ; k < m->ij->len ; k ++ ) {
    i = g_array_index(m->ij,gint,k)/(m->size[SISL_COL_NUMBER_OFFSET]) ;
    j = g_array_index(m->ij,gint,k) % (m->size[SISL_COL_NUMBER_OFFSET]) ;
    g_assert(i*m->size[SISL_COL_NUMBER_OFFSET] + j == 
	     g_array_index(m->ij,gint,k)) ;
/*     fprintf(stderr, "%d %d %d\n", k, i, j) ; */
    w[i] += v[j]*g_array_index(m->x,gdouble,k) ;
  }

  return 0 ;
}

gint sisl_sparse_mul_complex(sisl_matrix_sparse_t *m, gdouble *v,
			     gdouble *w)

{
  gint i, j, k ;

  for ( k = 0 ; k < m->ij->len ; k ++ ) {
    i = g_array_index(m->ij,gint,k)/(m->size[SISL_COL_NUMBER_OFFSET]) ;
    j = g_array_index(m->ij,gint,k) % (m->size[SISL_COL_NUMBER_OFFSET]) ;
    g_assert(i*m->size[SISL_COL_NUMBER_OFFSET] + j == 
	     g_array_index(m->ij,gint,k)) ;
    /* fprintf(stderr, "%d %d %d %lg %lg\n",  */
    /* 	    g_array_index(m->ij,gint,k), i, j, */
    /* 	    g_array_index(m->x,gdouble,2*k+0), */
    /* 	    g_array_index(m->x,gdouble,2*k+1)) ; */
    w[2*i+0] += 
      v[2*j+0]*g_array_index(m->x,gdouble,2*k+0) -
      v[2*j+1]*g_array_index(m->x,gdouble,2*k+1) ;      
    w[2*i+1] += 
      v[2*j+1]*g_array_index(m->x,gdouble,2*k+0) +
      v[2*j+0]*g_array_index(m->x,gdouble,2*k+1) ;      
  }

  return 0 ;
}

gint sisl_dense_sparse_mul_w(sisl_matrix_dense_t *A,
			     sisl_matrix_sparse_t *B,
			     gdouble a, gdouble c,
			     sisl_matrix_dense_t *C)

{
  gint i, j, k, n, nra, nca, /* nrb, */ ncb, /* nrc, */ ncc ;
  gdouble *xa, *xc ;

  nra = A->size[SISL_ROW_NUMBER_OFFSET] ;
  nca = A->size[SISL_COL_NUMBER_OFFSET] ;
  /* nrb = B->size[SISL_ROW_NUMBER_OFFSET] ; */
  ncb = B->size[SISL_COL_NUMBER_OFFSET] ;
  /* nrc = C->size[SISL_ROW_NUMBER_OFFSET] ; */
  ncc = C->size[SISL_COL_NUMBER_OFFSET] ;
  
  xc = (gdouble *)(C->x->data) ;
  for ( n = 0 ; n < C->x->len ; n ++ ) xc[n] *= c ;
  
  xa = (gdouble *)(A->x->data) ;
  for ( n = 0 ; n < B->ij->len ; n ++ ) {
    k = g_array_index(B->ij,gint,n)/ncb ;
    j = g_array_index(B->ij,gint,n) % ncb ;
    for ( i = 0 ; i < nra ; i ++ ) {
      xc[i*ncc+j] += a*xa[i*nca+k]*g_array_index(B->x,gdouble,n) ;
    }
  }
  
  return 0 ;
}

gint sisl_dense_sparse_mul_w_complex(sisl_matrix_dense_t *A,
				     sisl_matrix_sparse_t *B,
				     gsl_complex a, gsl_complex c,
				     sisl_matrix_dense_t *C)

{
  gint i, j, k, n, nra, nca, /* nrb, */ ncb, /* nrc, */ ncc ;
  gsl_complex *xa, *xc, z ;
  
  nra = A->size[SISL_ROW_NUMBER_OFFSET] ;
  nca = A->size[SISL_COL_NUMBER_OFFSET] ;
  /* nrb = B->size[SISL_ROW_NUMBER_OFFSET] ; */
  ncb = B->size[SISL_COL_NUMBER_OFFSET] ;
  /* nrc = C->size[SISL_ROW_NUMBER_OFFSET] ; */
  ncc = C->size[SISL_COL_NUMBER_OFFSET] ;
  
  xc = (gsl_complex *)(C->x->data) ;
  for ( n = 0 ; n < C->x->len/2 ; n ++ ) {
    xc[n] = gsl_complex_mul(xc[n], c) ;
  }
  
  xa = (gsl_complex *)(A->x->data) ;
  for ( n = 0 ; n < B->ij->len ; n ++ ) {
    k = g_array_index(B->ij,gint,n)/ncb ;
    j = g_array_index(B->ij,gint,n) % ncb ;
    for ( i = 0 ; i < nra ; i ++ ) {
      GSL_SET_COMPLEX(&z,
		      g_array_index(B->x,gdouble,2*n+0),
		      g_array_index(B->x,gdouble,2*n+1)) ;
      xc[i*ncc+j] =
	gsl_complex_add(xc[i*ncc+j],
			gsl_complex_mul(a,
					gsl_complex_mul(xa[i*nca+k],z))) ;
    }
  }
  
  return 0 ;
}
