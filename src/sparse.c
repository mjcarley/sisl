/* sparse.h
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

#include "sisl_private.h"

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

gint sisl_sparse_mul_c(sisl_matrix_sparse_t *m, gdouble *v,
		       gdouble *w)

{
  gint i, j, k ;

  for ( k = 0 ; k < m->ij->len ; k ++ ) {
    i = g_array_index(m->ij,gint,k)/(m->size[SISL_COL_NUMBER_OFFSET]) ;
    j = g_array_index(m->ij,gint,k) % (m->size[SISL_COL_NUMBER_OFFSET]) ;
    g_assert(i*m->size[SISL_COL_NUMBER_OFFSET] + j == 
	     g_array_index(m->ij,gint,k)) ;
    fprintf(stderr, "%d %d %d %lg %lg\n", 
	    g_array_index(m->ij,gint,k), i, j,
	    g_array_index(m->x,gdouble,2*k+0),
	    g_array_index(m->x,gdouble,2*k+1)) ;
    w[2*i+0] += 
      v[2*j+0]*g_array_index(m->x,gdouble,2*k+0) -
      v[2*j+1]*g_array_index(m->x,gdouble,2*k+1) ;      
    w[2*i+1] += 
      v[2*j+1]*g_array_index(m->x,gdouble,2*k+0) +
      v[2*j+0]*g_array_index(m->x,gdouble,2*k+1) ;      
  }

  return 0 ;
}
