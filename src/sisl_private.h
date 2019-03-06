/* sisl_private.h
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

#ifndef SISL_PRIVATE_H_INCLUDED
#define SISL_PRIVATE_H_INCLUDED

typedef struct _sisl_matrix_dense_t sisl_matrix_dense_t ;

struct _sisl_matrix_dense_t {
  gint size[4] ;
  GArray *x ;
} ;

#define SISL_MATRIX_DENSE_DATA(A) ((sisl_matrix_dense_t *)(A->m))

#define SISL_WORKSPACE_SIZE(w) (w->v->len)
#define SISL_WORKSPACE_VECTOR(w,i) (g_ptr_array_index(w->v,i))

#define SISL_CHECK_ARGUMENT_NULL(a)		\
  if ( a == NULL )				\
    g_error("%s: argument %s may not be NULL", __FUNCTION__, #a) 

#endif /*SISL_PRIVATE_H_INCLUDED*/
