/* iter_solver.h
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

#ifndef SISL_ITER_SOLVER_H_INCLUDED
#define SISL_ITER_SOLVER_H_INCLUDED

#include <glib.h>

#include "sisl.h"

gint sisl_solve_bicgstab(sisl_matrix_t *A, sisl_vector_t *x, 
			 sisl_vector_t *b,
			 gdouble tol, gint niter,
			 sisl_solver_workspace_t *w,
			 sisl_solver_performance_t *perf) ;
gint sisl_solve_bicgstab_c(sisl_matrix_t *A, sisl_vector_t *x, 
			   sisl_vector_t *b,
			   gdouble tol, gint niter,
			   sisl_solver_workspace_t *w,
			   sisl_solver_performance_t *perf) ;

gint sisl_solve_bicg(sisl_matrix_t *A, sisl_vector_t *x, 
		     sisl_vector_t *b,
		     gdouble tol, gint niter,
		     sisl_solver_workspace_t *w,
		     sisl_solver_performance_t *perf) ;
gint sisl_solve_cgs(sisl_matrix_t *A, sisl_vector_t *x, 
		    sisl_vector_t *b,
		    gdouble tol, gint niter,
		    sisl_solver_workspace_t *w,
		    sisl_solver_performance_t *perf) ;
gint sisl_solve_cgs_c(sisl_matrix_t *A, sisl_vector_t *x, 
		    sisl_vector_t *b,
		    gdouble tol, gint niter,
		    sisl_solver_workspace_t *w,
		      sisl_solver_performance_t *perf) ;

gint sisl_solve_cg(sisl_matrix_t *A, sisl_vector_t *x, sisl_vector_t *b,
		   gdouble tol, gint niter,
		   sisl_solver_workspace_t *w,
		   sisl_solver_performance_t *perf) ;

#endif /*SISL_ITER_SOLVER_H_INCLUDED*/
