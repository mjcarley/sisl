#include <stdio.h>
#include <glib.h>

#include <sisl.h>
#include <gsl/gsl_complex.h>

gint main(gint argc, gchar **argv)

{
  sisl_vector_t *v, *w ;
  sisl_matrix_t *A ;
  gint i, j, lrow, lcol ;
  gdouble d, dc[2] ;
  gsl_complex x ;
  sisl_solver_workspace_t *ws ;
  sisl_solver_performance_t p ;

  wmpi_initialize(&argc, &argv) ;

/*   v = sisl_vector_new(SISL_REAL) ; */
/*   w = sisl_vector_new(SISL_REAL) ; */
/*   A = sisl_matrix_new(SISL_REAL, SISL_MATRIX_DENSE) ; */
  v = sisl_vector_new(SISL_COMPLEX) ;
  w = sisl_vector_new(SISL_COMPLEX) ;
  A = sisl_matrix_new(SISL_COMPLEX, SISL_MATRIX_DENSE) ;

  ws = sisl_solver_workspace_new() ;

  lrow = 14 ; lcol = 14 ;
  sisl_vector_set_length(v, lrow) ;
  sisl_vector_set_length(w, lcol) ;
  sisl_matrix_set_block_size(A, lrow, lcol) ;
  sisl_matrix_row_number(A) = lrow ;
  sisl_matrix_column_number(A) = lcol ;
  sisl_matrix_local_row_start(A) = 0 ;
  sisl_matrix_local_row_end(A) = lrow ;

  for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
    GSL_SET_COMPLEX(&x,1+(gdouble)i,(gdouble)i-0.5) ;
    sisl_vector_set_complex(v,i,x) ;
    GSL_SET_COMPLEX(&x,1+(gdouble)i,(gdouble)i+4.5) ;
    sisl_vector_set_complex(w,i,x) ;
/*     sisl_vector_set(v,i,(gdouble)i-0.5) ; */
  }
  
  for ( i = 0 ; i < lrow ; i ++ ) {
    for ( j = 0 ; j < lcol ; j ++ ) {
      GSL_SET_COMPLEX(&x,1+(gdouble)i*(gdouble)j,(gdouble)(j*j)) ;
      sisl_matrix_set_complex(A, i, j, x) ;
/*       sisl_matrix_set(A, i, j, 0.1+(gdouble)i*(gdouble)j) ; */
    }  
    GSL_SET_COMPLEX(&x,2+(gdouble)i*(gdouble)i,(gdouble)(i*i)) ;
    sisl_matrix_set_complex(A, i, i, x) ;
/*     sisl_matrix_set(A, i, i, 1+(gdouble)i*(gdouble)i) ; */
  }

  sisl_solve(SISL_SOLVER_BICGSTAB, A, w, v, 1e-6, 128, ws, &p) ;
/*   sisl_solve(SISL_SOLVER_CGS, A, w, v, 1e-6, 128, ws, &p) ; */

/*   sisl_vector_clear(w) ; */
/*   sisl_matrix_vector_mul(A,v,w) ; */

/*   for ( i = 0 ; i < sisl_matrix_row_number(A) ; i ++ ) { */
/*     for ( j = 0 ; j < sisl_matrix_row_number(A) ; j ++ ) */
/*       fprintf(stderr, "%lg+j%lg ", */
/* 	      GSL_REAL(sisl_matrix_get_complex(A,i,j)), */
/* 	      GSL_IMAG(sisl_matrix_get_complex(A,i,j))) ; */
/*     fprintf(stderr, "\n") ; */
/*   } */

  for ( i = 0 ; i < sisl_vector_length(w) ; i ++ )
    fprintf(stderr, "%lg %lg\n",
	    GSL_REAL(sisl_vector_get_complex(w,i)),
	    GSL_IMAG(sisl_vector_get_complex(w,i))) ;

/*     fprintf(stderr, "%lg\n", sisl_vector_get(w,i)) ; */

  return 0 ;
}
