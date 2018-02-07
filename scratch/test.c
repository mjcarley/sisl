#include <stdio.h>

#include <glib.h>

#include <sisl.h>

gint main()

{
  sisl_matrix_t *A ;
  sisl_vector_t *v, *w ;
  gdouble x, y ;
  gint i, j, k ;
  gsl_complex z ;
  gint nr, nc ;

  nr = 64 ; nc = 128 ;
  A = sisl_matrix_new(SISL_COMPLEX, SISL_MATRIX_SPARSE) ;
  sisl_matrix_column_number(A) = nc ; sisl_matrix_row_number(A) = nr ;
  sisl_matrix_local_row_start(A) = 0 ; sisl_matrix_local_row_end(A) = nr ;

  v = sisl_vector_new(SISL_COMPLEX) ;
  w = sisl_vector_new(SISL_COMPLEX) ;

  sisl_vector_set_length(v, nc) ;

  for ( k = 0 ; k < nr ; k ++ ) {
    fscanf(stdin, "%d %d %lg %lg", &i, &j, &x, &y) ;
    GSL_SET_COMPLEX(&z,x,y) ;
    sisl_matrix_set_complex(A, i, j, z) ;
/*     fscanf(stdin, "%d %d %lg", &i, &j, &x) ; */
/*     sisl_matrix_set(A, i, j, x) ;     */
  }

  for ( k = 0 ; k < nc ; k ++ ) {
    fscanf(stdin, "%d %d %lg %lg", &i, &j, &x, &y) ;
    GSL_SET_COMPLEX(&z,x,y) ;
    sisl_vector_set_complex(v, k, z) ;
/*     fscanf(stdin, "%d %d %lg", &i, &j, &x) ; */
/*     sisl_vector_set(v, k, x) ; */
  }

  sisl_matrix_vector_mul(A, v, w) ;

  for ( k = 0 ; k < sisl_vector_length(w) ; k ++ ) {
    z = sisl_vector_get_complex(w, k) ;
    fprintf(stdout, "%1.16e %1.16e\n", GSL_REAL(z), GSL_IMAG(z)) ;
/*     fprintf(stdout, "%1.16e\n", sisl_vector_get(w, k)) ; */
  }

  return 0 ;

/*   for ( i = 0 ; i < sisl_matrix_row_number(A) ; i ++ )  */
/*     for ( j = sisl_matrix_column_number(A)-1 ; j >= 0 ; j -- )  */
/*       sisl_matrix_set(A, i, j, i*10+j) ; */
/*   for ( j = sisl_matrix_column_number(A)-1 ; j >= 0 ; j -- )  */
/*       sisl_matrix_set(A, j, j, 3.0*j) ; */

/*   i = 51 ; j = 35 ; */
/*   x = sisl_matrix_get(A, i, j) ; */
/*   fprintf(stdout, "%1.16e ", x) ; */
  
/*   return 0 ; */

/*   i = 23 ; j = 41 ;  */
/*   sisl_matrix_set(A, i, j, 3.7) ;     */
/*   x = sisl_matrix_get(A, i, j) ; */
/*   fprintf(stderr, "%1.16e ", x) ; */

  for ( i = 0 ; i < sisl_matrix_row_number(A) ; i ++ ) {
    for ( j = 0 ; j < sisl_matrix_column_number(A) ; j ++ ) {
      x = sisl_matrix_get(A, i, j) ;
      fprintf(stdout, "%1.16e ", x) ;
/*       z = sisl_matrix_get_complex(A, i, j) ; */
/*       fprintf(stdout, "%1.16e %1.16e ", GSL_REAL(z), GSL_IMAG(z)) ; */
    }
    fprintf(stdout, "\n") ;
  }
/*   i = 0 ; j = 3 ; */
/*   x = sisl_matrix_get(A, i, j) ; */
/*   fprintf(stderr, "(%d, %d, %lg)\n", i, j, x) ; */
/*   sisl_matrix_set(A, i, j, 1.7) ; */
/*   x = sisl_matrix_get(A, i, j) ; */
/*   fprintf(stderr, "(%d, %d, %lg)\n", i, j, x) ; */

  return 0 ;
}
