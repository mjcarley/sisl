#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include <sisl.h>

#include <wmpi.h>

static gdouble matrix_multiply_test(sisl_matrix_t *A, 
				    sisl_vector_t *b,
				    sisl_vector_t *w,
				    sisl_vector_t *r,
				    sisl_norm_t n)

{
  sisl_matrix_vector_mul(A, b, w) ;

  sisl_vector_sub(w, r) ;

  return sisl_vector_norm(w, n) ;
}

static void matrix_read(sisl_matrix_t **A, FILE *f)

{
  gint nc, nr, r, i, j ;
  guint rmin, rmax ;
  sisl_complex_t rc ;
  gsl_complex z ;
  gdouble x ;

  fscanf(f, "%d", &nr) ; fscanf(f, "%d", &nc) ;
  fscanf(f, "%d", &r) ;

  switch ( r ) { 
  case 1: rc = SISL_REAL ; break ;
  case 2: rc = SISL_COMPLEX ; break ;
  default: g_assert_not_reached() ; break ;
  }

  *A = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  sisl_matrix_distribution(*A) = SISL_DISTRIBUTED ;

  wmpi_split_range(0, nr, &rmin, &rmax) ; 

  sisl_matrix_set_block_size(*A, rmax-rmin, nc) ; 
  sisl_matrix_row_number(*A) = nr ;
  sisl_matrix_column_number(*A) = nc ;
  sisl_matrix_local_row_start(*A) = rmin ; 
  sisl_matrix_local_row_end(*A) = rmax ;

  switch ( rc ) {
  case SISL_REAL:
    for ( i = 0 ; i < nr ; i ++ )
      for ( j = 0 ; j < nc ; j ++ ) {
	fscanf(f, "%lg", &x) ;
	sisl_matrix_set(*A, i, j, x) ;
      }
    break ;
  case SISL_COMPLEX:
    for ( i = 0 ; i < nr ; i ++ )
      for ( j = 0 ; j < nc ; j ++ ) {
	fscanf(f, "%lg", &(GSL_REAL(z))) ;
	fscanf(f, "%lg", &(GSL_IMAG(z))) ;
	sisl_matrix_set_complex(*A, i, j, z) ;
      }
    break ;
  default: g_assert_not_reached() ; break ;
  }

  return ;
}

static void vector_read(sisl_vector_t **v, FILE *f)

{
  gint ne, r, i ;
  sisl_complex_t rc ;
  gsl_complex z ;
  gdouble x ;

  fscanf(f, "%d", &ne) ;
  fscanf(f, "%d", &r) ;

  switch ( r ) { 
  case 1: rc = SISL_REAL ; break ;
  case 2: rc = SISL_COMPLEX ; break ;
  default: g_assert_not_reached() ; break ;
  }

  *v = sisl_vector_new(rc) ;

  sisl_vector_set_length(*v, ne) ; 

  switch ( rc ) {
  case SISL_REAL:
    for ( i = 0 ; i < ne ; i ++ ) {
	fscanf(f, "%lg", &x) ;
	sisl_vector_set(*v, i, x) ;
    }
    break ;
  case SISL_COMPLEX:
    for ( i = 0 ; i < ne ; i ++ ) {
      fscanf(f, "%lg", &(GSL_REAL(z))) ;
      fscanf(f, "%lg", &(GSL_IMAG(z))) ;
      sisl_vector_set_complex(*v, i, z) ;
    }
    break ;
  default: g_assert_not_reached() ; break ;
  }

  return ;
}

gint main(gint argc, gchar **argv)

{
  sisl_matrix_t *A ;
  sisl_vector_t *v, *w, *r ;
  gdouble res ;
  FILE *f ;
  gchar pstr[16] ;

  wmpi_initialize(&argc, &argv) ;
  sprintf(pstr, "P%04d: ", wmpi_rank()) ;
  wmpi_log_status_set(wmpi_rank(), TRUE) ;
  wmpi_logging_init(stderr, pstr, G_LOG_LEVEL_DEBUG, wmpi_shutdown) ;
  sisl_logging_init(stderr, pstr, G_LOG_LEVEL_DEBUG, wmpi_shutdown) ;

  g_debug("%s: %d processes, rank %d", __FUNCTION__, 
	  wmpi_process_number(), wmpi_rank()) ;

  if ( argc < 2 ) {
    fprintf(stderr, "Usage: %s <input test file name>\n", argv[0]) ;
    wmpi_shutdown() ;
    exit(1) ;
  }

  g_debug("%s: opening input file %s", __FUNCTION__, argv[1]) ;

  f = fopen(argv[1], "r") ;
  if ( f == NULL ) {
    fprintf(stderr, "%s: cannot open file %s\n", argv[0], argv[1]) ;
    wmpi_shutdown() ;
    exit(1) ;
  }

  g_debug("%s: reading matrix", __FUNCTION__) ;

  matrix_read(&A, f) ;
  g_debug("%s: reading vector", __FUNCTION__) ;
  vector_read(&v, f) ;
  g_debug("%s: reading vector", __FUNCTION__) ;
  vector_read(&r, f) ;
  if ( sisl_is_real(v) )
    w = sisl_vector_new(SISL_REAL) ;
  else
    w = sisl_vector_new(SISL_COMPLEX) ;

  g_debug("%s: multiplication test", __FUNCTION__) ;
  res = matrix_multiply_test(A, v, w, r, SISL_NORM_1) ;

  if ( wmpi_rank() == 0 ) 
    fprintf(stderr, "residual: %lg\n", res) ;

  wmpi_shutdown() ;

  return 0 ;
}
