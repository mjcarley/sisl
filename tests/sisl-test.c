#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <glib.h>

#include <sisl.h>

#include <wmpi.h>

#include <gsl/gsl_complex_math.h>

gchar *cases[] = {"matrix_vector_mul",
		  "inverse",
		  "matrix_matrix_mul",
		  "matrix_vector_mul_diag",
		  "matrix_matrix_mul_weighted",
		  "matrix_triple_mul",
		  "matrix_dense_sparse_mul",
		  ""} ;

gint parse_case(gchar *arg) ;

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

static void matrix_inverse_test(gint nr, sisl_complex_t rc)

{
  sisl_matrix_t *A, *Ai, *C ;
  gint i, j ;
  gdouble x, offmax, onmax ;
  gsl_complex z ;
  
  A = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  Ai = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  C = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;

  sisl_matrix_random(A, nr, nr) ;

  sisl_matrix_copy(Ai, A) ;

  sisl_matrix_invert(Ai) ;

  if ( sisl_is_real(A) ) {
    /* for ( i = 0 ; i < sisl_matrix_row_number(Ai) ; i ++ ) { */
    /*   for ( j = 0 ; j < sisl_matrix_column_number(Ai) ; j ++ ) { */
    /* 	x = sisl_matrix_get(Ai, i, j) ; */
    /* 	fprintf(stderr, "%lg ", x) ; */
    /*   } */
    /*   fprintf(stderr, "\n") ; */
    /* } */

    /* fprintf(stderr, "\n") ; */

    sisl_matrix_matrix_mul(Ai, A, C) ;
    onmax = offmax = 0.0 ;
  
    for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
      for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
	x = sisl_matrix_get(C, i, j) ;
	if ( i != j ) 
	  offmax = MAX(offmax, fabs(x)) ;
	else
	  onmax = MAX(onmax, fabs(x-1.0)) ;
      }
    }
    
    fprintf(stderr, "Maximum error off diagonal: %lg; on diagonal: %lg\n",
	    offmax, onmax) ;
    
    sisl_matrix_matrix_mul(A, Ai, C) ;
    
    onmax = offmax = 0.0 ;
    
    for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
      for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
	x = sisl_matrix_get(C, i, j) ;
	if ( i != j ) 
	  offmax = MAX(offmax, fabs(x)) ;
	else
	  onmax = MAX(onmax, fabs(x-1.0)) ;
      }
    }
    
    fprintf(stderr, "Maximum error off diagonal: %lg; on diagonal: %lg\n",
	    offmax, onmax) ;
    return ;
  }

  /* for ( i = 0 ; i < sisl_matrix_row_number(Ai) ; i ++ ) { */
  /*   for ( j = 0 ; j < sisl_matrix_column_number(Ai) ; j ++ ) { */
  /*     z = sisl_matrix_get_complex(Ai, i, j) ; */
  /*     fprintf(stderr, "%lf+j%lf ", GSL_REAL(z), GSL_IMAG(z)) ; */
  /*   } */
  /*   fprintf(stderr, "\n") ; */
  /* } */

  /* fprintf(stderr, "\n") ; */

  sisl_matrix_matrix_mul(Ai, A, C) ;
  onmax = offmax = 0.0 ;
  
  for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
    for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
      z = sisl_matrix_get_complex(C, i, j) ;
      /* fprintf(stderr, "%lf+j%lf ", GSL_REAL(z), GSL_IMAG(z)) ; */
      if ( i != j ) 
	offmax = MAX(offmax, gsl_complex_abs(z)) ;
      else
	onmax = MAX(onmax,
		    gsl_complex_abs(gsl_complex_sub_real(z,1.0))) ;
    }
    /* fprintf(stderr, "\n") ; */
  }
  
  /* fprintf(stderr, "\n") ; */
  fprintf(stderr, "Maximum error off diagonal: %lg; on diagonal: %lg\n",
	  offmax, onmax) ;

    return ;
}

static void matrix_matrix_mul_test(sisl_matrix_t *A, sisl_matrix_t *B,
				   sisl_matrix_t *D)

{
  sisl_matrix_t *C ;
  gint i, j ;
  gdouble d, norm ;
  gsl_complex z ;
  
  C = sisl_matrix_new(A->c, SISL_MATRIX_DENSE) ;

  sisl_matrix_matrix_mul(A, B, C) ;

  norm = 0.0 ;
  
  if ( sisl_is_real(C) ) {
    for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
      for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
	d = sisl_matrix_get(C, i, j) - sisl_matrix_get(D, i, j) ;
	norm = MAX(norm, fabs(d)) ;
	/* d = sisl_matrix_get(C, i, j) ; */
	/* fprintf(stderr, "%lg ", d) ; */
      }
      /* fprintf(stderr, "\n") ; */
    }
    fprintf(stderr, "Norm: %lg\n", norm) ;

    return ;
  }

  for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
    for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
      z = sisl_matrix_get_complex(C, i, j) ;
      z = gsl_complex_sub(z, sisl_matrix_get_complex(D, i, j)) ;
      norm = MAX(norm, gsl_complex_abs(z)) ;
      /* fprintf(stderr, "%lg + j%lg ", GSL_REAL(z), GSL_IMAG(z)) ; */
    }
    /* fprintf(stderr, "\n") ; */
  }

  fprintf(stderr, "Norm: %lg\n", norm) ;

  return ;
}

static void matrix_vector_mul_diagonal_test(gint nr, gint nc,
					    sisl_complex_t rc) 

{
  sisl_matrix_t *A, *B ;
  sisl_vector_t *x, *y, *z ;
  gsl_complex a, b ;
  gint i ;
  
  A = sisl_matrix_new(rc, SISL_MATRIX_DIAGONAL) ;
  B = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;

  sisl_matrix_random(A, nr, nc) ;
  sisl_matrix_copy(B, A) ;
  
  x = sisl_vector_new(rc) ;
  y = sisl_vector_new(rc) ;
  z = sisl_vector_new(rc) ;

  sisl_vector_random(x, nc) ;
  
  sisl_matrix_vector_mul(A, x, y) ;
  sisl_matrix_vector_mul(B, x, z) ;

  if ( rc == SISL_REAL ) {
    for ( i = 0 ; i < sisl_vector_length(y) ; i ++ ) {
      fprintf(stderr, "%lg %lg %lg\n",
	      sisl_vector_get(y, i),
	      sisl_vector_get(z, i),
	      sisl_vector_get(y, i) - sisl_vector_get(z, i)) ;
    }
  } else {
    for ( i = 0 ; i < sisl_vector_length(y) ; i ++ ) {
      a = sisl_vector_get_complex(y, i) ;
      b = sisl_vector_get_complex(z, i) ;
      fprintf(stderr, "%lg+j%lg %lg+j%lg %lg+j%lg\n",
	      GSL_REAL(a), GSL_IMAG(a), 
	      GSL_REAL(b), GSL_IMAG(b), 
	      GSL_REAL(a) - GSL_REAL(b), 
	      GSL_IMAG(a) - GSL_IMAG(b)) ;
    }
  }

  sisl_vector_sub(z, y) ;

  fprintf(stderr, "Norm: %lg\n", sisl_vector_norm(z, SISL_NORM_1)) ;

  return ;
}

static void matrix_matrix_mul_weighted_test(gint nr, gint nc,
					    sisl_complex_t rc) 

{
  sisl_matrix_t *A, *B, *C, *D, *E ;
  gdouble wa, wc, norm, d ;
  gsl_complex wac, wcc, z1, z2 ;
  gint i, j ;
  
  A = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  B = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  C = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  D = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  E = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;

  sisl_matrix_random(A, nr, nc) ;
  sisl_matrix_random(B, nc, nr) ;
  sisl_matrix_random(C, nr, nr) ;
  sisl_matrix_copy(D, C) ;
  
  wa = g_random_double() ; wc = g_random_double() ;
  GSL_SET_COMPLEX(&wac, g_random_double(), g_random_double()) ;
  GSL_SET_COMPLEX(&wcc, g_random_double(), g_random_double()) ;
  
  sisl_matrix_matrix_mul(A, B, E) ;
  
  if ( rc == SISL_REAL ) {
    sisl_matrix_matrix_mul_w(A, B, C, wa, wc) ;
    norm = 0.0 ;
    for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
      for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
	d = wc*sisl_matrix_get(D, i, j) + wa*sisl_matrix_get(E, i, j) ;
	d -= sisl_matrix_get(C, i, j) ;
	norm = MAX(norm, fabs(d)) ;
      }
    }

    fprintf(stderr, "Norm: %lg\n", norm) ;
    
    return ;
  }

  sisl_matrix_matrix_mul_w_complex(A, B, C, wac, wcc) ;
  norm = 0.0 ;
  for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
    for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
      z1 = sisl_matrix_get_complex(D, i, j) ;
      z1 = gsl_complex_mul(z1, wcc) ;
      z2 = sisl_matrix_get_complex(E, i, j) ;
      z1 = gsl_complex_add(z1, gsl_complex_mul(z2, wac)) ;
      
      z2 = sisl_matrix_get_complex(C, i, j) ;
      z1 = gsl_complex_sub(z1, z2) ;
      norm = MAX(norm, gsl_complex_abs(z1)) ;
    }
  }
  
  fprintf(stderr, "Norm: %lg\n", norm) ;
    

  return ;
}

static void matrix_triple_mul_general_test(gint nr, gint nc,
					   sisl_matrix_density_t dA,
					   sisl_matrix_density_t dB,
					   sisl_matrix_density_t dC,
					   sisl_complex_t rc) 

{
  sisl_matrix_t *A, *B, *C, *D, *E, *F, *G ;
  gdouble wa, wd, norm, d ;
  gsl_complex wac, wdc, z1, z2 ;
  gint i, j, nra, nca, nrb, ncb, nrc, ncc ;
  
  fprintf(stderr, "weighted %s (%s)x(%s)x(%s) multiplication\n",
	  sisl_complex_string(rc),
	  sisl_density_string(dA),
	  sisl_density_string(dB),
	  sisl_density_string(dC)) ;

  A = sisl_matrix_new(rc, dA) ;
  B = sisl_matrix_new(rc, dB) ;
  C = sisl_matrix_new(rc, dC) ;
  D = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  E = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  F = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  G = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;

  if ( dB != SISL_MATRIX_DIAGONAL ) {
    nra = nr  ; nca = g_random_int_range(nc, nc*3) ;
    nrb = nca ; ncb = g_random_int_range(nc, nc*3) ;
    nrc = ncb ; ncc = nc ;
  } else {
    nra = nr  ; nca = g_random_int_range(nc, nc*3) ;
    nrb = nca ; ncb = nrb ;
    nrc = ncb ; ncc = nc ; 
  }
  if ( dC == SISL_MATRIX_DIAGONAL ) ncc = nrc ;
  
  sisl_matrix_random(A, nra, nca) ;
  sisl_matrix_random(B, nrb, ncb) ;
  sisl_matrix_random(C, nrc, ncc) ;
  sisl_matrix_random(D, nra, ncc) ;

  fprintf(stderr,
	  "  A: %dx%d\n"
	  "  B: %dx%d\n"
	  "  C: %dx%d\n",
	  nra, nca, nrb, ncb, nrc, ncc) ;
  
  wa = g_random_double() ; wd = g_random_double() ;
  
  GSL_SET_COMPLEX(&wac, g_random_double(), g_random_double()) ;
  GSL_SET_COMPLEX(&wdc, g_random_double(), g_random_double()) ;

  /*E = ABC*/
  sisl_matrix_copy(E, B) ;
  sisl_matrix_copy(G, C) ;
  sisl_matrix_matrix_mul(E, G, F) ;
  sisl_matrix_matrix_mul(A, F, E) ;
  /*F = D*/
  sisl_matrix_copy(F, D) ;
  
  if ( rc == SISL_REAL ) {
    sisl_matrix_triple_mul_w(A, B, C, D, wa, wd) ;
    fprintf(stderr, "  multiplication complete\n") ;

    norm = 0.0 ;
    for ( i = 0 ; i < sisl_matrix_row_number(D) ; i ++ ) {
      for ( j = 0 ; j < sisl_matrix_column_number(D) ; j ++ ) {
  	d = wd*sisl_matrix_get(F, i, j) + wa*sisl_matrix_get(E, i, j) ;
  	d -= sisl_matrix_get(D, i, j) ;
  	norm = MAX(norm, fabs(d)) ;
      }
    }

    fprintf(stderr, "Norm: %lg\n", norm) ;
    
    return ;
  }

  sisl_matrix_triple_mul_w_complex(A, B, C, D, wac, wdc) ;
  fprintf(stderr, "  multiplication complete\n") ;
  norm = 0.0 ;
  for ( i = 0 ; i < sisl_matrix_row_number(D) ; i ++ ) {
    for ( j = 0 ; j < sisl_matrix_column_number(D) ; j ++ ) {
      z1 = sisl_matrix_get_complex(F, i, j) ;
      z1 = gsl_complex_mul(z1, wdc) ;
      z2 = sisl_matrix_get_complex(E, i, j) ;
      z1 = gsl_complex_add(z1, gsl_complex_mul(z2, wac)) ;
      
      z2 = sisl_matrix_get_complex(D, i, j) ;
      z1 = gsl_complex_sub(z1, z2) ;
      norm = MAX(norm, gsl_complex_abs(z1)) ;
    }
  }
  
  fprintf(stderr, "Norm: %lg\n", norm) ;
    

  return ;
}


static void matrix_triple_mul_test(gint nr, gint nc, sisl_complex_t rc) 

{
  matrix_triple_mul_general_test(nr, nc,
				 SISL_MATRIX_DENSE,
				 SISL_MATRIX_DIAGONAL,
				 SISL_MATRIX_DENSE,
				 rc) ;

  matrix_triple_mul_general_test(nr, nc,
				 SISL_MATRIX_DENSE,
				 SISL_MATRIX_DIAGONAL,
				 SISL_MATRIX_DIAGONAL,
				 rc) ;

  matrix_triple_mul_general_test(nr, nc,
				 SISL_MATRIX_DENSE,
				 SISL_MATRIX_DIAGONAL,
				 SISL_MATRIX_SPARSE,
				 rc) ;

  return ;
}

static void matrix_dense_sparse_mul_test(gint nr, gint nc, sisl_complex_t rc)

{
  sisl_matrix_t *A, *B, *C, *D, *E, *F ;
  gdouble wa, wc, norm, d ;
  gsl_complex wac, wcc, z1, z2 ;
  gint i, j ;

  fprintf(stderr, "weighted (dense)x(sparse) %s multiplication\n",
	  sisl_complex_string(rc)) ;

  A = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  B = sisl_matrix_new(rc, SISL_MATRIX_SPARSE) ;
  C = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;

  D = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  E = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
  F = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;

  sisl_matrix_random(A, nr, nr) ;
  sisl_matrix_random(B, nr, nr) ;
  sisl_matrix_random(C, nr, nr) ;

  wa = g_random_double() ; wc = g_random_double() ;
  
  GSL_SET_COMPLEX(&wac, g_random_double(), g_random_double()) ;
  GSL_SET_COMPLEX(&wcc, g_random_double(), g_random_double()) ;

  /*E = B*/
  sisl_matrix_copy(D, B) ;
  sisl_matrix_copy(F, C) ;
  
  if ( rc == SISL_REAL ) {
    sisl_matrix_matrix_mul_w(A, B, C, wa, wc) ;

    sisl_matrix_matrix_mul(A, D, E) ;

    norm = 0.0 ;
    for ( i = 0 ; i < sisl_matrix_row_number(D) ; i ++ ) {
      for ( j = 0 ; j < sisl_matrix_column_number(D) ; j ++ ) {
  	d = wc*sisl_matrix_get(F, i, j) + wa*sisl_matrix_get(E, i, j) ;
  	d -= sisl_matrix_get(C, i, j) ;
  	norm = MAX(norm, fabs(d)) ;
      }
    }

    fprintf(stderr, "Norm: %lg\n", norm) ;
    
    return ;
  }

  sisl_matrix_matrix_mul_w_complex(A, B, C, wac, wcc) ;

  sisl_matrix_matrix_mul(A, D, E) ;

  norm = 0.0 ;
  for ( i = 0 ; i < sisl_matrix_row_number(C) ; i ++ ) {
    for ( j = 0 ; j < sisl_matrix_column_number(C) ; j ++ ) {
      z1 = sisl_matrix_get_complex(F, i, j) ;
      z1 = gsl_complex_mul(z1, wcc) ;
      z2 = sisl_matrix_get_complex(E, i, j) ;
      z1 = gsl_complex_add(z1, gsl_complex_mul(z2, wac)) ;
      
      z2 = sisl_matrix_get_complex(C, i, j) ;
      z1 = gsl_complex_sub(z1, z2) ;
      norm = MAX(norm, gsl_complex_abs(z1)) ;
    }
  }
  
  fprintf(stderr, "Norm: %lg\n", norm) ;
  
  return ;
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
  sisl_matrix_distribution(*A) = SISL_SINGLE ;
    /* DISTRIBUTED ; */

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

gint parse_case(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(cases[i]) != 0 ) {
    if ( !strcmp(cases[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}

gint main(gint argc, gchar **argv)

{
  sisl_matrix_t *A, *B, *D ;
  sisl_vector_t *v, *w, *r ;
  sisl_complex_t rc ;
  gdouble res ;
  FILE *f ;
  gint test, nr, nc, i ;
  gchar pstr[16] ;
  gchar ch, *ipfile ;
  
  wmpi_initialize(&argc, &argv) ;
  sprintf(pstr, "P%04d: ", wmpi_rank()) ;
  wmpi_log_status_set(wmpi_rank(), TRUE) ;
  /* wmpi_logging_init(stderr, pstr, G_LOG_LEVEL_DEBUG, wmpi_shutdown) ; */
  sisl_logging_init(stderr, pstr, G_LOG_LEVEL_DEBUG, wmpi_shutdown) ;

  g_debug("%s: %d processes, rank %d", __FUNCTION__, 
	  wmpi_process_number(), wmpi_rank()) ;
  test = 0 ;

  nr = 8 ; nc = 8 ;

  rc = SISL_REAL ;
  while ( (ch = getopt(argc, argv, "Cc:r:Tt:")) != EOF ) {
    switch (ch) {
    default:
    case 'h':
      fprintf(stderr, "Usage: %s (options) <input test file name>\n", argv[0]) ;
      fprintf(stderr, "Options:\n"
	      "  -C complex checks\n"
	      "  -c # number of columns in test matrices\n"
	      "  -r # number of rows in test matrices\n"
	      "  -T list tests and exit\n"
	      "  -t <test> test to run\n"
	      ) ;
      return 1 ;
      break ;
    case 'C': rc = SISL_COMPLEX ; break ;
    case 'c': nc = atoi(optarg) ; break ;
    case 'r': nr = atoi(optarg) ; break ;
    case 'T':
      fprintf(stderr, "Tests available:\n") ;
      i = 0 ;
      while ( strlen(cases[i]) != 0 ) {
	fprintf(stderr, "  %s\n", cases[i]) ;
	i ++ ;
      }
      return 0 ;
      break ;
    case 't':
      test = parse_case(optarg) ;
      break ;
    }
  }

  if ( test == 0 || test == 2 ) {
    /*these tests need external data*/
    ipfile = g_strdup(argv[argc-1]) ;
  
    g_debug("%s: opening input file %s", __FUNCTION__, ipfile) ;

    f = fopen(ipfile, "r") ;
    if ( f == NULL ) {
      fprintf(stderr, "%s: cannot open file %s\n", argv[0], ipfile) ;
      wmpi_shutdown() ;
      exit(1) ;
    }
  }
  
  if ( test == 0 ) {  
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

  if ( test == 1 ) {
    matrix_inverse_test(nr, rc) ;

    return 0 ;
  }
  
  if ( test == 2 ) {
    g_debug("%s: reading matrix", __FUNCTION__) ;

    matrix_read(&A, f) ;

    g_debug("%s: reading matrix", __FUNCTION__) ;

    matrix_read(&B, f) ;

    g_debug("%s: reading matrix", __FUNCTION__) ;

    matrix_read(&D, f) ;
    
    matrix_matrix_mul_test(A, B, D) ;

    return 0 ;
  }

  if ( test == 3 ) {
    matrix_vector_mul_diagonal_test(nr, nc, rc) ;
  }

  if ( test == 4 ) {
    matrix_matrix_mul_weighted_test(nr, nc, rc) ;
  }

  if ( test == 5 ) {
    matrix_triple_mul_test(nr, nc, rc) ;
  }

  if ( test == 6 ) {
    matrix_dense_sparse_mul_test(nr, nc, rc) ;
  }


  return 0 ;
}
