#ifndef _MATRIX_H_INCLUDED_
#define _MATRIX_H_INCLUDED_

#include <stdio.h>

#include <glib.h>

#include "vector.h"

typedef enum {
  MATRIX_SPARSE_ROWS,
  MATRIX_DENSE_ROWS,
  MATRIX_DENSE_BLOCK
} sisl_mat_layout_t ;

typedef struct {
  sisl_mat_layout_t layout ;
  sisl_dist_t dist ;
  data_complex_t rc ;
  gpointer m ;
} sisl_matrix_t ;

#define SISL_HEADER_LENGTH   64
#define SISL_HEADER_STRING   "MTXFILE"
#define SISL_HEADER_VERSION  "01.00"
#define SISL_HEADER_SIZE_FMT "%08dX%08d"

#define sisl_matrix_is_complex(m) ((m->rc)==DATA_COMPLEX)

sisl_matrix_t *sisl_mat_new(guint nrow, guint ncol, 
			    sisl_mat_layout_t layout, 
			    sisl_vector_density_t density,
			    data_complex_t rc, sisl_dist_t dist) ;
gint sisl_mat_free(sisl_matrix_t *m) ;
gint sisl_mat_add_element(sisl_matrix_t *m, guint i, guint j, gdouble x) ;
gint sisl_mat_write(sisl_matrix_t *m, FILE *f) ;
gint sisl_mat_write_sparse(sisl_matrix_t *m, FILE *f) ;
gint sisl_mat_set_size(sisl_matrix_t *m, guint rows, guint cols) ;
gint sisl_mat_clear(sisl_matrix_t *m) ;
gint sisl_mat_vector_multiply(sisl_matrix_t *m, sisl_vector_t *v, 
			      sisl_vector_t *w) ;
gint sisl_mat_trans_vector_multiply(sisl_matrix_t *m, 
				    sisl_vector_t *v, 
				    sisl_vector_t *w) ;
gint sisl_mat_set_element(sisl_matrix_t *m, guint i, guint j, gdouble x) ;
gint sisl_mat_size(sisl_matrix_t *m, guint *rows, guint *cols) ;
gdouble sisl_mat_get_element(sisl_matrix_t *m, guint i, guint j) ;
gint sisl_mat_compact(sisl_matrix_t *m) ;
gint sisl_mat_set_distribution(sisl_matrix_t *m, sisl_dist_t dist) ;
gint sisl_mat_addto_element(sisl_matrix_t *m, guint i, guint j, gdouble x) ;
gboolean sisl_mat_has_row(sisl_matrix_t *m, guint i) ;
sisl_dist_t sisl_mat_distribution(sisl_matrix_t *m) ;
sisl_vector_t *sisl_mat_get_row(sisl_matrix_t *m, guint i) ;
gint sisl_mat_set_all(sisl_matrix_t *m, gdouble x) ;
gint sisl_mat_split_chunks(sisl_matrix_t *m) ;
gint sisl_mat_dump_write(sisl_matrix_t *m, gchar fmt, gchar *file) ;
gint sisl_mat_dump_read(sisl_matrix_t *m, FILE *f) ;
gint sisl_mat_new_dump_read(sisl_matrix_t **m, FILE *f) ;
gchar *sisl_mat_file_header_string(gint rows, gint cols, gchar fmt) ;
gint sisl_mat_get_local_rows(sisl_matrix_t *m, gint *i, gint *j) ;
gint sisl_mat_set_local_rows(sisl_matrix_t *m, gint i, gint j) ;
#endif /* _MATRIX_H_INCLUDED_*/
