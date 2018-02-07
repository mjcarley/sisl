#ifndef _MATRIX_PRIVATE_H_INCLUDED_
#define _MATRIX_PRIVATE_H_INCLUDED_

#include <glib.h>

#include "vector.h"

typedef struct {
  guint nrows, nrows_max, ncols, ncols_max, *i, n ;
  guint row_last, i_last ; /*caching data*/
  sisl_vector_t **rows ;
} sisl_mat_sparse_rows_t ;

typedef struct {
  guint 
  nrows, /*number of rows in matrix*/
    nrows_max, /*number of rows allocated on this processor*/
    ncols, /*number of columns in matrix*/
    ncols_max, /*number of columns allocated in each row*/    
    *row_indices ; /*indices of first and last rows held on 
		     each processor*/
  gdouble *x ; /*dense block matrix will have its data available here*/
  sisl_vector_t **rows ;
} sisl_mat_dense_rows_t ;

sisl_mat_dense_rows_t *_pm_new_dense_rows(guint nrows, guint ncol,
					    data_complex_t rc,
					    sisl_vector_density_t density) ;
sisl_mat_sparse_rows_t *_pm_new_sparse_rows(guint nrows, guint ncol,
					      data_complex_t rc,
					      sisl_vector_density_t density) ;
sisl_mat_dense_rows_t *_pm_new_dense_block(guint nrows, guint ncol,
					   sisl_vector_density_t density,
					   data_complex_t rc) ;
gint _pm_clear_dense_rows(sisl_mat_dense_rows_t *m) ;
gint _pm_clear_sparse_rows(sisl_mat_sparse_rows_t *m) ;
gint _pm_set_size_dense_rows(sisl_mat_dense_rows_t *m, guint rows,
				 guint cols) ;
gint _pm_set_size_sparse_rows(sisl_mat_sparse_rows_t *m, guint rows,
				  guint cols) ;
gint _pm_write_sparse_rows(sisl_mat_sparse_rows_t *m, FILE *f) ;
gint _pm_write_dense_rows(sisl_mat_dense_rows_t *m, FILE *f) ;
gint _pm_write_sparse_sparse_rows(sisl_mat_sparse_rows_t *m, FILE *f) ;
gint _pm_write_sparse_dense_rows(sisl_mat_dense_rows_t *m, FILE *f) ;

gint _pm_add_element_dense_rows(sisl_mat_dense_rows_t *m,
				    guint i, guint j, gdouble x) ;
gint _pm_add_element_sparse_rows(sisl_mat_sparse_rows_t *m,
				     guint i, guint j, gdouble x) ;
gint _pm_vector_multiply_dense_rows(sisl_mat_dense_rows_t *m, 
				    sisl_vector_t *v, sisl_vector_t *w) ;
gint _pm_vector_multiply_sparse_rows(sisl_mat_sparse_rows_t *m, 
					  sisl_vector_t *v, sisl_vector_t *w) ;
gint _pm_trans_vector_multiply_sparse_rows(sisl_mat_sparse_rows_t *m, 
					   sisl_vector_t *v, 
					   sisl_vector_t *w) ;
gint _pm_trans_vector_multiply_dense_rows(sisl_mat_dense_rows_t *m, 
					  sisl_vector_t *v, 
					  sisl_vector_t *w) ;
gint _pm_set_element_dense_rows(sisl_mat_dense_rows_t *m,
				    guint i, guint j, gdouble x) ;
gint _pm_set_element_sparse_rows(sisl_mat_sparse_rows_t *m,
				     guint i, guint j, gdouble x) ;
gint _pm_size_sparse_rows(sisl_mat_sparse_rows_t *m,
			      guint *rows, guint *cols) ;
gint _pm_size_dense_rows(sisl_mat_dense_rows_t *m,
			     guint *rows, guint *cols) ;
gdouble _pm_get_element_dense_rows(sisl_mat_dense_rows_t *m, 
				       guint i, guint j) ;
gdouble _pm_get_element_sparse_rows(sisl_mat_sparse_rows_t *m, 
					guint i, guint j) ;
gint _pm_compact_dense_rows(sisl_mat_dense_rows_t *m) ;
gint _pm_compact_sparse_rows(sisl_mat_sparse_rows_t *m) ;
gint _pm_addto_element_dense_rows(sisl_mat_dense_rows_t *m,
				      guint i, guint j, gdouble x) ;
gint _pm_addto_element_sparse_rows(sisl_mat_sparse_rows_t *m,
				       guint i, guint j, gdouble x) ;
gboolean _pm_has_row_sparse_rows(sisl_mat_sparse_rows_t *m, guint i) ;
gboolean _pm_has_row_dense_rows(sisl_mat_dense_rows_t *m, guint i) ;
sisl_vector_t *_pm_get_row_dense_rows(sisl_mat_dense_rows_t *m, guint i) ;
sisl_vector_t *_pm_get_row_sparse_rows(sisl_mat_sparse_rows_t *m, guint i) ;
gint _pm_set_all_sparse_rows(sisl_mat_sparse_rows_t *m,
				 gdouble x) ;
gint _pm_set_all_dense_rows(sisl_mat_dense_rows_t *m,
				gdouble x) ;
gint _pm_split_chunks_dense_rows(sisl_mat_dense_rows_t *m) ;
gint _pm_get_local_rows_dense_rows(sisl_mat_dense_rows_t *m,
				   gint *i0, gint *i1) ;
gint _pm_set_local_rows_dense_rows(sisl_mat_dense_rows_t *m,
				   gint i0, gint i1) ;
gint _pm_dump_write_dense_rows(sisl_mat_dense_rows_t *m, gchar fmt, 
				   gchar *file) ;
gint _pm_dump_read_dense_rows(sisl_mat_dense_rows_t *m, 
				  gchar fmt, FILE *f) ;
gint _pm_dump_read_sparse_rows(sisl_mat_sparse_rows_t *m, 
				   gchar fmt, FILE *f) ;
gint _pm_dump_read_new_dense_rows(sisl_mat_dense_rows_t *m, 
				      gchar fmt, FILE *f) ;
gint _pm_dump_read_new_sparse_rows(sisl_mat_sparse_rows_t *m, 
				       gchar fmt, FILE *f) ;
gint _pm_free_sparse_rows(sisl_mat_sparse_rows_t *m) ;
gint _pm_free_dense_rows(sisl_mat_dense_rows_t *m) ;
#endif /* _MATRIX_PRIVATE_H_INCLUDED_*/

