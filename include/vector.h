#ifndef _VECTORS_H_INCLUDED_
#define _VECTORS_H_INCLUDED_

#include <stdio.h>

#include <glib.h>

typedef enum {
  DATA_REAL,
  DATA_COMPLEX
} data_complex_t ;

typedef enum {
  VECTOR_DENSE,
  VECTOR_SPARSE,
  VECTOR_NULL    /*use with care*/
} sisl_vector_density_t ;

typedef enum {
  SISL_NORM_INFINITY,
  SISL_NORM_1,
  SISL_NORM_2
} sisl_norm_t ;

typedef enum {
  SISL_SINGLE,
  SISL_MULTI
} sisl_dist_t ;

typedef struct {
  sisl_vector_density_t density ;
  sisl_dist_t dist ;
  data_complex_t rc ;
  gpointer v ;
} sisl_vector_t ;

#define sisl_vector_density(v) ((v)->density)
#define sisl_vector_is_dense(v) (v->density==VECTOR_DENSE)
#define sisl_vector_distribution(v) (v->dist)
#define sisl_vector_is_complex(v) (((v)->rc==DATA_COMPLEX))

sisl_vector_t *sisl_vector_new(guint nmax, sisl_vector_density_t density, 
			       data_complex_t rc, sisl_dist_t dist) ;
gint sisl_vector_clear(sisl_vector_t *v) ;
gint sisl_vector_inner_product(sisl_vector_t *v, sisl_vector_t *w, 
			       gdouble *ip) ;
gint sisl_vector_vector_multiply(sisl_vector_t *v, sisl_vector_t *w, 
				 sisl_vector_t *p) ;
gint sisl_vector_add_element(sisl_vector_t *v, guint i, gdouble x) ;
gint sisl_vector_set_length(sisl_vector_t *v, guint len) ;
gint sisl_vector_write(sisl_vector_t *v, FILE *f) ;
gint sisl_vector_write_sparse(sisl_vector_t *v, guint i, FILE *f) ;
gint sisl_vector_set_all(sisl_vector_t *v, gdouble x) ;
gint sisl_vector_scale(sisl_vector_t *v, gdouble x) ;
gint sisl_vector_copy(sisl_vector_t *v, sisl_vector_t *w) ;
gint sisl_vector_sub(sisl_vector_t *v, sisl_vector_t *w) ;
gint sisl_vector_add(sisl_vector_t *v, sisl_vector_t *w) ;
gint sisl_vector_norm(sisl_vector_t *v, sisl_norm_t nt, 
		      gdouble *n) ;
gint sisl_vector_set_element(sisl_vector_t *v, guint i, gdouble x) ;
gint sisl_vector_set_element_indexed(sisl_vector_t *v, guint i, guint idx,
				   gdouble x) ;
gint sisl_vector_add_element_indexed(sisl_vector_t *v, guint i,
				guint idx, gdouble x) ;
gint sisl_vector_addto_element_indexed(sisl_vector_t *v, guint i, 
				  guint idx, gdouble x) ;
gint sisl_vector_addto_block_indexed_weighted(sisl_vector_t *v, guint i,
					 guint idx, gdouble *x,
					 guint nb, gdouble w) ;
gdouble sisl_vector_get_element(sisl_vector_t *v, guint i) ;
guint sisl_vector_length(sisl_vector_t *v) ;
gint sisl_vector_compact(sisl_vector_t *v) ;
gint sisl_vector_addto_element(sisl_vector_t *v, guint i, gdouble x) ;
guint sisl_vector_longest(guint n) ;
gpointer sisl_vector_operation_buffer(guint n) ;
gint array_lookup_uint(guint *i, guint ni, guint n) ;
gint sisl_vector_wrap_array(sisl_vector_t *v, gdouble *x, guint n) ;
gint sisl_vector_unwrap_array(sisl_vector_t *v, gdouble **x, guint *n) ;
gint sisl_vector_set_distribution(sisl_vector_t *v, sisl_dist_t d) ;
gdouble sisl_vector_get_element_indexed(sisl_vector_t *v, guint i, guint idx) ;
gint sisl_vector_multiprocessor_combine(sisl_vector_t *v, guint *indices) ;
gint sisl_vector_multiprocessor_sum(sisl_vector_t *v) ;
gint sisl_vector_dump_read(sisl_vector_t *v, gchar fmt, FILE *f) ;
gint sisl_vector_dump_read_new(sisl_vector_t **v, gchar fmt, FILE *f) ;
gint sisl_vector_dump_write(sisl_vector_t *v, gchar fmt, FILE *f) ;
gint sisl_vector_dump_skip(gchar fmt, FILE *f) ;
gint sisl_vector_dump_array(gdouble *x, gint n, sisl_vector_density_t density,
			    data_complex_t rc, gchar fmt, FILE *f) ;
gint sisl_vector_add_weighted(sisl_vector_t *v, sisl_vector_t *w, gdouble wt) ;
gint sisl_vector_free(sisl_vector_t *v) ;
#endif /* _VECTORS_H_INCLUDED_*/
