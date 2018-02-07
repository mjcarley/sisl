#ifndef _VECTOR_PRIVATE_INCLUDED_
#define _VECTOR_PRIVATE_INCLUDED_

#include <stdio.h>
#include <glib.h>

#include "vector.h"

typedef struct {
  guint n,      /*number of elements*/
    nmax,       /*maximum number of elements*/
    *i,         /*element indices*/
    len ;       /*vector length*/
  gdouble *x ;  /*element values*/
} sisl_vector_sparse_t ;

typedef struct {
  guint n,       /*number of elements (equal to vector length)*/
    nmax ;       /*maximum number of elements*/
  gdouble *x ;   /*element values*/
} sisl_vector_dense_t ;

#define sisl_vector_dense_length(v) ((v)->n)
#define sisl_vector_sparse_length(v) ((v)->len)
#define sisl_vector_sparse_nelem(v) ((v)->n)

sisl_vector_sparse_t *_pvr_sparse_new(gint nmax, data_complex_t rc,
				    sisl_dist_t dist) ;
sisl_vector_dense_t *_pvr_dense_new(gint nmax, data_complex_t rc,
				  sisl_dist_t dist) ;

/*arithmetical operations*/
gint _pvr_inner_product_dense_dense(sisl_vector_dense_t *v,
				       sisl_vector_dense_t *w, gdouble *ip) ;
gint _pvr_inner_product_dense_sparse(sisl_vector_dense_t *v,
					   sisl_vector_sparse_t *w,
					   gdouble *ip) ;
gint _pvr_inner_product_sparse_sparse(sisl_vector_sparse_t *v,
					    sisl_vector_sparse_t *w,
					    gdouble *ip) ;
gint _pvr_scale_sparse(sisl_vector_sparse_t *v, gdouble x) ;
gint _pvr_scale_dense(sisl_vector_dense_t *v, gdouble x) ;
gint _pvr_sub_dense_sparse(sisl_vector_dense_t *v, sisl_vector_sparse_t *w) ;
gint _pvr_sub_sparse_dense(sisl_vector_sparse_t *v, sisl_vector_dense_t *w) ;
gint _pvr_sub_dense_dense(sisl_vector_dense_t *v, sisl_vector_dense_t *w) ;
gint _pvr_sub_sparse_sparse(sisl_vector_sparse_t *v, sisl_vector_sparse_t *w) ;
gint _pvr_add_dense_sparse(sisl_vector_dense_t *v, sisl_vector_sparse_t *w) ;
gint _pvr_add_dense_dense(sisl_vector_dense_t *v, sisl_vector_dense_t *w) ;
gint _pvr_add_sparse_sparse(sisl_vector_sparse_t *v, sisl_vector_sparse_t *w) ;
gint _pvr_norm_dense(sisl_vector_dense_t *v, sisl_norm_t nt, 
		    gdouble *n) ;
gint _pvr_norm_sparse(sisl_vector_sparse_t *v, sisl_norm_t nt,
		     gdouble *n) ;
gint _pvr_add_element_indexed_sparse(sisl_vector_sparse_t *v, guint i, 
				    guint idx, gdouble x) ;
gint _pvr_add_element_indexed_dense(sisl_vector_dense_t *v, guint i, 
				      guint idx, gdouble x) ;
gint _pvr_addto_element_indexed_sparse(sisl_vector_sparse_t *v, guint i, 
					  guint idx, gdouble x) ;
gint _pvr_addto_element_indexed_dense(sisl_vector_dense_t *v, guint i, 
					 guint idx, gdouble x) ;
gint _pvr_addto_block_indexed_weighted_dense(sisl_vector_dense_t *v, guint i,
						guint idx, gdouble *x,
						guint nb, gdouble w) ;
gint _pvr_addto_block_indexed_weighted_sparse(sisl_vector_sparse_t *v, guint i,
						 guint idx, gdouble *x,
						 guint nb, gdouble w) ;

gint _pvr_clear_dense(sisl_vector_dense_t *v) ;
gint _pvr_clear_sparse(sisl_vector_sparse_t *v) ;
gint _pvr_add_element_dense(sisl_vector_dense_t *v, guint i, gdouble x) ;
gint _pvr_add_element_sparse(sisl_vector_sparse_t *v, guint i, gdouble x) ;
gint _pvr_set_element_dense(sisl_vector_dense_t *v, guint i, gdouble x) ;
gint _pvr_set_element_sparse(sisl_vector_sparse_t *v, guint i, gdouble x) ;
gint _pvr_set_element_indexed_dense(sisl_vector_dense_t *v, 
				       guint i, guint idx, gdouble x) ;
gint _pvr_set_element_indexed_sparse(sisl_vector_sparse_t *v, 
					guint i, guint idx, gdouble x) ;

gint _pvr_set_length_dense(sisl_vector_dense_t *v, guint len) ;
gint _pvr_set_length_sparse(sisl_vector_sparse_t *v, guint len) ;
gint _pvr_write_dense(sisl_vector_dense_t *v,  data_complex_t rc,
		      FILE *f) ;
gint _pvr_write_sparse(sisl_vector_sparse_t *v, data_complex_t rc,
		       FILE *f) ;
gint _pvr_write_sparse_dense(sisl_vector_dense_t *v, guint row, FILE *f) ;
gint _pvr_write_sparse_sparse(sisl_vector_sparse_t *v, guint row, FILE *f) ;
gint _pvr_set_all_sparse(sisl_vector_sparse_t *v, gdouble x) ;
gint _pvr_set_all_dense(sisl_vector_dense_t *v, gdouble x) ;

gint _pvr_copy_dense_dense(sisl_vector_dense_t *v, sisl_vector_dense_t *w) ;
gint _pvr_copy_sparse_sparse(sisl_vector_sparse_t *v, sisl_vector_sparse_t *w) ;
gint _pvr_copy_dense_sparse(sisl_vector_dense_t *v, sisl_vector_sparse_t *w) ;

gdouble _pvr_get_element_dense(sisl_vector_dense_t *v, guint i) ;
gdouble _pvr_get_element_sparse(sisl_vector_sparse_t *v, guint i) ;
gdouble _pvr_get_element_indexed_sparse(sisl_vector_sparse_t *v, 
					   guint i, guint idx) ;
gdouble _pvr_get_element_indexed_dense(sisl_vector_dense_t *v, 
					  guint i, guint idx) ;
gint _pvr_compact_dense(sisl_vector_dense_t *v) ;
gint _pvr_compact_sparse(sisl_vector_sparse_t *v) ;
gint _pvr_addto_element_sparse(sisl_vector_sparse_t *v, guint i, gdouble x) ;
gint _pvr_addto_element_dense(sisl_vector_dense_t *v, guint i, gdouble x) ;
gint _pvr_wrap_array_dense(sisl_vector_dense_t *v, gdouble *x, guint n) ;
gint _pvr_unwrap_array_dense(sisl_vector_dense_t *v, gdouble **x, guint *n) ;
gint _pvr_multiprocessor_combine_dense(sisl_vector_dense_t *v, 
					  guint *indices) ;
gint _pvr_multiprocessor_combine_sparse(sisl_vector_dense_t *v, 
					   guint *indices) ;
gint _pvr_multiprocessor_sum_sparse(sisl_vector_sparse_t *v) ;
gint _pvr_multiprocessor_sum_dense(sisl_vector_dense_t *v) ;

gint _pvr_dump_write_dense(sisl_vector_dense_t *v, gchar rc,
			      gchar fmt, FILE *f) ;
gint _pvr_dump_write_sparse(sisl_vector_sparse_t *v, gchar rc,
			       gchar fmt, FILE *f) ;
gint _pvr_dump_read_sparse(sisl_vector_sparse_t *v, gchar fmt, gint n, FILE *f) ;
gint _pvr_dump_read_dense(sisl_vector_dense_t *v, gchar fmt, gint n, FILE *f) ;

gint _pvr_multiply_dense_sparse(sisl_vector_dense_t *v,
					   sisl_vector_sparse_t *w,
					   sisl_vector_t *p) ;
gint _pvr_multiply_dense_dense(sisl_vector_dense_t *v,
					  sisl_vector_dense_t *w,
					  sisl_vector_t *p) ;
gint _pvr_multiply_sparse_sparse(sisl_vector_sparse_t *v,
				       sisl_vector_sparse_t *w,
				       sisl_vector_t *p) ;

gint _pvr_add_weighted_dense_dense(sisl_vector_dense_t *v, 
				  sisl_vector_dense_t *w,
				  gdouble wt) ;
gint _pvr_add_weighted_dense_sparse(sisl_vector_dense_t *v, 
				   sisl_vector_sparse_t *w,
				   gdouble wt) ;
gint _pvr_add_weighted_sparse_sparse(sisl_vector_sparse_t *v, 
				    sisl_vector_sparse_t *w,
				    gdouble wt) ;
gint _pvr_free_dense(sisl_vector_dense_t *v) ;
gint _pvr_free_sparse(sisl_vector_sparse_t *v) ;

#endif /* _VECTOR_PRIVATE_INCLUDED_*/
