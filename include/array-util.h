#ifndef _ARRAY_UTIL_H_INCLUDED_
#define _ARRAY_UTIL_H_INCLUDED_

#include <glib.h>

gint array_util_linear_search_uint(guint *i, guint ni, guint n) ;
gint array_util_lookup_uint(guint *i, guint ni, guint n) ;
gint array_util_unique_uint(guint *x, guint *n) ;
gint array_util_remove_uint(guint *x, guint *n, guint k) ;

gdouble array_util_norm2(gdouble *x, guint n) ;

gint array_util_unique_pointer(gpointer *x, guint *n) ;
#endif /* _ARRAY_UTIL_H_INCLUDED_*/
