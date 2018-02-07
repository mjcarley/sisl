#ifndef _MPI_WRAPPER_H_INCLUDED_
#define _MPI_WRAPPER_H_INCLUDED_

#include <glib.h>
#include <unistd.h>

#ifdef WRAPPER_USE_MPI

#include "mpi.h"

#endif /*WRAPPER_USE_MPI*/

#ifndef MPI_ANY
#define MPI_ANY -1
#endif /*MPI_ANY*/

void mpi_log_func(const gchar *log_domain,
		  GLogLevelFlags log_level,
		  const gchar *message,
		  gpointer user_data) ;
gint mpi_initialize(gint *argc, gchar ***argv) ;
gint mpi_shutdown(void) ;
guint mpi_rank() ;
guint mpi_n_processes() ;
gint mpi_pause(void) ;
gint mpi_split_range(guint i, guint j, guint *pmin, guint *pmax) ;

gint mpi_log_status_set(guint rank, gboolean status) ;
gboolean mpi_log_status_check(guint rank) ;

gint mpi_sum_int(gint *a, gint *b, gint n) ;
gint mpi_sum_uint(guint *a, guint *b, gint n) ;
gint mpi_sum_double(gdouble *a, gdouble *b, gint n) ;

gint mpi_sum_all_int(gint *a, gint *b, guint n) ;
gint mpi_sum_all_uint(guint *a, guint *b, guint n) ;
gint mpi_sum_all_double(gdouble *a, gdouble *b, guint n) ;

gint mpi_checkpoint(guint c) ;

gint mpi_recv_int(gint *i, guint n, gint tag, gint sender) ;
gint mpi_send_int(gint *i, guint n, gint tag, gint recv) ;
gint mpi_recv_uint(guint *i, guint n, gint tag, gint sender) ;
gint mpi_send_uint(guint *i, guint n, gint tag, gint recv) ;
gint mpi_recv_double(gdouble *i, gdouble n, gint tag, gint sender) ;
gint mpi_send_double(gdouble *i, gdouble n, gint tag, gint recv) ;

gint mpi_broadcast_int(gint *data, gint n, gint proc) ;
gint mpi_broadcast_uint(guint *data, gint n, gint proc) ;
gint mpi_broadcast_double(gdouble *data, gint n, gint proc) ;
#endif /*_MPI_WRAPPER_H_INCLUDED_*/

