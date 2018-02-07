/* sisl-logging.c
 * 
 * Copyright (C) 2006 Michael Carley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/**
 * @defgroup logging Logging functions
 * 
 * @brief  Logging functions for use with GLIB logging facilities.
 * 
 * The logging functions provide control over the messages logged by
 * SISL, to ease debugging and to give access to messages from
 * specific processors during parallel execution.
 *
 * @{
 */

#include <stdio.h>

#include <glib.h>

#include "sisl.h"

#define SISL_LOGGING_DATA_WIDTH     4
#define SISL_LOGGING_DATA_FID       0
#define SISL_LOGGING_DATA_PREFIX    1
#define SISL_LOGGING_DATA_LEVEL     2
#define SISL_LOGGING_DATA_EXIT_FUNC 3

void sisl_logging_func(const gchar *log_domain,
		       GLogLevelFlags log_level,
		       const gchar *message,
		       gpointer data[]) ;

const gchar *sisl_logging_string(GLogLevelFlags level) ;

/** 
 * Return a string describing the log level of the message.
 * 
 * @param level log level.
 * 
 * @return string describing level.
 */

const gchar *sisl_logging_string(GLogLevelFlags level)

{
  const gchar *strings[] = {"RECURSION", 
			    "FATAL",
			    "ERROR",
			    "CRITICAL",
			    "WARNING",
			    "MESSAGE",
			    "INFO",
			    "DEBUG"} ;

  if ( G_LOG_LEVEL_ERROR & level) return strings[2] ; 
  if ( G_LOG_LEVEL_CRITICAL & level) return strings[3] ; 
  if ( G_LOG_LEVEL_WARNING & level) return strings[4] ; 
  if ( G_LOG_LEVEL_MESSAGE & level) return strings[5] ; 
  if ( G_LOG_LEVEL_INFO & level) return strings[6] ; 
  if ( G_LOG_LEVEL_DEBUG & level) return strings[7] ; 

  return NULL ;
}

/** 
 * Logging function for SISL. 
 * 
 * @param log_domain default log domain (see glib documentation);
 * @param log_level logging level (see glib documentation);
 * @param message logging message;
 * @param data array containing logging data, set by ::sisl_logging_init. 
 */

void sisl_logging_func(const gchar *log_domain,
		       GLogLevelFlags log_level,
		       const gchar *message,
		       gpointer data[])

{
  FILE *f = (FILE *)data[SISL_LOGGING_DATA_FID] ;
  gchar *p = (gchar *)data[SISL_LOGGING_DATA_PREFIX] ;
  GLogLevelFlags level = *(GLogLevelFlags *)data[SISL_LOGGING_DATA_LEVEL] ;
  gint (*exit_func)(void) = data[SISL_LOGGING_DATA_EXIT_FUNC] ;

  if ( log_level > level ) return ;

  fprintf(f, "%s%s-%s: %s\n", p, 
	  G_LOG_DOMAIN, sisl_logging_string(log_level),
	  message) ;

  if ( log_level <= G_LOG_LEVEL_ERROR ) {
    if ( exit_func != NULL ) exit_func() ;
  }

  return ;
}

/** 
 * Initialize SISL logging
 * 
 * @param f file stream for messages;
 * @param p string to prepend to messages;
 * @param log_level maximum logging level to handle (see gts_log);
 * @param exit_func function to call if exiting on an error.
 * 
 * @return 0 on success
 */

gint sisl_logging_init(FILE *f, gchar *p, 
		      GLogLevelFlags log_level,
		      gpointer exit_func)

{
  static gpointer data[SISL_LOGGING_DATA_WIDTH] ;
  static GLogLevelFlags level ;

  if ( f != NULL ) 
    data[SISL_LOGGING_DATA_FID] = f ;
  else
    data[SISL_LOGGING_DATA_FID] = stderr ;    
  if ( p != NULL ) 
    data[SISL_LOGGING_DATA_PREFIX] = g_strdup(p) ;
  else
    data[SISL_LOGGING_DATA_PREFIX] = g_strdup("") ;

  level = log_level ;
  data[SISL_LOGGING_DATA_LEVEL] = &level ;    
    
  g_log_set_handler (G_LOG_DOMAIN, 
		     G_LOG_FLAG_RECURSION |
		     G_LOG_FLAG_FATAL |   
		     G_LOG_LEVEL_ERROR |
		     G_LOG_LEVEL_CRITICAL |
		     G_LOG_LEVEL_WARNING |
		     G_LOG_LEVEL_MESSAGE |
		     G_LOG_LEVEL_INFO |
		     G_LOG_LEVEL_DEBUG,
		     (GLogFunc)sisl_logging_func, data);

  return 0 ;
}

/**
 * @}
 * 
 */
