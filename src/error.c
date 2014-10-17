/* error.c - Error functions */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#include <errno.h>

#include "error.h"

#ifndef HAVE_STRERROR
char *strerror(int errnum) {
  extern char *sys_errlist[];
  extern int sys_nerr;

  if (errnum > 0 && errnum < sys_nerr) {
    return sys_errlist[errnum]; }

  return (char *)"Unknown error type"; }
#endif /* HAVE_STRERROR */


/* Abort on fatal error */
void error_fatal(const char *str, const char *err) {

  if (err == NULL) { err = strerror(errno); }
  (void)fprintf(stderr, "Fatal: %s: %s\n", str, err);

  exit(EXIT_FAILURE); }

/* Warn for non fatal error */
void error_warn(const char *str, const char *err) {

  if (err == NULL) { err = strerror(errno); }
  (void)fprintf(stderr, "Warning: %s: %s\n", str, err);

  return; }
