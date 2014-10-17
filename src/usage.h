/*   File:                /home/edeveaud/Work/toppred/src/usage.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */


#ifndef __USAGE_H_
#define __USAGE_H_

#include "params.h"

#define PS "ps"
#define PNG "png"
#define PPM "ppm"
#define NONE "none"
#define X11 "x11"

/* output format */
#define NEW 1
#define OLD 2
#define HTML 3

/* usage display */
void usage(char *prog);

/* check plot format */
int check_output_format(char *prog, char *format);
void check_plot_format(char *prog, param_t *params);
void check_topo_format(char *prog, param_t *params);

#endif
