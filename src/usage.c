/*   File:                /home/edeveaud/Work/toppred/src/usage.c
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#include <stdio.h>

#include "usage.h"
#include "error.h"


void usage(char *prog) {
  FILE *PERR = stderr;
  (void)fprintf(PERR, "usage: %s [options] <file>\n", prog);

  (void)fprintf(PERR, "  -c <val>    ... Use <val> as certain cut-off.\n");

  (void)fprintf(PERR, "  -d <val>    ... Use <val> as critical distance between 2 transmembrane\n");
  (void)fprintf(PERR, "                  segments.\n");

  (void)fprintf(PERR, "  -e          ... For use with Eucaryotes.\n");

#ifdef HAVE_GNUPLOT
  (void)fprintf(PERR, "  -g <format> ... Display or produce Hydropphobic profile, in the specified\n");
  (void)fprintf(PERR, "                  <format> (ps, png, ppm, x11 and none).\n");
#endif /* HAVE_GNUPLOT */

  (void)fprintf(PERR, "  -h          ... Print this message and exit.\n");

  (void)fprintf(PERR, "  -H <file>   ... Use Hydrophobycitie values from <file>.\n");

  (void)fprintf(PERR, "  -n <val>    ... Use <val> as core window length.\n");

  (void)fprintf(PERR, "  -o <file>   ... Place the output into <file>.\n");

  (void)fprintf(PERR, "  -O <format> ... Print output in the specified\n");
  (void)fprintf(PERR, "                  <format> (old, new (default), html).\n");
  (void)fprintf(PERR, "  -p <val>    ... Use <val> as putative cut-off.\n");

  (void)fprintf(PERR, "  -q <val>    ... Use <val> as wedge window length.\n");

  (void)fprintf(PERR, "  -s <val>    ... Use <val> as critical loop length.\n");

#ifdef HAVE_LIBGD
  (void)fprintf(PERR, "  -t <format> ... Produce images of the topologies in the specified\n");
  (void)fprintf(PERR, "                  <format> (png and none).\n");
#endif

  (void)fprintf(PERR, "  -v          ... Print version number and exit.\n");
}

int check_output_format(char *prog, char *format) {

  if (strcmp (format, "html") == 0) { return HTML; }
  if (strcmp (format, "old") == 0) { return OLD; }
  if (strcmp (format, "new") == 0) { return NEW; }

  error_fatal(prog, "supported format are currently new, old and html");

  return 1;
}


void check_plot_format(char *prog, param_t *params) {
  char *format;
  format = params->plot_format;

  if(strcmp (format, PS) == 0) {
    params->plot_format = "set term postscript color solid\n";
    params->plot_outfile = "ps";
    return;
  }
  else if (strcmp (format, PNG) == 0) {
    params->plot_format= "set term png medium";
    params->plot_outfile = "png";
    return;
  }
  else if(strcmp (format, NONE) == 0) {
    params->gplot = FALSE;
    params->plot_format = NULL;
    return;
  }
  else 	if(strcmp (format, X11) == 0) {
    params->plot_format = "" ;
    params->plot_outfile = NULL;
    params->plot_pause = "pause -1";
    return;
  }
  else if(strcmp (format, PPM) == 0) {
    params->plot_format = "set term pbm medium color";
    params->plot_outfile = "ppm";
    return;
  }
  else{error_fatal(prog,
		   "supported format are currently ps, png, ppm, x11 and none");
  }
}

void check_topo_format(char *prog, param_t *params) {
  char *format;
  format = params->topo_format;

  if((strcmp (format, PNG) != 0) && (strcmp (format, NONE) !=0 ))  {
    error_fatal(prog,
		"supported format are currently png and none");
  }
}
