/* -----------------------------------------------------------------
 file      : /home/schuerer/toppred/src/mloutput.c

 author    : Schuerer
 description :


-------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "config.h"
#include "usage.h"

#include "output.h"
#include "mloutput.h"

/* internal macros */

/* internal prototypes */

FILE *init_html(char *name, char *dir) {
  char *filename;
  FILE *OUT;
  int len;

  len = strlen(dir) + 1 + strlen(name) + 5;

  if ((filename = malloc((size_t)len+1)) == NULL) {
    error_fatal("html output file", NULL); }

  (void)sprintf(filename, "%s/%s.html", dir, name);

  if ((OUT = fopen(filename, "w")) == NULL) {
    error_fatal(filename, NULL); }

  free(filename);

  return OUT; }


void html_header (FILE *OUT, char *name) {

  (void)fprintf(OUT, "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\"\n");
  (void)fprintf(OUT, " \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n");
  (void)fprintf(OUT, "<HTML>\n");

  (void)fprintf(OUT, "<HEAD>\n");
  (void)fprintf(OUT, "<META http-equiv=\"Content-Type\" content=\"text/html;");
  (void)fprintf(OUT, " charset=iso-8859-1\">\n");
  (void)fprintf(OUT, "<TITLE>Toppred prediction %s</TITLE>\n", name);
  (void)fprintf(OUT, "</HEAD>\n");

  (void)fprintf(OUT, "<BODY>\n");
  (void)fprintf(OUT, "<H1><CENTER>Toppred prediction %s</CENTER></H1>\n", name);
  (void)fprintf(OUT, "<PRE>\n");

  return; }


void html_parameters (FILE *OUT, param_t *p) {

  (void) fprintf (OUT, "<H4><CENTER>Algorithm specific parameters</CENTER></H4>\n\n");
  (void) fprintf (OUT, "<PRE>\n");
  new_print_parameters (p);
  (void) fprintf (OUT, "</PRE>\n");

}

void html_sequence (FILE *OUT, seq_t *seq) {

  (void)fprintf(OUT, "<H4><CENTER>Sequence</CENTER></H4>\n");
  start_phrase();
  (void) fprintf(OUT, "Sequence : %s  (%d res)\n", seq->id, seq->size);
  print_sequence(OUT, seq->seq);
  end_phrase();

}

void html_plot (FILE *OUT, char *seqname, param_t *p) {

  (void)fprintf(OUT, "<H4><CENTER>Hydrophobicity plot</CENTER></H4>\n");
  (void)fprintf(OUT, "<P><CENTER>\n");
#ifdef HAVE_GNUPLOT
  if (p->gplot && p->plot_outfile && !strcmp(p->plot_outfile, PNG)) {
    (void)fprintf(OUT, "<IMG SRC=\"./%s.png\">\n", seqname);
  }
#endif
  (void)fprintf(OUT, "<PRE>\n");
  (void)fprintf(OUT, "<P>\n");
  (void)fprintf(OUT, "<A HREF=\"./%s.hydro\">view hydrophobic values</A>\n",
		seqname);
  (void)fprintf(OUT, "</P>\n");
  (void)fprintf(OUT, "</CENTER></P>\n");

}

void html_topology (topoprint_t *topo, int nel, param_t *para) {

  FILE *OUT = para->OUT;
  int clen = para->seg_len;
  elprint_t *el = topo->elps;
  int i;

  /* header */
  (void) fprintf (OUT, "%8s %6s %6s %6s %4s %4s %8s %8s %8s %8s\n",
		  "HEADER  ", "START", "STOP", "LEN",
		  "PROB", "HP",
		  "DARGLYS", "DCYTEXT",  "DNCHARGE", "DNNEGPOS");


  /* topo summary */
  (void) fprintf (OUT, "%8s ", "TOPOLOGY");
#ifdef  HAVE_LIBGD
  if (!strcmp(para->topo_format, PNG)) {
    (void) fprintf (OUT, "<A HREF=\"./%s\">%3d</A>", topo->image, topo->nr);
  }
  else {
    (void) fprintf (OUT, "%3d", topo->nr);
  }
#else
  (void) fprintf (OUT, "%3d", topo->nr);
#endif

  (void) fprintf (OUT, "                  %3.2f       %6.2f   %6.2f  %8.2f %8.2f\n",
		  topo->prob, (double) topo->kr,
		  topo->cytext, (double) topo->nterm, topo->negpos);

  (void) fprintf (OUT, "%8s                                %7s  %7s  %8s\n",
		  "TOPOLOGY",
  		  new_orientation((double) (topo->kr+topo->ncharge)),
  		  new_orientation(topo->cytext * (-1.0)),
		  new_orientation((double) topo->nterm));

  /* topo elements */
  for (i=0; i<nel; i++) {
      if (i%2) { /* tmsegment (index impair) */
	(void) fprintf(OUT, "%8s %6d %6d %6d %4.2f %4.2f\n",
		       el[i].tm.type, el[i].tm.start+1, el[i].tm.stop,
		       el[i].tm.len, el[i].tm.prob, el[i].tm.H);
      }
      else { /* loop (index pair) */
	if (el[i].lo.len > clen) { /* cytext value is significant */
	  (void) fprintf(OUT, "%8s %6d %6d %6d           (%6.2f)  %6.2f\n",
			 el[i].lo.type, el[i].lo.start+1, el[i].lo.stop,
			 el[i].lo.len, (double) el[i].lo.kr, el[i].lo.cytext);
	}
	else if (el[i].lo.len != 0) { /* kr value is significant */
	  (void) fprintf(OUT, "%8s %6d %6d %6d            %6.2f  (%6.2f)\n",
			 el[i].lo.type, el[i].lo.start+1, el[i].lo.stop,
			 el[i].lo.len, (double) el[i].lo.kr, el[i].lo.cytext);
	}
      }
  }

  (void) fprintf(OUT, "\n");

  return;
}
