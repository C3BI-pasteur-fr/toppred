/* -----------------------------------------------------------------
 file      : /home/schuerer/toppred/src/topoprint.c

 author    : Schuerer
 description :


-------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#include "topoprint.h"
#include "loop.h"


/* init topology_t object with topo_t informations */

/* transforme topologie as topo_t object to topology_t with
   informations of all elements  */
int tp_decode (topoprint_t *topo, elem_t *elems, seq_t *sq, param_t *para)
{
  int putatives = topo->putatives;
  seg_t *segments = elems->segs;
  int nsegments = elems->nsegs;

  scale_t *scales = &(para->scales);
  int clen = para->seg_len;

  char *seq = sq->seq;
  size_t len;

  int krstart, krstop;
  double pos , neg;
  int type_loop, triangle, last_stop, first;
  int integrate, k, i, nel, side, kingdom;

  elprint_t *el = topo->elps;
  loprint_t *lo;
  tmprint_t *tm;

  len = strlen(seq);

  triangle = para->q;
  type_loop = (topo->kr > 0) ? 1 : (topo->kr < 0) ? -1 : 0;
  last_stop = 0;
  k=0;
  kingdom = para->kingdom;
  for(i=0; i<nsegments; i++) {

    /* segment integrated ? */
    integrate = 1;
    if (segments[i].kind == PUTATIVE)
      { integrate = putatives % 2 ; putatives /= 2; }

    /* */
    if (integrate)
      {
	/* loop */
	lo = &(el[k].lo);
	lo->type = (type_loop == 1) ? "CYT_LOOP" : (type_loop == -1) ? "EXT_LOOP" : "LOOP";
	lo->start = last_stop;
	lo->stop = segments[i].start;
	krstart = (lo->start - triangle < 0) ? 0 :  lo->start - triangle;
	krstop =  (lo->stop + triangle > len) ? len : lo->stop + triangle;
	lo->len = lo->stop - lo->start;
	if (lo->len != 0) {
	  lo->kr = countkr (seq, krstart, krstop);
	  lo->cytext = distance (seq, lo->start, lo->stop, scales);
	}
	else { lo->kr =0; lo->cytext = 0.0; }

	/* segment */
	k++;
	tm = &(el[k].tm);
	tm->type = "TRANSMEM";
	tm->nr = i+1;
	tm->start = segments[i].start;
	tm->stop = segments[i].stop;
	tm->len = tm->stop - tm->start;
	tm->prob = segments[i].probTM;
	/* tm->type_seg = (segments[i].kind == CERTAIN) ? "certain" : "putative"; */
	tm->H = segments[i].H;

	/* init next loop and segment */
	k++;
	type_loop *= -1;
	last_stop = tm->stop;
      } /* end if */
  } /* end for */

  /* last loop */
  lo = &(el[k].lo);
  lo->type = (type_loop == 1) ? "CYT_LOOP" : (type_loop == -1) ? "EXT_LOOP" : "LOOP";
  lo->start = last_stop;
  lo->stop = len;
  krstart = (lo->start - triangle < 0) ? 0 :  lo->start - triangle;
  krstop =  (lo->stop + triangle > len) ? len : lo->stop + triangle;
  lo->len = lo->stop - lo->start;
  if (lo->len != 0) {
    lo->kr = countkr (seq, krstart, krstop);
    lo->cytext = distance (seq, lo->start, lo->stop, scales);
  }
  else { lo->kr =0; lo->cytext = 0.0; }

  /* topology summary */
  nel = k+1; side = 1;
  topo->prob = 1.0;
  topo->cytext = 0.0;
  for (i=0; i<nel; i++) {
    if (i%2) {/* tm segment (index impair) */
      topo->prob = (el[i].tm.prob < topo->prob) ? el[i].tm.prob : topo->prob;
      side *= -1;
    }
    else { /* loop (index pair) */
      if (el[i].lo.len > clen) topo->cytext += el[i].lo.cytext * (double) side;
    }
  }

  first = el[1].tm.nr - 1;
  topo->nterm = nterminus(seq, segments[first].start, segments[first].stop, para->q);

  topo->ncharge = (el[0].lo.start != -1 && el[0].lo.len <= clen) ?
    ncharge(seq, kingdom) : 0;
  topo->pos = pos = el[0].lo.kr;
  topo->neg = neg = countneg(seq, el[0].lo.start, el[0].lo.stop);
  if (neg+pos == 0) topo->negpos = 0.0;
  else topo->negpos = ((double) (neg - pos))/((double) (neg+pos));

  return nel;
}


void tp_new_print (topoprint_t *topo, int nel, param_t *para) {

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
  if (para->web == 1)
    {
    (void) fprintf (OUT, "%8s <A HREF=\"./%s-%d.png\">%3d</A>                  %3.2f       %6.2f   %6.2f  %8.2f %8.2f\n",
		  "TOPOLOGY", topo->image, topo->nr, topo->nr, topo->prob,
		    (double) topo->kr,
		    topo->cytext, (double) topo->nterm, topo->negpos);
  }
  else {
    (void) fprintf (OUT, "%8s %3d                  %3.2f       %6.2f   %6.2f  %8.2f %8.2f\n",
		  "TOPOLOGY", topo->nr, topo->prob,
		  (double) topo->kr,
		  topo->cytext, (double) topo->nterm, topo->negpos);
  }
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

  /* end topology */
  (void) fprintf(OUT, "//\n");

  return;
}


void tp_toppred_fprintf (topoprint_t *topo, int nel, param_t *para) {

  FILE *OUT = para->OUT;
  int clen = para->seg_len;
  elprint_t *el = topo->elps;
  int i;

  (void) fprintf(OUT, "Structure %d\n\n", topo->nr);

  /* element information of the structure */
  (void) fprintf(OUT, "Transmembrane segments included in this structure:\n");

  /* number of included segments */
  (void) fprintf(OUT, "     Segment  ");
  for (i=0; i<nel; i++) {
    if (i%2) { /* tmsegment (index impair) */
      (void) fprintf(OUT, "%6d", el[i].tm.nr);
    }
  }
  (void) fprintf(OUT, "\n");

  /* loop length */
  (void) fprintf(OUT, " Loop length");
  for(i=0; i<nel; i++) {
     if (!(i%2)) { /* loop (index pair) */
       (void) fprintf(OUT, "%6d", el[i].lo.len);
     }
  }
  (void) fprintf(OUT, "\n");

  /* k+r profile */
  (void) fprintf(OUT, " K+R profile");

  for(i=0; i<nel; i++) {
    if ((i%4) == 0) { /* loop (index pair) */
      if (el[i].lo.len <= clen) { (void) fprintf(OUT, "%6.2f", (!i) ?  (double) (el[i].lo.kr+topo->ncharge) : (double) el[i].lo.kr); }
      else { (void) fprintf(OUT, "%6s", "+"); }
    }
    if ((i%4) == 2) { (void) fprintf(OUT, "      "); }
  }
  (void) fprintf(OUT, "\n");

  (void) fprintf(OUT, "            ");
  for(i=0; i<nel; i++) {
    if ((i%4) == 0) { (void) fprintf(OUT, "      "); }
    else if ((i%4) == 2) {
      if (el[i].lo.len <= clen) { (void) fprintf(OUT, "%6.2f", (double) el[i].lo.kr); }
      else { (void) fprintf(OUT, "%6s", "+"); }
    }
  }
  (void) fprintf(OUT, "\n");

  /* cyt-ext profile */
  (void) fprintf(OUT, "CYT-EXT prof");

  for(i=0; i<nel; i++) {
    if ((i%4) == 0) { /* loop (index pair) */
      if (el[i].lo.len > clen) { (void) fprintf(OUT, "%6.2f", el[i].lo.cytext); }
      else { (void) fprintf(OUT, "%6s", "-"); }
    }
    if ((i%4) == 2) { (void) fprintf(OUT, "      "); }
  }
  (void) fprintf(OUT, "\n");

  (void) fprintf(OUT, "            ");
  for(i=0; i<nel; i++) {
    if ((i%4) == 0) { (void) fprintf(OUT, "      "); }
    else if ((i%4) == 2) {
      if (el[i].lo.len > clen) { (void) fprintf(OUT, "%6.2f", el[i].lo.cytext); }
      else { (void) fprintf(OUT, "%6s", "-"); }
    }
  }
  (void) fprintf(OUT, "\n");

  (void) fprintf(OUT, "For CYT-EXT profile neg. values indicate cytoplasmic preference.\n\n");

  /* topology summary */

  (void) fprintf(OUT, "\nK+R difference: %.2f\n", (double) (topo->kr));
  (void) fprintf(OUT, "Tm probability: %.2f\n", topo->prob);
  (void) fprintf(OUT, "-> Orientation: %s\n", orientation((double) (topo->kr)));
  (void) fprintf(OUT, "\nCharge-difference over N-terminal Tm (+-15 residues): %.2f\n", (double) topo->nterm);
  (void) fprintf (OUT, "%20s: %.4f\n", "(NEG-POS)/(NEG+POS)", topo->negpos);

  (void) fprintf (OUT, "%20s: %.4f\n", "NEG", (double) topo->neg);
  (void) fprintf (OUT, "%20s: %.4f\n", "POS", (double) topo->pos);

  (void) fprintf(OUT, "-> Orientation: %s\n", orientation((double)topo->nterm));
  (void) fprintf(OUT, "\nCYT-EXT difference: %6.2f\n", topo->cytext);
  (void) fprintf(OUT, "-> Orientation: %s\n", orientation(topo->cytext * (-1.0)));

  return;
}
