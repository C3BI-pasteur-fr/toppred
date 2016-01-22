/* -----------------------------------------------------------------
 file      : /home/schuerer/toppred/src/baseprint.c

 author    : Schuerer
 description :


-------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "output.h"

#include "error.h"

/* internal macros */

/* internal prototypes */

void print_sequence (FILE *OUT, char *seq) {

  int k;
  char *sq;

  sq = seq;
  k = 1;
  while(*sq){
    (void)fputc(*sq, OUT);
    sq++;
    k++;
    if(k>60) {
      (void)fputc('\n', OUT);
      k=1;
    }
  }
  (void) fputc('\n', OUT);
}

void print_tmsummary (FILE *OUT, int nseg, seg_t *KSseg) {

  char *type;
  int i;

  (void)fprintf(OUT, "Candidate membrane-spanning segments:\n\n");
  (void)fprintf(OUT, " Helix Begin - End   Score Certainty\n");
  for(i=0; i< nseg; i++) {
    type = NULL;
    (void)fprintf(OUT, "%6d %5i - %-5i %-6.3f",
                  i+1, KSseg[i].start+1, KSseg[i].stop, KSseg[i].H);
    if(KSseg[i].kind == 0) type = "Putative";
    if(KSseg[i].kind == 1) type = "Certain";
    (void)fprintf(OUT, "%s\n", type);
  }
}

void new_print_parameters (param_t *p) {

  FILE *OUT = p->OUT;
  char *scale;
  char *cytext;

  scale = p->data_file;
  cytext = p->ce_file;

  (void) fprintf (OUT, "Algorithm specific parameters: \n\n");

  (void) fprintf (OUT, "Full window size : %d\n", 2*p->q+p->n);
  (void) fprintf (OUT, "Core window size : %d\n", p->n);
  (void) fprintf (OUT, "Wedge window size: %d\n", p->q);
  (void) fprintf (OUT, "Using hydrophobicity file: %s\n\n", scale);

  (void) fprintf (OUT, "Cutoff for certain transmembrane segments: %.2f\n",
		  p->c_cut);
  (void) fprintf (OUT, "Cutoff for putative transmembrane segments: %.2f\n",
		  p->p_cut);
  (void) fprintf (OUT, "Critical distance between 2 transmembrane segments: %d\n\n",
		  p->tmspacer);

  (void) fprintf (OUT, "Critical loop length:  %d\n\n", p->seg_len);
  (void) fprintf (OUT, "Kingdom: %s\n\n",
		  (p->kingdom == PROKARYOTE) ? "procaryote" : "eucaryote");
  (void) fprintf (OUT, "Using cyt/ext file: %s\n\n", cytext);
  /*  (void) fprintf (OUT, "\nCharge-pair energy: 0\n");  mettre en variable */

}

void old_print_firstparameters (param_t *p) {

  FILE *OUT = p->OUT;
  char *scale, *cytext;

  (void) fprintf (OUT, "\n");

  cytext = p->ce_file;
  scale = p->data_file;
  (void) fprintf (OUT, "Using hydrophobicity file: %s\n", scale);
  (void) fprintf (OUT, "Using cyt/ext file: %s\n", cytext);
}

void old_print_tmsummary (FILE *OUT, int nseg, seg_t *KSseg, char *seq) {

  char *type;
  int i;

  (void)fprintf(OUT, "Candidate membrane-spanning segments:\n\n");
  (void)fprintf(OUT, " Helix Begin   End   Score Certainty\n");
  for(i=0; i< nseg; i++) {
    type = NULL;
    (void)fprintf(OUT, "%6d %5i   %-5i %-6.3f",
                  i+1, KSseg[i].start+1, KSseg[i].stop, KSseg[i].H);
    if(KSseg[i].kind == 0) type = "Putative";
    if(KSseg[i].kind == 1) type = "Certain";
    (void)fprintf(OUT, "%s\n", type);
  }
}

void old_print_secondparameters (param_t *p) {

  FILE *OUT = p->OUT;

  (void) fprintf (OUT, "\n(p)rokaryotic or (e)ukaryotic: %c\n\n",
		  (p->kingdom == PROKARYOTE) ? 'p' : 'e');

  (void) fprintf (OUT, "\nCharge-pair energy: 0\n\n"); /* mettre en variable */
  (void) fprintf (OUT, "Length of full window (odd number!): %d\n\n", 2*p->q+p->n);
  (void) fprintf (OUT, "Length of core window (odd number!): %d\n\n", p->n);
  (void) fprintf (OUT, "Number of residues to add to each end of helix: 1\n\n");
  (void) fprintf (OUT, "Critical length: %d\n\n", p->seg_len);
  (void) fprintf (OUT, "Upper cutoff for candidates: %.2f\n\n", p->c_cut);
  (void) fprintf (OUT, "Lower cutoff for candidates: %.2f", p->p_cut);

}
