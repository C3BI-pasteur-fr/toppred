/*   File:                /home/edeveaud/Work/toppred/src/params.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */


#ifndef __PARAMS_H_
#define __PARAMS_H_

typedef struct scale_S {

  double cyt[26];
  double ext[26];
  double sd[26];

} scale_t;

typedef struct param_S {
  double Hdatas[26];
  double c_cut;	/* certain cut-off */
  double p_cut;	/* putative cut-off */
  int n;		/* core window */
  int q;		/* wedge window */
  int seg_len;	/* segment lentgh */
  int kingdom;       /* kingdom of species to which the sequence belongs */
  int n_topos;	/* numbers of topologies to calculate/print */
  int gplot;		/* to produce the ps Hphobe profile file or not */
  int hydro_file;	/* to produce the hydro file or not */
  char *out_dir;	/* where to store results files */
  FILE *OUT;		/* where to display the results */
  scale_t scales;       /* cyt-ext-sd scale file */

  char *plot_pause;
  char *topo_format;/* wich format to produce for the topo images */
  char *data_file;	/* hydrophobicity scale */
  char *ce_file;     /* cyt-ext-sd scale file */
  char *plot_format;	/* wich format to produce */
  char *plot_outfile;	/* prefix for the plot file */

  int output;        /* output format */
  int tmspacer;	/* critical tm distance */
  int web;		/* for web output format */
} param_t;

#define BUFFLEN 100
#define TRUE 1
#define FALSE 0


#endif
