/*   File:                /home/edeveaud/Work/toppred/src/graph.c
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#include <graph.h>
#include <gd.h>
#include <gdfontl.h>
#include <gdfontg.h>


#include "error.h"
#include "topology.h"


static int TM_COLOR;
static int TM_BORDER;
static int MBR_COLOR;
static int PUT_COLOR;
static int MBR_BORDER;
static int FONT_COLOR;
static int LOOP_COLOR;

static float COEFF;

static FILE *file_maker(char *title, char *dir, int num) ;

static gdImagePtr init_img(char *title, int num) ;

static int *loop_tweaker (elprint_t *element, int n) ;

static void draw_segment(gdImagePtr img, int pos, int n, int kind) ;

static void draw_loop(gdImagePtr img, int x1,int x2, int length, int real,
		      int kr, int kind, int first_last) ;

static void print_legend(gdImagePtr img,char *legend) ;

void topo_graph_print (topoprint_t *topo, int n, param_t *params, seq_t *seq) {

  char *seq_name;
  int  i,j;
  int *loop_length;
  int ll, real, kr, loc,first_last;
  int num;
  int  x1, x2;
  char *legend, *p;
  char *legend_text = "Segments included:";
  char *dir_name;
  FILE *pngout;
  gdImagePtr img;
  elprint_t *element;

  /* could we draw */
  if((n-1)/2 >MAX_SEGMENTS){
    error_warn("graphic_topo", "too many segments to be represented, skipped");
    return;
  }

  /* tweaking the values to avoid drawing overlaping segments */
  element = topo->elps;
  if ((loop_length = loop_tweaker(element, n)) == NULL) {
    error_warn("graphic_topo", "could not display topologies");
    return;
  }

  /* generating image holder file */
  dir_name = params->out_dir;
  seq_name=seq->id;
  num =   topo->nr;
  pngout = file_maker(seq_name, dir_name, num);

  /* allocate image */
  img = init_img(seq_name, num);

  /* there is at least n/2 -1 segment included
   * each one is at max 2 char long + 1 for space
   * thus generous legend size is  */
  i = strlen(legend_text) + n*3;
  if((legend = (char *)malloc(sizeof(char) * (i+1))) == NULL){
    error_fatal("memory", NULL);
  }
  p =legend;
  (void)sprintf(p, "%s", legend_text);
  p+=strlen(legend_text);


  /* draw topologie segments and loops  */
  if ( topo->kr <= 0) loc = CYT;
  else loc = EXT;

  x1 = x2 = 0;
  j = 0;
  for (i=0; i<n; i++) {
    if (i%2) {
      draw_segment(img, x2, element[i].tm.nr, (int)element[i].tm.prob);
      (void)sprintf(p, "%3d", element[i].tm.nr);
      p+=3;
    }
    else {
      x2 = x1 + loop_length[j];
      ll = element[i].lo.len;
      real = loop_length[j];
      kr = element[i].lo.kr;
      if (i == 0) first_last = FIRST;
      else if (i == n-1) first_last = LAST;
      else first_last = INNER;
      draw_loop(img, x1, x2, ll, real, kr, loc,first_last ) ;
      x1 = x2;
      loc *= -1;
      j++;
    }
  }

  print_legend(img, legend);

  /* dropping image to file */
  gdImagePng(img, pngout);
  gdImageDestroy(img);
  
  if(fclose(pngout) == EOF) {
    error_fatal("closing file", NULL);
  }

  /* memory cleaning */
  free(legend);
  free(loop_length);

  return;
}


static void print_legend(gdImagePtr img,char *legend) {
  gdFontPtr font;
  font = gdFontGiant;
  gdImageString(img, font, 10, 50, (unsigned char *)legend, FONT_COLOR);

  return;
}

static void draw_loop(gdImagePtr img, int x1, int x2, int length, int real, int kr, int loc, int first_last) {

  /* loop position */
  int xc , yc;
  int l, start, stop;
  /* legend related */
  gdFontPtr font;
  char *legend1, *legend2;
  int xf1, xf2, yf1, yf2;

  font = gdFontLarge;

  /* positioning the loop params */
  l = (int)(COEFF * real);
  yc = Y_CENTER - ((SEG_H /2) * loc);
  if (first_last == INNER) {
    xc = (int)(((x1 + x2)/2.0) * COEFF) + MARGIN;
    start = 180 - ((1 - loc) * 90);
    stop  = start + 180;
  }
  else if (first_last == FIRST) {
    xc = MARGIN;
    start = 270 - ((1 - loc) * 135);
    stop  = start + 90;
    l *= 2;
    gdImageChar(img, font, xc , yc, 'N', FONT_COLOR) ;
  }
  else { /* first_last == LAST */
    xc = (int)(((x1 + x2)/2.0) * COEFF) + l/2 + MARGIN;
    start = 180 - ((1 - loc) * 45);
    stop  = start + 90;
    l *= 2;
    gdImageChar(img, font, xc , yc, 'C', FONT_COLOR) ;
  }

  /* draw the loop */
  gdImageArc(img, xc, yc, l,LOOP_H, start, stop, LOOP_COLOR);


  /* now deal with the loop legend */
  if((legend1 = (char *)malloc(10)) == NULL){ error_fatal("memory", NULL); }
  if((legend2 = (char *)malloc(10)) == NULL){ error_fatal("memory", NULL); }

  (void)sprintf(legend1, "Ll = %d", length);
  (void)sprintf(legend2, "KR = %d", kr);

  /* legend position */
  if (first_last == INNER) {
    xf1 = xc - strlen(legend1)*font->w/2;
    xf2 = xc - strlen(legend2)*font->w/2;
  }
  else if (first_last == FIRST) {
    xf1 = xc ;
    xf2 = xc ;
  }
  else {
    xf1 = xc - strlen(legend1)*font->w;
    xf2 = xc - strlen(legend2)*font->w;
  }
  yf1 = Y_CENTER - (loc * (LOOP_H + SEG_H - font->h/2));
  yf2 = yf1 + 2*font->h/2;

  /* plot the loop legend */
  gdImageString(img, font, xf1, yf1, (unsigned char *)legend1, FONT_COLOR);
  gdImageString(img, font, xf2, yf2, (unsigned char *)legend2, FONT_COLOR);

  free(legend1);
  free(legend2);

  return;
}

/* segment representation routine */
static void draw_segment(gdImagePtr img, int pos, int n, int kind) {

  /* segment position */
  int x;
  int x1, x2, y1, y2;
  int w, l;
  int color;

  /* legend related */
  int xf, yf;
  int num_l, num_h;
  char *tag;
  gdFontPtr font;
  font = gdFontLarge;

  x = (int)(pos * COEFF) + MARGIN;

  /* segment coordinate */
  x1 = x - (SEG_L/2);
  x2 = x + (SEG_L/2);
  y1 = Y_CENTER - (int)((SEG_H - SEG_L) / 2);
  y2 = Y_CENTER + (int)((SEG_H - SEG_L) / 2);
  w = l = SEG_L + 2;

  /* segment color attribution */
  if (kind == 1){ color = PUT_COLOR; }
  else { color = TM_COLOR; }

  /* draw the segment border and fill it*/
  gdImageLine(img, x1, y1, x1, y2, TM_BORDER);
  gdImageLine(img, x2, y1, x2, y2, TM_BORDER);
  gdImageArc(img, x, y1, w, l, 180, 0, TM_BORDER);
  gdImageArc(img, x, y2, w, l, 0, 180, TM_BORDER);
  gdImageFillToBorder(img, x, Y_CENTER,  TM_BORDER, color);

  /* tag segment with is number */
  if((tag = (char *)malloc(3*sizeof(char))) == NULL) {
    error_fatal("memory", NULL);
  }
  (void)sprintf(tag, "%d", n);

  num_l=(strlen(tag) * font->w) /2 ;
  num_h = font->h / 2;
  xf = x - num_l;
  yf = Y_CENTER - num_h;

  gdImageString(img, font, xf, yf, (unsigned char *)tag, FONT_COLOR);
  free(tag);
}

/* tweaking the loop length values in order to avoid some overlaping segments
 * drawing
 * segment are drawn with a fixed width, we must check that the loop
 * representation is long enough to avoid some segment collision.
 * eg: shorter loop lenght are forced to be represented at least with
 * a lenght compatible with segment width
 * see below
 *          _                    ____
 *         | |                  |    |
 *        _____                ___   ___
 *       / / \ \              /   \ /   \
 *      | |   | |            |     |     |
 *      | |   | |            |     |     |
 *       \_\_/_/              \___/ \___/
 *
 * loop length < SEG_L   loop_lenght forced to SEG_L
 */

static int *loop_tweaker (elprint_t *element, int n){

  int check , cont, graph_len;
  int i, j, ll;
  int *loop_length;
  i = n - ((n-1)/2);
  if ((loop_length = malloc(i * sizeof(int))) == NULL) {
    error_fatal("memory", NULL);
  }
  for (i=0, j=0; i<n; i+=2, j++  ) {
    loop_length[j] = element[i].lo.len;
  }


  check = 10;
  j = n - ((n-1)/2);
  while (check > 0) {
    cont = 1;
    check--;
    graph_len = 0;

    for (i = 0; i< j; i++) { graph_len += loop_length[i]; }
    if ( graph_len == 0) graph_len = 1;
    COEFF = (float)(800 - 40)/ graph_len;
    for (i = 0; i<j; i++) {
      ll = loop_length[i]*COEFF;
      if (ll < 20) {
	cont = 0;
	loop_length[i] = (int)(20/COEFF)+1;
      }
    }
    if (cont == 1) {break;}
  }

  if (check <= 0) {
    free (loop_length);
    return NULL;
  }

  return loop_length;
}

/* prepare the png file wich will host the image */
static FILE *file_maker(char *title, char *dir, int num) {

  char *file_name, *p;
  int path_len;
  FILE *pngout;

  /* preparing filename  and allocating the result file ptr*/
  path_len =sizeof(char) * (strlen(dir) + 1 + strlen(title) + 6 + 5 + 2 + 1 );
  if((file_name=(char *)malloc(path_len)) == NULL){
    error_fatal("memory", NULL);
  }
  p = file_name;
  (void)sprintf(p, "%s/%s-%d.png",dir, title, num);

  if ((pngout = fopen(file_name, "wb")) == NULL) {
    error_fatal("topos", NULL);
  }

  free(file_name);
  return pngout;
}


/* allocate the image holder, and plot the legends and common drawing */
static gdImagePtr init_img(char *title, int num) {
  int y1, y2;
  int xf, yf, text_l, text_h;

  char  *cyto, *extra;
  char *struct_num;
  char *certain, *put;
  char *loop, *kr;

  gdImagePtr img;
  gdFontPtr font, sfont;

  img = gdImageCreate(IMAGE_LENGTH, IMAGE_HEIGTH);
  font =  gdFontGiant;
  sfont = gdFontLarge;


  /* colors are given in RGB see /usr/lib/X11/rgb.txt */
  TM_COLOR = gdImageColorAllocate(img, 255, 255, 255);            /* white  */
  TM_BORDER = FONT_COLOR = gdImageColorAllocate(img, 0, 0, 0);    /* black  */
  MBR_COLOR = gdImageColorAllocate(img, 219, 219, 219);           /* grey86 */
  MBR_BORDER = LOOP_COLOR = gdImageColorAllocate(img, 3, 3, 3);   /* grey1  */
  PUT_COLOR = gdImageColorAllocate(img, 161, 161, 161);           /* grey63 */

  y1 = Y_CENTER - (int)(SEG_H / 4);
  y2 = Y_CENTER + (int)(SEG_H / 4);

  /* plot the membrane */
  gdImageFilledRectangle(img, 0, y1, IMAGE_LENGTH , y2, MBR_COLOR);
  gdImageLine(img, 0, y1, IMAGE_LENGTH, y1, MBR_BORDER);
  gdImageLine(img, 0, y2, IMAGE_LENGTH, y2, MBR_BORDER);

  /* title and legend */

  /* localisation legend */
  cyto = "CYTOPLASM";
  extra = "EXTRACELLULAR";
  text_h = font->h / 2;

  /* cytoplasm legend */
  text_l = (strlen(cyto) * font->w) /2 ;
  xf = X_CENTER - text_l;
  yf = (int)(IMAGE_HEIGTH * 2.0/9.0) - text_h;
  gdImageString(img, font, xf, yf, (unsigned char *)cyto, FONT_COLOR);

  /* extracellular legend */
  text_l = (strlen(extra) * font->w) /2 ;
  xf = X_CENTER - text_l;
  yf = (int)(IMAGE_HEIGTH * 8.0/9.0) - text_h;
  gdImageString(img, font, xf, yf, (unsigned char *)extra, FONT_COLOR);

  /* title generation */
  gdImageString(img, font, 10, 10, (unsigned char *)title, FONT_COLOR);
  if((struct_num = (char *)malloc(14+3)) == NULL){
    error_fatal("memory", NULL);
  }
  (void)sprintf(struct_num, "Structure no. %d", num);

  gdImageString(img, font, 10, 30, (unsigned char *)struct_num, FONT_COLOR);
  free(struct_num);

  /* color meaning legend */
  gdImageRectangle(img, 600, 30, 620, 50, TM_BORDER);
  gdImageFilledRectangle(img, 600, 60, 620, 80, PUT_COLOR);
  gdImageRectangle(img, 600, 60, 620, 80, TM_BORDER);

  certain = "Segment Certain";
  put = "Segment Putative";
  gdImageString(img, sfont, 630, 32, (unsigned char *)put, FONT_COLOR);
  gdImageString(img, sfont, 630, 62, (unsigned char *)certain, FONT_COLOR);

  /* cellular localisation legend */
  loop = "Ll: Loop length";
  kr   = "KR: Number of Lys and Arg";
  gdImageString(img, font, 10, 350, (unsigned char *)loop, FONT_COLOR);
  gdImageString(img, font, 10, 370, (unsigned char *)kr, FONT_COLOR);

  /* positionnning the Margin */
  return img;

}
