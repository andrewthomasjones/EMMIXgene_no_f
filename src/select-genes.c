/* 
    $RCSfile: select-genes.c,v $
    
    $Revision: 1.7 $

    Change Log
    $Log: select-genes.c,v $
    Revision 1.7  2002/04/02 11:24:01  default
    Fixed apparent error in writing to the wrong file descriptor

    Revision 1.6  2002/03/10 00:04:41  default
    Better error handling

    Revision 1.5  2002/03/07 01:05:35  default
    Output file called *.cut.stats instead of *.stats

    Revision 1.4  2002/02/18 13:41:07  default

    Revision 1.3  2002/02/12 14:25:26  default
    Used emmix_win32.h header

    Revision 1.2  2002/02/10 14:08:54  default
    Conditionally #define Win32isms for system()
    Do not allocate big character array on the stack

*/

#define LINE_LENGTH 1000	/* max length of each line in microarray data */
#define MAX_GENES 3800		/* max number of genes in microarray data */
#define MAX_TISSUES 75		/* max number of tissues in microarray data */
#define MAX_SELECTED_GENES 2200	/* max number of selected genes */

#include <unistd.h>		    /* getpid, unlink */
#include <stdlib.h>		    /* exit, system */
#include <stdio.h>		    /* printf, scanf */
#include <string.h>	        /* strstr */
//#include "emmix_win32.h"

void call_emmix (char **argv, char *t, float *tl, int *sm, int *er, int col,
		 int g);
void call_emd (char **argv, int col, int g, int clus, int genecount);
void strip_cr (char s[LINE_LENGTH], char t[LINE_LENGTH]);
void removeTempFiles(void);

struct tab
{
  int nu;
  int sm;
  float tl;
};

int r, k, s1, s2, s3;		/* random starts, k-means starts, seeds 1 to 3 */

static char datafile[255];  /* name of input file, used for cleanup routine */

static char lpszErrMsg[1024]; /* holds error messages */

int reord;
int
tabcompare (const void *p1, const void *p2)
{
  float i = (*((struct tab *) p1)).tl;
  float j = (*((struct tab *) p2)).tl;
  if (i < j)
    return 1;
  if (i > j)
    return -1;
  return 0;
}

int main (int argc, char **argv)
{
  int row;			/* 2nd argument */
  int col;			/* 3rd argument */
  int g;			/* 4th argument */
  int x;			/* counter */
  int i;			/* another counter */
  int rc, b1, b2, r1, k1, clus;
  int genecount = 0;
  char s[LINE_LENGTH];
  struct tab gtab[MAX_SELECTED_GENES];	/* for our cut-down genes and -2 log \lambdas */

  char* sg[MAX_GENES];
  char* buf;
  FILE *f, *f2, *f3, *f4;

  /* CGT commented out 10/2/2002
  ** Could end up trying to allocate 3Mb of space
  ** on the stack. Will overflow typical Win32 box
  **
  char sg[MAX_GENES][LINE_LENGTH]; */	/* each line up to LINE_LENGTH chars 
					   long, max MAX_GENES lines */


  if (argc != 4)
    {
      printf
	("usage: select-genes microarray-filename number-of-rows number-of-columns\n");
      exit (1);
    }
  col = atoi (argv[3]);
  row = atoi (argv[2]);

  /* Allocate string array on the heap */
  /* Free()'d on exit                  */
  buf = (char*) calloc(MAX_GENES * LINE_LENGTH, sizeof(char));
  if (buf == NULL)
  {
    perror ("Could not allocate memory for SG");
    exit (1);
  }
  for (i=0;i<MAX_GENES;i++) sg[i]=(char*) (buf + i * LINE_LENGTH);

  strcpy(datafile, argv[1]);
  g = 1;
  /* rm -f error.$1.$$ log.$1.$$ smaller.$1.$$ input1.$$ */
  /* sprintf(s,"error.%d",getpid());
     unlink(s);
     sprintf(s,"log.%d",getpid());
     unlink(s);
     sprintf(s,"smaller.%d",getpid());
     unlink(s); */

  /* Comment out when debugging */
  /* atexit(removeTempFiles); */

  printf
    ("How many random starts for the fitting of t components to individual genes?\n");
  scanf ("%d", &r);
  printf
    ("How many k-means starts for the fitting of t components to individual genes?\n");
  scanf ("%d", &k);
  printf ("Enter random seed 1:\n");
  scanf ("%d", &s1);
  printf ("Enter random seed 2:\n");
  scanf ("%d", &s2);
  printf ("Enter random seed 3:\n");
  scanf ("%d", &s3);
  printf ("Enter threshold for likelihood ratio statistic:\n");	/* -2 log \lambda cutoff */
  scanf ("%d", &b1);
  printf ("Enter threshold for minimum cluster size:\n");	/* s_{min} cutoff */
  scanf ("%d", &b2);
  sprintf (s, "input1.%d", getpid ());
  f = fopen (s, "w");
  if (!f)
  {
    perror (s);
    exit(1);
  }
  sprintf (s, "input3.%d", getpid ());
  f2 = fopen (s, "w");
  if (!f2)
  {
    perror (s);
    exit(1);
  }
  fprintf (f, "2\n");		/* echo 2 >> input1.$$ */
  fprintf (f, "tmp.%d\n", getpid ());	/* echo tmp.$$ >> input1.$$ */
  fprintf (f, "out1.%d\n", getpid ());	/* echo out1.$$ >> input1.$$ */
  fprintf (f, "%d\n", col);	/* echo $3 >> input1.$$ */
  fprintf (f, "1\n");		/* echo 1 >> input1.$$ */
  fprintf (f, "1\n");		/* echo 1 >> input1.$$ */

  fprintf (f2, "2\n");		/* echo 2 >> input1.$$ */
  fprintf (f2, "tmp.%d\n", getpid ());	/* echo tmp.$$ >> input1.$$ */
  fprintf (f2, "out1.%d\n", getpid ());	/* echo out1.$$ >> input1.$$ */
  fprintf (f2, "%d\n", col);	/* echo $3 >> input1.$$ */
  fprintf (f2, "1\n");		/* echo 1 >> input1.$$ */
  fprintf (f2, "1\n");		/* echo 1 >> input1.$$ */

  fprintf (f, "1\ny\n9\n2\n.5\n0\n");

  fprintf (f2, "%d\n2\n3\n%d\n100\n%d\ny\n9\n2\n", 2, r, k);
  x = 0;
  while (x++ < 2)
    fprintf (f2, "1\n");
  fprintf (f2, "0\n%d\n%d\n%d\n", s1, s2, s3);

  fclose (f);
  fclose (f2);

  sprintf (s, "input2.%d", getpid ());
  f = fopen (s, "w");
  if (!f)
    perror (s);
  sprintf (s, "input4.%d", getpid ());
  f2 = fopen (s, "w");
  if (!f2)
    perror (s);

  fprintf (f, "2\n");
  fprintf (f, "tmp.%d\n", getpid ());
  fprintf (f, "out2.%d\n", getpid ());
  fprintf (f, "%d\n", col);
  fprintf (f, "1\n");
  fprintf (f, "1\n");
  fprintf (f, "%d\n", 1 + 1);
  fprintf (f, "2\n");
  fprintf (f, "3\n");
  fprintf (f, "%d\n", r);
  fprintf (f, "100\n");
  fprintf (f, "%d\n", k);
  fprintf (f, "y\n");
  fprintf (f, "9\n");
  fprintf (f, "2\n");

  fprintf (f2, "2\n");
  fprintf (f2, "tmp.%d\n", getpid ());
  fprintf (f2, "out2.%d\n", getpid ());
  fprintf (f2, "%d\n", col);
  fprintf (f2, "1\n");
  fprintf (f2, "1\n");
  fprintf (f2, "%d\n", 2 + 1);
  fprintf (f2, "2\n");
  fprintf (f2, "3\n");
  fprintf (f2, "%d\n", r);
  fprintf (f2, "100\n");
  fprintf (f2, "%d\n", k);
  fprintf (f2, "y\n");
  fprintf (f2, "9\n");
  fprintf (f2, "2\n");

  x = 0;
  while (x++ <= 1)
    {
      fprintf (f, "1\n");
    }
  x = 0;
  while (x++ <= 2)
    {
      fprintf (f2, "1\n");
    }

  fprintf (f, "0\n");
  fprintf (f, "%d\n%d\n%d\n", s1, s2, s3);

  fprintf (f2, "0\n");
  fprintf (f2, "%d\n%d\n%d\n", s1, s2, s3);
  fclose (f);
  fclose (f2);

  rc = 0;
  sprintf (s, "%s.cut.stats", argv[1]);
  f3 = fopen (s, "w");
  if (!f3)
  { 
    sprintf(lpszErrMsg, "ERROR: Cannot open %s", s);
    perror (lpszErrMsg);
    exit(1);
  }

  f = fopen (argv[1], "r");
  if (!f)
  {
    sprintf(lpszErrMsg, "ERROR: Cannot open %s", argv[1]);
    perror (lpszErrMsg);
    exit(1);
  }
  
  while (rc < row)
  {
      char t[LINE_LENGTH], last;
      int sm, er;
      float tl;

      /* Always clear memory */
      memset(t, 0, LINE_LENGTH);

      /* grab line i, change space to CR, remove blank lines */
      /* sed -n $i'p' $1 | tr ' ' '\n' | sed '/^$/d' > tmp.$$ */
      /* write output to tmp.$$ */
      
      fgets (s, LINE_LENGTH, f);

      strip_cr (s, t);
      rc++;

      /* Announce where we are at */
      fprintf(stderr, "ANNOUNCE: step %d of %d\n", rc, row);
      fflush(stderr);

      call_emmix (argv, t, &tl, &sm, &er, col, g);	/* get -2 log \lambda, smaller, error */

      if (tl > (float) b1)		/* if -2 log \lambda exceeds b_1 cutoff */
      {
        if (sm >= b2)		/* if smaller groups exceeds b_2 cutoff */
        {
            fprintf (f3, "%d %f %d\n", rc, tl, sm);
            fflush (f3);
            gtab[genecount].nu = genecount;
            gtab[genecount].tl = tl;
            gtab[genecount].sm = sm;
            strcpy (sg[genecount], s);
            genecount++;
            
        } else {
        
            /* g is the number of groups */
            call_emmix (argv, t, &tl, &sm, &er, col, g + 1);
            if (tl > b1)
            {
              fprintf (f3, "%d %f %d\n", rc, tl, sm);
              fflush (f3);
              gtab[genecount].nu = genecount;
              gtab[genecount].tl = tl;
              gtab[genecount].sm = sm;
              strcpy (sg[genecount], s);
              genecount++;
            }
        }

     } else {

        fprintf(stderr, "DEBUG: !(tl>b1) !(%f>%d)\n", tl, b1);
        fflush(stderr);

     } /* if (t1 > b1) */

  } /* while */
  
  fclose (f);
  fclose (f3);
  qsort ((void *) gtab, genecount, sizeof (struct tab), tabcompare);

  sprintf (s, "%s.cut", argv[1]);
  f4 = fopen (s, "w");
  if (!f4)
    {
    sprintf(lpszErrMsg, "ERROR: cannot open %s", s);
    perror (lpszErrMsg);
    exit(1);
    }

  fprintf(stderr, "DEBUG: genecount = %d\n", genecount);
  fflush(stderr);
  
  rc = 0;
  while (rc < genecount)
    fputs (sg[gtab[rc++].nu], f4);
  fclose (f4);

  sprintf (s, "%s.cut.sstats", argv[1]); /* sorted statistics */
  f3 = fopen (s, "w");
  if (!f3)
  {
    sprintf(lpszErrMsg, "ERROR: Cannot open %s", s);
    perror (lpszErrMsg);
    exit(1);
  }
  rc = 0;
  while (rc < genecount)
  {
    fprintf (f3,"%d %f %d\n",gtab[rc].nu+1, gtab[rc].tl, gtab[rc].sm);
    rc++;
  }
  fclose (f3);

  free(buf);
  return 0;
}

void removeTempFiles()
{

    char s[1024];

    /* clean up */
    sprintf (s, "input1.%d", getpid ());
    unlink (s);
    sprintf (s, "input2.%d", getpid ());
    unlink (s);
    sprintf (s, "input3.%d", getpid ());
    unlink (s);
    sprintf (s, "input4.%d", getpid ());
    unlink (s);
    sprintf (s, "%s.tmp", datafile);
    unlink (s);
    sprintf (s, "out1.%d", getpid ());
    unlink (s);
    sprintf (s, "out2.%d", getpid ());
    unlink (s);
    sprintf (s, "%sm", datafile);
    unlink (s);
    sprintf (s, "tmp.%d", getpid ());
    unlink (s);

}

void
call_emmix (char **argv, char *t, float *tl, int *sm, int *er, int col, int g)
{
  FILE *f2;
  float fl1, fl2;
  int group1[MAX_TISSUES], group2[MAX_TISSUES];
  int comp[MAX_TISSUES], x, y = 0, a1, a2, i;
  char s[LINE_LENGTH];

  sprintf (s, "tmp.%d", getpid ());
  f2 = fopen (s, "w");
  if (!f2)
  {
    sprintf(lpszErrMsg, "ERROR: Cannot open %s", s);
    perror (lpszErrMsg);
    exit(1);
  }
  fputs (t, f2);
  fprintf (f2, "\n");
  fclose (f2);

  /* flush input and output streams before calling system() */
  fflush(stdout);
  fflush(stdin);

  //~ if (g == 1)
    //~ {
      //~ sprintf (s, "%s < input1.%d > %s ", EMMIXT_EXE, getpid (), NUL);
      //~ system (s);
      //~ sprintf (s, "%s < input2.%d > %s ", EMMIXT_EXE, getpid (), NUL);
      //~ system (s);
    //~ }
  //~ else
    //~ {
      //~ sprintf (s, "%s < input3.%d > %s ", EMMIXT_EXE, getpid (), NUL);
      //~ system (s);
      //~ sprintf (s, "%s < input4.%d > %s ", EMMIXT_EXE, getpid (), NUL);
      //~ system (s);
    //~ }

  sprintf (s, "out1.%d", getpid ());
  f2 = fopen (s, "r");
  if (!f2)
  {
    perror (s);
    exit(1);
  }
  while (!strstr (s, "Log-L") && !feof (f2))
    fgets (s, LINE_LENGTH, f2);
  a1 = sscanf (s, " Final Log-Likelihood is %f", &fl1);
  fclose (f2);

  sprintf (s, "out2.%d", getpid ());
  f2 = fopen (s, "r");
  if (!f2)
  {
    perror (s);
    exit(1);
  }
  while (!strstr (s, "Log-L") && !feof (f2))
    fgets (s, LINE_LENGTH, f2);
  a2 = sscanf (s, " Final Log-Likelihood is %f", &fl2);
  fclose (f2);
  if (a1 && a2)
    *tl = 2.0 * (fl2 - fl1);
  else
    {
      *tl = -100;
      *sm = -100;
      *er = -100;
      return;
    }

  /* calculate error rate */
  /*sprintf (s, "%s-group", argv[1]);
     f2 = fopen (s, "r");
     if (!f2) perror(s);
     for (i = 0; i < col; i++)
     fscanf (f2, "%d", &group1[i]);
     fclose (f2); */

  /* The following reads out2.--- and looks for 
  ** line matching "Implied" and reads in 
  ** col number of values into group2[]
  */
  
  sprintf (s, "out2.%d", getpid ());
  f2 = fopen (s, "r");
  if (!f2)
  {
    sprintf(lpszErrMsg, "ERROR: Cannot open %s", s);
    perror (lpszErrMsg);
    exit(1);
  }

  while (!strstr (s, "Implied") && !feof (f2))
    fgets (s, LINE_LENGTH, f2);
  for (i = 0; i < col; i++)
    fscanf (f2, "%d", &group2[i]);

  fclose (f2);
  *er = 0;
  /* for (i = 0; i < col; i++)
     if (group1[i] != group2[i])
     (*er)++;
     if (col - *er < *er)
     *er = col - *er; */

  memset (comp, 0, sizeof (comp));
  for (i = 0; i < col; i++)
    comp[group2[i]]++;
  *sm = 9999;
  for (i = 1; i <= g + 1; i++)
    if (comp[i] < *sm)
      *sm = comp[i];

  /* 06/04/2002 CGT GJM 
  ** If minimum size is 1, the set -2 log lambda t = -100
  */
  if (*sm == 1)
  {
      *tl = -100;
      *er = -100;
  }

  /* For 2 vs 3, count the number of groups which
     sample size (sm) >= the threshold sample size
     eg. 1 1 1 2 2 2 2 3 3 3 3 3 ->
         3     4       5
     if threshold is 4, then number of values >= 4 is 2.
     Condition if 2 >= 2.
  */
}

/* strip_cr
** Replaces spaces and tabs with newline
** Parameters:
**  [in]    char* s
**  [out]   char* t
*/
void
strip_cr (char s[LINE_LENGTH], char t[LINE_LENGTH])
{
  int x, y;
  char last;
  y = 0;
  for (x = 0; x < strlen (s); x++)
    {
      last = s[x];
      if (s[x] == ' ' || s[x] == '\t')
        {
          t[y++] = '\n';
          continue;
        }
      if (s[x] == '\n' && last == '\n')
    	continue;
      t[y++] = s[x];
    }
  t[y] = 0;			/* NUL character to end string */
}

