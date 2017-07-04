/*
   Log: 
   
   $Log: cluster-tissues.c,v $
   Revision 1.4  2002/03/07 06:27:08  default
   Use factor analyzers

   Revision 1.3  2002/02/18 13:41:06  default
   Initial checkin

   Revision 1.2  2002/02/10 00:47:12  default
   Removed code which allocates memory for the string on the stack, and
   allocate it on the heap instead.
   
   Fix up some minor compiler warnings.

   Notes:

   The following parameters worked for me:
   $head -10 alon_r.dat > test.dat
   $./cluster-tissues.exe test.dat 10 62
   factors=4
   components=4
   random starts=15
   kmeans starts=15
   random seed 1=1234
   random seed 2=2345
   random seed 3=5432
   
*/
/* cluster-tissues which uses EMMIX-t rather than factor analyzers */

/* functionality is the same, but */
/* don't give it more than 10 genes at a time! */
#define LINE_LENGTH 5000	/* max length of each line in microarray data */
#define MAX_GENES 2200		/* max number of genes in microarray data */
#define MAX_TISSUES 200 	/* max number of tissues in microarray data */
#define MAX_ENTRY_LENGTH 20	/* max length of any single entry in the array */

#include <unistd.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>         /* calloc */
#include <string.h>
//~ #include "emmix_win32.h"

/* filename, number of clusters, number of rows, number of
   vars/dimensions, number of factor
   dimensions, number of groups, number of random & k-means starts,
   random seeds */
int f1, g, r, k, s1, s2, s3, clus;

FILE* fOUTPUT;

void call_emfac (char *data_file, int n, int p);
void call_emface (char *data_file, int n, int p);
void call_emmix  (char *data_file, int n, int p);

char lpszErrMsg[1024];

FILE* open_file_or_exit(char* filename, char* mode);

/* ---------------------------------- */
/* COMMON ROUTINES WITH CLUSTER-GENES */
/* ---------------------------------- */

FILE* open_file_or_exit(char* filename, char* mode)
{
    FILE* tmp;
    char  s[255];
    
    tmp = fopen(filename, mode);
    if (!tmp)
    {
        sprintf(s, "Cannot open %s", filename);
        perror(s);
        exit(1);
    }

    return tmp;
}

int interactive_input(int* opt, int* f1, int* g, 
    int* r, int* k, 
    int* s1, int* s2, int* s3, char* output_filename)
{
  printf ("Enter 1 to use factor analyzers, 2 to use EMMIX\n");
  scanf ("%d", opt);
  if (*opt == 1) { 
    printf ("How many factors?\n");
    scanf ("%d", f1); 
  }
  printf ("How many components?\n");
  scanf ("%d", g);
  printf ("How many random starts?\n");
  scanf ("%d", r);
  printf ("How many k-means starts?\n");
  scanf ("%d", k);
  printf ("Random Seed 1?\n");
  scanf ("%d", s1);
  printf ("Random Seed 2?\n");
  scanf ("%d", s2);
  printf ("Random Seed 3?\n");
  scanf ("%d", s3);
  printf ("Specify output filename?\n");
  scanf ("%s", output_filename);
  
  return 1;
}

void reorder(char* data_file, char* sorted_file)
/* Reorders a data file according to to the cluster-tissues output */
{
    char cmdline[1024];
    char output_file[1024];
    char current_dir[1024];
    char* output_dir;
    char* output_basename;
    FILE* batch_file;

    /*
    output_dir = dirname(data_file);
    output_basename = basename(data_file);
    */
    
    getcwd(current_dir, 1024);

    sprintf(output_file, "%s.reorder.dat", data_file); 
    sprintf(cmdline, ".\\bin\\python.exe EmmixGene\\reorder.py %s %s > %s", 
        data_file,
        sorted_file,
        output_file);

    /* This bit works
    ** sprintf(cmdline, "C:\\python21\\python.exe Hello.py > tmp.out");
    */

    /* Write out to a batch file */
    batch_file = fopen("reorder.bat","w");
    if (batch_file)
    {
        fputs(cmdline, batch_file);
        fputs("\n", batch_file);
        fclose(batch_file);
    }
    
    // system("reorder.bat");
    system(cmdline);
    fprintf(stderr, "DEBUG: Done system cmdline");
    fflush(stderr);

}

int transpose_file(char* in_file, char* out_file, int n, int p)
{
    FILE *f;
    int   i,j;
    char* buf;
    char* s[MAX_GENES][MAX_TISSUES];
    
    /* allocate character array on the heap */
    buf = (char*) calloc(MAX_GENES*MAX_TISSUES*MAX_ENTRY_LENGTH,sizeof(char));
    if (buf==NULL){
      perror("Could not allocate memory for input buffer\n");
      exit (1);
    }
    for (i=0; i<MAX_GENES;i++){
      for (j=0; j<MAX_TISSUES;j++){
          s[i][j]=(char*) (buf + (i*MAX_TISSUES+j) * MAX_ENTRY_LENGTH);
      }
    }

    /* transpose the input file */
    f = open_file_or_exit(in_file, "r");
    for (i = 0; i < n; i++)
        for (j = 0; j < p; j++)
          fscanf (f, "%s", s[i][j]);
    fclose (f);

    f = open_file_or_exit(out_file, "w");
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < n; j++)
            fprintf (f, "%s ", s[j][i]);
        fprintf (f, "\n");
    }
    fclose (f);

    free(buf);
    return 1;
}

int main (int argc, char **argv)
{
  int n;			/* number of rows - 2nd argument */
  int p;			/* number of cols - 3rd argument */
  int tmp, opt;

  /*
  We get stack fault on Win98 with GCC2.95.3
  Trying to allocate 168 Mb on the stack.
  char s[MAX_GENES][MAX_TISSUES][MAX_ENTRY_LENGTH];
  */
  
  char data_file[1024];
  char trans_file[1024];
  char output_filename[1024];   /* output filename */
  
  if (argc != 4)
    {
      printf
        ("usage: cluster-tissues microarray-filename number-of-rows number-of-columns\n");
        exit (1);
    }
  
  
  strcpy(data_file, argv[1]);
  n = atoi (argv[2]);
  p = atoi (argv[3]);
  
  /* do not change the .trans file name */
  /* HARD WIRED into other routines     */
  sprintf(trans_file, "%s.trans", data_file);
  
  /* Get input from user */
  interactive_input( &opt, &f1, &g, &r, &k, &s1, &s2, &s3, output_filename);

  transpose_file(data_file, trans_file, n,p);

  tmp = n;
  n = p;
  p = tmp;

  /* open the output file */
  fOUTPUT = open_file_or_exit(output_filename, "w");
  fprintf(stderr, "DEBUG: Output file opened\n");
  fflush(stderr);

  if (opt == 1)
  { 
    if (n > p)
        call_emfac (data_file, n, p);
    else
        call_emface (data_file, n, p); 
  }
  else
  {
    call_emmix (data_file, n, p);
  }

  /* close the output file */
  fclose(fOUTPUT);
  
  fprintf(stderr, "DEBUG: Output file closed\n");
  fflush(stderr);
  
  fprintf(stderr, "DEBUG: Reordering input file\n");
  fflush(stderr);
  reorder(data_file, output_filename);
  fprintf(stderr, "DEBUG: Finished reordering input file\n");
  fflush(stderr);
  
  /* unlink(trans_file); */
    
  return 0;
}

void
call_emfac (char* data_file, int n, int p)
{
  FILE *f;
  int i, j;
  char s[LINE_LENGTH];
  sprintf (s, "tmp.%d", getpid ());
  f = open_file_or_exit(s, "w");
  fprintf (f, "%s.trans\n", data_file);
  fprintf (f, "%s.out\n", data_file);
  fprintf (f, "%s.m\n", data_file);
  fprintf (f, "%s.mod\n", data_file);
  fprintf (f, "%d\n", n);
  fprintf (f, "%d\n", p);
  fprintf (f, "%d\n", p);
  fprintf (f, "0\n");		/* DON'T CLUSTER COLUMNS */
  fprintf (f, "%d\n", f1);
  fprintf (f, "0\n");
  fprintf (f, "%d\n", g);
  fprintf (f, "2\n3\n");
  fprintf (f, "%d\n100\n%d\n", r, k);
  fprintf (f, "1\n0\n%d\n%d\n%d\n", s1, s2, s3);
  fclose (f);
  //sprintf (s, "%s < tmp.%d > %s", EMMIXF1_EXE, getpid(), NUL);
  sprintf (s, "%s < tmp.%d > %s", EMMIXF1_EXE, getpid(), "EMMIX-f1.out");
  system (s);
  sprintf (s, "%s.out", data_file);
  f = open_file_or_exit(s, "r");
  while (!strstr (s, "Implied Grouping") && !feof(f))
    fgets (s, LINE_LENGTH, f);
  if (!strstr(s,"Implied Grouping"))
  {
    perror("ERROR: Could not determine implied grouping");
    goto exit_call_emfac;
  }
  
  /* Output results */
  for (i = 0; i < n; i++)
    {
      fscanf (f, "%d", &j);
      fprintf ( fOUTPUT, "%d  ", j);
      if (i > 0 && (i % 10) == 9)
	fprintf (fOUTPUT, "\n");
    }
  fprintf (fOUTPUT, "\n");

  exit_call_emfac:
  sprintf (s, "%s.out", data_file);
  unlink (s);
  sprintf (s, "%s.m", data_file);
  unlink (s);
  sprintf (s, "%s.mod", data_file);
  unlink (s);
  sprintf (s, "tmp.%d", getpid ());
  unlink (s);
}

void
call_emface (char* data_file, int n, int p)
{
  FILE *f;
  int i, j;
  char s[LINE_LENGTH];
  sprintf (s, "tmp.%d", getpid ());
  f = open_file_or_exit(s, "w");
  fprintf (f, "%s.trans\n", data_file);
  fprintf (f, "%s.out\n", data_file);
  fprintf (f, "%s.m\n", data_file);
  /*fprintf(f,"%s.mod\n",argv[1]); */
  fprintf (f, "%d\n", n);
  fprintf (f, "%d\n", p);
  fprintf (f, "%d\n", p);
  /*fprintf(f,"0\n"); */
  fprintf (f, "%d\n", f1);
  fprintf (f, "0\n");
  fprintf (f, "%d\n", g);
  fprintf (f, "2\n3\n");
  fprintf (f, "%d\n100\n%d\n", r, k);
  fprintf (f, "0\n0\n%d\n%d\n%d\n", s1, s2, s3);
  fclose (f);
  //sprintf (s, "%s < tmp.%d > %s", EMMIXF2_EXE, getpid(), NUL);
  sprintf (s, "%s < tmp.%d > %s", EMMIXF2_EXE, getpid(), "EMMIX-f2.out");
  system (s);
  sprintf (s, "%s.out", data_file);
  f = open_file_or_exit(s, "r");
  while (!strstr (s, "Final Allocation") && !feof(f))
    fgets (s, LINE_LENGTH, f);
  if (!strstr(s, "Final Allocation"))
  {
    perror("Could not determine final allocation");
    goto exit_call_emface;
  }
  fgets (s, LINE_LENGTH, f);
  for (i = 0; i < n; i++)
    {
      fscanf (f, "%d", &j);
      fprintf (fOUTPUT, "%d  ", j);
      if (i > 0 && (i % 10) == 9)
	fprintf (fOUTPUT, "\n");
    }
  fprintf (fOUTPUT, "\n");

  exit_call_emface:
  sprintf (s, "%s.out", data_file);
  unlink (s);
  sprintf (s, "%s.m", data_file);
  unlink (s);
  sprintf (s, "tmp.%d", getpid ());
  unlink (s);
}

void call_emmix (char* data_file, int n, int p)
{
  FILE *f;
  int i, j;
  char s[LINE_LENGTH];
  char emmix_in[LINE_LENGTH];
  char data_out[LINE_LENGTH];

  sprintf (emmix_in, "tmp.%d", getpid ());  /* stdin */
  sprintf (data_out, "%s.out", data_file); /* emmix results */

  /* Prepare stdin to emmix.exe */
  f = open_file_or_exit(emmix_in, "w");
  fprintf (f, "2\n");
  fprintf (f, "%s.trans\n", data_file);
  fprintf (f, "%s\n", data_out);
  fprintf (f, "%d\n", n);
  fprintf (f, "%d\n", p);
  fprintf (f, "%d\n", p);
  fprintf (f, "%d\n", g);
  fprintf (f, "2\n3\n");
  fprintf (f, "%d\n100\n%d\nn\n", r, k);
  fprintf (f, "%d\n%d\n%d\n", s1, s2, s3);
  fclose (f);
  
  //sprintf (s, "%s < tmp.%d > %s", EMMIX_EXE, getpid (), NUL);
  sprintf (s, "%s < %s > %s", EMMIX_EXE, emmix_in, "EMMIX.out");
  system (s);
  
  f = open_file_or_exit(data_out, "r");
  while (!strstr (s, "Implied grouping") & !feof(f))
    fgets (s, LINE_LENGTH, f);

  if (!strstr(s, "Implied grouping"))
  {
    perror("ERROR: Cannot determine implied grouping");
    goto exit_call_emmix;
  }

  for (i = 0; i < n; i++)
    {
      fscanf (f, "%d", &j);
      fprintf (fOUTPUT, "%d  ", j);
      if (i > 0 && (i % 10) == 9)
	fprintf (fOUTPUT, "\n");
    }
  fprintf (fOUTPUT, "\n");

  exit_call_emmix:
  fclose(f);
  /*
  unlink (data_out);
  unlink (emmix_in);
  */

}

