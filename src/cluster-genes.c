/*

    Department of Mathematics
    University Of Queensland
    Copyright (c) 2001-2002 All Rights Reserved
    
    Revision: $Revision: 1.11 $

    Change Log:
    $Log: cluster-genes.c,v $
    Revision 1.11  2002/05/06 01:46:01  default
    Bug fix: Deletes existing _groupmeans and _list files before appending

    Revision 1.10  2002/05/04 05:15:01  default
    Big refactoring ... getFinalLogLikehood, open_or_fail, temp filenames

    Revision 1.9  2002/04/16 11:25:01  default
    Send output to EMMIX-t.out, EMMIX-spher.out etc.
    Replace references to argv[1] with data_file

    Revision 1.8  2002/04/13 10:22:20  default
    Refactored makeInput1

    Revision 1.7  2002/04/13 10:01:38  default
    Refacted user_input
    Cleaned up error handling routines

    Revision 1.6  2002/03/07 06:30:12  default
    Removed test code

    Revision 1.5  2002/03/07 02:56:26  default
    Check all input files present before continuing
    More robust error handling

    Revision 1.4  2002/02/12 14:25:10  default
    Used emmix_win32.h header

    Revision 1.3  2002/02/10 00:53:35  default
    Used WIN32 style commandline for shelling out to EMMIX
    Removed the code which allocated a large string on the stack


*/

/* vim:set nowrap: */

#define LINE_LENGTH 5000	/* max length of each line in microarray data */
#define MAX_GENES 5000		/* max number of genes in microarray data */
#define MAX_TISSUES 200     /* max number of tissues in microarray data */
#define MAX_GROUPS 100		/* max number of groups to split into */
#define MAX_ENTRY_LENGTH 20	/* max length of any single entry in the array */
#define ENDL 10             /* new line character */
#define USE_TEMP_FILES 0    /* send groupmeans etc to temp files first */

#include <unistd.h>		/* getpid, unlink */
#include <stdlib.h>		/* exit, system */
#include <stdio.h>		/* printf, scanf */
#include <string.h>		/* strtok */
//~ #include "emmix_win32.h"
//~ #include <windows.h>

typedef char FILENAME[255];
FILENAME* group_filename;
FILENAME* list_filename;
FILENAME groupmean_filename;

void  call_emmix (char *data_file, char *t, float *tl, int *sm, int *er, int col,
		 int g);
void call_emd (char *data_file, int col, int g, int clus, int row, 
          int r, int k, int r1, int k1, float b2, 
          int nu, int s1, int s2, int s3, 
          int reord);
void  strip_cr (char s[LINE_LENGTH], char t[LINE_LENGTH]);
FILE* open_file_or_exit(char* filename, char* mode);
void  makeInput1(int,int,int,int,int,float,int,int,int,int,int,int);
void  initFilenames(char*, int);
void  makeTempFilesPermanent(char* data_file, int nGroups);

int errno;

/* void makeInput2(); */

struct tab
{
  int nu;
  int count;
  float tl;
  char s[LINE_LENGTH];
};


char datafile[255];

int debug_trace(char* msg)
{
    printf("XXX %s\n", msg);
    return 1;
}

/* initFilenames
** Initializes group_filename[], list_filename[], groupmean_filename[]
** 
*/
void initFilenames(char* data_file, int nGroups)
{

    int i;
    char lpstrTempPath[MAX_PATH];
    
    group_filename = calloc(nGroups+1, sizeof(FILENAME));
    list_filename  = calloc(nGroups+1, sizeof(FILENAME));
    
    GetTempPath(MAX_PATH,lpstrTempPath);

#if USE_TEMP_FILES 
    for (i=1; i<=nGroups; i++)
    {

        GetTempFileName(lpstrTempPath, "em", 0, group_filename[i]);
        GetTempFileName(lpstrTempPath, "em", 0, list_filename[i]);
        
    }

    GetTempFileName(lpstrTempPath, "em", 0, groupmean_filename);
#else
    for (i=1; i<=nGroups; i++)
    {
        sprintf(group_filename[i], "%s_group%d", data_file, i);
        sprintf(list_filename[i], "%s_list%d", data_file, i);
    }
    sprintf(groupmean_filename, "%s_groupmeans", data_file); 
#endif

}


/* unlinkOutputFiles
** clears the output files 
*/
void unlinkOutputFiles(int nGroups)
{
    int i;
    
    for (i=1; i<=nGroups; i++)
    {
        unlink(group_filename[i]);
        unlink(list_filename[i]);
    }

    unlink(groupmean_filename);

}


/* makeTempFilesPermanent
** moves the temporary files to their permanent names
*/
void makeTempFilesPermanent(char* data_file, int nGroups)
{
    
    int i;
    FILENAME one_group_filename;
    FILENAME one_list_filename;
    FILENAME new_groupmean_filename;
    
    for (i=1; i<=nGroups; i++)
    {
        /* 1. Define the new filename 
           2. Delete existing one if it exists
           3. Rename the original filename into the new filename
        */
        sprintf(one_group_filename, "%s_group%d", data_file, i);
        if (strcasecmp(group_filename[i], one_group_filename) != 0)
        {
            unlink(one_group_filename);
            if (rename(group_filename[i], one_group_filename) == -1)
                perror(strerror(errno));
        }

        sprintf(one_list_filename, "%s_list%d", data_file, i);
        if (strcasecmp(list_filename[i], one_list_filename) != 0)
        {
            unlink(one_list_filename);
            rename(list_filename[i], one_list_filename);
        }

    }

    sprintf(new_groupmean_filename, "%s_groupmeans", data_file); 
    if (strcasecmp(groupmean_filename, new_groupmean_filename) != 0)
    {
        unlink(new_groupmean_filename);
        rename(groupmean_filename, new_groupmean_filename);
    }

    free(group_filename);
    free(list_filename);
}

/* Gets cluster-genes parameters from the stdin */
int user_input(r,k,r1,k1,clus,b2,s1,s2,s3,reord,nu)
int *r,*k,*r1,*k1,*clus,*s1,*s2,*s3,*nu,*reord;
float *b2;
{

  char s[1024];
  printf
    ("How many random starts for the fitting of normal components to group means?\n");
  scanf ("%d", r);
  printf
    ("How many k-means starts for the fitting of normal components to group means?\n");
  scanf ("%d", k);
  printf
    ("How many random starts for the fitting of common spherical components to the selected genes?\n");
  scanf ("%d", r1);
  printf
    ("How many k-means starts for the fitting of common spherical components to the selected genes?\n");
  scanf ("%d", k1);
  printf ("How many groups into which to cluster the selected genes?\n");
  scanf ("%d", clus);

  /* printf ("Enter threshold for likelihood ratio statistic:\n"); 
  scanf ("%d", &B1); */
  printf ("Enter threshold for minimum cluster size:\n");       /* s_{min} cutoff */
  scanf ("%f", b2);


  printf ("Enter random seed 1:\n");
  scanf ("%d", s1);
  printf ("Enter random seed 2:\n");
  scanf ("%d", s2);
  printf ("Enter random seed 3:\n");
  scanf ("%d", s3);
  printf
    ("Do you wish to reorder the tissues as clustered on\nbasis of fitted group means?\n");
  scanf ("%s", s);
  if (s[0] == 'y' || s[0] == 'Y')
    *reord = 1;
  else
    *reord = 0;

  printf ("Do you wish use the robust version?\n");
  scanf ("%s", s);
  *nu = 0;
  if (s[0] == 'y' || s[0] == 'Y')
    {
      printf ("Enter common value of nu:\n");
      scanf ("%d", nu);
    }

  return 1; /* success */
}

void makeInput1(r,k,r1,k1,clus,b2,s1,s2,s3,reord,nu,col)
int  r,k,r1,k1,clus,s1,s2,s3,nu,reord,col;
float b2;
{

    FILE *f, *f2;
    int  x;
    char s[1024];
    
    sprintf (s, "input1.%d", getpid ());
    f = open_file_or_exit (s, "w");

    sprintf (s, "input3.%d", getpid ());
    f2 = open_file_or_exit (s, "w");
    fprintf (f, "2\n");        /* echo 2 >> input1.$$ */
    fprintf (f, "tmp.%d\n", getpid ());    /* echo tmp.$$ >> input1.$$ */
    fprintf (f, "out1.%d\n", getpid ());    /* echo out1.$$ >> input1.$$ */
    fprintf (f, "%d\n", col);    /* echo $3 >> input1.$$ */
    fprintf (f, "1\n");        /* echo 1 >> input1.$$ */
    fprintf (f, "1\n");        /* echo 1 >> input1.$$ */

    fprintf (f2, "2\n");        /* echo 2 >> input1.$$ */
    fprintf (f2, "tmp.%d\n", getpid ());    /* echo tmp.$$ >> input1.$$ */
    fprintf (f2, "out1.%d\n", getpid ());    /* echo out1.$$ >> input1.$$ */
    fprintf (f2, "%d\n", col);    /* echo $3 >> input1.$$ */
    fprintf (f2, "1\n");        /* echo 1 >> input1.$$ */
    fprintf (f2, "1\n");        /* echo 1 >> input1.$$ */

    fprintf (f, "1\ny\n9\n2\n.5\n0\n");

    fprintf (f2, "%d\n2\n3\n%d\n100\n%d\ny\n9\n2\n", 2, r, k);
    x = 0;
    while (x++ < 2)
        fprintf (f2, "1\n");
    fprintf (f2, "0\n%d\n%d\n%d\n", s1, s2, s3);

    fclose (f);
    fclose (f2);

    sprintf (s, "input2.%d", getpid ());
    f = open_file_or_exit(s, "w");
    sprintf (s, "input4.%d", getpid ());
    f2 = open_file_or_exit(s, "w");

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

}

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

int check_input_files(const char *filename)
{
    /* cluster-genes reads xxx.cut and xxx.cut.sstats */
    FILE* tmp;
    char  input_filename[255];
    int   success;

    success = 1;
    
    sprintf(input_filename, "%s.sstats", filename);    
    tmp = fopen(input_filename,"r");
    if (!tmp)
    {
        fprintf(stderr, "ERROR: %s does not exist\n", input_filename);
        fflush(stderr);
        success = 0;
    } else {
        fclose(tmp);
    }

    strcpy(input_filename, filename);
    tmp = fopen(input_filename,"r");
    if (!tmp)
    {
        fprintf(stderr, "ERROR: %s does not exist\n", input_filename);
        fflush(stderr);
        success = 0;
    } else {
        fclose(tmp);
    }

    return success;
}

void removeTempFiles()
{

    /* clean up */
    /*
    char s[1024];

    sprintf (s, "input1.%d", getpid ());
    unlink (s);
    sprintf (s, "input2.%d", getpid ());
    unlink (s);
    sprintf (s, "input3.%d", getpid ());
    unlink (s);
    sprintf (s, "input4.%d", getpid ());
    unlink (s);
    sprintf (s, "%s.emd", datafile);
    unlink (s);
    sprintf (s, "out1.%d", getpid ());
    unlink (s);
    sprintf (s, "out2.%d", getpid ());
    unlink (s);
    sprintf (s, "%sm", datafile);
    unlink (s);
    sprintf (s, "tmp.%d", getpid ());
    unlink (s);
    */

}

int resultsNotConverge(char* filename)
{
    FILE* fp;
    char  s[LINE_LENGTH];
    int   iFound = 0;
    
    fp = open_file_or_exit(filename, "r");

    while (!feof(fp))
    {
        fgets(s, LINE_LENGTH, fp);
        if (strstr(s, "did not converge"))
        {
            iFound = 1;
            break;
        }
    }
    
    fclose(fp);

    return iFound;
}

int copyEmmixTResults(char *data_file, int nGroup)
{
    char src[255];
    char dst[255];
    int i;

    for (i=1;i<=2;i++)
    {
        sprintf(src, "out%d.%d", i, getpid());
        sprintf(dst, "%s.emmixt.%d.%d", data_file, i, nGroup);
        file_copy(src, dst); 
    }

    /* also copy emmix-t inputs */
    sprintf(src, "tmp.%d", getpid());
    sprintf(dst, "%s.emmixt.in.%d", data_file, nGroup);
    file_copy(src, dst);


    /* if the result file did not converge, then warn user */
    sprintf(src, "%s.emmixt.%d.%d", data_file, 2, nGroup); 
    if (resultsNotConverge(src))
    {
        fprintf(stderr, "ANNOUNCE: %s did not converge\n", src);
        fflush (stderr);
    }

    return 1;

}

int main (int argc, char **argv)
{
  int row;			/* 2nd argument */
  int col;			/* 3rd argument */
  int g;			/* 4th argument */
  int clus;         /* number of groups */

  int r, k, s1, s2, s3;		/* random starts, k-means starts, seeds 1 to 3 */
  int r1, k1, nu;             /* sph random starts, sph k-means starts       */

  int reord;
  float b2;


  if (argc != 4)
    {
      printf
	("usage: cluster-genes microarray-filename number-of-rows number-of-columns\n");
      exit (1);
    }

  /* datafile is global and is used by removeTempFiles() */
  strcpy(datafile, argv[1]);
  if (!check_input_files(datafile))
  {
    exit (1);
  };
  
  col = atoi (argv[3]);
  row = atoi (argv[2]);

  if (col > MAX_TISSUES)
  {
    printf("Recompile with MAX_TISSUES > %d\n", col);
    exit (1);
  }
  
  g = 1;
  /* rm -f error.$1.$$ log.$1.$$ smaller.$1.$$ input1.$$ */
  /* sprintf(s,"error.%d",getpid());
     unlink(s);
     sprintf(s,"log.%d",getpid());
     unlink(s);
     sprintf(s,"smaller.%d",getpid());
     unlink(s); */

  atexit(removeTempFiles); 
  
  user_input(&r,&k,&r1,&k1,&clus,&b2,&s1,&s2,&s3,&reord,&nu);

  call_emd (datafile,col,g,clus,row,r,k,r1,k1,b2, 
            nu,s2,s2,s3,reord);

  if (sm < b2 )
  {
      call_emd (datafile,col,g+1,clus,row,r,k,r1,k1,b2, 
                nu,s2,s2,s3,reord);
  }

  return 0;
}

void
call_emd (char *data_file, int col, int g, int clus, int row, 
          int r, int k, int r1, int k1, float b2, 
          int nu, int s1, int s2, int s3, 
          int reord)
{
  char s[LINE_LENGTH];
  char t[LINE_LENGTH];
  int grea[MAX_GROUPS][MAX_TISSUES];	/* up to MAX_GROUPS group means, up to MAX_TISSUES samples per mean */

  struct tab gl[MAX_GROUPS];
  char   emd_input[255];        /* name of file to pipe stdin to Emmix-spher */
  int    gmap[MAX_GROUPS];		/* mapping between group numbers */
  FILE   *f, *f2, *f4;
  int    i, j, sm, er;
  float  tl;

  /* initializes all the temp. filenames */  
  initFilenames(data_file, clus);
  unlinkOutputFiles(clus);

  /* prepare input file for emmix-spher */
  sprintf (emd_input, "emd-input1.%d", getpid());
  f = open_file_or_exit(emd_input, "w");
  fprintf (f, "2\n");
  fprintf (f, "%s\n", data_file);	    /* input file  */
  fprintf (f, "%s.emd\n", data_file);   /* output file */
  fprintf (f, "%sm\n", data_file);	    /* matlab file */
  fprintf (f, "%d\n", row);
  fprintf (f, "%d\n", col);
  fprintf (f, "0\n");
  fprintf (f, "%d\n", clus);
  fprintf (f, "5\n3\n");	/* option 5, eq, option 3, random starts */
  fprintf (f, "%d\n100\n%d\n", r1, k1);	/* rand starts, 100%, k-means */
  if (nu)
  {
    fprintf (f, "y\n");
    fprintf (f, "11\n");
    fprintf (f, "%d\n", nu);
    fprintf (f, "0\n");
  } else {
    fprintf (f, "n\n");
  }
  fprintf (f, "%d\n%d\n%d\n", s1, s2, s3);	/* rand seeds */
  fclose (f);
  
  fprintf(stderr, "ANNOUNCE: EMMIX-spher started\n");
  fflush(stderr);

  sprintf (s, "%s < %s > %s", EMMIXS_EXE, emd_input, "EMMIX-spher.out");
  /* XXX CGT commented out to speed up testing */   
  system (s);

  fprintf(stderr, "ANNOUNCE: EMMIX-spher completed\n");
  fflush(stderr);

  sprintf (s, "%s.emd", data_file);
  f = open_file_or_exit(s, "r");

  /* go down to "Estimated mean...for each component" to get group means */
  while ((!strstr (s, "Estimated mean") || 
          !strstr (s, "for each component")) && !feof (f))
    fgets (s, LINE_LENGTH, f);
  if (feof (f))
    {
      perror("ERROR: Failed to find start\n");
      exit (1);
    }
    
  /* write out the group means to the groupmean file */
  /* note this file gets overwritten later           */
  f2 = open_file_or_exit (groupmean_filename, "w");
  for (j = 0; j < clus; j++)
    {
      for (i = 0; i < col; i++)
    	{
	      fscanf (f, "%s", s);
    	  fprintf (f2, "%s ", s);
    	}
      fprintf (f2, "\n");
    }
  fclose (f2);
  fclose (f);

  /* populate gl[].s with contents of groupmean file */
  f = open_file_or_exit (groupmean_filename, "r");
  for (i = 0; i < clus; i++)
    fgets (gl[i].s, LINE_LENGTH, f);
  fclose (f);

  /* calculate -2 log \lambda for group means */
  f4 = open_file_or_exit (groupmean_filename, "r");
  for (i = 0; i < clus; i++)
    {
      fgets (s, LINE_LENGTH, f4);
      if ( s[strlen(s)-1] != ENDL)
      {
        fprintf (stderr, "ERROR: The text in the groupmean file exceeds"
        " the maximum length of %d\n.", LINE_LENGTH);
        fflush (stderr);
        exit(1);
      }
      strip_cr (s, t);

      makeInput1(r,k,r1,k1,clus,b2,s1,s2,s3,reord,nu,col);
      /* makeInput2(); */
      
      fprintf(stderr, "ANNOUNCE: EMMIX-T.exe step %d of %d\n", i+1, clus);
      fflush(stderr);
      call_emmix (data_file, t, &tl, &sm, &er, col, g);

      if (sm < b2)
      {
	      call_emmix(data_file, t, &tl, &sm, &er, col, g+1);
      }
      
      /* make a permanent record of out1 and out2 */
      copyEmmixTResults(data_file, i+1);

      
      /* extract grouping for possible rearrangement */
      sprintf (s, "out2.%d", getpid ());
      f2 = open_file_or_exit(s, "r");
      while (!strstr (s, "Implied") && !feof (f2))
        fgets (s, LINE_LENGTH, f2);
      if (!feof(f2))
      {
          for (j = 0; j < col; j++)
            fscanf (f2, "%d", &(grea[i][j]));
      } else {
          fprintf(stderr, "ANNOUNCE: Could not determine "
                          "the implied grouping\n");
          fflush(stderr);
      }
      fclose (f2);

      gl[i].nu = i;
      gl[i].tl = tl;

      /*fgets (s, LINE_LENGTH, f); */
    }

  /* CGT group4 is failing regression test */
  printf ("CGT watch group4 carefully...\n");
  for (i=0;i<clus;i++)
  {
    for (j=0;j<col;j++)
        printf ("%d ", grea[i][j]);
    printf("\n");
  }
        
  fclose (f4);
  
  /* sort gl[] by gl[].tl */
  qsort ((void *) gl, clus, sizeof (struct tab), tabcompare);

  for (i = 0; i < clus; i++)
    gmap[gl[i].nu] = i;

  /* open .cut file */
  /*sprintf (s, "%s.cut", data_file); */
  f4 = open_file_or_exit (data_file, "r");

  /* open matlab file */
  sprintf (s, "%sm", data_file);
  f2 = open_file_or_exit(s, "r");

  for (i = 0; i < clus; i++)
    gl[i].count = 0;
  /* go through each a line at a time, writing lines to _group%d files */
  for (i = 0; i < row; i++)
    {
      FILE *f3, *f5, *f6;
      float f;
      int j,tmp,tmp2,x;
      fscanf (f2, "%d", &j);
      f5 = open_file_or_exit (list_filename[gmap[j-1]+1], "a");

      sprintf (s, "%s.sstats", data_file);
      f6 = open_file_or_exit(s, "r");
      /* scan through i lines */
      for (x = 0; x < i; x++)
          /*fscanf(f6,"%s",s);*/ fgets(s, LINE_LENGTH, f6);
      fscanf(f6,"%d %f %d",&tmp,&f,&tmp2);
      fclose(f6);

      fprintf (f5, "%d %f\n", i + 1, f); /* write number and -2ll */
      
      fclose (f5);
      gl[gmap[j - 1]].count++;
      f3 = open_file_or_exit(group_filename[gmap[j-1]+1], "a");
      fgets (s, LINE_LENGTH, f4);
      if (!reord) 
      {
      
    	fputs (s, f3);
        
      } else {
      
          /* reorder routine */
          char u[MAX_TISSUES][MAX_ENTRY_LENGTH],
               v[MAX_TISSUES][MAX_ENTRY_LENGTH];
          int k, m, l = 0;
          strcpy (u[0], strtok (s, " "));
          for (m = 1; m < col; m++)
            strcpy (u[m], strtok (0, " "));
          /* remove '\n' from u[col-1] */
          u[col - 1][strlen (u[col - 1]) - 1] = 0;
          for (k = 1; k <= g + 1; k++)
            for (m = 0; m < col; m++)
              if (grea[j - 1][m] == k)
                 strcpy (v[l++], u[m]);
          for (m = 0; m < col; m++)
            fprintf (f3, "%s ", v[m]);
          fprintf (f3, "\n");
      }
      fclose (f3);

    } /* for (i<row) */
    
  fclose (f4);
  fclose (f2);

  /* Write out gstats file */
  sprintf (s, "%s.gstats", data_file);
  f2 = open_file_or_exit (s, "w");
  for (i = 0; i < clus; i++)
    fprintf (f2, "%d %d %f\n", i+1, gl[i].count, gl[i].tl);
  fclose (f2);

  /* Write out groupmean file again, but this time
  ** it has been sorted 
  */
  f2 = open_file_or_exit (groupmean_filename, "w");
  for (i = 0; i < clus; i++)
    fprintf (f2, "%s", gl[i].s);
  fclose (f2);

  /* Rename the temp files to permanent names */
  makeTempFilesPermanent(data_file, clus);

}

int getFinalLogLikelihood(char* filename, float *result)
{

    FILE *fd;
    int  intSuccess;
    char s[LINE_LENGTH];

    fd = open_file_or_exit (filename, "r");
    fgets (s, LINE_LENGTH, fd);
    while (!strstr (s, "Log-Likelihood") && !feof (fd))
        fgets (s, LINE_LENGTH, fd);
    if (!feof(fd))  
    {
        intSuccess = sscanf (s, " Final Log-Likelihood is %f", result);
    }
    fclose (fd);

    return intSuccess;
}

void
call_emmix (char *data_file, char *t, float *tl, int *sm, int *er, int col, int g)
{
  FILE *f2;
  float fl1, fl2;   /* Final Log-Likelihood */
  int  group2[MAX_TISSUES];
  int  comp[MAX_TISSUES], a1, a2, i;
  char s[LINE_LENGTH];
  char filename[LINE_LENGTH];

  sprintf (s, "tmp.%d", getpid ());
  f2 = open_file_or_exit (s, "w");
  fputs (t, f2);
  fprintf (f2, "\n");
  fclose (f2);

  /* Only run for 1 vs 2 ie. g=1 */ 
  sprintf (s, "%s < input1.%d > %s", EMMIXT_EXE, getpid (), "EMMIX-T.out");
  system (s);
  sprintf (s, "%s < input2.%d > %s", EMMIXT_EXE, getpid (), "EMMIX-T.out");
  system (s);

  sprintf (filename, "out1.%d", getpid());
  a1 = getFinalLogLikelihood(filename, &fl1);

  sprintf (filename, "out2.%d", getpid());
  a2 = getFinalLogLikelihood(filename, &fl2);

  if (a1 && a2)
  {
    *tl = 2.0 * (fl2 - fl1);
  }
  else
    {

      *tl = -100;
      *sm = -100;
      *er = -100;
      
      return;
    }

  /* calculate error rate */
  /*sprintf (s, "%s-group", data_file);
     f2 = fopen (s, "r");
     if (!f2) perror(s);
     for (i = 0; i < col; i++)
     fscanf (f2, "%d", &group1[i]);
     fclose (f2); */

  sprintf (s, "out2.%d", getpid ());
  f2 = open_file_or_exit(s, "r");
  
  /* Search for the word "Implied" */
  while (!strstr (s, "Implied") && !feof (f2))
    fgets (s, LINE_LENGTH, f2);
  if (!feof(f2))
  {
      for (i = 0; i < col; i++)
      {
        fscanf (f2, "%d", &group2[i]);
      }
  }
  fclose (f2);

  *er = 0;
  /* for (i = 0; i < col; i++)
     if (group1[i] != group2[i])
     (*er)++;
     if (col - *er < *er)
     *er = col - *er; */

  memset (comp, 0, sizeof (comp));
  for (i = 0; i < col; i++)
  {
    comp[group2[i]]++;
  }

  /* Retrieve the smallest value of comp[]
  ** into sm 
  */
  *sm = 9999;
  for (i = 1; i <= g + 1; i++)
    if (comp[i] < *sm)
      *sm = comp[i];
}

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
