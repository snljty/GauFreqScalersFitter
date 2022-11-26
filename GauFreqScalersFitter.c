# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <stdbool.h>

/*
 * This program uses FtotalFrequencies/10 and ZPVE15/10 in JCTC,6,2872(2010)       *
 * to fit fundamental frequency scaler and zero-point energy scaler. *
 * The 15 molecules used are: HF, H2, N2, F2, CO, OH, Cl2, CO2, H2O, *
 * N2O, HCN, C2H2, H2CO, NH3, CH4.                                   *
 * The result will be saved to scalers_result.txt.                   *
 */

typedef struct {double slope, r_squared;} FitResult;

# define totalMolecules 15
# define totalFrequencies 38

# define maxNameLength 5
# define maxExtNameLength 5

# define freqThreshhold 3.

# define generateInputError -1
# define calculateError -2

const char * fileNames[] = {
  "HF" , "H2"  , "N2"  , "F2" , "CO" ,
  "OH" , "Cl2" , "CO2" , "H2O", "N2O",
  "HCN", "C2H2", "H2CO", "NH3", "CH4" };
/*  ZPE in the order listed above  */
const double ZPE_expe[] = {
   5.86353,
   6.2310,
   3.3618,
   1.3021,
   3.0929144,
   5.2915,
   0.7983,
   7.3,
  13.26,
   6.770,
  10.0,
  16.49,
  16.1,
  21.200,
  27.710 };
double ZPE_comp[totalMolecules] = {0.};
/*  Fundamental Frequency, as above order, but each individuals  *
 *  have these amounts of frequencies to be considered:          *
 *  1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 5, 6, 4, 4                  */
const double fundFreq_expe[] = {
  3959,
  4159,
  2330,
   894,
  2143,
  3568,
   554,
   667, 1333, 2349,
  1595, 3657, 3756,
   589, 1285, 2224,
   712, 2089, 3312,
   612,  730, 1974, 3289, 3374,
  1167, 1249, 1500, 1746, 2782, 2843,
   950, 1627, 3337, 3444,
  1306, 1534, 2917, 3019 };
double fundFreq_comp[totalFrequencies] = {0.};
unsigned int freqDataPos = 0;
char gauExe[] = "g16";
char const templateName[] = "template.gjf";
char const CMGDetector[] = "[CHARGE_MULTIPLICITY_AND_GEOMETRY]";
char const nameField[] = "[NAME]";

void Pauser(const char *prompt);
void PrintChargeMultiplicityAndGeometry(FILE *oflp, unsigned int i);
int PrepareInput(const char *level, unsigned int i);
int CalcFile(unsigned int i);
int ReadFreq(unsigned int i);
FitResult NoInterceptLinearFit(const double *x, const double *y, unsigned int n);

int main(int argc, const char *argv[])
{
  int iarg = 0;
  char buf[BUFSIZ + 1] = "";
  unsigned int i = 0;
  char level[BUFSIZ + 1] = "";
  FitResult fitFundFreq = {1., 1.};
  FitResult fitZPE = {1., 1.};
  FILE *iflp = NULL, *oflp = NULL;
  time_t startTime = 0, endTime = 0;
  bool isCMGFound = false; /* only will be used if file templateName exists */

  /* quick check help */
  for (iarg = 1; iarg < argc; ++ iarg)
  {
    if (! (strcmp(argv[iarg], "--help") && strcmp(argv[iarg], "-h") && strcmp(argv[iarg], "/?")))
    {
      fprintf(stderr, "Usage: %s [LEVEL]\n", argv[0]);
      fprintf(stderr, "\n");
      fprintf(stderr, "If file \"%s\" exists, then the level input will be ignored, \n", templateName);
      fprintf(stderr, "and level from \"%s\" will be used. In this file, \n", templateName);
      fprintf(stderr, "the CHARGE, MULTIPLICILY and GEOMETRY should be replaced by \"%s\".\n", CMGDetector);
      fprintf(stderr, "If \"%s\" exists in \"%s\", it will be replaced with the filename without the suffix.\n", \
        nameField, templateName);
      fprintf(stderr, "If file \"%s\" does not exist, the level from the first command argument \n", templateName);
      fprintf(stderr, "will be used. If it is blank, then the interactive mode will hint to input the level.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "Exiting now.\n");
      exit(EXIT_SUCCESS);
    }
  }

  /*  get level  */
  iflp = fopen(templateName, "rt");
  if (iflp)
  {
    while (fgets(buf, BUFSIZ, iflp))
    {
      if (! strncmp(buf, CMGDetector, strlen(CMGDetector)))
      {
        isCMGFound = true;
        break;
      }
    }
    fclose(iflp);
    iflp = NULL;
    if (! isCMGFound)
    {
      fprintf(stderr, "Error! Albeit file \"%s\" found, cannot locate \"%s\" in this file.\n", \
        templateName, CMGDetector);
      fprintf(stderr, "Please either writes \"%s\" in \"%s\" or remove this file.\n", \
        CMGDetector, templateName);
      fprintf(stderr, "Exiting abnormally.\n");
      exit(EXIT_FAILURE);
    }
    strcpy(level, "");
    printf("Level from file \"%s\" will be used.\n", templateName);
  }
  else
  {
    if (argc == 1)
    {
      printf("Tests %d small molecules and generates frequency scalers\n", \
        (int)(sizeof(ZPE_expe) / sizeof(double)));
      puts("for fundamental frequency and zero-point energy using Gaussian.\n");
      puts("Before using, you need to set system environment variable");
      puts("GAUSS_EXEDIR to the directory containing your gaussian executable.");
      puts("Also, this path should be added to environment variable PATH.\n");
      puts("Input the level and related keywords, e.g. ");
      puts("M062X/def2SVP EmpiricalDispersion=GD3 Integral=UltrafineGrid");
      puts("If you press ENTER directly, B3LYP/6-31G* will be used:");
      fgets(level, BUFSIZ, stdin);
      if (level[strlen(level) - 1] == '\n')
      {
        level[strlen(level) - 1] = '\0';
      }
      if (! strcmp(level, ""))
      {
        strcpy(level, "B3LYP/6-31G*");
      }
      puts("");
    }
    else
    {
      strcpy(level, argv[1]);
      printf("Level \"%s\" from the first command argument will be used.\n", level);
    }
  }

  /*  invoke Gaussian  */
  time(& startTime);
  for (i = 0; i < totalMolecules; ++ i)
  {
    /*  prepare input files  */
    if (PrepareInput(level, i))
    {
      exit(generateInputError);
    }
    /*  call Gaussian  */
    printf("Calculating %s...\n", fileNames[i]);
    if (CalcFile(i))
    {
      exit(calculateError);
    }
    ReadFreq(i);
  }
  time(& endTime);

  /*  Print Values  */
  if (argc == 1)
  {
    for (i = 0; i < totalFrequencies; ++ i)
    {
      printf("Calculated   frequency %2d: %9.4lf\n", i + 1, fundFreq_comp[i]);
      printf("Experimental frequency %2d: %9.4lf\n", i + 1, fundFreq_expe[i]);
    }
    puts("");
    for (i = 0; i < totalMolecules; ++ i)
    {
      printf("Calculated   ZPE %2d: %8.5lf\n", i + 1, ZPE_comp[i]);
      printf("Experimental ZPE %2d: %8.5lf\n", i + 1, ZPE_expe[i]);
    }
  }

  /*  fit values  */
  fitFundFreq = NoInterceptLinearFit(fundFreq_comp, fundFreq_expe, totalFrequencies);
  fitZPE = NoInterceptLinearFit(ZPE_comp, ZPE_expe, totalMolecules);

  /*  output result  */
  puts("");
  oflp = fopen("scalers_result.txt", "wt");
  if (* level)
  {
    fprintf(oflp, "level: %s\n", level);
  }
  else
  {
    fprintf(oflp, "level: see file \"%s\"\n", templateName);
  }
  puts("");
  printf("Scaler for fundamental frequency at this level: %6.4lf\n", fitFundFreq.slope);
  fprintf(oflp, "Scaler for fundamental frequency at this level: %6.4lf\n", fitFundFreq.slope);
  printf("Squared statistical significant linear coefficient of correlation: %6.4lf\n", fitFundFreq.r_squared);
  fprintf(oflp, "Squared statistical significant linear coefficient of correlation: %6.4lf\n", fitFundFreq.r_squared);
  printf("\n");
  fprintf(oflp, "\n");
  printf("Scaler for ZPE at this level: %6.4lf\n", fitZPE.slope);
  fprintf(oflp, "Scaler for ZPE at this level: %6.4lf\n", fitZPE.slope);
  printf("Squared statistical significant linear coefficient of correlation: %6.4lf\n", fitZPE.r_squared);
  fprintf(oflp, "Squared statistical significant linear coefficient of correlation: %6.4lf\n", fitZPE.r_squared);
  fclose(oflp);
  oflp = NULL;

  /*  final work  */
  printf("Time Elapsed: %d s.\n", (int)(endTime) - (int)startTime);
  puts("\nDone! The result has been saved to \"scalers_result.txt\".");
  if (argc == 1)
  {
    puts("Press ENTER directly will remove temporary files and exit.");
    puts("If anything is input, temporary files will be kept.");
    fgets(buf, BUFSIZ, stdin);
    if (! strcmp(buf, "\n"))
    {
      for (i = 0; i < totalMolecules; ++ i)
      {
        printf("Removing %s...\n", fileNames[i]);
        strncpy(buf, fileNames[i], maxNameLength);
        strcat(buf, ".gjf");
        remove(buf);
        strncpy(buf, fileNames[i], maxNameLength);
        strcat(buf, ".out");
        remove(buf);
        if (* level)
        {
          strncpy(buf, fileNames[i], maxNameLength);
          strcat(buf, ".chk");
          remove(buf);
        }
      }
      /* remove("fort.6"); */
      /* remove("fort.7"); */
    }
  }
  else
  {
    for (i = 0; i < totalMolecules; ++ i)
    {
      strncpy(buf, fileNames[i], maxNameLength);
      strcat(buf, ".gjf");
      remove(buf);
      strncpy(buf, fileNames[i], maxNameLength);
      strcat(buf, ".out");
      remove(buf);
    }
  }

  return 0;
}


void Pauser(const char *prompt)
{
  char c = '\0';
  puts(prompt);
  while ((c = getchar()) != '\n' && c != EOF)
  {
    ;
  }

  return;
}

void PrintChargeMultiplicityAndGeometry(FILE *oflp, unsigned int i)
{
  switch (i)
  {
  case 0:
    /*  HF  */
    fputs("0 1\n", oflp);
    fputs(" F                  0.00000000    0.00000000    0.09338300\n", oflp);
    fputs(" H                  0.00000000    0.00000000   -0.84044600\n", oflp);
    break;
  case 1:
    /*  H2  */
    fputs("0 1\n", oflp);
    fputs(" H                  0.00000000    0.00000000    0.37139400\n", oflp);
    fputs(" H                  0.00000000    0.00000000   -0.37139400\n", oflp);
    break;
  case 2:
    /*  N2  */
    fputs("0 1\n", oflp);
    fputs(" N                  0.00000000    0.00000000    0.55274900\n", oflp);
    fputs(" N                  0.00000000    0.00000000   -0.55274900\n", oflp);
    break;
  case 3:
    /*  F2  */
    fputs("0 1\n", oflp);
    fputs(" F                  0.00000000    0.00000000    0.70140600\n", oflp);
    fputs(" F                  0.00000000    0.00000000   -0.70140600\n", oflp);
    break;
  case 4:
    /*  CO  */
    fputs("0 1\n", oflp);
    fputs(" O                  0.00000000    0.00000000    0.48768200\n", oflp);
    fputs(" C                  0.00000000    0.00000000   -0.65024300\n", oflp);
    break;
  case 5:
    /*  OH  */
    fputs("0 2\n", oflp);
    fputs(" O                  0.00000000    0.00000000    0.10921300\n", oflp);
    fputs(" H                  0.00000000    0.00000000   -0.87370400\n", oflp);
    break;
  case 6:
    /*  Cl2  */
    fputs("0 1\n", oflp);
    fputs(" Cl                 0.00000000    0.00000000    1.02086800\n", oflp);
    fputs(" Cl                 0.00000000    0.00000000   -1.02086800\n", oflp);
    break;
  case 7:
    /*  CO2  */
    fputs("0 1\n", oflp);
    fputs(" C                  0.00000000    0.00000000    0.00000000\n", oflp);
    fputs(" O                  0.00000000    0.00000000    1.16914800\n", oflp);
    fputs(" O                  0.00000000    0.00000000   -1.16914800\n", oflp);
    break;
  case 8:
    /*  H2O  */
    fputs("0 1\n", oflp);
    fputs(" O                  0.00000000    0.00000000    0.11971800\n", oflp);
    fputs(" H                  0.00000000    0.76155000   -0.47887400\n", oflp);
    fputs(" H                  0.00000000   -0.76155000   -0.47887400\n", oflp);
    break;
  case 9:
    /*  N2O  */
    fputs("0 1\n", oflp);
    fputs(" N                  0.00000000    0.00000000   -0.07271700\n", oflp);
    fputs(" N                  0.00000000    0.00000000   -1.20710100\n", oflp);
    fputs(" O                  0.00000000    0.00000000    1.11984000\n", oflp);
    break;
  case 10:
    /*  HCN  */
    fputs("0 1\n", oflp);
    fputs(" C                  0.00000000    0.00000000   -0.50200400\n", oflp);
    fputs(" H                  0.00000000    0.00000000   -1.57213400\n", oflp);
    fputs(" N                  0.00000000    0.00000000    0.65488000\n", oflp);
    break;
  case 11:
    /*  C2H2  */
    fputs("0 1\n", oflp);
    fputs(" C                  0.00000000    0.00000000    0.60248600\n", oflp);
    fputs(" H                  0.00000000    0.00000000    1.66912800\n", oflp);
    fputs(" C                  0.00000000    0.00000000   -0.60248600\n", oflp);
    fputs(" H                  0.00000000    0.00000000   -1.66912800\n", oflp);
    break;
  case 12:
    /*  H2CO  */
    fputs("0 1\n", oflp);
    fputs(" C                  0.00000000    0.00000000   -0.52895600\n", oflp);
    fputs(" H                  0.00000000    0.93789200   -1.12366200\n", oflp);
    fputs(" H                  0.00000000   -0.93789200   -1.12366200\n", oflp);
    fputs(" O                  0.00000000    0.00000000    0.67763300\n", oflp);
    break;
  case 13:
    /*  NH3  */
    fputs("0 1\n", oflp);
    fputs(" N                  0.00000000    0.00000000    0.11935700\n", oflp);
    fputs(" H                  0.00000000    0.93857800   -0.27849900\n", oflp);
    fputs(" H                 -0.81283200   -0.46928900   -0.27849900\n", oflp);
    fputs(" H                  0.81283200   -0.46928900   -0.27849900\n", oflp);
    break;
  case 14:
    /*  CH4  */
    fputs("0 1\n", oflp);
    fputs(" C                  0.00000000    0.00000000    0.00000000\n", oflp);
    fputs(" H                  0.63133900    0.63133900    0.63133900\n", oflp);
    fputs(" H                 -0.63133900   -0.63133900    0.63133900\n", oflp);
    fputs(" H                 -0.63133900    0.63133900   -0.63133900\n", oflp);
    fputs(" H                  0.63133900   -0.63133900   -0.63133900\n", oflp);
    break;
  default:
    break;
  }
  return;
}

int PrepareInput(char const *level, unsigned int i)
{
  char fn[BUFSIZ + 1] = "";
  FILE *iflp = NULL, *oflp = NULL;
  char buf[BUFSIZ + 1] = "";
  char *tok = NULL;

  strncpy(fn, fileNames[i], maxNameLength);
  strncat(fn, ".gjf", maxExtNameLength);
  oflp = fopen(fn, "wt");
  if (! oflp)
  {
    Pauser("Error! Cannot create Gaussian input file!");
    return generateInputError;
  }
  printf("Progress %2d/%2d...\n", i + 1, totalMolecules);
  printf("Generating input file for %s...\n", fileNames[i]);
  iflp = fopen(templateName, "rt");
  if (iflp)
  {
    while (fgets(buf, BUFSIZ, iflp))
    {
      if (strncmp(buf, CMGDetector, strlen(CMGDetector)))
      {
        /* todo */
        tok = buf;
        while (* tok != '\0')
        {
          if (strncmp(tok, nameField, strlen(nameField)))
          {
            putc(* tok, oflp);
            ++ tok;
          }
          else
          {
            fputs(fileNames[i], oflp);
            tok += strlen(nameField);
          }
        }
      }
      else
      {
        PrintChargeMultiplicityAndGeometry(oflp, i);
      }
    }
    fclose(iflp);
    iflp = NULL;
  }
  else
  {
    fprintf(oflp, "%%chk=%s.chk\n", fileNames[i]);
    fprintf(oflp, "# Opt %s\n", level);
    fprintf(oflp, "\n");
    fprintf(oflp, "%s\n", fileNames[i]);
    fprintf(oflp, "\n");
    PrintChargeMultiplicityAndGeometry(oflp, i);
    fputs("\n", oflp);
    fputs("--link1--\n", oflp);
    fprintf(oflp, "%%chk=%s.chk\n", fileNames[i]);
    fprintf(oflp, "# Freq %s Geom=AllCheck Guess=Read\n", level);
    fputs("\n", oflp);
  }
  fclose(oflp);
  oflp = NULL;

  return 0;
}

int CalcFile(unsigned int i)
{
  char cmnd[BUFSIZ + 1] = "";

  for (;;)
  {
    strncpy(cmnd, gauExe, 8);
    strcat(cmnd, " ");
    strncat(cmnd, fileNames[i], maxNameLength);
    strcat(cmnd, ".gjf ");
    strncat(cmnd, fileNames[i], maxNameLength);
    strcat(cmnd, ".out");
    if (system(cmnd))
    {
      if (! i && ! strcmp(gauExe, "g16"))
      {
        /*  test whether to use g16 or g09  */
        puts("Calling g16 failed, try calling g09...this may NOT be an error.");
        strcpy(gauExe, "g09");
        continue;
      }
      else
      {
        printf("Error! Calculation of %s failed. Check your PATH and GAUSS_EXEDIR.", fileNames[i]);
        Pauser("Press ENTER to exit abnormally...");
        return calculateError;
      }
    }
    else
    {
      break;
    }
  }

  return 0;
}

int ReadFreq(unsigned int i)
{
  FILE *iflp = NULL;
  char flname[maxNameLength + maxExtNameLength] = "";
  char line[BUFSIZ + 1] = "";
  unsigned int t = 0;
  double tmpVal = 0.;
  unsigned int tmpValAmount = 0;
  char *readPos = NULL;

  strncpy(flname, fileNames[i], maxNameLength);
  strcat(flname, ".out");
  iflp = fopen(flname, "rt");
  /*  read ZPE  */
  while (fgets(line, BUFSIZ, iflp))
  {
    if (strstr(line, "Zero-point vibrational energy"))
    {
      fgets(line, BUFSIZ, iflp);
      sscanf(line, "%lf", & ZPE_comp[i]);
      break;
    }
  }
  fseek(iflp, 0, SEEK_SET);
  /*  read fundamental frequency  */
  while (fgets(line, BUFSIZ, iflp))
  {
    if (strstr(line, "and normal coordinates:"))
    {
      break;
    }
  }
  while (fgets(line, BUFSIZ, iflp))
  {
    if (strstr(line, "Thermochemistry"))
    {
      break;
    }
    if (readPos = strstr(line, "Frequencies"))
    {
      readPos += strlen("Frequencies --");
      while (sscanf(readPos, "%lf", & tmpVal))
      {
        while (* readPos == ' ')
        {
          ++ readPos;
        }
        readPos = strpbrk(readPos, " \n");
        if (tmpValAmount && fabs(fundFreq_comp[freqDataPos] / tmpValAmount - tmpVal) >= freqThreshhold)
        {
          /*  save old value  */
          fundFreq_comp[freqDataPos] /= tmpValAmount;
          tmpValAmount = 0;
          ++ freqDataPos;
        }
        /*  add this number as a fundFreq_comp[freqDataPos]  */
        fundFreq_comp[freqDataPos] += tmpVal;
        ++ tmpValAmount;
        if (* readPos == '\n')
        {
          break;
        }
      }
    }
  }
  /*  reading this file done, force to save the last value.  */
  fundFreq_comp[freqDataPos] /= tmpValAmount;
  tmpValAmount = 0;
  ++ freqDataPos;
  fclose(iflp);
  iflp = NULL;

  return 0;
}

FitResult NoInterceptLinearFit(const double *x, const double *y, unsigned int n)
{
  FitResult ans = {1., 1.};
  double Ex = 0., Ey = 0., Dx =0., Dy = 0.;
  double Exx = 0., Eyy = 0., Exy = 0., Cov = 0.;
  double t = 0.;
  unsigned int i = 0;

  /*  Ex = sum(x)/n                              *
   *  Ey = sum(y)/n                              *
   *  Exx = sum(x^2)/n                           *
   *  Eyy = sum(y^2)/n                           *
   *  Exy = sum(x*y)/n                           *
   *  Dx = Exx-(Ex)^2                            *
   *  Dy = Eyy-(Ey)^2                            *
   *  Cov(x, y) = Exy-Ex*Ey                      *
   *  slope = sum(x*y)/sum(x^2) = Exy / Exx      *
   *  r^2(slope * x + intercept, y) = r^2(x. y)  *
   *  = Cov(x, y)^2/(Dx*Dy)                      */

  /*  Average of x, y, x*y, x^2 and y^2  */
  for (i = 0; i < n; ++ i)
  {
    Ex += x[i];
    Ey += y[i];
    Exy += x[i] * y[i];
    Exx += x[i] * x[i];
    Eyy += y[i] * y[i];
  }
  Ex /= n;
  Ey /= n;
  Exy /= n;
  Exx /= n;
  Eyy /= n;

  /*  deviation var of x and y  */
  Dx = Exx - Ex * Ex;
  Dy = Eyy - Ey * Ey;
  Cov = Exy - Ex * Ey;

  /*  result  */
  ans.slope = Exy / Exx;
  ans.r_squared = Cov * Cov / (Dx * Dy);

  return ans;
}

