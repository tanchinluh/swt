/*
 * -------------------------------------------------------------------------
 * dwt1d.c -- 1-D signal decomposition and reconstruction
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------------
 */
#include "swt_common.h"
#include "dwt.h"

wavelet_identity wi[] = {
  {"haar",ORTH, HAAR, 0, haar_analysis_initialize , haar_synthesis_initialize}, 
  {"db1", ORTH, DAUBECHIES, 1, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db2", ORTH, DAUBECHIES, 2, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db3", ORTH, DAUBECHIES, 3, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db4", ORTH, DAUBECHIES, 4, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db5", ORTH, DAUBECHIES, 5, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db6", ORTH, DAUBECHIES, 6, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db7", ORTH, DAUBECHIES, 7, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db8", ORTH, DAUBECHIES, 8, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db9", ORTH, DAUBECHIES, 9, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db10", ORTH, DAUBECHIES, 10, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db11", ORTH, DAUBECHIES, 11, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db12", ORTH, DAUBECHIES, 12, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db13", ORTH, DAUBECHIES, 13, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db14", ORTH, DAUBECHIES, 14, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db15", ORTH, DAUBECHIES, 15, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db16", ORTH, DAUBECHIES, 16, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db17", ORTH, DAUBECHIES, 17, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db18", ORTH, DAUBECHIES, 18, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"db19", ORTH, DAUBECHIES, 19, daubechies_analysis_initialize, daubechies_synthesis_initialize},
  {"db20", ORTH, DAUBECHIES, 20, daubechies_analysis_initialize, daubechies_synthesis_initialize}, 
  {"coif1", ORTH, COIFLETS, 1, coiflets_analysis_initialize, coiflets_synthesis_initialize},
  {"coif2", ORTH, COIFLETS, 2, coiflets_analysis_initialize, coiflets_synthesis_initialize}, 
  {"coif3", ORTH, COIFLETS, 3, coiflets_analysis_initialize, coiflets_synthesis_initialize},
  {"coif4", ORTH, COIFLETS, 4, coiflets_analysis_initialize, coiflets_synthesis_initialize}, 
  {"coif5", ORTH, COIFLETS, 5, coiflets_analysis_initialize, coiflets_synthesis_initialize},
  {"sym2", ORTH, SYMLETS, 2, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym3", ORTH, SYMLETS, 3, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym4", ORTH, SYMLETS, 4, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym5", ORTH, SYMLETS, 5, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym6", ORTH, SYMLETS, 6, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym7", ORTH, SYMLETS, 7, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym8", ORTH, SYMLETS, 8, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym9", ORTH, SYMLETS, 9, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym10", ORTH, SYMLETS, 10, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym11", ORTH, SYMLETS, 11, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym12", ORTH, SYMLETS, 12, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym13", ORTH, SYMLETS, 13, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym14", ORTH, SYMLETS, 14, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym15", ORTH, SYMLETS, 15, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym16", ORTH, SYMLETS, 16, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym17", ORTH, SYMLETS, 17, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym18", ORTH, SYMLETS, 18, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"sym19", ORTH, SYMLETS, 19, symlets_analysis_initialize, symlets_synthesis_initialize},
  {"sym20", ORTH, SYMLETS, 20, symlets_analysis_initialize, symlets_synthesis_initialize}, 
  {"bior1.1", BIORTH, SPLINE_BIORTH, 11, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior1.3", BIORTH,SPLINE_BIORTH, 13, sp_bior_analysis_initialize, sp_bior_synthesis_initialize}, 
  {"bior1.5", BIORTH,SPLINE_BIORTH, 15, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior2.2", BIORTH,SPLINE_BIORTH, 22, sp_bior_analysis_initialize, sp_bior_synthesis_initialize}, 
  {"bior2.4", BIORTH,SPLINE_BIORTH, 24, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior2.6", BIORTH,SPLINE_BIORTH, 26, sp_bior_analysis_initialize, sp_bior_synthesis_initialize}, 
  {"bior2.8", BIORTH,SPLINE_BIORTH, 28, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior3.1", BIORTH,SPLINE_BIORTH, 31, sp_bior_analysis_initialize, sp_bior_synthesis_initialize}, 
  {"bior3.3", BIORTH,SPLINE_BIORTH, 33, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior3.5", BIORTH,SPLINE_BIORTH, 35, sp_bior_analysis_initialize, sp_bior_synthesis_initialize}, 
  {"bior3.7", BIORTH,SPLINE_BIORTH, 37, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior3.9", BIORTH,SPLINE_BIORTH, 39, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior4.4", BIORTH,SPLINE_BIORTH, 44, sp_bior_analysis_initialize, sp_bior_synthesis_initialize}, 
  {"bior5.5", BIORTH,SPLINE_BIORTH, 55, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"bior6.8", BIORTH,SPLINE_BIORTH, 68, sp_bior_analysis_initialize, sp_bior_synthesis_initialize},
  {"rbior1.1", BIORTH,SPLINE_RBIORTH, 11, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior1.3", BIORTH,SPLINE_RBIORTH, 13, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize}, 
  {"rbior1.5", BIORTH,SPLINE_RBIORTH, 15, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior2.2", BIORTH,SPLINE_RBIORTH, 22, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize}, 
  {"rbior2.4", BIORTH,SPLINE_RBIORTH, 24, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior2.6", BIORTH,SPLINE_RBIORTH, 26, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize}, 
  {"rbior2.8", BIORTH,SPLINE_RBIORTH, 28, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior3.1", BIORTH,SPLINE_RBIORTH, 31, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize}, 
  {"rbior3.3", BIORTH,SPLINE_RBIORTH, 33, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior3.5", BIORTH,SPLINE_RBIORTH, 35, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize}, 
  {"rbior3.7", BIORTH,SPLINE_RBIORTH, 37, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior3.9", BIORTH,SPLINE_RBIORTH, 39, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior4.4", BIORTH,SPLINE_RBIORTH, 44, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize}, 
  {"rbior5.5", BIORTH,SPLINE_RBIORTH, 55, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"rbior6.8", BIORTH,SPLINE_RBIORTH, 68, sp_rbior_analysis_initialize, sp_rbior_synthesis_initialize},
  {"beylkin", ORTH, BEYLKIN, 0, beylkin_analysis_initialize, beylkin_synthesis_initialize},
  {"vaidyanathan", ORTH, VAIDYANATHAN, 0, vaidyanathan_analysis_initialize, vaidyanathan_synthesis_initialize},
  {"dmey", ORTH, DMEY, 0, dmey_analysis_initialize, dmey_synthesis_initialize},
  {"bath4.0", ORTH, BATHLETS, 40, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.1", ORTH, BATHLETS, 41, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.2", ORTH, BATHLETS, 42, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.3", ORTH, BATHLETS, 43, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.4", ORTH, BATHLETS, 44, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.5", ORTH, BATHLETS, 45, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.6", ORTH, BATHLETS, 46, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.7", ORTH, BATHLETS, 47, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.8", ORTH, BATHLETS, 48, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.9", ORTH, BATHLETS, 49, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.10", ORTH, BATHLETS, 410, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.11", ORTH, BATHLETS, 411, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.12", ORTH, BATHLETS, 412, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.13", ORTH, BATHLETS, 413, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.14", ORTH, BATHLETS, 414, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath4.15", ORTH, BATHLETS, 415, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.0", ORTH, BATHLETS, 60, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.1", ORTH, BATHLETS, 61, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.2", ORTH, BATHLETS, 62, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.3", ORTH, BATHLETS, 63, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.4", ORTH, BATHLETS, 64, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.5", ORTH, BATHLETS, 65, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.6", ORTH, BATHLETS, 66, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.7", ORTH, BATHLETS, 67, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.8", ORTH, BATHLETS, 68, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.9", ORTH, BATHLETS, 69, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.10", ORTH, BATHLETS, 610, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.11", ORTH, BATHLETS, 611, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.12", ORTH, BATHLETS, 612, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.13", ORTH, BATHLETS, 613, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.14", ORTH, BATHLETS, 614, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"bath6.15", ORTH, BATHLETS, 615, bathlets_analysis_initialize, bathlets_synthesis_initialize},
  {"legd1", ORTH, LEGENDRE, 1, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd2", ORTH, LEGENDRE, 2, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd3", ORTH, LEGENDRE, 3, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd4", ORTH, LEGENDRE, 4, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd5", ORTH, LEGENDRE, 5, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd6", ORTH, LEGENDRE, 6, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd7", ORTH, LEGENDRE, 7, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd8", ORTH, LEGENDRE, 8, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"legd9", ORTH, LEGENDRE, 9, legendre_analysis_initialize, legendre_synthesis_initialize},
  {"fa1", ORTH, FARRAS, 1, farras_analysis_initialize, farras_synthesis_initialize},
  {"fa2", ORTH, FARRAS, 2, farras_analysis_initialize, farras_synthesis_initialize},
  {"ksq1", ORTH, KINGSBURYQ, 1, kingsburyq_analysis_initialize, kingsburyq_synthesis_initialize},
  {"ksq2", ORTH, KINGSBURYQ, 2, kingsburyq_analysis_initialize, kingsburyq_synthesis_initialize}
};

int waveletIdentityNum = sizeof(wi)/sizeof(wavelet_identity);

extend_method dwtMode = SYMH;

/*-------------------------------------------
 * wavelet name parser
 *-----------------------------------------*/

void
filter_clear ()
{
  int count;
  for(count=0;count<30;count++)
    {
      LowDecomFilCoef[count] = 0;
      LowReconFilCoef[count] = 0;
      HiDecomFilCoef[count] = 0;
      HiReconFilCoef[count] = 0;
    }
  return;
}

/*-------------------------------------------
 * Orthfilter Group
 *-----------------------------------------*/

void
orth_filt_group (double *filterIn, int sigInLength,
		 double *filterLowRec, double *filterLowDec,
		 double *filterHiRec, double *filterHiDec)
{
  int count;

  for (count = 0; count < sigInLength; count++)
    filterLowRec[count] = filterIn[count];
  wrev (filterLowRec, sigInLength, filterLowDec, sigInLength);
  qmf_even (filterLowRec, sigInLength, filterHiRec, sigInLength);
  wrev (filterHiRec, sigInLength, filterHiDec, sigInLength);
  return;
}

void 
bior_filt_group (double *f1, int sigInLength1, 
		 double *f2, int sigInLength2, 
		 double *lowDecom, int sigOutLength1, 
		 double *hiDecom, int sigOutLength2,
		 double *lowRecon, int sigOutLength3, 
		 double *hiRecon, int sigOutLength4)
{
  verbatim_copy (f2, sigInLength2, lowRecon, sigOutLength3);
  wrev (f1, sigInLength1, lowDecom, sigOutLength1);
  qmf_even (f1, sigInLength1, hiRecon, sigOutLength4);
  //qmf_odd (f1, sigInLength1, hiRecon, sigOutLength4);
  qmf_wrev (f2, sigInLength2, hiDecom, sigOutLength2);
  return;
}

/*-------------------------------------------
 * wavelet name parser
 *-----------------------------------------*/

void
wavelet_parser (char *wname, int *family, int *member)
{
  int count;

  *family = NOT_DEFINED;
  *member = NOT_DEFINED;
  for(count=0;count<waveletIdentityNum;count++)
    {
      if (strcmp(wname,wi[count].wname) == 0)
	{
	  *family = wi[count].family;
	  *member = wi[count].member;
	  break;
	}
    }
  return;
}

void
wavelet_fun_parser (char *wname, int *ii)
{
  int count;

  *ii = -1;
  for(count=0;count<waveletIdentityNum;count++)
    {
      if (strcmp(wname,wi[count].wname) == 0)
	{
	  *ii = count;
	  break;
	}
    }
  return;
}

/*void
wave_len_validate (int sigInLen, int waveLength, int *lev, int *val)
{
  int n;
  int m;
  int n1;
  int m1;
  float di;

  *val = 0;
  di = (float) sigInLen / (float) waveLength;
  if (di < 1)
    {
      *lev = 0;
      *val = 0;
      return;
    }
  else
    {
      n = (int) floor (log (di) / log ((float) 2));
      m = (int) ceil (log (di) / log ((float) 2));
      if ((((long) 1 << n) * waveLength == sigInLen)
	  || (((long) 1 << m) * waveLength == sigInLen))
	*lev = m + 1;
      else
	*lev = n + 1;
      *val = 1;
      //	  n1 = (int) floor (log (waveLength) / log ((float) 2));
      //m1 = (int) ceil (log (waveLength) / log ((float) 2));
      //if (n1 != m1)
      //    	  *lev = (*lev) - 1;
      if (di == 2)
	*lev -= 1;
      return;
    }
    }*/

void
wave_len_validate (int sigInLen, int waveLength, int *lev, int *val)
{
  int n;
  int m;

  *val = 0;
  if (sigInLen < 2*waveLength)
    {
      *lev = 0;
      *val = 0;
      return;
    }
  else
    {
      *val = 1;
      *lev = 0;
      n = sigInLen;
      do
	{
	  m = (int)floor((n + waveLength - 1)/2);
	  *lev += 1;
	  n = m;
	}
      while (m>=2*waveLength);
      return;
    }
}


void
dwt_print ()
{
  sciprint("\n**********************************************\n");
  switch (dwtMode) {
  case ZPD:
    {
      sciprint("**     DWT Extension Mode: Zero Padding     **\n");
      break;
    }
  case SYMH:
    {
      sciprint("** DWT Extension Mode: Half Symmetrization  **\n");
      break;
    }
  case SYMW:
    {
      sciprint("** DWT Extension Mode: Whole Symmetrization **\n");
      break;
    }
  case ASYMH:
    {
      sciprint("** DWT Extension Mode: Half Asymmetrization **\n");
      break;
    }
  case ASYMW:
    {
      sciprint("** DWT Extension Mode: Whole Asymmetrization**\n");
      break;
    }
  case SP0:
    {
      sciprint("** DWT Extension Mode: order 0 smooth padding*\n");
      break;
    }
  case SP1:
    {
      sciprint("** DWT Extension Mode: order 1 smooth padding*\n");
      break;
    }
  case PPD:
    {
      sciprint("**    DWT Extension Mode: Periodic Padding  **\n");
      break;
    }
  case PER:
    {
      sciprint("**    DWT Extension Mode: Periodization     **\n");
      break;
    }
  default:
    break;
  }
  sciprint("**********************************************\n");
  return;
}

void
dwt_write (char *mode, int *errCode)
{
  int count;
  *errCode = UNKNOWN_INPUT_ERR;
  
  for (count=0;count<extensionIdentityNum;count++)
    {
      if (strcmp(mode,ei[count].extMethodName) == 0)
	{
	  dwtMode = ei[count].extMethod;
	  *errCode = SUCCESS;
	  break;
	}
    }
  return;
}

void
dwt_parse(char **str1)
{
  int count;
  for (count=0;count<extensionIdentityNum;count++)
    {
      if (ei[count].extMethod == dwtMode)
	{
	  strcpy(*str1,ei[count].extMethodName);
	  break;
	}
    }
  return;
}

void
dwt_nex (double *sigIn, int sigInLength, double *lowDe, 
	 double *hiDe, int filterLen, double *approx, 
	 double *detail, int sigOutLength)
{
  //int sigInLengthTemp, 
  int sigOutLengthTemp, sigOutLengthPre;
  double *approxTemp, *approxPre; 
  double *detailTemp, *detailPre;

  //sigInLengthTemp = sigInLength + 2 * (filterLen - 1);
  //sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  //wextend_1D_center (sigIn, sigInLength, sigInTemp, 
  //     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLength + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, approxTemp, 
	sigOutLengthTemp, lowDe, filterLen);
  sigOutLengthPre = sigOutLengthTemp/2;
  approxPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 approxPre, sigOutLengthPre);
  wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);
  
  free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, detailTemp, 
	sigOutLengthTemp, hiDe, filterLen);
  detailPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 detailPre, sigOutLengthPre);
  wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);
  
  free(detailPre);
  free(detailTemp);
  //free(sigInTemp);
  return;
}


void
dwt (double *sigIn, int sigInLength, double *lowDe, 
     double *hiDe, int filterLen, double *approx, 
     double *detail, int sigOutLength, extend_method extMethod)
{
  int sigInLengthTemp, sigOutLengthTemp, sigOutLengthPre;
  double *sigInTemp, *approxTemp, *approxPre; 
  double *detailTemp, *detailPre;

  sigInLengthTemp = sigInLength + 2 * (filterLen - 1);
  sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  wextend_1D_center (sigIn, sigInLength, sigInTemp, 
		     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, approxTemp, 
	sigOutLengthTemp, lowDe, filterLen);
  sigOutLengthPre = sigOutLengthTemp/2;
  approxPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 approxPre, sigOutLengthPre);
  wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);
  
  free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, detailTemp, 
	sigOutLengthTemp, hiDe, filterLen);
  detailPre = malloc (sigOutLengthPre * sizeof(double));
  dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 detailPre, sigOutLengthPre);
  wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);
  
  free(detailPre);
  free(detailTemp);
  free(sigInTemp);
  return;
}

void
dwt_neo (double *sigIn, int sigInLength, double *lowDe, 
     double *hiDe, int filterLen, double *approx, 
     double *detail, int sigOutLength, extend_method extMethod)
{
  int sigInLengthTemp, sigOutLengthTemp, sigOutLengthPre;
  double *sigInTemp, *approxTemp, *approxPre; 
  double *detailTemp, *detailPre;

  sigInLengthTemp = sigInLength + 2 * filterLen;
  if ((extMethod == PER)&&(sigInLength%2 != 0))
	  sigInLengthTemp = sigInLength + 1 + 2 * filterLen;
  //if ((extMethod == PER)&&(sigInLength%2 == 0))
	//  sigInLengthTemp = sigInLength + 2 * filterLen;
  sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  wextend_1D_center (sigIn, sigInLength, sigInTemp, 
		     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, approxTemp, 
	sigOutLengthTemp, lowDe, filterLen);
  //sigOutLengthPre = sigOutLengthTemp/2;
  sigOutLengthPre = sigInLength + filterLen - 1;
  if ((extMethod==PER)&&(sigInLength%2 == 0))
	  sigOutLengthPre = sigInLength;
  if ((extMethod==PER)&&(sigInLength%2 !=0))
	  sigOutLengthPre = sigInLength + 1;
  approxPre = malloc (sigOutLengthPre * sizeof(double));
  wkeep_1D_center (approxTemp, sigOutLengthTemp, approxPre, sigOutLengthPre);
  dyaddown_1D_keep_even (approxPre, sigOutLengthPre,approx,sigOutLength);
  //dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 //approxPre, sigOutLengthPre);
  //wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);
  
  free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigInTemp, sigInLengthTemp, detailTemp, 
	sigOutLengthTemp, hiDe, filterLen);
  detailPre = malloc (sigOutLengthPre * sizeof(double));
  wkeep_1D_center (detailTemp, sigOutLengthTemp, detailPre, sigOutLengthPre);
  dyaddown_1D_keep_even (detailPre, sigOutLengthPre,detail,sigOutLength);
  //dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 //detailPre, sigOutLengthPre);
  //wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);
  
  free(detailPre);
  free(detailTemp);
  free(sigInTemp);
  return;
}


void
dwt_conv (double *sigIn, int sigInLength, double *lowDe, 
     double *hiDe, int filterLen, double *approx, 
     double *detail, int sigOutLength)
{
  
  //int sigOutLengthTemp;
  //double *approxTemp;
  //double *detailTemp;

  
  //sigOutLengthTemp = sigInLength + filterLen - 1;
  //approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  //conv (sigIn, sigInLength, approxTemp, 
	//sigOutLengthTemp, lowDe, filterLen);
  conv (sigIn, sigInLength, approx, 
	sigOutLength, lowDe, filterLen);
  
  //dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,approx,sigOutLength);
  
  //free(approxTemp);
  //detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  //conv (sigIn, sigInLength, detailTemp, 
	//sigOutLengthTemp, hiDe, filterLen);
  conv (sigIn, sigInLength, detail, 
	sigOutLength, hiDe, filterLen);
  
  //dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,detail,sigOutLength);
  
  //free(detailTemp);
  
  return;
}

void
dwt_no_extension (double *sigIn, int sigInLength, double *lowDe, 
     double *hiDe, int filterLen, double *approx, 
     double *detail, int sigOutLength)
{
  //int sigInLengthTemp, 
  int sigOutLengthTemp;
  //sigOutLengthPre;
  //double *sigInTemp, 
  double *approxTemp;//, *approxPre; 
  double *detailTemp;//, *detailPre;

  //sigInLengthTemp = sigInLength + 2 * filterLen;
  //if ((extMethod == PER)||(sigInLength%2 != 0))
	//  sigInLengthTemp = sigInLength + 1 + 2 * filterLen;
  //sigInTemp = malloc (sigInLengthTemp * sizeof (double));
  //wextend_1D_center (sigIn, sigInLength, sigInTemp, 
	//	     sigInLengthTemp, extMethod);
  sigOutLengthTemp = sigInLength + filterLen - 1;
  approxTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, approxTemp, 
	sigOutLengthTemp, lowDe, filterLen);
  //sigOutLengthPre = sigOutLengthTemp/2;
  //sigOutLengthPre = sigInLength + filterLen - 1;
  //approxPre = malloc (sigOutLengthPre * sizeof(double));
  //wkeep_1D_center (approxTemp, sigOutLengthTemp, approxPre, sigOutLengthPre);
  dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,approx,sigOutLength);
  //dyaddown_1D_keep_even (approxTemp, sigOutLengthTemp,
			 //approxPre, sigOutLengthPre);
  //wkeep_1D_center (approxPre, sigOutLengthPre, approx, sigOutLength);
  
  //free(approxPre);
  free(approxTemp);
  detailTemp = malloc (sigOutLengthTemp * sizeof(double));
  conv (sigIn, sigInLength, detailTemp, 
	sigOutLengthTemp, hiDe, filterLen);
  //detailPre = malloc (sigOutLengthPre * sizeof(double));
  //wkeep_1D_center (detailTemp, sigOutLengthTemp, detailPre, sigOutLengthPre);
  dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,detail,sigOutLength);
  //dyaddown_1D_keep_even (detailTemp, sigOutLengthTemp,
			 //detailPre, sigOutLengthPre);
  //wkeep_1D_center (detailPre, sigOutLengthPre, detail, sigOutLength);
  
  //free(detailPre);
  free(detailTemp);
  //free(sigInTemp);
  return;
}


void
idwt_complete_ex (double *approx, double *detail, int sigInLength, 
		  double *lowRe, double *hiRe, int filterLen, 
		  double *sigOut, int sigOutLength, 
		  extend_method extMethod) 
{
  int sigInLengthTemp, sigOutLengthTemp, count, ind, sigInLen;
  double *approxTemp, *detailTemp, *approxPre, *detailPre;
  double *sigOutPre, *approxEx, *detailEx;

  sigInLen = sigInLength + 2 * (filterLen - 1);
  approxEx = malloc(sigInLen*sizeof(double));
  detailEx = malloc(sigInLen*sizeof(double));

  wextend_1D_center (approx, sigInLength, approxEx, sigInLen,
		     extMethod);
  wextend_1D_center (detail, sigInLength, detailEx, sigInLen,
		     extMethod);
  
  sigInLengthTemp = 2 * sigInLen - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approxEx, sigInLen, 
		      approxTemp, sigInLengthTemp);
  dyadup_1D_feed_odd (detailEx, sigInLen, 
		      detailTemp, sigInLengthTemp);
  free(approxEx);
  free(detailEx);
  //printf("after dyadup!\n");
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  //printf("before conv!\n");
  //if ((!approxPre) || (!detailPre))
  //printf("Out of memory!\n");
  conv (approxTemp, sigInLengthTemp, approxPre, 
	sigOutLengthTemp, lowRe, filterLen);
  conv (detailTemp, sigInLengthTemp, detailPre, 
	sigOutLengthTemp, hiRe, filterLen);
  //printf("conv fin!\n");
  free(approxTemp);
  free(detailTemp);
  //printf("after conv!\n");
  sigOutPre = malloc(sigOutLengthTemp * sizeof(double));
  for(count=0;count<sigOutLengthTemp;count++)
    sigOutPre[count] = approxPre[count] + detailPre[count];
  free(approxPre);
  free(detailPre);
  //printf("before wkeep!\n");
  //wkeep_1D_center (sigOutPre, sigOutLengthTemp, 
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  wkeep_1D_index (sigOutPre, sigOutLengthTemp, 
  	  sigOut, sigOutLength, ind);
  //printf("sigLeng=%d,ind=%d\n",sigOutLengthTemp,ind);
  free(sigOutPre);
  //printf("leave idwt!\n");
  return;
}


void
idwt_complete (double *approx, double *detail, int sigInLength, 
	       double *lowRe, double *hiRe, int filterLen, 
	       double *sigOut, int sigOutLength) 
{
  int sigInLengthTemp, sigOutLengthTemp, count, ind;
  double *approxTemp, *detailTemp, *approxPre, *detailPre;
  double *sigOutPre;
  
  //printf("enter idwt!\n");
  sigInLengthTemp = 2 * sigInLength - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approx, sigInLength, 
		      approxTemp, sigInLengthTemp);
  dyadup_1D_feed_odd (detail, sigInLength, 
		      detailTemp, sigInLengthTemp);
  //printf("after dyadup!\n");
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  //printf("before conv!\n");
  //if ((!approxPre) || (!detailPre))
  //printf("Out of memory!\n");
  conv (approxTemp, sigInLengthTemp, approxPre, 
	sigOutLengthTemp, lowRe, filterLen);
  conv (detailTemp, sigInLengthTemp, detailPre, 
	sigOutLengthTemp, hiRe, filterLen);
  //printf("conv fin!\n");
  free(approxTemp);
  free(detailTemp);
  //printf("after conv!\n");
  sigOutPre = malloc(sigOutLengthTemp * sizeof(double));
  for(count=0;count<sigOutLengthTemp;count++)
    sigOutPre[count] = approxPre[count] + detailPre[count];
  free(approxPre);
  free(detailPre);
  //printf("before wkeep!\n");
  //wkeep_1D_center (sigOutPre, sigOutLengthTemp, 
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (sigOutPre, sigOutLengthTemp, 
  	  sigOut, sigOutLength, ind);
  //printf("sigLeng=%d,ind=%d\n",sigOutLengthTemp,ind);
  free(sigOutPre);
  //printf("leave idwt!\n");
  return;
}



void
idwt_neo (double *approx, double *detail, int sigInLength, 
	       double *lowRe, double *hiRe, int filterLen, 
	       double *sigOut, int sigOutLength) 
{
  int sigInLengthTemp, sigOutLengthTemp, count;//, ind;
  double *approxTemp, *detailTemp, *approxPre, *detailPre;
  double *sigOutPre;
  
  //printf("enter idwt!\n");
  sigInLengthTemp = 2 * sigInLength + 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_even (approx, sigInLength, 
		      approxTemp, sigInLengthTemp);
  dyadup_1D_feed_even (detail, sigInLength, 
		      detailTemp, sigInLengthTemp);
  //printf("after dyadup!\n");
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  //printf("before conv!\n");
  //if ((!approxPre) || (!detailPre))
  //printf("Out of memory!\n");
  conv (approxTemp, sigInLengthTemp, approxPre, 
	sigOutLengthTemp, lowRe, filterLen);
  conv (detailTemp, sigInLengthTemp, detailPre, 
	sigOutLengthTemp, hiRe, filterLen);
  //printf("conv fin!\n");
  free(approxTemp);
  free(detailTemp);
  //printf("after conv!\n");
  sigOutPre = malloc(sigOutLengthTemp * sizeof(double));
  for(count=0;count<sigOutLengthTemp;count++)
    sigOutPre[count] = approxPre[count] + detailPre[count];
  free(approxPre);
  free(detailPre);
  //printf("before wkeep!\n");
  wkeep_1D_center (sigOutPre, sigOutLengthTemp, 
     sigOut, sigOutLength);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  //wkeep_1D_index (sigOutPre, sigOutLengthTemp, 
  	  //sigOut, sigOutLength, ind);
  //printf("sigLeng=%d,ind=%d\n",sigOutLengthTemp,ind);
  free(sigOutPre);
  //printf("leave idwt!\n");
  return;
}


void
idwt_approx_ex (double *approx, int sigInLength, 
		double *lowRe, int filterLen, 
		double *sigOut, int sigOutLength,
		extend_method extMethod) 
{
  int sigInLengthTemp, sigOutLengthTemp, ind, sigInLen;
  double *approxTemp, *approxPre, *approxEx;

  sigInLen = sigInLength + 2 * (filterLen - 1);
  approxEx = malloc(sigInLen * sizeof(double));
  wextend_1D_center (approx, sigInLength, approxEx, sigInLen,
		     extMethod);

  sigInLengthTemp = 2 * sigInLen - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approxEx, sigInLen, 
		      approxTemp, sigInLengthTemp);
  free(approxEx);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (approxTemp, sigInLengthTemp, approxPre, 
	sigOutLengthTemp, lowRe, filterLen);
  free(approxTemp);

  //wkeep_1D_center (approxPre, sigOutLengthTemp, 
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (approxPre, sigOutLengthTemp, 
  	  sigOut, sigOutLength, ind);
  free(approxPre);
  return;
}



void
idwt_approx (double *approx, int sigInLength, 
	     double *lowRe, int filterLen, 
	     double *sigOut, int sigOutLength) 
{
  int sigInLengthTemp, sigOutLengthTemp, ind;
  double *approxTemp, *approxPre;

  sigInLengthTemp = 2 * sigInLength - 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (approx, sigInLength, 
		       approxTemp, sigInLengthTemp);
  
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (approxTemp, sigInLengthTemp, approxPre, 
	sigOutLengthTemp, lowRe, filterLen);
  free(approxTemp);

  //wkeep_1D_center (approxPre, sigOutLengthTemp, 
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (approxPre, sigOutLengthTemp, 
  	  sigOut, sigOutLength, ind);
  free(approxPre);
  return;
}

void
idwt_approx_neo (double *approx, int sigInLength, 
	     double *lowRe, int filterLen, 
	     double *sigOut, int sigOutLength) 
{
  int sigInLengthTemp, sigOutLengthTemp;//, ind;
  double *approxTemp, *approxPre;

  sigInLengthTemp = 2 * sigInLength + 1;
  approxTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_even (approx, sigInLength, 
		       approxTemp, sigInLengthTemp);
  
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  approxPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (approxTemp, sigInLengthTemp, approxPre, 
	sigOutLengthTemp, lowRe, filterLen);
  free(approxTemp);

  wkeep_1D_center (approxPre, sigOutLengthTemp, 
     sigOut, sigOutLength);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  //wkeep_1D_index (approxPre, sigOutLengthTemp, 
  	//  sigOut, sigOutLength, ind);
  free(approxPre);
  return;
}

void
idwt_detail (double *detail, int sigInLength, 
	     double *hiRe, int filterLen, 
	     double *sigOut, int sigOutLength) 
{
  int sigInLengthTemp, sigOutLengthTemp, ind;
  double *detailTemp, *detailPre;

  sigInLengthTemp = 2 * sigInLength - 1;
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (detail, sigInLength, 
		       detailTemp, sigInLengthTemp);
  
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (detailTemp, sigInLengthTemp, detailPre, 
	sigOutLengthTemp, hiRe, filterLen);
  free(detailTemp);

  //wkeep_1D_center (detailPre, sigOutLengthTemp, 
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (detailPre, sigOutLengthTemp, 
		  sigOut, sigOutLength, ind);
  free(detailPre);
  return;
}

void
idwt_detail_ex (double *detail, int sigInLength, 
		double *hiRe, int filterLen, 
		double *sigOut, int sigOutLength,
		extend_method extMethod) 
{
  int sigInLengthTemp, sigOutLengthTemp, ind, sigInLen;
  double *detailTemp, *detailPre, *detailEx;

  sigInLen = sigInLength + 2 * (filterLen - 1);
  detailEx = malloc(sigInLen * sizeof(double));
  wextend_1D_center (detail, sigInLength, detailEx, sigInLen,
		     extMethod);

  sigInLengthTemp = 2 * sigInLen - 1;
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_odd (detailEx, sigInLen, 
		      detailTemp, sigInLengthTemp);
  free(detailEx);
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (detailTemp, sigInLengthTemp, detailPre, 
	sigOutLengthTemp, hiRe, filterLen);
  free(detailTemp);

  //wkeep_1D_center (detailPre, sigOutLengthTemp, 
  //   sigOut, sigOutLength);
  ind = (int)(2 + (sigOutLengthTemp-sigOutLength)/2.0);
  wkeep_1D_index (detailPre, sigOutLengthTemp, 
		  sigOut, sigOutLength, ind);
  free(detailPre);
  return;
}


void
idwt_detail_neo (double *detail, int sigInLength, 
	     double *hiRe, int filterLen, 
	     double *sigOut, int sigOutLength) 
{
  int sigInLengthTemp, sigOutLengthTemp;//, ind;
  double *detailTemp, *detailPre;

  sigInLengthTemp = 2 * sigInLength + 1;
  detailTemp = malloc(sigInLengthTemp * sizeof(double));
  dyadup_1D_feed_even (detail, sigInLength, 
		       detailTemp, sigInLengthTemp);
  
  sigOutLengthTemp = sigInLengthTemp + filterLen - 1;
  detailPre = malloc (sigOutLengthTemp * sizeof(double));
  conv (detailTemp, sigInLengthTemp, detailPre, 
	sigOutLengthTemp, hiRe, filterLen);
  free(detailTemp);

  wkeep_1D_center (detailPre, sigOutLengthTemp, 
     sigOut, sigOutLength);
  //ind = 2 + (sigOutLengthTemp-sigOutLength)/2.0;
  //wkeep_1D_index (detailPre, sigOutLengthTemp, 
	//	  sigOut, sigOutLength, ind);
  free(detailPre);
  return;
}


void
wave_dec_len_cal (int filterLen, int sigLength,
		  int stride, int *waveDecLengthArray)
{
  int count = 0;
  int calLen;
  waveDecLengthArray[stride + 1] = sigLength;
  if (dwtMode!=PER)
    {
      calLen = sigLength;
      for (count = 0; count < stride; count++)
	{
	  calLen += (filterLen - 1);
	  waveDecLengthArray[stride-count]=(int)(floor(calLen/2));
	  calLen = *(waveDecLengthArray + stride - count);
	}
      waveDecLengthArray[0] = waveDecLengthArray[1];
    }
  else
    {
      for (count=stride; count > 0; count--)
	waveDecLengthArray[count] = 
	  (int)ceil(((double)(waveDecLengthArray[count+1]))/2.0);
     waveDecLengthArray[0] = waveDecLengthArray[1];
    }
  return;
}

void
upcoef_len_cal (int sigInLength, int filterLen, int stride, 
		int *sigOutLength, int *sigOutLengthDefault)
{
  int count;
  *sigOutLength = sigInLength;
  *sigOutLengthDefault = sigInLength;
  for(count=0;count<stride;count++)
    {
      *sigOutLengthDefault = 2*(*sigOutLengthDefault) + filterLen - 1;
      *sigOutLength = 2*(*sigOutLength) + filterLen - 2;
    }
  return;
}

void
upcoef (double *sigIn, int sigInLength, double *lowRe,double *hiRe, 
	int filterLen, double *sigOut, int sigOutLength, 
	int defaultLength, char *coefType, int step)
{
  int count, sigInLengthTemp, leng;
  double *sigInTemp, *sigOutTemp;

  sigInLengthTemp = 2 * sigInLength + filterLen - 2;
  //sigInLengthTemp = 2 * sigInLength + filterLen - 1;
  sigInTemp = malloc(defaultLength*sizeof(double));
  
  if (strcmp(coefType,"a")==0)
  {
	  //sciprint("recognized\n");
	  idwt_approx_neo (sigIn, sigInLength, lowRe, filterLen, 
		 sigInTemp, sigInLengthTemp);
  }
  else
    idwt_detail_neo (sigIn, sigInLength, hiRe, filterLen, 
		 sigInTemp, sigInLengthTemp); 

  if (step > 1)
    {
      sigOutTemp = malloc(defaultLength*sizeof(double));
      for(count=0;count<defaultLength;count++)
	sigOutTemp[count] = 0;
      leng = sigInLengthTemp;
      for(count=0;count<(step-1);count++)
	{
	  idwt_approx_neo (sigInTemp, leng, lowRe, filterLen,
		       sigOutTemp, leng*2+filterLen-2);
	  leng = leng*2+filterLen-2;
	  verbatim_copy (sigOutTemp, leng, sigInTemp, leng);
	}
      sigInLengthTemp = leng;
      free(sigOutTemp);
    }

 
  wkeep_1D_center (sigInTemp, sigInLengthTemp, sigOut, sigOutLength);
  free(sigInTemp);
  return;
}

void
wavedec (double *sigIn, int sigInLength, double *sigOut, 
	 int sigOutLength, double *lowDe, double *hiDe,
	 int filterLen, int *waveDecLengthArray, 
	 int lengthArrayLengh, int stride, extend_method extMethod)
{
  int count, pos, countt;
  int sigInLen;
  //int filterLen;
  double *pApprox, *pDetail, *pApproxTemp;
  double *pSig;

  pSig = sigIn;
  sigInLen = sigInLength;

  pApprox = malloc (sigInLength * sizeof (double));
  pApproxTemp = malloc (sigInLength * sizeof (double));
  for (count = 0; count < sigInLength; count++)
    {
      pApprox[count] = 0;
      pApproxTemp[count] = 0;
    }
  pos = sigOutLength - waveDecLengthArray[stride];
  pDetail = sigOut + pos;
  for (count = 0; count < stride; count++)
    {
      dwt_neo (pSig, sigInLen, lowDe, hiDe, filterLen, pApprox, 
	   pDetail, waveDecLengthArray[stride - count], extMethod);
      for (countt = 0; countt < waveDecLengthArray[stride - count]; countt++)
	pApproxTemp[countt] = pApprox[countt];
      pSig = pApproxTemp;
      sigInLen = waveDecLengthArray[stride - count];
      pos = pos - waveDecLengthArray[stride - count - 1];
      pDetail = (sigOut + pos);
    }
  for (count = 0; count < sigInLen; count++)
    sigOut[count] = pApprox[count];

  free (pApprox);
  free (pApproxTemp);
  return;
}


void
waverec (double *sigIn, int sigInLength, double *sigOut, 
	 int sigOutLength, double *lowRe, double *hiRe, 
	 int filterLen, int *waveDecLengthArray,
	 int lengthArraylength, int stride, 
	 extend_method extMethod)
{
  int count, pos, countt;
  int sigInLen;
  double *pApprox, *pDetail, *pApproxTemp;

  sigInLen = waveDecLengthArray[1];
  
  pApprox = malloc (sigOutLength * sizeof (double));
  pApproxTemp = malloc (sigOutLength * sizeof (double));
  for (count = 0; count < sigOutLength; count++)
    {
      pApprox[count] = 0;
      pApproxTemp[count] = 0;
    }
  pos = waveDecLengthArray[0];
  pDetail = sigIn + pos;
  for (count = 0; count < waveDecLengthArray[1]; count++)
    pApprox[count] = sigIn[count];

  for (count = 0; count < stride; count++)
    {
      idwt_neo (pApprox, pDetail, sigInLen, lowRe, hiRe, 
		     filterLen, pApproxTemp, 
		     waveDecLengthArray[count + 2]);
      for (countt = 0; countt < waveDecLengthArray[count + 2]; countt++)
	pApprox[countt] = pApproxTemp[countt];
      sigInLen = waveDecLengthArray[count + 2];
      pos += waveDecLengthArray[count + 1];
      pDetail = sigIn + pos;
    }
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = pApprox[count];

  free (pApprox);
  free (pApproxTemp);
  return;
}

void
wenergy (double *coef, int coefLen, int *lenArray, int arrayLen,
	 double *aE, int aELen, double *dE, int dELen)
{
  int count, countt, *pos;
  double ene;

  ene = 0;
  for(count=0;count<coefLen;count++)
    ene += coef[count]*coef[count];
  *aE = 0;
  for(count=0;count<lenArray[0];count++)
    *aE += coef[count]*coef[count];
  *aE = (*aE)*100/ene;

  pos = malloc (dELen*sizeof(int));
  for(count=0;count<dELen;count++)
    pos[count] = 0;
  pos[0] = lenArray[0];
  for(count=1;count<dELen;count++)
    pos[count] += (lenArray[count] + pos[count-1]); 
  for(count=0;count<dELen;count++)
    {
      dE[count] = 0;
      for(countt=0;countt<lenArray[count+1];countt++)
	dE[count] += coef[pos[count]+countt]*coef[pos[count]+countt];
      dE[count] = dE[count]*100/ene;
    }
  free(pos);
  return;
}

void
detcoef (double *sigIn, int sigInLength, int *waveDecLengthArray, 
	 int arrayLen, double *sigOut, int sigOutLength, 
	 int stride, int level)
{
  int leng, startCount, count;
  if (level != 0)
    {
      leng = 0;
      for (count = 0; count < level; count++)
	leng += waveDecLengthArray[stride - count];
    }

  startCount = sigInLength - leng;
  for (count = startCount; count <= (startCount + sigOutLength - 1); count++)
    sigOut[count - startCount] = sigIn[count];

  return;
}

void
appcoef (double *sigIn, int sigInLength, double *sigOut, 
	 int sigOutLength, double *lowRe, double *hiRe, 
	 int filterLen, int *waveDecLengthArray,
	 int lengthArraylength, int stride, int level, 
	 extend_method extMethod)
{
  int count, pos, countt;
  int sigInLen;
  double *pApprox, *pDetail, *pApproxTemp;

  if (level == stride)
    {
      for (count = 0; count < waveDecLengthArray[stride - level + 1]; count++)
	sigOut[count] = sigIn[count];
      return;
    }

  sigInLen = waveDecLengthArray[1];
  
  pApprox = malloc (sigOutLength * sizeof (double));
  pApproxTemp = malloc (sigOutLength * sizeof (double));
  for (count = 0; count < sigOutLength; count++)
    {
      pApprox[count] = 0;
      pApproxTemp[count] = 0;
    }
  pos = waveDecLengthArray[0];
  pDetail = sigIn + pos;
  for (count = 0; count < waveDecLengthArray[1]; count++)
    pApprox[count] = sigIn[count];

  for (count = 0; count < (stride - level); count++)
    {
      idwt_neo (pApprox, pDetail, sigInLen, lowRe, hiRe, 
		     filterLen, pApproxTemp, 
		     waveDecLengthArray[count + 2]);
      for (countt = 0; countt < waveDecLengthArray[count + 2]; countt++)
	pApprox[countt] = pApproxTemp[countt];
      sigInLen = waveDecLengthArray[count + 2];
      pos += waveDecLengthArray[count + 1];
      pDetail = sigIn + pos;
    }
  for (count = 0; count < sigOutLength; count++)
    sigOut[count] = pApprox[count];

  free (pApprox);
  free (pApproxTemp);
  return;
}

void
wrcoef (double *sigIn, int sigInLength, double *lowRe, double *hiRe,
	int filterLen, int *waveDecLengthArray, int arrayLen,
	double *sigOut, int sigOutLength, char *coefType, 
	int stride, int level, extend_method extMethod)
{

  int count = 0;
  int startCount, endCount, leng;
  double *sigOutTemp;

  sigOutTemp = malloc (sigInLength * sizeof (double));

  if (level != 0)
    {
      leng = 0;
      for (count = 0; count < level; count++)
	leng += waveDecLengthArray[stride - count];
    }

  if (strcmp (coefType, "d") == 0)
    {
      for (count = 0; count < sigInLength; count++)
	sigOutTemp[count] = 0;
      if (level != 0)
	{
	  startCount = sigInLength - leng;
	  endCount = startCount + waveDecLengthArray[stride - level + 1] - 1;
	  for (count = startCount; count <= endCount; count++)
	    sigOutTemp[count] = sigIn[count];
	}
    }
  else
    {
      for (count = 0; count < sigInLength; count++)
	sigOutTemp[count] = sigIn[count];
      if (level != 0)
	{
	  endCount = sigInLength - 1;
	  startCount = endCount - leng + 1;
	  for (count = startCount; count <= endCount; count++)
	    sigOutTemp[count] = 0;
	}
    }
  waverec (sigOutTemp, sigInLength, sigOut, sigOutLength, 
	   lowRe, hiRe, filterLen, waveDecLengthArray, arrayLen, 
	   stride, extMethod);
  //waverec (sigInLength, sigOutTemp, sigOutLength, sigOut,
  //   stride, waveDecLengthArray, waveType, mode);
  free (sigOutTemp);
  return;
}

void 
upwlev (double *coefArray, int coefLen, int *waveDecLengthArray,
	int arrayLen, double *lowRe, double *hiRe, int filterLen,
	double *newCoefArray, int newCoefLen, int *newLenArray, 
	int newArrayLen, double *approx, int approxLen, 
	int stride, extend_method extMethod)
{
  int count, pos1;
  char c='a';
  double *app, *det;

  //printf("enter upwlev!\n");
  for(count=0;count<approxLen;count++)
    approx[count]=coefArray[count];
  for(count=arrayLen-1;count>1;count--)
    newLenArray[count-1]=waveDecLengthArray[count];
  newLenArray[0]=newLenArray[1];

  pos1 = waveDecLengthArray[0] + waveDecLengthArray[1];
  for(count=coefLen-1;count>=pos1;count--)
    newCoefArray[count-coefLen+newCoefLen]=coefArray[count];

  app = malloc(waveDecLengthArray[1]*sizeof(double));
  det = malloc(waveDecLengthArray[1]*sizeof(double));
  for(count=0;count<waveDecLengthArray[1];count++)
    {
      app[count]=coefArray[count];
      det[count]=coefArray[count+waveDecLengthArray[1]];
    }
  idwt_neo (app, det, waveDecLengthArray[1], lowRe, hiRe, 
		 filterLen, newCoefArray, waveDecLengthArray[2]); 
  free(app);
  free(det);
  return;
}
