/*
 * -------------------------------------------------------------------------
 * dwt.h -- Declarations for wavelet functions.
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

/*********************************************
 * Macros
 ********************************************/

#define HAAR           0
#define DAUBECHIES     1
#define COIFLETS       2
#define SYMLETS        3
#define SPLINE_BIORTH  4
#define BEYLKIN        5
#define VAIDYANATHAN   6
#define DMEY           7
#define BATHLETS       8
#define LEGENDRE       9
#define SPLINE_RBIORTH 10
#define FARRAS         11
#define KINGSBURYQ     12
#define NOT_DEFINED    99

#define ORTH       0
#define BIORTH     1


/*********************************************
 * Wavelet Structure Declarations
 ********************************************/

typedef struct {
  int     length;
  double  *pLowPass;
  double  *pHiPass;
} swt_wavelet;

typedef void(*Func)(int member, swt_wavelet *pWaveStruct);

typedef struct {
  char  wname[20];
  int   rOrB;
  int   family;
  int   member;
  Func  analysis;
  Func  synthesis;
} wavelet_identity;



/*********************************************
 * Global Variable Declaration
 ********************************************/

extern double LowDecomFilCoef[80];
extern double LowReconFilCoef[80];
extern double HiDecomFilCoef[80];
extern double HiReconFilCoef[80];
extern wavelet_identity wi[];
extern int waveletIdentityNum;
extern extend_method dwtMode;
extern extension_identity ei[];
extern int extensionIdentityNum;
extern str_error_notification strErrNoti[];


//extern wavelet_identity wi[];


/*********************************************
 * Function Prototype
 ********************************************/

extern void haar_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void haar_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void daubechies_analysis_initialize (int memeber, swt_wavelet *pWaveStruct);
extern void daubechies_synthesis_initialize (int memeber, swt_wavelet *pWaveStruct);
extern void symlets_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void symlets_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void coiflets_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void coiflets_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void sp_bior_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void sp_bior_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void sp_rbior_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void sp_rbior_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void beylkin_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void beylkin_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void vaidyanathan_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void vaidyanathan_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void dmey_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void dmey_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void bathlets_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void bathlets_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void legendre_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void legendre_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void farras_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void farras_synthesis_initialize (int member, swt_wavelet *pWaveStruct);
extern void kingsburyq_analysis_initialize (int member, swt_wavelet *pWaveStruct);
extern void kingsburyq_synthesis_initialize (int member, swt_wavelet *pWaveStruct);

/*------------------------------------------*/
/* Wavelet Family Function                  */
/* -----------------------------------------*/
extern void filter_clear ();
extern void orth_filt_group (double *filterIn, int sigInLength,
			     double *filterLowRec, 
			     double *filterLowDec,
			     double *filterHiRec, 
			     double *filterHiDec);
extern void bior_filt_group (double *f1, int sigInLength1, 
			     double *f2, int sigInLength2, 
			     double *lowDecom, int sigOutLength1, 
			     double *hiDecom, int sigOutLength2,
			     double *lowRecon, int sigOutLength3, 
			     double *hiRecon, int sigOutLength4);
extern void wavelet_parser (char *wname, int *family, int *member);
extern void wavelet_fun_parser (char *wname, int *ii);
extern void wave_len_validate (int sigInLen, int waveLength, int *lev, int *val);
extern void dwt_print();
extern void dwt_write (char *mode, int *errCode);
extern void dwt_parse(char **strr);
extern void dwt (double *sigIn, int sigInLength, double *lowDe, 
		 double *hiDe, int filterLen, double *approx, 
		 double *detail, int sigOutLength, 
		 extend_method extMethod);
extern void dwt_neo (double *sigIn, int sigInLength, double *lowDe, 
		 double *hiDe, int filterLen, double *approx, 
		 double *detail, int sigOutLength, 
		 extend_method extMethod);
extern void dwt_nex (double *sigIn, int sigInLength, double *lowDe, 
		     double *hiDe, int filterLen, double *approx, 
		     double *detail, int sigOutLength);
extern void dwt_no_extension (double *sigIn, int sigInLength, double *lowDe, 
     double *hiDe, int filterLen, double *approx, 
     double *detail, int sigOutLength);
extern void dwt_conv (double *sigIn, int sigInLength, double *lowDe, 
     double *hiDe, int filterLen, double *approx, 
     double *detail, int sigOutLength);
extern void idwt_complete (double *approx, double *detail, 
			   int sigInLength, double *lowRe, 
			   double *hiRe, int filterLen, 
			   double *sigOut, int sigOutLength);
extern void idwt_neo (double *approx, double *detail, 
			   int sigInLength, double *lowRe, 
			   double *hiRe, int filterLen, 
			   double *sigOut, int sigOutLength);
extern void idwt_complete_ex (double *approx, double *detail, 
			      int sigInLength, double *lowRe, 
			      double *hiRe, int filterLen, 
			      double *sigOut, int sigOutLength, 
			      extend_method extMethod);
extern void idwt_approx (double *approx, int sigInLength, 
			 double *lowRe, int filterLen, 
			 double *sigOut, int sigOutLength);
extern void idwt_approx_ex (double *approx, int sigInLength, 
			    double *lowRe, int filterLen, 
			    double *sigOut, int sigOutLength,
			    extend_method extMethod); 
extern void idwt_approx_neo (double *approx, int sigInLength, 
			    double *lowRe, int filterLen, 
			    double *sigOut, int sigOutLength);
extern void idwt_detail (double *detail, int sigInLength, 
			 double *hiRe, int filterLen, 
			 double *sigOut, int sigOutLength);
extern void idwt_detail_ex (double *detail, int sigInLength, 
			    double *hiRe, int filterLen, 
			    double *sigOut, int sigOutLength,
			    extend_method extMethod); 
extern void idwt_detail_neo (double *detail, int sigInLength, 
			    double *hiRe, int filterLen, 
			    double *sigOut, int sigOutLength); 
extern void wave_dec_len_cal (int filterLen, int sigLength,
			      int stride, int *waveDecLengthArray);
extern void wavedec (double *sigIn, int sigInLength, double *sigOut, 
		     int sigOutLength, double *lowDe, double *hiDe,
		     int filterLen, int *waveDecLengthArray, 
		     int lengthArrayLengh, int stride, 
		     extend_method extMethod);
extern void waverec (double *sigIn, int sigInLength, double *sigOut, 
		     int sigOutLength, double *lowRe, double *hiRe, 
		     int filterLen, int *waveDecLengthArray,
		     int lengthArraylength, int stride, 
		     extend_method extMethod);
extern void wenergy (double *coef, int coefLen, int *lenArray, 
		     int arrayLen, double *aE, int aELen, 
		     double *dE, int dELen);
extern void detcoef (double *sigIn, int sigInLength, 
		     int *waveDecLengthArray, int arrayLen, 
		     double *sigOut, int sigOutLength, 
		     int stride, int level);
extern void appcoef (double *sigIn, int sigInLength, double *sigOut, 
		     int sigOutLength, double *lowRe, double *hiRe, 
		     int filterLen, int *waveDecLengthArray,
		     int lengthArraylength, int stride, int level, 
		     extend_method extMethod);
extern void wrcoef (double *sigIn, int sigInLength, double *lowRe, 
		    double *hiRe, int filterLen, 
		    int *waveDecLengthArray, int arrayLen,
		    double *sigOut, int sigOutLength, 
		    char *coefType, int stride, int level, 
		    extend_method extMethod);
extern void upcoef_len_cal (int sigInLength, int filterLen, 
			    int stride, int *sigOutLength,
			    int *sigOutLengthDefault);
extern void upwlev (double *coefArray, int coefLen, 
		    int *waveDecLengthArray,	int arrayLen, 
		    double *lowRe, double *hiRe, int filterLen,
		    double *newCoefArray, int newCoefLen, 
		    int *newLenArray, int newArrayLen, 
		    double *approx, int approxLen, int stride,
		    extend_method extMethod);
extern void upcoef (double *sigIn, int sigInLength, double *lowRe,
		    double *hiRe, int filterLen, double *sigOut, 
		    int sigOutLength, int defaultLength, 
		    char *coefType, int step);
extern void dwt2D (double *matrixIn, int matrixInRow, 
		   int matrixInCol, double *matrixOutApprox, 
		   double *matrixOutColDetail, 
		   double *matrixOutRowDetail, 
		   double *matrixOutDetail, int matrixOutRow, 
		   int matrixOutCol, double *lowDe, double *hiDe, 
		   int filterLen, extend_method extMethod);
extern void
dwt2D_neo_a (double *matrixIn, int matrixInRow, int matrixInCol,
	     double *matrixOutApprox, double *matrixOutColDetail,
	     double *matrixOutRowDetail, double *matrixOutDetail,
	     int matrixOutRow, int matrixOutCol, double *lowDeR, 
	     double *hiDeR, double *lowDeC, double *hiDeC,
	     int filterLen, extend_method extMethod);
extern void dwt2D_neo (double *matrixIn, int matrixInRow, 
		   int matrixInCol, double *matrixOutApprox, 
		   double *matrixOutColDetail, 
		   double *matrixOutRowDetail, 
		   double *matrixOutDetail, int matrixOutRow, 
		   int matrixOutCol, double *lowDe, double *hiDe, 
		   int filterLen, extend_method extMethod);
extern void idwt2D (double *matrixInApprox, 
		    double *matrixInColDetail,
		    double *matrixInRowDetail, 
		    double *matrixInDetail,
		    int matrixInRow, int matrixInCol, double *lowRe,
		    double *hiRe, int filterLen, double *matrixOut,
		    int matrixOutRow, int matrixOutCol, 
		    extend_method extMethod);
extern void idwt2D_neo (double *matrixInApprox, double *matrixInColDetail,
	double *matrixInRowDetail, double *matrixInDetail,
	int matrixInRow, int matrixInCol, double *lowRe,
	double *hiRe, int filterLen, double *matrixOut,
	int matrixOutRow, int matrixOutCol);
extern void
idwt2D_neo_a (double *matrixInApprox, double *matrixInColDetail,
	      double *matrixInRowDetail, double *matrixInDetail,
	      int matrixInRow, int matrixInCol, double *lowReR,
	      double *hiReR, double *lowReC, double *hiReC,
	      int filterLen, double *matrixOut,
	      int matrixOutRow, int matrixOutCol);

extern void wave_mem_cal (int *pLen, int stride, int *total);
extern void matrix_wavedec_len_cal (int matrixInRow, int matrixInCol,
				    int stride, int filterLen, 
				    int *pLen);
extern void matrix_locate (int stride, int *pLen, int *pH, 
			   int *pV, int *pD);
extern void wavedec2 (double *matrixIn, int matrixInRow, 
		      int matrixInCol, double *lowDe, double *hiDe, 
		      int filterLen, int *pLen, double *coef, 
		      int sigOutLength, int stride, 
		      extend_method extMethod);
extern void
wavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,
	   double *lowDeR, double *hiDeR, double *lowDeC,
	   double *hiDeC, int filterLen, int *pLen, 
	   double *coef, int sigOutLength, int stride, 
	   extend_method extMethod);
extern void waverec2 (double *coef, int sigInLength, double *lowRe, 
		      double *hiRe, int filterLen, double *matrixOut,
		      int matrixOutRow, int matrixOutCol, int *pLen, 
		      int stride, extend_method extMethod);
extern void
waverec2a (double *coef, int sigInLength, double *lowReR, 
	   double *hiReR, double *lowReC, double *hiReC, 
	   int filterLen, double *matrixOut, int matrixOutRow, 
	   int matrixOutCol, int *pLen, int stride, 
	   extend_method extMethod);

extern void wenergy_2output (double *coef, int sigInLength, 
			    int *pLen, double *ae, double *de, 
			    int deLength, int stride);
extern void wenergy_4output (double *coef, int sigInLength, 
			     int *pLen, double *ae, double *he, 
			     double *ve, double *de, int deLength, 
			     int stride);
extern void detcoef2 (double *coef, int sigInLength, double *coefOut,
		      int sigOutLength, int *pLen, int stride, 
		      int level, char *coefType);
extern void appcoef2 (double *coef, int sigInLength, double *lowRe, 
		      double *hiRe, int filterLen, double *coefOut, 
		      int matrixOutRow, int matrixOutCol, int *pLen, 
		      int stride, int level, extend_method extMethod);
extern void wrcoef2 (double *coef, int sigInLength, double *lowRe, 
		     double *hiRe, int filterLen, double *matrixOut, 
		     int matrixOutRow, int matrixOutCol, int *pLen, 
		     int stride, int level, char *type, 
		     extend_method extMethod);
extern void upwlev2 (double *coef, int sigInLength, double *lowRe, 
		     double *hiRe,
	 int filterLen, int *pLen, int matrixRow, int matrixCol,
	 double *approx, int approxLen, double *newCoef, 
	 int newCoefLen, int *newLenMatrix, int lenMatrixRow, 
	 int lenMatrixCol, int stride, extend_method extMethod);
extern void upcoef2 (double *matrixIn, int matrixInRow, 
		     int matrixInCol, double *lowRe, double *hiRe, 
		     int filterLen, double *matrixOut, 
		     int matrixOutRow, int matrixOutCol,
		     int matrixOutDefaultRow, 
		     int matrixOutDefaultCol,
		     int step, char *type);//, extend_method extMethod);

extern void dwt3d_tran(double *mat3DIn, int row1, int col1, int sli1,
		       double *mat3DOut, int row2, int col2, int sli2);

extern void dwt3d_line_forward(double *mat3DIn, int row1, int col1, int sli1,
			double *mat3DOutApp, double *mat3DOutDet,
			int row2, int col2, int sli2, 
			double *loDe, double *hiDe, int filterLen,
			extend_method extMethod);

extern void dwt3d_tran_z(double *mat3DIn, int row1, int col1, int sli1,
			 double *mat3DOut, int row2, int col2, int sli2);

extern void dwt3d_tran_z_inv(double *mat3DIn, int row1, int col1, int sli1,
		      double *mat3DOut, int row2, int col2, int sli2);

extern void dwt3d_combine(double *mat1, double *mat2, double *mat3,
		   double *mat4, double *mat5, double *mat6,
		   double *mat7, double *mat8, int rowIn, 
		   int colIn, int sliIn, double *matOut,
		   int rowOut, int colOut, int sliOut);
extern void dwt3d_line_reverse(double *mat3DInApp, double *mat3DInDet,
			int row1, int col1, int sli1, 
			double *mat3DOut, int row2, int col2, 
			int sli2,
			       double *loDe, double *hiDe, int filterLen);


extern void dwt3d_split(double *matIn, int rowIn, int colIn, int sliIn,
		 double *mat1, double *mat2, double *mat3,
		 double *mat4, double *mat5, double *mat6,
		 double *mat7, double *mat8, int rowOut, 
			int colOut, int sliOut);
extern void dwt3(double *mat3DIn, int row, int col, int sli,
	  double *mat3DOut, int row2, int col2, int sli2,
	  int r, int c, int s, double *Lo1, double *Hi1, 
	  double *Lo2, double *Hi2, double *Lo3, double *Hi3,
		 int fLen1, int fLen2, int fLen3, extend_method extMethod);
