/*
 * -------------------------------------------------------------------------
 * swt_common.h -- Declarations for wavelet functions.
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2007  Roger Liu
 * Copyright (C) 20010-2012  Holger Nahrstaedt
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

#ifndef __SWT_H__
#define __SWT_H__ 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define __USE_DEPRECATED_STACK_FUNCTIONS__
#include <stack-c.h>
#include <api_scilab.h>
#include <sciprint.h>
#include <MALLOC.h>
#include <Scierror.h>

/*********************************************
 * Macro
 ********************************************/

#define SUCCESS      0
#define DIM_ERR_ONE  1
#define DIM_ERR_VEC  2
#define DIM_ERR_MAT  3


#define POSITIVE_INTEGER_ONLY                               1
#define LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION          2
#define SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION            3
#define OPT_CHAR_NOT_VALID                                  4
#define EXTENSION_OPT_NOT_VALID                             5
#define WAVELET_NAME_NOT_VALID                              6
#define DECOMPOSITION_LEVEL_NOT_VALID                       7
#define MULTI_DECOM_LEVEL_LESS_THAN_TWO                     8
#define WRONG_LHS                                           9
#define UNKNOWN_INPUT_ERR                                   20

#define PI   3.1415926535897931159980

/*********************************************
 * Macros CWT
 ********************************************/

#define REAL    0
#define COMPLEX 1

#define PHI_ONLY     0
#define PSI_ONLY     1
#define PHI_PSI_BOTH 2

#define SINUS           0
#define POISSON         1
#define MEXICAN_HAT     2
#define MORLET          3
#define DOGAUSS         4
#define CMORLET         5
#define SHANNON         6
#define FBSP            7
#define CAUCHY          8
#define GAUSS           9
#define CGAUSS          10


/*********************************************
 * Extension Type
 ********************************************/

typedef enum {
  ZPD, SYMH, SYMW, ASYMH, ASYMW, 
  SP0, SP1, PPD, PER} extend_method;

/*********************************************
 * Structure Declarations
 ********************************************/
#ifndef __USE_DEPRECATED_STACK_FUNCTIONS__
typedef struct sciintmat {
         int m,n;
         int it ; 
         int l;   
         void *D;     
} SciIntMat ;
#endif

typedef struct {
  int     sigInLength;
  int     sigOutLength;
  double  *sigIn;
  double  *sigOut;
  void    (*func)(int sigInLength, int sigOutLenght,
		  double *sigIn, double *sigOut);
  struct sio *link;
  } swt_sio;

typedef struct {
  char extMethodName[6];
  extend_method extMethod;
} extension_identity;


typedef struct {
int   errorNumber;
char  message[150];
} str_error_notification;

typedef struct hypermat {
  SciIntMat sc; /* coding informations */
  int it;
  int size;
  double *R;
  double *I;
} HyperMat;


/*********************************************
 * Structures CWT
 ********************************************/

typedef void(*WScaleFunc)(double *x, int sigInLength, double *psi, int sigOutLength, double ys);


typedef struct {
	char wname[20];
	int     realOrComplex;
	int     family;
	int     phipsi;
	double  lb;
	double  ub;
	double cpsi;
	WScaleFunc scalef;
} cwt_identity;

typedef struct {
	char wname[20];
	char     realOrComplex[20];
	char     family[20];
} cwt_family;


/*********************************************
 * Global Variables CWT
 ********************************************/

extern cwt_identity ci[];
extern cwt_family cif[];
extern int cwtFamilyNum;
extern int cwtIdentityNum;


/*********************************************
 * Function Prototype
 ********************************************/

/*------------------------------------------*/
/* Utility Function                         */
/* -----------------------------------------*/
extern void matrix_tran (double *matrixIn, int matrixInRow, 
			 int matrixInCol, double *matrixOut, 
			 int matrixOutRow, int matrixOutCol);
extern void wrev (const double *sigIn, int sigInLength, 
		  double *sigOut, int sigOutLength);
extern void qmf_even (const double *sigIn, int sigInLength, 
		      double *sigOut, int sigOutLength);
extern void qmf_odd (double *sigIn, int sigInLength, 
		     double *sigOut, int sigOutLength);
extern void qmf_wrev (const double *sigIn, int sigInLength, 
		      double *sigOut, int sigOutLength);
extern void verbatim_copy (const double *sigIn, int sigInLength, 
			   double *sigOut, int sigOutLength);
extern void dyaddown_1D_keep_odd (double *sigIn, int sigInLength, 
			       double *sigOut, int sigOutLength);
extern void dyaddown_1D_keep_even (double *sigIn, int sigInLength,
				double *sigOut, int sigOutLength);
extern void dyaddown_2D_keep_odd_row (double *matrixIn, 
				      int matrixInRow, 
				      int matrixInCol, 
				      double *matrixOut, 
				      int matrixOutRow, 
				      int matrixOutCol);
extern void dyaddown_2D_keep_odd_col (double *matrixIn, 
				      int matrixInRow, 
				      int matrixInCol, 
				      double *matrixOut, 
				      int matrixOutRow, 
				      int matrixOutCol);
extern void dyaddown_2D_keep_even_row (double *matrixIn, 
				       int matrixInRow, 
				       int matrixInCol, 
				       double *matrixOut, 
				       int matrixOutRow, 
				       int matrixOutCol);
extern void dyaddown_2D_keep_even_col (double *matrixIn, 
				       int matrixInRow, 
				       int matrixInCol, 
				       double *matrixOut, 
				       int matrixOutRow, 
				       int matrixOutCol);
extern void dyaddown_2D_keep_odd (double *matrixIn, 
				  int matrixInRow, 
				  int matrixInCol, 
				  double *matrixOut, 
				  int matrixOutRow, 
				  int matrixOutCol);
extern void dyaddown_2D_keep_even (double *matrixIn, 
				   int matrixInRow, 
				   int matrixInCol, 
				   double *matrixOut, 
				   int matrixOutRow, 
				   int matrixOutCol);
extern void dyadup_1D_feed_odd (double *sigIn, int sigInLength,
				double *sigOut, int sigOutLength);
extern void dyadup_1D_feed_even (double *sigIn, int sigInLength,
				 double *sigOut, int sigOutLength);
extern void dyadup_2D_feed_odd_row (double *matrixIn, 
				    int matrixInRow, 
				    int matrixInCol, 
				    double *matrixOut, 
				    int matrixOutRow, 
				    int matrixOutCol);
extern void dyadup_2D_feed_odd_col (double *matrixIn, 
				    int matrixInRow, 
				    int matrixInCol, 
				    double *matrixOut, 
				    int matrixOutRow, 
				    int matrixOutCol);
extern void dyadup_2D_feed_even_row (double *matrixIn, 
				     int matrixInRow, 
				     int matrixInCol, 
				     double *matrixOut, 
				     int matrixOutRow, 
				     int matrixOutCol);
extern void dyadup_2D_feed_even_col (double *matrixIn, 
				     int matrixInRow, 
				     int matrixInCol, 
				     double *matrixOut, 
				     int matrixOutRow, 
				     int matrixOutCol);
extern void dyadup_2D_feed_odd (double *matrixIn, 
				int matrixInRow, 
				int matrixInCol, 
				double *matrixOut, 
				int matrixOutRow, 
				int matrixOutCol);
extern void dyadup_2D_feed_even (double *matrixIn, 
				 int matrixInRow, 
				 int matrixInCol, 
				 double *matrixOut, 
				 int matrixOutRow, 
				 int matrixOutCol);
extern void extend_method_parse (char *mode, extend_method *extMethod);
extern void wextend_1D_center (double *sigIn, int sigInLength, 
			       double *sigOut, int sigOutLength, 
			       extend_method method);
extern void wextend_1D_left (double *sigIn, int sigInLength, 
			       double *sigOut, int sigOutLength, 
			       extend_method method);
extern void wextend_1D_right (double *sigIn, int sigInLength, 
			       double *sigOut, int sigOutLength, 
			       extend_method method);
extern void wextend_2D (double *matrixIn, int matrixInRow, 
			int matrixInCol, double *matrixOut, 
			int matrixOutRow, int matrixOutCol, 
			extend_method extMethod, char *rowOpt, 
			char *colOpt);
extern void wextend_2D_row (double *matrixIn, int matrixInRow, 
			    int matrixInCol, double *matrixOut, 
			    int matrixOutRow, int matrixOutCol, 
			    extend_method extMethod, char *Opt);
extern void wextend_2D_col (double *matrixIn, int matrixInRow, 
			    int matrixInCol, double *matrixOut, 
			    int matrixOutRow, int matrixOutCol, 
			    extend_method extMethod, char *Opt);
extern void wkeep_1D_center (double *sigIn, int sigInLength,
			     double *sigOut, int sigOutLength);
extern void wkeep_1D_left (double *sigIn, int sigInLength,
			   double *sigOut, int sigOutLength);
extern void wkeep_1D_right (double *sigIn, int sigInLength,
			    double *sigOut, int sigOutLength);
extern void wkeep_1D_index (double *sigIn, int sigInLength,
			    double *sigOut, int sigOutLength, 
			    int first);
extern void wkeep_2D_center (double *matrixIn, int matrixInRow, 
			     int matrixInCol, double *matrixOut,
			     int matrixOutRow, int matrixOutCol);
extern void wkeep_2D_index (double *matrixIn, int matrixInRow, 
			    int matrixInCol, double *matrixOut,
			    int matrixOutRow, int matrixOutCol,
			    int rowFirst, int colFirst);
extern void conv (double *sigIn, int sigInLength,
		  double *sigOut, int sigOutLength,
		  double *fiter, int filterLength);
extern void i_conv (double *sigIn, int sigInLength,
		  double *sigOut, int sigOutLength,
		  double *fiter, int filterLength);
 extern void swt_exp2(int lev, int *outputV);
 extern void linspace(double lb, double ub, int n, double *sigOut, int sigOutLength);
//extern void ocumsum (double *sigIn, int sigInLength);

 //extern int GetRhsVarDouble(int pos,  int *m1, int *n1, double *input);
/*------------------------------------------*/
/* Validation Function                      */
/* -----------------------------------------*/
//extern int is_scalar (int row, int col);
//extern int is_vector (int row, int col);
//extern int is_matrix (int row, int col);
extern void void_check (int number, int *type);
extern void scalar_check (int number, int *type);
extern void vector_check (int number, int *type);
extern void matrix_check (int number, int *type);
extern void real_or_complex (int number, int *type);
extern int sci_matrix_vector_real (int number);
extern int sci_matrix_vector_complex (int number);
extern int sci_matrix_matrix_complex (int number);
extern int sci_matrix_scalar_real (int number);
extern int sci_matrix_matrix_real (int number);
extern int sci_matrix_void (int number);
extern int sci_strings_scalar (int number);
extern int sci_strings_vector (int number);
extern int sci_mlist_check (int number);
//extern int sci_mlist_real (int number);
extern int scalar_string_check(char *l, char c);
extern int length_check(int number, int leng);
extern int vector_length_check(int number1, int number2);
extern void vector_length_compare(int number, int leng, int *res);
extern void matrix_length_compare(int number, int rowLeng, 
				  int colLeng, int *resRow, 
				  int *resCol);
extern int matrix_length_check (int number1, int number2);
extern int matrix_col_length_check(int number, int leng);
extern int matrix_row_length_check(int number, int leng);
extern void extension_check(char *mode, int *type);
extern void wavelet_family_check(char *wname, int wf, int *type);
extern void validate_print (int errCode);
extern void swt_max(double *sigIn, int sigInLength, double *sigMax);
extern void swt_min(double *sigIn, int sigInLength, double *sigMin);
extern void wcodemat_abs(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, int minv, int maxv);
extern void swt_max_abs(double *sigIn, int sigInLength, double *sigMax);
extern void swt_min_abs(double *sigIn, int sigInLength, double *sigMin);
extern double swt_abs(double sigIn);
extern void wcodemat(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, int minv, int maxv);
extern void wcodematd(double *sigIn, int sigInLength, double *sigOut, int sigOutLength, double minv, double maxv);
extern void wcodemat_matrix (double *matrixIn, int matrixInRow, int matrixInCol,
	                   double *matrixOut, int matrixOutRow, int matrixOutCol, 
					   int minv, int maxv, int abso);
extern void wcodemat_matrix_col (double *matrixIn, int matrixInRow, int matrixInCol,
	                       double *matrixOut, int matrixOutRow, int matrixOutCol,
						   int minv, int maxv, int abso);
extern void wcodemat_matrix_row (double *matrixIn, int matrixInRow, int matrixInCol,
	                       double *matrixOut, int matrixOutRow, int matrixOutCol,
						   int minv, int maxv, int abso);

//

extern void wrev_validate (int *errCode);
extern void wrev2_form_validate (int *errCode);
extern void qmf_validate (int *flow, int *errCode);
extern void conv_validate (int *errCode);
extern void dyaddown_form_validate (int *flow, int *errCode);
extern void dyaddown_content_validate (char *l, int *errCode);
extern void dyadup_form_validate (int *flow, int *errCode);
extern void dyadup_content_validate (char *l, int *errCode);
extern void wkeep_form_validate (int *flow, int *errCode);
extern void wkeep_content_validate (int flow, int *errCode, 
				    int l1, int l2, int l3);
extern void wextend_form_validate (int *flow, int *errCode, int l1);
extern void wextend_content_validate (int flow, int *errCode, int l2,
				      int l3, int l4, int l5, 
				      char **str);
extern void wcodemat_form_validate (int *flow, int *errCode);
extern void wcodemat_content_validate (int *errCode, int flow, int l1, int l2, int l3, int l4);
extern void wnorm_form_validate (int *flow, int *errCode);
extern void orthfilt_form_validate (int *errCode);
extern void biorfilt_form_validate (int *errCode);
extern void dbwavf_form_validate (int *errCode);
extern void dbwavf_content_validate (int *errCode, char *wname);
extern void coifwavf_form_validate (int *errCode);
extern void coifwavf_content_validate (int *errCode, char *wname);
extern void symwavf_form_validate (int *errCode);
extern void symwavf_content_validate (int *errCode, char *wname);
extern void legdwavf_form_validate (int *errCode);
extern void legdwavf_content_validate (int *errCode, char *wname);
extern void biorwavf_form_validate (int *errCode);
extern void biorwavf_content_validate (int *errCode, char *wname);
extern void rbiorwavf_form_validate (int *errCode);
extern void rbiorwavf_content_validate (int *errCode, char *wname);
extern void wfilters_form_validate (int *errCode, int *flow, int l2);
extern void wfilters_content_validate(int *errCode, char *wname);
extern void wmaxlev_form_validate (int *errCode);
extern void dwt_form_validate(int *errCode, int *flow);
extern void dwt_content_validate(int *errCode, int flow, int l1, 
				 int l2, int l3, int l4, int l5);
extern void idwt_form_validate (int *errCode, int *flow);
extern void idwt_content_validate (int *errCode, int flow, int l1, 
				   int l2, int l3, int l4, int l5, 
				   int l6, int l7);
extern void wavedec_form_validate(int *errCode, int *flow);
extern void wavedec_content_validate(int *errCode, int flow, int l1,
				     int l2, int l3, int l4);
extern void waverec_form_validate(int *errCode, int *flow);
extern void waverec_content_validate(int *errCode, int flow, int l1,
				     int l2, int l3, int l4);
extern void wrcoef_form_validate (int *errCode, int *flow);
extern void appcoef_form_validate (int *errCode, int *flow);
extern void detcoef_form_validate (int *errCode, int *flow);
extern void wenergy_form_validate (int *errCode);
extern void appcoef_form_validate (int *errCode, int *flow);
extern void appcoef_content_validate (int *errCode, int flow, 
				      int l1, int l2,
				      int l3, int l4, int l5);
extern void wrcoef_content_validate (int *errCode, int flow, int l1, 
				     int l2, int l3, int l4, int l5, 
				     int l6);
extern void upcoef_form_validate (int *errCode, int *flow);
extern void uplwev_form_validate (int *errCode, int *flow);
extern void upcoef_content_validate (int *errCode, int flow, int l1, 
				     int l2, int l3, int l4, int l5, 
				     int l6);
extern void upwlev_form_validate (int *errCode, int *flow);
extern void upwlev_content_validate (int *errCode, int flow, int l1, 
				     int l2, int l3, int l4);
extern void dwt2_form_validate (int *errCode, int *flow);
extern void idwt2_form_validate (int *errCode, int *flow);
extern void wavedec2_form_validate (int *errCode, int *flow);
extern void wavedec2_content_validate (int *errCode, int flow, 
				       int l1, int l2, int l3, 
				       int l4);
extern void waverec2_form_validate (int *errCode, int *flow);
extern void wenergy2_form_validate (int *errCode, int *flow);
//extern void dwt2_form_validate (int *errCode, int *flow);
extern void dwt2_content_validate (int *errCode, int flow, int l1, 
				   int l2, int l3, int l4, int l5);
//extern void idwt2_form_validate (int *errCode, int *flow);
extern void idwt2_content_validate (int *errCode, int flow, int l1, 
				    int l2, int l3, int l4, int l5, 
				    int l6, int l7, int l8, int l9);
extern void waverec2_content_validate (int *errCode, int flow,int l1,
				       int l2, int l3, int l4);

extern void detcoef2_form_validate (int *errCode, int *flow);
extern void detcoef2_content_validate (int *errCode, int flow, 
				       int l1, int l2, int l3, 
				       int l4);
extern void appcoef2_form_validate (int *errCode, int *flow);
extern void appcoef2_content_validate (int *errCode, int flow, 
				       int l1, int l2, int l3, 
				       int l4, int l5);
extern void wrcoef2_form_validate (int *errCode, int *flow);
extern void wrcoef2_content_validate (int *errCode, int flow, 
				      int l1, int l2, int l3, 
				      int l4, int l5, int l6);
extern void upwlev2_form_validate (int *errCode, int *flow);
extern void upwlev2_content_validate (int *errCode, int flow, 
				      int l1, int l2,
				      int l3, int l4);
extern void upcoef2_form_validate (int *errCode, int *flow);
extern void upcoef2_content_validate (int *errCode, int flow, 
				      int l1, int l2, int l3, 
				      int l4, int l5, int l6);
extern void swt_form_validate (int *errCode, int *flow);
extern void swt_content_validate (int *errCode, int flow, 
				      int l1, int l2, int l3, int l4);
extern void iswt_form_validate (int *errCode, int *flow);
extern void iswt_content_validate (int *errCode, int flow, 
				      int l1, int l2, int l3, int l4);
extern void swt2_form_validate (int *errCode, int *flow);
extern void swt2_content_validate (int *errCode, int flow, 
				      int l1, int l2, int l3, int l4);
extern void iswt2_form_validate (int *errCode, int *flow);
extern void iswt2_content_validate (int *errCode, int flow, int l1, int l2, int l3, int l4, int l5, int l6);

extern void dwt3_form_validate (int *errCode, int *flow);
extern void dwt3_content_validate (int *errCode, int flow, int l1, int l2,
		       int l3, int l4, int l5, int l6, int l7,
		       int l8, int l9);
extern void idwt3_form_validate (int *errCode, int *flow);
extern void idwt3_content_validate (int *errCode, int flow, int l1, int l2,
		       int l3, int l4, int l5, int l6, int l7,
		       int l8);
extern void dualtree_form_validate (int *errCode, int *flow);
extern void dualtree_content_validate (int *errCode, int flow, int l1, int l2,
			   int l3, int l4, int l5, int l6);
extern void idualtree_form_validate (int *errCode, int *flow);
extern void dualtree2D_form_validate (int *errCode, int *flow);
extern void idualtree2D_form_validate (int *errCode, int *flow);
extern void icplxdual2D_form_validate (int *errCode, int *flow);

#endif //__SWT_H__
