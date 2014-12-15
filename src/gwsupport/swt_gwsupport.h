

// Copyright (C) 2012 - Michael Baudin
// Copyright (C) 2008-2009 - INRIA - Michael Baudin
// Copyright (C) 2009 - Digiteo - Michael Baudin
//
// This file must be used under the terms of the GNU Lesser General Public License license :
// http://www.gnu.org/copyleft/lesser.html

//
// gwsupport.h
//   Header for the C gateway support functions for DISTFUN
//
#ifndef __SCI_SWT_GWSUPPORT_H__
#define __SCI_SWT_GWSUPPORT_H__

#ifdef _MSC_VER
#if LIBSWTGWSUPPORT_EXPORTS
#define SWT_GWSUPPORT_IMPORTEXPORT __declspec (dllexport)
#else
#define SWT_GWSUPPORT_IMPORTEXPORT __declspec (dllimport)
#endif
#else
#define SWT_GWSUPPORT_IMPORTEXPORT
#endif

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS  "C" {
  # define __END_DECLS }
  #else
  # define __BEGIN_DECLS /* empty */
  # define __END_DECLS /* empty */
  #endif

  __BEGIN_DECLS
  //#define __USE_DEPRECATED_STACK_FUNCTIONS__
  static int SWT_GWSUPPORT_OK = 1;
  static int SWT_GWSUPPORT_ERROR = 0;
  static str_error_notification strErrNoti[] ={
    {POSITIVE_INTEGER_ONLY, "Input integer must be positive!\n"},
    {LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION, "Length Parameter is not valid for input vector!\n"},
    {SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION, "Size Parameter is not valid for input matrix!\n"},
    {OPT_CHAR_NOT_VALID, "Optional Charactor Parameter is not valid!\n"},
    {EXTENSION_OPT_NOT_VALID, "Extension Method is not valid!\n"},
    {WAVELET_NAME_NOT_VALID, "Wavelet Name is not valid!\n"},
    {DECOMPOSITION_LEVEL_NOT_VALID, "Input signal is not valid for selected decompostion level and wavelets!\n"},
    {MULTI_DECOM_LEVEL_LESS_THAN_TWO ,"Decomposition level must be no less than 2!\n"},
    {UNKNOWN_INPUT_ERR, "Unrecognized Input and Output Numbers or parameter not valid for the algorithm!\nPlease refer to help pages!\n"},
    {WRONG_LHS, "Output parameter quantity is inappropriate!\n"}
  };
  static int errorNum = sizeof(strErrNoti)/sizeof(str_error_notification);
  // typedef struct hypermat {
  //   SciIntMat sc; /* coding informations */
  //   int it;
  //   int size;
  //   double *R;
  //   double *I;
  // } HyperMat;
  // cwt_validate.c
  void sinus_content_validate(int *errCode, double* input1, double* input2, int* input3);
  void poisson_form_validate(int *errCode);
  void poisson_content_validate(int *errCode, double* input1, double* input2, int* input3);
  void Gauss_content_validate(int *errCode, double* input1, double* input2, int* input3, int* input4);
  void shanwavf_content_validate(int *errCode, double* input1, double* input2, int* input3, double* input4 , double* input5);
  void DOGauss_content_validate(int *errCode, double* input1, double* input2, int* input3);
  void mexihat_content_validate(int *errCode, double* input1, double* input2, int* input3);
  void cmorlet_content_validate(int *errCode, double* input1, double* input2, int* input3, double* input4 , double* input5);
  void cauchy_content_validate(int *errCode, double* input1, double* input2, int* input3);
  void fbspwavf_content_validate(int *errCode, double* input1, double* input2, int* input3, int* input4, double* input5, double* input6);
  void morlet_content_validate(int *errCode, double* input1, double* input2, int* input3);
  void wavefun_content_validate(int *errCode, char *input_string, int* input2);
  void wavefun2_content_validate(int *errCode, char *input_string, int* input2);


  void dwt_print();

  int void_check (int number, int *type);
  int scalar_check (int number, int *type);
  int vector_check (int number, int *type);
  int matrix_check (int number, int *type);
  int real_or_complex (int number, int *type);
  int sci_matrix_vector_real (int number);
  int sci_matrix_vector_complex (int number);
  int sci_matrix_matrix_complex (int number);
  int sci_matrix_scalar_real (int number);
  int sci_matrix_matrix_real (int number);
  int sci_matrix_void (int number);
  int sci_strings_scalar (int number);
  int sci_strings_vector (int number);
  int sci_mlist_check (int number);
  // int sci_mlist_real (int number);
  int scalar_string_check(char *l, char c);
  int length_check(int number, int leng);
  int vector_length_check(int number1, int number2);
  int vector_length_compare(int number, int leng, int *res);
  int matrix_length_compare(int number, int rowLeng, int colLeng, int *resRow, int *resCol);
  int matrix_length_check (int number1, int number2);
  int matrix_col_length_check(int number, int leng);
  int matrix_row_length_check(int number, int leng);
  void extension_check(char *mode, int *type);
  void wavelet_family_check(char *wname, int wf, int *type);
  void validate_print (int errCode);

   void wrev_validate (int *errCode);
   void wrev2_form_validate (int *errCode);
   void qmf_validate (int *flow, int *errCode);
   void conv_validate (int *errCode);
   void dyaddown_form_validate (int *flow, int *errCode);
   void dyaddown_content_validate (char *l, int *errCode);
   void dyadup_form_validate (int *flow, int *errCode);
   void dyadup_content_validate (char *l, int *errCode);
   void wkeep_form_validate (int *flow, int *errCode);
   void   wkeep_content_validate (int flow, int *errCode,
   double *input1, int *int_input2, int *int_input3);
   void wkeep_content_validate_string (int flow, int *errCode,
   double *input1, int *int_input2, char *input_string);
   void wextend_form_validate (int *flow, int *errCode, char* l1,int* int_l1);

void  wextend_content_validate (int flow, int *errCode, char *input_string1,
                double *input3, int* int_intput4, char *input_string2, char **str);
   void wcodemat_form_validate (int *flow, int *errCode);
   void wcodemat_content_validate (int *errCode, int flow, int* int_input2);
   void wnorm_form_validate (int *flow, int *errCode);
   void orthfilt_form_validate (int *errCode);
   void biorfilt_form_validate (int *errCode);
   void dbwavf_form_validate (int *errCode);
   void dbwavf_content_validate (int *errCode, char *wname);
   void coifwavf_form_validate (int *errCode);
   void coifwavf_content_validate (int *errCode, char *wname);
   void symwavf_form_validate (int *errCode);
   void symwavf_content_validate (int *errCode, char *wname);
   void legdwavf_form_validate (int *errCode);
   void legdwavf_content_validate (int *errCode, char *wname);
   void biorwavf_form_validate (int *errCode);
   void biorwavf_content_validate (int *errCode, char *wname);
   void rbiorwavf_form_validate (int *errCode);
   void rbiorwavf_content_validate (int *errCode, char *wname);
   void wfilters_form_validate (int *errCode, int *flow, char *input_string1);
   void wfilters_content_validate(int *errCode, char *wname);
   void wmaxlev_form_validate (int *errCode);
   void dwt_form_validate(int *errCode, int *flow);
   void dwt_content_validate(int *errCode, int flow, char* string1,char* string2,char* string3);
   void idwt_form_validate (int *errCode, int *flow);

  void idwt_content_validate (int *errCode, int flow,
             char * l3, int* l4, int* l5, char* l4_string, char* l5_string, char* l6, char* l7);
   void wavedec_form_validate(int *errCode, int *flow);
   void wavedec_content_validate(int *errCode, int flow, int* int1, char* string1);
   void waverec_form_validate(int *errCode, int *flow);
   void waverec_content_validate(int *errCode, int flow,char* l3);
   void wrcoef_form_validate (int *errCode, int *flow);
   void appcoef_form_validate (int *errCode, int *flow);
   void detcoef_form_validate (int *errCode, int *flow);
   void wenergy_form_validate (int *errCode);
   void appcoef_form_validate (int *errCode, int *flow);

   void      appcoef_content_validate (int *errCode, int flow, char* l3);
    void wrcoef_content_validate (int *errCode, int flow, char *l1,  char* l4, int *l5, int* l6);
   void upcoef_form_validate (int *errCode, int *flow);
   void uplwev_form_validate (int *errCode, int *flow);

   void   upcoef_content_validate (int *errCode, int flow, char* l1,
               char* l3, int* l4, int* l5, int* l6);
   void upwlev_form_validate (int *errCode, int *flow);

   void         upwlev_content_validate (int *errCode, int flow,char* l3);
   void dwt2_form_validate (int *errCode, int *flow);
   void idwt2_form_validate (int *errCode, int *flow);
   void wavedec2_form_validate (int *errCode, int *flow);
   void   wavedec2_content_validate (int *errCode, int flow, int* l2,   char* l3);

   void waverec2_form_validate (int *errCode, int *flow);
   void wenergy2_form_validate (int *errCode, int *flow);
  // void dwt2_form_validate (int *errCode, int *flow);
  void  dwt2_content_validate (int *errCode, int flow,  char* l2,
  char*  l3, char*  l4, char*  l5);
  // void idwt2_form_validate (int *errCode, int *flow);
   void idwt2_content_validate (int *errCode, int flow,char* l5, int* l6,char* l6_char, int* l7,char* l7_char,char* l8,char* l9);
   void waverec2_content_validate (int *errCode, int flow,char *l3);

   void detcoef2_form_validate (int *errCode, int *flow);
   void detcoef2_content_validate (int *errCode, int flow,  char* l1);
   void appcoef2_form_validate (int *errCode, int *flow);
   void appcoef2_content_validate (int *errCode, int flow,char* l3, int* l4, int* l5);
   void wrcoef2_form_validate (int *errCode, int *flow);
   void wrcoef2_content_validate (int *errCode, int flow, char* l1,  char* l4, int* l5, int* l6);
   void upwlev2_form_validate (int *errCode, int *flow);
   void upwlev2_content_validate (int *errCode, int flow,char *l3);
   void upcoef2_form_validate (int *errCode, int *flow);
   void upcoef2_content_validate (int *errCode, int flow,char* l1,  char* l3, int* l4, int* l5, int* l6);
   void swt_form_validate (int *errCode, int *flow);
   void swt_content_validate (int *errCode, int flow,int* l2, char* l3);
   void iswt_form_validate (int *errCode, int *flow);
   void iswt_content_validate (int *errCode, int flow,char* l2, char* l3);
   void swt2_form_validate (int *errCode, int *flow);
   void swt2_content_validate (int *errCode, int flow,int* l2, char* l3);
   void iswt2_form_validate (int *errCode, int *flow);
   void iswt2_content_validate (int *errCode, int flow, char* l2, char* l5);

   void dwt3_form_validate (int *errCode, int *flow);
   void dwt3_content_validate (int *errCode, int flow, char* l2,  char* l3, char* l4, char* l5, char* l6,      char* l8, char* l9);
   void idwt3_form_validate (int *errCode, int *flow);
   void idwt3_content_validate (int *errCode, int flow, char* l2,	char* char_l3,int* l3, char* char_l4,int* l4, int* l5,
   int* l8);
   void dualtree_form_validate (int *errCode, int *flow);
   void dualtree_content_validate (int *errCode, int flow, int l1, int l2,
  			   int l3, int l4, int l5, int l6);
 void idualtree_form_validate (int *errCode, int *flow);
void dualtree2D_form_validate (int *errCode, int *flow);
 void idualtree2D_form_validate (int *errCode, int *flow);
void icplxdual2D_form_validate (int *errCode, int *flow);

int swt_gwsupport_GetRealMatrixOfDoubles( char * fname, int ivar, int * _piRows , int * _piCols, double** _pdblReal );
int swt_gwsupport_GetCompexMatrixOfDoubles( char * fname, int ivar, int * _piRows , int * _piCols, double** _pdblReal , double** _pdblImg);
int swt_gwsupport_GetRealMatrixOfDoublesAsInteger( char * fname, int ivar, int * _piRows , int * _piCols, int** _pdblReal );
int swt_gwsupport_GetScalarString( char * fname, int ivar , char** mystring );
int swt_gwsupport_GetRealHypermatofdouble( char * fname, int ivar, int ** _dims , int *_ndims, double** _pdblReal );
int swt_gwsupport_AllocMatrixOfDoubles ( char * fname, int ovar , int _piRows , int _piCols , double ** _pdblReal );
int swt_gwsupport_AllocMatrixOfDoublesAsInteger ( char * fname, int ovar , int _piRows , int _piCols , int ** _pdblReal );
int swt_gwsupport_CreateMatrixOfString ( char * fname, int ovar , int _piRows , int _piCols , char ** pstData );
int swt_gwsupport_CreateHypermatOfDouble ( char * fname, int ovar , int* dims , int ndims, double* pdblReal);
int swt_gwsupport_AllocComplexMatrixOfDoubles ( char * fname, int ovar , int _piRows , int _piCols , double** _pdblReal, double** _pdblImg );
char** swt_gwsupport_GetMatrixOfString( char * fname, int ivar, int * _piRows , int * _piCols);
//int swt_gwsupport_GetRealMatrixOfUnsignedInteger16( char * fname, int ivar, int * _piRows , int * _piCols, unsigned short** _pusData16 );
int swt_gwsupport_AllocMatrixOfUnsignedInteger16 ( char * fname, int ovar , int _piRows , int _piCols , unsigned short** _pusData16);
int swt_gwsupport_AllocMatrixOfInteger32 ( char * fname, int ovar , int _piRows , int _piCols , int** _piData32);
int swt_gwsupport_GetType(int ivar);
int swt_gwsupport_GetMatrixdims(int ivar, int * _piRows , int * _piCols);
int swt_gwsupport_IsVarComplex(int ivar);
  __END_DECLS


  /* ==================================================================== */

  #endif /* __SCI_GWSUPPORT_H__ */
