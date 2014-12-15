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
#include <stdarg.h>
#include <math.h>
#include <string.h>
//#define __USE_DEPRECATED_STACK_FUNCTIONS__
//#include <stack-c.h>
#include "api_scilab.h"
#include <sciprint.h>
#include "MALLOC.h"
#include "Scierror.h"
#include "swtlib.h"
#include "swt_gwsupport.h"

/*********************************************
 * Global Variables CWT
 ********************************************/
//
extern cwt_identity ci[];
extern cwt_family cif[];
extern int cwtFamilyNum;
extern int cwtIdentityNum;


/*********************************************
 * Function Prototype
 ********************************************/


 //extern int GetRhsVarDouble(int pos,  int *m1, int *n1, double *input);
/*------------------------------------------*/
/* Validation Function                      */
/* -----------------------------------------*/
//extern int is_scalar (int row, int col);
//extern int is_vector (int row, int col);
//extern int is_matrix (int row, int col);
// extern int void_check (int number, int *type);
// extern int scalar_check (int number, int *type);
// extern int vector_check (int number, int *type);
// extern int matrix_check (int number, int *type);
// extern int real_or_complex (int number, int *type);
// extern int sci_matrix_vector_real (int number);
// extern int sci_matrix_vector_complex (int number);
// extern int sci_matrix_matrix_complex (int number);
// extern int sci_matrix_scalar_real (int number);
// extern int sci_matrix_matrix_real (int number);
// extern int sci_matrix_void (int number);
// extern int sci_strings_scalar (int number);
// extern int sci_strings_vector (int number);
// extern int sci_mlist_check (int number);
// //extern int sci_mlist_real (int number);
// extern int scalar_string_check(char *l, char c);
// extern int length_check(int number, int leng);
// extern int vector_length_check(int number1, int number2);
// extern int vector_length_compare(int number, int leng, int *res);
// extern int matrix_length_compare(int number, int rowLeng,
// 				  int colLeng, int *resRow,
// 				  int *resCol);
// extern int matrix_length_check (int number1, int number2);
// extern int matrix_col_length_check(int number, int leng);
// extern int matrix_row_length_check(int number, int leng);
// extern void extension_check(char *mode, int *type);
// extern void wavelet_family_check(char *wname, int wf, int *type);
// extern void validate_print (int errCode);



//
// //
//
#endif //__SWT_H__
