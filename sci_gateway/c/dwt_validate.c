/*
 * -------------------------------------------------------------------------
 * dwt_validate.c -- DWT validation
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


/*-------------------------------------------
 * orthfilt validation 
 *-----------------------------------------*/

void
orthfilt_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if (!sci_matrix_vector_real(1))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

/*-------------------------------------------
 * biorthfilt validation 
 *-----------------------------------------*/

void
biorfilt_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if ((!sci_matrix_vector_real(1)) || (!sci_matrix_vector_real(2)))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

/*-------------------------------------------
 * dbwavf validation 
 *-----------------------------------------*/

void
dbwavf_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if (!sci_strings_scalar(1))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
dbwavf_content_validate (int *errCode, char *wname)
{
  int type;
  *errCode = SUCCESS;
  wavelet_family_check(wname,DAUBECHIES,&type);
  if (!type)
    *errCode = WAVELET_NAME_NOT_VALID;
  return;
}


/*-------------------------------------------
 * coifwavf validation 
 *-----------------------------------------*/

void
coifwavf_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if (!sci_strings_scalar(1))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
coifwavf_content_validate (int *errCode, char *wname)
{
  int type;
  *errCode = SUCCESS;
  wavelet_family_check(wname,COIFLETS,&type);
  if (!type)
    *errCode = WAVELET_NAME_NOT_VALID;
  return;
}



/*-------------------------------------------
 * symwavf validation 
 *-----------------------------------------*/

void
symwavf_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if (!sci_strings_scalar(1))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
symwavf_content_validate (int *errCode, char *wname)
{
  int type;
  *errCode = SUCCESS;
  wavelet_family_check(wname,SYMLETS,&type);
  if (!type)
    *errCode = WAVELET_NAME_NOT_VALID;
  return;
}

/*-------------------------------------------
 * legdwavf validation 
 *-----------------------------------------*/

void
legdwavf_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if (!sci_strings_scalar(1))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
legdwavf_content_validate (int *errCode, char *wname)
{
  int type;
  *errCode = SUCCESS;
  wavelet_family_check(wname,LEGENDRE,&type);
  if (!type)
    *errCode = WAVELET_NAME_NOT_VALID;
  return;
}


/*-------------------------------------------
 * biorwavf validation 
 *-----------------------------------------*/

void
biorwavf_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if (!sci_strings_scalar(1))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
biorwavf_content_validate (int *errCode, char *wname)
{
  int type;
  *errCode = SUCCESS;
  wavelet_family_check(wname,SPLINE_BIORTH,&type);
  if (!type)
    *errCode = WAVELET_NAME_NOT_VALID;
  return;
}

/*-------------------------------------------
 * rbiorwavf validation 
 *-----------------------------------------*/

void
rbiorwavf_form_validate (int *errCode)
{
  *errCode = SUCCESS;
  if (!sci_strings_scalar(1))
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
rbiorwavf_content_validate (int *errCode, char *wname)
{
  int type;
  *errCode = SUCCESS;
  wavelet_family_check(wname,SPLINE_RBIORTH,&type);
  if (!type)
    *errCode = WAVELET_NAME_NOT_VALID;
  return;
}


/*-------------------------------------------
 * wfilters validation 
 *-----------------------------------------*/

void
wfilters_form_validate (int *errCode, int *flow, int l2)
{

  *errCode = SUCCESS;
  if ((Rhs==2) && (!sci_strings_scalar(2)))
    {
      *errCode = UNKNOWN_INPUT_ERR;
      return;
    }
  if ((Rhs==1) && sci_strings_scalar(1) && (Lhs==4))
    *flow = 1;
  else if ((Rhs==2) && sci_strings_scalar(1) && (Lhs==2) &&
	   cstk(l2)[0]=='d')
    *flow = 2;
  else if ((Rhs==2) && sci_strings_scalar(1) && (Lhs==2) &&
	   cstk(l2)[0]=='r')
    *flow = 3;
  else if ((Rhs==2) && sci_strings_scalar(1) && (Lhs==2) &&
	   cstk(l2)[0]=='l')
    *flow = 4;
  else if ((Rhs==2) && sci_strings_scalar(1) && (Lhs==2) &&
	   cstk(l2)[0]=='h')
    *flow = 5;
  else
	  *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void 
wfilters_content_validate(int *errCode, char *wname)
{
  int typ1, typ2, typ3, typ4, typ5, typ6, typ7, typ8, typ9, typ10, typ11;
  *errCode = SUCCESS;
  wavelet_family_check(wname,DAUBECHIES,&typ1);
  wavelet_family_check(wname,COIFLETS,&typ2);
  wavelet_family_check(wname,SYMLETS,&typ3);
  wavelet_family_check(wname,SPLINE_BIORTH,&typ4);
  wavelet_family_check(wname,HAAR,&typ5);
  wavelet_family_check(wname,BEYLKIN,&typ6);
  wavelet_family_check(wname,VAIDYANATHAN,&typ7);
  wavelet_family_check(wname,DMEY,&typ8);
  wavelet_family_check(wname,BATHLETS,&typ9);
  wavelet_family_check(wname,LEGENDRE,&typ10);
  wavelet_family_check(wname,SPLINE_RBIORTH,&typ11);
  //wavelet_family_check(wname,FARRAS,&typ12);
  //wavelet_family_check(wname,KINGSBURYQ,&typ13);
  if ((!typ1) && (!typ2) && (!typ3) && (!typ4) && (!typ5) && (!typ6) && (!typ7) && (!typ8) && (!typ9) && (!typ10) && (!typ11))
    *errCode = WAVELET_NAME_NOT_VALID;
  return;
}

/*-------------------------------------------
 * wmaxlev validation 
 *-----------------------------------------*/

void 
wmaxlev_form_validate(int *errCode)
{
  *errCode = UNKNOWN_INPUT_ERR;
  if (sci_matrix_scalar_real(1) && sci_strings_scalar(2))
    *errCode = SUCCESS;
  else if (sci_matrix_vector_real(1) && 
	   sci_strings_scalar(2) && length_check(1,2))
    *errCode = SUCCESS;
  return;
}

/*-------------------------------------------
 * dwt validation 
 *-----------------------------------------*/

void
dwt_form_validate(int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==2) && sci_matrix_vector_real(1) &&
      sci_strings_scalar(2))
    *flow = 1;
  else if ((Rhs==3) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && 
	   sci_matrix_vector_real(3) && vector_length_check(2,3))
    *flow = 2;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_strings_scalar(2) && sci_strings_scalar(3) && 
	   sci_strings_scalar(4))
    *flow = 3;
  else if ((Rhs==5) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   vector_length_check(2,3) &&
	   sci_strings_scalar(4) && sci_strings_scalar(5))
    *flow = 4;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
dwt_content_validate(int *errCode, int flow, int l1, 
		     int l2, int l3, int l4, int l5)
{
  int type;

  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode, cstk(l2));
      break;
    }
  case 2:
    {
      break;
    }
  case 3:
    {
      wfilters_content_validate(errCode, cstk(l2));
      extension_check(cstk(l4),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      if (strcmp(cstk(l3),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      break;
    }
  case 4:
    {
      extension_check(cstk(l5),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      if (strcmp(cstk(l4),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      break;
    }
  default:
    break;
  }
  return;
}


/*-------------------------------------------
 * idwt validation 
 *-----------------------------------------*/

void
idwt_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
 
  if ((Rhs==3) && sci_matrix_vector_real(1) &&
      sci_matrix_vector_real(2) && sci_strings_scalar(3) &&
      vector_length_check(1,2))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && vector_length_check(1,2) &&
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4))
    *flow = 2;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && vector_length_check(1,2) &&
	   sci_strings_scalar(3) && sci_matrix_scalar_real(4))
    *flow = 3;
  else if ((Rhs==5) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && vector_length_check(1,2) &&
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4) && sci_matrix_scalar_real(5))
    *flow = 4;
  else if ((Rhs==5) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && vector_length_check(1,2) &&
	   sci_strings_scalar(3) && sci_strings_scalar(4) &&
	   sci_strings_scalar(5))
    *flow = 5;
  else if ((Rhs==6) && sci_matrix_vector_real(1) && 
	   sci_matrix_vector_real(2) && vector_length_check(1,2) &&
	   sci_strings_scalar(3) && sci_matrix_scalar_real(4) &&
	   sci_strings_scalar(5) && sci_strings_scalar(6))
    *flow = 6;
  else if ((Rhs==6) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && vector_length_check(1,2) &&
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4) && sci_strings_scalar(5) &&
	   sci_strings_scalar(6))
    *flow = 7;
  else if ((Rhs==7) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && vector_length_check(1,2) &&
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4) && sci_matrix_scalar_real(5) &&
	   sci_strings_scalar(6) && sci_strings_scalar(7))
    *flow = 8;
  else if ((Rhs==3) && 
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_strings_scalar(3))
    *flow = 1;
  else if ((Rhs==4) &&
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4))
    *flow = 2;
  else if ((Rhs==4) && 
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_strings_scalar(3) && sci_matrix_scalar_real(4))
    *flow = 3;
  else if ((Rhs==5) && 
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4) && sci_matrix_scalar_real(5))
    *flow = 4;
  else if ((Rhs==5) && 
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_strings_scalar(3) && sci_strings_scalar(4) &&
	   sci_strings_scalar(5))
    *flow = 5;
  else if ((Rhs==6) && 
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_strings_scalar(3) && sci_matrix_scalar_real(4) &&
	   sci_strings_scalar(5) && sci_strings_scalar(6))
    *flow = 6;
  else if ((Rhs==6) && 
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4) && sci_strings_scalar(5) &&
	   sci_strings_scalar(6))
    *flow = 7;
  else if ((Rhs==7) && 
	   ((sci_matrix_vector_real(1) && sci_matrix_void(2)) || 
	    (sci_matrix_void(1) && sci_matrix_vector_real(2))) && 
	   sci_matrix_vector_real(3) && sci_matrix_vector_real(4) &&
	   vector_length_check(3,4) && sci_matrix_scalar_real(5) &&
	   sci_strings_scalar(6) && sci_strings_scalar(7))
    *flow = 8;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
idwt_content_validate (int *errCode, int flow, int l1, int l2,
		       int l3, int l4, int l5, int l6, int l7)
{
  int type;

  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode, cstk(l3));
      break;
    }
  case 2:
    {
      break;
    }
  case 3:
    {
      wfilters_content_validate(errCode, cstk(l3));
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 4:
    {
      if (istk(l5)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 5:
    {
      wfilters_content_validate(errCode, cstk(l3));
      if (strcmp(cstk(l4),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l5),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 6:
    {
      wfilters_content_validate(errCode, cstk(l3));
      if (strcmp(cstk(l5),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l6),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 7:
    {
      if (strcmp(cstk(l5),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l6),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 8:
    {
      if (strcmp(cstk(l6),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l7),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      if (istk(l5)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * wavedec validation 
 *-----------------------------------------*/

void
wavedec_form_validate(int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && (sci_matrix_vector_real(1))&& 
      sci_matrix_scalar_real(2) && sci_strings_scalar(3))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_scalar_real(2) && 
	   sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) &&
	   vector_length_check(3,4))
    *flow = 2;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
wavedec_content_validate(int *errCode, int flow, int l1,
			 int l2, int l3, int l4)
{
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      if (istk(l2)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 2:
    {
      if (istk(l2)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * waverec validation 
 *-----------------------------------------*/

void
waverec_form_validate(int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && (sci_matrix_vector_real(1))&& 
      sci_matrix_vector_real(2) && sci_strings_scalar(3))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && 
	   sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) &&
	   vector_length_check(3,4))
    *flow = 2;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
waverec_content_validate(int *errCode, int flow, int l1,
			 int l2, int l3, int l4)
{
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 2:
    {
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * wrcoef validation 
 *-----------------------------------------*/

void
wrcoef_form_validate(int *errCode, int *flow)
{
  *errCode = SUCCESS;
  
  if ((Rhs==4) && sci_strings_scalar(1) &&
      sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
      sci_strings_scalar(4))
    *flow = 1;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   sci_strings_scalar(4) && sci_matrix_scalar_real(5))
    *flow = 2;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && sci_matrix_vector_real(5) &&
	   vector_length_check(4,5))
    *flow = 3;
  else if ((Rhs==6) && sci_strings_scalar(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && sci_matrix_vector_real(5) &&
	   sci_matrix_scalar_real(6) && vector_length_check(4,5))
    *flow = 4;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
wrcoef_content_validate (int *errCode, int flow, int l1, int l2,
			 int l3, int l4, int l5, int l6)
{
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l4));
      if (scalar_string_check(cstk(l1),'a') || 
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = UNKNOWN_INPUT_ERR;
      break;
    }
  case 2:
    {
      wfilters_content_validate(errCode,cstk(l4));
      if (scalar_string_check(cstk(l1),'a') || 
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = UNKNOWN_INPUT_ERR;
      if ((scalar_string_check(cstk(l1),'a') && (istk(l5)[0] >= 0)) ||
	  (scalar_string_check(cstk(l1),'d') && (istk(l5)[0] > 0) ))
	*errCode = SUCCESS;
      else
	*errCode = UNKNOWN_INPUT_ERR;
      break;
    }
  case 3:
    {
      if (scalar_string_check(cstk(l1),'a') || 
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      break;
    }
  case 4:
    {
      if (scalar_string_check(cstk(l1),'a') || 
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = UNKNOWN_INPUT_ERR;
      if ((scalar_string_check(cstk(l1),'a') && (istk(l6)[0] >= 0)) ||
	  (scalar_string_check(cstk(l1),'d') && (istk(l6)[0] > 0) ))
	*errCode = SUCCESS;
      else
	*errCode = UNKNOWN_INPUT_ERR;
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * appcoef validation 
 *-----------------------------------------*/

void
appcoef_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_matrix_vector_real(1) &&
      sci_matrix_vector_real(2) && sci_strings_scalar(3))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && sci_strings_scalar(3) &&
	   sci_matrix_scalar_real(4) )
    *flow = 2;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4))
    *flow = 3;
  else if ((Rhs==5) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4) &&
	   sci_matrix_scalar_real(5))
    *flow = 4;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
appcoef_content_validate (int *errCode, int flow, int l1, int l2,
			  int l3, int l4, int l5)
{
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 2:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 3:
    {
      break;
    }
  case 4:
    {
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * detcoef validation 
 *-----------------------------------------*/

void
detcoef_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==2) && sci_matrix_vector_real(1) && 
      sci_matrix_vector_real(2))
    *flow = 1;
  else if ((Rhs==3) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_scalar_real(3))
    *flow = 2;
  else 
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

/*-------------------------------------------
 * wenergy validation 
 *-----------------------------------------*/

void
wenergy_form_validate (int *errCode)
{
  *errCode = UNKNOWN_INPUT_ERR;
  if ((Rhs==2) && sci_matrix_vector_real(1) && 
      sci_matrix_vector_real(2))
    *errCode = SUCCESS;
  return;
}

/*-------------------------------------------
 * upcoef validation 
 *-----------------------------------------*/

void
upcoef_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_strings_scalar(1) && 
      (sci_matrix_vector_real(2) || sci_matrix_scalar_real(2)) && sci_strings_scalar(3))
    *flow = 5;
  else if ((Rhs==4) && sci_strings_scalar(1) &&
	   (sci_matrix_vector_real(2) || sci_matrix_scalar_real(2)) && sci_strings_scalar(3) &&
	   sci_matrix_scalar_real(4))
    *flow = 1;
  else if ((Rhs==4) && sci_strings_scalar(1) && 
	   (sci_matrix_vector_real(2) || sci_matrix_scalar_real(2)) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4))
    *flow = 6;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   (sci_matrix_vector_real(2) || sci_matrix_scalar_real(2)) && sci_strings_scalar(3) &&
	   sci_matrix_scalar_real(4) && sci_matrix_scalar_real(5))
    *flow = 2;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   (sci_matrix_vector_real(2) || sci_matrix_scalar_real(2)) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4) &&
	   sci_matrix_scalar_real(5))
    *flow = 3;
  else if ((Rhs==6) && sci_strings_scalar(1) &&
	   (sci_matrix_vector_real(2) || sci_matrix_scalar_real(2)) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4) &&
	   sci_matrix_scalar_real(5) && sci_matrix_scalar_real(6))
    *flow = 4;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
upcoef_content_validate (int *errCode, int flow, int l1, int l2,
			 int l3, int l4, int l5, int l6)
{
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      if (scalar_string_check(cstk(l1),'a') ||
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      wfilters_content_validate(errCode,cstk(l3));
      if (istk(l4)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 2:
    {
      if (scalar_string_check(cstk(l1),'a') ||
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      wfilters_content_validate(errCode,cstk(l3));
      if (istk(l4)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      if (istk(l5)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 3:
    {
      if (scalar_string_check(cstk(l1),'a') ||
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      if (istk(l5)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 4:
    {
      if (scalar_string_check(cstk(l1),'a') ||
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      if (istk(l5)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      if (istk(l6)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 5:
    {
      if (scalar_string_check(cstk(l1),'a') ||
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 6:
    {
      if (scalar_string_check(cstk(l1),'a') ||
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * upwlev validation 
 *-----------------------------------------*/

void
upwlev_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_matrix_vector_real(1) && 
      sci_matrix_vector_real(2) && sci_strings_scalar(3))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4))
    *flow = 2;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void 
upwlev_content_validate (int *errCode, int flow, int l1, int l2,
			 int l3, int l4)
{
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 2:
    {
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * dwt2 validation 
 *-----------------------------------------*/

void
dwt2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==2) && sci_matrix_matrix_real(1) && 
      sci_strings_scalar(2))
    *flow = 1;
  else if ((Rhs==3) && sci_matrix_matrix_real(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   vector_length_check(2,3))
    *flow = 2;
  else if ((Rhs==4) && sci_matrix_matrix_real(1) &&
	   sci_strings_scalar(2) && sci_strings_scalar(3) &&
	   sci_strings_scalar(4))
    *flow = 3;
  else if ((Rhs==5) && sci_matrix_matrix_real(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_vector_real(3) &&
	   vector_length_check(2,3) && sci_strings_scalar(4) &&
	   sci_strings_scalar(5))
    *flow = 4;
  else 
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
dwt2_content_validate (int *errCode, int flow, int l1, int l2,
		       int l3, int l4, int l5)
{
  int type;
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l2));
      break;
    }
  case 2:
    {
      break;
    }
  case 3:
    {
      wfilters_content_validate(errCode,cstk(l2));
      if (strcmp(cstk(l3),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l4),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 4:
    {
      if (strcmp(cstk(l4),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l5),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * idwt2 validation 
 *-----------------------------------------*/

void
idwt2_form_validate (int *errCode, int *flow)
{
  if (sci_matrix_matrix_real(1) && sci_matrix_matrix_real(2) &&
      sci_matrix_matrix_real(3) && sci_matrix_matrix_real(4) &&
      matrix_length_check(1,2) && matrix_length_check(1,3) &&
      matrix_length_check(1,4))
    *errCode = SUCCESS;
  else if (sci_matrix_matrix_real(1) && sci_matrix_void(2) &&
	   sci_matrix_void(3) && sci_matrix_void(4))
    *errCode = SUCCESS;
  else if (sci_matrix_void(1) && sci_matrix_matrix_real(2) &&
	   sci_matrix_void(3) && sci_matrix_void(4))
    *errCode = SUCCESS;
  else if (sci_matrix_void(1) && sci_matrix_void(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_void(4))
    *errCode = SUCCESS;
  else if (sci_matrix_void(1) && sci_matrix_void(2) &&
	   sci_matrix_void(3) && sci_matrix_matrix_real(4))
    *errCode = SUCCESS;
  else if (sci_matrix_matrix_real(1) && sci_matrix_matrix_real(2) &&
	   sci_matrix_void(3) && sci_matrix_void(4) && 
	   matrix_length_check(1,2))
    *errCode = SUCCESS;
  else if (sci_matrix_matrix_real(1) && sci_matrix_void(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_void(4) &&
	   matrix_length_check(1,3))
    *errCode = SUCCESS;
  else if (sci_matrix_matrix_real(1) && sci_matrix_void(2) &&
	   sci_matrix_void(3) && sci_matrix_matrix_real(4) &&
	   matrix_length_check(1,4))
    *errCode = SUCCESS;
  else if (sci_matrix_void(1) && sci_matrix_matrix_real(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_void(4) &&
	   matrix_length_check(2,3))
    *errCode = SUCCESS;
  else if (sci_matrix_void(1) && sci_matrix_matrix_real(2) &&
	   sci_matrix_void(3) && sci_matrix_matrix_real(4) &&
	   matrix_length_check(2,4))
    *errCode = SUCCESS;
  else if (sci_matrix_void(1) && sci_matrix_void(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_matrix_real(4) &&
	   matrix_length_check(3,4))
    *errCode = SUCCESS;
  else if (sci_matrix_matrix_real(1) && sci_matrix_matrix_real(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_void(4) &&
	   matrix_length_check(1,2) && matrix_length_check(1,3))
    *errCode = SUCCESS;
  else if (sci_matrix_matrix_real(1) && sci_matrix_matrix_real(2) &&
	   sci_matrix_void(3) && sci_matrix_matrix_real(4) &&
	   matrix_length_check(1,2) && matrix_length_check(1,4))
    *errCode = SUCCESS;
  else if (sci_matrix_matrix_real(1) && sci_matrix_void(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_matrix_real(4) &&
	   matrix_length_check(1,3) && matrix_length_check(1,4))
    *errCode = SUCCESS;
  else if (sci_matrix_void(1) && sci_matrix_matrix_real(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_matrix_real(4) &&
	   matrix_length_check(2,3) && matrix_length_check(2,4))
    *errCode = SUCCESS;
  else
    {
      *errCode = UNKNOWN_INPUT_ERR;
      return;
    }
  
  if ((Rhs==5) && sci_strings_scalar(5))
    *flow = 1;
  else if ((Rhs==6) && sci_matrix_vector_real(5) && 
	   sci_matrix_vector_real(6) && vector_length_check(5,6))
    *flow = 2;
  else if ((Rhs==6) && sci_strings_scalar(5) && 
	   sci_matrix_vector_real(6) && length_check(6,2))
    *flow = 3;
  else if ((Rhs==7) && sci_matrix_vector_real(5) && 
	   sci_matrix_vector_real(6) && vector_length_check(5,6) && 
	   sci_matrix_vector_real(7) && length_check(7,2))
    *flow = 4;
  else if ((Rhs==7) && sci_strings_scalar(5) && 
	   sci_strings_scalar(6) && sci_strings_scalar(7))
    *flow = 5;
  else if ((Rhs==8) && sci_matrix_vector_real(5) && 
	   sci_matrix_vector_real(6) && vector_length_check(5,6) && 
	   sci_strings_scalar(7) && sci_strings_scalar(8))
    *flow = 6;
  else if ((Rhs==8) && sci_strings_scalar(5) && 
	   sci_matrix_vector_real(6) && length_check(6,2) && 
	   sci_strings_scalar(7) && sci_strings_scalar(8))
    *flow = 7;
  else if ((Rhs==9) && sci_matrix_vector_real(5) && 
	   sci_matrix_vector_real(6) && vector_length_check(5,6) && 
	   sci_matrix_vector_real(7) && length_check(7,2) && 
	   sci_strings_scalar(8) && sci_strings_scalar(9))
    *flow = 8;
  else
    *errCode = UNKNOWN_INPUT_ERR;

  return;
}

/*-------------------------------------------
 * idwt2 validation 
 *-----------------------------------------*/

void
idwt2_content_validate (int *errCode, int flow, int l1, int l2,
			int l3, int l4, int l5, int l6, int l7, 
			int l8, int l9)
{
  int type;
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l5));
      break;
    }
  case 2:
    {
      break;
    }
  case 3:
    {
      wfilters_content_validate(errCode,cstk(l5));
      if ((istk(l6)[0]<=0) || (istk(l6)[1]<=0))
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 4:
    {
      if ((istk(l7)[0]<=0) || (istk(l7)[1]<=0))
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 5:
    {
      wfilters_content_validate(errCode,cstk(l5));
      if (strcmp(cstk(l6),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l7),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 6:
    {
      if (strcmp(cstk(l7),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l8),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 7:
    {
      wfilters_content_validate(errCode,cstk(l5));
      if (strcmp(cstk(l7),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l8),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      if ((istk(l6)[0]<=0) || (istk(l6)[1]<=0))
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 8:
    {
      if (strcmp(cstk(l8),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l9),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      if ((istk(l7)[0]<=0) || (istk(l7)[1]<=0))
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * wavedec2 validation 
 *-----------------------------------------*/

void
wavedec2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_matrix_matrix_real(1) &&
      sci_matrix_scalar_real(2) && sci_strings_scalar(3))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_matrix_real(1) &&
	   sci_matrix_scalar_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4))
    *flow = 2;
  else 
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
wavedec2_content_validate (int *errCode, int flow, int l1, int l2,
			   int l3, int l4)
{
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      if (istk(l2)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 2:
    {
      if (istk(l2)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * waverec2 validation 
 *-----------------------------------------*/

void
waverec2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_matrix_vector_real(1) &&
      sci_matrix_matrix_real(2) && sci_strings_scalar(3) &&
      matrix_col_length_check(2,2))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_matrix_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && matrix_length_check(3,4) &&
	   matrix_col_length_check(2,2))
    *flow = 2;
  else
    *errCode = UNKNOWN_INPUT_ERR; 
  return;
}

void
waverec2_content_validate (int *errCode, int flow, int l1, int l2,
			   int l3, int l4)
{
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 2:
    {
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * wenergy2 validation 
 *-----------------------------------------*/

void
wenergy2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==2) && (Lhs==4) && sci_matrix_vector_real(1) &&
      sci_matrix_matrix_real(2) && matrix_col_length_check(2,2))
    *flow = 1;
  else if ((Rhs==2) && (Lhs==2) && sci_matrix_vector_real(1) &&
	   sci_matrix_matrix_real(2) && matrix_col_length_check(2,2))
    *flow = 2;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

/*-------------------------------------------
 * detcoef2 validation 
 *-----------------------------------------*/

void
detcoef2_form_validate (int *errCode, int *flow)
{
  if ((Rhs==4) && sci_strings_scalar(1) && 
      sci_matrix_vector_real(2) && sci_matrix_matrix_real(3) &&
      sci_matrix_scalar_real(4) && matrix_col_length_check(3,2))
    *errCode = SUCCESS;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}


void
detcoef2_content_validate (int *errCode, int flow, int l1, int l2,
			   int l3, int l4)
{
  if ((!strcmp(cstk(l1),"a")) || (!strcmp(cstk(l1),"h")) ||
      (!strcmp(cstk(l1),"v")) || (!strcmp(cstk(l1),"d")) ||
      (!strcmp(cstk(l1),"c")) || (!strcmp(cstk(l1),"all")) ||
      (!strcmp(cstk(l1),"compact")))
    *errCode = SUCCESS;
  else 
    *errCode = OPT_CHAR_NOT_VALID;
  return;
}


/*-------------------------------------------
 * appcoef2 validation 
 *-----------------------------------------*/

void
appcoef2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_matrix_vector_real(1) && 
      sci_matrix_matrix_real(2) && sci_strings_scalar(3) &&
      matrix_col_length_check(2,2))
    *flow = 2;
  else if ((Rhs==4) && sci_matrix_vector_real(1) && 
	   sci_matrix_matrix_real(2) && sci_strings_scalar(3) &&
	   sci_matrix_scalar_real(4) && matrix_col_length_check(2,2))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_matrix_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && 
	   matrix_col_length_check(2,2) && vector_length_check(3,4))
    *flow = 3;
  else if ((Rhs==5) && sci_matrix_vector_real(1) &&
	   sci_matrix_matrix_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && sci_matrix_scalar_real(5) &&
	   vector_length_check(3,4) && matrix_col_length_check(2,2))
    *flow = 4;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}


void
appcoef2_content_validate (int *errCode, int flow, int l1, int l2,
			   int l3, int l4, int l5)
{
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      if (istk(l4)[0]<0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 2:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 3:
    {
      break;
    }
  case 4:
    {
      if (istk(l5)[0]<0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  default:
    break;
  }
  return;
}


/*-------------------------------------------
 * wrcoef2 validation 
 *-----------------------------------------*/

void
wrcoef2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==4) && sci_strings_scalar(1) && 
      sci_matrix_vector_real(2) && sci_matrix_matrix_real(3) &&
      sci_strings_scalar(4) && matrix_col_length_check(3,2))
    *flow = 3;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_matrix_real(3) &&
	   sci_strings_scalar(4) && sci_matrix_scalar_real(5) &&
	   matrix_col_length_check(3,2))
    *flow = 1;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_matrix_real(3) &&
	   sci_matrix_vector_real(4) && sci_matrix_vector_real(5) &&
	   vector_length_check(4,5) && matrix_col_length_check(3,2))
    *flow = 4;
  else if ((Rhs==6) && sci_strings_scalar(1) &&
	   sci_matrix_vector_real(2) && sci_matrix_matrix_real(3) &&
	   sci_matrix_vector_real(4) && sci_matrix_vector_real(5) &&
	   vector_length_check(4,5) && 
	   matrix_col_length_check(3,2) && sci_matrix_scalar_real(6))
    *flow = 2;
  return;
}

void
wrcoef2_content_validate (int *errCode, int flow, int l1, int l2,
			  int l3, int l4, int l5, int l6)
{
  if (scalar_string_check(cstk(l1),'a') || 
      scalar_string_check(cstk(l1),'h') ||
      scalar_string_check(cstk(l1),'v') ||
      scalar_string_check(cstk(l1),'d'))
    *errCode = SUCCESS;
  else
    {
      *errCode = OPT_CHAR_NOT_VALID;
      return;
    }
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l4));
      if (istk(l5)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      if (scalar_string_check(cstk(l1),'a') || 
	  scalar_string_check(cstk(l1),'h') ||
	  scalar_string_check(cstk(l1),'v') ||
	  scalar_string_check(cstk(l1),'d'))
	*errCode = SUCCESS;
      else
	*errCode = OPT_CHAR_NOT_VALID;
      break;
    }
  case 2:
    {
      if (istk(l6)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 3:
    {
      wfilters_content_validate(errCode,cstk(l4));
      break;
    }
  case 4:
    {
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * upwlev2 validation 
 *-----------------------------------------*/

void
upwlev2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_matrix_vector_real(1) && 
      sci_matrix_matrix_real(2) && sci_strings_scalar(3) &&
      matrix_col_length_check(2,2))
    *flow = 1;
  else if ((Rhs==4) && sci_matrix_vector_real(1) &&
	   sci_matrix_matrix_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4) &&
	   matrix_col_length_check(2,2))
    *flow = 2;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
upwlev2_content_validate (int *errCode, int flow, int l1, int l2,
			  int l3, int l4)
{
  *errCode = SUCCESS;
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 2:
    {
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * upcoef2 validation 
 *-----------------------------------------*/

void
upcoef2_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==3) && sci_strings_scalar(1) && 
      sci_matrix_matrix_real(2) && sci_strings_scalar(3))
    *flow = 5;
  else if ((Rhs==4) && sci_strings_scalar(1) &&
	   sci_matrix_matrix_real(2) && sci_strings_scalar(3) &&
	   sci_matrix_scalar_real(4))
    *flow = 3;
  else if ((Rhs==4) && sci_strings_scalar(1) &&
	   sci_matrix_matrix_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4))
    *flow = 6;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   sci_matrix_matrix_real(2) && sci_strings_scalar(3) &&
	   sci_matrix_scalar_real(4) && sci_matrix_vector_real(5) &&
	   length_check(5,2))
    *flow = 1;
  else if ((Rhs==5) && sci_strings_scalar(1) &&
	   sci_matrix_matrix_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4) &&
	   sci_matrix_scalar_real(5))
    *flow = 4;
  else if ((Rhs==6) && sci_strings_scalar(1) &&
	   sci_matrix_matrix_real(2) && sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) && vector_length_check(3,4) &&
	   sci_matrix_scalar_real(5) && sci_matrix_vector_real(6) &&
	   length_check(6,2))
    *flow = 2;
  else 
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
upcoef2_content_validate (int *errCode, int flow, int l1, int l2,
			  int l3, int l4, int l5, int l6)
{
  if ((!strcmp(cstk(l1),"a")) || (!strcmp(cstk(l1),"h")) ||
      (!strcmp(cstk(l1),"v")) || (!strcmp(cstk(l1),"d")))
    *errCode = SUCCESS;
  else
    {
      *errCode = OPT_CHAR_NOT_VALID;
      return;
    }
  switch (flow) {
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l3));
      if ((istk(l4)[0]<=0) || (istk(l5)[0]<=0) ||
	  (istk(l5)[1]<=0))
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 2:
    {
      if ((istk(l5)[0]<=0) || (istk(l6)[0]<=0) ||
	  (istk(l6)[1]<=0))
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 3:
    {
      wfilters_content_validate(errCode,cstk(l3));
      if (istk(l4)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 4:
    {
      if (istk(l5)[0]<=0)
	*errCode = POSITIVE_INTEGER_ONLY;
      break;
    }
  case 5:
    {
      wfilters_content_validate(errCode,cstk(l3));
      break;
    }
  case 6:
    {
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * dwt3 validation 
 *-----------------------------------------*/
void
dwt3_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==2) && sci_mlist_check(1) && 
      sci_strings_scalar(2))
    *flow = 1;
  else if ((Rhs==3) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) && 
	   sci_matrix_vector_real(3) &&
	   vector_length_check(2,3))
    *flow = 2;
  else if ((Rhs==4) && sci_mlist_check(1) &&
	   sci_strings_scalar(2) &&
	   sci_strings_scalar(3) &&
	   sci_strings_scalar(4))
    *flow = 3;
  else if ((Rhs==5) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) &&
           sci_matrix_vector_real(3) &&
	   vector_length_check(2,3) &&
	   sci_strings_scalar(4) &&
	   sci_strings_scalar(5))
    *flow = 4;
  else if ((Rhs==4) && sci_mlist_check(1) &&
	   sci_strings_scalar(2) &&
	   sci_strings_scalar(3) &&
	   sci_strings_scalar(4))
    *flow = 5;
  else if ((Rhs==7) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) &&
	   sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) &&
	   sci_matrix_vector_real(5) &&
	   sci_matrix_vector_real(6) &&
	   sci_matrix_vector_real(7) &&
	   vector_length_check(2,3) &&
	   vector_length_check(4,5) &&
	   vector_length_check(6,7))
    *flow = 6;
  else if ((Rhs==6) && sci_mlist_check(1) &&
	   sci_strings_scalar(2) &&
	   sci_strings_scalar(3) &&
	   sci_strings_scalar(4) &&
	   sci_strings_scalar(5) &&
	   sci_strings_scalar(6) )
    *flow = 7;
  else if ((Rhs==9) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) &&
	   sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) &&
	   sci_matrix_vector_real(5) &&
	   sci_matrix_vector_real(6) &&
	   sci_matrix_vector_real(7) &&
	   vector_length_check(2,3) &&
	   vector_length_check(4,5) &&
	   vector_length_check(6,7) &&
	   sci_strings_scalar(8) &&
	   sci_strings_scalar(9) )
    *flow = 8;
  else 
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
dwt3_content_validate (int *errCode, int flow, int l1, int l2,
		       int l3, int l4, int l5, int l6, int l7,
		       int l8, int l9)
{
  int type, err1, err2, err3;
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l2));
      break;
    }
  case 2:
    {
      break;
    }
  case 3:
    {
      wfilters_content_validate(errCode,cstk(l2));
      if (strcmp(cstk(l3),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l4),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 4:
    {
      if (strcmp(cstk(l4),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l5),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 5:
    {
      wfilters_content_validate(&err1,cstk(l2));
      wfilters_content_validate(&err2,cstk(l3));
      wfilters_content_validate(&err3,cstk(l4));
      if ((err1 != SUCCESS) || (err1 != SUCCESS) || 
	  (err1 != SUCCESS))
	*errCode = WAVELET_NAME_NOT_VALID;
      break;
    }
  case 6:
    {
      
      break;
    }
  case 7:
    {
      wfilters_content_validate(&err1,cstk(l2));
      wfilters_content_validate(&err2,cstk(l3));
      wfilters_content_validate(&err3,cstk(l4));
      if ((err1 != SUCCESS) || (err1 != SUCCESS) || 
	  (err1 != SUCCESS))
	*errCode = WAVELET_NAME_NOT_VALID;
      if (strcmp(cstk(l5),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l6),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  case 8:
    {
      if (strcmp(cstk(l8),"mode"))
	*errCode = UNKNOWN_INPUT_ERR;
      extension_check(cstk(l9),&type);
      if (!type)
	*errCode = EXTENSION_OPT_NOT_VALID;
      break;
    }
  default:
    break;
  }
  return;
}


/*-------------------------------------------
 * idwt3 validation 
 *-----------------------------------------*/
void
idwt3_form_validate (int *errCode, int *flow)
{
  *errCode = SUCCESS;
  if ((Rhs==2) && sci_mlist_check(1) && 
      sci_strings_scalar(2))
    *flow = 1;
  else if ((Rhs==3) && sci_mlist_check(1) &&
	   sci_strings_scalar(2) && 
	   sci_matrix_vector_real(3) &&
	   length_check(3,3))
    *flow = 2;
  else if ((Rhs==3) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) &&
	   sci_matrix_vector_real(3) &&
	   vector_length_check(2,3))
    *flow = 3;
  else if ((Rhs==4) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) &&
           sci_matrix_vector_real(3) &&
	   vector_length_check(2,3) &&
	   sci_matrix_vector_real(4) &&
	   length_check(4,3))
    *flow = 4;
  else if ((Rhs==4) && sci_mlist_check(1) &&
	   sci_strings_scalar(2) &&
	   sci_strings_scalar(3) &&
	   sci_strings_scalar(4))
    *flow = 5;
  else if ((Rhs==5) && sci_mlist_check(1) &&
	   sci_strings_scalar(2) &&
	   sci_strings_scalar(3) &&
	   sci_strings_scalar(4) &&
	   sci_matrix_vector_real(5) &&
	   length_check(5,3))
    *flow = 6;
  else if ((Rhs==7) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) &&
	   sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) &&
	   sci_matrix_vector_real(5) &&
	   sci_matrix_vector_real(6) &&
	   sci_matrix_vector_real(7) &&
	   vector_length_check(2,3) &&
	   vector_length_check(4,5) &&
	   vector_length_check(6,7))
    *flow = 7;
  else if ((Rhs==8) && sci_mlist_check(1) &&
	   sci_matrix_vector_real(2) &&
	   sci_matrix_vector_real(3) &&
	   sci_matrix_vector_real(4) &&
	   sci_matrix_vector_real(5) &&
	   sci_matrix_vector_real(6) &&
	   sci_matrix_vector_real(7) &&
	   vector_length_check(2,3) &&
	   vector_length_check(4,5) &&
	   vector_length_check(6,7) &&
	   sci_matrix_vector_real(8) &&
	   length_check(8,3))
    *flow = 8;
  else 
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
idwt3_content_validate (int *errCode, int flow, int l1, 
			int l2,	int l3, int l4, int l5, 
			int l6, int l7, int l8)
{
  int type, err1, err2, err3;
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      wfilters_content_validate(errCode,cstk(l2));
      break;
    }
  case 2:
    {
      wfilters_content_validate(errCode,cstk(l2));
      if ((istk(l3)[0]<=0) || (istk(l3)[1]<=0) || (istk(l3)[2]<=0))
	{
	  *errCode = UNKNOWN_INPUT_ERR;
	}
      break;
    }
  case 3:
    {
      break;
    }
  case 4:
    {
      if ((istk(l4)[0]<=0) || (istk(l4)[1]<=0) || (istk(l4)[2]<=0))
	{
	  *errCode = UNKNOWN_INPUT_ERR;
	}
      break;
    }
  case 5:
    {
      wfilters_content_validate(&err1,cstk(l2));
      wfilters_content_validate(&err2,cstk(l3));
      wfilters_content_validate(&err3,cstk(l4));
      if ((err1 != SUCCESS) || (err1 != SUCCESS) || 
	  (err1 != SUCCESS))
	*errCode = WAVELET_NAME_NOT_VALID;
      break;
    }
  case 6:
    {
      wfilters_content_validate(&err1,cstk(l2));
      wfilters_content_validate(&err2,cstk(l3));
      wfilters_content_validate(&err3,cstk(l4));
      if ((err1 != SUCCESS) || (err1 != SUCCESS) || 
	  (err1 != SUCCESS))
	*errCode = WAVELET_NAME_NOT_VALID;
      if ((istk(l5)[0]<=0) || (istk(l5)[1]<=0) || (istk(l5)[2]<=0))
	{
	  *errCode = UNKNOWN_INPUT_ERR;
	}
      break;
    }
  case 7:
    {
      break;
    }
  case 8:
    {
      if ((istk(l8)[0]<=0) || (istk(l8)[1]<=0) || (istk(l8)[2]<=0))
	{
	  *errCode = UNKNOWN_INPUT_ERR;
	}
      break;
    }
  default:
    break;
  }
  return;
}
