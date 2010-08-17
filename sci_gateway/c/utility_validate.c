/*
 * -------------------------------------------------------------------------
 * utility_validate.c -- Utility Function Validation
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
 * convolution validation 
 *-----------------------------------------*/

void
conv_validate (int *errCode)
{
  if ((sci_matrix_vector_real(1) || sci_matrix_scalar_real(1)) && 
      (sci_matrix_vector_real(2) || sci_matrix_scalar_real(2)))
    *errCode = SUCCESS;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

/*-------------------------------------------
 * wrev validation 
 *-----------------------------------------*/

void
wrev_validate (int *errCode)
{
  if (sci_matrix_vector_real(1))
    *errCode = SUCCESS;
  else
    *errCode = UNKNOWN_INPUT_ERR;
   return;
}

/*-------------------------------------------
 * qmf validation 
 *-----------------------------------------*/

void
qmf_validate (int *flow, int *errCode)
{
  *errCode = SUCCESS;
  if ((Rhs == 1) && (sci_matrix_vector_real(1)))
    *flow = 1;
  else if ((Rhs == 2) && (sci_matrix_vector_real(1)) && 
		       (sci_matrix_scalar_real(2)))
    *flow = 2;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}


/*-------------------------------------------
 * dyaddown validation 
 *-----------------------------------------*/

void
dyaddown_form_validate (int *flow, int *errCode)
{
  *errCode = SUCCESS;
  if ((Rhs == 1) && (sci_matrix_vector_real(1)))
    *flow = 1;
  else if ((Rhs == 1) && (sci_matrix_matrix_real(1)))
    *flow = 3;
  else if ((Rhs == 2) && (sci_matrix_vector_real(1)) &&
	   (sci_matrix_scalar_real(2)))
    *flow = 2;
  else if ((Rhs == 2) && (sci_matrix_matrix_real(1)) &&
	   (sci_matrix_scalar_real(2)))
    *flow = 5;
  else if ((Rhs == 2) && (sci_matrix_matrix_real(1)) &&
	   (sci_strings_scalar(2)))
    *flow = 4;
  else if ((Rhs == 3) && (sci_matrix_matrix_real(1)) &&
	   (sci_matrix_scalar_real(2)) && (sci_strings_scalar(3)))
    *flow = 6;
  else if ((Rhs == 3) && (sci_matrix_matrix_real(1)) &&
	   (sci_matrix_scalar_real(3)) && (sci_strings_scalar(2)))
    *flow = 7;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
dyaddown_content_validate (char *l, int *errCode)
{
  if (scalar_string_check(l,'r') ||
       scalar_string_check(l,'c') ||
       scalar_string_check(l,'m'))
    *errCode = SUCCESS;
  else
    *errCode = OPT_CHAR_NOT_VALID;
  return;
}

/*-------------------------------------------
 * dyadup validation 
 *-----------------------------------------*/

void
dyadup_form_validate (int *flow, int *errCode)
{
  *errCode = SUCCESS;
  if ((Rhs == 1) && (sci_matrix_vector_real(1)))
    *flow = 1;
  else if ((Rhs == 1) && (sci_matrix_matrix_real(1)))
    *flow = 3;
  else if ((Rhs == 2) && (sci_matrix_vector_real(1)) &&
	   (sci_matrix_scalar_real(2)))
    *flow = 2;
  else if ((Rhs == 2) && (sci_matrix_matrix_real(1)) &&
	   (sci_matrix_scalar_real(2)))
    *flow = 5;
  else if ((Rhs == 2) && (sci_matrix_matrix_real(1)) &&
	   (sci_strings_scalar(2)))
    *flow = 4;
  else if ((Rhs == 3) && (sci_matrix_matrix_real(1)) &&
	   (sci_matrix_scalar_real(2)) && (sci_strings_scalar(3)))
    *flow = 6;
  else if ((Rhs == 3) && (sci_matrix_matrix_real(1)) &&
	   (sci_matrix_scalar_real(3)) && (sci_strings_scalar(2)))
    *flow = 7;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
dyadup_content_validate (char *l, int *errCode)
{
  if (scalar_string_check(l,'r') ||
       scalar_string_check(l,'c') ||
       scalar_string_check(l,'m'))
    *errCode = SUCCESS;
  else
    *errCode = OPT_CHAR_NOT_VALID;
  return;
}


/*-------------------------------------------
 * wkeep validation 
 *-----------------------------------------*/

void wkeep_form_validate (int *flow, int *errCode)
{
  *errCode = SUCCESS;
  if ((Rhs==2) && sci_matrix_vector_real(1) && 
      sci_matrix_scalar_real(2))
    *flow = 1;
  else if ((Rhs==2) && sci_matrix_matrix_real(1) && 
	   sci_matrix_vector_real(2) && length_check(2,2))
    *flow = 2;
  else if ((Rhs==3) && sci_matrix_vector_real(1) &&
	   sci_matrix_scalar_real(2) && sci_strings_scalar(3))
    *flow = 3;
  else if ((Rhs==3) && sci_matrix_vector_real(1) &&
	   sci_matrix_scalar_real(2) && sci_matrix_scalar_real(3))
    *flow = 4;
  else if ((Rhs==3) && sci_matrix_matrix_real(1) &&
	   sci_matrix_vector_real(2) && length_check(2,2) &&
	   sci_matrix_vector_real(3) && length_check(3,2))
    *flow = 5;
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
wkeep_content_validate (int flow, int *errCode, 
			int l1, int l2, int l3)
{
  int res, resRow, resCol;
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      if (istk(l2)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(1, istk(l2)[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      break;
    }
  case 2:
    {
      if ((istk(l2)[0] <= 0) || (istk(l2)[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(1, istk(l2)[0], istk(l2)[1], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	*errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
      break;
    }
  case 3:
    {
      if (istk(l2)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(1, istk(l2)[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      if ((cstk(l3)[0]!='r') && (cstk(l3)[0]!='l') &&
	  (cstk(l3)[0]!='c'))
	*errCode = OPT_CHAR_NOT_VALID;
      break;
    }
  case 4:
    {
      if (istk(l2)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      if (istk(l3)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(1, istk(l2)[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      vector_length_compare(1, istk(l3)[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      vector_length_compare(1, istk(l2)[0]+istk(l3)[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      break;
    }
  case 5:
    {
      if ((istk(l2)[0] <= 0) || (istk(l2)[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      if ((istk(l3)[0] <= 0) || (istk(l3)[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(1, istk(l2)[0], istk(l2)[1], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	*errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
      matrix_length_compare(1, istk(l3)[0], istk(l3)[1], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	*errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
      matrix_length_compare(1, istk(l2)[0]+istk(l3)[0], 
			    istk(l2)[1]+istk(l3)[1], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	*errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
      break;
    }
  default:
    break;
  }
  return;
}


/*-------------------------------------------
 * wextend validation 
 *-----------------------------------------*/

void 
wextend_form_validate (int *flow, int *errCode, int l1)
{
  *errCode = SUCCESS;
  if ((Rhs==4) && sci_strings_scalar(2) &&
      sci_matrix_vector_real(3) && sci_matrix_scalar_real(4))
    {
      if (sci_matrix_scalar_real(1))
	  {
	     if (istk(l1)[0]==1)
	        *flow = 2;
	     else
	        *errCode = UNKNOWN_INPUT_ERR;
	   }
      else if (sci_strings_scalar(1))
	  {
	     if ((!strcmp(cstk(l1),"1"))  ||
	         (!strcmp(cstk(l1),"1d")) ||
	         (!strcmp(cstk(l1),"1D")))
	        *flow = 2;
	     else
	        *errCode = UNKNOWN_INPUT_ERR;
	   }
      else
	    *errCode = UNKNOWN_INPUT_ERR;
      }
  else if ((Rhs==4) && sci_strings_scalar(2) &&
	   sci_matrix_matrix_real(3) && sci_matrix_scalar_real(4))
     {
        if (sci_matrix_scalar_real(1))
	      {
	         if (istk(l1)[0]==2)
	           *flow = 4;
	         else
	           *errCode = UNKNOWN_INPUT_ERR;
	       }
         else if (sci_strings_scalar(1))
	      {
	         if ((!strcmp(cstk(l1),"2"))  ||
	             (!strcmp(cstk(l1),"2d")) ||
	             (!strcmp(cstk(l1),"2D")))
	           *flow = 4;
	         else if ((!strcmp(cstk(l1),"ar")) ||
	                  (!strcmp(cstk(l1),"addrow")))
	           *flow = 7;
             else if ((!strcmp(cstk(l1),"ac")) ||
	                  (!strcmp(cstk(l1),"addcol")))
	           *flow = 9;
	  	     else
	           *errCode = UNKNOWN_INPUT_ERR;
	       }
          else
	           *errCode = UNKNOWN_INPUT_ERR;
       }
  else if ((Rhs==4) && sci_strings_scalar(2) &&
	   sci_matrix_matrix_real(3) &&
	   sci_matrix_vector_real(4) && length_check(4,2))
    {
      if (sci_matrix_scalar_real(1))
	{
	  if (istk(l1)[0]==2)
	    *flow = 5;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp(cstk(l1),"2"))  ||
	      (!strcmp(cstk(l1),"2d")) ||
	      (!strcmp(cstk(l1),"2D")))
	    *flow = 5;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else
	*errCode = UNKNOWN_INPUT_ERR;
    }
  else if ((Rhs==4) && sci_strings_scalar(2) &&
	   sci_matrix_matrix_real(3) &&
	   sci_matrix_scalar_real(4) && sci_strings_scalar(1))
    {
		if (sci_strings_scalar(1))
		{
           if ((!strcmp(cstk(l1),"ar")) ||
	           (!strcmp(cstk(l1),"addrow")))
	         *flow = 7;
           else if ((!strcmp(cstk(l1),"ac")) ||
	         (!strcmp(cstk(l1),"addcol")))
	         *flow = 9;
           else
	         *errCode = UNKNOWN_INPUT_ERR;
		}
		else
			*errCode = UNKNOWN_INPUT_ERR;
    }
  else if ((Rhs==5) && sci_strings_scalar(2) &&
	   sci_matrix_vector_real(3) && 
	   sci_matrix_scalar_real(4) && sci_strings_scalar(5))
    {
      if (sci_matrix_scalar_real(1))
	{
	  if (istk(l1)[0]==1)
	    *flow = 1;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp(cstk(l1),"1"))  ||
	      (!strcmp(cstk(l1),"1d")) ||
	      (!strcmp(cstk(l1),"1D")))
	    *flow = 1;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else
	*errCode = UNKNOWN_INPUT_ERR;
    }
  else if ((Rhs==5) && sci_strings_scalar(2) &&
	   sci_matrix_matrix_real(3) &&
	   sci_matrix_vector_real(4) && length_check(4,2) &&
	   sci_strings_vector(5) && length_check(5,2))
    {
      if (sci_matrix_scalar_real(1))
	{
	  if (istk(l1)[0]==2)
	    *flow = 3;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp(cstk(l1),"2"))  ||
	      (!strcmp(cstk(l1),"2d")) ||
	      (!strcmp(cstk(l1),"2D")))
	    *flow = 3;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else
	*errCode = UNKNOWN_INPUT_ERR;
    }
  else if ((Rhs==5) && sci_strings_scalar(2) &&
	   sci_matrix_matrix_real(3) &&
	   sci_matrix_vector_real(4) && length_check(4,2) &&
	   sci_strings_scalar(5))
    {
      if (sci_matrix_scalar_real(1))
	{
	  if (istk(l1)[0]==2)
	    *flow = 10;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp(cstk(l1),"2"))  ||
	      (!strcmp(cstk(l1),"2d")) ||
	      (!strcmp(cstk(l1),"2D")))
	    *flow = 10;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else
	*errCode = UNKNOWN_INPUT_ERR;
    }
  else if ((Rhs==5) && sci_strings_scalar(2) &&
	   sci_matrix_matrix_real(3) && 
	   sci_matrix_scalar_real(4) &&
	   sci_strings_scalar(5))
    {
	  //sciprint("enter Rhs==5\n");
      if (sci_strings_scalar(1))
	  {
          if ((!strcmp(cstk(l1),"ar")) ||
	          (!strcmp(cstk(l1),"addrow")))
	         *flow = 6;
          else if ((!strcmp(cstk(l1),"ac")) ||
	          (!strcmp(cstk(l1),"addcol")))
	         *flow = 8;
          else
	         *errCode = UNKNOWN_INPUT_ERR;
       }
	  else
		  *errCode = UNKNOWN_INPUT_ERR;
     }
  else
    *errCode = UNKNOWN_INPUT_ERR;
  return;
}

void
wextend_content_validate (int flow, int *errCode, int l2,
			  int l3, int l4, int l5, char **str)
{
  int type, res, resRow, resCol;
  *errCode = SUCCESS;
  extension_check(cstk(l2),&type);
  if (!type)
    {
      *errCode = EXTENSION_OPT_NOT_VALID;
      return;
    }
  switch (flow) {
  case 1:
    {
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(3, istk(l4)[0], &res);
      if (res == -1)
	{
	  *errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
	  return;
	}
      if ((cstk(l5)[0]!='b') && (cstk(l5)[0]!='l') &&
	  (cstk(l5)[0]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  case 2:
    {
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(3, istk(l4)[0], &res);
      if (res == -1)
	{
	  *errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
	  return;
	}
      break;
    }
  case 3:
    {
      if ((istk(l4)[0] <= 0) || (istk(l4)[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[1], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      if ((str[0][0]!='b') && (str[0][0]!='l') &&
	  (str[0][0]!='r') && (str[1][0]!='b') &&
	  (str[1][0]!='l') && (str[1][0]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  case 4:
    {
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[0], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      break;
    }
  case 5:
    {
      if ((istk(l4)[0] <= 0) || (istk(l4)[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[1], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      break;
    }
  case 6:
    {
      if (istk(l4)[0] <= 0)
	   *errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[0], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      if ((cstk(l5)[0]!='b') && (cstk(l5)[0]!='l') &&
	  (cstk(l5)[0]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  case 7:
    {
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[0], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      break;
    }
  case 8:
    {
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[0], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      if ((cstk(l5)[0]!='b') && (cstk(l5)[0]!='l') &&
	  (cstk(l5)[0]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  case 9:
    {
      if (istk(l4)[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[0], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      break;
    }
  case 10:
    {
      if ((istk(l4)[0] <= 0) || (istk(l4)[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, istk(l4)[0], istk(l4)[1], 
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      if ((cstk(l5)[0]!='b') && (cstk(l5)[0]!='l') &&
	  (cstk(l5)[0]!='r') && (cstk(l5)[1]!='b') &&
	  (cstk(l5)[1]!='l') && (cstk(l5)[1]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  default:
    break;
  }
  return;
}

/*-------------------------------------------
 * wcodemat validation 
 *-----------------------------------------*/
void wcodemat_form_validate (int *flow, int *errCode)
{
    *errCode = SUCCESS;
    if ((Rhs==1) && (sci_matrix_matrix_real(1)))
		*flow = 1;
	else if ((Rhs==2) && (sci_matrix_matrix_real(1)) && sci_matrix_scalar_real(2))
		*flow = 2;
	else if ((Rhs==3) && sci_matrix_matrix_real(1) && sci_matrix_scalar_real(2) && sci_strings_scalar(3)) 
		*flow = 3;
	else if ((Rhs==4) && sci_matrix_matrix_real(1) && sci_matrix_scalar_real(2) && sci_strings_scalar(3) && 
		      sci_matrix_scalar_real(4))
		*flow = 4;
	else
		*errCode = UNKNOWN_INPUT_ERR;
	return;
}

void wcodemat_content_validate (int *errCode, int flow, int l1, int l2, int l3, int l4)
{

	*errCode = SUCCESS;
	switch (flow){
		case 1:
			{
				break;
			}
		case 2:
			{
                if ((istk(l2)[0]<1) || (istk(l2)[0]>512))
					*errCode = UNKNOWN_INPUT_ERR;
				break;
			}
		case 3:
			{
				if ((istk(l2)[0]<1) || (istk(l2)[0]>512))
					*errCode = UNKNOWN_INPUT_ERR;
                
				break;
			}
		case 4:
			{
                if ((istk(l2)[0]<1) || (istk(l2)[0]>512))
					*errCode = UNKNOWN_INPUT_ERR;
              	break;
			}
		default:
			{
				*errCode  = UNKNOWN_INPUT_ERR;
				break;
			}
	}

	return;
}

/*-------------------------------------------
 * wrev2 validation 
 *-----------------------------------------*/
void wrev2_form_validate (int *errCode)
{
    *errCode = UNKNOWN_INPUT_ERR;
    if ((Rhs==2) && (sci_matrix_matrix_real(1)) &&
	(sci_matrix_scalar_real(2)))
      *errCode = SUCCESS;
    return;
}


void wnorm_form_validate (int *flow, int *errCode)
{
  *errCode = SUCCESS;
  if ((sci_matrix_matrix_real(1)) || (sci_matrix_vector_real(1)))
    {
      if (Rhs==1)
	*flow = 1;
      else if ((Rhs==3) && (sci_matrix_scalar_real(2)) &&
			 (sci_matrix_scalar_real(3)))
	*flow=2;
      else
	*errCode = UNKNOWN_INPUT_ERR;  
    }
  else
    *errCode = UNKNOWN_INPUT_ERR;
  
  return;
}


