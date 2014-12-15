/*
 * -------------------------------------------------------------------------
 * utility_validate.c -- Utility Function Validation
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
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
 //#include "swt_common.h"
 //#include "dwt.h"
 #include "swtlib.h"
 #include "swt_gwsupport.h"
 #include "api_scilab.h"
 #include "Scierror.h"
 //#include "localization.h"
 //#include "warningmode.h"
 #include "sciprint.h"

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
			double *input1, int *int_input2, int *int_input3)
{
  int res, resRow, resCol;
  *errCode = SUCCESS;
  switch (flow){
  case 1:
    {
      if (int_input2[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(1, int_input2[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      break;
    }
  case 2:
    {
      if ((int_input2[0] <= 0) || (int_input2[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(1, int_input2[0], int_input2[1],
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	*errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
      break;
    }
  case 4:
    {
      if (int_input2[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      if (int_input3[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(1, int_input2[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      vector_length_compare(1, int_input3[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      vector_length_compare(1, int_input2[0]+int_input3[0], &res);
      if (res == -1)
	*errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
      break;
    }
  case 5:
    {
      if ((int_input2[0] <= 0) || (int_input2[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      if ((int_input3[0] <= 0) || (int_input3[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(1, int_input2[0], int_input2[1],
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	*errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
      matrix_length_compare(1, int_input3[0], int_input3[1],
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	*errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
      matrix_length_compare(1, int_input2[0]+int_input3[0],
			    int_input2[1]+int_input3[1],
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

void
wkeep_content_validate_string (int flow, int *errCode,
double *input1, int *int_input2, char *input_string)
{
  int res, resRow, resCol;
  *errCode = SUCCESS;
  switch (flow){
            case 3:
            {
              if (int_input2[0] <= 0)
                *errCode = POSITIVE_INTEGER_ONLY;
                vector_length_compare(1, int_input2[0], &res);
                if (res == -1)
                  *errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
                  if ((input_string[0]!='r') && (input_string[0]!='l') &&
                    (input_string[0]!='c'))
                    *errCode = OPT_CHAR_NOT_VALID;
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
wextend_form_validate (int *flow, int *errCode, char* l1,int* int_l1)
{
  *errCode = SUCCESS;
  if ((Rhs==4) && sci_strings_scalar(2) &&
      sci_matrix_vector_real(3) && sci_matrix_scalar_real(4))
    {
      if (sci_matrix_scalar_real(1))
	  {
	     if (int_l1[0]==1)
	        *flow = 2;
	     else
	        *errCode = UNKNOWN_INPUT_ERR;
	   }
      else if (sci_strings_scalar(1))
	  {
	     if ((!strcmp((l1),"1"))  ||
	         (!strcmp((l1),"1d")) ||
	         (!strcmp((l1),"1D")))
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
	         if (int_l1[0]==2)
	           *flow = 4;
	         else
	           *errCode = UNKNOWN_INPUT_ERR;
	       }
         else if (sci_strings_scalar(1))
	      {
	         if ((!strcmp((l1),"2"))  ||
	             (!strcmp((l1),"2d")) ||
	             (!strcmp((l1),"2D")))
	           *flow = 4;
	         else if ((!strcmp((l1),"ar")) ||
	                  (!strcmp((l1),"addrow")))
	           *flow = 7;
             else if ((!strcmp((l1),"ac")) ||
	                  (!strcmp((l1),"addcol")))
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
	  if (int_l1[0]==2)
	    *flow = 5;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp((l1),"2"))  ||
	      (!strcmp((l1),"2d")) ||
	      (!strcmp((l1),"2D")))
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
           if ((!strcmp((l1),"ar")) ||
	           (!strcmp((l1),"addrow")))
	         *flow = 7;
           else if ((!strcmp((l1),"ac")) ||
	         (!strcmp((l1),"addcol")))
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
	  if (int_l1[0]==1)
	    *flow = 1;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp((l1),"1"))  ||
	      (!strcmp((l1),"1d")) ||
	      (!strcmp((l1),"1D")))
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
	  if (int_l1[0]==2)
	    *flow = 3;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp((l1),"2"))  ||
	      (!strcmp((l1),"2d")) ||
	      (!strcmp((l1),"2D")))
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
	  if (int_l1[0]==2)
	    *flow = 10;
	  else
	    *errCode = UNKNOWN_INPUT_ERR;
	}
      else if (sci_strings_scalar(1))
	{
	  if ((!strcmp((l1),"2"))  ||
	      (!strcmp((l1),"2d")) ||
	      (!strcmp((l1),"2D")))
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
          if ((!strcmp((l1),"ar")) ||
	          (!strcmp((l1),"addrow")))
	         *flow = 6;
          else if ((!strcmp((l1),"ac")) ||
	          (!strcmp((l1),"addcol")))
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
wextend_content_validate (int flow, int *errCode, char *input_string1,
			  double *input3, int* int_intput4, char *input_string2, char **str)
{
  int type, res, resRow, resCol;
  *errCode = SUCCESS;
  extension_check(input_string1,&type);
  if (!type)
    {
      *errCode = EXTENSION_OPT_NOT_VALID;
      return;
    }
  switch (flow) {
  case 1:
    {
      if (int_intput4[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(3, int_intput4[0], &res);
      if (res == -1)
	{
	  *errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
	  return;
	}
      if ((input_string2[0]!='b') && (input_string2[0]!='l') &&
	  (input_string2[0]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  case 2:
    {
      if (int_intput4[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      vector_length_compare(3, int_intput4[0], &res);
      if (res == -1)
	{
	  *errCode = LENGTH_DATA_NOT_VALID_FOR_VECTOR_DIMENSION;
	  return;
	}
      break;
    }
  case 3:
    {
      if ((int_intput4[0] <= 0) || (int_intput4[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[1],
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
      if (int_intput4[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[0],
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
      if ((int_intput4[0] <= 0) || (int_intput4[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[1],
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
      if (int_intput4[0] <= 0)
	   *errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[0],
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      if ((input_string2[0]!='b') && (input_string2[0]!='l') &&
	  (input_string2[0]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  case 7:
    {
      if (int_intput4[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[0],
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
      if (int_intput4[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[0],
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      if ((input_string2[0]!='b') && (input_string2[0]!='l') &&
	  (input_string2[0]!='r'))
	{
	  *errCode = OPT_CHAR_NOT_VALID;
	  return;
	}
      break;
    }
  case 9:
    {
      if (int_intput4[0] <= 0)
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[0],
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
      if ((int_intput4[0] <= 0) || (int_intput4[1] <= 0))
	*errCode = POSITIVE_INTEGER_ONLY;
      matrix_length_compare(3, int_intput4[0], int_intput4[1],
			    &resRow, &resCol);
      if ((resRow == -1 )||( resCol == -1))
	{
	  *errCode = SIZE_DATA_NOT_VALID_FOR_MATRIX_DIMENSION;
	  return;
	}
      if ((input_string2[0]!='b') && (input_string2[0]!='l') &&
	  (input_string2[0]!='r') && (input_string2[1]!='b') &&
	  (input_string2[1]!='l') && (input_string2[1]!='r'))
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

void wcodemat_content_validate (int *errCode, int flow, int* int_input2)
{

	*errCode = SUCCESS;
	switch (flow){
		case 1:
			{
				break;
			}
		case 2:
			{
                if ((int_input2[0]<1) || (int_input2[0]>512))
					*errCode = UNKNOWN_INPUT_ERR;
				break;
			}
		case 3:
			{
				if ((int_input2[0]<1) || (int_input2[0]>512))
					*errCode = UNKNOWN_INPUT_ERR;

				break;
			}
		case 4:
			{
                if ((int_input2[0]<1) || (int_input2[0]>512))
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
