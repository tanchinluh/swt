/*
 * -------------------------------------------------------------------------
 * cwt_int.c -- continuous wavelet transform interface
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
#include "swt_common.h"
#include "cwt.h"
// #include "stack-c.h"
// #define __USE_DEPRECATED_STACK_FUNCTIONS__

/*-------------------------------------------
 * Haar Scale Filter Generation
 *-----------------------------------------*/

/*int
int_haar (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
  
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  haar_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);

  haar_content_validate(&errCode, l1, l2, l3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
 
  m4 = 1;
  n4 = istk(l3)[0];
  m5 = 1;
  n5 = n4;

  CreateVar (4, "d", &m4, &n4, &l4);
  CreateVar (5, "d", &m5, &n5, &l5);
  linspace(stk(l1)[0], stk(l2)[0], istk(l3)[0], stk(l4), n4);
  haar(stk(l4), n4, stk(l5), n5, 1);

  LhsVar(1) = 5;
  LhsVar(2) = 4;
  return 0;
}*/

/*-------------------------------------------
 * Sinus Scale Filter Generation
 *-----------------------------------------*/

int
int_sinus (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5;
   int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
   int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    double *input1;
  double *input2;
  double *input3;
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2;
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  sinus_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}	
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  sinus_content_validate(&errCode, input1, input2, input3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
 
  m4 = 1;
  n4 = input3[0];
  m5 = 1;
  n5 = n4;
//   CreateVar (5, "d", &m5, &n5, &l5);
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (4, "d", &m4, &n4, &l4);
    _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  linspace(input1[0], input2[0], input3[0], output1, n4);
  sinus(output1, n4, output2, n5, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}


/*-------------------------------------------
 * Poisson Scale Filter Generation
 *-----------------------------------------*/

int
int_poisson (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5;
   int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
     int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    double *input1;
  double *input2;
  double *input3;
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2;
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  poisson_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}	
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  poisson_content_validate(&errCode, input1, input2, input3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
 
  m4 = 1;
  n4 = input3[0];
  m5 = 1;
  n5 = n4;
//   CreateVar (5, "d", &m5, &n5, &l5);
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (4, "d", &m4, &n4, &l4);
      _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  linspace(input1[0], input2[0], input3[0], output1, n4);
  poisson(output1, n4, output2, n5, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}

/*-------------------------------------------
 * Mexican Hat Scale Filter Generation
 *-----------------------------------------*/

int
int_mexihat (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5;
   int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
     int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    double *input1;
  double *input2;
  double *input3;
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2;
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  mexihat_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}	
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  mexihat_content_validate(&errCode, input1, input2, input3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
 
  m4 = 1;
  n4 = input3[0];
  m5 = 1;
  n5 = n4;
//   CreateVar (5, "d", &m5, &n5, &l5);
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (4, "d", &m4, &n4, &l4);
      _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  linspace(input1[0], input2[0], input3[0], output1, n4);
  mexihat(output1, n4, output2, n5, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}

/*-------------------------------------------
 * Morlet Scale Filter Generation
 *-----------------------------------------*/

int
int_morlet (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5;
   int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
       int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    double *input1;
  double *input2;
  double *input3;
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2;
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  morlet_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  morlet_content_validate(&errCode, input1, input2, input3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m4 = 1;
  n4 = input3[0];
  m5 = 1;
  n5 = n4;
//   CreateVar (5, "d", &m5, &n5, &l5);
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (4, "d", &m4, &n4, &l4);
        _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  linspace(input1[0], input2[0], input3[0], output1, n4);
  morlet(output1, n4, output2, n5, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}

/*-------------------------------------------
 * DoG Scale Filter Generation
 *-----------------------------------------*/

int
int_DOGauss (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5;
   int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
       int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    double *input1;
  double *input2;
  double *input3;
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2;
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  DOGauss_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
	_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

 DOGauss_content_validate(&errCode, input1, input2, input3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m4 = 1;
  n4 = input3[0];
  m5 = 1;
  n5 = n4;
//   CreateVar (5, "d", &m5, &n5, &l5);
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (4, "d", &m4, &n4, &l4);
        _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  linspace(input1[0], input2[0], input3[0], output1, n4);
  DOGauss(output1, n4, output2, n5, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}

/*-------------------------------------------
 * Gauss Scale Filter Generation
 *-----------------------------------------*/

int
int_Gauswavf (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5, l6, m6, n6;
   int minlhs = 2, maxlhs = 2, minrhs = 4, maxrhs = 4;
    int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    int * p_input_vector4 = NULL;
    double *input1;
  double *input2;
  double *input3;
   double *input4; 
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2;
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  Gauss_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
					_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (4, "i", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

 Gauss_content_validate(&errCode, input1, input2, input3, input4);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m5 = 1;
  n5 = input3[0];
  m6 = 1;
  n6 = n5;

  //   CreateVar (6, "d", &m6, &n6, &l6);
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m6, n6, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (5, "d", &m5, &n5, &l5);
        _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m5, n5, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  linspace(input1[0], input2[0], input3[0], output1, n5);
  Gauss(output1, n5, output2, n6, input4[0], 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}


/*-------------------------------------------
 * Shannon Scale Filter Generation
 *-----------------------------------------*/

int
int_shanwavf (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5, l6, m6, n6;
   int l7r, l7c, m7, n7;
   int minlhs = 2, maxlhs = 2, minrhs = 5, maxrhs = 5;
      int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    int * p_input_vector4 = NULL;
    int * p_input_vector5 = NULL;
    double *input1;
  double *input2;
  double *input3;
   double *input4; 
   double *input5; 
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2_r, *output2_i;
  int errCode;
  int it;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  shanwavf_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
					_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (4, "d", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (5, "d", &m5, &n5, &l5);
_SciErr = getVarAddressFromPosition(pvApiCtx,5, &p_input_vector5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector5 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 5. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector5, &m5, &n5, &input5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  shanwavf_content_validate(&errCode, input1, input2, input3, input4, input5);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m6 = 1;
  n6 = input3[0];
  m7 = 1;
  n7 = n6;
  it = 1;
//   CreateCVar (7, "d",&it, &m7, &n7, &l7r, &l7c);
  _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m7, n7, &output2_r, &output2_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

//   CreateVar (6, "d", &m6, &n6, &l6);
          _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m6, n6, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  linspace(input1[0], input2[0], input3[0], output1, n6);
  shanwavf(output1, n6, input5[0], input4[0], output2_r, output2_i, n7,1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;

  return 0;
}

/*-------------------------------------------
 * Complex Morlet Scale Filter Generation
 *-----------------------------------------*/

int
int_cmorlet (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5, l6, m6, n6;
   int l7r, l7c, m7, n7;
   int minlhs = 2, maxlhs = 2, minrhs = 5, maxrhs = 5;
        int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    int * p_input_vector4 = NULL;
    int * p_input_vector5 = NULL;
    double *input1;
  double *input2;
  double *input3;
   double *input4; 
   double *input5; 
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2_r, *output2_i;
  int errCode;
  int it;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  cmorlet_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
					_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (4, "d", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (5, "d", &m5, &n5, &l5);
_SciErr = getVarAddressFromPosition(pvApiCtx,5, &p_input_vector5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector5 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 5. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector5, &m5, &n5, &input5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  cmorlet_content_validate(&errCode, input1, input2, input3, input4, input5);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m6 = 1;
  n6 = input3[0];
  m7 = 1;
  n7 = n6;
  it = 1;
//   CreateCVar (7, "d",&it, &m7, &n7, &l7r, &l7c);
  _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m7, n7, &output2_r, &output2_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (6, "d", &m6, &n6, &l6);
            _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m6, n6, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}



  linspace(input1[0], input2[0], input3[0], output1, n6);
  cmorlet(output1, n6, input5[0], input4[0], output2_r, output2_i, n7, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}

/*-------------------------------------------
 * Frequency B-Spline Scale Filter Generation
 *-----------------------------------------*/

int
int_fbspwavf (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5, l6, m6, n6;
   int l7, m7, n7, l8r, l8c, m8, n8;
   int minlhs = 2, maxlhs = 2, minrhs = 6, maxrhs = 6;
       int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    int * p_input_vector4 = NULL;
    int * p_input_vector5 = NULL;
       int * p_input_vector6 = NULL; 
    double *input1;
  double *input2;
  double *input3;
   double *input4; 
   double *input5; 
     double *input6; 
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2_r, *output2_i;
  int errCode;
  int it;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  fbspwavf_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
					_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (4, "i", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (5, "d", &m5, &n5, &l5);
_SciErr = getVarAddressFromPosition(pvApiCtx,5, &p_input_vector5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector5 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 5. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector5, &m5, &n5, &input5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (6, "d", &m6, &n6, &l6);
_SciErr = getVarAddressFromPosition(pvApiCtx,6, &p_input_vector6);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector6 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 6. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector6, &m6, &n6, &input6);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  fbspwavf_content_validate(&errCode,  input1, input2, input3, input4, input5, input6);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m7 = 1;
  n7 = input3[0];
  m8 = 1;
  n8 = n7;
  it = 1;
//   CreateCVar (8, "d",&it, &m8, &n8, &l8r, &l8c);
  _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m8, n8, &output2_r, &output2_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (7, "d", &m7, &n7, &l7);
              _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m7, n7, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}



  linspace(input1[0], input2[0], input3[0], output1, n7);
  fbspwavf(output1, n7, input4[0], input6[0], input5[0], output2_r, output2_i, n8, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}


/*-------------------------------------------
 * Complex Cauchy Scale Filter Generation
 *-----------------------------------------*/

int
int_cauchy (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5r, m5, n5, l5i;
  //static int l7r, l7c, m7, n7;
   int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
   int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    double *input1;
  double *input2;
  double *input3;
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2_r, *output2_i;
  int errCode;
  int it;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  cauchy_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
  //GetRhsVar (4, "d", &m4, &n4, &l4);
  //GetRhsVar (5, "d", &m5, &n5, &l5);

  cauchy_content_validate(&errCode, input1, input2, input3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m4 = 1;
  n4 = input3[0];
  m5 = 1;
  n5 = n4;
  it = 1;
//   CreateCVar (5, "d",&it, &m5, &n5, &l5r, &l5i);
    _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output2_r, &output2_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
  
//   CreateVar (4, "d", &m4, &n4, &l4);
      _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}



  linspace(input1[0], input2[0], input3[0], output1, n4);
  cauchy_neo(output1, n4, output2_r, output2_i, n5, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}

/*-------------------------------------------
 * Complex Gauss Scale Filter Generation
 *-----------------------------------------*/

int
int_cgauss (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l5, m5, n5, l6r, l6i, m6, n6;
  
   int minlhs = 2, maxlhs = 2, minrhs = 4, maxrhs = 4;
      int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
    int * p_input_vector4 = NULL;
    double *input1;
  double *input2;
  double *input3;
   double *input4; 
    SciErr _SciErr;
  int type;
   double *output1;
   double *output2_r, *output2_i;
  int errCode;
  int it;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  Gauss_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
					_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (4, "i", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4 ,&type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
  

  Gauss_content_validate(&errCode, input1, input2, input3, input4);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m5 = 1;
  n5 = input3[0];
  m6 = 1;
  n6 = n5;
  it = 1;
//   CreateCVar (6, "d",&it, &m6, &n6, &l6r, &l6i);
    _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m6, n6, &output2_r, &output2_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar (5, "d", &m5, &n5, &l5);
          _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m5, n5, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}


  linspace(input1[0], input2[0], input3[0], output1, n5);
  cgauss(output1, n5, input4[0], output2_r, output2_i, n6, 1);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  return 0;
}

/*-------------------------------------------
 * Meyer Filter Generation
 *-----------------------------------------*/

int
int_meyeraux (char *fname)
{
   int l1, m1, n1, l2, m2, n2;
   int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 1;
  int * p_input_vector1 = NULL;
  double *input1;
  SciErr _SciErr;
  int type;
  double *output1;
//   GetRhsVar (1, "d", &m1, &n1, &l1);
      _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
  m2 = 1;
  n2 = 1;
//   CreateVar (2, "d", &m2, &n2, &l2);
          _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m2, n2, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  meyeraux(input1[0],output1);
  LhsVar(1) = Rhs + 1;
  return 0;
}

/*-------------------------------------------
 * Scale and Wavelet Function
 *-----------------------------------------*/

int
int_wavefun (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l6, m6, n6, l7, m7, n7, l3r, l3i;
   int minlhs = 2, maxlhs = 5, minrhs = 2, maxrhs = 2;
  int ind, family, member, s1, s2, leng, count, level, it, l, errCode;
  double *lowfltr, *hifltr, *buff;
  double one, *phi, *psi;
  char a[2]="a",d[2]="d";
  Func syn_fun, ana_fun;
  swt_wavelet pWaveStruct;
  int * p_input_string = NULL;
   char * input_string = NULL;
  int * p_input_vector2 = NULL;
  double *input1;
  double *input2;
  SciErr _SciErr;
  int type;
  double *output1;
  double *output1_r, *output1_i;
  double *output2;
  double *output3;
  double *output4;
  double *output5;
  
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  wavefun_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "c", &m1, &n1, &l1);
_SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_string);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
     getAllocatedSingleString(pvApiCtx, p_input_string, &input_string);	 
//   GetRhsVar (2, "i", &m2, &n2, &l2);
_SciErr = getVarAddressFromPosition(pvApiCtx,2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: second input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
  
  wavefun_content_validate(&errCode, input_string, input2);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  level = (int)input2[0];
  if (level==0)
	  level = 8;
  wavelet_fun_parser (input_string, &ind);
  one = 1.0;
  
  if ((ind!=-1) && (wi[ind].rOrB==ORTH))
  {
	  if (Lhs!=3)
	  {
          Scierror(999,"Output arguments number should be 3!\n");
		  return 0;
	  }
      wavelet_parser(input_string,&family,&member);
	  syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  cwt_upcoef_len_cal (1, pWaveStruct.length, level, 
	       &s1, &s2);
	  swt_exp2(level, &leng);
      m3 = 1;
      n3 = leng*(pWaveStruct.length-1)+1;
	  m4 = 1;
	  n4 = n3;
	  m5 = 1;
	  n5 = n3;
// 	  CreateVar(3, "d", &m3, &n3, &l3);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m3, n3, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(4, "d", &m4, &n4, &l4);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(5, "d", &m5, &n5, &l5);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 3, m5, n5, &output3);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
      phi=(double *) malloc(n3*sizeof(double));
      psi=(double *) malloc(n3*sizeof(double));

	  for(count=0;count<n3;count++)
	  {
		  phi[count] = 0;
		  psi[count] = 0;
	  }
	  l=1;
      cwt_upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, phi+l, 
	      s1, s1, a, level);
	  cwt_upcoef (&one, 1, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, psi+l, 
	      s1, s1, d, level);
	  if ((family==COIFLETS) || (family==SYMLETS) || (family==DMEY))
	  {

		  for(count=0;count<n3;count++)
		  {
		      psi[count] = -1*psi[count];
		 }
	  }
	  for(count=0;count<n3;count++)
	  {
		  output1[count] = phi[count]*pow(sqrt(2),level);
		  output2[count] = psi[count]*pow(sqrt(2),level);
	  }
	  linspace(0.0, (double)(pWaveStruct.length-1), n3, output3, n3);
	  free(phi);
	  free(psi);
          LhsVar(1) = Rhs + 1;
	  LhsVar(2) = Rhs + 2;
	  LhsVar(3) = Rhs + 3;
	    if (input_string != NULL)
    freeAllocatedSingleString(input_string);
	  filter_clear();
      return 0;
  }

  if ((ind!=-1) && (wi[ind].rOrB==BIORTH))
  {
	  if (Lhs!=5)
	  {
          Scierror(999,"Output arguments number should be 5!\n");
		  return 0;
	  }
      wavelet_parser(input_string,&family,&member);
	  ana_fun = wi[ind].analysis;
      (*ana_fun)(member, &pWaveStruct);
	  cwt_upcoef_len_cal (1, pWaveStruct.length, level, 
	       &s1, &s2);
      swt_exp2(level, &leng);
      m3 = 1;
      n3 = leng*(pWaveStruct.length-1)+1;
	  m4 = 1;
	  n4 = n3;
	  m5 = 1;
	  n5 = n3;
	  m6 = 1;
	  n6 = n3;
	  m7 = 1;
	  n7 = n3;
// 	  CreateVar(3, "d", &m3, &n3, &l3);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m3, n3, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(4, "d", &m4, &n4, &l4);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(5, "d", &m5, &n5, &l5);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 3, m5, n5, &output3);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(6, "d", &m6, &n6, &l6);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 4, m6, n6, &output4);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(7, "d", &m7, &n7, &l7);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 5, m7, n7, &output5);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  for(count=0;count<n3;count++)
	  {
		  output1[count] = 0;
		  output2[count] = 0;
		  output3[count] = 0;
		  output4[count] = 0;
	  }
	  lowfltr = (double *)  malloc(pWaveStruct.length*sizeof(double));
	  hifltr =(double *)  malloc(pWaveStruct.length*sizeof(double));
      wrev(pWaveStruct.pLowPass, pWaveStruct.length, lowfltr, pWaveStruct.length);
	  qmf_wrev(lowfltr,pWaveStruct.length,hifltr,pWaveStruct.length);
	  l=1;
      cwt_upcoef (&one, 1, lowfltr, hifltr, pWaveStruct.length, output1+l, 
	      s1, s1, a, level);
       cwt_upcoef (&one, 1, lowfltr, hifltr, pWaveStruct.length, output2+l, 
	      s1, s1, d, level);
	  free(lowfltr);
	  free(hifltr);
	  filter_clear();
      syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
       cwt_upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, output3+l, 
	      s1, s1, a, level);
	   cwt_upcoef (&one, 1, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, output4+l, 
	      s1, s1, d, level);
	   for(count=0;count<n3;count++)
	     {
			  output1[count] = output1[count]*pow(sqrt(2),level);
		      output2[count] = -1*output2[count]*pow(sqrt(2),level);
			  output3[count] = output3[count]*pow(sqrt(2),level);
		     output4[count] = -1*output4[count]*pow(sqrt(2),level);
		 }
	  linspace(0.0, (double)(pWaveStruct.length-1), n3, output5, n3);
	  LhsVar(1) = Rhs + 1;
	  LhsVar(2) = Rhs + 2;
	  LhsVar(3) = Rhs + 3;
	  LhsVar(4) = Rhs + 4;
	  LhsVar(5) = Rhs + 5;
	    if (input_string != NULL)
    freeAllocatedSingleString(input_string);
	  filter_clear();
	  return 0;
  }

  cwt_fun_parser(input_string, &ind);
  if ((ind!=-1) && (ci[ind].phipsi==PSI_ONLY) && (ci[ind].realOrComplex==REAL) )
  {
	  if (Lhs!=2)
	  {
          Scierror(999,"Output arguments number should be 2!\n");
		  return 0;
	  }
	  swt_exp2(level, &leng);
      m3 = 1;
      n3 = leng;
	  m4 = 1;
	  n4 = n3;
// 	  CreateVar(3, "d", &m3, &n3, &l3);
	  	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m3, n3, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(4, "d", &m4, &n4, &l4);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  linspace(ci[ind].lb, ci[ind].ub, n3, output2, n3);
	  (*(ci[ind].scalef))(output2,n3,output1,n3,ci[ind].cpsi);
	  LhsVar(1) = Rhs + 1;
	  LhsVar(2) = Rhs + 2;
	    if (input_string != NULL)
    freeAllocatedSingleString(input_string);
	  return 0;
  }

  if ((ind!=-1) && (ci[ind].phipsi==PSI_ONLY) && (ci[ind].realOrComplex==COMPLEX) )
  {
	  if (Lhs!=2)
	  {
          Scierror(999,"Output arguments number should be 2!\n");
		  return 0;
	  }
	  swt_exp2(level, &leng);
      m3 = 1;
      n3 = leng;
	  m4 = 1;
	  n4 = n3;
	  it = 1;
//       CreateCVar(3, "d", &it, &m3, &n3, &l3r, &l3i);
	 _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m3, n3, &output1_r,&output1_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(4, "d", &m4, &n4, &l4);
	  	  	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  buff = (double *) malloc(n3*2*sizeof(double));
	  linspace(ci[ind].lb, ci[ind].ub, n3, output2, n3);
	  (*(ci[ind].scalef))(output2,n3,buff,n3,ci[ind].cpsi);
	  for(count=0;count<n3;count++)
	  {
        output1_r[count]=buff[count];
		output1_i[count]=buff[count+n3];
	  }
	  free(buff);
	  LhsVar(1) = Rhs + 1;
	  LhsVar(2) = Rhs + 2;
	    if (input_string != NULL)
    freeAllocatedSingleString(input_string);
	  return 0;
  }
     return 0;
}

/*-------------------------------------------
 * 2D Scale and Wavelet Function 
 *-----------------------------------------*/

int
int_wavefun2 (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l6, m6, n6, l7, m7, n7;
   int minlhs = 5, maxlhs = 5, minrhs = 2, maxrhs = 2;
     int * p_input_string = NULL;
   char * input_string = NULL;
   int * p_input_vector2 = NULL;
   double *input2;
   SciErr _SciErr;
  int type;
  double *output1;
  double *output2;
  double *output3;
  double *output4;
  double *output5;
  
  int ind, family, member, s1, s2, leng, count, level, row, col, l, errCode;
  double one, *phi, *psi, *xval;
  char a[2]="a",d[2]="d";
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  wavefun2_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

//   GetRhsVar (1, "c", &m1, &n1, &l1);
_SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_string);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
     getAllocatedSingleString(pvApiCtx, p_input_string, &input_string);	 
//   GetRhsVar (2, "i", &m2, &n2, &l2);
     _SciErr = getVarAddressFromPosition(pvApiCtx,2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: second input vector must be double or int\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
  
  wavefun2_content_validate(&errCode, input_string, input2);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  level = input2[0];
  if (level==0)
	  level = 4;

  wavelet_fun_parser (input_string, &ind);
  one = 1.0;

 if ((ind!=-1) && (wi[ind].rOrB==ORTH))
  {
      wavelet_parser(input_string,&family,&member);
	  syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  cwt_upcoef_len_cal (1, pWaveStruct.length, level, 
	       &s1, &s2);
	  
	  swt_exp2(level, &leng);
      m3 = leng*(pWaveStruct.length-1)+1;
      n3 = leng*(pWaveStruct.length-1)+1;
	  m4 = m3;
	  n4 = n3;
	  m5 = m3;
	  n5 = n3;
	  m6 = m3;
	  n6 = n3;
	  m7 = m3;
	  n7 = n3;
// 	  CreateVar(3, "d", &m3, &n3, &l3);
	 _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m3, n3, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(4, "d", &m4, &n4, &l4);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m4, n4, &output2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(5, "d", &m5, &n5, &l5);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 3, m5, n5, &output3);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(6, "d", &m6, &n6, &l6);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 4, m6, n6, &output4);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
// 	  CreateVar(7, "d", &m7, &n7, &l7);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 5, m7, n7, &output5);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

	  phi=(double *) malloc(n3*sizeof(double));
	  psi=(double *) malloc(n3*sizeof(double));
      xval=(double *) malloc(n3*sizeof(double));
	  for(count=0;count<n3;count++)
	  {
		  phi[count] = 0;
		  psi[count] = 0;
		
	  }
	  l=(int)(floor((n3-s1)/2));
      cwt_upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, phi+l, 
	      s1, s1, a, level);
	  cwt_upcoef (&one, 1, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, psi+l, 
	      s1, s1, d, level);
	  linspace(0.0, (double)(pWaveStruct.length), n3, xval, n3);
	  if ((family==COIFLETS) || (family==SYMLETS) || (family==DMEY))
	  {

		  for(count=0;count<n3;count++)
	     {
		      psi[count] = -1*(psi[count]);
		 }
	  }
	  for (col=0;col<n3;col++)
	  {
		  for(row=0;row<n3;row++)
		  {
			  output1[row+col*n3]=phi[col]*phi[row]*pow(sqrt(2),level*2);
			  output2[row+col*n3]=psi[col]*phi[row]*pow(sqrt(2),level*2);
              output3[row+col*n3]=phi[col]*psi[row]*pow(sqrt(2),level*2);
			  output4[row+col*n3]=psi[col]*psi[row]*pow(sqrt(2),level*2);
			  output5[row+col*n3]=xval[col]*xval[row];
		  }
	  }
	  free(phi);
	  free(psi);
	  free(xval);
	  LhsVar(1) = Rhs + 1;
	  LhsVar(2) = Rhs + 2;
	  LhsVar(3) = Rhs + 3;
	  LhsVar(4) = Rhs + 4;
	  LhsVar(5) = Rhs + 5;
	  if (input_string != NULL)
    freeAllocatedSingleString(input_string);
	  filter_clear();
      return 0;
  }
  if ((ind!=-1) && (wi[ind].rOrB==BIORTH))
  {

          Scierror(999,"wavefun2 does not work with BIORTH splines!\n");
		  return 0;
  }
  cwt_fun_parser(input_string, &ind);
  if ((ind!=-1) && (ci[ind].phipsi==PSI_ONLY) && (ci[ind].realOrComplex==REAL) )
  {
    Scierror(999,"wavefun2 does not work with REAL splines!\n");
		  return 0;
  }
  
  if ((ind!=-1) && (ci[ind].phipsi==PSI_ONLY) && (ci[ind].realOrComplex==COMPLEX) )
  {
    Scierror(999,"wavefun2 does not work with COMPLEX splines!\n");
		  return 0;
  }
   return 0;
}

/*-------------------------------------------
 * CWT Utility
 *-----------------------------------------*/
/*int
int_wpsi (char *fname)
{
    static int l1, m1, n1, l2, m2, n2, l2r, l2i;
	static int minlhs=1, maxlhs = 1, minrhs = 1, maxrhs = 1;
	int ind1, ind2, it, count, family, member;
	double *fbuf;
	Func syn_fun;
    swt_wavelet pWaveStruct;

    CheckRhs (minrhs, maxrhs);
    CheckLhs (minlhs, maxlhs);

    GetRhsVar (1, "c", &m1, &n1, &l1);
    
	m2 = 1;
	it = 1;

	wavelet_fun_parser (cstk(l1), &ind1);
	cwt_fun_parser(cstk(l1), &ind2);

    if (ind1!=-1)
	{
       wavelet_parser(cstk(l1),&family,&member);
	   syn_fun = wi[ind1].synthesis;
       (*syn_fun)(member, &pWaveStruct);
	   n2=1024*(pWaveStruct.length-1)+1;
	   filter_clear();
	   CreateVar(2, "d", &m2, &n2, &l2);
	   full_range_scalef (cstk(l1), stk(l2), n2);
	}
	else if ((ind2!=-1) && (ci[ind2].realOrComplex==REAL))
	{
	   n2 = 1024;
       CreateVar(2, "d", &m2, &n2, &l2);
	   full_range_scalef (cstk(l1), stk(l2), n2);
	}
	else if ((ind2!=-1) && (ci[ind2].realOrComplex==COMPLEX))
	{
		n2 = 1024;
        CreateCVar(2, "d", &it, &m2, &n2, &l2r, &l2i);
        fbuf=malloc(n2*2*sizeof(double));
		full_range_scalef (cstk(l1), fbuf, 2*n2);
		for(count=0;count<n2;count++)
		{
			stk(l2r)[count]=fbuf[count];
			stk(l2i)[count]=fbuf[count+n2];
		}
		free(fbuf);
	}
	else
		return 0;

    LhsVar(1) = 2;
	return 0;
}

int
int_cwtscale (char *fname)
{
	static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
	static int minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
    double delta;

	GetRhsVar (1, "d", &m1, &n1, &l1);
    GetRhsVar (2, "i", &m2, &n2, &l2);

	m3 = 1;
	cwt_len_cal (m1*n1, istk(l2)[0], &n3, &delta);
	CreateVar(3, "d", &m3, &n3, &l3);
    scale_real (stk(l1), m1*n1, delta, stk(l3), n3);
    LhsVar(1) = 3;
	return 0;
}

int
int_cwtconv (char *fname)
{
    static int l1, m1, n1, l2, m2, n2, l2r, l2i, l3, m3, n3, l3r, l3i;
	static int minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	int it;

    CheckRhs (minrhs, maxrhs);
    CheckLhs (minlhs, maxlhs);

    GetRhsVar (1, "d", &m1, &n1, &l1);
    GetRhsCVar (2, "d", &it, &m2, &n2, &l2r, &l2i);

	m3 = 1;
	n3 = m1*n1;
	if (it!=1)
	{
        CreateVar(3, "d", &m3, &n3, &l3);
		cwt_conv_real (stk(l1), m1*n1, stk(l2), m2*n2, stk(l3), n3);
	}
	else
	{
		CreateCVar(3, "d", &it, &m3, &n3, &l3r, &l3i);
		cwt_conv_complex (stk(l1), m1*n1, stk(l2r), stk(l2i), m2*n2, 
					   stk(l3r), stk(l3i), m3*n3);
	}

	LhsVar(1) = 3; 
	return 0;
}*/

/*-------------------------------------------
 * Continuous Wavelet Transform
 *-----------------------------------------*/

int
int_cwt (char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   int l4, m4, n4, l4r, l4i, l1r, l1i;
   int minlhs = 1, maxlhs = 1, minrhs = 3, maxrhs = 3;
     int * p_input_vector1 = NULL;
    int * p_input_vector2 = NULL;
    int * p_input_string = NULL;
   char * input_string = NULL;
    double *input1,*input1_r,*input1_i;
  double *input2;
    SciErr _SciErr;
  int type;
   double *output1, *output1_r, *output1_i;
  int ind1, ind2, count, len, plen, family, member, i, it, scale, errCode, flow;
  double delta;
  double *psi, *f, *temp, *tempi;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  cwt_form_validate(&errCode, &flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
 if (flow==1)
 {
//   GetRhsVar (1, "d", &m1, &n1, &l1);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//   GetRhsVar (2, "d", &m2, &n2, &l2);
					_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar (3, "c", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_input_string);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
     getAllocatedSingleString(pvApiCtx, p_input_string, &input_string);	

  cwt_content_validate(&errCode, input_string);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m4 = m2*n2;
  n4 = m1*n1;
  it = 1;
  
  wavelet_fun_parser (input_string, &ind1);
  cwt_fun_parser(input_string, &ind2);

  if (ind1!=-1) 
  {
	  
	  wavelet_parser(input_string,&family,&member);
	  syn_fun = wi[ind1].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  plen=1024*(pWaveStruct.length-1)+1;
	  filter_clear();
	  psi=(double *) malloc(plen*sizeof(double));
	  full_range_scalef (input_string, psi, plen);
	  //ocumsum(psi,plen);
// 	  CreateVar(4, "d", &m4, &n4, &l4);
	  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  temp = (double *) malloc(m4*n4*sizeof(double));
	  //tim = (int)(floor(len/n4));
	  //scale = (int)(floor(stk(l2)[0]));
	  for (count=0;count<m4;count++)
	  {
		  scale = (int)(floor(input2[count]));
		  if (scale<1)
			  scale = 1;
		  cwt_len_cal (plen, (scale*pWaveStruct.length+1), &len, &delta);
		  f = (double *) malloc(len*sizeof(double));
		  scale_real (psi, plen, delta, f, len);
		  
		  for(i=0;i<len;i++)
			f[i]=f[i]/sqrt(scale);
		  cwt_conv_real (input1, n4, f, len, temp+count*n4, n4);
		  free(f);
	  }
	  
	  free(psi);
	  matrix_tran(temp,m4,n4,output1,n4,m4);
	  free(temp);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==REAL))
  {
	 
	  
// 	  CreateVar(4, "d", &m4, &n4, &l4);
	_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m4, n4, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  temp = (double *) malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  //scale = (int)(floor(stk(l2)[count]));
		  plen=(int)((ci[ind2].ub-ci[ind2].lb)*(input2[count]))+1;
		  if (plen<3)
			  plen = 3;
	      psi=(double *) malloc(plen*sizeof(double));
	      full_range_scalef (input_string, psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(input2[count]);
		  cwt_conv_real (input1, n4, psi, plen, temp+count*n4, n4);
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,output1,n4,m4);
	  free(temp);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==COMPLEX))
  {
	 
	  
// 	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	_SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m4, n4, &output1_r,&output1_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  temp = (double *) malloc(m4*n4*sizeof(double));
	  tempi = (double *) malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  plen=((int)((ci[ind2].ub-ci[ind2].lb)*(input2[count])+1))*2;//(int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count]+1))*2;
		  if (plen<3)
			  plen = 3;
	      psi=(double *) malloc(plen*sizeof(double));
	      full_range_scalef (input_string, psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(input2[count]);
		  cwt_conv_complex (input1, n4, psi, psi+plen/2, plen/2, temp+count*n4, tempi+count*n4, n4); 
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,output1_r,n4,m4);
	  matrix_tran(tempi,m4,n4,output1_i,n4,m4);
	  free(temp);
	  free(tempi);
  }
  LhsVar(1) = Rhs + 1;
  if (input_string != NULL)
    freeAllocatedSingleString(input_string);
  return 0;
 }
 
 if (flow==2)
 {
// 	  GetRhsCVar (1, "d", &it, &m1, &n1, &l1r, &l1i);
    _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_vector1);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector1, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: first input vector must be double\n",fname);	
		  return -1;
		}
		_SciErr = getComplexMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1_r, &input1_i);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
		}
//       GetRhsVar (2, "d", &m2, &n2, &l2);
					_SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 2. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//       GetRhsVar (3, "c", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_input_string);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
 getAllocatedSingleString(pvApiCtx, p_input_string, &input_string);
 
	   cwt_content_validate(&errCode, input_string);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m4 = m2*n2;
  n4 = m1*n1;
  it = 1;
  
  wavelet_fun_parser (input_string, &ind1);
  cwt_fun_parser(input_string, &ind2);

  if (ind1!=-1) 
  {
	  
	  wavelet_parser(input_string,&family,&member);
	  syn_fun = wi[ind1].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  plen=1024*(pWaveStruct.length-1)+1;
	  filter_clear();
	  psi=(double *) malloc(plen*sizeof(double));
	  full_range_scalef (input_string, psi, plen);

// 	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	 _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m4, n4, &output1_r,&output1_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  temp = (double *) malloc(m4*n4*sizeof(double));
	  tempi = (double *) malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  scale = (int)(floor(input2[count]));
		  if (scale<1)
			  scale = 1;
		  cwt_len_cal (plen, (scale*pWaveStruct.length+1), &len, &delta);
		  f = (double *) malloc(len*sizeof(double));
		  scale_real (psi, plen, delta, f, len);
		  
		  for(i=0;i<len;i++)
			f[i]=f[i]/sqrt(scale);
		  cwt_conv_real (input1_r, n4, f, len, temp+count*n4, n4);
		  cwt_conv_real (input1_i, n4, f, len, tempi+count*n4, n4);
	  }
	  free(f);
	  free(psi);
	  matrix_tran(temp,m4,n4,output1_r,n4,m4);
      matrix_tran(tempi,m4,n4,output1_i,n4,m4);
	  free(temp);
	  free(tempi);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==REAL))
  {
	 
	  
// 	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	_SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m4, n4, &output1_r,&output1_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  temp = (double *) malloc(m4*n4*sizeof(double));
	  tempi = (double *) malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  //scale = (int)(floor(stk(l2)[count]));
		  plen=(int)((ci[ind2].ub-ci[ind2].lb)*(input2[count]))+1;
		  if (plen<3)
			  plen = 3;
	      psi=(double *) malloc(plen*sizeof(double));
	      full_range_scalef (input_string, psi, plen);
          for(i=0;i<plen;i++)
		  psi[i]=psi[i]/sqrt(input2[count]);
		  cwt_conv_real (input1_r, n4, psi, plen, temp+count*n4, n4);
		  cwt_conv_real (input1_i, n4, psi, plen, tempi+count*n4, n4);
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,output1_r,n4,m4);
	  matrix_tran(tempi,m4,n4,output1_i,n4,m4);
	  free(temp);
	  free(tempi);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==COMPLEX))
  {
	 
	  
// 	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	_SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m4, n4, &output1_r,&output1_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	  temp = (double *) malloc(m4*n4*sizeof(double));
	  tempi = (double *) malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  plen=((int)((ci[ind2].ub-ci[ind2].lb)*(input2[count])+1))*2;//(int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count]+1))*2;
		  if (plen<3)
			  plen = 3;
	      psi=(double *) malloc(plen*sizeof(double));
	      full_range_scalef (input_string, psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(input2[count]);
		  cwt_conv_complex_complex (input1_r, input1_i, n4, psi, psi+plen/2, plen/2, temp+count*n4, tempi+count*n4, n4); 
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,output1_r,n4,m4);
	  matrix_tran(tempi,m4,n4,output1_i,n4,m4);
	  free(temp);
	  free(tempi);
  }

      LhsVar(1) = Rhs + 1;
      if (input_string != NULL)
	freeAllocatedSingleString(input_string);
      return 0;
 }
}

/*-------------------------------------------
 * Continuous Wavelet Transform Type II
 *-----------------------------------------*/

/*int
int_cwtzwei (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, lrx, lcx;
  static int minlhs = 1, maxlhs = 1, minrhs = 3, maxrhs = 3;
  int ind2, count, len, plen, i;
  double *psi, *temp;


  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "i", &m2, &n2, &l2);
  GetRhsVar (3, "c", &m3, &n3, &l3);


  m4 = m2*n2;
  n4 = m1*n1;
  

  cwt_fun_parser(cstk(l3), &ind2);

  if ((ind2!=-1) && (ci[ind2].realOrComplex==REAL))
  {
	 
	  
	  CreateVar(4, "d", &m4, &n4, &l4);
	  temp = malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  plen=(ci[ind2].ub-ci[ind2].lb)*istk(l2)[count]+1;
	      psi=malloc(plen*sizeof(double));
	      full_range_scalef (cstk(l3), psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(istk(l2)[count]);
		  cwt_conv_real (stk(l1), n4, psi, plen, temp+count*n4, n4);
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,stk(l4),n4,m4);
	  free(temp);
  }

  LhsVar(1) = 4;

  return 0;
}*/