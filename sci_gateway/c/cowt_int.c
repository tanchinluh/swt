/*
 * -------------------------------------------------------------------------
 * cowt_int.c -- dual tree complex wavelet transform interface
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
#include "dwt.h"
// #define __USE_DEPRECATED_STACK_FUNCTIONS__
// #include "stack-c.h"

int
int_FSfarras(char *fname)
{
  int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int minlhs=1, maxlhs=2, minrhs=1, maxrhs=1;
  swt_wavelet pWaveStruct;
  Func ana_fun,syn_fun;
  int family, member, ii, errCode;
  char s1[]={"fa1"};
  char s2[]={"fa2"};
  double *var1,*var2;
  int type1;
  int * p_input_string = NULL;
  char * input_string = NULL;
  SciErr _SciErr;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

//   GetRhsVar(1, "c", &m1, &n1, &l1);

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_string);
   if(_SciErr.iErr)
    {
	   printError(&_SciErr, 0);
	   return 0;
    }

    getAllocatedSingleString(pvApiCtx, p_input_string, &input_string);	    

  var1 = (double *)malloc(40*sizeof(double));
  var2 = (double *)malloc(40*sizeof(double));

  wavelet_parser(s1,&family,&member);
  wavelet_fun_parser (s1, &ii);
  ana_fun = wi[ii].analysis;
  (*ana_fun)(member, &pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var1, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var1+pWaveStruct.length, pWaveStruct.length);

  syn_fun = wi[ii].synthesis;
  (*syn_fun)(member,&pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var2, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var2+pWaveStruct.length, pWaveStruct.length);

  wavelet_parser(s2,&family,&member);
  wavelet_fun_parser (s2, &ii);
  ana_fun = wi[ii].analysis;
  (*ana_fun)(member, &pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var1+20, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var1+30, pWaveStruct.length);

  syn_fun = wi[ii].synthesis;
  (*syn_fun)(member,&pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var2+20, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var2+30, pWaveStruct.length);

  //if ((cstk(l1)[0] == 'f') || (cstk(l1)[0] == 'F'))
  if ((input_string[0] == 'f') || (input_string[0] == 'F'))  
    {
      double *ptr,*ptr2;
      m2 = 4;
      n2 = 10;
      m3 = 4;
      n3 = 10;
//       CreateVar(2, "d", &m2, &n2, &l2);
//       CreateVar(3, "d", &m3, &n3, &l3);
      _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m2, n2, &ptr);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }
       _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m3, n3, &ptr2);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }	         
      
      matrix_tran(var1,m2,n2,ptr,n2,m2);
      matrix_tran(var2,m3,n3,ptr2,n3,m3);
      LhsVar(1) = Rhs+1;
      LhsVar(2) = Rhs+2;
    }
//   else if ((cstk(l1)[0] == 'a') || (cstk(l1)[0] == 'A'))
  else if ((input_string[0] == 'a') || (input_string[0] == 'A'))
    {
      double *ptr;
      m2 = 4;
      n2 = 10;
//       CreateVar(2, "d", &m2, &n2, &l2);
       _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m2, n2, &ptr);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }
//     matrix_tran(var1,m2,n2,stk(l2),n2,m2);
      matrix_tran(var1,m2,n2,ptr,n2,m2);
      LhsVar(1) = Rhs + 1;
    }
//   else if ((cstk(l1)[0] == 's') || (cstk(l1)[0] == 'S'))
  else if ((input_string[0] == 's') || (input_string[0] == 'S')) 
    {
      double *ptr;
      m2 = 4;
      n2 = 10;
//       CreateVar(2, "d", &m2, &n2, &l2);
      _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m2, n2, &ptr);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }
//       matrix_tran(var2,m2,n2,stk(l2),n2,m2);
      matrix_tran(var2,m2,n2,ptr,n2,m2);
      LhsVar(1) = Rhs + 1;
    }
  else
    {
      errCode = UNKNOWN_INPUT_ERR;
      validate_print (errCode);
    }
  if (input_string != NULL)
    freeAllocatedSingleString(input_string);
  free(var1);
  free(var2);
  return 0;
}


int
int_dualfilt1(char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int minlhs=1, maxlhs=2, minrhs=1, maxrhs=1;
  swt_wavelet pWaveStruct;
  Func ana_fun,syn_fun;
  int family, member, ii, errCode;
  char s1[]={"ksq1"};
  char s2[]={"ksq2"};
  double *var1,*var2;
  int type1;
  int * p_input_string = NULL;
  char * input_string = NULL;
  SciErr _SciErr;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

//   GetRhsVar(1, "c", &m1, &n1, &l1);
   _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_input_string);
   if(_SciErr.iErr)
    {
	   printError(&_SciErr, 0);
	   return 0;
    }

    getAllocatedSingleString(pvApiCtx, p_input_string, &input_string);	    


  var1 = (double *)malloc(40*sizeof(double));
  var2 = (double *)malloc(40*sizeof(double));

  wavelet_parser(s1,&family,&member);
  wavelet_fun_parser (s1, &ii);
  ana_fun = wi[ii].analysis;
  (*ana_fun)(member, &pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var1, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var1+pWaveStruct.length, pWaveStruct.length);

  syn_fun = wi[ii].synthesis;
  (*syn_fun)(member,&pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var2, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var2+pWaveStruct.length, pWaveStruct.length);

  wavelet_parser(s2,&family,&member);
  wavelet_fun_parser (s2, &ii);
  ana_fun = wi[ii].analysis;
  (*ana_fun)(member, &pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var1+20, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var1+30, pWaveStruct.length);

  syn_fun = wi[ii].synthesis;
  (*syn_fun)(member,&pWaveStruct);

  verbatim_copy(pWaveStruct.pLowPass, pWaveStruct.length, var2+20, pWaveStruct.length);
  verbatim_copy(pWaveStruct.pHiPass, pWaveStruct.length, var2+30, pWaveStruct.length);

//   if (cstk(l1)[0] == 'f')
  if (input_string[0] == 'f')
    {
      double *ptr,*ptr2;
      m2 = 4;
      n2 = 10;
      m3 = 4;
      n3 = 10;
//       CreateVar(2, "d", &m2, &n2, &l2);
//       CreateVar(3, "d", &m3, &n3, &l3);
      _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m2, n2, &ptr);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }
       _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m3, n3, &ptr2);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }	         
      
      matrix_tran(var1,m2,n2,ptr,n2,m2);
      matrix_tran(var2,m3,n3,ptr2,n3,m3);
      LhsVar(1) = Rhs +1;
      LhsVar(2) = Rhs +2;
    }
//   else if (cstk(l1)[0] == 'a')
  else if (input_string[0] == 'a')
    {
      double *ptr;
      m2 = 4;
      n2 = 10;
//       CreateVar(2, "d", &m2, &n2, &l2);
            _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m2, n2, &ptr);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }
      matrix_tran(var1,m2,n2,ptr,n2,m2);
      LhsVar(1) = Rhs + 1;
    }
//   else if (cstk(l1)[0] == 's')
  else if (input_string[0] == 's')
    {
       double *ptr;
      m2 = 4;
      n2 = 10;
//       CreateVar(2, "d", &m2, &n2, &l2);
            _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m2, n2, &ptr);
			if(_SciErr.iErr)
		        {
			    printError(&_SciErr, 0);
			    return 0;
		         }
      matrix_tran(var2,m2,n2,ptr,n2,m2);
      LhsVar(1) = Rhs + 1;
    }
  else
    {
      errCode = UNKNOWN_INPUT_ERR;
      validate_print (errCode);
    }
    
  if (input_string != NULL)
    freeAllocatedSingleString(input_string);
  free(var1);
  free(var2);
  return 0;
}

int
int_dualtree(char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l6, m6, n6, l7, m7, n7, l8, m8, n8;
   int l5r, l5c;
  static int minlhs=2, maxlhs=2, minrhs=4, maxrhs=6;
  int errCode, flow, calLen, temLen, count, ln, it, stride, val;
  double *f1, *f2;
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
   
  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  dualtree_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
  {
    validate_print (errCode);
    return 0;			
  }

//   GetRhsVar(1, "d", &m1, &n1, &l1);
// GetRhsVarDouble(1,  &m1, &n1, input1);
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
			return -1;
		}
//   GetRhsVar(2, "i", &m2, &n2, &l2);
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
//   GetRhsVar(3, "d", &m3, &n3, &l3);
    		_SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(4, "d", &m4, &n4, &l4);
    		_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  wave_len_validate (m1*n1, n3, &stride, &val);
  if ((!val) || (stride<input2[0]))
    {
      sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
      return 0;
    }

  if (dwtMode==PER)
    {
      ln = 0;
      calLen = n1 * m1;
      for (count = 0; count < input2[0]; count++)
	{
	  calLen = (int)ceil(((double)(calLen))/2.0);
	  temLen = calLen;
	  ln += temLen;
	}
      ln += temLen;
    }
  else
    {
      ln = 0;
      calLen = n1 * m1;
      for (count = 0; count < input2[0]; count++)
	{
	  calLen += n3 - 1;
	  temLen = calLen/2;
	  ln += temLen;
	  calLen = temLen;
	}
      ln += temLen;
    }

  f1 = (double *)malloc(m3*n3*sizeof(double));
  f2 = (double *)malloc(m3*n3*sizeof(double));
  matrix_tran(input3,n3,m3,f1,m3,n3);
  matrix_tran(input4,n3,m3,f2,m3,n3);
  
  switch (flow) {
  case 1:
    {
      double *ptr_r,*ptr_i,*ptr2;
      int *ptr2_int;
      int i,j;
      //sciprint("flow 1\n");
      m5 = 1;
      n5 = ln;
      m6 = 1;
      n6 = input2[0] + 2;
      it = 1;
//       CreateCVar (5, "d",&it, &m5, &n5, &l5r, &l5c);

	_SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &ptr_r,&ptr_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//       CreateVar(6, "i", &m6, &n6, &l6);
      	_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m6, n6, &ptr2);
	//_SciErr = allocMatrixOfInteger32(pvApiCtx, Rhs + 2, m6, n6, &ptr2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	ptr2_int = (int*)malloc(sizeof(int) * m6 * n6);
      wave_dec_len_cal (n3, m1*n1, (int)input2[0], ptr2_int);
      cowavedec (input1, m1*n1, ptr_r, ptr_i, n5, 
		 f1, f1+n3, f1+n3*2, f1+n3*3,
		 f2, f2+n3, f2+n3*2, f2+n3*3, 
		 n3, ptr2_int, n6, (int)input2[0], dwtMode);
      
      for(i = 0 ; i < m6 ; i++)
        for(j = 0 ; j < n6 ; j++)
                    ptr2[i + m6 * j] = (double)ptr2_int[i + m6 * j];
        
    
    free(ptr2_int);
      LhsVar(1) = Rhs + 1;
      LhsVar(2) = Rhs + 2;
      break;
    }
  default:
    sciprint("input not valid\n");
    break;
  }

  free(f1);
  free(f2);

  return 0;
}

int
int_idualtree(char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l1r, l1c;
   int minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
  int errCode, flow, it, val, count, len;
  double *f1, *f2;
  
  int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
  int * p_input_vector4 = NULL;
  double *input1_r,*input1_i;
  double *input2;
  int *input2_int;
  double *input3;
  double *input4;
  SciErr _SciErr;
  int type;
   double *output1;
   int i,j;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  idualtree_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
  {
    validate_print (errCode);
    return 0;			
  }

//   GetRhsCVar(1, "d", &it, &m1, &n1, &l1r, &l1c);
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
		_SciErr = getComplexMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1_r,&input1_i);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(2, "i", &m2, &n2, &l2);
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
//   GetRhsVar(3, "d", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(4, "d", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  len = 0;
  for (count = 0; count < (m2 * n2 - 1); count++)
    len += input2[count];
  if (len != m1 * n1)
    {
      sciprint ("Inputs are not coef and length array!\n");
      return 0;
    }
  val = 0;
  for (count = 0; count < (m2 * n2 - 1); count++)
    {
      if (input2[count] > input2[count + 1])
	{
	  val = 1;
	  break;
	}
    }
  if (val != 0)
    {
      sciprint ("Inputs are not coef and length array!\n");
      return 0;
    }

  if (input2[0] < n3)
    {
      sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
      return 0;
    }

  m5 = 1;
  n5 = input2[m2*n2-1];
//   CreateVar(5, "d", &m5, &n5, &l5);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
  f1 = (double *)malloc(m3*n3*sizeof(double));
  f2 = (double *)malloc(m3*n3*sizeof(double));
  matrix_tran(input3,n3,m3,f1,m3,n3);
  matrix_tran(input4,n3,m3,f2,m3,n3);
  input2_int = (int*)malloc(sizeof(int) * m2 * n2);
        for(i = 0 ; i < m2 ; i++)
        for(j = 0 ; j < n2 ; j++)
                    input2_int[i + m2 * j] = (int)input2[i + m2 * j];

  cowaverec (input1_r, input1_i, m1*n1, output1, m5*n5, 
	     f1, f1+n3, f1+n3*2, f1+n3*3,
	     f2, f2+n3, f2+n3*2, f2+n3*3, 
	     n3, input2_int, m2*n2, m2*n2-2, dwtMode);
  
  free(input2_int);
  free(f1);
  free(f2);
  LhsVar(1) = Rhs + 1;

  return 0;
}


int
int_dualtree2D(char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l6, m6, n6, l7, m7, n7, l8, m8, n8;
   int l5r, l5c;
   int minlhs=2, maxlhs=2, minrhs=4, maxrhs=4;
  int errCode, flow, calLen, temLen, count, ln, it, stride, val;
  int val1, val2, stride1, stride2, row, col, total, *pLen; 
  double *f1, *f2, *mr, *mi;
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
      double *output1_r,*output1_i,*output2;

      
  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  dualtree2D_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
  {
    validate_print (errCode);
    return 0;			
  }

//   GetRhsVar(1, "d", &m1, &n1, &l1);
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
			return -1;
		}
//   GetRhsVar(2, "i", &m2, &n2, &l2);
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
//   GetRhsVar(3, "d", &m3, &n3, &l3);
    		_SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(4, "d", &m4, &n4, &l4);
    		_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}


  wave_len_validate (m1, n3, &stride1, &val1);
  wave_len_validate (n1, n3, &stride2, &val2);
  if ((val1 == 0) || (val2 == 0))
    {
      sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
      return 0;
    }
  stride = (stride1 > stride2) ? stride2 : stride1;
  if ((input2[0] < 1) || (input2[0] > stride))
    {
      sciprint ("Level Parameter is not valid for input matrix!\n");
      return 0;
    }
  pLen = (int *)malloc ((input2[0] + 2) * 2 * sizeof (int));
  matrix_wavedec_len_cal (m1, n1, (int)input2[0], n3, pLen);
  wave_mem_cal (pLen, (int)input2[0], &total);

  f1 = (double *)malloc(m3*n3*sizeof(double));
  f2 = (double *)malloc(m3*n3*sizeof(double));
  
  matrix_tran(input3,n3,m3,f1,m3,n3);
  matrix_tran(input4,n3,m3,f2,m3,n3);

  it = 1;
  m5 = 1;
  n5 = total;
  m6 = input2[0] + 2;
  n6 = 2;
//   CreateCVar (5, "d",&it, &m5, &n5, &l5r, &l5c);
  _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output1_r,&output1_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar(6, "i", &m6, &n6, &l6);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 2, m6, n6, &output2);
	//_SciErr = allocMatrixOfInteger32(pvApiCtx, Rhs + 2, m6, n6, &ptr2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

	
  mr = (double *)malloc(m5*n5*sizeof(double));
  mi = (double *)malloc(m5*n5*sizeof(double));
  for (row = 0; row < m6; row++)
    {
      for (col = 0; col < n6; col++)
	output2[row + col * m6] = (double)pLen[col + row * n6];
    }

  cowavedec2 (input1, m1, n1, f1, f1+n3, f2, f2+n3,
	      n3, pLen, mr, total, (int)input2[0], dwtMode);

  cowavedec2 (input1, m1, n1, f1+n3*2, f1+n3*3, f2+n3*2, f2+n3*3,
	      n3, pLen, mi, total, (int)input2[0], dwtMode);

  copmd (mr,mi,total,pLen[0],pLen[1],output1_r,output1_i);

  
  free(pLen);
  free(f1);
  free(f2);
  free(mr);
  free(mi);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;

  return 0;
}


int
int_idualtree2D(char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l1r, l1c;
   int minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
  int errCode, flow, it, val, *pLen, size, row, col, i;
  double *f1, *f2, *mr, *mi, *maxR, *maxI;
  int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
  int * p_input_vector4 = NULL;
  double *input1_r,*input1_i;
  double *input2;
  double *input3;
  double *input4;
  SciErr _SciErr;
  int type;
   double *output1;
   
  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  idualtree2D_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
  {
    validate_print (errCode);
    return 0;			
  }

  it = 1;

//   GetRhsCVar(1, "d", &it, &m1, &n1, &l1r, &l1c);
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
		_SciErr = getComplexMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1_r,&input1_i);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(2, "i", &m2, &n2, &l2);
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
//   GetRhsVar(3, "d", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(4, "d", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  if ((input2[0] < n3) || (input2[m2] < n3))
    {
      sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
      return 0;
    }
  size = 0;
  for (row = 1; row < (m2 - 1); row++)
    {
      size += (input2[row]) * (input2[row + m2]);
    }
  size = size * 3 + (input2[0]) * (input2[m2]);
  if (m1 * n1 != size)
    {
      sciprint ("Inputs are not size array and coefs!\n");
      return 0;
    }
  val = 0;
  if ((input2[0] != input2[1]) ||
      (input2[m2] != input2[m2+1]))
    val += 1;
  for (row = 1; row < (m2 - 1); row++)
    {
      if (input2[row] >= input2[row + 1])
	val += 1;
      if (input2[row + m2] >= input2[row + m2 + 1])
	val += 1;
    }
  if (val != 0)
    {
      sciprint ("Inputs are not size array!\n");
      return 0;
    }
  pLen = (int *)malloc (m2 * n2 * sizeof (int));
  for (row = 0; row < n2; row++)
    {
      for (col = 0; col < m2; col++)
	pLen[row + col * n2] = (int)input2[col + row * m2];
    }

  m5 = pLen[(m2 - 1) * n2];
  n5 = pLen[(m2 - 1) * n2 + 1];

//   CreateVar(5, "d", &m5, &n5, &l5);
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  f1 = (double *)malloc(m3*n3*sizeof(double));
  f2 = (double *)malloc(m3*n3*sizeof(double));
  matrix_tran(input3,n3,m3,f1,m3,n3);
  matrix_tran(input4,n3,m3,f2,m3,n3);
  mr = (double *)malloc(m1*n1*sizeof(double));
  mi = (double *)malloc(m1*n1*sizeof(double));
  maxR = (double *)malloc(m5*n5*sizeof(double));
  maxI = (double *)malloc(m5*n5*sizeof(double));

  copmr (input1_r,input1_i,m1*n1,pLen[0],pLen[1],mr,mi);

  cowaverec2 (mr, m1*n1, f1, f1+n3, f2, f2+n3,
        n3, maxR, m5, n5, pLen, m2-2, dwtMode);
  cowaverec2 (mi, m1*n1, f1+n3*2, f1+n3*3, f2+n3*2, f2+n3*3,
        n3, maxI, m5, n5, pLen, m2-2, dwtMode);

  for(i=0;i<m5*n5;i++)
    output1[i] = (maxR[i]+maxI[i])/2;

  free(pLen);
  free(mr);
  free(mi);
  free(maxR);
  free(maxI);
  free(f1);
  free(f2);
  LhsVar(1) = Rhs + 1;

  return 0;
}


int
int_cplxdual2D(char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l6, m6, n6, l7, m7, n7, l8, m8, n8;
   int l5r, l5c, l7r, l7c;
   int minlhs=3, maxlhs=3, minrhs=4, maxrhs=4;
  int errCode, flow, calLen, temLen, count, ln, it, stride, val;
  int val1, val2, stride1, stride2, row, col, total, *pLen; 
  double *f1, *f2, *mr, *mi, *mrr, *mii;
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
      double *output1_r,*output1_i,*output2;
      int *output2_int;
      double *output3_r,*output3_i;
      int i,j;
  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  dualtree2D_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
  {
    validate_print (errCode);
    return 0;			
  }

//   GetRhsVar(1, "d", &m1, &n1, &l1);
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
			return -1;
		}
//   GetRhsVar(2, "i", &m2, &n2, &l2);
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
//   GetRhsVar(3, "d", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 3. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector3, &m3, &n3, &input3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(4, "d", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  wave_len_validate (m1, n3, &stride1, &val1);
  wave_len_validate (n1, n3, &stride2, &val2);
  if ((val1 == 0) || (val2 == 0))
    {
      sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
      return 0;
    }
  stride = (stride1 > stride2) ? stride2 : stride1;
  if ((input2[0] < 1) || (input2[0] > stride))
    {
      sciprint ("Level Parameter is not valid for input matrix!\n");
      return 0;
    }
  pLen = malloc ((input2[0] + 2) * 2 * sizeof (int));
  matrix_wavedec_len_cal (m1, n1, input2[0], n3, pLen);
  wave_mem_cal (pLen, input2[0], &total);

  f1 = (double *)malloc(m3*n3*sizeof(double));
  f2 = (double *)malloc(m3*n3*sizeof(double));
  
  matrix_tran(input3,n3,m3,f1,m3,n3);
  matrix_tran(input4,n3,m3,f2,m3,n3);

  it = 1;
  m5 = 1;
  n5 = total;
  m6 = input2[0] + 2;
  n6 = 2;
  m7 = m5;
  n7 = n5;
//   CreateCVar (5, "d",&it, &m5, &n5, &l5r, &l5c);
    _SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 1, m5, n5, &output1_r,&output1_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
	_SciErr = allocComplexMatrixOfDouble(pvApiCtx, Rhs + 2, m7, n7, &output3_r,&output3_i);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}
//   CreateVar(6, "i", &m6, &n6, &l6);
_SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 3, m6, n6, &output2);
	//_SciErr = allocMatrixOfInteger32(pvApiCtx, Rhs + 2, m6, n6, &ptr2);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  mr = (double *)malloc(m5*n5*sizeof(double));
  mi = (double *)malloc(m5*n5*sizeof(double));
//   CreateCVar(7, "d", &it, &m7, &n7, &l7r, &l7c);

  mrr = (double *)malloc(m5*n5*sizeof(double));
  mii = (double *)malloc(m5*n5*sizeof(double));
  for (row = 0; row < m6; row++)
    {
      for (col = 0; col < n6; col++)
	output2[row + col * m6] = (double)pLen[col + row * n6];
    }

  /* (1,1) */
  cowavedec2a (input1, m1, n1, f1, f1+n3, f1, f1+n3, f2, f2+n3,
	       f2, f2+n3, n3, pLen, mr, total, (int)input2[0], dwtMode);
  /* (2,2) */
  cowavedec2a (input1, m1, n1, f1+n3*2, f1+n3*3, f1+n3*2, f1+n3*3,
	       f2+n3*2, f2+n3*3, f2+n3*2, f2+n3*3,
	       n3, pLen, mi, total,  (int)input2[0], dwtMode);
  /* (1,2) */
  cowavedec2a (input1, m1, n1, f1, f1+n3, f1+n3*2, f1+n3*3, 
	       f2, f2+n3, f2+n3*2, f2+n3*3, n3, pLen, mrr, total, 
	       (int)input2[0], dwtMode);
  /* (2,1) */
  cowavedec2a (input1, m1, n1, f1+n3*2, f1+n3*3, f1, f1+n3,  
	       f2+n3*2, f2+n3*3, f2, f2+n3, n3, pLen, mii, total, 
	       (int)input2[0], dwtMode);


  copmd (mr,mi,total,pLen[0],pLen[1],output1_r,output1_i);

  copmd (mrr,mii,total,pLen[0],pLen[1],output3_r,output3_i);

	
  free(pLen);
  free(f1);
  free(f2);
  free(mr);
  free(mi);
  free(mrr);
  free(mii);

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  LhsVar(3) = Rhs + 3;

  return 0;
}

int
int_icplxdual2D(char *fname)
{
   int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
   int l5, m5, n5, l1r, l1c, l6, m6, n6, l2r, l2c;
   int minlhs=1, maxlhs=1, minrhs=5, maxrhs=5;
  int errCode, flow, it, val, *pLen, size, row, col, i;
  double *f1, *f2, *mr, *mi, *maxR, *maxI, *mrr, *mii, *maxRR, *maxII;
 int * p_input_vector1 = NULL;
  int * p_input_vector2 = NULL;
  int * p_input_vector3 = NULL;
  int * p_input_vector4 = NULL;
  int * p_input_vector5 = NULL;
  double *input1_r,*input1_i;
  double *input2_r,*input2_i;
  double *input3;
  double *input4;
    double *input5;
  SciErr _SciErr;
  int type;
   double *output1;
  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  icplxdual2D_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
  {
    validate_print (errCode);
    return 0;			
  }

  it = 1;

//   GetRhsCVar(1, "d", &it, &m1, &n1, &l1r, &l1c);
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
		_SciErr = getComplexMatrixOfDouble(pvApiCtx, p_input_vector1, &m1, &n1, &input1_r,&input1_i);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsCVar(2, "d", &it, &m2, &n2, &l2r, &l2c);
 _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_input_vector2);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return 0;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector2, &type);
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
		_SciErr = getComplexMatrixOfDouble(pvApiCtx, p_input_vector2, &m2, &n2, &input2_r,&input2_i);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(3, "i", &m3, &n3, &l3);
_SciErr = getVarAddressFromPosition(pvApiCtx,3, &p_input_vector3);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector3, &type);
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
//   GetRhsVar(4, "d", &m4, &n4, &l4);
_SciErr = getVarAddressFromPosition(pvApiCtx,4, &p_input_vector4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector4, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 4. input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector4, &m4, &n4, &input4);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
//   GetRhsVar(5, "d", &m5, &n5, &l5);
_SciErr = getVarAddressFromPosition(pvApiCtx,5, &p_input_vector5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
                _SciErr = getVarType(pvApiCtx, p_input_vector5, &type);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}
		if (type!=sci_matrix)
		{
		  Scierror (999,"%s: 5 input vector must be double\n",fname);	
		  return 0;
		}
		_SciErr = getMatrixOfDouble(pvApiCtx, p_input_vector5, &m5, &n5, &input5);
		if(_SciErr.iErr)
		{
			printError(&_SciErr, 0);
			return -1;
		}

  if ((input3[0] < n4) || (input3[m3] < n4))
    {
      sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
      return 0;
    }
  size = 0;
  for (row = 1; row < (m3 - 1); row++)
    {
      size += (input3[row]) * (input3[row + m3]);
    }
  size = size * 3 + (input3[0]) * (input3[m3]);
  if ((m1 * n1 != size) || (m2 * n2 != size))
    {
      sciprint ("Inputs are not size array and coefs!\n");
      return 0;
    }
  val = 0;
  if ((input3[0] != input3[1]) ||
      (input3[m3] != input3[m3+1]))
    val += 1;
  for (row = 1; row < (m3 - 1); row++)
    {
      if (input3[row] >= input3[row + 1])
	val += 1;
      if (input3[row + m3] >= input3[row + m3 + 1])
	val += 1;
    }
  if (val != 0)
    {
      sciprint ("Inputs are not size array!\n");
      return 0;
    }
  pLen = (int *)malloc (m3 * n3 * sizeof (int));
  for (row = 0; row < n3; row++)
    {
      for (col = 0; col < m3; col++)
	pLen[row + col * n3] = (int)input3[col + row * m3];
    }

  m6 = pLen[(m3 - 1) * n3];
  n6 = pLen[(m3 - 1) * n3 + 1];

//   CreateVar(6, "d", &m6, &n6, &l6);
    _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m6, n6, &output1);
	if(_SciErr.iErr)
	 {
		printError(&_SciErr, 0);
		return -1;
	}

  f1 = (double *)malloc(m4*n4*sizeof(double));
  f2 = (double *)malloc(m4*n4*sizeof(double));
  matrix_tran(input4,n4,m4,f1,m4,n4);
  matrix_tran(input5,n4,m4,f2,m4,n4);
  mr = (double *)malloc(m1*n1*sizeof(double));
  mi = (double *)malloc(m1*n1*sizeof(double));
  mrr = (double *)malloc(m1*n1*sizeof(double));
  mii = (double *)malloc(m1*n1*sizeof(double));
  maxR = (double *)malloc(m6*n6*sizeof(double));
  maxI = (double *)malloc(m6*n6*sizeof(double));
  maxRR = (double *)malloc(m6*n6*sizeof(double));
  maxII = (double *)malloc(m6*n6*sizeof(double));

  copmr (input1_r,input1_i,m1*n1,pLen[0],pLen[1],mr,mi);
  copmr (input2_r,input2_i,m1*n1,pLen[0],pLen[1],mrr,mii);

  /* (1,1) */
  cowaverec2a (mr, m1*n1, f1, f1+n4, f1, f1+n4, 
	       f2, f2+n4, f2, f2+n4, n4, maxR, 
	       m6, n6, pLen, m3-2, dwtMode);
  /* (2,2) */
  cowaverec2a (mi, m1*n1, f1+n4*2, f1+n4*3, 
	       f1+n4*2, f1+n4*3, f2+n4*2, f2+n4*3,
	       f2+n4*2, f2+n4*3, n4, maxI, m6, n6, 
	       pLen, m3-2, dwtMode);
  /* (1,2) */
  cowaverec2a (mrr, m1*n1, f1, f1+n4, f1+n4*2, f1+n4*3,
	       f2, f2+n4, f2+n4*2, f2+n4*3,
	       n4, maxRR, m6, n6, pLen, m3-2, dwtMode);
  /* (2,2) */
  cowaverec2a (mii, m1*n1, f1+n4*2, f1+n4*3, f1, f1+n4,
	       f2+n4*2, f2+n4*3, f2, f2+n4,
	       n4, maxII, m6, n6, pLen, m3-2, dwtMode);

  for(i=0;i<m6*n6;i++)
    output1[i] = (maxR[i]+maxI[i]+maxRR[i]+maxII[i])/4;

  free(pLen);
  free(mr);
  free(mi);
  free(mrr);
  free(mii);
  free(maxR);
  free(maxI);
  free(maxRR);
  free(maxII);
  free(f1);
  free(f2);
  LhsVar(1) = Rhs + 1;

  return 0;
}
