/*
 * -------------------------------------------------------------------------
 * dwt2d_int.c -- 2-D signal decomposition and reconstruction interface
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
#include "stack-c.h"

int
int_dwt2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7, l8, m8, n8;
  static int l9, m9, n9;
  static int minrhs=2, maxrhs=5, minlhs=4, maxlhs=4;
  int errCode, flow, family, member, ii; 
  int stride1, val1, stride2, val2;
  Func ana_fun;
  swt_wavelet pWaveStruct;
  extend_method extMethod;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  dwt2_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
    {
	  //sciprint("enter here!\n");
      validate_print (errCode);
      return 0;			
    }

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;
  l5 = 0;

  switch (flow) {
  case 1:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "c", &m2, &n2, &l2);
      dwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l2),&family,&member);
      wavelet_fun_parser (cstk(l2), &ii);
      ana_fun = wi[ii].analysis;
      (*ana_fun)(member, &pWaveStruct);
      wave_len_validate (m1, pWaveStruct.length, &stride1, &val1);
      wave_len_validate (n1, pWaveStruct.length, &stride2, &val2);
      if ((val1 == 0) || (val2 == 0))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m3 = (int)floor((m1 + pWaveStruct.length - 1)/2);
      n3 = (int)floor((n1 + pWaveStruct.length - 1)/2);
	  if (dwtMode==PER)
	  {
		  m3 = (int)ceil(((double)m1)/2.0);
		  n3 = (int)ceil(((double)n1)/2.0);
	  }
      m4 = m3;
      n4 = n3;
      m5 = m3;
      n5 = n3;
      m6 = m3;
      n6 = n3;
      CreateVar(3, "d", &m3, &n3, &l3);
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      dwt2D_neo (stk(l1), m1, n1, stk(l3), stk(l4), stk(l5), stk(l6),
	     m3, n3, pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
	     pWaveStruct.length, dwtMode);
      LhsVar(1) = 3;
      LhsVar(2) = 4;
      LhsVar(3) = 5;
      LhsVar(4) = 6;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      dwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((m1<m2*n2) || (n1<m2*n2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m4 = (int)floor((m1 + m2*n2 - 1)/2);
      n4 = (int)floor((n1 + m2*n2 - 1)/2);
	  if (dwtMode==PER)
	  {
		  m4 = (int)ceil(((double)m1)/2.0);
		  n4 = (int)ceil(((double)n1)/2.0);
	  }
      m5 = m4;
      m6 = m4;
      m7 = m4;
      n5 = n4;
      n6 = n4;
      n7 = n4;
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      CreateVar(7, "d", &m7, &n7, &l7);
      dwt2D_neo (stk(l1), m1, n1, stk(l4), stk(l5), stk(l6), stk(l7),
	     m4, n4, stk(l2), stk(l3), m2*n2, dwtMode);
      LhsVar(1) = 4;
      LhsVar(2) = 5;
      LhsVar(3) = 6;
      LhsVar(4) = 7;
      break;
    }
  case 3:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "c", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "c", &m4, &n4, &l4);
      dwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l2),&family,&member);
      wavelet_fun_parser (cstk(l2), &ii);
      ana_fun = wi[ii].analysis;
      (*ana_fun)(member, &pWaveStruct);
      wave_len_validate (m1, pWaveStruct.length, &stride1, &val1);
      wave_len_validate (n1, pWaveStruct.length, &stride2, &val2);
      if ((val1 == 0) || (val2 == 0))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      //dwt_write(cstk(l4),&errCode);
      extend_method_parse (cstk(l4), &extMethod);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      m5 = (int)floor((m1 + pWaveStruct.length - 1)/2);
      n5 = (int)floor((n1 + pWaveStruct.length - 1)/2);
	  if (extMethod==PER)
	  {
		  m5 = (int)ceil(((double)m1)/2.0);
		  n5 = (int)ceil(((double)n1)/2.0);
	  }
      m6 = m5;
      m7 = m5;
      m8 = m5;
      n6 = n5;
      n7 = n5;
      n8 = n5;
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      CreateVar(7, "d", &m7, &n7, &l7);
      CreateVar(8, "d", &m8, &n8, &l8);
      dwt2D_neo (stk(l1), m1, n1, stk(l5), stk(l6), stk(l7), stk(l8),
	     m5, n5, pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
	     pWaveStruct.length, extMethod);
      LhsVar(1) = 5;
      LhsVar(2) = 6;
      LhsVar(3) = 7;
      LhsVar(4) = 8;
      filter_clear();
      break;
    }
  case 4:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "c", &m4, &n4, &l4);
      GetRhsVar(5, "c", &m5, &n5, &l5);
      dwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((m1<m2*n2) || (n1<m2*n2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      //dwt_write(cstk(l5),&errCode);
      extend_method_parse (cstk(l5), &extMethod);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      m6 = (int)floor((m1 + m2*n2 - 1)/2);
      n6 = (int)floor((n1 + m2*n2 - 1)/2);
	  if (extMethod==PER)
	  {
		  m6 = (int)ceil(((double)m1)/2.0);
		  n6 = (int)ceil(((double)n1)/2.0);
	  }
      m7 = m6;
      m8 = m6;
      m9 = m6;
      n7 = n6;
      n8 = n6;
      n9 = n6;
      CreateVar(6, "d", &m6, &n6, &l6);
      CreateVar(7, "d", &m7, &n7, &l7);
      CreateVar(8, "d", &m8, &n8, &l8);
      CreateVar(9, "d", &m9, &n9, &l9);
      dwt2D_neo (stk(l1), m1, n1, stk(l6), stk(l7), stk(l8), stk(l9),
	     m6, n6, stk(l2), stk(l3), m2*n2, extMethod);
      LhsVar(1) = 6;
      LhsVar(2) = 7;
      LhsVar(3) = 8;
      LhsVar(4) = 9;
      break;
    }
  default:
    break;
  }
  return 0;
}

int
int_idwt2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7, l8, m8, n8;
  static int l9, m9, n9, l10, m10, n10;
  static int minrhs=5, maxrhs=9, minlhs=1, maxlhs=1;
  int errCode, flow, family, member, ii, stride, val1, val2;
  int m, n, count;
  double *ca, *ch, *cv, *cd, *vo;
  Func syn_fun;
  swt_wavelet pWaveStruct;
  extend_method extMethod;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  idwt2_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;
  l5 = 0;
  l6 = 0;
  l7 = 0;
  l8 = 0;
  l9 = 0;

  GetRhsVar(1, "d", &m1, &n1, &l1);
  GetRhsVar(2, "d", &m2, &n2, &l2);
  GetRhsVar(3, "d", &m3, &n3, &l3);
  GetRhsVar(4, "d", &m4, &n4, &l4);

  m = m1>0?m1:(m2>0?m2:(m3>0?m3:m4));
  n = n1>0?n1:(n2>0?n2:(n3>0?n3:n4));

  if (!m1 || !m2 || !m3 || !m4)
    {
      vo = malloc(m*n*sizeof(double));
      for(count=0;count<m*n;count++)
	vo[count] = 0;
    }
  if (!m1)
    ca = vo;
  else
    ca = stk(l1);
  if (!m2)
    ch = vo;
  else
    ch = stk(l2);
  if (!m3)
    cv = vo;
  else
    cv = stk(l3);
  if (!m4)
    cd = vo;
  else
    cd = stk(l4);

  switch (flow) {
  case 1:
    {
      GetRhsVar(5, "c", &m5, &n5, &l5);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l5),&family,&member);
      wavelet_fun_parser (cstk(l5), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate (m, pWaveStruct.length, &stride, &val1);
      wave_len_validate (n, pWaveStruct.length, &stride, &val2);
      if ((!val1) || (!val2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m6 = m*2 - pWaveStruct.length + 2;
      n6 = n*2 - pWaveStruct.length + 2;
	  if (dwtMode==PER)
	  {
        m6 = 2*m;
		n6 = 2*n;
	  }
      CreateVar(6, "d", &m6, &n6, &l6);
      //idwt2D (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	    //  pWaveStruct.pHiPass, pWaveStruct.length, 
	      //stk(l6), m6, n6, dwtMode);
	  idwt2D_neo (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, 
	      stk(l6), m6, n6);
      LhsVar(1) = 6;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(5, "d", &m5, &n5, &l5);
      GetRhsVar(6, "d", &m6, &n6, &l6);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((m<m5*n5) || (n<m5*n5))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m7 = m*2 - m5*n5 + 2;
      n7 = n*2 - m5*n5 + 2;
	  if (dwtMode==PER)
	  {
        m7 = 2*m;
		n7 = 2*n;
	  }
      CreateVar(7, "d", &m7, &n7, &l7);
      //idwt2D (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	    //  m5*n5, stk(l7), m7, n7, dwtMode);
	  idwt2D_neo (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	      m5*n5, stk(l7), m7, n7);
      LhsVar(1) = 7;
      break;
    }
  case 3:
    {
      GetRhsVar(5, "c", &m5, &n5, &l5);
      GetRhsVar(6, "i", &m6, &n6, &l6);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l5),&family,&member);
      wavelet_fun_parser (cstk(l5), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate (m, pWaveStruct.length, &stride, &val1);
      wave_len_validate (n, pWaveStruct.length, &stride, &val2);
      if ((!val1) || (!val2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      /*if ((istk(l6)[0]>(2*m-pWaveStruct.length+2)) ||
	  (istk(l6)[1]>(2*n-pWaveStruct.length+2)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}*/

	  if ((istk(l6)[0]>(2*m+pWaveStruct.length)) ||
	  (istk(l6)[1]>(2*n+pWaveStruct.length)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}
      m7 = istk(l6)[0];
      n7 = istk(l6)[1];
	  //sciprint("m=%d\n",m7);
      CreateVar(7, "d", &m7, &n7, &l7);
      //idwt2D (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	    //  pWaveStruct.pHiPass, pWaveStruct.length, 
	      //stk(l7), m7, n7, dwtMode);
	  idwt2D_neo (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, 
	      stk(l7), m7, n7);
      LhsVar(1) = 7;
      filter_clear();
      break;
    }
  case 4:
    {
      GetRhsVar(5, "d", &m5, &n5, &l5);
      GetRhsVar(6, "d", &m6, &n6, &l6);
      GetRhsVar(7, "i", &m7, &n7, &l7);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((m<m5*n5) || (n<m5*n5))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      /*if ((istk(l7)[0]>(2*m-m5*n5+2)) ||
	  (istk(l7)[1]>(2*n-m5*n5+2)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}*/

	  if ((istk(l7)[0]>(2*m+m5*n5)) ||
	  (istk(l7)[1]>(2*n+m5*n5)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}
      m8 = istk(l7)[0];
      n8 = istk(l7)[1];
      CreateVar(8, "d", &m8, &n8, &l8);
      //idwt2D (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	    //  m5*n5, stk(l8), m8, n8, dwtMode);
	  idwt2D_neo (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	      m5*n5, stk(l8), m8, n8);
      LhsVar(1) = 8;
      break;
    }
  case 5:
    {
      GetRhsVar(5, "c", &m5, &n5, &l5);
      GetRhsVar(6, "c", &m6, &n6, &l6);
      GetRhsVar(7, "c", &m7, &n7, &l7);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l5),&family,&member);
      wavelet_fun_parser (cstk(l5), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate (m, pWaveStruct.length, &stride, &val1);
      wave_len_validate (n, pWaveStruct.length, &stride, &val2);
      if ((!val1) || (!val2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      extend_method_parse (cstk(l7), &extMethod);
      m8 = m*2 - pWaveStruct.length + 2;
      n8 = n*2 - pWaveStruct.length + 2;
	  if (extMethod==PER)
	  {
        m8 = 2*m;
		n8 = 2*n;
	  }
      CreateVar(8, "d", &m8, &n8, &l8);
      //idwt2D (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	    //  pWaveStruct.pHiPass, pWaveStruct.length, 
	      //stk(l8), m8, n8, extMethod);
	  idwt2D_neo (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, 
	      stk(l8), m8, n8);
      LhsVar(1) = 8;
      filter_clear();
      break;
    }
  case 6:
    {
      GetRhsVar(5, "d", &m5, &n5, &l5);
      GetRhsVar(6, "d", &m6, &n6, &l6);
      GetRhsVar(7, "c", &m7, &n7, &l7);
      GetRhsVar(8, "c", &m8, &n8, &l8);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((m<m5*n5) || (n<m5*n5))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      extend_method_parse (cstk(l8), &extMethod);
      m9 = m*2 - m5*n5 + 2;
      n9 = n*2 - m5*n5 + 2;
	  if (extMethod==PER)
	  {
        m9 = 2*m;
		n9 = 2*n;
	  }
      CreateVar(9, "d", &m9, &n9, &l9);
      //idwt2D (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	    //  m5*n5, stk(l9), m9, n9, extMethod);
	  idwt2D_neo (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	      m5*n5, stk(l9), m9, n9);
      LhsVar(1) = 9;
      break;
    }
  case 7:
    {
      GetRhsVar(5, "c", &m5, &n5, &l5);
      GetRhsVar(6, "i", &m6, &n6, &l6);
      GetRhsVar(7, "c", &m7, &n7, &l7);
      GetRhsVar(8, "c", &m8, &n8, &l8);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l5),&family,&member);
      wavelet_fun_parser (cstk(l5), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate (m, pWaveStruct.length, &stride, &val1);
      wave_len_validate (n, pWaveStruct.length, &stride, &val2);
      if ((!val1) || (!val2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      /*if ((istk(l6)[0]>(2*m-pWaveStruct.length+2)) ||
	  (istk(l6)[1]>(2*n-pWaveStruct.length+2)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}*/

	  if ((istk(l6)[0]>(2*m+pWaveStruct.length)) ||
	  (istk(l6)[1]>(2*n+pWaveStruct.length)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}
      extend_method_parse (cstk(l8), &extMethod);
      m9 = istk(l6)[0];
      n9 = istk(l6)[1];
      CreateVar(9, "d", &m9, &n9, &l9);
      //idwt2D (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	    //  pWaveStruct.pHiPass, pWaveStruct.length, 
	      //stk(l9), m9, n9, extMethod);
	  idwt2D_neo (ca, ch, cv, cd, m, n, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, 
	      stk(l9), m9, n9);//, extMethod);
      LhsVar(1) = 9;
      filter_clear();
      break;
    }
  case 8:
    {
      GetRhsVar(5, "d", &m5, &n5, &l5);
      GetRhsVar(6, "d", &m6, &n6, &l6);
      GetRhsVar(7, "i", &m7, &n7, &l7);
      GetRhsVar(8, "c", &m8, &n8, &l8);
      GetRhsVar(9, "c", &m9, &n9, &l9);
      idwt2_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7,l8,l9);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((m<m5*n5) || (n<m5*n5))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      /*if ((istk(l7)[0]>(2*m-m5*n5+2)) ||
	  (istk(l7)[1]>(2*n-m5*n5+2)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}*/

	  if ((istk(l7)[0]>(2*m+m5*n5)) ||
	  (istk(l7)[1]>(2*n+m5*n5)))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}
      extend_method_parse (cstk(l9), &extMethod);
      m10 = istk(l7)[0];
      n10 = istk(l7)[1];
      CreateVar(10, "d", &m10, &n10, &l10);
      //idwt2D (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	    //  m5*n5, stk(l10), m10, n10, extMethod);
	  idwt2D_neo (ca, ch, cv, cd, m, n, stk(l5), stk(l6),
	      m5*n5, stk(l10), m10, n10);
      LhsVar(1) = 10;
      break;
    }
  default:
    break;
  }

  //sciprint("flow = %d\n",flow);
  return 0;
}

int
int_wavedec2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6;
  static int minrhs=3, maxrhs=4, minlhs=2, maxlhs=2;
  int errCode, flow, family, member, ii, row, col, total, *pLen; 
  int stride1, val1, stride2, val2, stride;
  Func ana_fun;
  swt_wavelet pWaveStruct;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  //if (dwtMode==PER)
  //{
    //sciprint("PER extension mode is not supported for multiple level decomposition!\n");
	//sciprint("Please change the extension mode by dwtmode command.\n");
	//return 0;
  //}

  wavedec2_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;

  switch (flow) {
  case 1:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      wavedec2_content_validate(&errCode,flow,l1,l2,l3,l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      ana_fun = wi[ii].analysis;
      (*ana_fun)(member, &pWaveStruct);
      wave_len_validate (m1, pWaveStruct.length, &stride1, &val1);
      wave_len_validate (n1, pWaveStruct.length, &stride2, &val2);
      if ((val1 == 0) || (val2 == 0))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      stride = (stride1 > stride2) ? stride2 : stride1;
      if ((istk (l2)[0] < 1) || (istk (l2)[0] > stride))
	{
	  sciprint ("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      
      pLen = malloc ((istk (l2)[0] + 2) * 2 * sizeof (int));
      matrix_wavedec_len_cal (m1, n1, istk (l2)[0], 
			      pWaveStruct.length, pLen);
      //for(row=0;row<(istk(l2)[0]+2)*2;row++)
      //sciprint("%d\n",pLen[row]);

      wave_mem_cal (pLen, istk (l2)[0], &total);
      m4 = 1;
      n4 = total;
      m5 = istk (l2)[0] + 2;
      n5 = 2;
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "i", &m5, &n5, &l5);
      for (row = 0; row < m5; row++)
	{
	  for (col = 0; col < n5; col++)
	    istk (l5)[row + col * m5] = pLen[col + row * n5];
	}
      wavedec2 (stk(l1), m1, n1, pWaveStruct.pLowPass, 
		pWaveStruct.pHiPass, pWaveStruct.length,
		pLen, stk(l4), m4*n4, istk(l2)[0], dwtMode);
      LhsVar(1) = 4;
      LhsVar(2) = 5;
      filter_clear();
      free(pLen);
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      wavedec2_content_validate(&errCode,flow,l1,l2,l3,l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate (m1, m3*n3, &stride1, &val1);
      wave_len_validate (n1, m3*n3, &stride2, &val2);
      if ((val1 == 0) || (val2 == 0))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      stride = (stride1 > stride2) ? stride2 : stride1;
      if ((istk (l2)[0] < 1) || (istk (l2)[0] > stride))
	{
	  sciprint ("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      pLen = malloc ((istk (l2)[0] + 2) * 2 * sizeof (int));
      matrix_wavedec_len_cal (m1, n1, istk (l2)[0], m3*n3, pLen);
      wave_mem_cal (pLen, istk (l2)[0], &total);
      m5 = 1;
      n5 = total;
      m6 = istk (l2)[0] + 2;
      n6 = 2;
      CreateVar (5, "d", &m5, &n5, &l5);
      CreateVar (6, "i", &m6, &n6, &l6);
      for (row = 0; row < m6; row++)
	{
	  for (col = 0; col < n6; col++)
	    istk (l6)[row + col * m6] = pLen[col + row * n6];
	}
      wavedec2 (stk(l1), m1, n1, stk(l3), stk(l4), m3*n3,
		pLen, stk(l5), m5*n5, istk(l2)[0], dwtMode);
      LhsVar (1) = 5;
      LhsVar (2) = 6;
      free(pLen);
      break;
    }
  default:
    break;
  }
  return 0;
}

int
int_waverec2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5;
  static int minrhs=3, maxrhs=4, minlhs=1, maxlhs=1;
  int errCode, flow, family, member, ii, val, row, col, size, *pLen;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  waverec2_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;

  switch (flow) {
  case 1:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      waverec2_content_validate (&errCode, flow, l1, l2, l3, l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if ((istk(l2)[0] < pWaveStruct.length) || 
	  (istk(l2)[m2] < pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      size = 0;
      for (row = 1; row < (m2 - 1); row++)
	{
	  size += (istk (l2)[row]) * (istk (l2)[row + m2]);
	}
      size = size * 3 + (istk (l2)[0]) * (istk (l2)[m2]);
      if (m1 * n1 != size)
	{
	  sciprint ("Inputs are not size array and coefs!\n");
	  return 0;
	}
      val = 0;
      if ((istk(l2)[0] != istk(l2)[1]) ||
	  (istk(l2)[m2] != istk(l2)[m2+1]))
	val += 1;
      for (row = 1; row < (m2 - 1); row++)
	{
	  if (istk (l2)[row] >= istk (l2)[row + 1])
	    val += 1;
	  if (istk (l2)[row + m2] >= istk (l2)[row + m2 + 1])
	    val += 1;
	}
      if (val != 0)
	{
	  sciprint ("Inputs are not size array!\n");
	  return 0;
	}
      pLen = malloc (m2 * n2 * sizeof (int));
      for (row = 0; row < n2; row++)
	{
	  for (col = 0; col < m2; col++)
	    pLen[row + col * n2] = istk (l2)[col + row * m2];
	}
      m4 = pLen[(m2 - 1) * n2];
      n4 = pLen[(m2 - 1) * n2 + 1];
      CreateVar (4, "d", &m4, &n4, &l4);
      waverec2 (stk(l1), m1*n1, pWaveStruct.pLowPass,
		pWaveStruct.pHiPass, pWaveStruct.length, 
		stk(l4), m4, n4, pLen, m2-2, dwtMode);
      LhsVar(1) = 4;
      filter_clear();
      free(pLen);
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      waverec2_content_validate (&errCode, flow, l1, l2, l3, l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l2)[0] < m3*n3) || 
	  (istk(l2)[m2] < m3*n3))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      size = 0;
      for (row = 1; row < (m2 - 1); row++)
	{
	  size += (istk (l2)[row]) * (istk (l2)[row + m2]);
	}
      size = size * 3 + (istk (l2)[0]) * (istk (l2)[m2]);
      if (m1 * n1 != size)
	{
	  sciprint ("Inputs are not size array and coefs!\n");
	  return 0;
	}
      val = 0;
      if ((istk(l2)[0] != istk(l2)[1]) ||
	  (istk(l2)[m2] != istk(l2)[m2+1]))
	val += 1;
      for (row = 1; row < (m2 - 1); row++)
	{
	  if (istk (l2)[row] >= istk (l2)[row + 1])
	    val += 1;
	  if (istk (l2)[row + m2] >= istk (l2)[row + m2 + 1])
	    val += 1;
	}
      if (val != 0)
	{
	  sciprint ("Inputs are not size array!\n");
	  return 0;
	}
      pLen = malloc (m2 * n2 * sizeof (int));
      for (row = 0; row < n2; row++)
	{
	  for (col = 0; col < m2; col++)
	    pLen[row + col * n2] = istk (l2)[col + row * m2];
	}
      m5 = pLen[(m2 - 1) * n2];
      n5 = pLen[(m2 - 1) * n2 + 1];
      CreateVar (5, "d", &m5, &n5, &l5);
      waverec2 (stk(l1), m1*n1, stk(l3), stk(l4),
		m3*n3, stk(l5), m5, n5, pLen, m2-2, dwtMode);
      LhsVar(1) = 5;
      free(pLen);
      break;
    }
  default:
    break;
  }

  return 0;
}

int
int_wenergy2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6;
  static int minrhs=2, maxrhs=2, minlhs=2, maxlhs=4;
  int errCode, flow, size, val, row, col, *pLen;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  wenergy2_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar(1, "d", &m1, &n1, &l1);
  GetRhsVar(2, "i", &m2, &n2, &l2);
  size = 0;
  for (row = 1; row < (m2 - 1); row++)
    {
      size += (istk (l2)[row]) * (istk (l2)[row + m2]);
    }
  size = size * 3 + (istk (l2)[0]) * (istk (l2)[m2]);
  if (m1 * n1 != size)
    {
      sciprint ("Inputs are not size array and coefs!\n");
      return 0;
    }
  val = 0;
  if ((istk(l2)[0] != istk(l2)[1]) ||
      (istk(l2)[m2] != istk(l2)[m2+1]))
    val += 1;
  for (row = 1; row < (m2 - 1); row++)
    {
      if (istk (l2)[row] >= istk (l2)[row + 1])
	val += 1;
      if (istk (l2)[row + m2] >= istk (l2)[row + m2 + 1])
	val += 1;
    }
  if (val != 0)
    {
      sciprint ("Inputs are not size array!\n");
      return 0;
    }
  pLen = malloc (m2 * n2 * sizeof (int));
  for (row = 0; row < n2; row++)
    {
      for (col = 0; col < m2; col++)
	pLen[row + col * n2] = istk (l2)[col + row * m2];
    }

  switch (flow) {
  case 1:
    {
      m3 = 1;
      m4 = 1;
      m5 = 1;
      m6 = 1;
      n3 = 1;
      n4 = m2 - 2;
      n5 = m2 - 2;
      n6 = m2 - 2;
      CreateVar(3, "d", &m3, &n3, &l3);
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      wenergy_4output (stk(l1), m1*n1, pLen, stk(l3), stk(l4), 
		       stk(l5), stk(l6), n4, m2 - 2);
      LhsVar(1) = 3;
      LhsVar(2) = 4;
      LhsVar(3) = 5;
      LhsVar(4) = 6;
      break;
    }
  case 2:
    {
      m3 = 1;
      n3 = 1;
      m4 = 1;
      n4 = m2 - 2;
      CreateVar(3, "d", &m3, &n3, &l3);
      CreateVar(4, "d", &m4, &n4, &l4);
      wenergy_2output (stk(l1), m1*n1, pLen, stk(l3), stk(l4), 
		      n4, m2-2);
      LhsVar(1) = 3;
      LhsVar(2) = 4;
      break;
    }
  default:
    break;
  }
  free(pLen);
  return 0;
}

int
int_detcoef2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7;
  static int minrhs=4, maxrhs=4, minlhs=1, maxlhs=3;
  int errCode, flow, row, size, val, *pLen, col;
  static char h[] = { "h" };
  static char v[] = { "v" };
  static char d[] = { "d" };
  
  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);
  
  detcoef2_form_validate(&errCode, &flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar(1, "c", &m1, &n1, &l1);
  GetRhsVar(2, "d", &m2, &n2, &l2);
  GetRhsVar(3, "i", &m3, &n3, &l3);
  GetRhsVar(4, "i", &m4, &n4, &l4);
  
  detcoef2_content_validate(&errCode,flow,l1,l2,l3,l4);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  size = 0;
  for (row = 1; row < (m3 - 1); row++)
    {
      size += (istk (l3)[row]) * (istk (l3)[row + m3]);
    }
  size = size * 3 + (istk (l3)[0]) * (istk (l3)[m3]);
  if (m2 * n2 != size)
    {
      sciprint ("Inputs are not size array and coefs!\n");
      return 0;
    }
  val = 0;
  if ((istk(l3)[0] != istk(l3)[1]) ||
      (istk(l3)[m3] != istk(l3)[m3+1]))
    val += 1;
  for (row = 1; row < (m3 - 1); row++)
    {
      if (istk (l3)[row] >= istk (l3)[row + 1])
	val += 1;
      if (istk (l3)[row + m3] >= istk (l3)[row + m3 + 1])
	val += 1;
    }
  if (val != 0)
    {
      sciprint ("Inputs are not size array!\n");
      return 0;
    }
  if ((istk(l4)[0]<=0) || (istk(l4)[0]>(m3-2)))
    {
      sciprint ("Level Parameter is not valid for input matrix!\n");
      return 0;
    }


  pLen = malloc (m3 * n3 * sizeof (int));
  for (row = 0; row < n3; row++)
    {
      for (col = 0; col < m3; col++)
	pLen[row + col * n3] = istk (l3)[col + row * m3];
    }

  if ((Lhs==1) && (!strcmp(cstk(l1),"h")))
    {
      m5 = pLen[(m3 - 1 - istk (l4)[0]) * n3];
      n5 = pLen[(m3 - 1 - istk (l4)[0]) * n3 + 1];
      CreateVar(5, "d", &m5, &n5, &l5);
      detcoef2 (stk(l2), m2*n2 , stk(l5), m5*n5, pLen, m3-2, 
		istk(l4)[0], cstk(l1));
      LhsVar(1) = 5;
    }
  else if ((Lhs==1) && (!strcmp(cstk(l1),"v")))
    {
      m5 = pLen[(m3 - 1 - istk (l4)[0]) * n3];
      n5 = pLen[(m3 - 1 - istk (l4)[0]) * n3 + 1];
      CreateVar(5, "d", &m5, &n5, &l5);
      detcoef2 (stk(l2), m2*n2 , stk(l5), m5*n5, pLen, m3-2, 
		istk(l4)[0], cstk(l1));
      LhsVar(1) = 5;
    }
  else if ((Lhs==1) && (!strcmp(cstk(l1),"d")))
    {
      m5 = pLen[(m3 - 1 - istk (l4)[0]) * n3];
      n5 = pLen[(m3 - 1 - istk (l4)[0]) * n3 + 1];
      CreateVar(5, "d", &m5, &n5, &l5);
      detcoef2 (stk(l2), m2*n2 , stk(l5), m5*n5, pLen, m3-2, 
		istk(l4)[0], cstk(l1));
      LhsVar(1) = 5;
    }
  else if ((Lhs==1) && ((!strcmp(cstk(l1),"c")) ||
			(!strcmp(cstk(l1),"compact"))))
    {
      m5 = pLen[(m3 - 1 - istk (l4)[0]) * n3];
      n5 = pLen[(m3 - 1 - istk (l4)[0]) * n3 + 1]*3;
      CreateVar(5, "d", &m5, &n5, &l5);
      detcoef2 (stk(l2), m2*n2 , stk(l5), m5*n5/3, pLen, m3-2, 
		istk(l4)[0], h);
      detcoef2 (stk(l2), m2*n2 , stk(l5)+m5*n5/3, m5*n5/3, pLen, 
		m3-2, istk(l4)[0], v);
      detcoef2 (stk(l2), m2*n2 , stk(l5)+m5*n5*2/3, m5*n5/3, pLen, m3-2, 
		istk(l4)[0], d);
      LhsVar(1) = 5;
    }

  else if ((Lhs==3) && ((!strcmp(cstk(l1),"a")) ||
			(!strcmp(cstk(l1),"all"))))
    {
      m5 = pLen[(m3 - 1 - istk (l4)[0]) * n3];
      n5 = pLen[(m3 - 1 - istk (l4)[0]) * n3 + 1];
      m6 = m5;
      m7 = m5;
      n6 = n5;
      n7 = n5;
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      CreateVar(7, "d", &m7, &n7, &l7);
      detcoef2 (stk(l2), m2*n2 , stk(l5), m5*n5, pLen, m3-2, 
		istk(l4)[0], h);
      detcoef2 (stk(l2), m2*n2 , stk(l6), m6*n6, pLen, m3-2, 
		istk(l4)[0], v);
      detcoef2 (stk(l2), m2*n2 , stk(l7), m7*n7, pLen, m3-2, 
		istk(l4)[0], d);
      LhsVar(1) = 5;
      LhsVar(2) = 6;
      LhsVar(3) = 7;
    }
  else 
    {
      sciprint("Unknown Input Error!\n");
      return 0;
    }
  free(pLen);
  return 0;
}


int
int_appcoef2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6;
  static int minrhs=3, maxrhs=5, minlhs=1, maxlhs=1;
  int errCode, flow, size, row, val, *pLen, family, member, ii, col;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  appcoef2_form_validate(&errCode, &flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;
  l5 = 0;

  //sciprint("after form check!\n");
  GetRhsVar(1, "d", &m1, &n1, &l1);
  GetRhsVar(2, "i", &m2, &n2, &l2);
  size = 0;
  for (row = 1; row < (m2 - 1); row++)
    {
      size += (istk (l2)[row]) * (istk (l2)[row + m2]);
    }
  size = size * 3 + (istk (l2)[0]) * (istk (l2)[m2]);
  if (m1 * n1 != size)
    {
      sciprint ("Inputs are not size array and coefs!\n");
      return 0;
    }
  val = 0;
  if ((istk(l2)[0] != istk(l2)[1]) ||
      (istk(l2)[m2] != istk(l2)[m2+1]))
    val += 1;
  for (row = 1; row < (m2 - 1); row++)
    {
      if (istk (l2)[row] >= istk (l2)[row + 1])
	val += 1;
      if (istk (l2)[row + m2] >= istk (l2)[row + m2 + 1])
	val += 1;
    }
  if (val != 0)
    {
      sciprint ("Inputs are not size array!\n");
      return 0;
    }

  pLen = malloc (m2 * n2 * sizeof (int));
  for (row = 0; row < n2; row++)
    {
      for (col = 0; col < m2; col++)
	pLen[row + col * n2] = istk (l2)[col + row * m2];
    }

  switch (flow) {
  case 1:
    {
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      appcoef2_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      //sciprint("after content check!\n");
      if (istk(l4)[0]>m2-2)
	{
	  sciprint("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if ((istk(l2)[0] < pWaveStruct.length) || 
	  (istk(l2)[m2] < pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = pLen[(m2 - 1 - istk (l4)[0]) * n2];
      n5 = pLen[(m2 - 1 - istk (l4)[0]) * n2 + 1];
      CreateVar (5, "d", &m5, &n5, &l5);
      appcoef2 (stk(l1), m1*n1, pWaveStruct.pLowPass, 
		pWaveStruct.pHiPass, pWaveStruct.length, 
		stk(l5), m5, n5, pLen, m2-2, istk(l4)[0], dwtMode);
      LhsVar(1) = 5;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(3, "c", &m3, &n3, &l3);
      appcoef2_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if ((istk(l2)[0] < pWaveStruct.length) || 
	  (istk(l2)[m2] < pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m4 = pLen[0];
      n4 = pLen[1];
      CreateVar(4, "d", &m4, &n4, &l4);
      appcoef2 (stk(l1), m1*n1, pWaveStruct.pLowPass, 
		pWaveStruct.pHiPass, pWaveStruct.length, 
		stk(l4), m4, n4, pLen, m2-2, m2-2, dwtMode);
      LhsVar(1) = 4;
      filter_clear();
      break;
    }
  case 3:
    {
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      appcoef2_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l2)[0] < m3*n3) || 
	  (istk(l2)[m2] < m3*n3))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = pLen[0];
      n5 = pLen[1];
      CreateVar(5, "d", &m5, &n5, &l5);
      appcoef2 (stk(l1), m1*n1, stk(l3), stk(l4), m3*n3, 
		stk(l5), m5, n5, pLen, m2-2, m2-2, dwtMode);
      LhsVar(1) = 5;
      break;
    }
  case 4:
    {
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      appcoef2_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l2)[0] < m3*n3) || 
	  (istk(l2)[m2] < m3*n3))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      if (istk(l5)[0]>m2-2)
	{
	  sciprint("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      m6 = pLen[(m2 - 1 - istk (l5)[0]) * n2];
      n6 = pLen[(m2 - 1 - istk (l5)[0]) * n2 + 1];
      CreateVar (6, "d", &m6, &n6, &l6);
      appcoef2 (stk(l1), m1*n1, stk(l3), stk(l4), m3*n3, 
		stk(l6), m6, n6, pLen, m2-2, istk(l5)[0], dwtMode);
      LhsVar(1) = 6;
      break;
    }
  default:
    break;
  }
  free(pLen);
  return 0;
}


int
int_wrcoef2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7;
  static int minrhs=4, maxrhs=6, minlhs=1, maxlhs=1;
  int errCode, flow, size, row, val, *pLen, family, member, ii, col;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  wrcoef2_form_validate (&errCode,&flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
  
  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;
  l5 = 0;
  l6 = 0;
  
  GetRhsVar(1, "c", &m1, &n1, &l1);
  GetRhsVar(2, "d", &m2, &n2, &l2);
  GetRhsVar(3, "i", &m3, &n3, &l3);

  size = 0;
  for (row = 1; row < (m3 - 1); row++)
    {
      size += (istk (l3)[row]) * (istk (l3)[row + m3]);
    }
  size = size * 3 + (istk (l3)[0]) * (istk (l3)[m3]);
  if (m2 * n2 != size)
    {
      sciprint ("Inputs are not size array and coefs!\n");
      return 0;
    }
  val = 0;
  if ((istk(l3)[0] != istk(l3)[1]) ||
      (istk(l3)[m3] != istk(l3)[m3+1]))
    val += 1;
  for (row = 1; row < (m3 - 1); row++)
    {
      if (istk (l3)[row] >= istk (l3)[row + 1])
	val += 1;
      if (istk (l3)[row + m3] >= istk (l3)[row + m3 + 1])
	val += 1;
    }
  if (val != 0)
    {
      sciprint ("Inputs are not size array!\n");
      return 0;
    }

  pLen = malloc (m3 * n3 * sizeof (int));
  for (row = 0; row < n3; row++)
    {
      for (col = 0; col < m3; col++)
	pLen[row + col * n3] = istk (l3)[col + row * m3];
    }

  switch (flow) {
  case 1:
    {
      GetRhsVar(4, "c", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      wrcoef2_content_validate (&errCode,flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l4),&family,&member);
      wavelet_fun_parser (cstk(l4), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if ((istk(l3)[0] < pWaveStruct.length) || 
	  (istk(l3)[m3] < pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      if ((istk(l5)[0]<=0) || (istk(l5)[0]>(m3-2)))
	{
	  sciprint ("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      m6 = pLen[(m3 - 1) * n3];
      n6 = pLen[(m3 - 1) * n3 + 1];
      CreateVar (6, "d", &m6, &n6, &l6);
      wrcoef2 (stk(l2), m2*n2, pWaveStruct.pLowPass, 
	       pWaveStruct.pHiPass, pWaveStruct.length, 
	       stk(l6), m6, n6, pLen, m3-2, istk(l5)[0], 
	       cstk(l1), dwtMode);
      LhsVar(1) = 6;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "d", &m5, &n5, &l5);
      GetRhsVar(6, "i", &m6, &n6, &l6);
      wrcoef2_content_validate (&errCode,flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l6)[0]<=0) || (istk(l6)[0]>(m3-2)))
	{
	  sciprint ("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      if ((istk(l3)[0] < m4*n4) || 
	  (istk(l3)[m3] < m4*n4))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m7 = pLen[(m3 - 1) * n3];
      n7 = pLen[(m3 - 1) * n3 + 1];
      CreateVar (7, "d", &m7, &n7, &l7);
      wrcoef2 (stk(l2), m2*n2, stk(l4), stk(l5), 
	       m4*n4, stk(l7), m7, n7, pLen, m3-2, istk(l6)[0], 
	       cstk(l1), dwtMode);
      LhsVar(1) = 7;
      break;
    }
  case 3:
    {
      GetRhsVar(4, "c", &m4, &n4, &l4);
      wrcoef2_content_validate (&errCode,flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l4),&family,&member);
      wavelet_fun_parser (cstk(l4), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if ((istk(l3)[0] < pWaveStruct.length) || 
	  (istk(l3)[m3] < pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = pLen[(m3 - 1) * n3];
      n5 = pLen[(m3 - 1) * n3 + 1];
      CreateVar (5, "d", &m5, &n5, &l5);
      wrcoef2 (stk(l2), m2*n2, pWaveStruct.pLowPass, 
	       pWaveStruct.pHiPass, pWaveStruct.length, 
	       stk(l5), m5, n5, pLen, m3-2, m3-2, 
	       cstk(l1), dwtMode);
      LhsVar(1) = 5;
      filter_clear();
      break;
    }
  case 4:
    {
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "d", &m5, &n5, &l5);
      wrcoef2_content_validate (&errCode,flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l3)[0] < m4*n4) || 
	  (istk(l3)[m3] < m4*n4))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m6 = pLen[(m3 - 1) * n3];
      n6 = pLen[(m3 - 1) * n3 + 1];
      CreateVar (6, "d", &m6, &n6, &l6);
      wrcoef2 (stk(l2), m2*n2, stk(l4), stk(l5), 
	       m4*n4, stk(l6), m6, n6, pLen, m3-2, m3-2, 
	       cstk(l1), dwtMode);
      LhsVar(1) = 6;
      break;
    }
  default:
    break;
  }
  free(pLen);
  return 0;
}


int
int_upwlev2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7;
  static int minrhs=3, maxrhs=4, minlhs=3, maxlhs=3;
  int errCode, flow, size, row, val, *pLen, family, member, ii, col;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  upwlev2_form_validate (&errCode, &flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;

  GetRhsVar(1, "d", &m1, &n1, &l1);
  GetRhsVar(2, "i", &m2, &n2, &l2);

  size = 0;
  for (row = 1; row < (m2 - 1); row++)
    {
      size += (istk (l2)[row]) * (istk (l2)[row + m2]);
    }
  size = size * 3 + (istk (l2)[0]) * (istk (l2)[m2]);
  if (m1 * n1 != size)
    {
      sciprint ("Inputs are not size array and coefs!\n");
      return 0;
    }
  val = 0;
  if ((istk(l2)[0] != istk(l2)[1]) ||
      (istk(l2)[m2] != istk(l2)[m2+1]))
    val += 1;
  for (row = 1; row < (m2 - 1); row++)
    {
      if (istk (l2)[row] >= istk (l2)[row + 1])
	val += 1;
      if (istk (l2)[row + m2] >= istk (l2)[row + m2 + 1])
	val += 1;
    }
  if (val != 0)
    {
      sciprint ("Inputs are not size array!\n");
      return 0;
    }

  pLen = malloc (m2 * n2 * sizeof (int));
  for (row = 0; row < n2; row++)
    {
      for (col = 0; col < m2; col++)
	pLen[row + col * n2] = istk (l2)[col + row * m2];
    }

  if (m2<=3)
    {
      sciprint("Inputs are not coef and length!\n");
      return 0;
    }

  if (m2<4)
  {
	  sciprint("Decomposition level less than 2 is not accepted!\n");
      return 0;
  }

  switch (flow) {
  case 1:
    {
      GetRhsVar(3, "c", &m3, &n3, &l3);
      upwlev2_content_validate (&errCode, flow, l1, l2, l3, l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if ((istk(l2)[0] < pWaveStruct.length) || 
	  (istk(l2)[m2] < pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      
      m4 = 1;
      m5 = m2 - 1;
      n4 = n1*m1 - 4*pLen[0]*pLen[1] + pLen[4]*pLen[5];
      n5 = 2;
      m6 = pLen[0];
      n6 = pLen[1];
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "i", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      upwlev2 (stk(l1), m1*n1, pWaveStruct.pLowPass, 
	       pWaveStruct.pHiPass,pWaveStruct.length,
	       pLen, m2, n2, stk(l6), m6*n6, stk(l4), m4*n4,
	       istk(l5), m5, n5, m2-2, dwtMode);
      LhsVar(1) = 4;
      LhsVar(2) = 5;
      LhsVar(3) = 6;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      upwlev2_content_validate (&errCode, flow, l1, l2, l3, l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l2)[0] < m3*n3) || 
	  (istk(l2)[m2] < m3*n3))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      m6 = m2 - 1;
      n5 = n1*m1 - 4*pLen[0]*pLen[1] + pLen[4]*pLen[5];
      n6 = 2;
      m7 = pLen[0];
      n7 = pLen[1];
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "i", &m6, &n6, &l6);
      CreateVar(7, "d", &m7, &n7, &l7);
      upwlev2 (stk(l1), m1*n1, stk(l3), stk(l4), m3*n3,
	       pLen, m2, n2, stk(l7), m7*n7, stk(l5), m5*n5,
	       istk(l6), m6, n6, m2-2, dwtMode);
      LhsVar(1) = 5;
      LhsVar(2) = 6;
      LhsVar(3) = 7;
      break;
    }
  default:
    break;
  }
  free(pLen);
  return 0;
}


int
int_upcoef2 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7;
  static int minrhs=3, maxrhs=6, minlhs=1, maxlhs=1;
  int errCode, flow, s2, s3, s4, s1, family, member, ii;
  Func syn_fun;
  swt_wavelet pWaveStruct;
  
  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);
  
  upcoef2_form_validate (&errCode, &flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;
  l5 = 0;
  l6 = 0;

  GetRhsVar(1, "c", &m1, &n1, &l1);
  GetRhsVar(2, "d", &m2, &n2, &l2);

  switch (flow) {
  case 1:
    {
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      upcoef2_content_validate (&errCode, flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if (istk(l4)[0]<1)
	{
	  sciprint("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      if ((m2<pWaveStruct.length) || (n2<pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      upcoef_len_cal (m2, pWaveStruct.length, istk(l4)[0], 
		      &s1, &s2);
      upcoef_len_cal (n2, pWaveStruct.length, istk(l4)[0], 
		      &s3, &s4);
      if ((istk(l5)[0]>s1) || (istk(l5)[1]>s3))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}
      m6 = istk(l5)[0];
      n6 = istk(l5)[1];
      CreateVar(6, "d", &m6, &n6, &l6);
      upcoef2 (stk(l2), m2, n2,	pWaveStruct.pLowPass, 
	       pWaveStruct.pHiPass, pWaveStruct.length,
	       stk(l6), m6, n6, s1, s3, istk(l4)[0], cstk(l1), 
	       dwtMode);
      LhsVar(1) = 6;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      GetRhsVar(6, "i", &m6, &n6, &l6);
      upcoef2_content_validate (&errCode, flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if (istk(l5)[0]<1)
	{
	  sciprint("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      if ((m2<m3*n3) || (n2<m3*n3))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      upcoef_len_cal (m2, m3*n3, istk(l5)[0], &s1, &s2);
      upcoef_len_cal (n2, m3*n3, istk(l5)[0], &s3, &s4);
      if ((istk(l6)[0]>s1) || (istk(l6)[1]>s3))
	{
	  sciprint("Size Parameter is not valid for input matrix!\n");
	  return 0;
	}
      m7 = istk(l6)[0];
      n7 = istk(l6)[1];
      CreateVar(7, "d", &m7, &n7, &l7);
      upcoef2 (stk(l2), m2, n2,	stk(l3), stk(l4), m3*n3,
	       stk(l7), m7, n7, s1, s3, istk(l5)[0], cstk(l1), 
	       dwtMode);
      LhsVar(1) = 7;
      break;
    }
  case 3:
    {
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      upcoef2_content_validate (&errCode, flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if (istk(l4)[0]<1)
	{
	  sciprint("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      if ((m2<pWaveStruct.length) || (n2<pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
	  //sciprint("before upcoef_len_cal!\n");
      upcoef_len_cal (m2, pWaveStruct.length, istk(l4)[0], 
		      &s1, &s2);
      upcoef_len_cal (n2, pWaveStruct.length, istk(l4)[0], 
		      &s3, &s4);
      m5 = s1;
      n5 = s3;
	  //sciprint("after upcoef_len_cal!\n");
      CreateVar(5, "d", &m5, &n5, &l5);
      upcoef2 (stk(l2), m2, n2,	pWaveStruct.pLowPass, 
	       pWaveStruct.pHiPass, pWaveStruct.length,
	       stk(l5), m5, n5, s1, s3, istk(l4)[0], cstk(l1), 
	       dwtMode);
      LhsVar(1) = 5;
      filter_clear();
      break;
    }
  case 4:
    {
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      upcoef2_content_validate (&errCode, flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if (istk(l5)[0]<1)
	{
	  sciprint("Level Parameter is not valid for input matrix!\n");
	  return 0;
	}
      if ((m2<m3*n3) || (n2<m3*n3))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      upcoef_len_cal (m2, m3*n3, istk(l5)[0], &s1, &s2);
      upcoef_len_cal (n2, m3*n3, istk(l5)[0], &s3, &s4);
      m6 = s1;
      n6 = s3;
      CreateVar(6, "d", &m6, &n6, &l6);
      upcoef2 (stk(l2), m2, n2,	stk(l3), stk(l4), m3*n3, 
	       stk(l6), m6, n6, s1, s3, istk(l5)[0], cstk(l1), 
	       dwtMode);
      LhsVar(1) = 6;
      break;
    }
  case 5:
    {
      GetRhsVar(3, "c", &m3, &n3, &l3);
      upcoef2_content_validate (&errCode, flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if ((m2<pWaveStruct.length) || (n2<pWaveStruct.length))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      upcoef_len_cal (m2, pWaveStruct.length, 1, &s1, &s2);
      upcoef_len_cal (n2, pWaveStruct.length, 1, &s3, &s4);
      m4 = s1;
      n4 = s3;
      CreateVar(4, "d", &m4, &n4, &l4);
      upcoef2 (stk(l2), m2, n2,	pWaveStruct.pLowPass, 
	       pWaveStruct.pHiPass, pWaveStruct.length, 
	       stk(l4), m4, n4, s1, s3, 1, cstk(l1), 
	       dwtMode);
      LhsVar(1) = 4;
      filter_clear();
      break;
    }
  case 6:
    {
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      upcoef2_content_validate (&errCode, flow,l1,l2,l3,l4,l5,l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((m2<m3*n3) || (n2<m3*n3))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      upcoef_len_cal (m2, m3*n3, 1, &s1, &s2);
      upcoef_len_cal (n2, m3*n3, 1, &s3, &s4);
      m5 = s1;
      n5 = s3;
      CreateVar(5, "d", &m5, &n5, &l5);
      upcoef2 (stk(l2), m2, n2,	stk(l3), stk(l4), m3*n3, 
	       stk(l5), m5, n5, s1, s3, 1, cstk(l1), 
	       dwtMode);
      LhsVar(1) = 5;
      break;
    }
  default:
    break;
  }
  
  return 0;
}
