/*
 * -------------------------------------------------------------------------
 * dwt1d_int.c -- 1-D signal decomposition and reconstruction interface
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
// #include <stack-c.h>

int 
int_dwt(char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int l7, m7, n7, minlhs=2, maxlhs=2, minrhs=2, maxrhs=5;
  int errCode, flow, family, member, ii, stride, val;
  Func ana_fun;
  swt_wavelet pWaveStruct;
  extend_method extMethod;

  CheckRhs (minrhs,maxrhs);
  CheckLhs (minlhs,maxlhs);

  //GetRhsVar(1, "d", &m1, &n1, &l1);

  dwt_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
	  //sciprint("right here!\n");
      return 0;			
    }

  //CheckLhs (minlhs,maxlhs);

  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;
  l5 = 0;

  switch (flow){
  case 1:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "c", &m2, &n2, &l2);
      dwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l2),&family,&member);
      wavelet_fun_parser (cstk(l2), &ii);
      ana_fun = wi[ii].analysis;
      (*ana_fun)(member, &pWaveStruct);
      wave_len_validate (m1*n1, pWaveStruct.length, 
			 &stride, &val);
	 
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m3 = 1;
      m4 = 1;
      n3 = (int)floor((m1*n1 + pWaveStruct.length - 1)/2);
      if (dwtMode==PER)
	     n3 = (int)ceil(((double)(m1*n1))/2.0);
      n4 = n3;
      CreateVar(3, "d", &m3, &n3, &l3);
      CreateVar(4, "d", &m4, &n4, &l4);
	  dwt_neo (stk(l1), m1*n1, pWaveStruct.pLowPass, 
	   pWaveStruct.pHiPass, pWaveStruct.length, 
	   stk(l3), stk(l4), n3, dwtMode);
      filter_clear();
      LhsVar(1) = 3;
      LhsVar(2) = 4;
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      dwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate (m1*n1, m3*n3, &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m4 = 1;
      m5 = 1;
      n4 = (int)floor((m1*n1 + m2*n2 - 1)/2);
      if (dwtMode==PER)
         n4 = (int)ceil(((double)(m1*n1))/2.0);
      n5 = n4;
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "d", &m5, &n5, &l5);
      dwt_neo (stk(l1), m1*n1, stk(l2), stk(l3), m2*n2, 
	   stk(l4), stk(l5), n4, dwtMode);
      LhsVar(1) = 4;
      LhsVar(2) = 5;
      break;
    }
  case 3:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "c", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "c", &m4, &n4, &l4);
      dwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      //dwt_write(cstk(l4),&errCode);
      extend_method_parse (cstk(l4), &extMethod);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l2),&family,&member);
      wavelet_fun_parser (cstk(l2), &ii);
      ana_fun = wi[ii].analysis;
      (*ana_fun)(member, &pWaveStruct);
      wave_len_validate (m1*n1, pWaveStruct.length, 
			 &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      m6 = 1;
      n5 = (int)floor((m1*n1 + pWaveStruct.length - 1)/2);
      if (extMethod==PER)
        n5 = (int)ceil(((double)(m1*n1))/2.0);
      n6 = n5;
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      dwt_neo (stk(l1), m1*n1, pWaveStruct.pLowPass, 
	   pWaveStruct.pHiPass, pWaveStruct.length, 
	   stk(l5), stk(l6), n5, extMethod);
      LhsVar(1) = 5;
      LhsVar(2) = 6;
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
      dwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      //dwt_write(cstk(l4),&errCode);
      extend_method_parse (cstk(l5), &extMethod);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate (m1*n1, m3*n3, &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m6 = 1;
      m7 = 1;
      n6 = (int)floor((m1*n1 + m2*n2 - 1)/2);
      if (extMethod==PER)
         n6 = (int)ceil(((double)(m1*n1))/2.0);
      n7 = n6;
      CreateVar(6, "d", &m6, &n6, &l6);
      CreateVar(7, "d", &m7, &n7, &l7);
      dwt_neo (stk(l1), m1*n1, stk(l2), stk(l3), m2*n2, 
	   stk(l6), stk(l7), n6, extMethod);
      LhsVar(1) = 6;
      LhsVar(2) = 7;
      break;
    }
  default:
    break;
  }

  //sciprint("flow=%d\n", flow);
  return 0;
}

int
int_idwt (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7, l8, m8, n8;
  static minlhs=1, maxlhs=1, minrhs=3, maxrhs=7;
  int errCode, flow, family, member, ii, stride, val;
  Func syn_fun;
  swt_wavelet pWaveStruct;
  extend_method extMethod;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  idwt_form_validate(&errCode,&flow);
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
 
  switch (flow) {
  case 1:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 pWaveStruct.length, &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m4 = 1;
      if (m1*n1 != 0)// && (dwtMode!=PER))
	n4 = m1*n1*2 - pWaveStruct.length + 2;
      if ((m1*n1 != 0) && (dwtMode==PER))
         n4 = m1*n1*2;
      if (m1*n1 == 0)// && (dwtMode!=PER))
	n4 = m2*n2*2 - pWaveStruct.length + 2;
      if ((m1*n1 == 0) && (dwtMode==PER))
         n4 = m2*n2*2;
      CreateVar(4, "d", &m4, &n4, &l4);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, pWaveStruct.pHiPass, 
			pWaveStruct.length, stk(l4), m4*n4); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, pWaveStruct.pLowPass, 
		     pWaveStruct.length, stk(l4), m4*n4); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
			  pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
			  pWaveStruct.length, stk(l4), m4*n4);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 4;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 (m3*n3), &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      if (m1*n1 != 0)// && (dwtMode!=PER))
	n5 = m1*n1*2 - m3*n3 + 2;
      if ((m1*n1 != 0) && (dwtMode==PER))
        n5 = m1*n1*2;
      if (m1*n1 == 0) //&& (dwtMode!=PER))
	n5 = m2*n2*2 - m3*n3 + 2;
      if ((m1*n1 == 0) && (dwtMode==PER))
        n5 = m2*n2*2;
      CreateVar(5, "d", &m5, &n5, &l5);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, stk(l4), 
		     n4*m4, stk(l5), m5*n5); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, stk(l3), 
			n3*m3, stk(l5), m5*n5); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
			  stk(l3), stk(l4), 
			  m3*n3, stk(l5), m5*n5);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 5;
      break;
    }
  case 3:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 pWaveStruct.length, &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      n5 = istk(l4)[0];
      /*if (n5 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		pWaveStruct.length - 1))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}*/
	  if (n5 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		pWaveStruct.length))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}
      CreateVar(5, "d", &m5, &n5, &l5);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, pWaveStruct.pHiPass, 
			pWaveStruct.length, stk(l5), m5*n5); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, pWaveStruct.pLowPass, 
		     pWaveStruct.length, stk(l5), m5*n5); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
		       pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
		       pWaveStruct.length, stk(l5), m5*n5);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 5;
      filter_clear();
      break;
    }
  case 4:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 (m3*n3), &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m6 = 1;
      n6 = istk(l5)[0];
      /*if (n6 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		m3*n3 - 1))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}*/
	  if (n6 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		m3*n3))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}
      CreateVar(6, "d", &m6, &n6, &l6);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, stk(l4), 
		     n4*m4, stk(l6), m6*n6); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, stk(l3), 
			n3*m3, stk(l6), m6*n6); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
		       stk(l3), stk(l4), 
		       m3*n3, stk(l6), m6*n6);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 6;
      break;
    }
  case 5:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "c", &m4, &n4, &l4);
      GetRhsVar(5, "c", &m5, &n5, &l5);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 pWaveStruct.length, &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      extend_method_parse (cstk(l5), &extMethod);
      m6 = 1;
      if (m1*n1 != 0)// && (dwtMode!=PER))
	n6 = m1*n1*2 - pWaveStruct.length + 2;
      if ((m1*n1 != 0) && (extMethod==PER))
        n6 = m1*n1*2;
      if (m1*n1 == 0) //&& (dwtMode!=PER))
	n6 = m2*n2*2 - pWaveStruct.length + 2;
      if ((m1*n1 == 0) && (extMethod==PER))
        n6 = m2*n2*2;
      CreateVar(6, "d", &m6, &n6, &l6);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, pWaveStruct.pHiPass, 
		     pWaveStruct.length, stk(l6), m6*n6); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, pWaveStruct.pLowPass, 
		     pWaveStruct.length, stk(l6), m6*n6); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
		       pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
		       pWaveStruct.length, stk(l6), m6*n6);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 6;
      filter_clear();
      break;
    }
  case 6:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      GetRhsVar(5, "c", &m5, &n5, &l5);
      GetRhsVar(6, "c", &m6, &n6, &l6);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 pWaveStruct.length, &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      extend_method_parse (cstk(l6), &extMethod);
      m7 = 1;
      n7 = istk(l4)[0];
      /*if (n7 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		pWaveStruct.length - 1))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}*/
	  if (n7 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		pWaveStruct.length))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}
      CreateVar(7, "d", &m7, &n7, &l7);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, pWaveStruct.pHiPass, 
			pWaveStruct.length, stk(l7), m7*n7); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, pWaveStruct.pLowPass, 
		     pWaveStruct.length, stk(l7), m7*n7); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
			  pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
			  pWaveStruct.length, stk(l7), m7*n7);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 7;
      filter_clear();
      break;
    }
  case 7:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "c", &m5, &n5, &l5);
      GetRhsVar(6, "c", &m6, &n6, &l6);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 (m3*n3), &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      extend_method_parse (cstk(l6), &extMethod);
      m7 = 1;
      if (m1*n1 != 0)// && (dwtMode!=PER))
	n7 = m1*n1*2 - m3*n3 + 2;
      if ((m1*n1 != 0) && (extMethod==PER))
          n7 = m1*n1*2;
      if (m1*n1 == 0)// && (dwtMode!=PER))
	n7 = m2*n2*2 - m3*n3 + 2;
      if ((m1*n1 == 0) && (extMethod==PER))
         n7 = m2*n2*2;
      CreateVar(7, "d", &m7, &n7, &l7);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, stk(l4), 
		     n4*m4, stk(l7), m7*n7); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, stk(l3), 
		     n3*m3, stk(l7), m7*n7); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
		       stk(l3), stk(l4), 
		       m3*n3, stk(l7), m7*n7);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 7;
      break;
    }
  case 8:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      GetRhsVar(6, "c", &m6, &n6, &l6);
      GetRhsVar(7, "c", &m7, &n7, &l7);
      idwt_content_validate(&errCode,flow,l1,l2,l3,l4,l5,l6,l7);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate ((m1*n1>0)?(m1*n1):(m2*n2), 
			 (m3*n3), &stride, &val);
      if (!val)
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      extend_method_parse (cstk(l7), &extMethod);
      m8 = 1;
      n8 = istk(l5)[0];
      /*if (n8 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		m3*n3 - 1))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}*/
	  if (n8 > (((m1*n1>0)?(m1*n1):(m2*n2))*2 + 
		m3*n3))
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}
      CreateVar(8, "d", &m8, &n8, &l8);
      if ((m1*n1==0) && (m2*n2!=0))
	idwt_detail_neo (stk(l2), m2*n2, stk(l4), 
		     n4*m4, stk(l8), m8*n8); 
      else if ((m1*n1!=0) && (m2*n2==0))
	idwt_approx_neo (stk(l1), m1*n1, stk(l3), 
		     n3*m3, stk(l8), m8*n8); 
      else if ((m1*n1!=0) && (m2*n2!=0))
	idwt_neo (stk(l1), stk(l2), m1*n1, 
		       stk(l3), stk(l4), 
		       m3*n3, stk(l8), m8*n8);//, extMethod);
      else
	{
	  errCode=UNKNOWN_INPUT_ERR; 
	  validate_print (errCode);
	  return 0;
	}
      LhsVar(1) = 8;
      break;
    }
  default:
    break;
  }

  return 0;
}

int
int_wavedec (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int minlhs=2, maxlhs=2, minrhs=3, maxrhs=4;
  int errCode, flow, stride, val, ii, family, member;
  int count, calLen, temLen;
  Func ana_fun;
  swt_wavelet pWaveStruct;

  CheckLhs(minlhs,maxlhs);
  CheckRhs(minrhs,maxrhs);

  


  wavedec_form_validate(&errCode,&flow);
  if (errCode != SUCCESS)
    {
      //sciprint("flow error!\n");
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
      wavedec_content_validate(&errCode,flow,l1,l2,l3,l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      ana_fun = wi[ii].analysis;
      (*ana_fun)(member, &pWaveStruct);
      wave_len_validate (m1*n1, pWaveStruct.length, &stride, &val);
      if ((!val) || (stride<istk(l2)[0]))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m4 = 1;
      m5 = 1;
      n4 = 0;
      calLen = n1 * m1;
      for (count = 0; count < istk (l2)[0]; count++)
	{
	  calLen += pWaveStruct.length - 1;
	  temLen = calLen/2;
	  n4 += temLen;
	  calLen = temLen;
	}
      n4 += temLen;
	  if (dwtMode==PER)
	  {
	    n4 = 0;
        calLen = n1 * m1;
        for (count = 0; count < istk (l2)[0]; count++)
	     {
	       //calLen += m3*n3 - 1;
		   calLen = (int)ceil(((double)(calLen))/2.0);
	       temLen = calLen;
	       n4 += temLen;
	       //calLen = temLen;
	     }
        n4 += temLen;
	  }
      n5 = istk (l2)[0] + 2;
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "i", &m5, &n5, &l5);
      wave_dec_len_cal (pWaveStruct.length, m1*n1, 
			istk(l2)[0], istk(l5));
      wavedec (stk(l1), m1*n1, stk(l4), m4*n4, 
	       pWaveStruct.pLowPass, pWaveStruct.pHiPass,
	       pWaveStruct.length, istk(l5), n5, 
	       istk(l2)[0], dwtMode);
      LhsVar(1) = 4;
      LhsVar(2) = 5;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
	  //sciprint("after get data!\n");
      wavedec_content_validate(&errCode,flow,l1,l2,l3,l4);
      if (errCode != SUCCESS)
	{
      //sciprint("enter here\n");
	  validate_print (errCode);
	  return 0;			
	}
      wave_len_validate (m1*n1, m3*n3, &stride, &val);
      if ((!val) || (stride<istk(l2)[0]))
	{
      //sciprint("enter here!\n");
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      m6 = 1;
      n5 = 0;
      calLen = n1 * m1;
      for (count = 0; count < istk (l2)[0]; count++)
	  {
	     calLen += m3*n3 - 1;
	     temLen = calLen/2;
	     n5 += temLen;
	     calLen = temLen;
	  }
      n5 += temLen;
	  if (dwtMode==PER)
	  {
	    n5 = 0;
        calLen = n1 * m1;
        for (count = 0; count < istk (l2)[0]; count++)
	     {
	       //calLen += m3*n3 - 1;
		   calLen = (int)ceil(((double)(calLen))/2.0);
	       temLen = calLen;
	       n5 += temLen;
	       //calLen = temLen;
	     }
        n5 += temLen;
	  }
      n6 = istk (l2)[0] + 2;
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "i", &m6, &n6, &l6);
      wave_dec_len_cal (m3*n3, m1*n1, 
			istk(l2)[0], istk(l6));
      wavedec (stk(l1), m1*n1, stk(l5), m5*n5, stk(l3), stk(l4), 
	       m3*n3, istk(l6), n6, istk(l2)[0], dwtMode);
      LhsVar(1) = 5;
      LhsVar(2) = 6;
      break;
    }
  default:
    break;
  }

  return 0;
}

int
int_waverec (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minrhs=3, maxrhs=4, minlhs=1, maxlhs=1;
  int errCode, flow, len, count, val, family, member, ii;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckLhs(minlhs,maxlhs);
  CheckRhs(minrhs,maxrhs);

  waverec_form_validate(&errCode,&flow);
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
      waverec_content_validate(&errCode,flow,l1,l2,l3,l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if (istk (l2)[0] < pWaveStruct.length)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      //for(count=0;count<pWaveStruct.length;count++)
      //	printf("%f\n",pWaveStruct.pLowPass[count]);
      m4 = 1;
      n4 = istk(l2)[m2*n2-1];
      CreateVar(4, "d", &m4, &n4, &l4);
      waverec (stk(l1), m1*n1, stk(l4), m4*n4, 
	       pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
	       pWaveStruct.length, istk(l2), m2*n2, 
	       m2*n2-2, dwtMode);
      LhsVar(1) = 4;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      waverec_content_validate(&errCode,flow,l1,l2,l3,l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      if (istk (l2)[0] < m3*n3)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      n5 = istk(l2)[m2*n2-1];
      CreateVar(5, "d", &m5, &n5, &l5);
      waverec (stk(l1), m1*n1, stk(l5), m5*n5, 
	       stk(l3), stk(l4), m3*n3, istk(l2), m2*n2, 
	       m2*n2-2, dwtMode);
      LhsVar(1) = 5;
      break;
    }
  default:
    break;
  }
  return 0;
}

int
int_wrcoef (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int l7, m7, n7;
  static int minrhs=4, maxrhs=6, minlhs=1, maxlhs=1;
  int errCode, flow, len, count, val, family, member, ii;
  Func syn_fun;
  swt_wavelet pWaveStruct;
  
  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  wrcoef_form_validate (&errCode, &flow);
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

  switch (flow) {
  case 1:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "i", &m3, &n3, &l3);
      GetRhsVar(4, "c", &m4, &n4, &l4);
      //printf("enter flow 1\n");
      len = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	len += istk (l3)[count];
      if (len != m2 * n2)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	{
	  if (istk (l3)[count] > istk (l3)[count + 1])
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
      //printf("before content check\n");
      wrcoef_content_validate(&errCode, flow, l1, l2, l3, 
			      l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      //printf("before wavelet parser\n");
      wavelet_parser(cstk(l4),&family,&member);
      wavelet_fun_parser (cstk(l4), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      //printf("after parser\n");
      if (istk (l3)[0] < pWaveStruct.length)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      n5 = istk(l3)[m3*n3-1];
      CreateVar(5, "d", &m5, &n5, &l5);
      //printf("before wrcoef\n");
      wrcoef (stk(l2), m2*n2, pWaveStruct.pLowPass, 
	      pWaveStruct.pHiPass, pWaveStruct.length, 
	      istk(l3), m3*n3, stk(l5), m5*n5, cstk(l1), 
	      m3*n3-2, m3*n3-2, dwtMode);
      LhsVar(1) = 5; 
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "i", &m3, &n3, &l3);
      GetRhsVar(4, "c", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      len = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	len += istk (l3)[count];
      if (len != m2 * n2)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	{
	  if (istk (l3)[count] > istk (l3)[count + 1])
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
      wrcoef_content_validate(&errCode, flow, l1, l2, l3, 
			      l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l4),&family,&member);
      wavelet_fun_parser (cstk(l4), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if (istk (l3)[0] < pWaveStruct.length)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      if (istk(l5)[0] > (m3*n3-2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m6 = 1;
      n6 = istk(l3)[m3*n3-1];
      CreateVar(6, "d", &m6, &n6, &l6);
      wrcoef (stk(l2), m2*n2, pWaveStruct.pLowPass, 
	      pWaveStruct.pHiPass, pWaveStruct.length, 
	      istk(l3), m3*n3, stk(l6), m6*n6, cstk(l1), 
	      m3*n3-2, istk(l5)[0], dwtMode);
      LhsVar(1) = 6;
      filter_clear();
      break;
    }
  case 3:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "i", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "d", &m5, &n5, &l5);
      len = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	len += istk (l3)[count];
      if (len != m2 * n2)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	{
	  if (istk (l3)[count] > istk (l3)[count + 1])
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
      wrcoef_content_validate(&errCode, flow, l1, l2, l3, 
			      l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if (istk (l3)[0] < m4*n4)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m6 = 1;
      n6 = istk(l3)[m3*n3-1];
      CreateVar(6, "d", &m6, &n6, &l6);
      wrcoef (stk(l2), m2*n2, stk(l4), stk(l5), m4*n4, 
	      istk(l3), m3*n3, stk(l6), m6*n6, cstk(l1), 
	      m3*n3-2, m3*n3-2, dwtMode);
      LhsVar(1) = 6;
      break;
    }
  case 4:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "i", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "d", &m5, &n5, &l5);
      GetRhsVar(6, "i", &m6, &n6, &l6);
      len = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	len += istk (l3)[count];
      if (len != m2 * n2)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m3 * n3 - 1); count++)
	{
	  if (istk (l3)[count] > istk (l3)[count + 1])
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
      wrcoef_content_validate(&errCode, flow, l1, l2, l3, 
			      l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if (istk(l6)[0] > (m3*n3-2))
	{
	  sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      if (istk (l3)[0] < m4*n4)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m7 = 1;
      n7 = istk(l3)[m3*n3-1];
      CreateVar(7, "d", &m7, &n7, &l7);
      wrcoef (stk(l2), m2*n2, stk(l4), stk(l5), m4*n4, 
	      istk(l3), m3*n3, stk(l7), m7*n7, cstk(l1), 
	      m3*n3-2, istk(l6)[0], dwtMode);
      LhsVar(1) = 7;
      break;
    }
  default:
    break;
  }

  return 0;
}

int
int_appcoef (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int minrhs=3, maxrhs=5, minlhs=1, maxlhs=1;
  int errCode, flow, len, count, val, family, member, ii;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  appcoef_form_validate (&errCode, &flow);
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

  switch (flow){
  case 1:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      appcoef_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if (istk (l2)[0] < pWaveStruct.length)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m4 = 1;
      n4 = istk(l2)[0];
      CreateVar(4, "d", &m4, &n4, &l4);
      appcoef (stk(l1), m1*n1, stk(l4), m4*n4,  
	       pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
	       pWaveStruct.length, istk(l2), m2*n2,
	       m2*n2-2, m2*n2-2, dwtMode);
      LhsVar(1) = 4;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      appcoef_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l4)[0]>(m2*n2-2)) || (istk(l4)[0]<0))
	{
	  sciprint ("Level Parameter is not valid for input vector!\n");
	  return 0;
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if (istk (l2)[0] < pWaveStruct.length)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      n5 = istk(l2)[n2 * m2 - 2 - istk(l4)[0] + 1];
      CreateVar(5, "d", &m5, &n5, &l5);
      appcoef (stk(l1), m1*n1, stk(l5), m5*n5, 
	       pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
	       pWaveStruct.length, istk(l2), m2*n2,
	       m2*n2-2, istk(l4)[0], dwtMode);
      LhsVar(1) = 5;
      filter_clear();
      break;
    }
  case 3:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      appcoef_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if (istk (l2)[0] < m3*n3)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m5 = 1;
      n5 = istk(l2)[0];
      CreateVar(5, "d", &m5, &n5, &l5);
      appcoef (stk(l1), m1*n1, stk(l5), m5*n5, stk(l3), stk(l4), 
	       m3*n3, istk(l2), m2*n2, m2*n2-2, m2*n2-2, dwtMode);
      LhsVar(1) = 5;
      break;
    }
  case 4:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      appcoef_content_validate(&errCode, flow, l1, l2, l3, l4, l5);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if ((istk(l5)[0]>(m2*n2-2)) || (istk(l5)[0]<0))
	{
	  sciprint ("Level Parameter is not valid for input vector!\n");
	  return 0;
	}
      if (istk (l2)[0] < m3*n3)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
      m6 = 1;
      n6 = istk (l2)[n2 * m2 - 2 - istk (l5)[0] + 1];
      CreateVar(6, "d", &m6, &n6, &l6);
      appcoef (stk(l1), m1*n1, stk(l6), m6*n6, stk(l3), stk(l4), 
	       m3*n3, istk(l2), m2*n2, m2*n2-2, istk(l5)[0], dwtMode);
      LhsVar(1) = 6;
      break;
    }
  default:
    break;
  }

  return 0;
}

int
int_detcoef (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int minrhs=2, maxrhs=3, minlhs=1, maxlhs=1;
  int errCode, flow, len, val, count;
  
  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  detcoef_form_validate (&errCode, &flow);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  switch (flow) {
  case 1:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      m3 = 1;
      n3 = istk(l2)[1];
      CreateVar(3, "d", &m3, &n3, &l3);
      detcoef(stk(l1), m1*n1, istk(l2), m2*n2, stk(l3), m3*n3,
	      m2*n2-2, m2*n2-2);
      LhsVar(1) = 3;
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "i", &m3, &n3, &l3);
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      if ((istk(l3)[0]>(m2*n2-2)) || (istk(l3)[0]<1))
	{
		sciprint ("Level Parameter is not valid for input vector!\n");
	  return 0;
	}
      m4 = 1;
      n4 = istk (l2)[n2 * m2 - istk (l3)[0] - 1];
      CreateVar(4, "d", &m4, &n4, &l4);
      detcoef(stk(l1), m1*n1, istk(l2), m2*n2, stk(l4), m4*n4,
	      m2*n2-2, istk(l3)[0]);
      LhsVar(1) = 4;
      break;
    }
  default:
    break;
  }

  return 0;
}

int
int_wenergy (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int minrhs=2, maxrhs=2, minlhs=2, maxlhs=2;
  int errCode, len, val, count;
  
  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);
  
  wenergy_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
  
  GetRhsVar(1, "d", &m1, &n1, &l1);
  GetRhsVar(2, "i", &m2, &n2, &l2);

  len = 0;
  for (count = 0; count < (m2 * n2 - 1); count++)
    len += istk (l2)[count];
  if (len != m1 * n1)
    {
      sciprint ("Inputs are not coef and length array!\n");
      return 0;
    }
  val = 0;
  for (count = 0; count < (m2 * n2 - 1); count++)
    {
      if (istk (l2)[count] > istk (l2)[count + 1])
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
      
  m3 = 1;
  n3 = 1;
  m4 = 1;
  n4 = m2 * n2 - 2;
  CreateVar(3, "d", &m3, &n3, &l3);
  CreateVar(4, "d", &m4, &n4, &l4);
  wenergy(stk(l1), m1 * n1, istk(l2), m2 * n2, 
	  stk(l3), m3 * n3, stk(l4), m4 * n4);
  LhsVar(1) = 3;
  LhsVar(2) = 4;
  return 0;
}

int
int_upcoef (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int l7, m7, n7;
  static int minrhs=3, maxrhs=6, minlhs=1, maxlhs=1;
  int errCode, flow, family, member, ii;
  int s1, s2;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  upcoef_form_validate (&errCode, &flow);
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

  switch (flow) {
  case 1:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      upcoef_content_validate(&errCode, flow, l1, l2, 
			      l3, l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
     // if (pWaveStruct.length>m2*n2)
	//{
	  //sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  //return 0;
	//}
      upcoef_len_cal (m2*n2, pWaveStruct.length, istk(l4)[0], 
		      &s1, &s2);
      m5 = 1;
      n5 = s1;
      CreateVar(5, "d", &m5, &n5, &l5);
      upcoef (stk(l2), m2*n2, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, stk(l5), 
	      m5*n5, s1, cstk(l1), istk(l4)[0]);
      LhsVar(1) = 5;
      break;
    }
  case 2:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      GetRhsVar(4, "i", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      upcoef_content_validate(&errCode, flow, l1, l2, 
			      l3, l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
     // if (pWaveStruct.length>m2*n2)
	//{
	 // sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	 // return 0;
	//}
      upcoef_len_cal (m2*n2, pWaveStruct.length, istk(l4)[0], 
		      &s1, &s2);
      if (istk(l5)[0]>s1)
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}
      m6 = 1;
      n6 = istk(l5)[0];
      CreateVar(6, "d", &m6, &n6, &l6);
      upcoef (stk(l2), m2*n2, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, stk(l6), 
	      m6*n6, s1, cstk(l1), istk(l4)[0]);
      LhsVar(1) = 6;
      filter_clear();
      break;
    }
  case 3:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      upcoef_content_validate(&errCode, flow, l1, l2, 
			      l3, l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      //if (m3*n3>m2*n2)
	//{
	  //sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  //return 0;
	//}
      upcoef_len_cal (m2*n2, m3*n3, istk(l5)[0], &s1, &s2);
      m6 = 1;
      n6 = s1;
      CreateVar(6, "d", &m6, &n6, &l6);
      upcoef (stk(l2), m2*n2, stk(l3), stk(l4), m3*n3, stk(l6), 
	      m6*n6, s1, cstk(l1), istk(l5)[0]);
      LhsVar(1) = 6;
      break;
    }
  case 4:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      GetRhsVar(5, "i", &m5, &n5, &l5);
      GetRhsVar(6, "i", &m6, &n6, &l6);
      upcoef_content_validate(&errCode, flow, l1, l2, 
			      l3, l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
     // if (m3*n3>m2*n2)
	//{
	  //sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  //return 0;
	//}
      upcoef_len_cal (m2*n2, m3*n3, istk(l5)[0], &s1, &s2);
      if (istk(l6)[0]>s1)
	{
	  sciprint("Length Parameter is not valid for input vector!\n");
	  return 0;
	}
      m7 = 1;
      n7 = istk(l6)[0];
      CreateVar(7, "d", &m7, &n7, &l7);
      upcoef (stk(l2), m2*n2, stk(l3), stk(l4), m3*n3, stk(l7), 
	      m7*n7, s1, cstk(l1), istk(l5)[0]);
      LhsVar(1) = 7;
      break;
    }
  case 5:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "c", &m3, &n3, &l3);
      upcoef_content_validate(&errCode, flow, l1, l2, 
			      l3, l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      //if (pWaveStruct.length>m2*n2)
	//{
	 // sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	 // return 0;
	//}
      upcoef_len_cal (m2*n2, pWaveStruct.length, 1, &s1, &s2);
      m4 = 1;
      n4 = s1;
      CreateVar(4, "d", &m4, &n4, &l4);
      upcoef (stk(l2), m2*n2, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, stk(l4), 
	      m4*n4, s1, cstk(l1), 1);
      LhsVar(1) = 4;
      filter_clear();
      break;
    }
  case 6:
    {
      GetRhsVar(1, "c", &m1, &n1, &l1);
      GetRhsVar(2, "d", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      upcoef_content_validate(&errCode, flow, l1, l2, 
			      l3, l4, l5, l6);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
     // if (m3*n3>m2*n2)
	//{
	 // sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  //return 0;
	//}
      upcoef_len_cal (m2*n2, m3*n3, 1, &s1, &s2);
      m5 = 1;
      n5 = s1;
      CreateVar(5, "d", &m5, &n5, &l5);
      upcoef (stk(l2), m2*n2, stk(l3), stk(l4), m3*n3,
	      stk(l5), m5*n5, s1, cstk(l1), 1);
      LhsVar(1) = 5;
      break;
    }
  default:
    break;
  }
  return 0;
}

int
int_upwlev (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int l7, m7, n7;
  static int minrhs=3, maxrhs=4, minlhs=3, maxlhs=3;
  int errCode, flow, len, count, val, family, member, ii;
  Func syn_fun;
  swt_wavelet pWaveStruct;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  upwlev_form_validate (&errCode, &flow);
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
      upwlev_content_validate(&errCode, flow, l1, l2, l3, l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if (m2*n2 <= 3)
	{
	  sciprint("Inputs are not coef and length array!\n");
	  return 0;
	}
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      wavelet_parser(cstk(l3),&family,&member);
      wavelet_fun_parser (cstk(l3), &ii);
      syn_fun = wi[ii].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      if (istk (l2)[0] < pWaveStruct.length)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
	  if (m2*n2<4)
	  {
         sciprint("Decomposition level less than 2 is not accepted!\n"); 
         return 0;
	  }
      m4 = 1;
      m5 = 1;
      m6 = 1;
      n4 = m1*n1 - 2*istk(l2)[0] + istk(l2)[2];
      n5 = m2 * n2 - 1;
      n6 = istk(l2)[0];
      CreateVar(4, "d", &m4, &n4, &l4);
      CreateVar(5, "i", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
      //printf("before upwlev!\n");
      upwlev (stk(l1), m1*n1, istk(l2), m2*n2, pWaveStruct.pLowPass, 
	      pWaveStruct.pHiPass, pWaveStruct.length,
	      stk(l4), m4*n4, istk(l5), m5*n5, stk(l6), m6*n6, 
	      m2*n2-2, dwtMode);
      //printf("after upwlev!\n");
      LhsVar(1) = 4;
      LhsVar(2) = 5;
      LhsVar(3) = 6;
      filter_clear();
      break;
    }
  case 2:
    {
      GetRhsVar(1, "d", &m1, &n1, &l1);
      GetRhsVar(2, "i", &m2, &n2, &l2);
      GetRhsVar(3, "d", &m3, &n3, &l3);
      GetRhsVar(4, "d", &m4, &n4, &l4);
      upwlev_content_validate(&errCode, flow, l1, l2, l3, l4);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      if (m2*n2 <= 3)
	{
	  sciprint("Inputs are not coef and length array!\n");
	  return 0;
	}
      len = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	len += istk (l2)[count];
      if (len != m1 * n1)
	{
	  sciprint ("Inputs are not coef and length array!\n");
	  return 0;
	}
      val = 0;
      for (count = 0; count < (m2 * n2 - 1); count++)
	{
	  if (istk (l2)[count] > istk (l2)[count + 1])
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
      if (istk(l2)[0]<m3*n3)
	{
	  sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	  return 0;
	}
	  if (m2*n2<4)
	  {
         sciprint("Decomposition level less than 2 is not accepted!\n"); 
         return 0;
	  }
      m5 = 1;
      m6 = 1;
      m7 = 1;
      n5 = m1*n1 - 2*istk(l2)[0] + istk(l2)[2];
      n6 = m2 * n2 - 1;
      n7 = istk(l2)[0];
      CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "i", &m6, &n6, &l6);
      CreateVar(7, "d", &m7, &n7, &l7);
      upwlev (stk(l1), m1*n1, istk(l2), m2*n2, stk(l3), 
	      stk(l4), m3*n3,
	      stk(l5), m5*n5, istk(l6), m6*n6, stk(l7), m7*n7, 
	      m2*n2-2, dwtMode);
      LhsVar(1) = 5;
      LhsVar(2) = 6;
      LhsVar(3) = 7;
      break;
    }
  default:
    break;
  }
  return 0;
}

