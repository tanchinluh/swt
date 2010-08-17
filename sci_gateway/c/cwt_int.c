/*
 * -------------------------------------------------------------------------
 * cwt_int.c -- continuous wavelet transform interface
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
#include "cwt.h"
#include "stack-c.h"

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
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
  
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  sinus_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);

  sinus_content_validate(&errCode, l1, l2, l3);
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
  sinus(stk(l4), n4, stk(l5), n5, 1);

  LhsVar(1) = 5;
  LhsVar(2) = 4;
  return 0;
}


/*-------------------------------------------
 * Poisson Scale Filter Generation
 *-----------------------------------------*/

int
int_poisson (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
  
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  poisson_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);

  poisson_content_validate(&errCode, l1, l2, l3);
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
  poisson(stk(l4), n4, stk(l5), n5, 1);

  LhsVar(1) = 5;
  LhsVar(2) = 4;
  return 0;
}

/*-------------------------------------------
 * Mexican Hat Scale Filter Generation
 *-----------------------------------------*/

int
int_mexihat (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
  
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  mexihat_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);

  mexihat_content_validate(&errCode, l1, l2, l3);
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
  mexihat(stk(l4), n4, stk(l5), n5, 1);

  LhsVar(1) = 5;
  LhsVar(2) = 4;
  return 0;
}

/*-------------------------------------------
 * Morlet Scale Filter Generation
 *-----------------------------------------*/

int
int_morlet (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
  
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  morlet_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);

  morlet_content_validate(&errCode, l1, l2, l3);
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
  morlet(stk(l4), n4, stk(l5), n5, 1);

  LhsVar(1) = 5;
  LhsVar(2) = 4;
  return 0;
}

/*-------------------------------------------
 * DoG Scale Filter Generation
 *-----------------------------------------*/

int
int_DOGauss (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
  
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  DOGauss_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);

 DOGauss_content_validate(&errCode, l1, l2, l3);
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
  DOGauss(stk(l4), n4, stk(l5), n5, 1);

  LhsVar(1) = 5;
  LhsVar(2) = 4;
  return 0;
}

/*-------------------------------------------
 * Gauss Scale Filter Generation
 *-----------------------------------------*/

int
int_Gauswavf (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int minlhs = 2, maxlhs = 2, minrhs = 4, maxrhs = 4;
  
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  Gauss_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);
  GetRhsVar (4, "i", &m4, &n4, &l4);

 Gauss_content_validate(&errCode, l1, l2, l3, l4);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m5 = 1;
  n5 = istk(l3)[0];
  m6 = 1;
  n6 = n5;

  
  CreateVar (5, "d", &m5, &n5, &l5);
  CreateVar (6, "d", &m6, &n6, &l6);
  linspace(stk(l1)[0], stk(l2)[0], istk(l3)[0], stk(l5), n5);
  Gauss(stk(l5), n5, stk(l6), n6, istk(l4)[0], 1);

  LhsVar(1) = 6;
  LhsVar(2) = 5;
  return 0;
}


/*-------------------------------------------
 * Shannon Scale Filter Generation
 *-----------------------------------------*/

int
int_shanwavf (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int l7r, l7c, m7, n7;
  static int minlhs = 2, maxlhs = 2, minrhs = 5, maxrhs = 5;
  
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

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);
  GetRhsVar (4, "d", &m4, &n4, &l4);
  GetRhsVar (5, "d", &m5, &n5, &l5);

  shanwavf_content_validate(&errCode, l1, l2, l3, l4, l5);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m6 = 1;
  n6 = istk(l3)[0];
  m7 = 1;
  n7 = n6;
  it = 1;

  CreateVar (6, "d", &m6, &n6, &l6);
  CreateCVar (7, "d",&it, &m7, &n7, &l7r, &l7c);

  linspace(stk(l1)[0], stk(l2)[0], istk(l3)[0], stk(l6), n6);
  shanwavf(stk(l6), n6, stk(l5)[0], stk(l4)[0], stk(l7r), stk(l7c), n7,1);

  LhsVar(1) = 7;
  LhsVar(2) = 6;

  return 0;
}

/*-------------------------------------------
 * Complex Morlet Scale Filter Generation
 *-----------------------------------------*/

int
int_cmorlet (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int l7r, l7c, m7, n7;
  static int minlhs = 2, maxlhs = 2, minrhs = 5, maxrhs = 5;
  
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

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);
  GetRhsVar (4, "d", &m4, &n4, &l4);
  GetRhsVar (5, "d", &m5, &n5, &l5);

  cmorlet_content_validate(&errCode, l1, l2, l3, l4, l5);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m6 = 1;
  n6 = istk(l3)[0];
  m7 = 1;
  n7 = n6;
  it = 1;

  CreateVar (6, "d", &m6, &n6, &l6);
  CreateCVar (7, "d",&it, &m7, &n7, &l7r, &l7c);

  linspace(stk(l1)[0], stk(l2)[0], istk(l3)[0], stk(l6), n6);
  cmorlet(stk(l6), n6, stk(l5)[0], stk(l4)[0], stk(l7r), stk(l7c), n7, 1);

  LhsVar(1) = 7;
  LhsVar(2) = 6;
  return 0;
}

/*-------------------------------------------
 * Frequency B-Spline Scale Filter Generation
 *-----------------------------------------*/

int
int_fbspwavf (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int l7, m7, n7, l8r, l8c, m8, n8;
  static int minlhs = 2, maxlhs = 2, minrhs = 6, maxrhs = 6;
  
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

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);
  GetRhsVar (4, "i", &m4, &n4, &l4);
  GetRhsVar (5, "d", &m5, &n5, &l5);
  GetRhsVar (6, "d", &m6, &n6, &l6);

  fbspwavf_content_validate(&errCode, l1, l2, l3, l4, l5, l6);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m7 = 1;
  n7 = istk(l3)[0];
  m8 = 1;
  n8 = n7;
  it = 1;

  CreateVar (7, "d", &m7, &n7, &l7);
  CreateCVar (8, "d",&it, &m8, &n8, &l8r, &l8c);

  linspace(stk(l1)[0], stk(l2)[0], istk(l3)[0], stk(l7), n7);
  fbspwavf(stk(l7), n7, istk(l4)[0], stk(l6)[0], stk(l5)[0], stk(l8r), stk(l8c), n8, 1);

  LhsVar(1) = 8;
  LhsVar(2) = 7;
  return 0;
}


/*-------------------------------------------
 * Complex Cauchy Scale Filter Generation
 *-----------------------------------------*/

int
int_cauchy (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5r, m5, n5, l5i;
  //static int l7r, l7c, m7, n7;
  static int minlhs = 2, maxlhs = 2, minrhs = 3, maxrhs = 3;
  
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

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);
  //GetRhsVar (4, "d", &m4, &n4, &l4);
  //GetRhsVar (5, "d", &m5, &n5, &l5);

  cauchy_content_validate(&errCode, l1, l2, l3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m4 = 1;
  n4 = istk(l3)[0];
  m5 = 1;
  n5 = n4;
  it = 1;

  CreateVar (4, "d", &m4, &n4, &l4);
  CreateCVar (5, "d",&it, &m5, &n5, &l5r, &l5i);

  linspace(stk(l1)[0], stk(l2)[0], istk(l3)[0], stk(l4), n4);
  cauchy_neo(stk(l4), n4, stk(l5r), stk(l5i), n5, 1);

  LhsVar(1) = 5;
  LhsVar(2) = 4;
  return 0;
}

/*-------------------------------------------
 * Complex Gauss Scale Filter Generation
 *-----------------------------------------*/

int
int_cgauss (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5, l6r, l6i, m6, n6;
  
  static int minlhs = 2, maxlhs = 2, minrhs = 4, maxrhs = 4;
  
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

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "i", &m3, &n3, &l3);
  GetRhsVar (4, "i", &m4, &n4, &l4);
  

  Gauss_content_validate(&errCode, l1, l2, l3, l4);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }


  m5 = 1;
  n5 = istk(l3)[0];
  m6 = 1;
  n6 = n5;
  it = 1;

  CreateVar (5, "d", &m5, &n5, &l5);
  CreateCVar (6, "d",&it, &m6, &n6, &l6r, &l6i);

  linspace(stk(l1)[0], stk(l2)[0], istk(l3)[0], stk(l5), n5);
  cgauss(stk(l5), n5, istk(l4)[0], stk(l6r), stk(l6i), n6, 1);

  LhsVar(1) = 6;
  LhsVar(2) = 5;
  return 0;
}

/*-------------------------------------------
 * Meyer Filter Generation
 *-----------------------------------------*/

int
int_meyeraux (char *fname)
{
  static int l1, m1, n1, l2, m2, n2;
  static int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 1;
  
  GetRhsVar (1, "d", &m1, &n1, &l1);
  m2 = 1;
  n2 = 1;
  CreateVar (2, "d", &m2, &n2, &l2);

  meyeraux(stk(l1)[0],stk(l2));
  LhsVar(1) = 2;
  return 0;
}

/*-------------------------------------------
 * Scale and Wavelet Function
 *-----------------------------------------*/

int
int_wavefun (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7, l3r, l3i;
  static int minlhs = 2, maxlhs = 5, minrhs = 2, maxrhs = 2;
  int ind, family, member, s1, s2, leng, count, level, it, l, errCode;
  double one, *lowfltr, *hifltr, *buff;
  char a[2]="a",d[2]="d";
  Func syn_fun, ana_fun;
  swt_wavelet pWaveStruct;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  wavefun_form_validate(&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  GetRhsVar (1, "c", &m1, &n1, &l1);
  GetRhsVar (2, "i", &m2, &n2, &l2);
  
  wavefun_content_validate(&errCode, l1, l2);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  level = istk(l2)[0];
  if (level==0)
	  level = 8;
  wavelet_fun_parser (cstk(l1), &ind);
  one = 1.0;
  
  if ((ind!=-1) && (wi[ind].rOrB==ORTH))
  {
	  if (Lhs!=3)
	  {
          sciprint("Output arguments number should be 3!\n");
		  return 0;
	  }
      wavelet_parser(cstk(l1),&family,&member);
	  syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  upcoef_len_cal (1, pWaveStruct.length, level, 
	       &s1, &s2);
	  swt_exp2(level, &leng);
      m3 = 1;
      n3 = leng*(pWaveStruct.length-1)+1;
	  m4 = 1;
	  n4 = n3;
	  m5 = 1;
	  n5 = n3;
	  CreateVar(3, "d", &m3, &n3, &l3);
	  CreateVar(4, "d", &m4, &n4, &l4);
	  CreateVar(5, "d", &m5, &n5, &l5);
	  for(count=0;count<n3;count++)
	  {
		  stk(l3)[count] = 0;
		  stk(l4)[count] = 0;
	  }
	  l=1;
      upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, &(stk(l3)[l]), 
	      s1, s1, a, level);
	  upcoef (&one, 1, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, &(stk(l4)[l]), 
	      s1, s1, d, level);
	  if ((family==COIFLETS) || (family==SYMLETS) || (family==DMEY))
	  {

		  for(count=0;count<n3;count++)
	     {
	         stk(l4)[count] = -1*stk(l4)[count];
		 }
	  }
	  for(count=0;count<n3;count++)
	  {
		  stk(l3)[count] = stk(l3)[count]*pow(sqrt(2),level);
		  stk(l4)[count] = stk(l4)[count]*pow(sqrt(2),level);
	  }
	  linspace(0.0, (double)(pWaveStruct.length-1), n3, stk(l5), n3);
      LhsVar(1) = 3;
	  LhsVar(2) = 4;
	  LhsVar(3) = 5;
	  filter_clear();
      return 0;
  }

  if ((ind!=-1) && (wi[ind].rOrB==BIORTH))
  {
	  if (Lhs!=5)
	  {
          sciprint("Output arguments number should be 5!\n");
		  return 0;
	  }
      wavelet_parser(cstk(l1),&family,&member);
	  ana_fun = wi[ind].analysis;
      (*ana_fun)(member, &pWaveStruct);
	  upcoef_len_cal (1, pWaveStruct.length, level, 
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
	  CreateVar(3, "d", &m3, &n3, &l3);
	  CreateVar(4, "d", &m4, &n4, &l4);
	  CreateVar(5, "d", &m5, &n5, &l5);
      CreateVar(6, "d", &m6, &n6, &l6);
	  CreateVar(7, "d", &m7, &n7, &l7);
	  for(count=0;count<n3;count++)
	  {
		  stk(l3)[count] = 0;
		  stk(l4)[count] = 0;
		  stk(l5)[count] = 0;
		  stk(l6)[count] = 0;
	  }
	  lowfltr = malloc(pWaveStruct.length*sizeof(double));
	  hifltr = malloc(pWaveStruct.length*sizeof(double));
      wrev(pWaveStruct.pLowPass, pWaveStruct.length, lowfltr, pWaveStruct.length);
	  qmf_wrev(lowfltr,pWaveStruct.length,hifltr,pWaveStruct.length);
	  l=1;
      upcoef (&one, 1, lowfltr, hifltr, pWaveStruct.length, &(stk(l3)[l]), 
	      s1, s1, a, level);
      upcoef (&one, 1, lowfltr, hifltr, pWaveStruct.length, &(stk(l4)[l]), 
	      s1, s1, d, level);
	  free(lowfltr);
	  free(hifltr);
	  filter_clear();
      syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
      upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, &(stk(l5)[l]), 
	      s1, s1, a, level);
	  upcoef (&one, 1, pWaveStruct.pLowPass,
	      pWaveStruct.pHiPass, pWaveStruct.length, &(stk(l6)[l]), 
	      s1, s1, d, level);
	   for(count=0;count<n3;count++)
	     {
			  stk(l3)[count] = stk(l3)[count]*pow(sqrt(2),level);
		      stk(l4)[count] = -1*stk(l4)[count]*pow(sqrt(2),level);
			  stk(l5)[count] = stk(l5)[count]*pow(sqrt(2),level);
		      stk(l6)[count] = -1*stk(l6)[count]*pow(sqrt(2),level);
		 }
	  linspace(0.0, (double)(pWaveStruct.length-1), n3, stk(l7), n3);
	  LhsVar(1) = 3;
	  LhsVar(2) = 4;
	  LhsVar(3) = 5;
	  LhsVar(4) = 6;
	  LhsVar(5) = 7;
	  filter_clear();
	  return 0;
  }

  cwt_fun_parser(cstk(l1), &ind);
  if ((ind!=-1) && (ci[ind].phipsi==PSI_ONLY) && (ci[ind].realOrComplex==REAL) )
  {
	  if (Lhs!=2)
	  {
          sciprint("Output arguments number should be 2!\n");
		  return 0;
	  }
	  swt_exp2(level, &leng);
      m3 = 1;
      n3 = leng;
	  m4 = 1;
	  n4 = n3;
      CreateVar(3, "d", &m3, &n3, &l3);
	  CreateVar(4, "d", &m4, &n4, &l4);
	  linspace(ci[ind].lb, ci[ind].ub, n3, stk(l4), n3);
	  (*(ci[ind].scalef))(stk(l4),n3,stk(l3),n3,ci[ind].cpsi);
	  LhsVar(1) = 3;
	  LhsVar(2) = 4;
	  return 0;
  }

  if ((ind!=-1) && (ci[ind].phipsi==PSI_ONLY) && (ci[ind].realOrComplex==COMPLEX) )
  {
	  if (Lhs!=2)
	  {
          sciprint("Output arguments number should be 2!\n");
		  return 0;
	  }
	  swt_exp2(level, &leng);
      m3 = 1;
      n3 = leng;
	  m4 = 1;
	  n4 = n3;
	  it = 1;
      CreateCVar(3, "d", &it, &m3, &n3, &l3r, &l3i);
	  CreateVar(4, "d", &m4, &n4, &l4);
	  buff = malloc(n3*2*sizeof(double));
	  linspace(ci[ind].lb, ci[ind].ub, n3, stk(l4), n3);
	  (*(ci[ind].scalef))(stk(l4),n3,buff,n3,ci[ind].cpsi);
	  for(count=0;count<n3;count++)
	  {
        stk(l3r)[count]=buff[count];
		stk(l3i)[count]=buff[count+n3];
	  }
	  free(buff);
	  LhsVar(1) = 3;
	  LhsVar(2) = 4;
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
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int l5, m5, n5, l6, m6, n6, l7, m7, n7;
  static int minlhs = 5, maxlhs = 5, minrhs = 2, maxrhs = 2;
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

  GetRhsVar (1, "c", &m1, &n1, &l1);
  GetRhsVar (2, "i", &m2, &n2, &l2);
  
  wavefun2_content_validate(&errCode, l1, l2);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  level = istk(l2)[0];
  if (level==0)
	  level = 4;

  wavelet_fun_parser (cstk(l1), &ind);
  one = 1.0;

 if ((ind!=-1) && (wi[ind].rOrB==ORTH))
  {
      wavelet_parser(cstk(l1),&family,&member);
	  syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  upcoef_len_cal (1, pWaveStruct.length, level, 
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
	  CreateVar(3, "d", &m3, &n3, &l3);
	  CreateVar(4, "d", &m4, &n4, &l4);
	  CreateVar(5, "d", &m5, &n5, &l5);
	  CreateVar(6, "d", &m6, &n6, &l6);
	  CreateVar(7, "d", &m7, &n7, &l7);

	  phi=malloc(n3*sizeof(double));
	  psi=malloc(n3*sizeof(double));
      xval=malloc(n3*sizeof(double));
	  for(count=0;count<n3;count++)
	  {
		  phi[count] = 0;
		  psi[count] = 0;
		
	  }
	  l=(int)(floor((n3-s1)/2));
      upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, phi+l, 
	      s1, s1, a, level);
	  upcoef (&one, 1, pWaveStruct.pLowPass,
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
			  stk(l3)[row+col*n3]=phi[col]*phi[row]*pow(sqrt(2),level*2);
			  stk(l4)[row+col*n3]=psi[col]*phi[row]*pow(sqrt(2),level*2);
              stk(l5)[row+col*n3]=phi[col]*psi[row]*pow(sqrt(2),level*2);
			  stk(l6)[row+col*n3]=psi[col]*psi[row]*pow(sqrt(2),level*2);
			  stk(l7)[row+col*n3]=xval[col]*xval[row];
		  }
	  }
	  free(phi);
	  free(psi);
	  free(xval);
      LhsVar(1) = 3;
	  LhsVar(2) = 4;
	  LhsVar(3) = 5;
	  LhsVar(4) = 6;
	  LhsVar(5) = 7;
	  filter_clear();
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
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l4r, l4i, l1r, l1i;
  static int minlhs = 1, maxlhs = 1, minrhs = 3, maxrhs = 3;
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
  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  GetRhsVar (3, "c", &m3, &n3, &l3);

  cwt_content_validate(&errCode, l1, l2, l3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m4 = m2*n2;
  n4 = m1*n1;
  it = 1;
  
  wavelet_fun_parser (cstk(l3), &ind1);
  cwt_fun_parser(cstk(l3), &ind2);

  if (ind1!=-1) 
  {
	  
	  wavelet_parser(cstk(l3),&family,&member);
	  syn_fun = wi[ind1].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  plen=1024*(pWaveStruct.length-1)+1;
	  filter_clear();
	  psi=malloc(plen*sizeof(double));
	  full_range_scalef (cstk(l3), psi, plen);
	  //ocumsum(psi,plen);
	  CreateVar(4, "d", &m4, &n4, &l4);
	  temp = malloc(m4*n4*sizeof(double));
	  //tim = (int)(floor(len/n4));
	  //scale = (int)(floor(stk(l2)[0]));
	  for (count=0;count<m4;count++)
	  {
		  scale = (int)(floor(stk(l2)[count]));
		  if (scale<1)
			  scale = 1;
		  cwt_len_cal (plen, (scale*pWaveStruct.length+1), &len, &delta);
		  f = malloc(len*sizeof(double));
		  scale_real (psi, plen, delta, f, len);
		  
		  for(i=0;i<len;i++)
			f[i]=f[i]/sqrt(scale);
		  cwt_conv_real (stk(l1), n4, f, len, temp+count*n4, n4);
	  }
	  free(f);
	  free(psi);
	  matrix_tran(temp,m4,n4,stk(l4),n4,m4);
	  free(temp);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==REAL))
  {
	 
	  
	  CreateVar(4, "d", &m4, &n4, &l4);
	  temp = malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  //scale = (int)(floor(stk(l2)[count]));
		  plen=(int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count]))+1;
		  if (plen<3)
			  plen = 3;
	      psi=malloc(plen*sizeof(double));
	      full_range_scalef (cstk(l3), psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(stk(l2)[count]);
		  cwt_conv_real (stk(l1), n4, psi, plen, temp+count*n4, n4);
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,stk(l4),n4,m4);
	  free(temp);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==COMPLEX))
  {
	 
	  
	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	  temp = malloc(m4*n4*sizeof(double));
	  tempi = malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  plen=((int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count])+1))*2;//(int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count]+1))*2;
		  if (plen<3)
			  plen = 3;
	      psi=malloc(plen*sizeof(double));
	      full_range_scalef (cstk(l3), psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(stk(l2)[count]);
		  cwt_conv_complex (stk(l1), n4, psi, psi+plen/2, plen/2, temp+count*n4, tempi+count*n4, n4); 
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,stk(l4r),n4,m4);
	  matrix_tran(tempi,m4,n4,stk(l4i),n4,m4);
	  free(temp);
	  free(tempi);
  }
  LhsVar(1) = 4;

  return 0;
 }
 
 if (flow==2)
 {
	  GetRhsCVar (1, "d", &it, &m1, &n1, &l1r, &l1i);
      GetRhsVar (2, "d", &m2, &n2, &l2);
      GetRhsVar (3, "c", &m3, &n3, &l3);

	   cwt_content_validate(&errCode, l1, l2, l3);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  m4 = m2*n2;
  n4 = m1*n1;
  it = 1;
  
  wavelet_fun_parser (cstk(l3), &ind1);
  cwt_fun_parser(cstk(l3), &ind2);

  if (ind1!=-1) 
  {
	  
	  wavelet_parser(cstk(l3),&family,&member);
	  syn_fun = wi[ind1].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  plen=1024*(pWaveStruct.length-1)+1;
	  filter_clear();
	  psi=malloc(plen*sizeof(double));
	  full_range_scalef (cstk(l3), psi, plen);

	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	  temp = malloc(m4*n4*sizeof(double));
	  tempi = malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  scale = (int)(floor(stk(l2)[count]));
		  if (scale<1)
			  scale = 1;
		  cwt_len_cal (plen, (scale*pWaveStruct.length+1), &len, &delta);
		  f = malloc(len*sizeof(double));
		  scale_real (psi, plen, delta, f, len);
		  
		  for(i=0;i<len;i++)
			f[i]=f[i]/sqrt(scale);
		  cwt_conv_real (stk(l1r), n4, f, len, temp+count*n4, n4);
		  cwt_conv_real (stk(l1i), n4, f, len, tempi+count*n4, n4);
	  }
	  free(f);
	  free(psi);
	  matrix_tran(temp,m4,n4,stk(l4r),n4,m4);
      matrix_tran(tempi,m4,n4,stk(l4i),n4,m4);
	  free(temp);
	  free(tempi);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==REAL))
  {
	 
	  
	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	  temp = malloc(m4*n4*sizeof(double));
	  tempi = malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  //scale = (int)(floor(stk(l2)[count]));
		  plen=(int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count]))+1;
		  if (plen<3)
			  plen = 3;
	      psi=malloc(plen*sizeof(double));
	      full_range_scalef (cstk(l3), psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(stk(l2)[count]);
		  cwt_conv_real (stk(l1r), n4, psi, plen, temp+count*n4, n4);
		  cwt_conv_real (stk(l1i), n4, psi, plen, tempi+count*n4, n4);
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,stk(l4r),n4,m4);
	  matrix_tran(tempi,m4,n4,stk(l4i),n4,m4);
	  free(temp);
	  free(tempi);
  }

  if ((ind2!=-1) && (ci[ind2].realOrComplex==COMPLEX))
  {
	 
	  
	  CreateCVar(4, "d", &it, &m4, &n4, &l4r, &l4i);
	  temp = malloc(m4*n4*sizeof(double));
	  tempi = malloc(m4*n4*sizeof(double));
	  for (count=0;count<m4;count++)
	  {
		  plen=((int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count])+1))*2;//(int)((ci[ind2].ub-ci[ind2].lb)*(stk(l2)[count]+1))*2;
		  if (plen<3)
			  plen = 3;
	      psi=malloc(plen*sizeof(double));
	      full_range_scalef (cstk(l3), psi, plen);
          for(i=0;i<plen;i++)
			psi[i]=psi[i]/sqrt(stk(l2)[count]);
		  cwt_conv_complex_complex (stk(l1r), stk(l1i), n4, psi, psi+plen/2, plen/2, temp+count*n4, tempi+count*n4, n4); 
		  free(psi);
	  }
	  matrix_tran(temp,m4,n4,stk(l4r),n4,m4);
	  matrix_tran(tempi,m4,n4,stk(l4i),n4,m4);
	  free(temp);
	  free(tempi);
  }

      LhsVar(1) = 4;

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