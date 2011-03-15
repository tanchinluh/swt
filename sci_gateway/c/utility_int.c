/*
 * -------------------------------------------------------------------------
 * utility_int.c -- utility function interface
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2007  Roger Liu
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

/*-------------------------------------------
 * convolution
 *-----------------------------------------*/

int
int_conv (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int minlhs = 1, maxlhs = 1, minrhs = 2, maxrhs = 2;
  int errCode;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);
  /* Get Input */
  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  conv_validate (&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);	/* Indicate error  */
      return 0;			/* Stop when validation fails */
    }
  /* Memory Management */
  m3 = 1;
  n3 = m1 * n1 + m2 * n2 - 1;
  CreateVar (3, "d", &m3, &n3, &l3);
  /* Actual Processing */
  conv (stk (l1), n1 * m1, stk (l3), n3, stk (l2), n2 * m2);
  /* Return Value */
  LhsVar (1) = 3;
  return 0;
}

/*-------------------------------------------
 * circular or periodic convolution
 *-----------------------------------------*/

int
int_iconv (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int minlhs = 1, maxlhs = 1, minrhs = 2, maxrhs = 2;
  int errCode;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);
  /* Get Input */
  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "d", &m2, &n2, &l2);
  conv_validate (&errCode);

  if (errCode != SUCCESS)
    {
      validate_print (errCode);	
      return 0;			
    }
  /* Memory Management */
  m3 = 1;
  if (m1*n1>m2*n2)
      n3 = m1 * n1;
  else
	  n3 = m2 * n2;
  CreateVar (3, "d", &m3, &n3, &l3);
  /* Actual Processing */
  if (m1*n1>m2*n2)
     i_conv (stk (l1), n1 * m1, stk (l3), n3, stk (l2), n2 * m2);
  else
      i_conv (stk (l2), n2 * m2, stk (l3), n3, stk (l1), n1 * m1);
  /* Return Value */
  LhsVar (1) = 3;
  return 0;
}


/*-------------------------------------------
 * wrev
 *-----------------------------------------*/

int
int_wrev (char *fname)
{
  static int l1, m1, n1, l2, m2, n2;
  static int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 1;
  int errCode;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  wrev_validate (&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
  GetRhsVar (1, "d", &m1, &n1, &l1);

  m2 = m1;
  n2 = n1;
  CreateVar (2, "d", &m2, &n2, &l2);

  wrev (stk (l1), n1 * m1, stk (l2), n1 * m1);

  LhsVar (1) = 2;
  return 0;
}

/*-------------------------------------------
 * qmf
 *-----------------------------------------*/

int
int_qmf (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 2;
  int errCode, flow;
  /* Input Validation  */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);
  
  qmf_validate (&flow, &errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);	
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);

  switch (flow){
  case 1:
    {
      m2 = m1;
      n2 = n1;
      CreateVar (2, "d", &m2, &n2, &l2);
      qmf_even (stk (l1), n1 * m1, stk (l2), n1 * m1);
      LhsVar (1) = 2;
      break;
    }
  case 2:
    {
      GetRhsVar (2, "i", &m2, &n2, &l2);
      m3 = m1;
      n3 = n1;
      CreateVar (3, "d", &m3, &n3, &l3);
      if ((istk(l2)[0] % 2) == 0) /* even */
	qmf_even (stk (l1), n1 * m1, stk (l3), n1 * m1);
      else
	qmf_odd (stk (l1), n1 * m1, stk (l3), n1 * m1);
      LhsVar (1) = 3;
      break;
    }
  default:
    break;
  }
  return 0;
}


/*-------------------------------------------
 * dyaddown
 *-----------------------------------------*/

int
int_dyaddown (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 3;
  int flow, errCode;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  flow = 0;
  dyaddown_form_validate (&flow, &errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  switch (flow){
  case 1:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      if (m1 > 1)
	{
	  m2 = m1 / 2;// + (m1%2==0)?0:1;
	  n2 = 1;
	}
      else
	{
	  m2 = 1;
	  n2 = n1 / 2;// + (n1%2==0)?0:1;
	}
      CreateVar (2, "d", &m2, &n2, &l2);
      dyaddown_1D_keep_even (stk(l1), n1 * m1, stk(l2), n2 * m2);
      LhsVar (1) = 2;
      break;
    }
  case 2:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
    //  if (m1 > 1)
	//{
	 // m3 = m1 / 2;
	 // n3 = 1;
	//}
    //  else
	//{
	//  m3 = 1;
	//  n3 = n1 / 2;
	//}
    //  CreateVar (3, "d", &m3, &n3, &l3);
      if (istk(l2)[0] % 2 == 0)
	  {
		if (m1 > 1)
		{
		 m3 = m1 / 2;
		 n3 = 1;
		}
		else
		{
		  m3 = 1;
		  n3 = n1 / 2;
		}
		CreateVar (3, "d", &m3, &n3, &l3);
		dyaddown_1D_keep_even (stk(l1), n1 * m1, stk(l3), n3 * m3);
	  }
      else
	  {
		if (m1 > 1)
		{
			m3 = m1 / 2;// + ((m1%2)==0)?0:1;
			n3 = 1;
			if (m1 % 2 != 0)
				m3 += 1;
		}
		else
		{
		  m3 = 1;
		  n3 = n1 / 2;// + ((n1%2)==0)?0:1;
		  if (n1 % 2 != 0)
			  n3 += 1;
		}
		//sprintf("%d\n",m3);
		//sprintf("%d\n",n3);
		CreateVar (3, "d", &m3, &n3, &l3);
		dyaddown_1D_keep_odd (stk(l1), n1 * m1, stk(l3), n3 * m3);
	  }
      LhsVar (1) = 3;
      break;
    }
  case 3:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      //m2 = m1 / 2;
	  m2 = m1;
      n2 = n1 / 2;
      CreateVar (2, "d", &m2, &n2, &l2);
      //dyaddown_2D_keep_even (stk(l1), m1, n1, stk(l2), m2, n2); 
	  dyaddown_2D_keep_even_col(stk(l1), m1, n1, stk(l2), m2, n2);
      LhsVar (1) = 2;
      break;
    }
  case 4:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "c", &m2, &n2, &l2);
      //printf("%c\n",cstk(l2)[0]);
      dyaddown_content_validate(cstk(l2), &errCode);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if (strcmp(cstk(l2),"r")==0)
	{
	  m3 = m1 / 2;
	  n3 = n1;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyaddown_2D_keep_even_row(stk(l1), m1, n1, 
				    stk(l3), m3, n3); 
	}
      if (strcmp(cstk(l2),"c")==0)
	{
	  m3 = m1;
	  n3 = n1 / 2;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyaddown_2D_keep_even_col(stk(l1), m1, n1, 
				    stk(l3), m3, n3); 
	}
      if (strcmp(cstk(l2),"m")==0)
	{
	  m3 = m1 / 2;
	  n3 = n1 / 2;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyaddown_2D_keep_even(stk(l1), m1, n1, stk(l3), m3, n3); 
	}
      LhsVar (1) = 3;
      break;
    }
  case 5:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      //m3 = m1 / 2;
	  
      if (istk(l2)[0]%2==0)
	  {
 	    m3 = m1;
        n3 = n1 / 2;
        CreateVar (3, "d", &m3, &n3, &l3);
		dyaddown_2D_keep_even_col(stk(l1), m1, n1, stk(l3), m3, n3); 
	  }
      else
	  {
        m3 = m1;
		n3 = n1 /2;
		if (n1 % 2 != 0)
			n3 += 1;
        CreateVar (3, "d", &m3, &n3, &l3);
	    dyaddown_2D_keep_odd_col(stk(l1), m1, n1, stk(l3), m3, n3); 
	  }
      LhsVar (1) = 3;
      break;
    }
  case 6:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      GetRhsVar (3, "c", &m3, &n3, &l3);
      dyaddown_content_validate(cstk(l3), &errCode);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if ((istk(l2)[0])%2==0)
	{
	  if (strcmp(cstk(l3),"r")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_even_row(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 / 2;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_even_col(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"m")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1 / 2;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_even(stk(l1), m1, n1, 
				    stk(l4), m4, n4); 
	    }
	}
      else
	{
	  if (strcmp(cstk(l3),"r")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1;
		  if (m1 % 2 != 0)
			  m4 += 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_odd_row(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 / 2;
		  if (n1 % 2 != 0)
			  n4 += 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_odd_col(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"m")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1 / 2;
		  if (m1 % 2 != 0)
			  m4 += 1;
		  if (n1 % 2 != 0)
			  n4 += 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_odd(stk(l1), m1, n1, 
				   stk(l4), m4, n4); 
	    }
	}
      LhsVar (1) = 4;
      break;
    }
  case 7:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "i", &m3, &n3, &l3);
      dyaddown_content_validate(cstk(l2), &errCode);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if ((istk(l3)[0])%2==0)
	{
	  if (strcmp(cstk(l2),"r")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_even_row(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 / 2;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_even_col(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"m")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1 / 2;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_even(stk(l1), m1, n1, 
				    stk(l4), m4, n4); 
	    }
	}
      else
	{
	  if (strcmp(cstk(l2),"r")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1;
		  if (m1 % 2 != 0)
			  m4 += 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_odd_row(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 / 2;
		  if (n1 % 2 != 0)
			  n4 += 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_odd_col(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"m")==0)
	    {
	      m4 = m1 / 2;
	      n4 = n1 / 2;
		  if (m1 % 2 != 0)
			  m4 += 1;
		  if (n1 % 2 != 0)
			  n4 += 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyaddown_2D_keep_odd(stk(l1), m1, n1, 
				   stk(l4), m4, n4); 
	    }
	}
      LhsVar (1) = 4;
      break;
    }
  default:
    break;
  }
  return 0;
}

/*-------------------------------------------
 * dyadup
 *-----------------------------------------*/

int
int_dyadup (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 3;
  int flow, errCode;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  flow = 0;
  dyadup_form_validate (&flow, &errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

  switch (flow){
  case 1:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      if (m1 > 1)
	{
	  m2 = m1 * 2 + 1;
	  n2 = 1;
	}
      else
	{
	  m2 = 1;
	  n2 = n1 * 2 + 1;
	}
      CreateVar (2, "d", &m2, &n2, &l2);
      dyadup_1D_feed_even (stk(l1), n1 * m1, stk(l2), n2 * m2);
      LhsVar (1) = 2;
      break;
    }
  case 2:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      if (istk(l2)[0] % 2 == 0)
	{
	  if (m1 > 1)
	    {
	      m3 = m1 * 2 - 1;
	      n3 = 1;
	    }
	  else
	    {
	      m3 = 1;
	      n3 = n1 * 2 - 1;
	    }
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyadup_1D_feed_odd (stk(l1), n1 * m1, stk(l3), n3 * m3);
	}
      else
	{
	  if (m1 > 1)
	    {
	      m3 = m1 * 2 + 1;
	      n3 = 1;
	    }
	  else
	    {
	      m3 = 1;
	      n3 = n1 * 2 + 1;
	    }
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyadup_1D_feed_even (stk(l1), n1 * m1, stk(l3), n3 * m3);
	}
      LhsVar (1) = 3;
      break;
    }
  case 3:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      //m2 = m1 * 2 + 1;
      //n2 = n1 * 2 + 1;
	  m2 = m1;
	  n2 = n1 * 2 + 1;
      CreateVar (2, "d", &m2, &n2, &l2);
	  dyadup_2D_feed_even_col (stk(l1), m1, n1, stk(l2), m2, n2);
      //dyadup_2D_feed_even (stk(l1), m1, n1, stk(l2), m2, n2); 
      LhsVar (1) = 2;
      break;
    }
  case 4:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "c", &m2, &n2, &l2);
      dyadup_content_validate(cstk(l2), &errCode);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if (strcmp(cstk(l2),"r")==0)
	{
	  m3 = m1 * 2 + 1;
	  n3 = n1;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyadup_2D_feed_even_row(stk(l1), m1, n1, 
				  stk(l3), m3, n3); 
	}
      if (strcmp(cstk(l2),"c")==0)
	{
	  m3 = m1;
	  n3 = n1 * 2 + 1;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyadup_2D_feed_even_col(stk(l1), m1, n1, 
				    stk(l3), m3, n3); 
	}
      if (strcmp(cstk(l2),"m")==0)
	{
	  m3 = m1 * 2 + 1;
	  n3 = n1 * 2 + 1;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  dyadup_2D_feed_even(stk(l1), m1, n1, stk(l3), m3, n3); 
	}
      LhsVar (1) = 3;
      break;
    }
  case 5:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      if (istk(l2)[0]%2==0)
	{
	  //m3 = m1 * 2 - 1;
	  //n3 = n1 * 2 - 1;
	  m3 = m1;
	  n3 = n1 * 2 - 1;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  //dyadup_2D_feed_odd(stk(l1), m1, n1, stk(l3), m3, n3); 
	  dyadup_2D_feed_odd_col(stk(l1), m1, n1, stk(l3), m3, n3);
	}
      else
	{
	  //m3 = m1 * 2 + 1;
	  //n3 = n1 * 2 + 1;
	  m3 = m1;
	  n3 = n1 * 2 + 1;
	  CreateVar (3, "d", &m3, &n3, &l3);
	  //dyadup_2D_feed_even(stk(l1), m1, n1, stk(l3), m3, n3); 
	  dyadup_2D_feed_even_col(stk(l1), m1, n1, stk(l3), m3, n3);
	}
      LhsVar (1) = 3;
      break;
    }
  case 6:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      GetRhsVar (3, "c", &m3, &n3, &l3);
      dyadup_content_validate(cstk(l3), &errCode);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if ((istk(l2)[0])%2==0)
	{
	  if (strcmp(cstk(l3),"r")==0)
	    {
	      m4 = m1 * 2 - 1;
	      n4 = n1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_odd_row(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 * 2 - 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_odd_col(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"m")==0)
	    {
	      m4 = m1 * 2 - 1;
	      n4 = n1 * 2 - 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_odd(stk(l1), m1, n1, 
				    stk(l4), m4, n4); 
	    }
	}
      else
	{
	  if (strcmp(cstk(l3),"r")==0)
	    {
	      m4 = m1 * 2 + 1;
	      n4 = n1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_even_row(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 * 2 + 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_even_col(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l3),"m")==0)
	    {
	      m4 = m1 * 2 + 1;
	      n4 = n1 * 2 + 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_even(stk(l1), m1, n1, 
				   stk(l4), m4, n4); 
	    }
	}
      LhsVar (1) = 4;
      break;
    }
  case 7:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "i", &m3, &n3, &l3);
      dyadup_content_validate(cstk(l2), &errCode);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if ((istk(l3)[0])%2==0)
	{
	  if (strcmp(cstk(l2),"r")==0)
	    {
	      m4 = m1 * 2 - 1;
	      n4 = n1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_odd_row(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 * 2 - 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_odd_col(stk(l1), m1, n1, 
					stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"m")==0)
	    {
	      m4 = m1 * 2 - 1;
	      n4 = n1 * 2 - 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_odd(stk(l1), m1, n1, 
				    stk(l4), m4, n4); 
	    }
	}
      else
	{
	  if (strcmp(cstk(l2),"r")==0)
	    {
	      m4 = m1 * 2 + 1;
	      n4 = n1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_even_row(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"c")==0)
	    {
	      m4 = m1;
	      n4 = n1 * 2 + 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_even_col(stk(l1), m1, n1, 
				       stk(l4), m4, n4); 
	    }
	  if (strcmp(cstk(l2),"m")==0)
	    {
	      m4 = m1 * 2 + 1;
	      n4 = n1 * 2 + 1;
	      CreateVar (4, "d", &m4, &n4, &l4);
	      dyadup_2D_feed_even(stk(l1), m1, n1, 
				   stk(l4), m4, n4); 
	    }
	}
      LhsVar (1) = 4;
      break;
    }
  default:
    break;
  }
  return 0;
}

/*-------------------------------------------
 * wkeep
 *-----------------------------------------*/

int
int_wkeep (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int minlhs = 1, maxlhs = 1, minrhs = 2, maxrhs = 3;
  int flow, errCode;//, row1, row2, col1, col2;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  flow = 0;
  wkeep_form_validate (&flow, &errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }
  l1 = 0;
  l2 = 0;
  l3 = 0;
  switch (flow){
  case 1:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      wkeep_content_validate(flow, &errCode, l1, l2, l3);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if (n1 > 1)
	{
	  m3 = 1;
	  n3 = istk (l2)[0];
	}
      else
	{
	  m3 = istk (l2)[0];
	  n3 = 1;
	}
      CreateVar (3, "d", &m3, &n3, &l3);
      wkeep_1D_center(stk(l1), n1 * m1, stk(l3), n3 * m3);
      LhsVar (1) = 3;
      break;
    }
  case 2:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      wkeep_content_validate(flow, &errCode, l1, l2, l3);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}            
      m3 = istk(l2)[0];
      n3 = istk(l2)[1];
      CreateVar (3, "d", &m3, &n3, &l3);
      wkeep_2D_center(stk(l1), m1, n1, stk(l3), m3, n3);
      LhsVar (1) = 3;
      break;
    }
  case 3:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      GetRhsVar (3, "c", &m3, &n3, &l3);
      wkeep_content_validate(flow, &errCode, l1, l2, l3);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if (n1 > 1)
	{
	  m4 = 1;
	  n4 = istk (l2)[0];
	}
      else
	{
	  m4 = istk (l2)[0];
	  n4 = 1;
	}
      CreateVar (4, "d", &m4, &n4, &l4);
      if (((cstk(l3)[0] - 'l') == 0)||((cstk(l3)[0] - 'L') == 0))
	wkeep_1D_left(stk(l1), n1 * m1, stk(l4), n4 * m4);
      else if (((cstk(l3)[0] - 'c') == 0)||((cstk(l3)[0] - 'C') == 0))
	wkeep_1D_center(stk(l1), n1 * m1, stk(l4), n4 * m4);
      else if (((cstk(l3)[0] - 'r') == 0)||((cstk(l3)[0] - 'R') == 0))
	wkeep_1D_right(stk(l1), n1 * m1, stk(l4), n4 * m4);
      LhsVar (1) = 4;
      break;
    }
  case 4:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      GetRhsVar (3, "i", &m3, &n3, &l3);
      wkeep_content_validate(flow, &errCode, l1, l2, l3);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      if (n1 > 1)
	{
	  m4 = 1;
	  n4 = istk (l2)[0];
	}
      else
	{
	  m4 = istk (l2)[0];
	  n4 = 1;
	}
      CreateVar (4, "d", &m4, &n4, &l4);
      wkeep_1D_index(stk(l1),n1 * m1,stk(l4),n4 * m4,istk(l3)[0]);
      LhsVar (1) = 4;
      break;
    }
  case 5:
    {
      GetRhsVar (1, "d", &m1, &n1, &l1);
      GetRhsVar (2, "i", &m2, &n2, &l2);
      GetRhsVar (3, "i", &m3, &n3, &l3);
      wkeep_content_validate(flow, &errCode, l1, l2, l3);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;	
	}
      m4 = istk(l2)[0];
      n4 = istk(l2)[1];
      CreateVar (4, "d", &m4, &n4, &l4);
      wkeep_2D_index(stk(l1), m1, n1, stk(l4), 
		     m4, n4, istk(l3)[0], istk(l3)[1]);
      LhsVar (1) = 4;
      break;
    }
  default:
    break;
  }
  return 0;
}

/*-------------------------------------------
 * wextend
 *-----------------------------------------*/

int
int_wextend (char *fname)
{

  static int l1, m1, n1, l2, m2, n2, m3, n3, l3;
  static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
  static int minlhs = 1, maxlhs = 1, minrhs = 4, maxrhs = 5;
  int flow, errCode;// adr, po, row;//col;
  char c = 'b';
  char cr[2] = "b";
  char **str;
  extend_method extMethod;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  if (GetType(1) == sci_matrix)
    GetRhsVar (1, "i", &m1, &n1, &l1);
  if (GetType(1) == sci_strings)
    GetRhsVar (1, "c", &m1, &n1, &l1);

  flow = 0;
  //sciprint("before form validate\n");
  wextend_form_validate (&flow, &errCode, l1);
  //sciprint("flow=%d\n",flow);
  if (errCode != SUCCESS)
    {
	  //sciprint("form validate error!\n");
      validate_print (errCode);
      return 0;			
    }
  l2 = 0;
  l3 = 0;
  l4 = 0;
  l5 = 0;
  
  
  switch(flow) {
  case 1:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      GetRhsVar (5, "c", &m5, &n5, &l5);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && ((m3*n3)%2 != 0))
	{
	  if (((cstk(l5)[0]-'l') == 0) || ((cstk(l5)[0]-'r') == 0))
	    {
	      if (m3 > 1)
		{
		  m6 = m3 + istk(l4)[0] + 1;
		  n6 = 1;
		}
	      else
		{
		  m6 = 1;
		  n6 = n3 + istk(l4)[0] + 1;
		}
	    }
	  if ((cstk(l5)[0]-'b') == 0)
	    {
	      if (m3 > 1)
		{
		  m6 = m3 + 2*istk(l4)[0] + 1;
		  n6 = 1;
		}
	      else
		{
		  m6 = 1;
		  n6 = n3 + 2*istk(l4)[0] + 1;
		}
	    }
	}
      else
	{
	  if (((cstk(l5)[0]-'l') == 0) || ((cstk(l5)[0]-'r') == 0))
	    {
	      if (m3 > 1)
		{
		  m6 = m3 + istk(l4)[0];
		  n6 = 1;
		}
	      else
		{
		  m6 = 1;
		  n6 = n3 + istk(l4)[0];
		}
	    }
	  if ((cstk(l5)[0]-'b') == 0)
	    {
	      if (m3 > 1)
		{
		  m6 = m3 + 2*istk(l4)[0];
		  n6 = 1;
		}
	      else
		{
		  m6 = 1;
		  n6 = n3 + 2*istk(l4)[0];
		}
	    }
	}
      CreateVar (6, "d", &m6, &n6, &l6);
      if (strcmp(cstk(l5),"l") == 0)
	wextend_1D_left(stk(l3), m3*n3, stk(l6), m6*n6, extMethod);
      else if (strcmp(cstk(l5),"r") == 0)
	wextend_1D_right(stk(l3), m3*n3, stk(l6), m6*n6, extMethod);
      else
	wextend_1D_center(stk(l3), m3*n3, stk(l6), m6*n6, extMethod);
      LhsVar (1) = 6;
      break;
    }
  case 2:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && ((m3*n3)%2 != 0))
	{
	  if (m3 > 1)
	    {
	      m5 = m3 + 2*istk(l4)[0] + 1;
	      n5 = 1;
	    }
	  else
	    {
	      m5 = 1;
	      n5 = n3 + 2*istk(l4)[0] + 1;
	    }
	}
      else
	{
	  if (m3 > 1)
	    {
	      m5 = m3 + 2*istk(l4)[0];
	      n5 = 1;
	    }
	  else
	    {
	      m5 = 1;
	      n5 = n3 + 2*istk(l4)[0];
	    }
	}
      CreateVar (5, "d", &m5, &n5, &l5);
      wextend_1D_center(stk(l3), m3*n3, stk(l5), m5*n5, extMethod);
      LhsVar (1) = 5;
      break;
    }
  case 3:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      GetRhsVar (5, "S", &m5, &n5, &str);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (m3%2 != 0))
	{
	  if (((str[0][0]-'l') == 0) || ((str[0][0]-'r') == 0))
	    m6 = m3 + istk(l4)[0] + 1;
	  if ((str[0][0]-'b') == 0)
	    m6 = m3 + 2*istk(l4)[0] + 1;
	}
      else
	{
	  if (((str[0][0]-'l') == 0) || ((str[0][0]-'r') == 0))
	    m6 = m3 + istk(l4)[0];
	  if ((str[0][0]-'b') == 0)
	    m6 = m3 + 2*istk(l4)[0];
	}
      if ((extMethod == PER) && (n3%2 != 0))
	{
	  if (((str[1][0]-'l') == 0) || ((str[1][0]-'r') == 0))
	    n6 = n3 + istk(l4)[1] + 1;
	  if ((str[1][0]-'b') == 0)
	    n6 = n3 + 2*istk(l4)[1] + 1;
	}
      else
	{
	  if (((str[1][0]-'l') == 0) || ((str[1][0]-'r') == 0))
	    n6 = n3 + istk(l4)[1];
	  if ((str[1][0]-'b') == 0)
	    n6 = n3 + 2*istk(l4)[1];
	}
      CreateVar (6, "d", &m6, &n6, &l6);
      //wextend_2D (stk(l3), m3, n3, stk(l6), m6, n6, 
		//  extMethod, &str[0][0], &str[1][0]);
	  wextend_2D (stk(l3), m3, n3, stk(l6), m6, n6, 
		  extMethod, &str[1][0], &str[0][0]);
      LhsVar (1) = 6;
      FreeRhsSVar(str);
      break;
    }
  case 4:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (m3%2 != 0))
	m5 = m3 + 2*istk(l4)[0] + 1;
      else
	m5 = m3 + 2*istk(l4)[0];
      if ((extMethod == PER) && (n3%2 != 0))
	n5 = n3 + 2*istk(l4)[0] + 1;
      else
	n5 = n3 + 2*istk(l4)[0];
      CreateVar (5, "d", &m5, &n5, &l5);
      wextend_2D (stk(l3), m3, n3, stk(l5), m5, n5, 
		  extMethod, &c, &c);
      LhsVar (1) = 5;
      break;
    }
  case 5:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (m3%2 != 0))
	m5 = m3 + 2*istk(l4)[0] + 1;
      else
	m5 = m3 + 2*istk(l4)[0];
      if ((extMethod == PER) && (n3%2 != 0))
	n5 = n3 + 2*istk(l4)[1] + 1;
      else
	n5 = n3 + 2*istk(l4)[1];
      CreateVar (5, "d", &m5, &n5, &l5);
      wextend_2D (stk(l3), m3, n3, stk(l5), m5, n5, 
		  extMethod, &c, &c);
      LhsVar (1) = 5;
      break;
    }
  case 6:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      GetRhsVar (5, "c", &m5, &n5, &l5);
	  //sciprint("before wextend content validate!\n");
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
	  //sciprint("after wextend content validate!\n");
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (m3%2 != 0))
	{
	  if ((strcmp(cstk(l5),"l")==0) || (strcmp(cstk(l5),"r")==0))
	    m6 = m3 + istk(l4)[0] + 1;
	  if ((strcmp(cstk(l5),"b") == 0))
	    m6 = m3 + 2*istk(l4)[0] + 1;
	}
      else
	{
	  if ((strcmp(cstk(l5),"l")==0) || (strcmp(cstk(l5),"r")==0))
	    m6 = m3 + istk(l4)[0];
	  if ((strcmp(cstk(l5),"b") == 0))
	    m6 = m3 + 2*istk(l4)[0];
	}
      n6 = n3;
      CreateVar (6, "d", &m6, &n6, &l6);
	  //sciprint("after creation!flow=6\n");
      wextend_2D_row (stk(l3), m3, n3, stk(l6), m6, n6, 
		  extMethod, cstk(l5));
      LhsVar (1) = 6;
      break;
    }
  case 7:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (m3%2 != 0))
	m5 = m3 + 2*istk(l4)[0] + 1;
      else
	m5 = m3 + 2*istk(l4)[0];
      n5 = n3;
      CreateVar (5, "d", &m5, &n5, &l5);
      wextend_2D_row (stk(l3), m3, n3, stk(l5), m5, n5, 
		  extMethod, &cr);
      LhsVar (1) = 5;
      break;
    }
  case 8:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      GetRhsVar (5, "c", &m5, &n5, &l5);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (n3%2 != 0))
	{
	  if ((strcmp(cstk(l5),"l")==0) || (strcmp(cstk(l5),"r")==0))
	    n6 = n3 + istk(l4)[0] + 1;
	  if ((strcmp(cstk(l5),"b") == 0))
	    n6 = n3 + 2*istk(l4)[0] + 1;
	}
      else
	{
	  if ((strcmp(cstk(l5),"l")==0) || (strcmp(cstk(l5),"r")==0))
	    n6 = n3 + istk(l4)[0];
	  if ((strcmp(cstk(l5),"b") == 0))
	    n6 = n3 + 2*istk(l4)[0];
	}
      m6 = m3;
      CreateVar (6, "d", &m6, &n6, &l6);
      wextend_2D_col (stk(l3), m3, n3, stk(l6), m6, n6, 
		  extMethod, cstk(l5));
      LhsVar (1) = 6;
      break;
    }
  case 9:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (n3%2 != 0))
	n5 = n3 + 2*istk(l4)[0] + 1;
      else
	n5 = n3 + 2*istk(l4)[0];
      m5 = m3;
      CreateVar (5, "d", &m5, &n5, &l5);
      wextend_2D_col (stk(l3), m3, n3, stk(l5), m5, n5, 
		  extMethod, &cr);
      LhsVar (1) = 5;
      break;
    }
  case 10:
    {
      GetRhsVar (2, "c", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      GetRhsVar (4, "i", &m4, &n4, &l4);
      GetRhsVar (5, "c", &m5, &n5, &l5);
      wextend_content_validate (flow, &errCode, l2, l3, l4, l5, str);
      if (errCode != SUCCESS)
	{
	  validate_print (errCode);
	  return 0;			
	}
      extend_method_parse (cstk(l2), &extMethod);
      if ((extMethod == PER) && (m3%2 != 0))
	{
	  if ((cstk(l5)[0]=='l') || (cstk(l5)[0]=='r'))
	    m6 = m3 + istk(l4)[0] + 1;
	  if (cstk(l5)[0]=='b')
	    m6 = m3 + 2*istk(l4)[0] + 1;
	}
      else
	{
	  if ((cstk(l5)[0]=='l') || (cstk(l5)[0]=='r'))
	    m6 = m3 + istk(l4)[0];
	  if (cstk(l5)[0]=='b')
	    m6 = m3 + 2*istk(l4)[0];
	}
      if ((extMethod == PER) && (n3%2 != 0))
	{
	  if ((cstk(l5)[1]=='l') || (cstk(l5)[1]=='r'))
	    n6 = n3 + istk(l4)[1] + 1;
	  if (cstk(l5)[1]=='b')
	    n6 = n3 + 2*istk(l4)[1] + 1;
	}
      else
	{
	  if ((cstk(l5)[1]=='l') || (cstk(l5)[1]=='r'))
	    n6 = n3 + istk(l4)[1];
	  if (cstk(l5)[1]=='b')
	    n6 = n3 + 2*istk(l4)[1];
	}

	
      CreateVar (6, "d", &m6, &n6, &l6);
      wextend_2D (stk(l3), m3, n3, stk(l6), m6, n6, 
		  extMethod, &cstk(l5)[1], &cstk(l5)[0]);
	  //wextend_2D (stk(l3), m3, n3, stk(l6), m6, n6, 
		//  extMethod, &cstk(l5)[0], &cstk(l5)[1]);
      LhsVar (1) = 6;
      break;
    }
  default:
    break;
  }
  //sciprint("flow=%d\n",flow);
  return 0;
}

/*-------------------------------------------
 * wcodemat
 *-----------------------------------------*/
int_wcodemat (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int l4, m4, n4, l5, m5, n5;
  static int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 4;
  int errCode, flow, op, zero, typ, mn, inc;
  double *var;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  wcodemat_form_validate (&flow, &errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);	
      return 0;			
    }
  l1 = 0;
  l2 = 0;
  l3 = 0;
  l4 = 0;
  typ = I_UINT16;
  zero = 0;
  inc = 1;

  switch (flow){
  case 1:
	  {
          GetRhsVar (1, "d", &m1, &n1, &l1);
		  wcodemat_content_validate (&errCode, flow, l1, l2, l3, l4);
          if (errCode != SUCCESS)
            {
              validate_print (errCode);	
              return 0;			
            }
		  m2 = m1;
		  n2 = n1;
		  var=malloc(m2*n2*sizeof(double));
		  l2 = typ = I_UINT16;
          CreateVar (2, "I", &m2, &n2, &l2);
		  wcodemat_matrix (stk(l1), m1, n1, var, m2, n2, 1, 64, 1);
		  mn = m2*n2;
		  C2F(tpconv)(&zero,&typ,&mn,var, &inc,istk(l2), &inc);
		  free(var);
		  LhsVar(1) = 2;
		  break;
	  }
  case 2:
	  {
          GetRhsVar (1, "d", &m1, &n1, &l1);
          GetRhsVar (2, "i", &m2, &n2, &l2);
		  wcodemat_content_validate (&errCode, flow, l1, l2, l3, l4);
          if (errCode != SUCCESS)
            {
              validate_print (errCode);	
              return 0;			
            }
		  m3 = m1;
		  n3 = n1;
		  l3 = typ = I_UINT16;
		  CreateVar (3, "I", &m3, &n3, &l3);
		  var = malloc(m3*n3*sizeof(double));
		  wcodemat_matrix (stk(l1), m1, n1, var, m3, n3, 1, istk(l2)[0], 1);
		  mn = m3*n3;
		  C2F(tpconv)(&zero,&typ,&mn,var, &inc,istk(l3), &inc);
		  free(var);
		  LhsVar(1) = 3;
		  break;
	  }
  case 3:
	  {
          GetRhsVar (1, "d", &m1, &n1, &l1);
          GetRhsVar (2, "i", &m2, &n2, &l2);
          GetRhsVar (3, "c", &m3, &n3, &l3);
		  wcodemat_content_validate (&errCode, flow, l1, l2, l3, l4);
          if (errCode != SUCCESS)
            {
              validate_print (errCode);	
              return 0;			
            }
		  if ((!strcmp(cstk(l3),"r"))  || (!strcmp(cstk(l3),"row"))||
	                (!strcmp(cstk(l3),"c"))  ||	(!strcmp(cstk(l3),"col"))||
	                (!strcmp(cstk(l3),"m"))  || (!strcmp(cstk(l3),"mat")))
		  {

			  m4 = m1;
		      n4 = n1;
			  var=malloc(m4*n4*sizeof(double));
			  l4 = typ = I_UINT16;
		      CreateVar (4, "I", &m4, &n4, &l4);
			  if ((!strcmp(cstk(l3),"c")) || (!strcmp(cstk(l3),"col")))
				  wcodemat_matrix_col (stk(l1), m1, n1, var, m4, n4, 1, istk(l2)[0], 1);
			  else if ((!strcmp(cstk(l3),"r")) || (!strcmp(cstk(l3),"row")))
				  wcodemat_matrix_row (stk(l1), m1, n1, var, m4, n4, 1, istk(l2)[0], 1);
			  else if ((!strcmp(cstk(l3),"m")) || (!strcmp(cstk(l3),"mat")))
				  wcodemat_matrix (stk(l1), m1, n1, var, m4, n4, 1, istk(l2)[0], 1);
			  mn = m4*n4;
              C2F(tpconv)(&zero,&typ,&mn,var, &inc,istk(l4), &inc);
			  free(var);
		      LhsVar(1) = 4;
		  }
		  else
		  {
			  sciprint("option argument not valid!\n");
			  return 0;
		  }
		  break;
	  }
  case 4:
	  {
		  GetRhsVar (1, "d", &m1, &n1, &l1);
          GetRhsVar (2, "i", &m2, &n2, &l2);
          GetRhsVar (3, "c", &m3, &n3, &l3);
          GetRhsVar (4, "i", &m4, &n4, &l4);
		  wcodemat_content_validate (&errCode, flow, l1, l2, l3, l4);
          if (errCode != SUCCESS)
            {
              validate_print (errCode);	
              return 0;			
            }
		  if ((!strcmp(cstk(l3),"r"))  || (!strcmp(cstk(l3),"row"))||
	                (!strcmp(cstk(l3),"c"))  ||	(!strcmp(cstk(l3),"col"))||
	                (!strcmp(cstk(l3),"m"))  || (!strcmp(cstk(l3),"mat")))
		  {
		      m5 = m1;
		      n5 = n1;
			  l5 = typ = I_UINT16;
		      CreateVar (5, "I", &m5, &n5, &l5);
			  if (istk(l4)[0] !=0)
				  op = 1;
			  else
				  op = 0;
			  var = malloc(m5*n5*sizeof(double));
			  if ((!strcmp(cstk(l3),"c")) || (!strcmp(cstk(l3),"col")))
				  wcodemat_matrix_col (stk(l1), m1, n1, var, m5, n5, 1, istk(l2)[0], op);
			  else if ((!strcmp(cstk(l3),"r")) || (!strcmp(cstk(l3),"row")))
				  wcodemat_matrix_row (stk(l1), m1, n1, var, m5, n5, 1, istk(l2)[0], op);
			  else if ((!strcmp(cstk(l3),"m")) || (!strcmp(cstk(l3),"mat")))
				  wcodemat_matrix (stk(l1), m1, n1, var, m5, n5, 1, istk(l2)[0], op);
			  mn = m5*n5;
			  C2F(tpconv)(&zero,&typ,&mn,var, &inc,istk(l5), &inc);
			  free(var);
		      LhsVar(1) = 5;
		  }
		  else
		  {
			  sciprint("option argument not valid!\n");
			  return 0;
		  }
		  break;
	  }
  default:
	  break;
    }
  
  return 0;
}

int_ind2rgb (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int minlhs = 1, maxlhs = 1, minrhs = 2, maxrhs = 2;
  int mL, nL, ms, ns, lc2, it, mi, ni, mn, inc, zero, count;
  int si[3];
  static char *Str[]= { "hm","dims","entries"};
  double *var, *temp;
  SciIntMat ssi, M;

  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  if ((GetType(1)!=sci_ints)|| (GetType(2)!=sci_matrix))
  {
	  sciprint("Argument 1 should be integer matrix and 2 should be Nx3 double matrix");
      return 0;
  }
  GetRhsVar (1, "I", &m1, &n1, &M);
  mn = m1*n1;
  inc  = 1;
  zero = 0;
  //C2F(tpconv)(&M.it,&zero,&mn, M.D, &inc, stk(l1), &inc);

  if (M.it==1)
  {
	  sciprint("index matrix should be real integer!\n");
	  return 0;
  }
    
  GetRhsCVar(2, "d", &it, &m2, &n2, &l2, &lc2);
  
  if (it==1)
  {
	  sciprint("colormap should be real matrix!\n");
	  return 0;
  }
  if ((n2!=3))
  {
	  sciprint("colormap should be Nx3 matrix!\n");
      return 0;
  }

  si[0] = m1;
  si[1] = n1;
  si[2] = 3;
  ssi.m = 1;
  ssi.n = 3;
  ssi.l = 100;
  ssi.it = 4;
  ssi.D = si;
  var = malloc(m1*n1*3*sizeof(double));
  temp = malloc(m1*n1*sizeof(double));


  C2F(tpconv)(&M.it,&zero,&mn, M.D, &inc, temp, &inc);
  for(count=0;count<m1*n1;count++)
  {
     if ((int)(temp[count]) < 1)
	 {
         var[count] = stk(l2)[0];
		 var[count+m1*n1] = stk(l2)[m2];
		 var[count+2*m1*n1] = stk(l2)[2*m2];
	 }
	 if ((int)(temp[count]) < m2)
	 {
		 var[count]=stk(l2)[(int)(temp[count])-1];
		 var[count+m1*n1]=stk(l2)[(int)(temp[count])+m2-1]; 
		 var[count+m1*n1*2]=stk(l2)[(int)(temp[count])+2*m2-1];
	 }
	 else
	 {
		 var[count] = stk(l2)[m2-1];
		 var[count+m1*n1] = stk(l2)[2*m2-1];
		 var[count+2*m1*n1] = stk(l2)[3*m2-1];
	 }
  }
  free(temp);
  mi = 1;
  ni = 3;
  ms = 1;
  ns = 3;
  mL = 3;
  nL = 1;
  m3 = m1*n1*3;
  n3 = 1;
  CreateVar(3, "m", &mL, &nL, &l3);
  CreateListVarFromPtr(3,1,"S",&ms,&ns,Str);
  CreateListVarFromPtr(3,2,"I",&mi,&ni,&ssi);
  CreateListVarFromPtr(3,3,"d",&m3, &n3, &var);
  free(var);
  LhsVar(1) = 3;
 
  return 0;
}


int
int_mat3Dtran (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, mL2, nL2, s4;
  static int it3, mL3, nL3, lL3, lcL3, mL1, nL1, s3, l4, m4, n4;
  static int minlhs = 1, maxlhs = 1, minrhs = 3, maxrhs = 3;
  char **Str2;
  static char *Str[]= { "hm","dims","entries"};
  int mL=3, nL=1, ms=1, ns=3, mi=1, ni=3;
  int r, c, s, zero, inc, mn, row, col, sli, r2;
  SciIntMat ssi, M;
  int si[3];
  double *temp, *var3, *var2;
  int errCode;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);
  /* Get Input */
  GetRhsVar(1,"m",&m1,&n1,&l1);
  CheckLength(1,m1,3);
  GetListRhsVar(1,1,"S",&mL1,&nL1,&Str2);

  if ( strcmp(Str2[0],"hm") != 0) 
    {
      Scierror(999,"Argument %d is not an hypermatrix\r\n",1);
      return 0;
    }
  FreeRhsSVar(Str2);
  GetListRhsVar(1,2,"I",&mL2,&nL2,&M);
  GetListRhsCVar(1,3,"d",&it3,&mL3,&nL3,&lL3,&lcL3);
  
  GetRhsVar (2, "i", &m2, &n2, &l2);
  if ((istk(l2)[0] != 1) && (istk(l2)[0] != 2) && (istk(l2)[0] != 3) &&
      (istk(l2)[0] != 4) && (istk(l2)[0] != 5) && (istk(l2)[0] != 0))
    {
      sciprint("the second argument should be integer from 1 to 6!\n");
      return 0;
    }

  GetRhsVar (3, "i", &m3, &n3, &l3);
  if ((istk(l3)[0] != 0) && (istk(l3)[0] != 1) )
    {
      sciprint("the second argument should be integer 1 or 0!\n");
      return 0;
    }

  if (it3 == 1)
    {
      Scierror(999,"Argument %d should be real hypermatrix\r\n",1);
      return 0;
    }
  if ((mL2 != 1) || (nL2 != 3))
    {
      Scierror(999,"Argument %d dimension error\r\n",1);
      return 0;
    }
  
  mn = mL2*nL2;
  inc  = 1;
  zero = 0;

  temp = malloc(m1*n1*sizeof(double));
  C2F(tpconv)(&M.it,&zero,&mn, M.D, &inc, temp, &inc);
  row = (int)temp[0];
  col = (int)temp[1];
  sli = (int)temp[2];
  if (istk(l3)[0]==0)
    {
      switch (istk(l2)[0])
	{
	case 0:
	  {
	    m4 = (int)temp[0];
	    n4 = (int)temp[1];
	    s4 = (int)temp[2];
	    break;
	  }
	case 1:
	  {
	    m4 = (int)temp[1];
	    n4 = (int)temp[0];
	    s4 = (int)temp[2];
	    break;
	  }
	case 2:
	  {
	    m4 = (int)temp[0];
	    n4 = (int)temp[2];
	    s4 = (int)temp[1];
	    break;
	  }
	case 3:
	  {
	    m4 = (int)temp[2];
	    n4 = (int)temp[0];
	    s4 = (int)temp[1];
	    break;
	  }
	case 4:
	  {
	    m4 = (int)temp[1];
	    n4 = (int)temp[2];
	    s4 = (int)temp[0];
	    break;
	  }
	case 5:
	  {
	    m4 = (int)temp[2];
	    n4 = (int)temp[1];
	    s4 = (int)temp[0];
	    break;
	  }
	default:break;
	}
    }
  else
    {

      switch (istk(l2)[0])
	{
	case 0:
	  {
	    m4 = (int)temp[0];
	    n4 = (int)temp[1];
	    s4 = (int)temp[2];
	    break;
	  }
	case 1:
	  {
	    m4 = (int)temp[1];
	    n4 = (int)temp[0];
	    s4 = (int)temp[2];
	    break;
	  }
	case 2:
	  {
	    m4 = (int)temp[0];
	    n4 = (int)temp[2];
	    s4 = (int)temp[1];
	    break;
	  }
	case 3:
	  {
	    m4 = (int)temp[1];
	    n4 = (int)temp[2];
	    s4 = (int)temp[0];
	    break;
	  }
	case 4:
	  {
	    m4 = (int)temp[2];
	    n4 = (int)temp[0];
	    s4 = (int)temp[1];
	    break;
	  }
	case 5:
	  {
	    m4 = (int)temp[2];
	    n4 = (int)temp[1];
	    s4 = (int)temp[0];
	    break;
	  }
	default:break;
	}
    }

  si[0] = m4;
  si[1] = n4;
  si[2] = s4;
  //si[3] = 8;
  ssi.m = 1;
  ssi.n = 3;
  ssi.it = 4;
  ssi.l = 100;
  ssi.D = si;
  m4 = m4*n4*s4;
  n4 = 1;
  var3 = malloc(m4*n4*sizeof(double));
  
  if (istk(l3)[0]==0)
    {
      switch (istk(l2)[0])
	{
	case 0:
	  {
	    verbatim_copy (stk(lL3), row*col*sli, var3, row*col*sli);
	    break;
	  }
	case 1:
	  {
	    dwt3d_tran(stk(lL3), col, row, sli, var3, row, col, sli);
	    break;
	  }
	case 2:
	  {
	    var2 = malloc(m4*n4*sizeof(double));
	    dwt3d_tran_z(stk(lL3), col, row, sli, var2, row, sli, col);
	    dwt3d_tran(var2, row, sli, col, var3, sli, row, col);
	    free(var2);
	    break;
	  }
	case 3:
	  {
	    dwt3d_tran_z(stk(lL3), col, row, sli, var3, row, sli, col);
	    break;
	  }
	case 4:
	  {
	    var2 = malloc(m4*n4*sizeof(double));
	    dwt3d_tran(stk(lL3), col, row, sli, var3, row, col, sli);
	    dwt3d_tran_z(var3, row, col, sli, var2, col, sli, row);
	    dwt3d_tran(var2, col, sli, row, var3, sli, col, row);
	    free(var2);
	    break;
	  }
	case 5:
	  {
	    var2 = malloc(m4*n4*sizeof(double));
	    dwt3d_tran(stk(lL3), col, row, sli, var2, row, col, sli);
	    dwt3d_tran_z(var2, row, col, sli, var3, col, sli, row);
	    free(var2);
	    break;
	  }
	default:break;
	}
    }
  else
    {
      switch (istk(l2)[0])
	{
	case 0:
	  {
	    verbatim_copy (stk(lL3), row*col*sli, var3, row*col*sli);
	    break;
	  }
	case 1:
	  {
	    dwt3d_tran(stk(lL3), col, row, sli, var3, row, col, sli);
	    break;
	  }
	case 2:
	  {
	    var2 = malloc(m4*n4*sizeof(double));
	    dwt3d_tran(stk(lL3), col, row, sli, var2, row, col, sli);
	    dwt3d_tran_z_inv(var2, row, col, sli, var3, row, sli, col);
	    free(var2);
	    break;
	  }
	case 3:
	  {
	    dwt3d_tran_z_inv(stk(lL3), col, row, sli, var3, col, sli, row);
	    break;
	  }
	case 4:
	  {
	    dwt3d_tran_z(stk(lL3), col, row, sli, var3, col, sli, row);
	    break;
	  }
	case 5:
	  {
	    var2 = malloc(m4*n4*sizeof(double));
	    dwt3d_tran(stk(lL3), col, row, sli, var2, row, col, sli);
	    dwt3d_tran_z(var2, row, col, sli, var3, row, sli, col);
	    free(var2);
	    break;
	  }
	default:break;
	}
    }   
  
 
  CreateVar(4, "m", &mL, &nL, &l4);
  CreateListVarFromPtr(4,1,"S",&ms,&ns,Str);
  CreateListVarFromPtr(4,2,"I",&mi,&ni,&ssi);
  CreateListVarFromPtr(4,3,"d",&m4, &n4, &var3);
  free(var3);
  LhsVar(1) = 4;

  return 0;
}

int
int_wrev3 (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, mL2, nL2, s4;
  static int it3, mL3, nL3, lL3, lcL3, mL1, nL1, s3, l4, m4, n4;
  static int minlhs = 1, maxlhs = 1, minrhs = 2, maxrhs = 2;
  char **Str2;
  static char *Str[]= { "hm","dims","entries"};
  int mL=3, nL=1, ms=1, ns=3, mi=1, ni=3;
  int r, c, s, zero, inc, mn, row, col, sli, r2, i;
  SciIntMat ssi, M;
  int si[3];
  double *temp, *var3, *var2;
  int errCode;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);
  /* Get Input */
  GetRhsVar(1,"m",&m1,&n1,&l1);
  CheckLength(1,m1,3);
  GetListRhsVar(1,1,"S",&mL1,&nL1,&Str2);

  if ( strcmp(Str2[0],"hm") != 0) 
    {
      Scierror(999,"Argument %d is not an hypermatrix\r\n",1);
      return 0;
    }
  FreeRhsSVar(Str2);
  GetListRhsVar(1,2,"I",&mL2,&nL2,&M);
  GetListRhsCVar(1,3,"d",&it3,&mL3,&nL3,&lL3,&lcL3);
  
  GetRhsVar (2, "i", &m2, &n2, &l2);
  if ((istk(l2)[0] != 1) && (istk(l2)[0] != 2) && (istk(l2)[0] != 3) &&
      (istk(l2)[0] != 4) && (istk(l2)[0] != 5) && (istk(l2)[0] != 6) && istk(l2)[0] != 7)
    {
      sciprint("the second argument should be integer from 1 to 7!\n");
      return 0;
    }

  if (it3 == 1)
    {
      Scierror(999,"Argument %d should be real hypermatrix\r\n",1);
      return 0;
    }
  if ((mL2 != 1) || (nL2 != 3))
    {
      Scierror(999,"Argument %d dimension error\r\n",1);
      return 0;
    }
  
  mn = mL2*nL2;
  inc  = 1;
  zero = 0;

  temp = malloc(m1*n1*sizeof(double));
  C2F(tpconv)(&M.it,&zero,&mn, M.D, &inc, temp, &inc);
  row = (int)temp[0];
  col = (int)temp[1];
  sli = (int)temp[2];
  free(temp);
  m3 = row;
  n3 = col;
  s3 = sli;

  si[0] = m3;
  si[1] = n3;
  si[2] = s3;

  ssi.m = 1;
  ssi.n = 3;
  ssi.it = 4;
  ssi.l = 100;
  ssi.D = si;
  m3 = m3*n3*s3;
  n3 = 1;
  var3 = malloc(m3*n3*sizeof(double));

  switch(istk(l2)[0])
    {
    case 1:
      {
	dwt3d_tran(stk(lL3), col, row, sli, var3, row, col, sli);
	var2 = malloc(m3*n3*sizeof(double));
	for(i=0;i<row*sli;i++)
	  wrev(var3+i*col,col,var2+i*col,col);
	dwt3d_tran(var2, row, col, sli, var3, col, row, sli);
	free(var2);
	break;
      }
    case 2:
      {
	for(i=0;i<col*sli;i++)
	  wrev(stk(lL3)+i*row,row,var3+i*row,row);
	break;
      }
    case 3:
      {
	dwt3d_tran_z(stk(lL3), col, row, sli, var3, row, sli, col);
	var2 = malloc(m3*n3*sizeof(double));
	for(i=0;i<row*col;i++)
	  wrev(var3+i*sli,sli,var2+i*sli,sli);
	dwt3d_tran_z_inv(var2, row, sli, col, var3, col, row, sli);
	free(var2);
	break;
      }
    case 4:
      {
	var2 = malloc(m3*n3*sizeof(double));
	dwt3d_tran(stk(lL3), col, row, sli, var2, row, col, sli);
	for(i=0;i<row*sli;i++)
	  wrev(var2+i*col,col,var3+i*col,col);
	dwt3d_tran(var3, row, col, sli, var2, col, row, sli);
	for(i=0;i<col*sli;i++)
	  wrev(var2+i*row,row,var3+i*row,row);
	free(var2);
	break;
      }
    case 5:
      {
	var2 = malloc(m3*n3*sizeof(double));
	dwt3d_tran(stk(lL3), col, row, sli, var2, row, col, sli);
	for(i=0;i<row*sli;i++)
	  wrev(var2+i*col,col,var3+i*col,col);
	dwt3d_tran(var3, row, col, sli, var2, col, row, sli);
	dwt3d_tran_z(var2, col, row, sli, var3, row, sli, col);
	for(i=0;i<row*col;i++)
	  wrev(var3+i*sli,sli,var2+i*sli,sli);
	dwt3d_tran_z_inv(var2, row, sli, col, var3, col, row, sli);
	free(var2);
	break;
      }
    case 6:
      {
	var2 = malloc(m3*n3*sizeof(double));
	for(i=0;i<col*sli;i++)
	  wrev(stk(lL3)+i*row,row,var2+i*row,row);
	dwt3d_tran_z(var2, col, row, sli, var3, row, sli, col);
	for(i=0;i<row*col;i++)
	  wrev(var3+i*sli,sli,var2+i*sli,sli);
	dwt3d_tran_z_inv(var2, row, sli, col, var3, col, row, sli);
	free(var2);
	break;
      }
    case 7:
      {
	var2 = malloc(m3*n3*sizeof(double));
	for(i=0;i<col*sli;i++)
	  wrev(stk(lL3)+i*row,row,var3+i*row,row);
	dwt3d_tran_z(var3, col, row, sli, var2, row, sli, col);
	for(i=0;i<row*col;i++)
	  wrev(var2+i*sli,sli,var3+i*sli,sli);
	dwt3d_tran_z_inv(var3, row, sli, col, var2, col, row, sli);
	dwt3d_tran(var2, col, row, sli, var3, row, col, sli);
	for(i=0;i<row*sli;i++)
	  wrev(var3+i*col,col,var2+i*col,col);
	dwt3d_tran(var2, row, col, sli, var3, col, row, sli);
	free(var2);
	break;
      }
    default:break;
    }

  CreateVar(3, "m", &mL, &nL, &l3);
  CreateListVarFromPtr(3,1,"S",&ms,&ns,Str);
  CreateListVarFromPtr(3,2,"I",&mi,&ni,&ssi);
  CreateListVarFromPtr(3,3,"d",&m3, &n3, &var3);
  free(var3);
  LhsVar(1) = 3;

  return 0;
} 

int
int_wrev2(char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
  static int minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  int i, errCode;
  double *var;

  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);

  wrev2_form_validate (&errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);	
      return 0;			
    }

  GetRhsVar (1, "d", &m1, &n1, &l1);
  GetRhsVar (2, "i", &m2, &n2, &l2);

  if ((istk(l2)[0] != 1) && (istk(l2)[0] != 2) && (istk(l2)[0] != 3))
    {
      sciprint("second argument should be integer from 1 to 3!\n");
      return 0;
    }


  m3 = m1;
  n3 = n1;

  CreateVar (3, "d", &m3, &n3, &l3);

  switch(istk(l2)[0])
    {
    case 1:
      {
	var = malloc(m1*n1*sizeof(double));
	matrix_tran(stk(l1),n1,m1,stk(l3),m1,n1);
	for(i=0;i<m1;i++)
	  wrev(stk(l3)+i*n1,n1,var+i*n1,n1);
	matrix_tran(var,m1,n1,stk(l3),n1,m1);
	free(var);
	break;
      }
    case 2:
      {
	for(i=0;i<n1;i++)
	  wrev(stk(l1)+i*m1,m1,stk(l3)+i*m1,m1);
	break;
      }
    case 3:
      {
	var = malloc(m1*n1*sizeof(double));
	for(i=0;i<n1;i++)
	  wrev(stk(l1)+i*m1,m1,var+i*m1,m1);
	matrix_tran(var,n1,m1,stk(l3),m1,n1);
	for(i=0;i<m1;i++)
	  wrev(stk(l3)+i*n1,n1,var+i*n1,n1);
	matrix_tran(var,m1,n1,stk(l3),n1,m1);
	free(var);
	break;
      }
    default: break;
    }

  LhsVar(1) = 3;

  return 0;
}

int
int_wnorm (char *fname)
{
  static int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
  static int minlhs = 1, maxlhs = 1, minrhs = 1, maxrhs = 3;
  int errCode, flow;
  /* Input Validation */
  CheckRhs (minrhs, maxrhs);
  CheckLhs (minlhs, maxlhs);
  
  /* Get Input */
  GetRhsVar (1, "d", &m1, &n1, &l1);
  flow = 1;
  wnorm_form_validate (&flow, &errCode);
  if (errCode != SUCCESS)
    {
      validate_print (errCode);	/* Indicate error  */
      return 0;			/* Stop when validation fails */
    }

  switch (flow){
  case 1:
    {
      m2 = m1;
      n2 = n1;
      CreateVar (2, "d", &m2, &n2, &l2);
      wcodematd(stk(l1), m1*n1, stk(l2), m2*n2, 0.0, 1.0);
      LhsVar(1) = 2;
      break;
    }
  case 2:
    {
      m4 = m1;
      n4 = n1;
      GetRhsVar (2, "d", &m2, &n2, &l2);
      GetRhsVar (3, "d", &m3, &n3, &l3);
      CreateVar (4, "d", &m4, &n4, &l4);
      if (stk(l2)[0] >= stk(l3)[0])
	{
	  Scierror(999,"min value must be smaller than max value!\n");
	  return 0;
	}
      wcodematd(stk(l1), m1*n1, stk(l4), m4*n4, stk(l2)[0], stk(l3)[0]);
      LhsVar(1) = 4;
      break;
    }
  default:break;
      }
  
  return 0;
}
