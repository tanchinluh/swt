/*
 * -------------------------------------------------------------------------
 * swt_int.c -- SWT interface
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
 #include "swt.h"
 #include "stack-c.h"

 int 
 int_swt(char *fname)
 {
   static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
   static int minlhs=1, maxlhs=2, minrhs=3, maxrhs=4;
   int errCode, flow, family, member, ii, ep,val,stride;
   Func ana_fun;
   swt_wavelet pWaveStruct; 

   CheckRhs (minrhs,maxrhs);
   CheckLhs (minlhs,maxlhs);

   //sprint("after check\n");
   swt_form_validate(&errCode,&flow);
   if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

   l1 = 0;
   l2 = 0;
   l3 = 0;
   l4 = 0;
   
   sciprint("flow=%d\n",flow);

   switch (flow){
   case 1:
	   {
           GetRhsVar(1, "d", &m1, &n1, &l1);
           GetRhsVar(2, "i", &m2, &n2, &l2);
           GetRhsVar(3, "c", &m3, &n3, &l3);
           swt_content_validate(&errCode,flow,l1,l2,l3,14);
           if (errCode != SUCCESS)
	        {
	           validate_print (errCode);
	           return 0;			
	        }
		    swt_exp2(istk(l2)[0],&ep);
			if (((m1*n1)%ep)!=0)
			{
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			}
			
            //sciprint("after content validateion!\n");
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
            //sciprint("after wavelet resolution!\n");
            m4 = istk(l2)[0] + 1;
            n4 = m1 * n1;
			CreateVar(4, "d", &m4, &n4, &l4);
			//sciprint("after creating!\n");
			swt_out1 (stk(l1), m1*n1, stk(l4), m4, n4,
				pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
				pWaveStruct.length, istk(l2)[0]);
            //sciprint("after out!\n");
			filter_clear();
            LhsVar(1) = 4;
		    break;
	   }
   case 2:
	   {
            GetRhsVar(1, "d", &m1, &n1, &l1);
            GetRhsVar(2, "i", &m2, &n2, &l2);
            GetRhsVar(3, "d", &m3, &n3, &l3);
            GetRhsVar(4, "d", &m4, &n4, &l4);
			swt_content_validate(&errCode,flow,l1,l2,l3,14);
            if (errCode != SUCCESS)
	         {
	            validate_print (errCode);
	            return 0;			
	         }
			 swt_exp2(istk(l2)[0],&ep);
			 if (((m1*n1)%ep)!=0)
			 {
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
			 wave_len_validate (m1*n1, m3*n3, &stride, &val);
            if ((!val) || (stride<istk(l2)[0]))
				{
					sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
					return 0;
				}
             m5 = istk(l2)[0] + 1;
             n5 = m1 * n1;
			 CreateVar(5, "d", &m5, &n5, &l5);
			 swt_out1 (stk(l1), m1*n1, stk(l5), m5, n5,
				stk(l3), stk(l4), m3*n3, istk(l2)[0]);
			 //filter_clear();
             LhsVar(1) = 5;
		     break;
	   }
   case 3:
	   {
		     GetRhsVar(1, "d", &m1, &n1, &l1);
             GetRhsVar(2, "i", &m2, &n2, &l2);
             GetRhsVar(3, "c", &m3, &n3, &l3);
             swt_content_validate(&errCode,flow,l1,l2,l3,14);
             if (errCode != SUCCESS)
	          {
	             validate_print (errCode);
	             return 0;			
	          }
			 swt_exp2(istk(l2)[0],&ep);
			 if (((m1*n1)%ep)!=0)
			 {
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
	     //sciprint("after content validateion!\n");
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
		//sciprint("after wavelet resolution!\n");
             
	     m4 = istk(l2)[0];
             n4 = m1 * n1;
			 m5 = m4;
			 n5 = n4;
			 CreateVar(4, "d", &m4, &n4, &l4);
			 CreateVar(5, "d", &m5, &n5, &l5);
		//sciprint("after creating!\n");
             swt_out2 (stk(l1), m1*n1, stk(l4), stk(l5), m4, n4, 
				       pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
				       pWaveStruct.length, istk(l2)[0]);
            //sciprint("after out!\n");
			 filter_clear();
             LhsVar(1) = 4;
             LhsVar(2) = 5;
		     break;
	   }
   case 4:
	   {
		     GetRhsVar(1, "d", &m1, &n1, &l1);
             GetRhsVar(2, "i", &m2, &n2, &l2);
             GetRhsVar(3, "d", &m3, &n3, &l3);
             GetRhsVar(4, "d", &m4, &n4, &l4);
			 swt_content_validate(&errCode,flow,l1,l2,l3,14);
             if (errCode != SUCCESS)
	          {
	            validate_print (errCode);
	            return 0;			
	          }
			 swt_exp2(istk(l2)[0],&ep);
			 if (((m1*n1)%ep)!=0)
			 {
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
			 wave_len_validate (m1*n1, m3*n3, &stride, &val);
             if ((!val) || (stride<istk(l2)[0]))
				{
					sciprint("Input signal is not valid for selected decompostion level and wavelets!\n");
					return 0;
				}
              m5 = istk(l2)[0];
              n5 = m1 * n1;
			  m6 = m5;
			  n6 = n5;
			  CreateVar(5, "d", &m5, &n5, &l5);
			  CreateVar(6, "d", &m6, &n6, &l6);
			  swt_out2 (stk(l1), m1*n1, stk(l5), stk(l6), m5, n5, 
				        stk(l3), stk(l4), m3*n3, istk(l2)[0]);
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
 int_iswt(char *fname)
 {
   static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   static int l4, m4, n4, l5, m5, n5;
   static int minlhs=1, maxlhs=1, minrhs=2, maxrhs=4;
   int flow, errCode, family, member, ii,ep;
   Func syn_fun;
   swt_wavelet pWaveStruct;;

   CheckRhs (minrhs,maxrhs);
   CheckLhs (minlhs,maxlhs);

   iswt_form_validate(&errCode,&flow);
   if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

   l1 = 0;
   l2 = 0;
   l3 = 0;
   l4 = 0;

   switch (flow){
   case 1:
	   {
		   GetRhsVar(1, "d", &m1, &n1, &l1);
           GetRhsVar(2, "c", &m2, &n2, &l2);
		   iswt_content_validate(&errCode,flow,l1,l2,l3,14);
           if (errCode != SUCCESS)
	        {
	           validate_print (errCode);
	           return 0;			
	        }
		   wavelet_parser(cstk(l2),&family,&member);
           wavelet_fun_parser (cstk(l2), &ii);
           syn_fun = wi[ii].synthesis;
          (*syn_fun)(member, &pWaveStruct);
		   swt_exp2(m1-1,&ep);
			 if ((n1%ep)!=0)
			 {
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
		   if (n1 < 2 * pWaveStruct.length)
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
		   m3 = 1;
		   n3 = n1;
           CreateVar(3, "d", &m3, &n3, &l3);
           iswt_input1 (stk(l1), m1, n1, stk(l3), n3, pWaveStruct.pLowPass, 
			            pWaveStruct.pHiPass, pWaveStruct.length);
		   LhsVar(1) = 3;
		   filter_clear();
		   break;
	   }
   case 2:
	   {
		   GetRhsVar(1, "d", &m1, &n1, &l1);
           GetRhsVar(2, "d", &m2, &n2, &l2);
           GetRhsVar(3, "c", &m3, &n3, &l3);
		   iswt_content_validate(&errCode,flow,l1,l2,l3,14);
           if (errCode != SUCCESS)
	        {
	           validate_print (errCode);
	           return 0;			
	        }
		   wavelet_parser(cstk(l3),&family,&member);
           wavelet_fun_parser (cstk(l3), &ii);
           syn_fun = wi[ii].synthesis;
          (*syn_fun)(member, &pWaveStruct);
		  swt_exp2(m1,&ep);
			 if ((n1%ep)!=0)
			 {
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
			 if (n1 < 2 * pWaveStruct.length)
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
		   m4 = 1;
		   n4 = n1;
           CreateVar(4, "d", &m4, &n4, &l4);
		   iswt_input2 (stk(l1), stk(l2), m1, n1, stk(l4), n4, pWaveStruct.pLowPass, 
			            pWaveStruct.pHiPass, pWaveStruct.length);
		   LhsVar(1) = 4;
		   filter_clear();
		   break;
	   }
   case 3:
	   {
		   GetRhsVar(1, "d", &m1, &n1, &l1);
           GetRhsVar(2, "d", &m2, &n2, &l2);
           GetRhsVar(3, "d", &m3, &n3, &l3);
		   iswt_content_validate(&errCode,flow,l1,l2,l3,14);
           if (errCode != SUCCESS)
	        {
	           validate_print (errCode);
	           return 0;			
	        }
		   swt_exp2(m1-1,&ep);
			 if ((n1%ep)!=0)
			 {
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
			 if (n1 < 2 * m2*n2)
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
		   m4 = 1;
		   n4 = n1;
           CreateVar(4, "d", &m4, &n4, &l4);
		   iswt_input1 (stk(l1), m1, n1, stk(l4), n4, stk(l2), 
			            stk(l3), m3*n3);
		   LhsVar(1) = 4;
		   break;
	   }
   case 4:
	   {
		   GetRhsVar(1, "d", &m1, &n1, &l1);
           GetRhsVar(2, "d", &m2, &n2, &l2);
           GetRhsVar(3, "d", &m3, &n3, &l3);
		   GetRhsVar(4, "d", &m4, &n4, &l4);
		   iswt_content_validate(&errCode,flow,l1,l2,l3,14);
           if (errCode != SUCCESS)
	        {
	           validate_print (errCode);
	           return 0;			
	        }
		   swt_exp2(m1,&ep);
			 if ((n1%ep)!=0)
			 {
				sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
			 if (n1 < 2 * m3*n3)
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
		   m5 = 1;
		   n5 = n1;
           CreateVar(5, "d", &m5, &n5, &l5);
		   iswt_input2 (stk(l1), stk(l2), m1, n1, stk(l5), n5, stk(l3), 
			            stk(l4), m3*n3);
		   LhsVar(1) = 5;
		   break;
	   }
   default:
	   break;
   }

   return 0;
 }

 int 
 int_swt2(char *fname)
 {
   static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
   static int l7, m7, n7, l8, m8, n8;
   static int minlhs=1, maxlhs=4, minrhs=3, maxrhs=4;
   int errCode, flow, family, member, ii, ep,val1,stride1;
   int val2, stride2,stride;
   Func ana_fun;
   swt_wavelet pWaveStruct;
   static char *Str[]= { "hm","dims","entries"};
   int mL=3, nL=1, ms=1, ns=3, mi=1, ni=3;
   int si[3];
   double *var4, *var5, *var6, *var7, *var8;
   SciIntMat ssi;

   CheckRhs (minrhs,maxrhs);
   CheckLhs (minlhs,maxlhs);

   swt2_form_validate(&errCode,&flow);
   if (errCode != SUCCESS)
    {
      validate_print (errCode);
      return 0;			
    }

   l1 = 0;
   l2 = 0;
   l3 = 0;
   l4 = 0;
   si[0] = 0;
   si[1] = 0;
   si[2] = 3;

   switch (flow){
	   case 1:
		   {
			   GetRhsVar(1, "d", &m1, &n1, &l1);
               GetRhsVar(2, "i", &m2, &n2, &l2);
               GetRhsVar(3, "c", &m3, &n3, &l3);
               swt2_content_validate(&errCode,flow,l1,l2,l3,14);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }
			   swt_exp2(istk(l2)[0],&ep);
			   if (((m1%ep)!=0) || ((n1%ep)!=0))
			     {
				     sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
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
			   
			   if (istk(l2)[0]==1)
			   {
				   m4 = m1*n1*4;
			       n4 = 1;
				   si[0] = m1;
				   si[1] = n1;
				   si[2] = 4;
				   ssi.m = 1;
				   ssi.n = 3;
				   ssi.l = 100;
				   ssi.it = 4;
				   ssi.D = si;
				   var4 = malloc(m4*n4*sizeof(double));
				   swt2_output1_step(stk(l1), m1, n1, var4, m1, n1,
				        pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
				              pWaveStruct.length, istk(l2)[0]);
				   CreateVar(4, "m", &mL, &nL, &l4);
				   CreateListVarFromPtr(4,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(4,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(4,3,"d",&m4, &n4, &var4);
				   free(var4);
			   }
			   else 
			   {
				   m4 = m1 * n1 * 3 * istk(l2)[0] + m1 * n1;
				   n4 = 1;
				   si[0] = m1;
				   si[1] = n1;
				   si[2] = 3 * istk(l2)[0] + 1;
				   ssi.m = 1;
				   ssi.n = 3;
				   ssi.l = 100;
				   ssi.it = 4;
				   ssi.D = si;
				   var4 = malloc(m4*n4*sizeof(double));
				   swt2_output1_step(stk(l1), m1, n1, var4, m1, n1,
				        pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
				              pWaveStruct.length, istk(l2)[0]);
				   CreateVar(4, "m", &mL, &nL, &l4);
				   CreateListVarFromPtr(4,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(4,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(4,3,"d",&m4, &n4, &var4);
				   free(var4);
			   }
			   filter_clear();
			   LhsVar(1) = 4;
			   break;
		   }
	   case 2:
		   {
			   GetRhsVar(1, "d", &m1, &n1, &l1);
               GetRhsVar(2, "i", &m2, &n2, &l2);
               GetRhsVar(3, "d", &m3, &n3, &l3);
               GetRhsVar(4, "d", &m4, &n4, &l4);
			   swt_content_validate(&errCode,flow,l1,l2,l3,14);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }
			   swt_exp2(istk(l2)[0],&ep);
			   if (((m1%ep)!=0) || ((n1%ep)!=0))
			     {
				    sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
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
			   
			   if (istk(l2)[0]==1)
			   {
				   m5 = m1 * n1 * 4;
			       n5 = 1;
				   si[0] = m1;
				   si[1] = n1;
				   si[2] = 4;
				   ssi.m = 1;
				   ssi.n = 3;
				   ssi.l = 100;
				   ssi.it = 4;
				   ssi.D = si;
				   var5 = malloc(m5*n5*sizeof(double));
				   swt2_output1_step(stk(l1), m1, n1, var5, m1, n1,
				        stk(l3), stk(l4), m3*n3, istk(l2)[0]);
				   CreateVar(5, "m", &mL, &nL, &l5);
				   CreateListVarFromPtr(5,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(5,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(5,3,"d",&m5, &n5, &var5);
				   free(var5);
				   
			   }
			   else
			   {
				   m5 = m1 * n1 * 3 * istk(l2)[0] + m1 * n1;
				   n5 = 1;
				   si[0] = m1;
				   si[1] = n1;
				   si[2] = 3 * istk(l2)[0] + 1;
				   ssi.m = 1;
				   ssi.n = 3;
				   ssi.l = 100;
				   ssi.it = 4;
				   ssi.D = si;
				   var5 = malloc(m5*n5*sizeof(double));
				   swt2_output1_step(stk(l1), m1, n1, var5, m1, n1,
				        stk(l3), stk(l4), m3*n3, istk(l2)[0]);
				   CreateVar(5, "m", &mL, &nL, &l5);
				   CreateListVarFromPtr(5,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(5,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(5,3,"d",&m5, &n5, &var5);
				   free(var5);
			   }
			   LhsVar(1) = 5;
			   
			   break;
		   }
	   case 3:
		   {
			   GetRhsVar(1, "d", &m1, &n1, &l1);
               GetRhsVar(2, "i", &m2, &n2, &l2);
               GetRhsVar(3, "c", &m3, &n3, &l3);
               swt_content_validate(&errCode,flow,l1,l2,l3,14);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }
			   swt_exp2(istk(l2)[0],&ep);
			   if (((m1%ep)!=0) || ((n1%ep)!=0))
			     {
				    sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
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
			   
			   if (istk(l2)[0]==1)
			   {
				   m4 = m1;
			       m5 = m1;
			       m6 = m1;
			       m7 = m1;
			       n4 = n1;
			       n5 = n1;
			       n6 = n1;
			       n7 = n1;
				   CreateVar(4, "d", &m4, &n4, &l4);
				   CreateVar(5, "d", &m5, &n5, &l5);
                   CreateVar(6, "d", &m6, &n6, &l6);
                   CreateVar(7, "d", &m7, &n7, &l7);
				   swt2_output4_step(stk(l1), m1, n1,
	                          stk(l4), stk(l5), stk(l6), stk(l7),
				              m1, n1, pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
				              pWaveStruct.length, istk(l2)[0]);
			   }
			   else
			   {
				   m4 = m1 * n1 * istk(l2)[0];
				   n4 = 1;
				   m5 = m4;
				   m6 = m4;
				   m7 = m4;
				   n5 = n4;
				   n6 = n4;
				   n7 = n4;
				   si[0] = m1;
				   si[1] = n1;
				   si[2] = istk(l2)[0];
				   ssi.m = 1;
				   ssi.n = 3;
				   ssi.l = 100;
				   ssi.it = 4;
				   ssi.D = si;
				   var4 = malloc(m4*n4*sizeof(double));
				   var5 = malloc(m5*n5*sizeof(double));
				   var6 = malloc(m6*n6*sizeof(double));
				   var7 = malloc(m7*n7*sizeof(double));
				   swt2_output4_step(stk(l1), m1, n1,
	                          var4, var5, var6, var7,
				              m1, n1, pWaveStruct.pLowPass, pWaveStruct.pHiPass, 
				              pWaveStruct.length, istk(l2)[0]);


				   CreateVar(4, "m", &mL, &nL, &l4);
				   CreateListVarFromPtr(4,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(4,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(4,3,"d",&m4, &n4, &var4);
				   CreateVar(5, "m", &mL, &nL, &l5);
				   CreateListVarFromPtr(5,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(5,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(5,3,"d",&m5, &n5, &var5);
				   CreateVar(6, "m", &mL, &nL, &l6);
				   CreateListVarFromPtr(6,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(6,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(6,3,"d",&m5, &n5, &var6);
				   CreateVar(7, "m", &mL, &nL, &l7);
				   CreateListVarFromPtr(7,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(7,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(7,3,"d",&m7, &n7, &var7);
				   free(var4);
				   free(var5);
				   free(var6);
				   free(var7);
			   }
			   filter_clear();
			   LhsVar(1) = 4;
			   LhsVar(2) = 5;
			   LhsVar(3) = 6;
			   LhsVar(4) = 7;
			   break;
		   }
	   case 4:
		   {
			   GetRhsVar(1, "d", &m1, &n1, &l1);
               GetRhsVar(2, "i", &m2, &n2, &l2);
               GetRhsVar(3, "d", &m3, &n3, &l3);
               GetRhsVar(4, "d", &m4, &n4, &l4);
			   swt_content_validate(&errCode,flow,l1,l2,l3,14);
               if (errCode != SUCCESS)
	            {
	              validate_print (errCode);
	              return 0;			
	            }
			   swt_exp2(istk(l2)[0],&ep);
			   if (((m1%ep)!=0) || ((n1%ep)!=0))
			     {
				    sciprint("Input length should be multiples of power of 2! Please extend the input!\n");
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
			   
			   if (istk(l2)[0]==1)
			   {
				   m5 = m1;
			       m6 = m1;
			       m7 = m1;
			       m8 = m1;
			       n5 = n1;
			       n6 = n1;
			       n7 = n1;
			       n8 = n1;
				   CreateVar(5, "d", &m5, &n5, &l5);
				   CreateVar(6, "d", &m6, &n6, &l6);
                   CreateVar(7, "d", &m7, &n7, &l7);
                   CreateVar(8, "d", &m8, &n8, &l8);
				   swt2_output4_step(stk(l1), m1, n1,
	                          stk(l5), stk(l6), stk(l7), stk(l8),
				              m1, n1, stk(l3), stk(l4), 
				              m3*n3, istk(l2)[0]);
			   }
			   else
			   {
				   m5 = m1 * n1 * istk(l2)[0];
				   n5 = 1;
				   m6 = m5;
				   m7 = m5;
				   m8 = m5;
				   n6 = n5;
				   n7 = n5;
				   n8 = n5;
				   si[0] = m1;
				   si[1] = n1;
				   si[2] = istk(l2)[0];
				   ssi.m = 1;
				   ssi.n = 3;
				   ssi.l = 100;
				   ssi.it = 4;
				   ssi.D = si;
				   var5 = malloc(m5*n5*sizeof(double));
				   var6 = malloc(m6*n6*sizeof(double));
				   var7 = malloc(m7*n7*sizeof(double));
				   var8 = malloc(m8*n8*sizeof(double));

				   swt2_output4_step(stk(l1), m1, n1,
	                          var5, var6, var7, var8,
				              m1, n1, stk(l3), stk(l4), 
				              m3*n3, istk(l2)[0]);
				   CreateVar(5, "m", &mL, &nL, &l5);
				   CreateListVarFromPtr(5,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(5,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(5,3,"d",&m5, &n5, &var5);
				   CreateVar(6, "m", &mL, &nL, &l6);
				   CreateListVarFromPtr(6,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(6,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(6,3,"d",&m5, &n5, &var6);
				   CreateVar(7, "m", &mL, &nL, &l7);
				   CreateListVarFromPtr(7,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(7,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(7,3,"d",&m7, &n7, &var7);
				   CreateVar(8, "m", &mL, &nL, &l8);
				   CreateListVarFromPtr(8,1,"S",&ms,&ns,Str);
				   CreateListVarFromPtr(8,2,"I",&mi,&ni,&ssi);
				   CreateListVarFromPtr(8,3,"d",&m8, &n8, &var8);
				   free(var5);
				   free(var6);
				   free(var7);
				   free(var8);

			   }
			   LhsVar(1) = 5;
			   LhsVar(2) = 6;
			   LhsVar(3) = 7;
			   LhsVar(4) = 8;
			   break;
		   }
	   default:
		   break;
   }

   return 0;
 }


 int 
 int_iswt2(char *fname)
 {
   static int l1, m1, n1, l2, m2, n2, l3, m3, n3;
   static int l4, m4, n4, l5, m5, n5, l6, m6, n6;
   static int l7, m7, n7, l, m, n;
   static int minlhs=1, maxlhs=1, minrhs=2, maxrhs=6;
   int flow, errCode, family, member, ii,ep, step;
   Func syn_fun;
   swt_wavelet pWaveStruct;
   char **Str2;
   int mL1, nL1, mL2, nL2, mL3, nL3, lL3, it3;
   SciIntMat lL2;
   int lcL3, lLA3, lLH3, lLV3, lLD3;
   int m2A, m2H, m2V, m2D, n2A, n2H, n2V, n2D, stepA, stepH, stepV, stepD;
   //HyperMat *H;

   CheckRhs (minrhs,maxrhs);
   CheckLhs (minlhs,maxlhs);

   iswt2_form_validate(&errCode,&flow);
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
   

   switch (flow)
   {
   case 1:
	   {
		   GetRhsVar(1,"m",&m1,&n1,&l1);
           CheckLength(1,m1,3);
		   GetListRhsVar(1,1,"S",&mL1,&nL1,&Str2);
           if ( strcmp(Str2[0],"hm") != 0) 
              {
                 Scierror(999,"Argument %d is not an hypermatrix\r\n",1);
                 return 0;
              }
		   FreeRhsSVar(Str2);
           GetListRhsVar(1,2,"I",&mL2,&nL2,&lL2);
           GetListRhsCVar(1,3,"d",&it3,&mL3,&nL3,&lL3,&lcL3);
           
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
           GetRhsVar(2,"c",&m2,&n2,&l2);

           iswt2_content_validate(&errCode,flow,l1,l2,l3,14,l5,l6);
           if (errCode != SUCCESS)
	        {
	           validate_print (errCode);
	           return 0;			
	        }
		   wavelet_parser(cstk(l2),&family,&member);
           wavelet_fun_parser (cstk(l2), &ii);
           syn_fun = wi[ii].synthesis;
          (*syn_fun)(member, &pWaveStruct);
		   

		   m3 = istk(lL2.l)[0];
		   n3 = istk(lL2.l)[1];
		   step = (istk(lL2.l)[2] - 1)/3;
		   if (step<1)
				   {
				      sciprint("Inputs are not swt2 results!\n");
				      return 0;
			       }
		   swt_exp2(step,&ep);
			 if (((m3%ep)!=0) || ((n3%ep)!=0))
			 {
				sciprint("Input size should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
		   if ((m3 < pWaveStruct.length * 2) || (n3 < pWaveStruct.length * 2))
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
		   
		   CreateVar(3, "d", &m3, &n3, &l3);
		   iswt2_input1_step(stk(lL3),  m3, n3, stk(l3), m3, n3,
		               pWaveStruct.pLowPass, pWaveStruct.pHiPass, pWaveStruct.length, step);

		   LhsVar(1) = 3;
		   filter_clear();
           break;
	   }
   case 2:
	   {
		   GetRhsVar(1,"m",&m1,&n1,&l1);
           CheckLength(1,m1,3);
		   GetListRhsVar(1,1,"S",&mL1,&nL1,&Str2);
           if ( strcmp(Str2[0],"hm") != 0) 
              {
                 Scierror(999,"Argument %d is not an hypermatrix\r\n",1);
                 return 0;
              }
		   FreeRhsSVar(Str2);
           GetListRhsVar(1,2,"I",&mL2,&nL2,&lL2);
           GetListRhsCVar(1,3,"d",&it3,&mL3,&nL3,&lL3,&lcL3);
           
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
           GetRhsVar(2,"d",&m2,&n2,&l2);
		   GetRhsVar(3,"d",&m3,&n3,&l3);
		   iswt2_content_validate(&errCode,flow,l1,l2,l3,14,l5,l6);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }

		   m4 = istk(lL2.l)[0];
		   n4 = istk(lL2.l)[1];
		   step = (istk(lL2.l)[2] - 1)/3;
		   if (step<1)
				   {
				      sciprint("Inputs are not swt2 results!\n");
				      return 0;
			       }
		   swt_exp2(step,&ep);
			 if (((m4%ep)!=0) || ((n4%ep)!=0))
			 {
				sciprint("Input size should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
		   if ((m4 < m3 * n3 * 2) || (n4 < m3 * n3 * 2))
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
		   
		   CreateVar(4, "d", &m4, &n4, &l4);
		   iswt2_input1_step(stk(lL3),  m4, n4, stk(l4), m4, n4,
		               stk(l2), stk(l3), m2*n2, step);

		   LhsVar(1) = 4;

		   break;
	   }
   case 3:
	   {
		   if (GetType(1)==sci_matrix)
		   {
			   GetRhsVar(1, "d", &m1, &n1, &l1);
               GetRhsVar(2, "d", &m2, &n2, &l2);
               GetRhsVar(3, "d", &m3, &n3, &l3);
               GetRhsVar(4, "d", &m4, &n4, &l4);
			   GetRhsVar(5, "c", &m5, &n5, &l5);
			   iswt2_content_validate(&errCode,flow,l1,l2,l3,14,l5,l6);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }
		       wavelet_parser(cstk(l5),&family,&member);
               wavelet_fun_parser (cstk(l5), &ii);
               syn_fun = wi[ii].synthesis;
              (*syn_fun)(member, &pWaveStruct);

			   m6 = m1;
			   n6 = n1;
			   step = 1;

			     if (step<1)
				   {
				      sciprint("Inputs are not swt2 results!\n");
				      return 0;
			       }
			   swt_exp2(step,&ep);
			 if (((m6%ep)!=0) || ((n6%ep)!=0))
			 {
				sciprint("Input size should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
		   if ((m6 < pWaveStruct.length * 2) || (n6 < pWaveStruct.length * 2))
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
			   CreateVar(6, "d", &m6, &n6, &l6);
			   iswt2_input4_step(stk(l1), stk(l2), stk(l3), stk(l4), m1, n1,
		                         stk(l6), m6, n6, pWaveStruct.pLowPass, 
								 pWaveStruct.pHiPass, pWaveStruct.length, step);

		   }
		   else
		   {
			   /* A */
			   GetRhsVar(1,"m",&m1,&n1,&l1);
               CheckLength(1,m1,3);
		       GetListRhsVar(1,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",1);
                     return 0;
                  }
               GetListRhsVar(1,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(1,3,"d",&it3,&mL3,&nL3,&lLA3,&lcL3);
           
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
			   m2A = istk(lL2.l)[0];
			   n2A = istk(lL2.l)[1];
			   stepA = istk(lL2.l)[2];

			   /* H */
			   GetRhsVar(2,"m",&m2,&n2,&l2);
               CheckLength(2,m2,3);
		       GetListRhsVar(2,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",2);
                     return 0;
                  }
               GetListRhsVar(2,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(2,3,"d",&it3,&mL3,&nL3,&lLH3,&lcL3);
           
		       if (it3 == 1)
			      {
                     Scierror(999,"Argument %d should be real hypermatrix\r\n",2);
                     return 0;
                  }
		       if ((mL2 != 1) || (nL2 != 3))
                  {
                     Scierror(999,"Argument %d dimension error\r\n",2);
                     return 0;
                  }
			   m2H = istk(lL2.l)[0];
			   n2H = istk(lL2.l)[1];
			   stepH = istk(lL2.l)[2];

			   /* V */
			   GetRhsVar(3,"m",&m3,&n3,&l3);
               CheckLength(3,m3,3);
		       GetListRhsVar(3,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",3);
                     return 0;
                  }
               GetListRhsVar(3,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(3,3,"d",&it3,&mL3,&nL3,&lLV3,&lcL3);
           
		       if (it3 == 1)
			      {
                     Scierror(999,"Argument %d should be real hypermatrix\r\n",3);
                     return 0;
                  }
		       if ((mL2 != 1) || (nL2 != 3))
                  {
                     Scierror(999,"Argument %d dimension error\r\n",3);
                     return 0;
                  }
			   m2V = istk(lL2.l)[0];
			   n2V = istk(lL2.l)[1];
			   stepV = istk(lL2.l)[2];

			   /* D */
			   GetRhsVar(4,"m",&m4,&n4,&l4);
               CheckLength(4,m4,3);
		       GetListRhsVar(4,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",4);
                     return 0;
                  }
               GetListRhsVar(4,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(4,3,"d",&it3,&mL3,&nL3,&lLD3,&lcL3);
           
		       if (it3 == 1)
			      {
                     Scierror(999,"Argument %d should be real hypermatrix\r\n",4);
                     return 0;
                  }
		       if ((mL2 != 1) || (nL2 != 3))
                  {
                     Scierror(999,"Argument %d dimension error\r\n",4);
                     return 0;
                  }
			   m2D = istk(lL2.l)[0];
			   n2D = istk(lL2.l)[1];
			   stepD = istk(lL2.l)[2];
			   FreeRhsSVar(Str2);

			   if ((m2A != m2H) || (m2A != m2V) || (m2A != m2D)) 
			   {
				   sciprint("Inputs are not swt2 results!\n");
				   return 0;
			   }
			   if ((n2A != n2H) || (n2A != n2V) || (n2A != n2D)) 
			   {
				   sciprint("Inputs are not swt2 results!\n");
				   return 0;
			   }
			   if ((stepA != stepH) || (stepA != stepV) || (stepA != stepD)) 
			   {
				   sciprint("Inputs are not swt2 results!\n");
				   return 0;
			   }

			   GetRhsVar(5, "c", &m5, &n5, &l5);
			   iswt2_content_validate(&errCode,flow,l1,l2,l3,14,l5,l6);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }
		       wavelet_parser(cstk(l5),&family,&member);
               wavelet_fun_parser (cstk(l5), &ii);
               syn_fun = wi[ii].synthesis;
              (*syn_fun)(member, &pWaveStruct);


			   m6 = istk(lL2.l)[0];
			   n6 = istk(lL2.l)[1];
			   step = istk(lL2.l)[2];

			     if (step<1)
				   {
				      sciprint("Inputs are not swt2 results!\n");
				      return 0;
			       }
			   swt_exp2(step,&ep);
			 if (((m6%ep)!=0) || ((n6%ep)!=0))
			 {
				sciprint("Input size should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
		   if ((m6 < pWaveStruct.length * 2) || (n6 < pWaveStruct.length * 2))
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
			   CreateVar(6, "d", &m6, &n6, &l6);
			   iswt2_input4_step(stk(lLA3), stk(lLH3), stk(lLV3), stk(lLD3), m6, n6,
		                         stk(l6), m6, n6, pWaveStruct.pLowPass, 
								 pWaveStruct.pHiPass, pWaveStruct.length, step);
		       

		   }

		   LhsVar(1) = 6;
		   filter_clear();
		   break;
	   }
   case 4:
	   {
		   if (GetType(1)==sci_matrix)
		   {
			   GetRhsVar(1, "d", &m1, &n1, &l1);
               GetRhsVar(2, "d", &m2, &n2, &l2);
               GetRhsVar(3, "d", &m3, &n3, &l3);
               GetRhsVar(4, "d", &m4, &n4, &l4);
			   GetRhsVar(5, "d", &m5, &n5, &l5);
			   GetRhsVar(6, "d", &m6, &n6, &l6);
			   iswt2_content_validate(&errCode,flow,l1,l2,l3,14,l5,l6);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }


			   m7 = m1;
			   n7 = n1;
			   step = 1;
			   swt_exp2(step,&ep);
			 if (((m7%ep)!=0) || ((n7%ep)!=0))
			 {
				sciprint("Input size should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
		   if ((m7 < m5 * n5 * 2) || (n7 < m5 * n5 * 2))
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
			   CreateVar(7, "d", &m7, &n7, &l7);
			   iswt2_input4_step(stk(l1), stk(l2), stk(l3), stk(l4), m7, n7,
		                         stk(l7), m7, n7, stk(l5), stk(l6), m5*n5, step);
			   


		   }
		   else
		   {
			   /* A */
			   GetRhsVar(1,"m",&m1,&n1,&l1);
               CheckLength(1,m1,3);
		       GetListRhsVar(1,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",1);
                     return 0;
                  }
               GetListRhsVar(1,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(1,3,"d",&it3,&mL3,&nL3,&lLA3,&lcL3);
           
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
			   m2A = istk(lL2.l)[0];
			   n2A = istk(lL2.l)[1];
			   stepA = istk(lL2.l)[2];

			   /* H */
			   GetRhsVar(2,"m",&m2,&n2,&l2);
               CheckLength(2,m2,3);
		       GetListRhsVar(2,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",2);
                     return 0;
                  }
               GetListRhsVar(2,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(2,3,"d",&it3,&mL3,&nL3,&lLH3,&lcL3);
           
		       if (it3 == 1)
			      {
                     Scierror(999,"Argument %d should be real hypermatrix\r\n",2);
                     return 0;
                  }
		       if ((mL2 != 1) || (nL2 != 3))
                  {
                     Scierror(999,"Argument %d dimension error\r\n",2);
                     return 0;
                  }
			   m2H = istk(lL2.l)[0];
			   n2H = istk(lL2.l)[1];
			   stepH = istk(lL2.l)[2];

			   /* V */
			   GetRhsVar(3,"m",&m3,&n3,&l3);
               CheckLength(3,m3,3);
		       GetListRhsVar(3,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",3);
                     return 0;
                  }
               GetListRhsVar(3,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(3,3,"d",&it3,&mL3,&nL3,&lLV3,&lcL3);
           
		       if (it3 == 1)
			      {
                     Scierror(999,"Argument %d should be real hypermatrix\r\n",3);
                     return 0;
                  }
		       if ((mL2 != 1) || (nL2 != 3))
                  {
                     Scierror(999,"Argument %d dimension error\r\n",3);
                     return 0;
                  }
			   m2V = istk(lL2.l)[0];
			   n2V = istk(lL2.l)[1];
			   stepV = istk(lL2.l)[2];

			   /* D */
			   GetRhsVar(4,"m",&m4,&n4,&l4);
               CheckLength(4,m4,3);
		       GetListRhsVar(4,1,"S",&mL1,&nL1,&Str2);
               if ( strcmp(Str2[0],"hm") != 0) 
                  {
                     Scierror(999,"Argument %d is not an hypermatrix\r\n",4);
                     return 0;
                  }
               GetListRhsVar(4,2,"I",&mL2,&nL2,&lL2);
               GetListRhsCVar(4,3,"d",&it3,&mL3,&nL3,&lLD3,&lcL3);
           
		       if (it3 == 1)
			      {
                     Scierror(999,"Argument %d should be real hypermatrix\r\n",4);
                     return 0;
                  }
		       if ((mL2 != 1) || (nL2 != 3))
                  {
                     Scierror(999,"Argument %d dimension error\r\n",4);
                     return 0;
                  }
			   m2D = istk(lL2.l)[0];
			   n2D = istk(lL2.l)[1];
			   stepD = istk(lL2.l)[2];
			   FreeRhsSVar(Str2);

			   if ((m2A != m2H) || (m2A != m2V) || (m2A != m2D)) 
			   {
				   sciprint("Inputs are not swt2 results!\n");
				   return 0;
			   }
			   if ((n2A != n2H) || (n2A != n2V) || (n2A != n2D)) 
			   {
				   sciprint("Inputs are not swt2 results!\n");
				   return 0;
			   }
			   if ((stepA != stepH) || (stepA != stepV) || (stepA != stepD)) 
			   {
				   sciprint("Inputs are not swt2 results!\n");
				   return 0;
			   }

			   GetRhsVar(5, "d", &m5, &n5, &l5);
               GetRhsVar(6, "d", &m6, &n6, &l6);
			   iswt2_content_validate(&errCode,flow,l1,l2,l3,14,l5,l6);
               if (errCode != SUCCESS)
	            {
	               validate_print (errCode);
	               return 0;			
	            }


			   m7 = istk(lL2.l)[0];
			   n7 = istk(lL2.l)[1];
			   step = istk(lL2.l)[2];
			   swt_exp2(step,&ep);
			 if (((m7%ep)!=0) || ((n7%ep)!=0))
			 {
				sciprint("Input size should be multiples of power of 2! Please extend the input!\n");
				return 0;
			 }
		   if ((m7 < m5 * n5 * 2) || (n7 < m5 * n5 * 2))
	         {
	            sciprint ("Input signal is not valid for selected decompostion level and wavelets!\n");
	            return 0;
	         }
			   CreateVar(7, "d", &m7, &n7, &l7);
			   iswt2_input4_step(stk(lLA3), stk(lLH3), stk(lLV3), stk(lLD3), m7, n7,
		                         stk(l7), m7, n7, stk(l5), stk(l6), m5*n5, step);


		   }
		   
		   LhsVar(1) = 7;
		   break;
	   }
   default:
	   break;
   }
   return 0;
 }
