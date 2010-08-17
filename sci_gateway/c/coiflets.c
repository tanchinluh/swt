/*
 * -------------------------------------------------------------------------
 * coiflet.c -- Coiflet wavelets coefficients.
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

/*********************************************
 * Local Variable (Filter Coefficent)
 ********************************************/

static const double coif1[6] = { 
	             -0.072732619512853897, 0.33789766245780922, 
				 0.85257202021225542, 0.38486484686420286, 
				 -0.072732619512853897, -0.01565572813546454};
	

static const double coif2[12] = { 
	             0.016387336463522112, -0.041464936781759151, 
				 -0.067372554721963018, 0.38611006682116222, 
				 0.81272363544554227, 0.41700518442169254, 
				 -0.076488599078306393, -0.059434418646456898, 
				 0.023680171946334084, 0.0056114348193944995, 
				 -0.0018232088707029932, -0.00072054944536451221};
	
static const double coif3[18] = { 
	            -0.0037935128644910141, 0.0077825964273254182, 
				0.023452696141836267, -0.0657719112818555, 
				-0.061123390002672869, 0.4051769024096169, 
				0.79377722262562056, 0.42848347637761874, 
				-0.071799821619312018, -0.082301927106885983, 
				0.034555027573061628, 0.015880544863615904, 
				-0.0090079761366615805, -0.0025745176887502236, 
				0.0011175187708906016, 0.00046621696011288631, 
				-7.0983303138141252e-005, -3.4599772836212559e-005};
	

static const double coif4[24] = { 
	            0.00089231366858231456, -0.0016294920126017326, 
				-0.0073461663276420935, 0.016068943964776348, 
				0.026682300156053072, -0.081266699680878754, 
				-0.056077313316754807, 0.41530840703043026, 
				0.78223893092049901, 0.4343860564914685, 
				-0.066627474263425038, -0.096220442033987982, 
				0.039334427123337491, 0.025082261844864097, 
				-0.015211731527946259, -0.0056582866866107199, 
				0.0037514361572784571, 0.0012665619292989445, 
				-0.00058902075624433831, -0.00025997455248771324, 
				6.2339034461007128e-005, 3.1229875865345646e-005, 
				-3.2596802368833675e-006, -1.7849850030882614e-006};
	
	

static const double coif5[30] = { 
	            -0.00021208083980379827, 0.00035858968789573785, 
				0.0021782363581090178, -0.004159358781386048, 
				-0.010131117519849788, 0.023408156785839195, 
				0.02816802897093635, -0.091920010559696244, 
				-0.052043163176243773, 0.42156620669085149, 
				0.77428960365295618, 0.43799162617183712, 
				-0.062035963962903569, -0.10557420870333893, 
				0.041289208750181702, 0.032683574267111833, 
				-0.019761778942572639, -0.0091642311624818458, 
				0.0067641854480530832, 0.0024333732126576722, 
				-0.0016628637020130838, -0.00063813134304511142, 
				0.00030225958181306315, 0.00014054114970203437, 
				-4.1340432272512511e-005, -2.1315026809955787e-005, 
				3.7346551751414047e-006, 2.0637618513646814e-006, 
				-1.6744288576823017e-007, -9.517657273819165e-008};
	
/*********************************************
 * Global Function
 ********************************************/

void
coiflets_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{

  double *pFilterCoef;

  pWaveStruct->length = 6 * member;

  switch (member)
    {
    case 1:
      pFilterCoef = coif1;
      break;
    case 2:
      pFilterCoef = coif2;
      break;
    case 3:
      pFilterCoef = coif3;
      break;
    case 4:
      pFilterCoef = coif4;
      break;
    case 5:
      pFilterCoef = coif5;
      break;
    default:
      printf("db%d is not available!\n",member);
      exit(0);
    }

  wrev(pFilterCoef, pWaveStruct->length, 
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(pFilterCoef, pWaveStruct->length,
	   HiDecomFilCoef, pWaveStruct->length);
  //for (count = 0; count < pWaveStruct->length; count++)
    //LowDecomFilCoef[count] *= sqrt(2.0);
  //for (count = 0; count < pWaveStruct->length; count++)
    //HiDecomFilCoef[count] *= sqrt(2.0);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;

  return;
}


void
coiflets_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{

  double *pFilterCoef;

  pWaveStruct->length = 6 * member;

  switch (member)
    {
    case 1:
      pFilterCoef = coif1;
      break;
    case 2:
      pFilterCoef = coif2;
      break;
    case 3:
      pFilterCoef = coif3;
      break;
    case 4:
      pFilterCoef = coif4;
      break;
    case 5:
      pFilterCoef = coif5;
      break;
    default:
      printf("db%d is not available!\n",member);
      exit(0);
    }

  verbatim_copy(pFilterCoef, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(pFilterCoef, pWaveStruct->length,
      HiReconFilCoef, pWaveStruct->length );
  //for (count = 0; count < pWaveStruct->length; count++)
    //LowReconFilCoef[count] *= sqrt(2.0);
  //for (count = 0; count < pWaveStruct->length; count++)
    //HiReconFilCoef[count] *= sqrt(2.0);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;

  return;
}
