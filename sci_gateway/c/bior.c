/*
 * -------------------------------------------------------------------------
 * bspline.c -- Biorthogonal wavelets coefficents.
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

static const double h1[10] = {0.0, 0.0, 0.0, 0.0, 
			      0.70710678118654752440084436210, 
			      0.70710678118654752440084436210, 
			      0.0, 0.0, 0.0, 0.0};

static const double hm1_11[2] = {0.70710678118654752440084436210,
				0.70710678118654752440084436210};

static const double hm1_13[6] = {-0.0883883476483184405501055452631,
				0.0883883476483184405501055452631,
				0.70710678118654752440084436210, 
				0.70710678118654752440084436210, 
				0.0883883476483184405501055452631,
				-0.0883883476483184405501055452631};

static const double hm1_15[10] = {0.0165728151840597076031447897368,
				 -0.0165728151840597076031447897368,
				 -0.1215339780164378557563951247368,
				 0.1215339780164378557563951247368,
				 0.70710678118654752440084436210, 
				 0.70710678118654752440084436210, 
				 0.1215339780164378557563951247368,
				 -0.1215339780164378557563951247368,
				 -0.0165728151840597076031447897368,
				 0.0165728151840597076031447897368};


static const double h2[18] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			      0.3535533905932737622004221810524,
			      0.7071067811865475244008443621048,
			      0.3535533905932737622004221810524,
			      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
			      0.0};

static const double hm2_22[6] = {-0.1767766952966368811002110905262,
				 0.3535533905932737622004221810524,
				 1.0606601717798212866012665431573,
				 0.3535533905932737622004221810524,
				 -0.1767766952966368811002110905262,
				 0.0};

static const double hm2_24[10] = {0.0331456303681194152062895794737,
				  -0.0662912607362388304125791589473,
				  -0.1767766952966368811002110905262,
				  0.4198446513295125926130013399998,
				  0.9943689110435824561886873842099,
				  0.4198446513295125926130013399998,
				  -0.1767766952966368811002110905262,
				  -0.0662912607362388304125791589473,
				  0.0331456303681194152062895794737,
				  0.0};

static const double hm2_26[14] = {-0.0069053396600248781679769957237,
				  0.0138106793200497563359539914474,
				  0.0469563096881691715422435709210,
				  -0.1077232986963880994204411332894,
				  -0.1698713556366120029322340948025,
				  0.4474660099696121052849093228945,
				  0.9667475524034829435167794013152,
				  0.4474660099696121052849093228945,
				  -0.1698713556366120029322340948025,
				  -0.1077232986963880994204411332894,
				  0.0469563096881691715422435709210,
				  0.0138106793200497563359539914474,
				  -0.0069053396600248781679769957237,
				  0.0};

static const double hm2_28[18] = {0.0015105430506304420992449678146,
				  -0.0030210861012608841984899356291,
				  -0.0129475118625466465649568669819,
				  0.0289161098263541773284036695929,
				  0.0529984818906909399392234421792,
				  -0.1349130736077360572068505539514,
				  -0.1638291834340902345352542235443,
				  0.4625714404759165262773590010400,
				  0.9516421218971785225243297231697,
				  0.4625714404759165262773590010400,
				  -0.1638291834340902345352542235443,
				  -0.1349130736077360572068505539514,
				  0.0529984818906909399392234421792,
				  0.0289161098263541773284036695929,
				  -0.0129475118625466465649568669819,
				  -0.0030210861012608841984899356291,
				  0.0015105430506304420992449678146,
				  0.0};

static const double h3[20] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			      0.1767766952966368811002110905262,
			      0.5303300858899106433006332715786,
			      0.5303300858899106433006332715786,
			      0.1767766952966368811002110905262,
			      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
static const double hm3_31[4] = {-0.3535533905932737622004221810524,
				 1.0606601717798212866012665431573,
				 1.0606601717798212866012665431573,
				 -0.3535533905932737622004221810524};

static const double hm3_33[8] = {0.0662912607362388304125791589473,
				 -0.1988737822087164912377374768420,
				 -0.1546796083845572709626847042104,
				 0.9943689110435824561886873842099,
				 0.9943689110435824561886873842099,
				 -0.1546796083845572709626847042104,
				 -0.1988737822087164912377374768420,
				 0.0662912607362388304125791589473};

static const double hm3_35[12] = {-0.0138106793200497563359539914474,
				  0.0414320379601492690078619743421,
				  0.0524805814161890740766251675000,
				  -0.2679271788089652729175074340788,
				  -0.0718155324642587329469607555263,
				  0.9667475524034829435167794013152,
				  0.9667475524034829435167794013152,
				  -0.0718155324642587329469607555263,
				  -0.2679271788089652729175074340788,
				  0.0524805814161890740766251675000,
				  0.0414320379601492690078619743421,
				  -0.0138106793200497563359539914474};

static const double hm3_37[16] = {0.0030210861012608841984899356291,
				  -0.0090632583037826525954698068873,
				  -0.0168317654213106405344439270765,
				  0.0746639850740189951912512662623,
				  0.0313329787073628846871956180962,
				  -0.3011591259228349991008967259990,
				  -0.0264992409453454699696117210896,
				  0.9516421218971785225243297231697,
				  0.9516421218971785225243297231697,
				  -0.0264992409453454699696117210896,
				  -0.3011591259228349991008967259990,
				  0.0313329787073628846871956180962,
				  0.0746639850740189951912512662623,
				  -0.0168317654213106405344439270765,
				  -0.0090632583037826525954698068873,
				  0.0030210861012608841984899356291};


static const double hm3_39[20] = {-0.0006797443727836989446602355165,
				  0.0020392331183510968339807065496,
				  0.0050603192196119810324706421788,
				  -0.0206189126411055346546938106687,
				  -0.0141127879301758447558029850103,
				  0.0991347824942321571990197448581,
				  0.0123001362694193142367090236328,
				  -0.3201919683607785695513833204624,
				  0.0020500227115698857061181706055,
				  0.9421257006782067372990864259380,
				  0.9421257006782067372990864259380,
				  0.0020500227115698857061181706055,
				  -0.3201919683607785695513833204624,
				  0.0123001362694193142367090236328,
				  0.0991347824942321571990197448581,
				  -0.0141127879301758447558029850103,
				  -0.0206189126411055346546938106687,
				  0.0050603192196119810324706421788,
				  0.0020392331183510968339807065496,
				  -0.0006797443727836989446602355165};


static const double h4[10] = {
                  0.0, -0.064538882628697058, 
				  -0.040689417609164058, 0.41809227322161724, 
				  0.7884856164055829, 0.41809227322161724, 
				  -0.040689417609164058, -0.064538882628697058, 
				  0.0, 0.0};

static const double hm4_44[10] = {
	              0.03782845550726404, -0.023849465019556843,
                  -0.11062440441843718,0.37740285561283066,
                  0.85269867900889385,0.37740285561283066,
                  -0.11062440441843718, -0.023849465019556843,
                  0.03782845550726404, 0.0};

static const double h5[12] = {
                  0.013456709459118716, -0.0026949668801115071, 
				  -0.13670658466432914, -0.093504697400938863, 
				  0.47680326579848425, 0.89950610974864842, 
				  0.47680326579848425, -0.093504697400938863, 
				  -0.13670658466432914, -0.0026949668801115071, 
				  0.013456709459118716, 0.0};

static const double hm5_55[12] = {
	              0.0, 0.03968708834740544,
				  0.0079481086372403219, -0.054463788468236907,
                  0.34560528195603346,0.73666018142821055,
                  0.34560528195603346, -0.054463788468236907,
                  0.0079481086372403219, 0.03968708834740544,
                  0.0, 0.0};

static const double h6[18] = {
                  0.0, 0.0, 
				  0.0, 0.014426282505624435, 
				  0.014467504896790148, -0.078722001062628819, 
				  -0.040367979030339923, 0.41784910915027457, 
				  0.75890772945365415, 0.41784910915027457, 
				  -0.040367979030339923, -0.078722001062628819, 
				  0.014467504896790148, 0.014426282505624435, 
				  0.0, 0.0, 
				  0.0, 0.0};

static const double hm6_68[18] = {
                 0.0019088317364812906, -0.0019142861290887667, 
				 -0.016990639867602342, 0.01193456527972926, 
				 0.04973290349094079, -0.077263173167204144, 
				 -0.09405920349573646, 0.42079628460982682, 
				 0.82592299745840225, 0.42079628460982682, 
				 -0.09405920349573646, -0.077263173167204144, 
				 0.04973290349094079, 0.01193456527972926, 
				 -0.016990639867602342, -0.0019142861290887667, 
				 0.0019088317364812906, 0.0};

/*********************************************
 * Global Function
 ********************************************/

void
sp_bior_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
  double *pFilterCoef, *pFilterCoefMirror;

  switch (member) {
  case 11:
    {
      pWaveStruct->length = 2;
      pFilterCoef = hm1_11;
      pFilterCoefMirror = h1+4;
      break;
    }
  case 13:
    {
      pWaveStruct->length = 6;
      pFilterCoef = hm1_13;
      pFilterCoefMirror = h1+2;
      break;
    }
  case 15:
    {
      pWaveStruct->length = 10;
      pFilterCoef = hm1_15;
      pFilterCoefMirror = h1;
      break;
    }
  case 22:
    {
      pWaveStruct->length = 6;
      pFilterCoef = hm2_22;
      pFilterCoefMirror = h2+6;
      break;
    }
  case 24:
    {
      pWaveStruct->length = 10;
      pFilterCoef = hm2_24;
      pFilterCoefMirror = h2+4;
      break;
    }
  case 26:
    {
      pWaveStruct->length = 14;
      pFilterCoef = hm2_26;
      pFilterCoefMirror = h2+2;
      break;
    }
  case 28:
    {
      pWaveStruct->length = 18;
      pFilterCoef = hm2_28;
      pFilterCoefMirror = h2;
      break;
    }
  case 31:
    {
      pWaveStruct->length = 4;
      pFilterCoef = hm3_31;
      pFilterCoefMirror = h3+8;
      break;
    }
  case 33:
    {
      pWaveStruct->length = 8;
      pFilterCoef = hm3_33;
      pFilterCoefMirror = h3+6;
      break;
    }
  case 35:
    {
      pWaveStruct->length = 12;
      pFilterCoef = hm3_35;
      pFilterCoefMirror = h3+4;
      break;
    }
  case 37:
    {
      pWaveStruct->length = 16;
      pFilterCoef = hm3_37;
      pFilterCoefMirror = h3+2;
      break;
    }
  case 39:
    {
      pWaveStruct->length = 20;
      pFilterCoef = hm3_39;
      pFilterCoefMirror = h3;
      break;
    }
  case 44:
	{
      pWaveStruct->length = 10;
      pFilterCoef = hm4_44;
      pFilterCoefMirror = h4;
      break;
	}
  case 55:
	{
      pWaveStruct->length = 12;
      pFilterCoef = hm5_55;
      pFilterCoefMirror = h5;
      break;
	}
  case 68:
	{
      pWaveStruct->length = 18;
      pFilterCoef = hm6_68;
      pFilterCoefMirror = h6;
      break;
	}
  default:
    break;
  }

  wrev(pFilterCoef, pWaveStruct->length, 
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(pFilterCoefMirror, pWaveStruct->length, 
	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;
  return;
}


void
sp_bior_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
  double *pFilterCoef, *pFilterCoefMirror;

  
  switch (member) {
  case 11:
    {
      pWaveStruct->length = 2;
      pFilterCoef = h1+4;
      pFilterCoefMirror = hm1_11;
      break;
    }
  case 13:
    {
      pWaveStruct->length = 6;
      pFilterCoef = h1+2;
      pFilterCoefMirror = hm1_13;
      break;
    }
  case 15:
    {
      pWaveStruct->length = 10;
      pFilterCoef = h1;
      pFilterCoefMirror = hm1_15;
      break;
    }
  case 22:
    {
      pWaveStruct->length = 6;
      pFilterCoef = h2+6;
      pFilterCoefMirror = hm2_22;
      break;
    }
  case 24:
    {
      pWaveStruct->length = 10;
      pFilterCoef = h2+4;
      pFilterCoefMirror = hm2_24;
      break;
    }
  case 26:
    {
      pWaveStruct->length = 14;
      pFilterCoef = h2+2;
      pFilterCoefMirror = hm2_26;
      break;
    }
  case 28:
    {
      pWaveStruct->length = 18;
      pFilterCoef = h2;
      pFilterCoefMirror = hm2_28;
      break;
    }
  case 31:
    {
      pWaveStruct->length = 4;
      pFilterCoef = h3+8;
      pFilterCoefMirror = hm3_31;
      break;
    }
  case 33:
    {
      pWaveStruct->length = 8;
      pFilterCoef = h3+6;
      pFilterCoefMirror = hm3_33;
      break;
    }
  case 35:
    {
      pWaveStruct->length = 12;
      pFilterCoef = h3+4;
      pFilterCoefMirror = hm3_35;
      break;
    }
  case 37:
    {
      pWaveStruct->length = 16;
      pFilterCoef = h3+2;
      pFilterCoefMirror = hm3_37;
      break;
    }
  case 39:
    {
      pWaveStruct->length = 20;
      pFilterCoef = h3;
      pFilterCoefMirror = hm3_39;
      break;
    }
  case 44:
	{
      pWaveStruct->length = 10;
      pFilterCoef = h4;
      pFilterCoefMirror = hm4_44;
      break;
	}
  case 55:
	{
      pWaveStruct->length = 12;
      pFilterCoef = h5;
      pFilterCoefMirror = hm5_55;
      break;
	}
  case 68:
	{
      pWaveStruct->length = 18;
      pFilterCoef = h6;
      pFilterCoefMirror = hm6_68;
      break;
	}
  default:
    break;
  };

  verbatim_copy(pFilterCoef, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(pFilterCoefMirror, pWaveStruct->length,
	      HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;
  return;
}


void
sp_rbior_analysis_initialize (int member, swt_wavelet *pWaveStruct)
{
  double *pFilterCoef, *pFilterCoefMirror;

  switch (member) {
  case 11:
    {
      pWaveStruct->length = 2;
      pFilterCoef = h1+4;
	  pFilterCoefMirror = hm1_11;
      break;
    }
  case 13:
    {
      pWaveStruct->length = 6;
      pFilterCoef = h1+2;
	  pFilterCoefMirror = hm1_13;
      break;
    }
  case 15:
    {
      pWaveStruct->length = 10;
      pFilterCoef = h1;
	  pFilterCoefMirror = hm1_15;
      break;
    }
  case 22:
    {
      pWaveStruct->length = 6;
      pFilterCoef = h2+6;
	  pFilterCoefMirror = hm2_22;
      break;
    }
  case 24:
    {
      pWaveStruct->length = 10;
      pFilterCoef = h2+4;
	  pFilterCoefMirror = hm2_24;
      break;
    }
  case 26:
    {
      pWaveStruct->length = 14;
      pFilterCoef = h2+2;
	  pFilterCoefMirror = hm2_26;
      break;
    }
  case 28:
    {
      pWaveStruct->length = 18;
      pFilterCoef = h2;
	  pFilterCoefMirror = hm2_28;
      break;
    }
  case 31:
    {
      pWaveStruct->length = 4;
      pFilterCoef = h3+8;
	  pFilterCoefMirror = hm3_31;
      break;
    }
  case 33:
    {
      pWaveStruct->length = 8;
      pFilterCoef = h3+6;
	  pFilterCoefMirror = hm3_33;
      break;
    }
  case 35:
    {
      pWaveStruct->length = 12;
      pFilterCoef = h3+4;
	  pFilterCoefMirror = hm3_35;
      break;
    }
  case 37:
    {
      pWaveStruct->length = 16;
      pFilterCoef = h3+2;
	  pFilterCoefMirror = hm3_37;
      break;
    }
  case 39:
    {
      pWaveStruct->length = 20;
      pFilterCoef = h3;
	  pFilterCoefMirror = hm3_39;
      break;
    }
  case 44:
	{
      pWaveStruct->length = 10;
      pFilterCoef = h4;
	  pFilterCoefMirror = hm4_44;
      break;
	}
  case 55:
	{
      pWaveStruct->length = 12;
      pFilterCoef = h5;
	  pFilterCoefMirror = hm5_55;
      break;
	}
  case 68:
	{
      pWaveStruct->length = 18;
      pFilterCoef = h6;
	  pFilterCoefMirror = hm6_68;
      break;
	}
  default:
    break;
  }

  wrev(pFilterCoef, pWaveStruct->length, 
       LowDecomFilCoef, pWaveStruct->length);
  qmf_wrev(pFilterCoefMirror, pWaveStruct->length, 
	   HiDecomFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowDecomFilCoef;
  pWaveStruct->pHiPass = HiDecomFilCoef;
  return;
}


void
sp_rbior_synthesis_initialize (int member, swt_wavelet *pWaveStruct)
{
  double *pFilterCoef, *pFilterCoefMirror;

  
  switch (member) {
  case 11:
    {
      pWaveStruct->length = 2;
      pFilterCoef = hm1_11;
	  pFilterCoefMirror = h1+4;
      break;
    }
  case 13:
    {
      pWaveStruct->length = 6;
      pFilterCoef = hm1_13;
	  pFilterCoefMirror = h1+2;
      break;
    }
  case 15:
    {
      pWaveStruct->length = 10;
      pFilterCoef = hm1_15;
	  pFilterCoefMirror = h1;
      break;
    }
  case 22:
    {
      pWaveStruct->length = 6;
      pFilterCoef = hm2_22;
	  pFilterCoefMirror = h2+6;
      break;
    }
  case 24:
    {
      pWaveStruct->length = 10;
      pFilterCoef = hm2_24;
	  pFilterCoefMirror = h2+4;
      break;
    }
  case 26:
    {
      pWaveStruct->length = 14;
      pFilterCoef = hm2_26;
	  pFilterCoefMirror = h2+2;
      break;
    }
  case 28:
    {
      pWaveStruct->length = 18;
      pFilterCoef = hm2_28;
	  pFilterCoefMirror = h2;
      break;
    }
  case 31:
    {
      pWaveStruct->length = 4;
      pFilterCoef = hm3_31;
	  pFilterCoefMirror = h3+8;
      break;
    }
  case 33:
    {
      pWaveStruct->length = 8;
      pFilterCoef = hm3_33;
	  pFilterCoefMirror = h3+6;
      break;
    }
  case 35:
    {
      pWaveStruct->length = 12;
      pFilterCoef = hm3_35;
	  pFilterCoefMirror = h3+4;
      break;
    }
  case 37:
    {
      pWaveStruct->length = 16;
      pFilterCoef = hm3_37;
	  pFilterCoefMirror = h3+2;
      break;
    }
  case 39:
    {
      pWaveStruct->length = 20;
      pFilterCoef = hm3_39;
	  pFilterCoefMirror = h3;
      break;
    }
  case 44:
	{
      pWaveStruct->length = 10;
      pFilterCoef = hm4_44;
	  pFilterCoefMirror = h4;
      break;
	}
  case 55:
	{
      pWaveStruct->length = 12;
      pFilterCoef = hm5_55;
	  pFilterCoefMirror = h5;
      break;
	}
  case 68:
	{
      pWaveStruct->length = 18;
      pFilterCoef = hm6_68;
	  pFilterCoefMirror = h6;
      break;
	}
  default:
    break;
  };

  verbatim_copy(pFilterCoef, pWaveStruct->length,
		LowReconFilCoef, pWaveStruct->length);
  qmf_even(pFilterCoefMirror, pWaveStruct->length,
	      HiReconFilCoef, pWaveStruct->length);
  pWaveStruct->pLowPass = LowReconFilCoef;
  pWaveStruct->pHiPass = HiReconFilCoef;
  return;
}
