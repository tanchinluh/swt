/*
 * -------------------------------------------------------------------------
 * swt.h -- Declarations for stationary wavelet functions.
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

/*********************************************
 * Function Prototype
 ********************************************/
extern void swt_conv(double *sigIn, int sigInLength,
	           double *approx, int approxLength,
			   double *detail, int detailLength,
			   double *filterLow, double *filterHi,
			   int filterLength);
extern void  swt_out1 (double *sigIn, int sigInLength,
	            double *sigOutMatrix, int rowLength, 
				int colLength, double *filterLow, 
				double *filterHi, int filterLength, int step);

extern void swt_out2 (double *sigIn, int sigInLength,
	            double *approxMatrix, double *detailMatrix,
				int rowLength, int colLength, double *filterLow,
				double *filterHi, int filterLength, int step);

extern void iswt_conv (double *approx, double *detail, int sigInLength,
	             double *sigOut, int sigOutLength, double *filterLow,
				 double *filterHi, int filterLength);
extern void iswt_conv_step (double *approx, double *detail, int sigInLength,
	                  double *sigOut, int sigOutLength, double *filterLow,
				      double *filterHi, int filterLength, int level);

extern void iswt_input1 (double *matrixIn, int rowLength, int colLength,
	              double *sigOut, int sigOutLength, double *filterLow,
				  double *filterHi, int filterLength);

extern void iswt_input2 (double *matrixApproxIn, double *matrixDetailIn,
	               int rowLength, int colLength,
	              double *sigOut, int sigOutLength, double *filterLow,
				  double *filterHi, int filterLength);
extern void swt2_output4(double *matrixIn, int matrixInRow, int matrixInCol,
	               double *matrixOutApprox, double *matrixOutColDetail,
				   double *matrixOutRowDetail, double *matrixOutDetail,
				   int matrixOutRow, int matrixOutCol,
				   double *filterLow, double *filterHi, 
				   int filterLength, int step);
extern void swt2_output4_step(double *matrixIn, int matrixInRow, int matrixInCol,
	               double *matrixOutApprox, double *matrixOutColDetail,
				   double *matrixOutRowDetail, double *matrixOutDetail,
				   int matrixOutRow, int matrixOutCol,
				   double *filterLow, double *filterHi, 
				   int filterLength, int step);
extern void swt2_output1_step(double *matrixIn, int matrixInRow, 
	                    int matrixInCol,  double *matrixOut,
				        int matrixOutRow, int matrixOutCol,
				        double *filterLow, double *filterHi, 
				        int filterLength, int step);
extern void iswt2(double *matrixInApprox, double *matrixInColDetail,
		   double *matrixInRowDetail, double *matrixInDetail,
		   int matrixInRow, int matrixInCol,
		   double *matrixOut, int matrixOutRow, int matrixOutCol,
		   double *filterLow, double *filterHi, 
		   int filterLength, int step);
extern void iswt2_input4_step(double *matrixInApprox, double *matrixInColDetail,
		               double *matrixInRowDetail, double *matrixInDetail,
		               int matrixInRow, int matrixInCol,
		               double *matrixOut, int matrixOutRow, int matrixOutCol,
		               double *filterLow, double *filterHi, 
		               int filterLength, int step);
extern void iswt2_input1_step(double *matrixIn,  int matrixInRow, int matrixInCol,
		               double *matrixOut, int matrixOutRow, int matrixOutCol,
		               double *filterLow, double *filterHi, 
		               int filterLength, int step);
 
 
