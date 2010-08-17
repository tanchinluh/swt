/*
 * -------------------------------------------------------------------------
 * cwt.h -- Continuous Wavelet Transform Declaration.
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
#include "dwt.h"


/*********************************************
 * Macros
 ********************************************/

#define REAL    0
#define COMPLEX 1

#define PHI_ONLY     0
#define PSI_ONLY     1
#define PHI_PSI_BOTH 2

#define SINUS           0
#define POISSON         1
#define MEXICAN_HAT     2
#define MORLET          3
#define DOGAUSS         4
#define CMORLET         5
#define SHANNON         6
#define FBSP            7
#define CAUCHY          8
#define GAUSS           9
#define CGAUSS          10

/*********************************************
 * Structures
 ********************************************/

typedef void(*WScaleFunc)(double *x, int sigInLength, double *psi, int sigOutLength, double ys);


typedef struct {
	char wname[20];
	int     realOrComplex;
	int     family;
	int     phipsi;
	double  lb;
	double  ub;
	double cpsi;
	WScaleFunc scalef;
} cwt_identity;

/*********************************************
 * Global Variables
 ********************************************/

extern cwt_identity ci[];

/*********************************************
 * Function Declaration
 ********************************************/

extern void sinus_form_validate(int *errCode);
extern void sinus_content_validate(int *errCode, int l1, int l2, int l3);
extern void sinus(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void poisson_form_validate(int *errCode);
extern void poisson_content_validate(int *errCode, int l1, int l2, int l3);
extern void poisson(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void mexihat(double *x, int sigInLength, double *psi, int sigOutLength, double ys);

extern void mexihat_form_validate(int *errCode);
extern void mexihat_content_validate(int *errCode, int l1, int l2, int l3);
extern void morlet(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void morlet_form_validate(int *errCode);
extern void morlet_content_validate(int *errCode, int l1, int l2, int l3);
extern void DOGauss(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void DOGauss_form_validate(int *errCode);
extern void DOGauss_content_validate(int *errCode, int l1, int l2, int l3);
extern void Gauss(double *x, int sigInLength, double *psi, int sigOutLength, int n, double ys);
extern void Gaus1(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gaus2(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gaus3(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gaus4(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gaus5(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gaus6(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gaus7(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gaus8(double *x, int sigInLength, double *psi, int sigOutLength, double ys);
extern void Gauss_form_validate(int *errCode);
extern void Gauss_content_validate(int *errCode, int l1, int l2, int l3, int l4);
extern void cgauss(double *x, int sigInLength, int p, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau1(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau1_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
extern void cgau2(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau2_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
extern void cgau3(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau3_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
extern void cgau4(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau4_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);

extern void cgau5(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau5_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
extern void cgau6(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau6_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
extern void cgau7(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau7_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);
extern void cgau8(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cgau8_packet(double *x, int sigInLength, double *f, int sigOutLength, double ys);

extern void cmorlet_form_validate(int *errCode);
extern void cmorlet_content_validate(int *errCode, int l1, int l2, int l3, int l4, int l5);
extern void cmorlet(double *x, int sigInLength, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys);
extern void cmorlet_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys);
extern void shanwavf_form_validate(int *errCode);
extern void shanwavf_content_validate(int *errCode, int l1, int l2, int l3, int l4, int l5);
extern void shanwavf(double *x, int sigInLength, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys);
extern void shanwavf_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys);
extern void fbspwavf(double *x, int sigInLength,int m, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys);
extern void fbspwavf_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys);
extern void fbspwavf_form_validate(int *errCode);
extern void fbspwavf_content_validate(int *errCode, int l1, int l2, int l3, int l4, int l5, int l6);
extern void cauchy(double *x, int sigInLength, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys);
extern void cauchy_neo(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys);
extern void cauchy_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys);
extern void cauchy_form_validate(int *errCode);
extern void cauchy_content_validate(int *errCode, int l1, int l2, int l3);

extern void meyeraux(double x, double *y);

extern void wavefun_form_validate(int *errCode);
extern void wavefun_content_validate(int *errCode, int l1, int l2);
extern void wavefun2_form_validate(int *errCode);
extern void wavefun2_content_validate(int *errCode, int l1, int l2);


extern void cwt_fun_parser(char *wname, int *ind);
extern void cwt_len_cal (int sigInLength, int scale, int *sigOutLength, double *delta);
extern void full_range_scalef (char *wname, double *f, int sigOutLength);
extern void scale_real (double *f, int sigInLength, double delta, double *fout, int sigOutLength);
//extern void scale_complex (double *f, int sigInLength, double delta, double *fout, int sigOutLength);
extern void cwt_conv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength);
extern void cwt_iconv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength);
extern void cwt_conv_complex (double *sigIn, int sigInLength, double *fr, double *fi, int filterLen, double *sigOutR, double *sigOutI, int sigOutLength);
extern void cwt_conv_complex_complex (double *a, double *b, int sigInLength,double *c, double *d,
							   int filterLen, double *sigOutR, double *sigOutI, int sigOutLength);

extern void cwt_form_validate(int *errCode, int *flow);
extern void cwt_content_validate(int *errCode, int l1, int l2, int l3);