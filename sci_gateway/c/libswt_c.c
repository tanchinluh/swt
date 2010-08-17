#include <mex.h> 
#include <sci_gateway.h>
#include <api_scilab.h>
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc int_wrev;
extern Gatefunc int_wrev2;
extern Gatefunc int_qmf;
extern Gatefunc int_conv;
extern Gatefunc int_iconv;
extern Gatefunc int_dyaddown;
extern Gatefunc int_dyadup;
extern Gatefunc int_wkeep;
extern Gatefunc int_wextend;
extern Gatefunc int_wcodemat;
extern Gatefunc int_ind2rgb;
extern Gatefunc int_mat3Dtran;
extern Gatefunc int_wrev3;
extern Gatefunc int_orthfilt;
extern Gatefunc int_dbwavf;
extern Gatefunc int_coifwavf;
extern Gatefunc int_symwavf;
extern Gatefunc int_legdwavf;
extern Gatefunc int_biorwavf;
extern Gatefunc int_rbiorwavf;
extern Gatefunc int_biorfilt;
extern Gatefunc int_wfilters;
extern Gatefunc int_wmaxlev;
extern Gatefunc int_dwtmode;
extern Gatefunc int_dwt;
extern Gatefunc int_idwt;
extern Gatefunc int_wavedec;
extern Gatefunc int_waverec;
extern Gatefunc int_wrcoef;
extern Gatefunc int_appcoef;
extern Gatefunc int_detcoef;
extern Gatefunc int_wenergy;
extern Gatefunc int_upcoef;
extern Gatefunc int_upwlev;
extern Gatefunc int_dwt2;
extern Gatefunc int_idwt2;
extern Gatefunc int_wavedec2;
extern Gatefunc int_waverec2;
extern Gatefunc int_wenergy2;
extern Gatefunc int_appcoef2;
extern Gatefunc int_detcoef2;
extern Gatefunc int_wrcoef2;
extern Gatefunc int_upwlev2;
extern Gatefunc int_upcoef2;
extern Gatefunc int_swt;
extern Gatefunc int_iswt;
extern Gatefunc int_swt2;
extern Gatefunc int_iswt2;
extern Gatefunc int_sinus;
extern Gatefunc int_poisson;
extern Gatefunc int_mexihat;
extern Gatefunc int_morlet;
extern Gatefunc int_DOGauss;
extern Gatefunc int_cmorlet;
extern Gatefunc int_shanwavf;
extern Gatefunc int_fbspwavf;
extern Gatefunc int_cauchy;
extern Gatefunc int_Gauswavf;
extern Gatefunc int_cgauss;
extern Gatefunc int_wavefun;
extern Gatefunc int_wavefun2;
extern Gatefunc int_cwt;
extern Gatefunc int_dwt3;
extern Gatefunc int_idwt3;
extern Gatefunc int_FSfarras;
extern Gatefunc int_dualfilt1;
extern Gatefunc int_dualtree;
extern Gatefunc int_idualtree;
extern Gatefunc int_dualtree2D;
extern Gatefunc int_idualtree2D;
extern Gatefunc int_cplxdual2D;
extern Gatefunc int_icplxdual2D;
extern Gatefunc int_wnorm;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,int_wrev,"wrev"},
  {(Myinterfun)sci_gateway,int_wrev2,"wrev2"},
  {(Myinterfun)sci_gateway,int_qmf,"qmf"},
  {(Myinterfun)sci_gateway,int_conv,"conv"},
  {(Myinterfun)sci_gateway,int_iconv,"iconv"},
  {(Myinterfun)sci_gateway,int_dyaddown,"dyaddown"},
  {(Myinterfun)sci_gateway,int_dyadup,"dyadup"},
  {(Myinterfun)sci_gateway,int_wkeep,"wkeep"},
  {(Myinterfun)sci_gateway,int_wextend,"wextend"},
  {(Myinterfun)sci_gateway,int_wcodemat,"wcodemat"},
  {(Myinterfun)sci_gateway,int_ind2rgb,"ind2rgb"},
  {(Myinterfun)sci_gateway,int_mat3Dtran,"wrot3"},
  {(Myinterfun)sci_gateway,int_wrev3,"wrev3"},
  {(Myinterfun)sci_gateway,int_orthfilt,"orthfilt"},
  {(Myinterfun)sci_gateway,int_dbwavf,"dbwavf"},
  {(Myinterfun)sci_gateway,int_coifwavf,"coifwavf"},
  {(Myinterfun)sci_gateway,int_symwavf,"symwavf"},
  {(Myinterfun)sci_gateway,int_legdwavf,"legdwavf"},
  {(Myinterfun)sci_gateway,int_biorwavf,"biorwavf"},
  {(Myinterfun)sci_gateway,int_rbiorwavf,"rbiorwavf"},
  {(Myinterfun)sci_gateway,int_biorfilt,"biorfilt"},
  {(Myinterfun)sci_gateway,int_wfilters,"wfilters"},
  {(Myinterfun)sci_gateway,int_wmaxlev,"wmaxlev"},
  {(Myinterfun)sci_gateway,int_dwtmode,"dwtmode"},
  {(Myinterfun)sci_gateway,int_dwt,"dwt"},
  {(Myinterfun)sci_gateway,int_idwt,"idwt"},
  {(Myinterfun)sci_gateway,int_wavedec,"wavedec"},
  {(Myinterfun)sci_gateway,int_waverec,"waverec"},
  {(Myinterfun)sci_gateway,int_wrcoef,"wrcoef"},
  {(Myinterfun)sci_gateway,int_appcoef,"appcoef"},
  {(Myinterfun)sci_gateway,int_detcoef,"detcoef"},
  {(Myinterfun)sci_gateway,int_wenergy,"wenergy"},
  {(Myinterfun)sci_gateway,int_upcoef,"upcoef"},
  {(Myinterfun)sci_gateway,int_upwlev,"upwlev"},
  {(Myinterfun)sci_gateway,int_dwt2,"dwt2"},
  {(Myinterfun)sci_gateway,int_idwt2,"idwt2"},
  {(Myinterfun)sci_gateway,int_wavedec2,"wavedec2"},
  {(Myinterfun)sci_gateway,int_waverec2,"waverec2"},
  {(Myinterfun)sci_gateway,int_wenergy2,"wenergy2"},
  {(Myinterfun)sci_gateway,int_appcoef2,"appcoef2"},
  {(Myinterfun)sci_gateway,int_detcoef2,"detcoef2"},
  {(Myinterfun)sci_gateway,int_wrcoef2,"wrcoef2"},
  {(Myinterfun)sci_gateway,int_upwlev2,"upwlev2"},
  {(Myinterfun)sci_gateway,int_upcoef2,"upcoef2"},
  {(Myinterfun)sci_gateway,int_swt,"swt"},
  {(Myinterfun)sci_gateway,int_iswt,"iswt"},
  {(Myinterfun)sci_gateway,int_swt2,"swt2"},
  {(Myinterfun)sci_gateway,int_iswt2,"iswt2"},
  {(Myinterfun)sci_gateway,int_sinus,"sinus"},
  {(Myinterfun)sci_gateway,int_poisson,"poisson"},
  {(Myinterfun)sci_gateway,int_mexihat,"mexihat"},
  {(Myinterfun)sci_gateway,int_morlet,"morlet"},
  {(Myinterfun)sci_gateway,int_DOGauss,"DOGauss"},
  {(Myinterfun)sci_gateway,int_cmorlet,"cmorwavf"},
  {(Myinterfun)sci_gateway,int_shanwavf,"shanwavf"},
  {(Myinterfun)sci_gateway,int_fbspwavf,"fbspwavf"},
  {(Myinterfun)sci_gateway,int_cauchy,"cauwavf"},
  {(Myinterfun)sci_gateway,int_Gauswavf,"gauswavf"},
  {(Myinterfun)sci_gateway,int_cgauss,"cgauswavf"},
  {(Myinterfun)sci_gateway,int_wavefun,"wavefun"},
  {(Myinterfun)sci_gateway,int_wavefun2,"wavefun2"},
  {(Myinterfun)sci_gateway,int_cwt,"cwt"},
  {(Myinterfun)sci_gateway,int_dwt3,"dwt3"},
  {(Myinterfun)sci_gateway,int_idwt3,"idwt3"},
  {(Myinterfun)sci_gateway,int_FSfarras,"FSfarras"},
  {(Myinterfun)sci_gateway,int_dualfilt1,"dualfilt1"},
  {(Myinterfun)sci_gateway,int_dualtree,"dualtree"},
  {(Myinterfun)sci_gateway,int_idualtree,"idualtree"},
  {(Myinterfun)sci_gateway,int_dualtree2D,"dualtree2D"},
  {(Myinterfun)sci_gateway,int_idualtree2D,"idualtree2D"},
  {(Myinterfun)sci_gateway,int_cplxdual2D,"cplxdual2D"},
  {(Myinterfun)sci_gateway,int_icplxdual2D,"icplxdual2D"},
  {(Myinterfun)sci_gateway,int_wnorm,"wnorm"},
};
 
int C2F(libswt_c)()
{
  Rhs = Max(0, Rhs);
  if (*(Tab[Fin-1].f) != NULL) 
  {
     pvApiCtx->pstName = (char*)Tab[Fin-1].name;
    (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  }
  return 0;
}
