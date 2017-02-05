// year to be computed
#define YEAR 1989

// constants
// reference density
#define RHO 1027
// constant horizontal scale
#define WLNGTH 400e3

// bounding box
#define LATMIN -65
#define LATMAX 65
#define LONMIN 0
#define LONMAX 360

// flags
// debug flag
#define DBGFLG 3
// stress correction flag
#define STRSCOR 1
// run the slab model

#define SLABFLG 1
// run the hybrid extension
#define HYBRIDFLG 1
// atocorrelation flag
#define ACFLG 0

// paths to dynamic data;
// years are replaced with %d
// months are replaced with %02d

//#define MLDPATH "/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
//#define NPATH "/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_N_month%02d_regrid.nc"
//#define STRSPATH "/run/media/georg/TRANSCEND/EIS2/NCEP-CFSR/stress/p2/wndstrs.gdas.%d%02d.fast.nc"

#define MLDPATH "../data/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
#define NPATH "../data/MIMOC/MIMOC_ML_v2.2_N_month%02d_regrid.nc"
#define STRSPATH "../data/NCEP-CFSR/stress/wndstrs.gdas.%d%02d.fast.nc"

//#define OUTPATH_N "/run/media/georg/4TB-ext/hybrid/eis3c_test/test_sym_20_n_%02d.nc"
//#define AUXPATH_N "/run/media/georg/4TB-ext/hybrid/eis3c_test/test_sym_20_aux_n_%02d.nc"
//#define OUTPATH_S "/run/media/georg/4TB-ext/hybrid/eis3c_test/test_sym_20_s_%02d.nc"
//#define AUXPATH_S "/run/media/georg/4TB-ext/hybrid/eis3c_test/test_sym_20_aux_s_%02d.nc"

#define OUTPATH_N "../results/%d_%02d_n.nc"
#define AUXPATH_N "../results/%d_%02d_aux_n.nc"
#define OUTPATH_S "../results/%d_%02d_s.nc"
#define AUXPATH_S "../results/%d_%02d_aux_s.nc"

//#define MLDPATH "/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
//#define MLDPATH "/cygdrive/e/Data/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
//#define STRSPATH "/cygdrice/e/Data/NCEP-CFSR/stress/p2/wndstrs.gdas.%d%02d.fast.nc"
//#define NPATH "/cygdrive/e/Data/MIMOC/MIMOC_ML_v2.2_N2_month%02d.nc"
//#define OUTPATH "/cygdrive/e/Data/results/eis3/test_lh_%02d.nc"
//#define AUXPATH "/cygdrive/e/Data/results/eis3/test_lh_aux_%02d.nc"

