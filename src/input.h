// year to be computed
#define YEAR 1989

// constants
#define RHO 1027

// bounding box
#define LATMIN -65
#define LATMAX 65
#define LONMIN 0
#define LONMAX 360

// flags
// debug flag
#define DBGFLG 1
// stress correction flag
#define STRSCOR 1

// paths to dynamic data;
// years are replaced with %d
// months are replaced with %02d

#define MLDPATH "/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
#define NPATH "/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_N_month%02d_regrid.nc"
#define STRSPATH "/run/media/georg/TRANSCEND/EIS2/NCEP-CFSR/stress/p2/wndstrs.gdas.%d%02d.fast.nc"
#define OUTPATH "/run/media/georg/TRANSCEND/EIS2/results/eis3_test/test_%02d.nc"
#define AUXPATH "/run/media/georg/TRANSCEND/EIS2/results/eis3_test/test_aux_%02d.nc"
//#define MLDPATH "/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
//#define MLDPATH "/cygdrive/e/Data/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
//#define STRSPATH "/cygdrice/e/Data/NCEP-CFSR/stress/p2/wndstrs.gdas.%d%02d.fast.nc"
//#define NPATH "/cygdrive/e/Data/MIMOC/MIMOC_ML_v2.2_N2_month%02d.nc"
//#define OUTPATH "/cygdrive/e/Data/results/eis3/test_lh_%02d.nc"
//#define AUXPATH "/cygdrive/e/Data/results/eis3/test_lh_aux_%02d.nc"

