// year to be computed
#define YEAR 1996

// constants
#define RHO 1027

// bounding box
#define LATMIN 20
#define LATMAX 21
#define LONMIN 320
#define LONMAX 321

// flags
#define DBGFLG 0
#define STRSCOR 1

// paths to dynamic data;
// years are replaced with %d
// months are replaced with %02d
#define MLDPATH "/run/media/georg/TRANSCEND/EIS2/MIMOC/MIMOC_ML_v2.2_CT_SA_MLP_month%02d_regrid.nc"
#define STRSPATH "/run/media/georg/TRANSCEND/EIS2/NCEP-CFSR/stress/p2/wndstrs.gdas.%d%02d.fast.nc"
#define OUTPATH "/run/media/georg/TRANSCEND/EIS2/results/eis3_test/test_lh_%02d.nc"