typedef struct RGBCLR {short r, g, b;} rgbclr;

extern void char2img(char c, short **cm);

extern void defcolmap(double *icmval,int colmap);
extern int getcolmapnumber();
extern void getcolmap(int colmap, short **cm);

extern rgbclr getcolour();

extern void Img2RGB(short **imat, int nx, int ny, int colmap, 
                    short ***rgbmat, int ni, int nj);

extern void writeGryPSHeader(FILE *fpw,char *sfpw, int nx, int ny);
extern void writeColPSHeader(FILE *fpw,char *sfpw, int nx, int ny);
extern void writePSFoot(FILE *fpw);
extern void writeColPS(char *sfpw, short ***rgbmat, int nx, int ny);
extern void writeGryPS(char *sfpw, short ***rgbmat, int nx, int ny);

extern void writeRGB(char *sfpw, short ***imat, int nx, int ny);
extern void writeJPG(char *sfpw, short ***imat, int nx, int ny);
