#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#define MAXCHARLEN 100

// the length of the DFTs is set to 2^5 * 7^3 in accordance with the maximum cutoff for long damping times.
// this value is set to preserve performance while making sure the cutoff is long enough
// edit only if you know what you do.
#define DFFT_LEN 10976