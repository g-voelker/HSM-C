#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "header.h"

int getstring(char* buffer, size_t buflen, FILE *fp) {
  // this is a wrapper for fgets
  if (fgets(buffer, buflen, fp) != 0) {
    // read length of string in buffer
    size_t len = strlen(buffer);

    // deal with windows line ends and endline characters
    char* winchar = strrchr(buffer, '\r');
    if(winchar && winchar[1]=='\n' && winchar[2]=='\0') {
      winchar[0] = '\0';
    }
    return 1;
  } else {
    return 0;
  }
}

double retval(char* buffer, char* expect){
  // allocate variables
  double value;
  size_t lenExpect, lenBuffer, lenValue;
  int ns, match = 1;

  // get length of expected variable
  lenExpect = strlen(expect);

  // get length of total buffer and value
  lenBuffer = strlen(buffer);
  lenValue = lenBuffer - lenExpect;

  // get length of variable name in setup file
  size_t lenVar;
  for (ns=0; ns<lenBuffer; ns++) {
    if (buffer[ns] == ' ') {
      lenVar = (size_t) ns;
    }
  }

  // check if the length of the variable in file matches length of expected variable
  if (lenVar==lenExpect) {
    // check if characters of the line match the expected variable
    for (ns = 0; ns < lenExpect; ns++) {
      // if any character does not match: set (int) match to 0 and exit loop
      if (buffer[ns] != expect[ns]) {
        match = 0;
        break;
      }
    }
  } else {
    // if the length is wrong we do not have a match.
    match = 0;
  }

  // read value if there is a match
  // there is no heuristic check for the success of reading the right value!
  if (match==1) {

    // allocate a buffer for the value and fill with characters
    char bufferValue[lenValue];
    for (ns = 0; ns<lenValue; ns++)
      bufferValue[ns] = buffer[ns + lenExpect];

    // get value as a double
    value = strtod(bufferValue, NULL);

    // return value
    return(value);

  } else {

    // if the variable is not the expected throw an error
    // replace this error with the value error macro in the main code
    VARERR(expect, buffer);fflush(NULL);
  }
}

void map(double* values, double *params, int valLen, int pathLen) {
  // map the values to the params array
  // if you add more initial parameters you will need to assign a slot in the parameters array here.
  // at the moment the params array declared in the main file has exactly the needed length
  // but leaves space to insert more constants in between.
  int nmap;

  // mapping from the first (values) to the second index (params)
  int mapping[VALLEN - PATHLEN][2] = {{0, 0},
                                      {1, 1},
                                      {2, 2},
                                      {3, 3},
                                      {4, 4},
                                      {5, 10},
                                      {6, 11},
                                      {7, 12},
                                      {8, 13},
                                      {9, 14},
                                      {10, 15},
                                      {11, 16}
  };

  // assign values to params (double)
  for (nmap=0; nmap < valLen - pathLen; nmap++) {
    params[mapping[nmap][1]] = values[mapping[nmap][0]];
  }
}

void readtxt(double* params, char **paths, int valLen, int pathLen) {
  FILE *ptrFile;
  ptrFile = fopen("input.txt", "r");

  // proceed only if file exists
  if (ptrFile != NULL) {

    // set array with expected parameter names and (double) values
    char *expect[valLen];
    double values[valLen];
    // flags
    expect[0] = "SLABFLG";
    expect[1] = "HYBRIDFLG";
    expect[2] = "DIVFLG";
    expect[3] = "ACFLG";
    expect[4] = "STRSCOR";

    // parameters
    expect[5] = "YEAR";
    expect[6] = "LATMIN";
    expect[7] = "LATMAX";
    expect[8] = "LONMIN";
    expect[9] = "LONMAX";

    // constants
    expect[10] = "RHO";
    expect[11] = "WLNGTH";

    // file paths
    expect[12] = "NPATH"; // paths[0]
    expect[13] = "STRSPATH"; // ...
    expect[14] = "MLDPATH";
    expect[15] = "AUXPATH_N"; // paths[3]
    expect[16] = "AUXPATH_S";
    expect[17] = "OUTPATH_N"; // ...
    expect[18] = "OUTPATH_S"; // paths[6]

    // one can assign more expected parameters here, just make sure the values and the expect arrays are long enough
    // also make sure you map new entries correctly in map()

    // declare the buffer the text lines will be read into
    size_t buflen = MAXCHARLEN;
    char buffer[MAXCHARLEN];

    // declare auxiliary integers
    int nread = 0, shift = 0;

    // read the whole file
    while (getstring(buffer, buflen, ptrFile) != 0) {
      // ignore comment lines beginning with # as well as blank lines
      if ((buffer[0] == '#') || (buffer[0] == '\0') || (buffer[0] == '\n')) {
        shift++;
      } else {
        // put the next parameter into the corresponding slot in the parameter array
        if (nread - shift < valLen - pathLen) {
          values[nread - shift] = retval(buffer, expect[nread - shift]);
        } else {
          paths[nread - shift + pathLen - valLen] = buffer;
        }
      }
      nread++;
    }

    // close the file
    fclose(ptrFile);

    // map the values to the params array
    map(values, params, valLen, pathLen);

  } else {
    // if file does not exists throw an error and quit.
    IOERR;
  }
}
