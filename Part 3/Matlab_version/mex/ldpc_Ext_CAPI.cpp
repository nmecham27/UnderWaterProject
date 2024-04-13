/*==========================================================
 * mexcpp.cpp - example in MATLAB External Interfaces
 *
 * Illustrates how to use some C++ language features in a MEX-file.
 *
 * The routine simply defines a class, constructs a simple object,
 * and displays the initial values of the internal variables. It
 * then sets the data members of the object based on the input given
 * to the MEX-file and displays the changed values.
 *
 * This file uses the extension .cpp. Other common C++ extensions such
 * as .C, .cc, and .cxx are also supported.
 *
 * The calling syntax is:
 *
 *              mexcpp( num1, num2 )
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2018 The MathWorks, Inc.
 *
 *========================================================*/

#include <iostream>
#include "mex.h"
#include "bLDPC.hpp"




/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    
     uint64_t* ptr = static_cast<uint64_t*>(mxGetData(prhs[0]));
     qcldpc* ldpc_ptr = reinterpret_cast<qcldpc*>(*ptr);
     
     double* doubleCode = mxGetPr(prhs[1]);
     int mesgLen = (*ldpc_ptr).infoLen;
     int codeLen = (*ldpc_ptr).codeLen;
     double *mesg = new double[mesgLen];
     double *code = new double[codeLen];
     for(int i = 0; i < codeLen; ++i)
     {
         code[i] = doubleCode[i];
     }

	 (*ldpc_ptr).extract_mesg(mesg, code);

     plhs[0] = mxCreateDoubleMatrix(mesgLen, 1, mxREAL);
     double *output_p = mxGetPr(plhs[0]);
     for(int i = 0; i < mesgLen; ++i)
         output_p[i] = (double)mesg[i];

     delete []mesg;
     delete []code;
}
