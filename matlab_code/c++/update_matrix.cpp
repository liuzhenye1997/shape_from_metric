#include "mex.h"
#include"math.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize actual_number_of_non_zeros;
	double *ptr;
	double *value;
#if MX_HAS_INTERLEAVED_COMPLEX
	mxComplexDouble *complexPtr;
#endif

	/* Check for proper number of input and output arguments */
	if (nrhs != 2 && nrhs != 3)
	{
		mexErrMsgIdAndTxt("MATLAB:updatevalue:invalidNumInputs",
			"Two or three input argument required.");
	}

	if (!mxIsSparse(prhs[0]))
	{
		mexErrMsgIdAndTxt("MATLAB:updatevalue:invalidInputSparisty",
			"Input argument must be sparse\n");
	}
	if (nrhs == 2)
	{
		ptr = mxGetPr(prhs[0]);
		value = mxGetPr(prhs[1]);
		for (int i = 0; i < mxGetNzmax(prhs[0]); i++)
		{
			ptr[i] = value[i];
		}
	}
	else
	{
		mxInt32 *index;
		ptr = mxGetPr(prhs[0]);
		value = mxGetPr(prhs[1]);
		index = mxGetInt32s(prhs[2]);
		for (int i = 0; i < mxGetNzmax(prhs[0]); i++)
		{
	/*		mxprintf("%ld %lf\n", index[i], index[i]);*/
		/*	if (index[i] >= mxGetNzmax(prhs[0]))
			{
				break;
			}*/
			
			ptr[i] = value[index[i]];
		}
	}

}
