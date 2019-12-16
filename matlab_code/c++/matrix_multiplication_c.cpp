
# include "mex.h"
void mexArrayAdd(const int Row, const int Col, const double *A, const double *B, double *C,bool transpose) 
{
	int i, j,k,l;
	if (transpose == false)
	{
		for (l = 0; l < Row; l++)
		{		
			for (j = 0; j < 4; j++)
			{
				for (k = 0; k < 4; k++)
				{
					C[l + j * Row + k * Row * 4] = 0;
					for (i = 0; i < 4; i++)
					{
						C[ j * Row + k * Row * 4+l] = C[j * Row + k * Row * 4 + l] + A[l + j * Row + i * Row * 4] * B[l + i * Row + k * Row * 4];
					}
				}
			}
		}
	}
	else
	{
		for (l = 0; l < Row; l++)
		{
			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 4; j++)
				{
					for (k = 0; k < 4; k++)
					{
						C[j * Row + k * Row * 4 + l] = C[j * Row + k * Row * 4 + l] + A[l + i * Row + j * Row * 4] * B[l + i * Row + k * Row * 4];
					}
				}
			}
		}
	}

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *A, *B, *C;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 1) {
		mexErrMsgTxt("One output required.");
	}
	else if (nrhs > 10) {
		mexErrMsgTxt("Two input required.");
	}
	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])) {
		mexErrMsgTxt("Input Array must be double.");
	}
	//获得矩阵的行数
	RowA = mxGetM(prhs[0]);
	//获得矩阵的列数
	ColA = mxGetN(prhs[0]);
	RowB = mxGetM(prhs[1]);
	ColB = mxGetN(prhs[1]);
	//判断行列是否相等
	if (RowA != RowB || ColA != ColB) {
		mexErrMsgTxt("Rows and Cols must be same.");
	}
	//获取指向输入参数的指针
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
	transpose=*mxGetPr(prhs[2]);
	//生成输出参量的mxArray
	plhs[0] = mxCreateDoubleMatrix(16, RowA, mxREAL);
	//获取指向输出参数的指针
	C = mxGetPr(plhs[0]);
	mexArrayAdd(RowA, ColA, A, B, C, transpose);

}

