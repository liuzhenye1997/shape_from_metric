
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
	//��þ��������
	RowA = mxGetM(prhs[0]);
	//��þ��������
	ColA = mxGetN(prhs[0]);
	RowB = mxGetM(prhs[1]);
	ColB = mxGetN(prhs[1]);
	//�ж������Ƿ����
	if (RowA != RowB || ColA != ColB) {
		mexErrMsgTxt("Rows and Cols must be same.");
	}
	//��ȡָ�����������ָ��
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
	transpose=*mxGetPr(prhs[2]);
	//�������������mxArray
	plhs[0] = mxCreateDoubleMatrix(16, RowA, mxREAL);
	//��ȡָ�����������ָ��
	C = mxGetPr(plhs[0]);
	mexArrayAdd(RowA, ColA, A, B, C, transpose);

}

