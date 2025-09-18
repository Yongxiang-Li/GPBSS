#include <stdio.h>
#include "mex.h"

/* Declare global variables */
double **a, *b;
mwIndex **each_ir, **each_jc, *ir, *each_rows, cnt;
mwSize N;

/* Define function which recursively computes column in output matrix */
void compute_output_column(mwIndex c, mwIndex n, double x, mwIndex ind) {
	double x_new;  // mexPrintf("%d\n", ind);
	mwIndex i, ind_new;
	for (i = each_jc[n][c]; i < each_jc[n][c + 1]; ++i) {
		x_new = x*a[n][i];
		ind_new = ind*each_rows[n] + each_ir[n][i];
		if (n < N - 1) {
			compute_output_column(c, n + 1, x_new, ind_new);
		}
		else {
			b[cnt] = x_new;
			ir[cnt] = ind_new;
			++cnt;
		}
	}
}

/* mex interface */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Declare other variables */
	mwSize c, n, b_nnz, no_cols;

	/* Get input variables */
	N = nrhs;
	a = malloc(N*sizeof(double *));
	each_ir = malloc(N*sizeof(mwIndex *));
	each_jc = malloc(N*sizeof(mwIndex *));
	each_rows = malloc(N*sizeof(mwIndex));
	for (n = 0; n < N; ++n) {
		a[n] = mxGetPr(prhs[n]);
		each_ir[n] = mxGetIr(prhs[n]);
		each_jc[n] = mxGetJc(prhs[n]);
		each_rows[n] = mxGetM(prhs[n]);
	}
	no_cols = mxGetN(prhs[0]); // mexPrintf("%d\n", no_rows);

	/* Compute no rows in output matrix */
	mwIndex no_rows = 1;
	for (n = 0; n < N; ++n) {
		no_rows *= mxGetM(prhs[n]);
	}

	/* Compute nnz in output matrix */
	b_nnz = 1;
	for (c = 0; c < no_cols; ++c) {
		mwIndex prod = 1;
		for (n = 0; n < N; ++n){
			prod *= each_jc[n][c + 1] - each_jc[n][c];
		}
		b_nnz += prod;
	}
	/* Create sparse output matrix */
	plhs[0] = mxCreateSparse(no_rows, no_cols, b_nnz, mxREAL);
	b = mxGetPr(plhs[0]);
	ir = mxGetIr(plhs[0]);
	mwIndex *jc = mxGetJc(plhs[0]);

	/* Compute jc for output matrix */
	jc[0] = 0;
	for (c = 0; c < no_cols; ++c) {
		mwIndex prod = 1;
		for (n = 0; n < N; ++n) {
			prod *= each_jc[n][c + 1] - each_jc[n][c];
		}
		jc[c + 1] = jc[c] + prod;
	}

	/* Compute non-zero elements and ir vector for output matrix */
	cnt = 0;
	for (c = 0; c < no_cols; ++c) {
		compute_output_column(c, 0, 1.0, 0);
	}

	/* Free dynamically allocated memory */
	free(each_rows);
	free(each_jc);
	free(each_ir);
	free(a);
}