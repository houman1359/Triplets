#include <math.h>
#include "mex.h"
#include "matrix.h"

/* Input Arguments */

#define	B_in          prhs[0]
#define	data_in       prhs[1]
#define	grid_in       prhs[2]

/* Output Arguments */

#define	Ker_grid_out  plhs[0]


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{ 
    const double *B, *data, *grid;
    double *Ker_grid;
    int i, j, k;
    int is_vec;
    mwSize n, m;
    mwSize g_m, g_n, k_m, d_m, d_n, b_m;
    double c1, c2, c3, b1, b2, b3, a;
 
    /* Check for proper number of arguments */
    if (nrhs != 3)
    {
        mexErrMsgTxt ("3 input arguments required.");
    }
    if (nlhs > 1) 
    {
        mexErrMsgTxt ("1 output argument only.");
    }

    /* Assign pointers to the various arguments */ 
    B = mxGetPr (B_in);
    data = mxGetPr (data_in);
    grid = mxGetPr (grid_in);

    g_m = mxGetM (grid_in);
    g_n = mxGetN (grid_in);
    d_m = mxGetM (data_in);
    d_n = mxGetN (data_in);
    b_m = mxGetM (B_in);
    k_m = 3;

    Ker_grid_out = mxCreateDoubleMatrix (3, g_n, mxREAL);
    Ker_grid = mxGetPr (Ker_grid_out);

    is_vec = mxGetN (B_in) == 1 || mxGetM (B_in) == 1;

    for (k = 0; k < g_n; k++)
    {
        for (i = 0; i < 3; i++)
            Ker_grid[i + k_m * k] = 0.0;

        for (i = 0; i < d_n; i++)
        {
            c1 = -data[0 + d_m * i] + grid[0 + g_m * k];
            
            if (is_vec)
            {
                b1 = B[0];
            }
            else
            {
                b1 = B[0 + b_m * i];
            }
            a = exp((-c1*c1)/(2*b1*b1))  / (sqrt(2*M_PI)*b1*d_n);

            Ker_grid[0 + k_m * k] += a;
            Ker_grid[1 + k_m * k] += a * c1;
            Ker_grid[2 + k_m * k] += a * c1 * c1;
        }
    }
}
