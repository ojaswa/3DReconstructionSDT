/*
Author: Ojaswa Sharma <ojaswa@iiitd.ac.in>
Original code  for MVC interpolation from VTK MVC Triangle mesh interpolator
Thu 28 May 2015
*/

#include "mex.h"
#include <math.h>
//#include <omp.h>

typedef unsigned int uint;
#define PI 3.141592653589793

extern "C" int mxUnshareArray(mxArray *array_ptr, int level);
void computeMVC(double* x, double * vertices, int n_vertices, uint * triangles, int n_triangles, double *weights, double *uVec, double *dist);



// Utility functions
inline double norm(const double v[3]){return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);}
inline double distance2BetweenPoints(const double x[3], const double y[3])
{
    return ( ( x[0] - y[0] ) * ( x[0] - y[0] )
    + ( x[1] - y[1] ) * ( x[1] - y[1] )
    + ( x[2] - y[2] ) * ( x[2] - y[2] ) );
}
inline double determinant3x3(const double c1[3], const double c2[3], const double c3[3])
{
    return c1[0] * c2[1] * c3[2] + c2[0] * c3[1] * c1[2] + c3[0] * c1[1] * c2[2] -
            c1[0] * c3[1] * c2[2] - c2[0] * c1[1] * c3[2] - c3[0] * c2[1] * c1[2];
}
 
// Entry point
// Input args: (x, vertices, triangles, weights, uvec, dist)
// For faster indexing arrays have these sizes - x: 3x1, points: 3xn, triangles: 3xm, weights: 1xn, uvec: 3xn, dist: 1xn
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if(nrhs != 3)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 3) // No output arguments
        mexErrMsgTxt("Too many output arguments.");
    
    //Get the point for which weighst are to be computed
    int nr_x = mxGetM(prhs[0]);
    int nc_x = mxGetN(prhs[0]);
    if(nc_x != 1 && nr_x != 3)
        mexErrMsgTxt("First argument is not a 3x1 array.");
    double *x = mxGetPr(prhs[0]);
    
    // Get vertices of the mesh
    int nr_vertices = mxGetM(prhs[1]);
    int nc_vertices = mxGetN(prhs[1]);
    if(nr_vertices != 3)
        mexErrMsgTxt("Second argument is not an nx3 array.");
    double *vertices = mxGetPr(prhs[1]);
    
    // Get the triangles
    int nr_triangles = mxGetM(prhs[2]);
    int nc_triangles = mxGetN(prhs[2]);
    if(nr_triangles != 3)
        mexErrMsgTxt("Third argument is not an nx3 array.");
    uint *triangles = (uint *) mxGetData(prhs[2]);
    
    // Get the weights vector --output array
    plhs[0] = mxCreateDoubleMatrix(1, nc_vertices, mxREAL);
    double *weights = mxGetPr(plhs[0]);
    
    // Get the uVec matrix --output array
    plhs[1] = mxCreateDoubleMatrix(3, nc_vertices, mxREAL);
    double *uVec = mxGetPr(plhs[1]);

    // Get the dist vector --output array
    plhs[2] = mxCreateDoubleMatrix(1, nc_vertices, mxREAL);
    double *dist = mxGetPr(plhs[2]);
    
    computeMVC(x, vertices, nc_vertices, triangles, nc_triangles, weights, uVec, dist);
}

void computeMVC(double* x, double * vertices, int n_vertices, uint * triangles, int n_triangles, double *weights, double *uVec, double *dist)
{
    //Initialize weights to zero
    for(int i=0; i<n_vertices; i++)
        weights[i] = 0.0;
    
    static const double eps = 0.000000001;
    bool exit0 = false;
    #pragma omp parallel for
    for(int i=0; i<n_vertices; i++)
    {
        // point-to-vertex vector
        uVec[3*i]   = vertices[3*i  ] - x[0];
        uVec[3*i+1] = vertices[3*i+1] - x[1];
        uVec[3*i+2] = vertices[3*i+2] - x[2];
        
        //distance
        dist[i] = norm(uVec + 3*i);
        
        // handle special case when the point is really close to a vertex
        if (dist[i] < eps)
        {
            exit0 = true;
            weights[i] = 1.0;
        } else {
            
            // project onto unit sphere
            uVec[3*i]   /= dist[i];
            uVec[3*i+1] /= dist[i];
            uVec[3*i+2] /= dist[i];
        }
    }
    if (exit0)
        return;

    // Now loop over all triangle to compute weights
    bool exit1 = false;
    #pragma omp parallel for
    for(int t = 0; t<n_triangles; t++)
    {
        // vertex id
        uint pid0 = triangles[3*t];
        uint pid1 = triangles[3*t+1];
        uint pid2 = triangles[3*t+2];
        
        // unit vector
        double *u0 = uVec + 3*pid0;
        double *u1 = uVec + 3*pid1;
        double *u2 = uVec + 3*pid2;
        
        // edge length
        double l0 = sqrt(distance2BetweenPoints(u1, u2));
        double l1 = sqrt(distance2BetweenPoints(u2, u0));
        double l2 = sqrt(distance2BetweenPoints(u0, u1));
        
        // angle
        double theta0 = 2.0*asin(l0/2.0);
        double theta1 = 2.0*asin(l1/2.0);
        double theta2 = 2.0*asin(l2/2.0);
        double halfSum = (theta0 + theta1 + theta2)/2.0;
        
        // special case when the point lies on the triangle
        if (PI - halfSum < eps)
        {
            for (int i=0; i < n_vertices; ++i)
            {
                weights[i] = 0.0;
            }
            
            weights[pid0] = sin(theta0)* dist[pid1]* dist[pid2];
            weights[pid1] = sin(theta1)* dist[pid2]* dist[pid0];
            weights[pid2] = sin(theta2)* dist[pid0]* dist[pid1];
            
            double sumWeight = weights[pid0] + weights[pid1] + weights[pid2];
            
            weights[pid0] /= sumWeight;
            weights[pid1] /= sumWeight;
            weights[pid2] /= sumWeight;
            
            exit1 = true;
        } else {
            
            // coefficient
            double sinHalfSum = sin(halfSum);
            double sinHalfSumSubTheta0 = sin(halfSum-theta0);
            double sinHalfSumSubTheta1 = sin(halfSum-theta1);
            double sinHalfSumSubTheta2 = sin(halfSum-theta2);
            double sinTheta0 = sin(theta0);
            double sinTheta1 = sin(theta1);
            double sinTheta2 = sin(theta2);
            
            double c0 = 2 * sinHalfSum * sinHalfSumSubTheta0 / sinTheta1 / sinTheta2 - 1;
            double c1 = 2 * sinHalfSum * sinHalfSumSubTheta1 / sinTheta2 / sinTheta0 - 1;
            double c2 = 2 * sinHalfSum * sinHalfSumSubTheta2 / sinTheta0 / sinTheta1 - 1;
            
            if (fabs(c0) > 1)
            {
                c0 = c0 > 0 ? 1 : -1;
            }
            if (fabs(c1) > 1)
            {
                c1 = c1 > 0 ? 1 : -1;
            }
            if (fabs(c2) > 1)
            {
                c2 = c2 > 0 ? 1 : -1;
            }
            
            // sign
            double det = determinant3x3(u0, u1, u2);
            
            if (fabs(det) < eps)
                continue;
            
            double detSign = det > 0 ? 1 : -1;
            double sign0 = detSign * sqrt(1 - c0*c0);
            double sign1 = detSign * sqrt(1 - c1*c1);
            double sign2 = detSign * sqrt(1 - c2*c2);
            
            // if x lies on the plane of current triangle but outside it, ignore
            // the current triangle.
            if (fabs(sign0) < eps || fabs(sign1) < eps || fabs(sign2) < eps)
                continue;
            
            // weight
            weights[pid0] += (theta0-c1*theta2-c2*theta1) / (dist[pid0]*sinTheta1*sign2);
            weights[pid1] += (theta1-c2*theta0-c0*theta2) / (dist[pid1]*sinTheta2*sign0);
            weights[pid2] += (theta2-c0*theta1-c1*theta0) / (dist[pid2]*sinTheta0*sign1);
        }
    }
    if(exit1)
        return;

    // normalize weights
    double sumWeight = 0.0;
    for (int i=0; i < n_vertices; ++i)
        sumWeight += weights[i];
    
    if (fabs(sumWeight) < eps)
        return;
    
    #pragma omp parallel for
    for (int i=0; i < n_vertices; ++i)
        weights[i] /= sumWeight;
}
