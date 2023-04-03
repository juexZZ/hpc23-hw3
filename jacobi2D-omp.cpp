// jacobian 2d omp
#include <stdio.h>
// for timing
#include "utils.h"

// map value, include the boundary
double map_value(double* u, long N, long i, long j){
    if(i<0 || i>=N || j<0 || j>=N){
        return 0;
    }
    else{
        return u[i*N+j];
    }
}

double get_residual(double* u, double f, long N, double h2){
    double res = 0;
        for(long i=0; i<N; i++){
            for(long j=0; j<N; j++){
                double lb = map_value(u, N, i-1, j-1);
                double lu = map_value(u, N, i-1, j+1);
                double rb = map_value(u, N, i+1, j-1);
                double ru = map_value(u, N, i+1, j+1);
                double diff = (-(lb+lu+rb+ru)+4*u[i*N+j])/h2 - f;
                res += diff * diff;
            }
        }
    return res;
}

//sequential version
double jacobi2d_seq(double* u,  double* u_prev, double f, long N, double h2, long max_iter, double tol=1e-6){
    double initial_res = get_residual(u_prev, f, N, h2);
    double res=0;
    for(long iter=0; iter<max_iter; iter++){
        for(long i=0; i<N; i++){
            for(long j=0; j<N; j++){
                double lb = map_value(u_prev, N, i-1, j-1);
                double lu = map_value(u_prev, N, i-1, j+1);
                double rb = map_value(u_prev, N, i+1, j-1);
                double ru = map_value(u_prev, N, i+1, j+1);
                u[i*N+j] = 0.25 * (h2*f + lb + lu + rb + ru);
            }
        }
        // get residual
        res = get_residual(u, f, N, h2);
        if (res < tol * initial_res){
            printf("Converged after %ld iterations, res= %f\n", iter, res);
            break;
        }
        u_prev = u;
    }
    return res;
}

#if defined(_OPENMP)
#include <omp.h>
// parallel version
double get_residual_omp(double* u, double f, long N, double h2){
    double res = 0;
    #pragma omp parallel for reduction(+:res) collapse(2) schedule(static)
    for(long i=0; i<N; i++){
        for(long j=0; j<N; j++){
            double lb = map_value(u, N, i-1, j-1);
            double lu = map_value(u, N, i-1, j+1);
            double rb = map_value(u, N, i+1, j-1);
            double ru = map_value(u, N, i+1, j+1);
            double diff = (-(lb+lu+rb+ru)+4*u[i*N+j])/h2 - f;
            res += diff * diff;
        }
    }
    return res;
}

double jacobi2d_omp(double* u, double* u_prev, double f, long N, double h2, long max_iter, double tol=1e-6){
    double initial_res = get_residual_omp(u_prev, f, N, h2);
    double res = 0;
    for(long iter=0; iter<max_iter; iter++){
        #pragma omp parallel for collapse(2) schedule(static)
        for(long i=0; i<N; i++){
            for(long j=0; j<N; j++){
                double lb = map_value(u_prev, N, i-1, j-1);
                double lu = map_value(u_prev, N, i-1, j+1);
                double rb = map_value(u_prev, N, i+1, j-1);
                double ru = map_value(u_prev, N, i+1, j+1);
                u[i*N+j] = 0.25 * (h2*f + lb + lu + rb + ru);
            }
        }
        // get residual
        res = get_residual_omp(u, f, N, h2);
        if (res < tol * initial_res){
            printf("Converged after %ld iterations, res= %f\n", iter, res);
            break;
        }
        u_prev = u;
    }
    return res;
}
#endif

int main(int argc, char** argv){
    long N = 1000;
    long max_iter = 10000;
    long N2 = N*N;
    int max_threads = 16;
    double f = 1.0; // f is set to 1 for simplicity, should be a vector of size N2 for general case
    double tol = 1e-6;
    double *u = (double*) malloc(N2 * sizeof(double));
    double *u_prev = (double*) malloc(N2 * sizeof(double));
    for(long i=0; i<N2; i++){
        u[i] = 0;
        u_prev[i] = 0;
    }
    double h = 1.0 / (N+1);
    double h2 = h*h;
    // sequential, if openmp is defined, run sequential with openmp timer first for reference
    #if defined(_OPENMP)
    printf("OpenMP is defined\n");
    double t = omp_get_wtime();
    double res_seq = jacobi2d_seq(u, u_prev, f, N, h2, max_iter, tol);
    double seq_time = omp_get_wtime() - t;
    printf("Sequential time: %f with residual %f\n", seq_time, res_seq);
    #else
    printf("OpenMP is not defined, only run the sequential version\n");
    Timer t;
    t.tic();
    double res = jacobi2d_seq(u, u_prev, f, N, h2, max_iter, tol);
    double seq_time = t.toc();
    printf("Sequential time: %f with residual %f\n", seq_time, res);
    #endif
    // parallel
    #if defined(_OPENMP)
    omp_set_num_threads(max_threads); // set number of threads
    double t_omp = omp_get_wtime();
    double res_omp = jacobi2d_omp(u, u_prev, f, N, h2, max_iter, tol);
    double omp_time = omp_get_wtime() - t_omp;
    printf("Parallel time: %f with residual %f\n", omp_time, res_omp);
    #endif
    return 0;
}