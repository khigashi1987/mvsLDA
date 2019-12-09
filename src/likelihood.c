/*
    likelihood.c
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double log_multi_beta_vector(int *vec, int length, double padding){
    double *new_vec;
    int i;
    double left;
    double right;
    
    if((new_vec = calloc(length,sizeof(double))) == NULL){
        fprintf(stderr,"lda_likelihood:: cannot allocate new_vec.\n");
        exit(1);
    }
    
    for(i = 0;i < length;i++){
        new_vec[i] = (double)vec[i] + padding;
    }
    
    left = 0.0;
    for(i = 0;i < length;i++){
        left += lgamma(new_vec[i]);
    }
    right = 0.0;
    for(i = 0;i < length;i++){
        right += new_vec[i];
    }
    right = lgamma(right);
    
    return (left - right);
}

static double log_multi_beta_scalar(double val, int K){
    return ((double)K * lgamma(val) - lgamma((double)K * val));
}

double loglikelihood(int **n_mz, int **n_zw, int *n_m, 
    int nclass, int nlex, int ndocs, int nresp, double **resp,
    double *alpha, double beta, double **eta, double **empirical_z, double nu2, double sigma2){
    double lik = 0.0;
    double sum_alpha = 0.0;
    double left_const = 0.0;
    int m, k, t;
    double prediction;
    
    // P(W | Z, beta)
    for(k = 0;k < nclass;k++){
        lik += log_multi_beta_vector(n_zw[k], nlex, beta);
        lik -= log_multi_beta_scalar(beta, nlex);
    }
    
    // P(Z | alpha)
    for(k = 0;k < nclass;k++)
        sum_alpha += alpha[k];
    left_const += lgamma(sum_alpha);
    for(k = 0;k < nclass;k++)
        left_const -= lgamma(alpha[k]);
    for(m = 0;m < ndocs;m++){
        lik += left_const;
        for(k = 0;k < nclass;k++)
            lik += lgamma((double)n_mz[m][k] + alpha[k]);
        lik -= lgamma((double)n_m[m] + sum_alpha);
    }

    // P(eta | mu=0, nu2)
    for(t = 0;t < nresp;t++){
        for(k = 0;k < nclass; k++){
            lik += -0.5 * (log(2.0 * M_PI * nu2) + (eta[k][t]*eta[k][t] / nu2));
        }
    }

    // P(Y | eta, Z, sigma2)
    for(t = 0;t < nresp;t++){
        for(m = 0;m < ndocs;m++){
            prediction = 0.0;
            for(k = 0;k < nclass;k++){
                prediction += empirical_z[m][k] * eta[k][t];
            }
            lik += -0.5 * (log(2.0 * M_PI * sigma2) + (pow((resp[m][t] - prediction), 2.0) / sigma2));
        }
    }
    
    return lik;
}