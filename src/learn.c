/*
    learn.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multifit.h>
#include "learn.h"
#include "feature.h"
#include "imatrix.h"
#include "dmatrix.h"
#include "util.h"
#include "likelihood.h"
#include "hyper.h"

void mvslda_learn(document *data, double **resp, double *alpha, double beta, double nu2, double sigma2, 
    int nclass, int nlex, int dlenmax, int nresp, int maxiter, double **phi, double **theta, double **eta, 
    int **n_mz, int **n_zw, FILE *likp, FILE *hyperp, unsigned long int random_seed){
    document *dp;
    int ndocs;
    int *n_m;
    int *n_z;
    int ***topics;
    int word_index;
    int word_num;
    double sum_alpha;
    double *left;
    double *center;
    double *log_right;
    double *p_z;
    double *log_p_z;
    double *cum_sum_p_z;
    double log_Z, sum_p_r, sum_empirical_z;
    double temp_prediction;
    double lik;
    double **temp_phi;
    double **temp_theta;
    double **empirical_z;
    int z;
    int it;
    int m, w, t, i, j, k;
    const gsl_rng_type *T;
    gsl_rng *r;
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    // count data length
    for(dp = data, ndocs = 0;(dp->len) != -1;dp++, ndocs++)
        ;
    
    // initialize buffers
    if((n_m = calloc(ndocs,sizeof(int))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate n_m.\n");
        return;
    }
    if((n_z = calloc(nclass,sizeof(int))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate n_z.\n");
        return;
    }
    if((left = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate left.\n");
        return;
    }
    if((center = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate center.\n");
        return;
    }
    if((log_right = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate log_cright.\n");
        return;
    }
    if((p_z = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate p_z.\n");
        return;
    }
    if((log_p_z = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate log_p_z.\n");
        return;
    }
    if((cum_sum_p_z = calloc((nclass+1),sizeof(double))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate cum_sum_p_z.\n");
        return;
    }
    if((topics = calloc(ndocs,sizeof(int **))) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate topics.\n");
        return;
    }
    if((empirical_z = dmatrix(ndocs, nclass)) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate empirical_z.\n");
        return;
    }
    if((temp_phi = dmatrix(nlex, nclass)) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate temp_phi.\n");
        exit(1);
    }
    if((temp_theta = dmatrix(ndocs, nclass)) == NULL){
        fprintf(stderr,"mvslda_learn:: cannot allocate temp_theta.\n");
        exit(1);
    }
    
    printf("Number of documents          = %d\n",ndocs);
    printf("Number of unique words       = %d\n",nlex);
    printf("Number of latent classes     = %d\n",nclass);
    printf("Number of responses          = %d\n",nresp);
    printf("Number of iteration          = %d\n",maxiter);
    
    // choose an arbitrary topic as first topic for word
    gsl_rng_set(r, random_seed);
    for(dp = data, m = 0;(dp->len) != -1;dp++, m++){
        if((topics[m] = calloc((dp->len), sizeof(int *))) == NULL){
            fprintf(stderr,"mvslda_learn:: cannot allocate topics[m].\n");
            return;
        }
        for(w = 0;w < (dp->len);w++){
            if((topics[m][w] = calloc((dp->cnt[w]), sizeof(int))) == NULL){
                fprintf(stderr,"mvslda_learn:: cannot allocate topics[m][w].\n");
                return;
            }
            word_index = dp->id[w];
            word_num = dp->cnt[w];
            for(i = 0;i < word_num;i++){
                z = (int)gsl_rng_uniform_int(r, nclass);
                n_mz[m][z] += 1;
                n_m[m] += 1;
                n_zw[z][word_index] += 1;
                n_z[z] += 1;
                topics[m][w][i] = z;
            }
        }
    }
    // initialize eta ([nclass, nresp])
    for(k = 0;k < nclass;k++)
        for(t = 0;t < nresp;t++)
            eta[k][t] = 0.0;
    
    // learning main
    for(it = 0;it < maxiter;it++){
        printf("iteration %2d/%3d..\n", it + 1, maxiter);
        fflush(stdout);
        sum_alpha = 0.0;
        for(k = 0;k < nclass;k++)
            sum_alpha += alpha[k];
        for (dp = data, m = 0; (dp->len) != -1; dp++, m++){
            // for words
            for(w = 0;w < (dp->len);w++){
                word_index = dp->id[w];
                word_num = dp->cnt[w];
                for(i = 0;i < word_num;i++){
                    z = topics[m][w][i];
                    n_mz[m][z] -= 1;
                    n_m[m] -= 1;
                    n_zw[z][word_index] -= 1;
                    n_z[z] -= 1;
                    
                    // compute conditional distribution log_p_z
                    // log_p_z left ... theta term
                    for(k = 0;k < nclass;k++){
                        left[k] = (double)n_mz[m][k] + alpha[k];
                        left[k] /= ((double)n_m[m] + sum_alpha);
                    }
                    // log_p_z center ... phi term
                    for(k = 0;k < nclass;k++){
                        center[k] = (double)n_zw[k][word_index] + beta;
                        center[k] /= ((double)n_z[k] + (double)nlex * beta);
                    }
                    // temporal log_p_z (left and center)
                    for(k = 0; k < nclass;k++){
                        log_p_z[k] = log(left[k]) + log(center[k]);
                    }
                    // p_z right ... eta term
                    sum_empirical_z = 0.0;
                    for(k = 0; k < nclass;k++){
                        empirical_z[m][k] = (double)n_mz[m][k];
                        sum_empirical_z += (double)n_mz[m][k];
                    }
                    for(k = 0; k < nclass;k++){
                        empirical_z[m][k] = empirical_z[m][k] / sum_empirical_z;
                    }
                    for(t = 0; t < nresp; t++){
                        temp_prediction = 0.0;
                        for(k = 0;k < nclass;k++){
                            temp_prediction += eta[k][t] * empirical_z[m][k]; // dot(eta, z_d)
                        }
                        for(k = 0;k < nclass;k++){
                            log_right[k] = 1.0;
                            log_right[k] *= 1.0 / (2 * sigma2);
                            log_right[k] *= (eta[k][t] / (double)n_m[m]);
                            log_right[k] *= (2 * (resp[m][t] - temp_prediction) - (eta[k][t] / (double)n_m[m]));
                            log_p_z[k] += log_right[k];
                        }
                    }
                    // conditional distribution log_p_z
                    // log_Z = logsumexp(logP_k1 + logP_k2 + ... logP_kK)
                    log_Z = logsumexp(log_p_z, nclass);
                    for(k = 0;k < nclass;k++){
                        p_z[k] = exp(log_p_z[k] - log_Z); // normalize to obtain probabilities
                    }
                    // random sampling from p_z
                    z = sampling_multinomial(r, p_z, cum_sum_p_z, nclass);
                    // update buffers
                    n_mz[m][z] += 1;
                    n_m[m] += 1;
                    n_zw[z][word_index] += 1;
                    n_z[z] += 1;
                    topics[m][w][i] = z;
                }
            }
        }

        // for eta update
        //least squares for dot(Z, eta) = resp
        for(m = 0; m < ndocs; m++){
            sum_empirical_z = 0.0;
            for(k = 0; k < nclass;k++){
                empirical_z[m][k] = (double)n_mz[m][k];
                sum_empirical_z += (double)n_mz[m][k];
            }
            for(k = 0; k < nclass;k++){
                empirical_z[m][k] = empirical_z[m][k] / sum_empirical_z;
            }
        }
        for(t = 0; t < nresp; t++){
            double chisq;
            gsl_matrix *Z, *cov;
            gsl_vector *y, *c;
            gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (ndocs, nclass);
            Z = gsl_matrix_alloc(ndocs, nclass);
            cov = gsl_matrix_alloc(nclass, nclass);
            y = gsl_vector_alloc(ndocs);
            c = gsl_vector_alloc(nclass);
            for(m = 0;m < ndocs;m++){
                gsl_vector_set(y, m, resp[m][t]);
                for(k = 0;k < nclass;k++){
                    gsl_matrix_set(Z, m, k, empirical_z[m][k]);
                }
            }
            gsl_multifit_linear(Z, y, c, cov, &chisq, work);
            for(k = 0;k < nclass;k++){
                eta[k][t] = gsl_vector_get(c, k);
            }
            gsl_multifit_linear_free(work);
            gsl_matrix_free(Z);
            gsl_matrix_free(cov);
            gsl_vector_free(y);
            gsl_vector_free(c);
        }
        
        // update hyperparameters.
        update_alpha(alpha, n_m, n_mz, ndocs, nclass);
        beta = update_beta(beta, n_z, n_zw, nclass, nlex);
        
        // compute likelihood.
        lik = loglikelihood(n_mz, n_zw, n_m, nclass, nlex, ndocs, nresp, resp, alpha, beta, eta, empirical_z, nu2, sigma2);
        printf("\tlikelihood ... %.8f\n",lik);
        printf("\talpha = \n\t");
        for(k = 0;k < nclass;k++)
            printf("%.8f ",alpha[k]);
        printf("\n\tbeta ... %.2f\n",beta);
        fprintf(likp,"%.8f\n",lik);
        for(k = 0;k < nclass;k++)
            fprintf(hyperp,"%.8f,",alpha[k]);
        fprintf(hyperp,"%.8f\n",beta);
    }
    
    // compute matrix phi ([nlex, nclass] matrix)
    for(w = 0;w < nlex;w++)
        for(k = 0;k < nclass;k++)
            temp_phi[w][k] = (double)n_zw[k][w] + beta;
    normalize_matrix_col(phi, temp_phi, nlex, nclass);
    
    // compute matrix theta ([ndocs, nclass])
    for(m = 0;m < ndocs;m++)
        for(k = 0;k < nclass;k++)
            temp_theta[m][k] = (double)n_mz[m][k] + alpha[k];
    normalize_matrix_row(theta, temp_theta, ndocs, nclass);
    
    free(n_m);
    free(n_z);
    free(left);
    free(center);
    free(log_right);
    free(p_z);
    free(log_p_z);
    free(cum_sum_p_z);
    
    for(dp = data, m = 0;(dp->len) != -1;dp++, m++){
        for(w = 0;w < (dp->len);w++){
            free(topics[m][w]);
        }
        free(topics[m]);
    }
    free(topics);
    free_dmatrix(temp_phi, nlex);
    free_dmatrix(temp_theta, ndocs);
    free_dmatrix(empirical_z, ndocs);
    
    return;
}

int sampling_multinomial(gsl_rng *r, double *p, double *cum_sum_p, int len_p){
    int k, z;
    double sampling;
    
    cum_sum_p[0] = 0.0;
    for(k = 0;k < len_p;k++){
        cum_sum_p[k+1] = cum_sum_p[k] + p[k];
    }
    sampling = gsl_rng_uniform(r);
    for(k = 0;k < len_p;k++){
        if((sampling >= cum_sum_p[k]) && (sampling < cum_sum_p[k+1])){
            z = k;
            break;
        }
    }
    return z;
}
