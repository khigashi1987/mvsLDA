/*
    mvslda.c
    Supervised Latent Dirichlet Allocation with Asymmetric Dirichlet prior, main driver.
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mvslda.h"
#include "learn.h"
#include "writer.h"
#include "feature.h"
#include "dmatrix.h"
#include "imatrix.h"
#include "util.h"

int main(int argc, char **argv){
    document *data;
    double **resp;
    FILE *pp, *tp, *likp, *hyperp; // for phi, theta, n_mz, n_zw, likelihood, hyperparameters
    FILE *ep; // for eta
    char c;
    int nlex, dlenmax;
    int ndoc;
    int nclass = NCLASS_DEFAULT;
    int nresp = NRESP_DEFAULT;
    int maxiter = MAXITER_DEFAULT;
    double alpha_init = ALPHA_DEFAULT;
    double *alpha;
    double beta = BETA_DEFAULT;
    double nu2 = NU2_DEFAULT;
    double sigma2 = SIGMA2_DEFAULT;
    unsigned long int random_seed = RANDOM_SEED_DEFAULT;
    double **phi;
    double **theta;
    double **eta;
    int **n_mz;
    int **n_zw;
    int k;
    
    while((c = getopt(argc, argv, "I:K:A:B:Y:S:h")) != -1){
        switch(c){
            case 'I': maxiter = atoi(optarg); break;
            case 'K': nclass = atoi(optarg); break;
            case 'A': alpha_init = atof(optarg); break;
            case 'B': beta = atof(optarg); break;
            case 'Y': nresp = atoi(optarg); break;
            case 'S': random_seed = atoi(optarg); break;
            case 'h': usage(); break;
            default: usage(); break;
        }
    }
    if(!(argc - optind == 3))
        usage();
    
    // open data
    if((data = feature_matrix(argv[optind], &nlex, &dlenmax, &ndoc)) == NULL){
        fprintf(stderr, "mvslda:: cannot open training data.\n");
        exit(1);
    }
    if((resp = load_dmatrix(argv[optind+1], ndoc, nresp)) == NULL){
        fprintf(stderr, "mvslda:: cannot open response data file.\n");
        exit(1);
    }
    
    // open model output
    if(((pp = fopen(strconcat(argv[optind+2], ".phi"),"w")) == NULL)
    || ((tp = fopen(strconcat(argv[optind+2], ".theta"),"w")) == NULL)
    || ((ep = fopen(strconcat(argv[optind+2], ".eta"),"w")) == NULL)
    || ((likp = fopen(strconcat(argv[optind+2], ".lik"),"w")) == NULL)
    || ((hyperp = fopen(strconcat(argv[optind+2], ".hyper"),"w")) == NULL)){
        fprintf(stderr, "mvslda:: cannot open model outputs.\n");
        exit(1);
    }
    
    // allocate parameters
    if((alpha = calloc(nclass, sizeof(double))) == NULL){
        fprintf(stderr, "corrlda:: cannot allocate alpha.\n");
        exit(1);
    }
    for(k = 0;k < nclass;k++)
        alpha[k] = alpha_init;
    
    if((phi = dmatrix(nlex, nclass)) == NULL){
        fprintf(stderr, "mvslda:: cannot allocate phi.\n");
        exit(1);
    }
    if((theta = dmatrix(ndoc, nclass)) == NULL){
        fprintf(stderr, "mvslda:: cannot allocate theta.\n");
        exit(1);
    }
    if((eta = dmatrix(nclass, nresp)) == NULL){
        fprintf(stderr, "mvslda:: cannot allocate eta.\n");
        exit(1);
    }

    // n_mz ... number of times document and topic z co-occur
    if((n_mz = imatrix(ndoc, nclass)) == NULL){
        fprintf(stderr, "mvslda:: cannot allocate n_mz.\n");
        exit(1);
    }
    // n_zw ... number of times topic and word w co-occur
    if((n_zw = imatrix(nclass, nlex)) == NULL){
        fprintf(stderr, "mvslda:: cannot allocate n_zw.\n");
        exit(1);
    }
    
    mvslda_learn(data, resp, alpha, beta, nu2, sigma2, nclass, nlex, dlenmax, nresp, maxiter, phi, theta, eta, n_mz, n_zw, likp, hyperp, random_seed);
    mvslda_write(pp, tp, ep, phi, theta, eta, nclass, nlex, nresp, ndoc);
    
    free_feature_matrix(data);
    free_dmatrix(resp, ndoc);
    free_dmatrix(phi, nlex);
    free_dmatrix(theta, ndoc);
    free_dmatrix(eta, nclass);
    free_imatrix(n_mz, ndoc);
    free_imatrix(n_zw, nclass);
    
    fclose(pp);
    fclose(tp);
    fclose(ep);
    fclose(likp);
    fclose(hyperp);
    
    exit(0);
}

void usage(void){
    printf("usage: %s [-I maxiter] [-K n_classes] [-A alpha] [-B beta] [-Y nresp] [-S random_seed] doc resp model\n", "mvslda");
    exit(0);
}
