/*///////////////////////////////// TO COMPILE /////////////////////////////////

gcc [Code File] -Wall -O3 -DHAVE_INLINE -lgsl -lgslcblas -o [Out File]

//////////////////////////////////// TO RUN ////////////////////////////////////

./[Out File] [matrix size (int)]       // number of species
             [input filename (string)] // file storing the matrix
             [seed (int)]              // seed for random number generator
             [maximize (1 or -1)]      // are we maximizing (1) or minimizing (-1)?
             [para1 (int)]             // first parameter for search algorithm
             [para2 (int)]             // second parameter for search algoritm
             [BvsC (0 or 1)]           // mutate offdiagonal (1) or diagonal (0)?
             [SearchAlg (1 or 2)]      // Genetic algoritm (1) or Hill climber (2)?

//////////////////////////////////// OUTPUT ////////////////////////////////////

A file containing a modified matrix. The name of the file is that of the input
file with a modified extension depending on the values of maximize and BvsC.

//////////////////////////////////////////////////////////////////////////////*/

// Standard Libraries
#include <stdio.h>
#include <string.h>
#include <math.h>

// Matrices, vectors, eigenvalues
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>

// Random number generation and random distributions
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Set multiplier for eigenvalues to increase difference
#define EIGMULT 10.0

// Global variables
int n;
int npairs;
gsl_vector * eval;
gsl_matrix * M;
gsl_matrix * A;
gsl_vector * diagonal;
gsl_matrix * pairs;
gsl_matrix * indices;
gsl_eigen_symm_workspace * w;


//////////////////// Functions for reading/writing matrices ////////////////////
// read a matrix stored in filename and assign to M
// M is initialized in the process
gsl_matrix * read_matrix(int n,
				char * filename){
	gsl_matrix * M = gsl_matrix_calloc(n, n);
	FILE * F;
	F = fopen(filename, "rb");
	gsl_matrix_fscanf(F, M);
	fclose(F);
	return M;
}

int initialize_structures(char * filename){
    M = gsl_matrix_calloc(n, n);
	FILE * F;
	F = fopen(filename, "rb");
	gsl_matrix_fscanf(F, M);
	fclose(F);
	
	// Now initialize all the variables needed for the eigenvalue search
    npairs = n * (n - 1) / 2;
    A = gsl_matrix_calloc(n, n);
    eval = gsl_vector_alloc(n);
    w = gsl_eigen_symm_alloc(n);
    pairs = gsl_matrix_calloc(npairs, 2);
    indices = gsl_matrix_calloc(npairs, 2);
    diagonal = gsl_vector_calloc(n);
    int i, j, k;
    k = 0;
    for (i = 0; i < n; i++){
        for (j = i; j < n; j++){
            if (j == i){
                gsl_vector_set(diagonal, j, gsl_matrix_get(M, j, j));
            }
            else {
                gsl_matrix_set(pairs, k, 0, gsl_matrix_get(M, i, j));
                gsl_matrix_set(pairs, k, 1, gsl_matrix_get(M, j, i));
                gsl_matrix_set(indices, k, 0, i);
                gsl_matrix_set(indices, k, 1, j);
                k++;
            }
        }
    }
	return 0;
}

int free_structures(){
    gsl_matrix_free(M);
    gsl_matrix_free(A);
    gsl_matrix_free(pairs);
    gsl_matrix_free(indices);
    gsl_vector_free(diagonal);
    gsl_vector_free(eval);
    gsl_eigen_symm_free(w);
	return 0;
}

// print a matrix either to a file or to a stream
int print_matrix(FILE * F,
				 gsl_matrix * M){
	int i, j;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			fprintf(F, "%f ", gsl_matrix_get(M, i, j));
		}
		fprintf(F, "\n");
	}
	return 0;
}

///////////////////////// Functions for random numbers /////////////////////////
int random_setup(gsl_rng ** r, int seed){
    const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	*r = gsl_rng_alloc(T);
	gsl_rng_set (*r, seed);
	return 0;
}

int random_free(gsl_rng ** r){
    gsl_rng_free(*r);
    return 0;
}

/////////////////////// Functions for mutating matrices ////////////////////////
// Mutate the matrix M by swapping two coefficients in the
// upper-triangular part and the corresponding lower-triangular
// coefficients as well, to keep symmetry
int mutate_off_diagonal(int n,
                        gsl_rng * r,
                        gsl_matrix * M){
    int i = gsl_rng_uniform_int(r, n);
    int j = gsl_rng_uniform_int(r, n);
    while(j == i){
        j = gsl_rng_uniform_int(r, n);
    }
    int k = gsl_rng_uniform_int(r, n);
    int l = gsl_rng_uniform_int(r, n);
    while(k == l){
        l = gsl_rng_uniform_int(r, n);
    }
    double tmp;
    tmp = gsl_matrix_get(M, i, j);
    gsl_matrix_set(M, i, j, gsl_matrix_get(M, k, l));
    gsl_matrix_set(M, k, l, tmp);

    tmp = gsl_matrix_get(M, j, i);
    gsl_matrix_set(M, j, i, gsl_matrix_get(M, l, k));
    gsl_matrix_set(M, l, k, tmp);
    return 0;
}

// Mutate by swapping two elements on the diagonal
int mutate_diagonal(int n,
                    gsl_rng * r,
                    gsl_matrix * M){
    
    int i = gsl_rng_uniform_int(r, n);
    int j = i;
    while(j == i){
        j = gsl_rng_uniform_int(r, n);
    }
    double tmp;
    tmp = gsl_matrix_get(M, i, i);
    gsl_matrix_set(M, i, i, gsl_matrix_get(M, j, j));
    gsl_matrix_set(M, j, j, tmp);
    return 0;
}

//////////////////// Functions for calculating eigenvalues /////////////////////
// Setup the environment for computing the eigenvalues
int eigenvalues_setup(int n,
					  gsl_matrix ** tmpM,
					  gsl_vector ** eval,
					  gsl_eigen_symm_workspace ** w){
	// Allocate temporary matrix for calculations
	*tmpM = gsl_matrix_calloc(n, n);
	// Allocate vector for storing eigenvalues
	*eval = gsl_vector_calloc(n);
	// Allocate workspace for eigenvalue calculation
	*w = gsl_eigen_symm_alloc(n);
	return 0;
}

// Free memory associated with eigenvalues
int eigenvalues_free(gsl_matrix ** tmpM,
					 gsl_vector ** eval,
					 gsl_eigen_symm_workspace ** w){
	gsl_matrix_free(*tmpM);
	gsl_vector_free(*eval);
	gsl_eigen_symm_free(*w);
	return 0;
}

// Find the largest eigenvalue of the symmetric matrix M
double find_max_eigen(gsl_matrix * M,
					  gsl_matrix * tmpM,
					  gsl_vector * eval,
					  gsl_eigen_symm_workspace * w){
	// Copy the matrix M, as it will be destroyed
	gsl_matrix_memcpy(tmpM, M);
	// Find the eigenvalues
	gsl_eigen_symm(tmpM, eval, w);
	// return the maximum eigenvalue
	return gsl_vector_max(eval);
}

int build_A(gsl_vector * solution, gsl_vector * diag){
    // build matrix A
    int i;
    gsl_matrix_set_zero(A);
    for (i = 0; i < npairs; i++){
        gsl_matrix_set(A, gsl_matrix_get(indices, gsl_vector_get(solution, i), 0), 
                          gsl_matrix_get(indices, gsl_vector_get(solution, i), 1),
                          gsl_matrix_get(pairs, i, 1));
        gsl_matrix_set(A, gsl_matrix_get(indices, gsl_vector_get(solution, i), 1), 
                          gsl_matrix_get(indices, gsl_vector_get(solution, i), 0),
                          gsl_matrix_get(pairs, i, 0));    
    }
    for (i = 0; i < n; i++){
        gsl_matrix_set(A, i, i, gsl_vector_get(diag, i));
    }
    return 0;
}

double find_leading_eigen(gsl_vector * solution, gsl_vector * diag){
    build_A(solution, diag);
    // always find the eigenvalues of A;
    // it is destroyed by the operation
    gsl_eigen_symm(A, eval, w);
	// return the maximum eigenvalue
	return gsl_vector_max(eval);
}

////////////////////// Functions to Run Search Algorithms //////////////////////
int GeneticAlgorithm (int n,            // number of species
                      char * filename,  // file storing the matrix
                      int seed,
                      int maximize,
                      int npop,
                      int ngen,
                      int BvsC){ // 1 for offdiagonal optimization, 0 for diagonal
    int i = 0;
    // initialize structures
    fprintf(stderr, "Initializing\n");
    i = initialize_structures(filename);

    gsl_rng * r = NULL;
    i = random_setup(&r, seed);

    if (maximize == 1) {
        fprintf(stderr, "GA to Maximize\n");
    }
    else if (maximize == -1) {
        fprintf(stderr, "GA to Minimize\n");
    }
    gsl_vector * bestsol = gsl_vector_calloc(npairs);
    gsl_vector * bestdiag = gsl_vector_calloc(n); 
    double bestfit = 0.0;

    for (i = 0; i < npairs; i++){
        gsl_vector_set(bestsol, i, i);
    }
    gsl_vector_memcpy(bestdiag, diagonal);
    
    bestfit = maximize * find_leading_eigen(bestsol, bestdiag);
    fprintf(stderr, "Original Matrix: %f\n", maximize * bestfit);

    gsl_vector * popsol[npop];
    gsl_vector * popdiag[npop];
    gsl_vector * pop2sol[npop];
    gsl_vector * pop2diag[npop];

    // initialize populations
    for (i = 0; i < npop; i++){
            popsol[i] = gsl_vector_calloc(npairs);
            popdiag[i] = gsl_vector_calloc(n);

            pop2sol[i] = gsl_vector_calloc(npairs);
            pop2diag[i] = gsl_vector_calloc(n);

            gsl_vector_memcpy(popsol[i], bestsol);
            gsl_vector_memcpy(popdiag[i], bestdiag);

            gsl_ran_shuffle(r, popsol[i]->data, n, sizeof (double));
            gsl_ran_shuffle(r, popdiag[i]->data, n, sizeof (double));
    }

    gsl_vector * fitness = gsl_vector_calloc(npop);
    int j, dad, mom, k1, k2;
    double tmp;
    for (i = 0; i < ngen; i++){
        // compute fitness and save best sol
        for (j = 0; j < npop; j++){
            tmp = maximize * find_leading_eigen(popsol[j], popdiag[j]);
            gsl_vector_set(fitness, j, tmp);
            if (tmp > bestfit){
                bestfit = tmp;
                if (BvsC == 1){
                    gsl_vector_memcpy(bestsol, popsol[j]);
                }
                else {
                    gsl_vector_memcpy(bestdiag, popdiag[j]);
                }
                fprintf(stderr, "%d %f\n", i, maximize * bestfit);
            }    
        }
        // reproduce
        for (j = 0; j < npop; j++){
            if (j < 50){
                if (BvsC == 1){
                    gsl_vector_memcpy(pop2sol[j], bestsol);
                }
                else {
                    gsl_vector_memcpy(pop2diag[j], bestdiag);
                }
            }
            else{
                dad = gsl_rng_uniform_int(r, npop);
                mom = gsl_rng_uniform_int(r, npop);
                if (gsl_vector_get(fitness, dad) > gsl_vector_get(fitness, mom)){
                    if (BvsC == 1){
                        gsl_vector_memcpy(pop2sol[j], popsol[dad]);
                    }
                    else {
                        gsl_vector_memcpy(pop2diag[j], popdiag[dad]);
                    }
                }
                else {
                    if (BvsC == 1){
                        gsl_vector_memcpy(pop2sol[j], popsol[mom]);
                    }
                    else {
                        gsl_vector_memcpy(pop2diag[j], popdiag[mom]);
                    }
                }
            }
            // mutation
            if (BvsC == 1){
                k1 = gsl_rng_uniform_int(r, npairs);
                k2 = gsl_rng_uniform_int(r, npairs);
                while (k1 == k2) {
                    k2 = gsl_rng_uniform_int(r, npairs);
                }
                tmp = gsl_vector_get(pop2sol[j], k1);
                gsl_vector_set(pop2sol[j], k1, gsl_vector_get(pop2sol[j], k2));
                gsl_vector_set(pop2sol[j], k2, tmp);
            }
            else {
                k1 = gsl_rng_uniform_int(r, n);
                k2 = gsl_rng_uniform_int(r, n);
                while (k1 == k2) {
                    k2 = gsl_rng_uniform_int(r, n);
                }
                tmp = gsl_vector_get(pop2diag[j], k1);
                gsl_vector_set(pop2diag[j], k1, gsl_vector_get(pop2diag[j], k2));
                gsl_vector_set(pop2diag[j], k2, tmp);
            }
        }
        // copy over
        for (j = 0; j < npop; j++){
            if (BvsC == 1){
                gsl_vector_memcpy(popsol[j], pop2sol[j]);
            }
            else {
                gsl_vector_memcpy(popdiag[j], pop2diag[j]);
            }
        }
        
    }

    FILE * F;
    char OutFileName[1000];
    if (BvsC == 1) {
        if (maximize == 1){
            sprintf(OutFileName,"%s.Bmax--%d", filename, seed);
        }
        if (maximize == -1){
            sprintf(OutFileName,"%s.Bmin--%d", filename, seed);
        }
    } else {
        if (maximize == 1){
            sprintf(OutFileName,"%s.max--%d", filename, seed);
        }
        if (maximize == -1){
            sprintf(OutFileName,"%s.min--%d", filename, seed);
        }
    }
    F = fopen(OutFileName, "wb");
    build_A(bestsol, bestdiag);
    print_matrix(F, A);
    fclose(F); 
    
    gsl_vector_free(bestsol);

    for (i = 0; i < npop; i++){
        gsl_vector_free(popsol[i]);
        gsl_vector_free(pop2sol[i]);
    }
    gsl_vector_free(fitness);
    
    free_structures();
    random_free(&r);
    return 0;
}

// Hill-climber with multiple sampling at each step to find the matrix structure
// that maximizes or minimizes the leading eigenvalue
int HillClimb (int n,           // number of species
               char * filename, // file storing the matrix
               int num_steps,   // number of steps without improvement before quitting
               int num_try,     // number of solutions to try at each step
               int seed,        // random seed
               int maximize,    // 1 to maximize, -1 to minimize
               int BvsC){       // mutate diagonal (0) or off-diagonal (1)
    int i = 0;
    double fit = 0.0;
    
    // read in the matrix
    gsl_matrix * M = NULL;
    M = read_matrix(n, filename);
    
    // initialize the random number generator
    gsl_rng * r = NULL;
    i = random_setup(&r, seed);
    
    // setup eigenvalues
    gsl_matrix * tmpM = NULL; // temporary matrix for eigenvalue calculation
    gsl_vector * eval = NULL; // vector for storing eigenvalues
    gsl_eigen_symm_workspace * w = NULL; // workspace for eigenvalues
    i = eigenvalues_setup(n, &tmpM, &eval, &w);
	
    // SPECIFIC TO HC
    fprintf(stderr, "starting HC of leading eigenvalue with %d steps,
                     %d tries per step\n", num_steps, num_try);
    
    gsl_matrix * M1 = gsl_matrix_calloc(n, n);
    gsl_matrix * M2 = gsl_matrix_calloc(n, n);
    

    fit = find_max_eigen(M, tmpM, eval, w);
   
    double fit1 = fit;
    double fit2 = fit;
    
    int how_many_steps = 0;
    int counter = 0;
    while(how_many_steps < num_steps){
        counter++;
        how_many_steps++;
        gsl_matrix_memcpy(M1, M);
        fit1 = fit;
        // now mutate several times and save the best
        for (i = 0; i < num_try; i++){
            gsl_matrix_memcpy(M2, M);
            if (BvsC == 1){
                mutate_off_diagonal(n, r, M2);
            }
            else {
                mutate_diagonal(n, r, M2);
            }
            fit2 = find_max_eigen(M2, tmpM, eval, w);
            if ((fit1 - fit2) * maximize < 0){
                // we have a better solution
                fit1 = fit2;
                gsl_matrix_memcpy(M1, M2);
            }
        }
        if ((fit - fit1) * maximize < 0){
            // accept the move
            fit = fit1;
            gsl_matrix_memcpy(M, M1);
            how_many_steps = 0;
            fprintf(stderr, "step %d -- new best solution: %.6f\n", counter, fit);
        }
        else if ((counter % 1000) == 0){
            fprintf(stderr, "step %d\n", counter);
        }
    }
    
    gsl_matrix_free(M1);
    gsl_matrix_free(M2);

    // END SPECIFIC TO HC    

    // Stuff for output
    char OutFileName[1000];
    char FileNameRoot[1000];
    char *dotBm;
    strcpy(FileNameRoot, filename);
    dotBm = strrchr(FileNameRoot, '.');
    *dotBm = '\0';
    // Save the results
    if (maximize == 1) {
        if (BvsC == 1) {
            sprintf(OutFileName,"%s.Bmax--%d", FileNameRoot, seed);
        }
        else{
            sprintf(OutFileName,"%s.max--%d", FileNameRoot, seed);
        }

    } else if (maximize == -1) {
        if (BvsC == 1) {
            sprintf(OutFileName,"%s.Bmin--%d", FileNameRoot, seed);
        }
        else{
            sprintf(OutFileName,"%s.min--%d", FileNameRoot, seed);
        }
    }

	int l, m;
    FILE * F = fopen(OutFileName,"w");
	for (l = 0; l < n; l++){
		for (m = 0; m < n; m++){
			fprintf(F, "%f ", gsl_matrix_get(M, l, m));
		}
		fprintf(F, "\n");
	}
    fclose(F);
    // free matrix M
    gsl_matrix_free(M);
    // free eigenvalue-related variables
    i = eigenvalues_free(&tmpM, &eval, &w);
    // free random number generator
    i = random_free(&r);
    return 0;
}


int main (int argc, char *argv[]){
    n = atoi(argv[1]);             // number of species
    char * filename = argv[2];     // file storing the matrix
    int seed = atoi(argv[3]);      // seed for random number generator
    int maximize = atoi(argv[4]);  // are we maximizing (1) or minimizing (-1)?
    int para1 = atoi(argv[5]);     // first parameter for search algorithm
    int para2 = atoi(argv[6]);     // second parameter for search algoritm
    int BvsC = atoi(argv[7]);      // what to mutate: 1 for offdiagonal, 0 for diagonal
    int SearchAlg = atoi(argv[8]); // Genetic algoritm (1) or Hill climber (2)?
    
    if (SearchAlg == 1) {
        GeneticAlgorithm (n, filename, seed, maximize, para1, para2, BvsC);
    }
    else if (SearchAlg == 2) {
        HillClimb (n, filename, para1, para2, seed, maximize, BvsC);
    }
    return(0);
}
