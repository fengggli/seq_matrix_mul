#include "seq_matrix_mul.h"


/*
 * generate a random matrix with double value -1~1
 * init_type
 *  1: init with random number (-1~1)
 *  0: init with 0
 */
int init_matrix(double ** m, int n, int init_type){

    

    int i, j;
    double tmp;
    srand(time(NULL));
    // allocate space
    *m = (double*)malloc(sizeof(double)*n*n);
    if(*m ==NULL){
        perror("allocate error");
        exit(-1);
    }

    // initiate with random number
    for( i = 0; i< n; i++){
        for(j = 0; j < n ; j++){

            if(init_type == 1){
                tmp = (double)rand()/RAND_MAX*2.0-1.0;
            }
            else{
                tmp = 0;
            }
            (*m)[i*n + j] = tmp;
        }
    }
    return 0;
}

void free_matrix(double *m){
    if(m != NULL){
        free(m);
    }
}

void print_matrix(double *m, int n){
    int i, j;
    for(i = 0 ; i < n ; i++){
        for(j = 0; j< n; j++){
            printf("%4.3lf ", m[i*n +j]);
        }
        printf("\n");
    }
}


// return 0 if two matrix are equal
/*
int compare_matrix(double *A, double *B, int n){
    int i, j;
    int equal =0;
    for(i = 0; i < n; i++){
        for(j = 0; j < n ; j++){
            if(A[i*n +j]!= B[i*n+j]){
                return -1;
            }
        }
    }
    return 0;
}
*/
void blas_cal(double *C, double* A, double *B, int n){
    int m, k;
    m = n;
    k = n;

    double  alpha =1;
    double  beta = 0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                            m, n, k, alpha, A, k, B, n, beta, C, n);
}




void seq_cal_ijk(double *C, double *A, double *B, int n){
    int i, j, k;
    double tmp;
    for(i = 0; i< n; i++){
        for(j = 0; j < n;j++){
            tmp = 0;
            for(k = 0; k<n;k++){
                tmp += A[i*n +k]*B[k*n+ j]; // cij +=A[i][k]+ B[k][j]
            }
            C[i*n +j] += tmp;
        }
    }
}

void seq_cal_ikj(double *C, double *A, double *B, int n){
    int i, j, k;
    //double cij;
    for(i = 0; i< n; i++){
        for(k = 0; k< n;k++){
            for(j = 0; j < n;j++){
                /*
                cij = C[i*n + j];
                cij += A[i*n +k]* B[k*n+ j]; // cij +=A[i][k]+ B[k][j]
                C[i*n +j] = cij;
                */

                C[i*n +j] += A[i*n +k]* B[k*n+ j];
            }
        }
    }
}

void seq_cal_jik(double *C, double *A, double *B, int n){
    int i, j, k;
    double tmp;

    for(j = 0; j < n;j++){
        for(i = 0; i< n; i++){
            tmp = 0;
            for(k = 0; k<n;k++){
                tmp += A[i*n +k]*B[k*n+ j]; // cij +=A[i][k]+ B[k][j]
            }
            C[i*n +j] += tmp;
        }
    }
}

void seq_cal_jki(double *C, double *A, double *B, int n){
    int i, j, k;
    //double cij;

    for(j = 0; j < n;j++){
        for(k = 0; k<n;k++){
            for(i = 0; i< n; i++){
                /*
                cij = C[i*n + j];
                cij += A[i*n +k]*B[k*n+ j]; // cij +=A[i][k]+ B[k][j]
                C[i*n +j] = cij;
                */
                C[i*n +j] += A[i*n +k]*B[k*n+ j];
            }
        }
    }
}

void seq_cal_kij(double *C, double *A, double *B, int n){
    int i, j, k;

    for(k = 0; k< n;k++){
        for(i = 0; i< n; i++){
            for(j = 0; j < n;j++){
                /*
                cij = C[i*n + j];
                cij += A[i*n +k]*B[k*n+ j]; // cij +=A[i][k]+ B[k][j]
                C[i*n +j] += cij;
                */
                C[i*n +j] += A[i*n +k]*B[k*n+ j];
            }
        }
    }
}

void seq_cal_kji(double *C, double *A, double *B, int n){
    int i, j, k;

    for(k = 0; k<n;k++){
        for(j = 0; j < n;j++){
            for(i = 0; i< n; i++){
                /*
                cij = C[i*n + j];
                cij += A[i*n +k]*B[k*n+ j]; // cij +=A[i][k]+ B[k][j]
                C[i*n +j] = cij;
                */

                C[i*n +j] += A[i*n +k]*B[k*n+ j];
            }
        }
    }
}

// verify results usig blas
int verify(double *C, double *C_blas, int n){
    int i, j;

    for(i = 0; i < n; i++){
        for(j = 0; j< n; j++){
            if(abs(C[i*n +j] -  C_blas[i*n+j]) > 0.000001){
                printf("%3.12lf != %3.12lf\n", C[i*n +j], C_blas[i*n +j]);
                return -1;
            }
        }
    }
    return 0;
}


/*
 * main routine
 * input:
 *  N:matrix size
 * return:
 *  0 if success
 */
int main(int argc, char * argv[]){
    int N, i, j, retval;
    double *C, *A, *B, *C_blas;

    double btime, etime;

    double total_time[6] = {0};

    // l1 data cache access and L1 data cache miss
    long long count_l1_a[6] = {0};
    long long count_l1_m[6] = {0};

    int Events[NUM_EVENTS]={PAPI_L1_DCA, PAPI_L1_DCM,PAPI_TLB_DM};
    int EventSet = PAPI_NULL;
    long long values[NUM_EVENTS];


    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) 
        exit(1);

    // create event set
    retval = PAPI_create_eventset (&EventSet);
    /* Add events to the eventset */
    retval = PAPI_add_events (EventSet, Events, NUM_EVENTS);


    char *method[6] = {
        "ijk",
        "ikj",
        "jik",
        "jki",
        "kij",
        "kji"
    };
    // matrix size
    if(argc != 2){
        printf("format: aprun seq_matrix_mul N, now exit\n");
        exit(-1);
    }

    N = atoi(argv[1]);
    if(N >0){
        printf("matrix all have dimension %d x %d\n", N, N);
    }

    // there will be six types of loops
    //int cal_types = 0;
    void (*seq_cal[6])(double*, double*, double*, int);
    seq_cal[0] = &seq_cal_ijk;
    seq_cal[1] = &seq_cal_ikj;
    seq_cal[2] = &seq_cal_jik;
    seq_cal[3] = &seq_cal_jki;
    seq_cal[4] = &seq_cal_kij;
    seq_cal[5] = &seq_cal_kji;


   
    // generate three input matrix
    C = NULL;
    B = NULL;
    A = NULL;
    C_blas =NULL;

    for(j = 0 ; j < NUM_EXP ; j ++){
        printf("****************************************************\n");
        printf("start No.%d experiment \n", j);

        if(init_matrix(&B, N, 1) == 0){
            printf("\trandom matrix B is generated \n");
        }
        if(init_matrix(&A, N, 1) == 0){
            printf("\trandom matrix A is generated \n");
        }

        // start blas_version

        if(init_matrix(&C_blas, N, 0) == 0){
            printf("\tstart blas  version...\n");
        }

        blas_cal(C_blas, A, B, N);
        //print_matrix(C_blas, N);
        //printf("\t\tElapsed Time: %16.9f second\n", etime - btime);
        
        //start sequtial execution
        for(i = 0; i<6; i++){

            printf("\tstart method %s...", method[i]);

            init_matrix(&C, N, 0);

            // init counter
            btime = get_cur_time(); 
            retval = PAPI_start (EventSet);

            (*seq_cal[i])(C, A, B, N);

            // get counter
            retval = PAPI_stop (EventSet, values);
            etime = get_cur_time(); 
            //print_matrix(C, N);

            // sum to get avg
            total_time[i] += etime-btime;
            count_l1_a[i] += values[0];
            count_l1_m[i] += values[1];


            // use mkl intel scientic library to verify
            // use cray libsci(blas to verify)
            if(verify(C, C_blas, N) == 0){

//#ifdef VERBOSE
                printf("\t\tElapsed Time: %16.9f second, %ld miss in %ld acess, mr=%lf\t, tlb = %ld\t", etime - btime, values[1],values[0], (values[1]+0.0)/values[0], values[2]);
                printf("correct\n");
//#endif
            }else{
                free_matrix(C);
                free_matrix(A);
                free_matrix(B);
                free_matrix(C_blas);

                printf("calculation not correct, now exit\n");
                return -1;
            }
            free_matrix(C);
        }
        free_matrix(A);
        free_matrix(B);
        free_matrix(C_blas);

        printf("\tall matrix are freed\n"); 
    }
    // print average
    printf("the average cacluation time of %d runs:\n", NUM_EXP);

    fprintf(stderr, "%d", N);

    for(i = 0; i <  6; i ++){
        double avg_time = total_time[i]/NUM_EXP;
        double avg_gflops = 2*N*N*N*1E-9/avg_time;
        double avg_l1_a = count_l1_a[i]/NUM_EXP;
        double avg_l1_m = count_l1_m[i]/NUM_EXP;
        double avg_l1_mr = avg_l1_m/avg_l1_a;

        printf("method %s:\t %lf s, %lf gflops, L1 data access %lf , L1 data miss %lf, miss rate %lf\n", method[i],avg_time, avg_gflops, avg_l1_a, avg_l1_m, avg_l1_mr);
        fprintf(stderr, "\t%lf\t%lf\t%lf\t%lf", avg_gflops, avg_l1_a, avg_l1_m, avg_l1_mr);
    }
    fprintf(stderr, "\n");
    return 0;
}


    


   

