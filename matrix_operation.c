#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define MIN 0.0
#define MAX 10.0
#define TOT_RUN 10

double w_end, w_start;
double stime_mpi, etime_mpi;


//global variable for the speedup and the efficency (glocal because I use it in two different function)
double avg_time_sequential_cs; //for the checksym

// all the function useful for the code
void initializeMatrix(float **mat, int _n);
void print(float **mat, int _n);

int checkSym(float **mat, int _n) ; 
void matTranspose(float **mat, float **transpose, int _n);


int  checkSymMPI(float **mat, int _n, int rank, int size);
void matTransposeMPI(float **mat, float **transpose, int _n, int rank, int size);

int checksquare(int num); 
void checktraspose(float **SEQ, float **B, int _n);

int main(int argc, char** argv){

  // Initialize MPI and get the total processes and the current process rank
  MPI_Init(&argc, &argv);
  
  int myrank, size;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  //Initialize the random number generator with the current time to have different random sequences on each program execution 
  srand(time(NULL)); 
    
  //take a value (given by the user) as input
  int n = atoi(argv[1]); 
   
  double total_time_sequential = 0.0;
  double time_sequential = 0.0; 
  double avg_time_sequential = 0.0;
  double total_time_mpi = 0.0;
  double time_mpi = 0.0; 
  double avg_time_mpi = 0.0;
  
  double avg_speedup = 0.0;
  double avg_efficiency = 0.0;
  double avg_scalability = 0.0;
  double bandwidth = 0.0;   
  
  if (myrank == 0) {
  if (checksquare(n) == 0) {
    printf("the dimention of the matrix is not a power of two\n"); 
    MPI_Finalize();
    return 1; 
  } }
  
  if ((n % size) != 0) {
        if (myrank == 0) {
            printf("the dimention of the matrix is not divisible by the number of processes\n");
        }
        MPI_Finalize();
        return 1;
    }
  
  // Allocate memory for matrix M
  float **M = (float **) malloc(n * sizeof(float *)); 
  float **TSEQ = (float **) malloc(n * sizeof(float *));
  float **TMPI = (float **) malloc(n * sizeof(float *));
  
  for (int i = 0; i < n; i++) {
        M[i] = (float *) malloc(n * sizeof(float));
        TSEQ[i] = (float *) malloc(n * sizeof(float));
        TMPI[i] = (float *) malloc(n * sizeof(float));
  } 
    
  //check for allocation success
  if ( M == NULL) {
    printf (" Memory allocation failed \n") ;
    MPI_Finalize();
    return 1; 
  } 
  
  if (myrank == 0) { initializeMatrix(M, n); }
   
  //print(M, n);

//sequential__
  
  //check for the allocation success
  if ( TSEQ == NULL) {
    printf (" Memory allocation failed \n") ;
    MPI_Finalize();
    return 1;
  }
  
  for (int i = 0; i < n; ++i) {
        MPI_Bcast(M[i], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }
  
  if (myrank == 0) {
    if (checkSym(M, n) == 1 ){   
      printf("the matrix M is symmetric \n");   
    } else { 
     
       for (int run = 0; run < TOT_RUN; run++) {
       
        double wstart = MPI_Wtime(); 
        matTranspose(M, TSEQ, n); 
        double wend = MPI_Wtime();
        
        time_sequential = wend - wstart;
        total_time_sequential += time_sequential;
        
      }
        avg_time_sequential = total_time_sequential / TOT_RUN;
        bandwidth = ((2 * (n *n) * sizeof(float))/avg_time_sequential)/ 1e9; 
        
        if (myrank == 0) {
        printf("\n____sequential____\n - transposition routine check: %.4g seconds\n - bandwidth: %.2f", avg_time_sequential, bandwidth);} 
        printf("\n"); 
              
        //checktraspose(M, TSEQ, n);    
      } 
  }
 
//Message Passing Interface__
 
  if (TMPI == NULL) {
    if (myrank == 0){
      printf (" Memory allocation failed \n");
      MPI_Finalize();
      return 1;
  }}
  
  if (checkSymMPI(M, n, myrank, size) == 1) {
  
    if (myrank == 0){
        printf("the matrix M is symmetric \n"); } 
  } else {
 
    for (int i = 0; i < n; ++i) {
          MPI_Bcast(M[i], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    
    for (int run = 0; run < TOT_RUN; run++) {
      
      stime_mpi = MPI_Wtime();
      matTransposeMPI(M, TMPI, n, myrank, size); 
      etime_mpi = MPI_Wtime();
      
      time_mpi = etime_mpi - stime_mpi;
      total_time_mpi += time_mpi;
    
  }
    avg_time_mpi = total_time_mpi / TOT_RUN; 
    avg_speedup = avg_time_sequential/ avg_time_mpi;
    avg_efficiency = avg_speedup / size;
    bandwidth = ((2 * (n *n) * sizeof(float))/avg_time_mpi)/ 1e9; 
    avg_scalability = (avg_time_sequential * n)/(avg_time_mpi * n * size); 
    
  if (myrank == 0) {
    printf("\n____Message Passing Interface (MPI)____\n - transposition routine check: %.4g seconds\n - speedup: %.2f \n - efficiency: %.2f%% \n - bandwidth: %.2f \n - scalability: %.2f\n ", avg_time_mpi,  avg_speedup, avg_efficiency * 100, bandwidth, avg_scalability);}
  }

  //________________________Free the allocated memory______________________________
  for (int i = 0; i < n; i++) {
    free(M[i]);
    free(TSEQ[i]);
    free(TMPI[i]);
  }
  free(M); 
  free(TSEQ);
  free(TMPI);
  
  MPI_Finalize();
 
  return 0;
}

//check if dimention of the matrix is a power of two
int checksquare(int num) {

  if (num == 0) { return 0;}
  
  while (num%2 == 0) { 
  num= num/2;} 
  return (num==1);
   
}

//check if the transpose of the sequential code is equal to the other 
void checktranspose(float **SEQ, float **B, int _n) {

  int check = 1; 
  
  for (int i = 0; i < _n; i++) {
 	  for(int j = 0; j < _n; j++ ) {
  		if (SEQ[i][j] != B[i][j]){ check = 0;} 
  	}
  }
  
  if (check == 1) { printf("\nequal\n"); 
  } else { printf("\nnot equal\n"); }

}


//inizializate the matrix
void initializeMatrix(float **mat, int _n) {
  
  for (int i = 0; i < _n; i++) {
  	for(int j = 0; j < _n; j++ ) {
     //insert randomic number in the matrix
  		mat[i][j] = MIN + (float)rand()/(float) (RAND_MAX / (MAX-MIN)); 
  	}
  }
  
}

//print the matrix
void print(float **mat, int _n) {

  for (int i = 0; i < _n; i++) {
  	for(int j = 0; j < _n; j++ ) {
      //.1 to specify only one digit after the decimal point 
 		  printf("%.1f  ", mat[i][j]);
  	}
  	printf("\n"); 
  }
  
}

//control if the matrix is symmetric 
int checkSym(float **mat, int _n) {

  double total_time_sequential_cs= 0.0;
  double time_sequential_cs = 0.0; 
  int is_symmetric= 1; 
  
  for (int run = 0; run < TOT_RUN; run++) {
  
  //Start time
    w_start = MPI_Wtime();
    
    for (int i = 0; i < _n; i++) {
      for(int j = 0; j < _n; j++ ) {
        if (mat[i][j] != mat[j][i]){ 
          is_symmetric = 0;
          } 
      }
    }
    
    //End time 
    w_end = MPI_Wtime();
    
    time_sequential_cs = w_end - w_start;
    total_time_sequential_cs += time_sequential_cs;
  }
  
    avg_time_sequential_cs = total_time_sequential_cs / TOT_RUN;
    //printf("Execution times of the checksym routine: %12.4g seconds\n", avg_time_sequential_cs); 
  
  
  
  return is_symmetric; 		
}

//calculate the traspose of the matrix 
void matTranspose(float **mat, float **transpose, int _n) {
   	
    for (int i = 0; i < _n; i++) {
    	for(int j = 0; j < _n; j++ ) { 
        transpose[i][j] = mat[j][i]; 
      }
    }
  
}

int checkSymMPI(float **mat, int _n, int rank, int size) {

  int is_symmetric = 1;
  
  int rows_per_process = _n / size; 
  
  int start_row = rank * rows_per_process;
  int end_row = 0; 
  
  if (rank == (size-1)) {
      end_row = _n;
  } else {
      end_row = start_row + rows_per_process;
  }
  
  float *recv = (float *) malloc(_n * rows_per_process * sizeof(float));
  float *send = (float *) malloc(_n * rows_per_process * sizeof(float));
  
  MPI_Alltoall(send, rows_per_process * rows_per_process, MPI_FLOAT, recv, rows_per_process * rows_per_process, MPI_FLOAT, MPI_COMM_WORLD);
    
  for (int k = 0; k < size; ++k) {
    for (int i = 0; i < rows_per_process; ++i) {
      for (int j = 0; j < rows_per_process; ++j) {
        int global_row = start_row + j;
        int global_col = k * rows_per_process + i;
          if (global_row < _n && global_col < _n) { 
            if (mat[global_row][global_col] != recv[k * rows_per_process * rows_per_process + i * rows_per_process + j]) {
              is_symmetric = 0;
            }
          }
        }
      }
    }
    
  free(recv);
  free(send);
  MPI_Barrier(MPI_COMM_WORLD);

  return is_symmetric;
}

//calculate the traspose of the matrix using MPI
void matTransposeMPI(float **mat, float **transpose, int _n, int rank, int size) {
  
  int rows_per_process = _n / size;
  int start_row = rank * rows_per_process;
  int end_row = 0; 
  
  if (rank == (size-1)) {
      end_row = _n;
  } else {
      end_row = start_row + rows_per_process;
  }

  float *send = (float *) malloc(_n * rows_per_process * sizeof(float));
  float *recv = (float *) malloc(_n * rows_per_process * sizeof(float));

  for (int i = 0; i < rows_per_process; ++i) {
      for (int j = 0; j < _n; ++j) {
          send[i * _n + j] = mat[start_row + i][j];
      }
  }

     
  MPI_Alltoall(send, rows_per_process * rows_per_process, MPI_FLOAT, recv, rows_per_process * rows_per_process, MPI_FLOAT, MPI_COMM_WORLD);

  for (int k = 0; k < size; ++k) {
    for (int i = 0; i < rows_per_process; ++i) {
        for (int j = 0; j < rows_per_process; ++j) {
            transpose[start_row + j][k * rows_per_process + i] = recv[k * rows_per_process * rows_per_process + i * rows_per_process + j];
        }
    }
  }
  //free the allocated memory 
  free(send);
  free(recv); 
  MPI_Barrier(MPI_COMM_WORLD);
  
}
