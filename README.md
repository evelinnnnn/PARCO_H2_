# Homework 2: Parallelizing matrix operations using MPI
### Evelin Begher - 235188 - evelin.begher@studenti.unitn.it <br>

I wrote the code in C, and to simplify the execution process, I created a [PBS script](matrix_operation.pbs) that allows me to run the instructions in a single step.<br>
<br>
I can compile my code in another way by directly using the processor and executing the following command: <br>
$ module load gcc91<br>
$ module load mpich-3.2.1--gcc-9.1.0<br>
$ mpicc -o test matrix_transpose.c<br>
$ mpirun -np ***k*** ./test ***n***<br>  
where ***n*** represents the matrix dimension (from $$2^4$$ to $$2^{12}$$) and ***k*** the processes using for MPI(2, 4, 16, 32, 64)<br>
This procedure takes more time, especially when changing the matrix dimensions for various configurations. <br>
<br>

In this [code](matrix_operation.c) we Initialize a random matrix *M* *n×n*. Then, we implemented a function to check if the matrix is symmetric (**checkSym**) . Additionally, we created a function to compute the transpose of *M* and store the result in a new matrix *T* (**matTranspose**). <br>
The assignment requires us to implement the function **checkSym** and the function **matTranspose** using two different methods: sequential implementation (already did it in the previous assignment) and parallelization using MPI . The goal of this assignment is to measure the time taken to perform both the symmetry check and the transpose operation, and to explore how we can optimize performance, especially when working with large datasets or complex calculations.

### Here there are the steps to follow to run my code: 

* Download the [MATRIX_TRANSPOSE_ H1_PBS](matrix_operation.pbs) and [MATRIX_TRANSPOSE_ H1.C](matrix_operation.c) files. These will be used for the job submission and running the program.

* Connect to the HPC Cluster:

    * Set up a secure VPN connection to the University of Trento network.
    * Open your SSH client and connect to the cluster.

*  Enter your university login credentials 
      * Upload the Files: Once connected, navigate to your desired directory on the cluster or create a new one. Use the file transfer feature in your SSH client to upload the [MATRIX_TRANSPOSE_ H1_PBS](matrix_operation.pbs) and [MATRIX_TRANSPOSE_ H1.C](matrix_operation.c) to that directory.

* Request an Interactive Session and Reserve a Node:
     * Navigate to the directory where the files are located: cd ./directory
     * Request an interactive session on a compute node with the desired specifications (ncpus=64:mpiprocs=64:mem=1mb) using the following command: qsub matrix_operation.pbs
     * Once the job completes, you can find the results in a file named job.o in the same directory where the job was submitted.
