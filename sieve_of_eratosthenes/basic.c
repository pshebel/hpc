#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) ( BLOCK_LOW((id)+1,p,n)-1 )
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW( (id)+1, p, n) - BLOCK_LOW( (id), p, n  ) )
#define BLOCK_OWNER(index,p,n) ( ( ((p)*(index)+1)-1 ) / (n) )
        
int main (int argc, char *argv[]) {
        int id, p;
        long i, n, size,proc0_size, low_value, high_value, first, index, prime, count, global_count;
        double elapsed_time;
        char *marked;

        MPI_Init (&argc, &argv); 
        MPI_Barrier(MPI_COMM_WORLD); 
        elapsed_time = -MPI_Wtime(); 
        MPI_Comm_rank (MPI_COMM_WORLD, &id); 
        MPI_Comm_size (MPI_COMM_WORLD, &p);
        if (argc != 2) {
                if (!id) printf ("Command line: %s <m>\n", argv[0]); 
                MPI_Finalize(); 
                return 0;
        }
        n = atol(argv[1]);


        // determine which section of the solution space we will work on
        // low = id*n/p
        // high = (id+1)*n/p - 1
        low_value = 2 + BLOCK_LOW(id,p,n-1);
        high_value = 2 + BLOCK_HIGH(id,p,n-1);
        size = BLOCK_SIZE(id,p,n-1);
        proc0_size = (n-1)/p;
        if ((2 + proc0_size) < (int) sqrt(n)) {
                if (!id) printf ("Too many processes\n"); 
                MPI_Finalize();
                return 0;
        }
        
        marked = (char *) malloc (size);
        if (marked == NULL) {
                printf ("Cannot allocate enough memory\n"); 
                MPI_Finalize();
                return 0;
        }

        for (i = 0; i < size; i++) marked[i] = 0; 

        if (!id) index = 0;
        prime = 2;
        do {
                // if (prime * prime > low_value)
                //         first = prime * prime - low_value;
                // else {
                //         if (!(low_value % prime)) first = 0;
                //         else first = prime - (low_value % prime);
                // }
                long start = (low_value / prime) * prime;
                if (start < prime * prime) start = prime * prime;
                if (start < low_value) start += prime;

                for (long j = start; j <= high_value; j += prime) {
                        if (j >= low_value) {
                                marked[j - low_value] = 1;
                        }
                }

                if (!id) {
                        // 0 only
                        // find the smallest unmarked number > prime
                        while (marked[++index]);
                        prime = index + 2;
                }
                MPI_Bcast (&prime, 1, MPI_INT, 0, MPI_COMM_WORLD); 
        } while (prime * prime <= n);

        // count the primes
        count = 0;
        for (i = 0; i < size; i++)
                if (!marked[i]) count++;
        MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        elapsed_time += MPI_Wtime();
        if (!id) {
                printf ("%d primes are less than or equal to %ld\n", global_count, n);
                printf ("Total elapsed time: %10.6f\n", elapsed_time); }
        MPI_Finalize ();
        return 0; 
}
