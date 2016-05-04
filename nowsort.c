
/*
 * CSC 258: NowSort Feroz RAuf
 *	Using MPI, openMP, mergeSort,hdfs to read files, and Simd vector operations
 *
 */



#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/time.h>
#include <string.h>
#include <limits.h>


#define BUCKET_SIZE 32000


int proc_id;
int num_procs;
struct timeval start, end;
long total;
long gig_value = 250000000;
long num_elements;
int num_buckets;
int **storage;



void parseArgs(int argc, char ** argv)
{
	char c;
	int num_gigs;
    while ((c = getopt (argc, argv, "n:")) != -1)
        switch (c)
        {
            case 'n':
                num_gigs = atoi(optarg);
                break;
            default:
                printf("Error: malformed arguments!\n");
                exit(-1);
				break;
        }
    num_elements = gig_value * num_gigs;
}
void allocate(int proc_id)
{
	int arr[4] = {3,3,4,10};
	int i,j;
	num_buckets = num_elements/BUCKET_SIZE + 1;
	//load balancing number of buckets across the number of cores 24 cycle 1, 24 cycle 2, 32 cycle 3, and 56 node2x14a
	//keeping number of buckets dynamic, but starting at half as many buckets for entire project
	if(proc_id!=0)
	{
		num_buckets/=2;
	}
	//alocate buckets for each processor
	storage = (int*)malloc(sizeof(int)*(num_buckets));
	//allocate memory for each bucket
	for(i = 0;i < num_buckets;i++)
		storage[i] = (int*)malloc(sizeof(int)*BUCKET_SIZE);
}
int getProc(int num)
{
	int i;
	int arr[4] = {3,6,10,17};
	for(i = 0; i < 4; i++)
	{
		if(num <= INT_MAX* (arr[i]/17))
			return i+1;
	}
	return 4;
}
void sendNumbers()
{
	int i,j,k;
	int load_balance[4] = {3,6,10,17};
	srand(time(NULL));
	int bucket_threshold = INT_MAX/num_buckets;
	int ** arr;
	int *loc;
	loc = (int*)malloc(sizeof(int)*4);
	arr = (int*)malloc(sizeof(int)*4);
	int *num_sent = (int*)malloc(sizeof(int)*4);
	for(i = 0;i < 4; i++)
	{
		loc[i] = 0;
		num_sent[i] = 0;
		arr[i] = (int*)malloc(sizeof(int)*BUCKET_SIZE);
	}
	for(i = 0;i< num_elements;i++)
	{
		int j = rand();
		int proc = getProc(j);
		arr[proc][loc[proc]++] = j;
		for(j = 0;j < 4; j++)
		{
			if(loc[j]==BUCKET_SIZE)
			{
				int orig = j;
				while(num_sent[j]!=((num_buckets/17)*load_balance[j])) j++;
				MPI_Send(&arr[orig],BUCKET_SIZE,MPI_INT,j+1,0,MPI_COMM_WORLD);
				arr[orig] = (int*)malloc(sizeof(int)*BUCKET_SIZE);
				loc[orig] = 0;
				num_sent[j]++;
				break;
			}
		}
	}
}

void receiveNumbers(int proc_id)
{
	int i,j;
	MPI_Status status;
	int load_balance[4] = {3,6,10,17};
	int loc = 0;
	/*for(i = 0;i < num_buckets; i++)
	{
		MPI_Recv(&(storage[i]),BUCKET_SIZE,MPI_INT,0,0,MPI_COMM_WORLD,&status);
	}*/
}
int partition( int* arr, int l, int r) {
   int pivot, i, j, t;
   pivot = arr[l];
   i = l; j = r+1;
		
   while( 1)
   {
   	do ++i; while( arr[i] <= pivot && i <= r );
   	do --j; while( arr[j] > pivot );
   	if( i >= j ) break;
   	t = arr[i]; arr[i] = arr[j]; arr[j] = t;
   }
   t = arr[l]; arr[l] = arr[j]; arr[j] = t;
   return j;
}
void quickSort( int* arr, int l, int r)
{
   int j;

   if( l < r ) 
   {
       j = partition( arr, l, r);
       quickSort( arr, l, j-1);
       quickSort( arr, j+1, r);
   }
	
}
void heapsort(int* arr, unsigned int N) 
{
    if(N==0) 
      return;
    int t; 
    unsigned int n = N, parent = N/2, index, child; 
    while (1) { 
        if (parent > 0) { 
            t = arr[--parent]; 
        } else {
            n--;              
            if (n == 0) {
                return; /
            }
            t = arr[n];        
            arr[n] = arr[0];  
        }
        index = parent; 
        child = index * 2 + 1; 
        while (child < n) {
            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {
                child++;
            }
            if (arr[child] > t) {
                arr[index] = arr[child];
                index = child; 
                child = index * 2 + 1;
            } else {
                break; 
            }
        }
        arr[index] = t; 
    }
}

void sort(int * arr)
{
	quickSort(arr,0,BUCKET_SIZE);
	heapsort(arr,BUCKET_SIZE);
}
void sortNumbers(int proc_id)
{
	int i,j;
	int load_balance[4] = {3,6,10,17};
	int num_buckets = (num_elements/BUCKET_SIZE) * (load_balance[proc_id-1]/17);
	//omp_set_num_threads();
	#pragma omp parallel for
	for (i = 0; i < num_buckets; i++)
	{
		sort(storage[i]);
	}
}

int main(int argc, char **argv)
{  
    long t;
    parseArgs(argc, argv);
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    
    int resultlen;
    char pname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(pname, &resultlen);
    printf("CSUG Computer %s, with proc_id %d\n", pname, proc_id);

    allocate(proc_id);
    //MPI_Barrier(MPI_COMM_WORLD);

    // Wait for all machines to finish allocating memory
    if (proc_id == 0) {
        gettimeofday(&start, NULL);
        sendNumbers();
    } else {
    	receiveNumbers(proc_id);
    	sortNumbers(proc_id);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    if (proc_id == 0) {
    	gettimeofday(&end, NULL);
        t = (end.tv_sec - start.tv_sec) * 1000000;
    	t += end.tv_usec - start.tv_usec;
    	total = ((double) t) / 1000000.0;
    	printf("Completed time: %lf seconds\n",total);
    }
    
    //MPI_Finalize();

    return 0;
}
