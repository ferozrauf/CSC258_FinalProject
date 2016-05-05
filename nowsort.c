
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


#define BUCKET_SIZE 64000


int proc_id;
int num_procs;
clock_t begin, end;
long total;
long gig_value = 250000000;
long num_elements;
int num_buckets;
int **storage;
int **buckets;
int last_bucket_length;
int *bucket_length;
int *load_balance;
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
	int i,j,k;
	num_buckets = num_elements/BUCKET_SIZE + 1;
	printf("Number of elements %d in bucket size %d, so num_buckets is %d \n",num_elements,BUCKET_SIZE,num_buckets);
	load_balance = malloc(sizeof(int)*4);
	for(i = 0; i< 4 ;i++)
		load_balance[i] = i==0 ? arr[0] : load_balance[i-1] + arr[i];
	//load balancing number of buckets across the number of cores 24 cycle 1, 24 cycle 2, 32 cycle 3, and 56 node2x14a
	//keeping number of buckets dynamic, but starting at half as many buckets for entire project
	storage = (int*)malloc(sizeof(int*)*num_buckets);
	for(i = 0; i< num_buckets; i++)
		storage[i] = (int*)malloc(sizeof(int)*BUCKET_SIZE);
	/*int alloc = num_buckets * arr[proc_id];
	alloc /= 17;
	alloc += 1;*/
	bucket_length = (int*)calloc(num_buckets,sizeof(int));
	for(i = 0;i < 4;i++)
	{	
		int temp = num_buckets/17;
		temp += 1;
		load_balance[i] *= temp;
	}
	load_balance[3] = num_buckets;
}

int getProc(int bucket_num)
{
	//int load_balance[4] = {3,6,10,17};
	int j;
	for(j = 0; j < 4; j++)
	{
		if(bucket_num<load_balance[j])
			return j;
	}
	return 4;
}
int getBucket(int num,int num_buckets,int *bucket_size)
{
	int i = INT_MAX/num_buckets;
	i = num/i;
	i = i % num_buckets;
	//printf("Initial bucket placement %d \n",i);
	while(bucket_size[i]==BUCKET_SIZE && i<num_buckets) 
	{
		i++;
		i = i % num_buckets;
	}
	//printf("Final bucket placement %d \n",i);

	return i;
}
void sendNumbers()
{
	int i,j,k,l;
	//int load_balance[4] = {3,6,10,17};
	//int bucket_threshold = INT_MAX/num_buckets;
	srand(time(NULL));
	printf("Doing bucket sort on random numbers, num_buckets %d \n",num_buckets);

	for(i = 0;i< num_elements;i++)
	{
		int temp = rand();
		j = getBucket(temp,num_buckets,bucket_length);
		//printf("Adding the %dith number: %d ,into bucket %d of length %d\n",i,temp,j,bucket_length[j]);
		storage[j][bucket_length[j]] = temp;
		bucket_length[j]++;
	}
	printf("Buckets allocated, num_buckets %d \n",num_buckets);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(bucket_length,num_buckets,MPI_INT,0,MPI_COMM_WORLD);
	printf("PROC_ID %d has it's buckets\n",0);
	//MPI_Bcast(load_balance,4,MPI_INT,0,MPI_COMM_WORLD);
	for(i = load_balance[0]; i< num_buckets; i++)
	{
		MPI_Send(storage[i],bucket_length[i],MPI_INT,getProc(i),i,MPI_COMM_WORLD);
	}
	printf("Buckets sent, sent %d of them \n",num_buckets);
}

void receiveNumbers(int proc_id)
{
	int i,j;
	MPI_Status status;
	printf("PROC_ID %d is waiting on process 0 \n",proc_id);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(bucket_length,num_buckets,MPI_INT,0,MPI_COMM_WORLD);
	printf("ID: %d received bucket lengths \n",proc_id);
	//buckets = malloc(sizeof(int*)*load_balance[proc_id]);
	printf("ID: %d getting buckets from %d to %d \n",proc_id,load_balance[proc_id-1],load_balance[proc_id]);

	for(i = load_balance[proc_id-1]; i < num_buckets && i < load_balance[proc_id] ; i++)
	{
		MPI_Recv(storage[i],bucket_length[i],MPI_INT,0,i,MPI_COMM_WORLD,&status);

	}

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
                return; 
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

void sort(int * arr,int length)
{
	quickSort(arr,0,length);
}
void sortNumbers(int proc_id)
{
	int i,j,k,start,end;
	start = proc_id==0 ? 0 : load_balance[proc_id-1];
	end = load_balance[proc_id];
	int threads[4] = {24,24,32,56};
	omp_set_num_threads(threads[proc_id]);
	#pragma omp parallel for
	for(i = start;i < end; i++)
		sort(storage[i],bucket_length[i]);
	printf("Process %d has finished sorting buckets\n",proc_id);
	/*int length = 0;
	for(i = start; i< end;i++)
		length += bucket_length[i];
	int* arr = malloc(sizeof(int)*length);
	
	k = 0;
	for(i = 0;i < end ; i++)
	{
		for(j = 0; j < bucket_length[i];j++)
			arr[k + j] = buckets[i][j];
		k+=bucket_length[i];
	}
	sort(arr,length);*/
}

int main(int argc, char **argv)
{  
	struct timeval  tv1, tv2;

    parseArgs(argc, argv);
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    
    int resultlen;
    char pname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(pname, &resultlen);
    printf("CSUG Computer %s, with proc_id %d\n", pname, proc_id);

    allocate(proc_id);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Memory allocation for %s, process_id %d is completed\n",pname,proc_id);
    // Wait for all machines to finish allocating memory
    if (proc_id == 0) {
		gettimeofday(&tv1, NULL);
        sendNumbers();
    } else {
    	receiveNumbers(proc_id);
    }
    printf("CSUG computer %s proc_id %d is starting sorting \n",pname,proc_id);
    sortNumbers(proc_id);
    MPI_Barrier(MPI_COMM_WORLD);
    if (proc_id == 0) {
    	gettimeofday(&tv2, NULL);

		printf ("Total time = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
    }
    
    MPI_Finalize();

    return 0;
}
