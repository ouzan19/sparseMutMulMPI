#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mmio.h"
#define SIZE 512

int findIndex(int* vectorIndex, int col){
	
	int i=0;
	for(;;i++)
		if(vectorIndex[i] == col)
			return i;
	
}

void sort(int* I,int* J,int* V,int n){

 int c,d,swap,position;
 for ( c = 0 ; c < ( n - 1 ) ; c++ )
   {
      position = c;
 
      for ( d = c + 1 ; d < n ; d++ )
      {
         if ( I[position] > I[d] )
            position = d;
      }
      if ( position != c )
      {
         swap = I[c];
         I[c] = I[position];
         I[position] = swap;

	
         swap = J[c];
         J[c] = J[position];
	 J[position] = swap;

	
         swap = V[c];
         V[c] = V[position];
	 V[position] = swap;
      }
   }

}



int main(int argc,char**argv){

	
    int *I,*J;
	int *matrix,*vector,*outvector,*out;
	int myrank,P;
	double tstart,tend;
	MPI_Status status;

	MPI_Init (&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&P);


	int perProcessor = SIZE/P;
	
	if(myrank == 0){
	
		//printf("Hello world\n");

		srand(time(0));
		int ret_code;
		MM_typecode matcode;
		FILE *f;
		int M, N, nz;   
		int i;
		

		if (argc < 2)
		{
		//	printf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
			return 1;
		}
		else    
		{ 
			if ((f = fopen(argv[1], "r")) == NULL) 
			{
		//		printf(stderr, "No file named %s \n", argv[1]);
				return 1;
			
			}
		}

		if (mm_read_banner(f, &matcode) != 0)
		{
			printf("Could not process Matrix Market banner.\n");
			return 1;
		}


	
	

		if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
				mm_is_sparse(matcode) )
		{
			printf("Sorry, this application does not support ");
			//printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
			return 1;
		}

	

		if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
			return 1;


		//printf("Vector:\n");
		vector = (int*) malloc(sizeof(int)*N);
		outvector = (int*) malloc(sizeof(int)*M);
		
		for(i=0;i<M;i++)
			outvector[i] = 0;
		for(i=0;i<N;i++){
			vector[i] =i+1;// rand() % 100;
			//printf("%d %d\n",i,vector[i]);
		}
		out    = (int*) malloc(sizeof(int)*M);
		I = (int *) malloc(nz * sizeof(int));
		J = (int *) malloc(nz * sizeof(int));
		matrix = (int *) malloc(nz * sizeof(int));


	

		for (i=0; i<nz; i++)
		{
			fscanf(f, "%d %d %d\n", &I[i], &J[i], &matrix[i]);
			I[i]--;  
			J[i]--;
		}

		if (f !=stdin) fclose(f);

		//sort(I,J,matrix,nz);

		mm_write_banner(stdout, matcode);
		mm_write_mtx_crd_size(stdout, M, N, nz);
		//for (i=0; i<nz; i++)
			//fprintf(stdout, "%d %d %d\n", I[i], J[i], matrix[i]);
		
		/*
		if(P==1){
			
			tstart = MPI_Wtime();
			for(i=0;i<nz;i++)
				outvector[ I[i] ] += matrix[i]*vector[J[i]];
			
		
		
		
			tend = MPI_Wtime();
			printf("\nFinal Result\n");
			for(i=0;i<M;i++)
				printf("%d ",outvector[i]);
			printf("\n");
		
			printf( "Elapsed time is %f\n", tend - tstart ); 
		
		
		
			MPI_Finalize();
			return 0;
			
			
			
		}*/
	
		int rowPerProcessor = M/P;
		//int **toProcessors = (int**) malloc(sizeof(int*)*P);
		//int **toProcessorsIndex = (int**) malloc(sizeof(int*)*P);
		//int **toProcessorsI = (int**) malloc(sizeof(int*)*P);
		//int* counts = (int*) malloc(sizeof(int)*P);
		//int* countsM = (int*)malloc(sizeof(int)*P);
	    int* selected = (int*) malloc(sizeof(int)*N);	
		//int **toProcessorsJ = (int**) malloc(sizeof(int*)*P);
		//int **toProcessorsM = (int**) malloc(sizeof(int*)*P);
		
		int* countsToSend = (int*) malloc(sizeof(int)*3*P);
		int** dataToSend = (int**) malloc(sizeof(int*)*P);
			
		for(i=0;i<P;i++){

			countsToSend[i*3] = rowPerProcessor;
		
			int start = i*rowPerProcessor;
			int end = (i+1)* rowPerProcessor;
			int count = 0;
			int j;
			

			
			for(j=0;j<nz;j++){	
				
				if(I[j] >= start && I[j]<end  ){

					count++;	

				}
			}
			
			//toProcessorsI[i] = (int*) malloc(sizeof(int)*count); 
			//toProcessorsJ[i] = (int*) malloc(sizeof(int)*count); 
			//toProcessorsM[i] = (int*) malloc(sizeof(int)*count); 
			//countsM[i] = count;
	
			countsToSend[3*i + 1] = count;
			
			count = 0;
			for(j=0;j<nz;j++)
				selected[J[j]] = 0;
			for(j=0;j<nz;j++){
				
				if(!selected[J[j]] && I[j] >= start && I[j]<end  )
				{
					count++;
					selected[J[j]] = 1;
				} 	 

			}
			//toProcessors[i] = (int*) malloc(sizeof(int)*count);
			//counts[i] = count;
			countsToSend[3*i + 2] = count;
			//toProcessorsIndex[i] = (int*) malloc(sizeof(int)*count);
			
			dataToSend[i] = (int*) malloc(sizeof(int)* (3*countsToSend[3*i + 1] + 2*countsToSend[3*i + 2]));
			
		}

		
		for(i=0;i<P;i++){
			int dataToSendIndex = 0;
			int start = i*rowPerProcessor;
			int end = (i+1)* rowPerProcessor;
			int count = 0;
			int j;
			
			/*
			for(j=0;j<nz;j++){	
				
				if(I[j] >= start && I[j]<end  ){

					toProcessorsI[i][count] = I[j];
					toProcessorsM[i][count] = matrix[j];
					toProcessorsJ[i][count] = J[j];					
					count++;
					
				}
				
			}*/
			
			for(j=0;j<nz;j++)
				if(I[j] >= start && I[j]<end ) 
					dataToSend[i][dataToSendIndex++] = I[j];
					
			for(j=0;j<nz;j++)		
				if(I[j] >= start && I[j]<end  )
					dataToSend[i][dataToSendIndex++] = J[j];
					
			for(j=0;j<nz;j++)
				if(I[j] >= start && I[j]<end  )
					dataToSend[i][dataToSendIndex++] = matrix[j];					
				
				
			/*	
			count = 0;
			for(j=0;j<nz;j++)
				selected[J[j]] = 0;
			for(j=0;j<nz;j++){
	
				if(!selected[J[j]] && I[j] >= start && I[j]<end  ){
					
					toProcessors[i][count] = vector[J[j]];
					
					toProcessorsIndex[i][count] = J[j];
					selected[J[j]] = 1;					
					count++;	
					

				}
			}
			*/
			
			for(j=0;j<nz;j++)
				selected[J[j]] = 0;
			
			for(j=0;j<nz;j++)
				if(!selected[J[j]] && I[j] >= start && I[j]<end  ){
					
					dataToSend[i][dataToSendIndex++] = J[j];
					selected[J[j]] = 1;	
					
				}
					
					
			for(j=0;j<nz;j++)
				selected[J[j]] = 0;
			
			for(j=0;j<nz;j++)
				if(!selected[J[j]] && I[j] >= start && I[j]<end  ){
					
					dataToSend[i][dataToSendIndex++] = vector[J[j]];	
					selected[J[j]] = 1;	

				}					
			
				
		
			
			
			
			
			
			//printf("--------------------------------------\n");
		}

		/*
		for(i=0;i<P;i++){
			
			printf("Processes %d: \n",i);
			
			
			printf("Vector:\n");
			int j;
			for(j=0;j<counts[i];j++){
				
				printf("%d %d\n",toProcessorsIndex[i][j],toProcessors[i][j]);


			}
			
			printf("Matrix:\n");
			for(j=0;j<countsM[i];j++){
				
				printf("%d %d %d\n",toProcessorsI[i][j],toProcessorsJ[i][j],toProcessorsM[i][j]);


			}
			printf("---------------------------------------------------\n");


		}
		*/
		
		MPI_Status *status = (MPI_Status*) malloc(sizeof(MPI_Status)*P);
        MPI_Request *reqs =  (MPI_Request*) malloc(sizeof(MPI_Request)*P);
		 MPI_Request *reqs2 =  (MPI_Request*) malloc(sizeof(MPI_Request)*P);

		
		for(i=1;i<P;i++){
			
			
			 MPI_Isend(&(countsToSend[3*i]),3,MPI_INT,i,0,MPI_COMM_WORLD,&(reqs[i]));
			
			MPI_Isend(dataToSend[i],3*countsToSend[3*i + 1] + 2*countsToSend[3*i + 2],MPI_INT,i,1,MPI_COMM_WORLD,&(reqs2[i]));

			
			//MPI_Send(&(countsToSend[3*i]),3,MPI_INT,i,0,MPI_COMM_WORLD);
			//MPI_Send(dataToSend[i],3*countsToSend[3*i + 1] + 2*countsToSend[3*i + 2],MPI_INT,i,1,MPI_COMM_WORLD);
			
			//MPI_Send(&rowPerProcessor,1,MPI_INT,i,0,MPI_COMM_WORLD);
			//MPI_Send(&(counts[i]),1,MPI_INT,i,1,MPI_COMM_WORLD);
			//MPI_Send(toProcessors[i],counts[i],MPI_INT,i,2,MPI_COMM_WORLD);	
			//MPI_Send(toProcessorsIndex[i],counts[i],MPI_INT,i,3,MPI_COMM_WORLD);
	
			//MPI_Send(&(countsM[i]),1,MPI_INT,i,4,MPI_COMM_WORLD);
			//MPI_Send(toProcessorsI[i],countsM[i],MPI_INT,i,5,MPI_COMM_WORLD);	
			//MPI_Send(toProcessorsJ[i],countsM[i],MPI_INT,i,6,MPI_COMM_WORLD);	
			//MPI_Send(toProcessorsM[i],countsM[i],MPI_INT,i,7,MPI_COMM_WORLD);	
	
		}
		tstart = MPI_Wtime();	
		int count,count2;
		rowPerProcessor = countsToSend[0];
		count2 = countsToSend[1];
		count = countsToSend[2];
		
		int *out = (int*) malloc(sizeof(int)*rowPerProcessor);
		int *outReceived = (int*) malloc(sizeof(int)*rowPerProcessor);
		for(i=0;i<rowPerProcessor;i++)
			out[i] = 0;
		
		
		for(i=0;i<count2;i++){
			
			int index = findIndex(dataToSend[0]+3*count2,dataToSend[0][count2+i]);
			
			
			outvector[ dataToSend[0][i] - myrank*rowPerProcessor ] += dataToSend[0][3*count2+count+index]*dataToSend[0][2*count2+i];
			
		}
		
	
			
		
		for(i=1;i<P;i++){
			
			MPI_Wait(&(reqs2[i]), &(status[i])); 
			MPI_Status stat;
			MPI_Request req;
			MPI_Irecv(out,rowPerProcessor,MPI_INT,i,8,MPI_COMM_WORLD,&req);
			
			//MPI_Recv(out,rowPerProcessor,MPI_INT,i,8,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Wait(&req,&stat); 
			
			int j;
			for(j=0;j<rowPerProcessor;j++)
				outvector[i*rowPerProcessor + j] = out[j];
				
			
		}
		
		
		tend = MPI_Wtime();
	/*	
		printf("\nFinal Result\n");
		for(i=0;i<M;i++)
			printf("%d ",outvector[i]);
		printf("\n");
	*/	
		printf( "Elapsed time is %f\n", tend - tstart ); 
		
	} 

	else {


		int count,count2,i,rowPerProcessor;
		//int *receivedVector, *receivedVectorIndex;
		//int *I, *J, *M;
		int* countsReceived = (int*) malloc(sizeof(int)*3);
		int* receivedData;
		
		MPI_Status stat;
		MPI_Request req;

		
		MPI_Irecv(countsReceived,3,MPI_INT,0,0,MPI_COMM_WORLD,&req);
		
		//MPI_Recv(countsReceived,3,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		MPI_Wait(&req,&stat); 
		
		rowPerProcessor = countsReceived[0];
		count2 = countsReceived[1];
		count = countsReceived[2];
		
		int *out = (int*) malloc(sizeof(int)*rowPerProcessor);
		for(i=0;i<rowPerProcessor;i++)
			out[i] = 0;
		
		//MPI_Recv(&count,1,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		//receivedVector = (int*) malloc(sizeof(int)*count);
		//receivedVectorIndex= (int*) malloc(sizeof(int)*count);
		
		
		//I = (int*) malloc(sizeof(int)*count2);
		//J = (int*) malloc(sizeof(int)*count2);
		//M = (int*) malloc(sizeof(int)*count2);
		
		receivedData = (int*) malloc(sizeof(int)* (2*count + 3*count2));
		
		//MPI_Recv(receivedData,2*count + 3*count2,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Irecv(receivedData,2*count + 3*count2,MPI_INT,0,1,MPI_COMM_WORLD,&req);
		MPI_Wait(&req,&stat); 
		/*
		printf("\n//////////////////////////////////////////////////////////////////////\n");
		for(i=0;i<2*count + 3*count2;i++)
			printf("%d ",receivedData[i]);
		
		printf("\n//////////////////////////////////////////////////////////////////////\n");
		*/
		//memcpy((void*)I,(void*)receivedData,count2);
		//memcpy((void*)J,(void*)(receivedData+count2),count2);
		//memcpy((void*)M,(void*)(receivedData+2*count2),count2);
		//memcpy((void*)receivedVectorIndex,(void*)(receivedData+3*count2),count);
		//memcpy((void*)receivedVector,(void*)(receivedData+3*count2+count),count);
		
		

		//MPI_Recv(receivedVector,count,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//MPI_Recv(receivedVectorIndex,count,MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		//printf("********************************************\nVector\n");
		//for(i=0;i<count;i++)	
			//printf("%d,%d - ",receivedVectorIndex[i],receivedVector[i]);
		

		//MPI_Recv(&count2,1,MPI_INT,0,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
				

		//MPI_Recv(I,count2,MPI_INT,0,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
		//MPI_Recv(J,count2,MPI_INT,0,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//MPI_Recv(M,count2,MPI_INT,0,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		//printf("\nMatrix\n");

		//for(i=0;i<count2;i++)
			//printf("%d, %d, %d - ",I[i],J[i],M[i]);
	



		for(i=0;i<count2;i++){
			
			int index = findIndex(receivedData+3*count2,receivedData[count2+i]);
			
			
			out[ receivedData[i] - myrank*rowPerProcessor ] += receivedData[3*count2+count+index]*receivedData[2*count2+i];
			
		}
		/*
		printf("\nResult\n");
		for(i=0;i<rowPerProcessor;i++)
			printf("%d ",out[i]);
		*/
		//printf("\n****************************************************\n");
		
		
		MPI_Isend(out,rowPerProcessor,MPI_INT,0,8,MPI_COMM_WORLD,&req);	

		//MPI_Send(out,rowPerProcessor,MPI_INT,0,8,MPI_COMM_WORLD);	






	}

	MPI_Finalize();
	return 0;

}
