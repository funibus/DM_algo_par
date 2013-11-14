#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 100
#define min(a,b) (a<=b?a:b)

int main(int argc, char** argv){
  double d=1/3;
  int i,j,p,l,s;
  int k;
  int nb_step =1;
  MPI_Init(NULL, NULL);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  p=world_size;
  k=n/p;
  int k_col=k, k_row=k;
  if(world_rank/k<n % p){k_row++;};
  if(world_rank % k<n % p){k_col++;};

  int data[n][n];
  srandom(clock());

  for(i=0;i<n;i++){
    for(i=0;i<n;i++){
      if((double)random()/(double)RAND_MAX<d){data[i][j]=1;}
      else{data[i][j]=0;}
    }
  }

  int start_i = world_rank/k * k + min(n% p,world_rank/k); 
  int start_j = world_rank %k * k + min(world_rank%k,n%p);
	printf("\n********\n proc %d :%d,%d\n**********\n", world_rank,start_i, start_j);  
	 int mydata[k_row+2][k_col+2];

    for( i=0; i<k_row; i++){
      for(j=0; j<k_col; j++){

	mydata[i+1][j+1]= data[start_i+i][start_j+j];
      }
    };
  int nord= world_rank - k +p % p;
  int sud = world_rank +k % p;
  int est = world_rank +1 % k;
  int ouest = world_rank -1 +k % k;

  //gérer les send/receive
  MPI_Send(mydata[1], k_col, MPI_INT, nord,  0, MPI_COMM_WORLD);
  MPI_Send(mydata[k_row], k_col, MPI_INT, sud,  0, MPI_COMM_WORLD);
  for(i=1;i<=k_row;i++){
    MPI_Send(&mydata[i][1], 1, MPI_INT, ouest,  i, MPI_COMM_WORLD);
    MPI_Send(&mydata[i][k_col], 1, MPI_INT, est,  i, MPI_COMM_WORLD);
  };
  MPI_Recv(mydata[0], k_col, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
  MPI_Recv(mydata[k_row+1], k_col, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    

  for(i=1;i<=k_row;i++){
    MPI_Recv(&mydata[i][1], 1, MPI_INT, ouest, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[i][k_col], 1, MPI_INT, est, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }


  //  int trnord[k_col], trsud[k_col], trest[k_row], trouest[k_row];//tableaux reçus 
  //int tsnord[k_col], tssud[k_col], tsest[k_row], tsouest[k_row];//tableaux à envoyer

  
  for( s=1;s<=nb_step;s++){
    //étape    
    int voisins[k_row][k_col];
    for( i=1; i<=k_row; i++){
      for(j=1; j<=k_col; j++){
	int compte=0;
	for(l=i-1;l<=i+1;l++){
	  if(mydata[l,j-1]){compte++;};
	  if(mydata[l,j+1]){compte++;}
	  ;}
	if(mydata[i-1+1,j+1]){compte++;};
	if(mydata[i+1,j]){compte++;};
	voisins[i][j]=compte;
	;}
      ;}

    for( i=1; i<=k_row; i++){
      for(j=1; j<=k_col; j++){
	if(voisins[i][j]==3){mydata[i][j]=1;}
	else{if(voisins[i][j]!=2){mydata[i][j]=0;};
	}
      }
    }
  }
  int compte=0;
for( i=1; i<=k_row; i++){
  for(j=1; j<=k_col; j++){
    if(mydata[i][j]){compte++;};
  }
 }
  
  MPI_Finalize();
}
