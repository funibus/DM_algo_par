#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>

#define min(a,b) (a<=b?a:b)

int sqrti(int x){
  float y =x;
  return floor(sqrt(y));

}

int main(int argc, char** argv){
  int n=100; 
  double d=0.3;
  int i,j,p,l,s;
  int k;
  int nb_step =0;


if (argc != 4)
    {
        printf("il n'y a pas le bon nombre d'arguments.\n");
        return 1;
    }
 n = atoi(argv[1]);
 nb_step=atoi(argv[2]);
 d=atof(argv[3]);

  MPI_Init(NULL, NULL);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  p=world_size;
  int kp=sqrti(p);//nb de processus par ligne/colonne

  int row= world_rank/kp;
  int col= world_rank%kp;
  
  k=sqrti(n*n/p); 
  int k_col=k, k_row=k;//nb de lignes/colonnes pour le processus
  
  if(row< n % kp){k_row++;};
  if(col< n % kp){k_col++;};

  int data[n][n];
  srand(time(NULL));
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if((double)rand()/(double)RAND_MAX<d){data[i][j]=1;}
      else{data[i][j]=0;}
    }
  }
  
  int start_i = row*k + min(n%kp, row);
  int start_j = col*k + min(n%kp,col);
  printf("\n********\n proc %d :%d,%d\n**********\n", world_rank,start_i, start_j);  
  int mydata[k_row+2][k_col+2];
  for( i=0; i<k_row+2; i++){
    for(j=0; j<k_col+2; j++){
      mydata[i][j]= 0;
    }
  };
  
  for( i=0; i<k_row; i++){
    for(j=0; j<k_col; j++){
      mydata[i+1][j+1]= data[start_i+i][start_j+j];
    }
  };

  int fnord(i){
    return (i - kp +p) % p;
  };
  int fsud(i){
    return (i +kp) % p;
  };
  int fest(i){
    return (i +1) % kp + kp*row;
  };
  int fouest(i){
    return (world_rank -1 +kp) % kp +kp*row;
  }

  int nord= fnord(world_rank);
  int sud= fsud(world_rank);
  int est= fest(world_rank);
  int ouest = fouest(world_rank);


  //pour communiquer directement avec les proc diagonaux (interdit à priori)
  /*  int ne= fnord(fest(world_rank));
  int no=fnord(fouest(world_rank));
  int se=fsud(fest(world_rank));
  int so= fsud(fouest(world_rank));*/


  //gérer les send/receive

  int sendw[k_row+1], sende[k_row+1], recw[k_row+1], rece[k_row+1];
  for(i=1;i<=k_row;i++){
    sendw[i]=mydata[i][1];
    sende[i]=mydata[i][k_col];
  };

  for( s=1;s<=nb_step;s++){

    MPI_Request request;
    MPI_Isend(mydata[1], k_col, MPI_INT, nord,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);
    MPI_Isend(mydata[k_row], k_col, MPI_INT, sud,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);
  
    //  printf("proc %d : %d,moi,%d\n",world_rank,ouest,est);

    MPI_Isend(sendw, k_row, MPI_INT, ouest,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);
    MPI_Isend(sende, k_row, MPI_INT, est,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);
  

    MPI_Recv(mydata[k_row+1], k_col, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
    MPI_Recv(mydata[0], k_col, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
  
    MPI_Recv(rece, k_row, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//&request);
    MPI_Recv(recw, k_row, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //printf("proc %d : message reçu\n",world_rank);

    /*    //communication avec les processus en diagonale
    MPI_Isend(&mydata[0][0],1,MPI_INT, no, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[0][k_col+1],1,MPI_INT, ne, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[k_row+1][0],1,MPI_INT, so, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[k_row+1][k_col+1],1,MPI_INT, se, 0, MPI_COMM_WORLD, &request);

    MPI_Recv(&mydata[k_row+1][k_col+1], 1, MPI_INT, se, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[k_row+1][0], 1, MPI_INT, so, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[0][k_col+1], 1, MPI_INT, ne, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[0][0], 1, MPI_INT, no, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    */


    MPI_Isend(&mydata[0][0],1,MPI_INT, nord, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[k_row+1][0],1,MPI_INT, ouest, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[k_row+1][k_col+1],1,MPI_INT, sud, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[0][k_col+1],1,MPI_INT, est, 0, MPI_COMM_WORLD, &request);

    MPI_Recv(&mydata[k_row+1][0], 1, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[k_row+1][k_col+1], 1, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[k_col+1][0], 1, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[0][0], 1, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


  
    //étape    
    int voisins[k_row][k_col];
    for( i=1; i<=k_row; i++){
      for(j=1; j<=k_col; j++){
	int compte=0;
	for(l=i-1;l<=i+1;l++){
	  if(mydata[l][j-1]){compte++;};
	  if(mydata[l][j+1]){compte++;}
	  ;}
	if(mydata[i-1+1][j+1]){compte++;};
	if(mydata[i+1][j]){compte++;};
	voisins[i][j]=compte;
	;}
      ;}
    
    for( i=1; i<=k_row; i++){
      for(j=1; j<=k_col; j++){

	if(voisins[i][j]){printf("proc %d : %d,%d : %d voisins\n", world_rank,i,j,voisins[i][j]);;}
	if(voisins[i][j]==3){mydata[i][j]=1;}
	else{if(voisins[i][j]!=2){mydata[i][j]=0;};
	}
      }
    }
  }
  int alive=0;
  for( i=1; i<=k_row; i++){
    for(j=1; j<=k_col; j++){
      if(mydata[i][j]){alive++;};
    }
  }
  printf("alive : %d\n", alive);

  int global_alive;
  MPI_Reduce(&alive,&global_alive,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if(world_rank==0){
    printf("global_alive : %d\n", global_alive);
    printf("densité : %f\n", (double)global_alive/(double)(n*n));

  }
  
  MPI_Finalize();
}
