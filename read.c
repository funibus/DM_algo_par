#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>

#define min(a,b) (a<=b?a:b)
#define max(a,b) (a<=b?b:a)

int sqrti(int x){
  float y =x;
  return floor(sqrt(y));

}

double puissance (int i,int j){ //renvoie i^j
  int k;
  double r=1.0;
  if(j>0){
  for (k=0; k<j;k++){r=i*r;}
  return r;}
  else{  
    for (k=0; k<-j;k++){r=r/i;}
    return r;
  }
}
  
 

int main(int argc, char** argv){
 
  int i,j,p,l,s;
  int k;
  int nb_step =3;
  char* file_name;
  FILE* file;
  if (argc != 3)
    {
      printf("il n'y a pas le bon nombre d'arguments.\n");
      return 1;
    }
  nb_step = atoi(argv[1]);
  file_name=argv[2];
  

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
  
  file= fopen (file_name, "rt");
 int mr, mc;
  fscanf(file,"%d",&mr);
  fscanf(file,"%d",&mc);
  //  printf("%d,%d\n",mr,mc);
  int n= max(mr,mc);  

  k=sqrti(n*n/p); 
  
  int k_col=k, k_row=k;//nb de lignes/colonnes pour le processus
  
  if(row< n % kp){k_row++;};
  if(col< n % kp){k_col++;};

  int start_i = row*k + min(n%kp, row);
  int start_j = col*k + min(n%kp,col);

  int my_rows[k_row+1];
  int tmp;
  double diviseur= puissance(10,mc-(k_col+start_j)); 
  //  if(world_rank==3) printf("%d,%f\n",k_col,diviseur);
    for(i=0;i<start_i; i++)//on ignore les start_i premières lignes
    {
      fscanf(file,"%d", &tmp);
    }

    for(i=0;i<k_row; i++)
    {//on lit les suivantes
      if(i+start_i<mr){ 
	fscanf(file,"%d", &tmp);
	my_rows[i]=floor(tmp/diviseur);
      }
	else{my_rows[i]=0;};
    }

    
  for(i=0;i<start_i; i++)
    {
      scanf("%d", &tmp);
    }
  
     int mydata[k_row+2][k_col+2];
  for( i=0; i<k_row+2; i++){
    for(j=0; j<k_col+2; j++){
      mydata[i][j]= 0;
    }
  };

  for( i=1; i<=k_row; i++){
    for(j=k_col; j>0; j--){
      mydata[i][j]= my_rows[i-1]%10;
      my_rows[i-1]/=10;
    }
  };

  fclose(file);
  MPI_Finalize();

  return 1;
}
