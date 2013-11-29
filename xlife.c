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

int maxi(a,b){
  if (a < b)
    return b;
  return a;
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
  int n; 
  int i,j,p,s;
  int k;
  int nb_step;
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
  n= maxi(mr,mc); 
  k=sqrti(n*n/p); 
  int k_col=k+2, k_row=k+2;//nb de lignes/colonnes pour le processus, on ajoute 2 pour les lignes et colonnes des voisins
  
  if(row< n % kp){k_row++;};
  if(col< n % kp){k_col++;};

  int prochain_expand_col = n%kp; //c'est la plus petite colonne qui a une colonne de moins que sa voisine de gauche
  int prochain_expand_row = n%kp; //idem pour les lignes

  
  int start_i = row*k + min(n%kp, row); //indice de la premiere ligne du "monde" dont je suis chargé
  int start_j = col*k + min(n%kp,col); //indice de la premiere colonne du "monde" dont je suis chargé
	
  int my_rows[k_row];
  int tmp;
  double diviseur= puissance(10,mc-(k_col-2+start_j)); 
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
  int** mydata =(int **) malloc((k_row)*sizeof(int*));
  for (i = 0; i<k_row; i++)
    {
    mydata[i] = (int*) calloc((k_col),sizeof(int));
    }
    for( i=0; i<k_row; i++){
    for(j=0; j<k_col; j++){
      mydata[i][j]= 0;
    }
  };

  for( i=1; i<=k_row-2; i++){
    for(j=k_col-2; j>0; j--){
      mydata[i][j]= my_rows[i-1]%10;
      my_rows[i-1]/=10;
    }
  };

  /*  printf("proc %d\n",world_rank);
  for( i=1; i<=k_row-2; i++){
    for(j=1; j<=k_col-2; j++){
      printf("%d ", mydata[i][j]);
    }
    printf("\n");
    };*/

  fclose(file);

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

for (i = 0; i<k_row; i++)
	{
		for (j = 0; j<k_col; j++)
			fprintf(stderr, "proc %d, case (%d,%d) il y a %d\n", world_rank, i, j, mydata[i][j]);
}

	MPI_Request request;
	int temp;
	for( s=1;s<=nb_step;s++){

		// voir s'il faudra étendre le tableau
		int expand_n = 0, expand_s = 0, expand_o = 0, expand_e = 0;
		if(row == 0)
		{
		  for(j=0; j<k_col-1;j++)
				{
					if (mydata[1][j] == 1)
						expand_n = 1;
				}
		}

		if(row == kp-1)
		{
		  for(j=0; j<k_col-1;j++)
		    {
					if (mydata[k_row-2][j] == 1)
						expand_s = 1;
				}
		}

		if(col == 0)
		{
		  for(i=0; i<k_row-1;i++)
		    {
					if (mydata[i][1] == 1)
						expand_o = 1;
				}
		}

		if(col == kp-1)
		{
		  for(i=0; i<k_row-1;i++)
		    {
					if (mydata[i][k_col-2] == 1)
						expand_e = 1;
				}
		}

//**********les lignes et colonnes du bord se mettent d'accord pour savoir s'il faut etendre ou pas**********

		if (col == 0) //les processeurs de la col 0 mettent en commun pour savoir s'ils doivent etendre
		{
			if (row != 0) //propagation des 1 vers le sud dans la 1ere colonne
			{
				MPI_Recv(&temp, 1, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  expand_o = maxi(expand_o,temp);
			}

			if (row != kp-1)
			{
				MPI_Isend(&expand_o, 1, MPI_INT, sud,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}

			if (row != kp-1) //on remonte
				MPI_Recv(&temp, 1, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			expand_o = maxi(expand_o,temp);
			if (row != 0)
			{
				MPI_Isend(&expand_o, 1, MPI_INT, nord,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}
		}

		if (col == kp-1) //processeurs de la col kp-1
		{
			if (row != 0) //propagation des 1 vers le sud dans la derniere colonne
				MPI_Recv(&temp, 1, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			expand_e = maxi(expand_e,temp);
			if (row != kp-1)
			{
				MPI_Isend(&expand_e, 1, MPI_INT, sud,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}

			if (row != kp-1) //on remonte
				MPI_Recv(&temp, 1, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			expand_e = maxi(expand_e,temp);
			if (row != 0)
			{
				MPI_Isend(&expand_e, 1, MPI_INT, nord,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}
		}

		if (row == 0) //les processeurs de la ligne 0 mettent en commun pour savoir s'ils doivent etendre
		{
			if (col != 0) //propagation des 1 vers l'est dans la 1ere ligne
				MPI_Recv(&temp, 1, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			expand_n = maxi(expand_n,temp);
			if (col != kp-1)
			{
				MPI_Isend(&expand_n, 1, MPI_INT, est,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}

			if (col != kp-1) //vers l'ouest
				MPI_Recv(&temp, 1, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			expand_n = maxi(expand_n,temp);
			if (col != 0)
			{
				MPI_Isend(&expand_n, 1, MPI_INT, ouest,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}
		}

		if (row == kp-1) //les processeurs de la derniere ligne mettent en commun pour savoir s'ils doivent etendre
		{
			if (col != 0) //propagation des 1 vers l'est dans la derniere ligne
				MPI_Recv(&temp, 1, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			expand_s = maxi(expand_s,temp);
			if (col != kp-1)
			{
				MPI_Isend(&expand_s, 1, MPI_INT, est,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}

			if (col != kp-1) //vers l'ouest
				MPI_Recv(&temp, 1, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			expand_s = maxi(expand_s,temp);
			if (col != 0)
			{
				MPI_Isend(&expand_s, 1, MPI_INT, ouest,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
			}
		}
	//*******les lignes et colonnes du bord envoient aux autres proc s'il faut etendre ou pas, et dans quelle direction******

		if (col != 0) //la colonne de gauche envoie ses infos
			MPI_Recv(&expand_o, 1, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (col != kp-1)
		{
				MPI_Isend(&expand_o, 1, MPI_INT, est,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
		}

		if (col != kp-1) //la colonne de droite envoie ses infos
			MPI_Recv(&expand_e, 1, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (col != 0)
		{
				MPI_Isend(&expand_e, 1, MPI_INT, ouest,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
		}

		if (row != 0) //la ligne du haut envoie ses infos
			MPI_Recv(&expand_n, 1, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (row != kp-1)
		{
				MPI_Isend(&expand_n, 1, MPI_INT, sud,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
		}

		if (row != kp-1) //la ligne du bas envoie ses infos
			MPI_Recv(&expand_s, 1, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (row != 0)
		{
				MPI_Isend(&expand_s, 1, MPI_INT, nord,  0, MPI_COMM_WORLD,&request);
		  	MPI_Request_free(&request);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	//******gerer l'expansion******

		if (expand_n == 1) //s'il faut etendre au nord
		{
			if (row == prochain_expand_row)
			{
				k_row++;
				mydata = (int**) realloc(mydata, k_row*sizeof(int*));
				for (i = k_row-1; i>0; i--)
					mydata[i] = mydata[i-1];
				mydata[0] = (int*) malloc(k_col*sizeof(int));
				if (row != 0)
					MPI_Recv(mydata[0], k_col, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				else
				{
					for (j = 0; j<k_col; j++)
						mydata[0][j] = 0;
				}
			}

			else if (row == 0 && row < prochain_expand_row)
			{
				MPI_Send(mydata[k_row-1], k_col, MPI_INT, sud,  0, MPI_COMM_WORLD);
				free(mydata[k_row-1]);
				for (i = k_row-1; i>0; i--)
					mydata[i] = mydata[i-1];
				mydata[0] = malloc(k_col*sizeof(int));
				for (j = 0; j<k_col; j++)
					mydata[0][j] = 0;
			}
			else if (row < prochain_expand_row)
			{
				int* realloc_nord = malloc(k_col*sizeof(int));
				MPI_Recv(realloc_nord, k_col, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(mydata[k_row-1], k_col, MPI_INT, sud,  0, MPI_COMM_WORLD);
				free(mydata[k_row-1]);
				for (i = k_row-1; i>0; i--)
					mydata[i] = mydata[i-1];
				mydata[0] = malloc(k_col*sizeof(int));
				for (j = 0; j<k_col; j++)
					mydata[0][j] = realloc_nord[j];
				free(realloc_nord);
			}
			prochain_expand_row = (prochain_expand_row +1)%kp;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (expand_s == 1) //s'il faut etendre au sud
		{
			if (row == prochain_expand_row)
			{
				k_row++;
				mydata = (int**) realloc(mydata, k_row*sizeof(int));
				mydata[k_row-1] = malloc(k_col*sizeof(int));
				if (row != kp-1)
					MPI_Recv(mydata[k_row-1], k_col, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				else
				{
					for (j = 0; j<k_col; j++)
						mydata[k_row-1][j] = 0;
				}
			}

			else if (row == kp-1 && row > prochain_expand_row)
			{
				MPI_Send(mydata[0], k_col, MPI_INT, nord,  0, MPI_COMM_WORLD);
				free(mydata[0]);
				for (i = 0; i<k_row-1; i++)
					mydata[i] = mydata[i+1];
				mydata[k_row-1] = malloc(k_col*sizeof(int));
				for (j = 0; j<k_col; j++)
					mydata[k_row-1][j] = 0;
			}
			else if (row > prochain_expand_row)
			{
				int* realloc_sud = malloc(k_col*sizeof(int));
				MPI_Recv(realloc_sud, k_col, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(mydata[0], k_col, MPI_INT, nord,  0, MPI_COMM_WORLD);
				free(mydata[0]);
				for (i =0; i<k_row-1; i++)
					mydata[i] = mydata[i+1];
				mydata[0] = malloc(k_col*sizeof(int));
				for (j = 0; j<k_col; j++)
					mydata[k_row-1][j] = realloc_sud[j];
				free(realloc_sud);
			}
			prochain_expand_row = (prochain_expand_row +1)%kp;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (expand_o == 1) //s'il faut etendre a l'ouest
		{
			if (col == prochain_expand_col)
			{
				k_col++;
				for (i = 0; i<k_row; i++)
					mydata[i] = (int*) realloc(mydata[i], k_col*sizeof(int));
				int* temp_rec = malloc(k_row*sizeof(int));
				for (i = 0; i<k_row; i++)
				{
					for (j = k_col-1; j>0; j--)
						mydata[i][j] = mydata[i][j-1];
				}
				if (col != 0)
				{
					MPI_Recv(temp_rec, k_row, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for (i = 0; i<k_row; i++)
						mydata[i][0] = temp_rec[i];
				}
				else
				{
					for (i = 0; i<k_row; i++)
						mydata[i][0] = 0;
				}
				free(temp_rec);
			}

			else if (col == 0 && col < prochain_expand_col)
			{
				int* temp_envoi = malloc(k_row*sizeof(int));
				for (i = 0; i< k_row; i++)
					temp_envoi[i] = mydata[i][k_col-1];
				MPI_Send(temp_envoi, k_row, MPI_INT, est,  0, MPI_COMM_WORLD);
				for (i = 0; i<k_row; i++)
				{
					for (j = k_col-1; j>0; j--)
						mydata[i][j] = mydata[i][j-1];
				}
				for (i = 0;i<k_row; i++)
					mydata[i][0] = 0;
				free(temp_envoi);
			}
			else if (col < prochain_expand_col)
			{
				int* realloc_ouest = malloc(k_row*sizeof(int));
				int* temp_envoi = malloc(k_row*sizeof(int));
				for (i = 0; i< k_row; i++)
					temp_envoi[i] = mydata[i][k_col-1];
				MPI_Recv(realloc_ouest, k_row, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(temp_envoi, k_row, MPI_INT,est,  0, MPI_COMM_WORLD);
				for (i = 0; i<k_row; i++)
				{
					for (j = k_col-1; j>0; j--)
						mydata[i][j] = mydata[i][j-1];
				}
				for (i = 0;i<k_row; i++)
					mydata[i][0] = realloc_ouest[i];
			
				free(temp_envoi);
				free(realloc_ouest);
			}
			prochain_expand_col = (prochain_expand_col +1)%kp;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (expand_e == 1) //s'il faut etendre a l'est
		{
			if (col == prochain_expand_col)
			{
				k_col++;
				for (i = 0; i<k_row; i++)
					mydata[i] = (int*) realloc(mydata[i], k_col*sizeof(int));
				int* temp_rec = malloc(k_row*sizeof(int));

				if (col != kp-1)
				{
					MPI_Recv(temp_rec, k_row, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for (i = 0; i<k_row; i++)
						mydata[i][k_col-1] = temp_rec[i];
				}
				else
				{
					for (i = 0; i<k_row; i++)
						mydata[i][k_col-1] = 0;
				}
				free(temp_rec);
			}

			else if (col == kp-1 && col > prochain_expand_col)
			{
				int* temp_envoi = malloc(k_row*sizeof(int));
				for (i = 0; i< k_row; i++)
					temp_envoi[i] = mydata[i][0];
				MPI_Send(temp_envoi, k_row, MPI_INT, ouest,  0, MPI_COMM_WORLD);
				for (i = 0; i<k_row; i++)
				{
					for (j = 0; j<k_col-1; j++)
						mydata[i][j] = mydata[i][j+1];
				}
				for (i = 0;i<k_row; i++)
					mydata[i][k_col-1] = 0;
				free(temp_envoi);
			}
			else if (col > prochain_expand_col)
			{
				int* realloc_est = malloc(k_row*sizeof(int));
				int* temp_envoi = malloc(k_row*sizeof(int));
				for (i = 0; i< k_row; i++)
					temp_envoi[i] = mydata[i][0];
				MPI_Recv(realloc_est, k_row, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(temp_envoi, k_row, MPI_INT,ouest,  0, MPI_COMM_WORLD);
				for (i = 0; i<k_row; i++)
				{
					for (j = 0; j<k_col-1; j++)
						mydata[i][j] = mydata[i][j+1];
				}
				for (i = 0;i<k_row; i++)
					mydata[i][k_col-1] = realloc_est[i];
			
				free(temp_envoi);
				free(realloc_est);
			}
			prochain_expand_col = (prochain_expand_col +1)%kp;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		//gérer les send/receive

    int sendw[k_row], sende[k_row], recw[k_row], rece[k_row];
		for(i=1;i<k_row-1;i++){
		  sendw[i]=mydata[i][1];
		  sende[i]=mydata[i][k_col-2];
		};

    MPI_Isend(mydata[1], k_col, MPI_INT, nord,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);
    MPI_Isend(mydata[k_row-2], k_col, MPI_INT, sud,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);

    MPI_Recv(mydata[k_row-1], k_col, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
    MPI_Recv(mydata[0], k_col, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

    MPI_Isend(sendw, k_row, MPI_INT, ouest,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);
    MPI_Isend(sende, k_row, MPI_INT, est,  0, MPI_COMM_WORLD,&request);
    MPI_Request_free(&request);
  
    MPI_Recv(rece, k_row, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(recw, k_row, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	for(i=1;i<k_row-1;i++){
		mydata[i][0] = recw[i];
		mydata[i][k_col-1] = rece[i];
	}

    MPI_Isend(&mydata[1][0],1,MPI_INT, nord, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[1][k_col-1],1,MPI_INT, nord, 1, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[k_row-2][0],1,MPI_INT, sud, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&mydata[k_row-2][k_col-1],1,MPI_INT, sud, 1, MPI_COMM_WORLD, &request);

    MPI_Recv(&mydata[k_row-1][0], 1, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[k_row-1][k_col-1], 1, MPI_INT, sud, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[0][0], 1, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&mydata[0][k_col-1], 1, MPI_INT, nord, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //étape    
    int** voisins = malloc(k_row*sizeof(int*));
		for (i = 0; i<k_row; i++)
			voisins[i] = calloc(k_col,sizeof(int));

    for(i=1; i<k_row-1; i++){
      for(j=1; j<k_col-1; j++){
				voisins[i][j] = mydata[i-1][j-1]+mydata[i-1][j]+mydata[i-1][j+1]+mydata[i][j-1]+mydata[i][j+1]+mydata[i+1][j-1]+mydata[i+1][j]+mydata[i+1][j+1];
			;}
    ;}

    
    for(i=1; i<k_row-1; i++){
      for(j=1; j<k_col-1; j++){
				if(voisins[i][j]==3){mydata[i][j]=1;}
				else if(voisins[i][j]!=2){mydata[i][j]=0;};
      }
    }
		for (i = 0; i<k_row; i++)
			free(voisins[i]);
		free(voisins);

  }
  int alive=0;
  for( i=1; i<k_row-1; i++){
    for(j=1; j<k_col-1; j++){
      if(mydata[i][j]){alive++;};
    }
  }

  int global_alive;
  MPI_Reduce(&alive,&global_alive,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if(world_rank==0){
    printf("global_alive : %d\n", global_alive);

  }
  
  MPI_Finalize();

	return 0;
}
