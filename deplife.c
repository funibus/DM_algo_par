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

  int** mydata =(int **) malloc(k_row*sizeof(int*));
  for (i = 0; i<k_row; i++)
      {
	  mydata[i] = (int*) malloc(k_col*sizeof(int));
      }

  srand(time(NULL));
  for(i=0;i<k_row;i++){
    for(j=0;j<k_col;j++){
      if((double)rand()/(double)RAND_MAX<d){mydata[i][j]=1;}
      else{mydata[i][j]=0;}
    }
  }

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
    return (i -1 +kp) % kp +kp*row;
  }

  int nord= fnord(world_rank);
  int sud= fsud(world_rank);
  int est= fest(world_rank);
  int ouest = fouest(world_rank);

  int ne= fnord(fest(world_rank));
  int no=fnord(fouest(world_rank));
  int se=fsud(fest(world_rank));
  int so= fsud(fouest(world_rank));
  //gérer les send/receive

  int* sendo = calloc(k_row, sizeof(int));
  int* sende = calloc(k_row, sizeof(int));
  int* reco = calloc(k_row, sizeof(int));
  int* rece = calloc(k_row, sizeof(int));
  int* recn = calloc(k_col, sizeof(int));
  int* recs = calloc(k_col, sizeof(int));
  int case_ne = 0, case_no = 0, case_se = 0, case_so = 0;
  for(i=0;i<k_row;i++){
    sendo[i]=mydata[i][0];
    sende[i]=mydata[i][k_col-1];
  };

//fprintf(stderr, "proc %d a %d lignes et %d colonnes\n", world_rank,k_row,k_col);

  for( s=1;s<=nb_step;s++){

    MPI_Request request;

    //echange avec les proc est et ouest
    if (col != 0)
    {
        MPI_Isend(sendo, k_row, MPI_INT, ouest,  0, MPI_COMM_WORLD,&request);
        MPI_Request_free(&request);
    }
    if (col != kp-1)
    {
        MPI_Isend(sende, k_row, MPI_INT, est,  0, MPI_COMM_WORLD,&request);
        MPI_Request_free(&request);
    }

    if (col != kp-1)
        MPI_Recv(rece, k_row, MPI_INT, est, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (col != 0)
        MPI_Recv(reco, k_row, MPI_INT, ouest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    //envoie au nord et recoit du sud
    if (row != 0 && s == 1) //a partir de l'etape 2, on envoie au nord des qu'on a fait le calcul (plus loin dans le code)
    {
	if (col != 0)
 	{
            MPI_Isend(&mydata[0][0],1,MPI_INT, no, 0, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
	}
	if (col != kp-1)
	{
            MPI_Isend(&mydata[0][k_col-1],1,MPI_INT, ne, 0, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
	}
        MPI_Isend(mydata[0], k_col, MPI_INT, nord,  0, MPI_COMM_WORLD,&request);
        MPI_Request_free(&request);
    }

    if (row != kp-1)
    {
	if (col != kp-1)
            MPI_Recv(&case_se, 1, MPI_INT, se, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if (col != 0)
            MPI_Recv(&case_so, 1, MPI_INT, so, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Recv(recs, k_col, MPI_INT, sud, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //recoit du nord quand c'est pret
    if (row != 0)
    {
	if (col != kp-1)
            MPI_Recv(&case_ne, 1, MPI_INT, ne, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if (col != 0)
            MPI_Recv(&case_no, 1, MPI_INT, no, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Recv(recn, k_col, MPI_INT, nord, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for(i=0; i<k_row; i++)
    {
        for (j = 0; j<k_col; j++)
        {
            int voisin_no, voisin_n,voisin_ne,voisin_o,voisin_e,voisin_so,voisin_s,voisin_se;
            if (i == 0 && j == 0)
                voisin_no = case_no;
	    else if (i == 0)
		voisin_no = recn[j-1];
	    else if (j == 0)
		voisin_no = reco[i-1];
            else
                voisin_no = mydata[i-1][j-1];

            if (i == 0 && j == k_col-1)
                voisin_ne = case_ne;
	    else if (i == 0)
		voisin_ne = recn[j+1];
	    else if (j == k_col-1)
		voisin_ne = rece[i-1];
            else
                voisin_ne = mydata[i-1][j+1];

            if (i == k_row-1 && j == 0)
                voisin_so = case_so;
	    else if (i == k_row-1)
		voisin_so = recs[j-1];
	    else if (j == 0)
		voisin_so = reco[i+1];
            else
                voisin_so = mydata[i+1][j-1];

            if (i == k_row-1 && j == k_col-1)
                voisin_se = case_se;
	    else if (i == k_row-1)
		voisin_se = recs[j+1];
	    else if (j == k_col-1)
		voisin_se = rece[i+1];
            else
                voisin_se = mydata[i+1][j+1];

            if (i == 0)
                voisin_n = recn[j];
            else
                voisin_n = mydata[i-1][j];

            if (i == k_row-1)
                voisin_s = recs[j];
            else
                voisin_s = mydata[i+1][j];

            if (j == 0)
                voisin_o = reco[i];
            else
                voisin_o = mydata[i][j-1];

            if (j == k_col-1)
                voisin_e = rece[i];
            else
                voisin_e = mydata[i][j+1];

            //mise a jour de la case i,j
            int nb_voisins = voisin_no+voisin_n+voisin_ne+voisin_o+voisin_e+voisin_so+voisin_s+voisin_se;

            if (mydata[i][j] == 1 && (nb_voisins <2 || nb_voisins >3) )
            {
                mydata[i][j] = 0;
            }
            else if (mydata[i][j] == 0 && nb_voisins == 3)
            {
                mydata[i][j] = 1;
            }

        }

        if (row != 0 && i == 0 && s < nb_step) //si on est pas dans la premiere ligne,
        {                       //on envoie notre premiere ligne au proc du dessus des qu'elle est calculee
	    if (col != 0)
	    {
                MPI_Isend(&mydata[0][0],1,MPI_INT, no, 0, MPI_COMM_WORLD, &request);
                MPI_Request_free(&request);
	    }
	    if (col != kp-1)
	    {
                MPI_Isend(&mydata[0][k_col-1],1,MPI_INT, ne, 0, MPI_COMM_WORLD, &request);
                MPI_Request_free(&request);
	    }
            MPI_Isend(mydata[0], k_col, MPI_INT, nord,  0, MPI_COMM_WORLD,&request);
            MPI_Request_free(&request);
        }
    }

    if (row != kp-1)
    {
	if (col != kp-1)
	{
            MPI_Isend(&mydata[k_row-1][k_col-1],1,MPI_INT, se, 0, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
	}
	if (col != 0)
	{
            MPI_Isend(&mydata[k_row-1][0],1,MPI_INT, so, 0, MPI_COMM_WORLD, &request);
            MPI_Request_free(&request);
	}
        MPI_Isend(mydata[k_row-1], k_col, MPI_INT, sud,  0, MPI_COMM_WORLD,&request);
        MPI_Request_free(&request);
    }
  }


  int alive=0;
  for( i=0; i<k_row; i++)
  {
    for(j=0; j<k_col; j++)
    {
      if(mydata[i][j] == 1)
        alive++;
    }
  }

  int global_alive;
  MPI_Reduce(&alive,&global_alive,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if(world_rank==0){
    printf("global_alive : %d\n", global_alive);
    printf("densité : %f\n", (double)global_alive/(double)(n*n));

  }

  MPI_Finalize();
}
