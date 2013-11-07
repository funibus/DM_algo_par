#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<mpi.h>

int max(int a, int b)
{
    if (a < b)
        return b;
    else
        return a;
}

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        printf("il n'y a pas le bon nombre d'arguments.\n");
        return 1;
    }

    int rule = atoi(argv[1]);
    int nb_steps = atoi(argv[2]);
    char* input = argv[3];
    int i,j;
    int total;

    /*on ecrit la loi sous forme d'un tableau (plus facile a utiliser*/
    int tab_loi[8];
    for (i = 0; i<8; i++)
    {
        if (rule%2 == 0)
        {
            tab_loi[i] = 0;
        }
        else
        {
            tab_loi[i] = 1;
        }
        rule = rule/2;
    }

    /*on convertit l'entree en tableau d'entiers*/
    int taille_entree = strlen(input);
    int* entree = malloc(sizeof(int)*taille_entree);
    for (i = 0; i<taille_entree; i++)
    {
        entree[i] = (int) (input[i]-'0');
    }

    int taille_finale = taille_entree + 2*nb_steps;

    /*debut du programme parallele*/
    MPI_Init(NULL, NULL);
    int nb_proc;
    int mon_rang;
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &mon_rang);

    /* creation de tableaux de la taille voulue pour stocker un bout du grand tableau dans chaque proc*/
    int taille_blocs = max(taille_entree, taille_finale/nb_proc +3);
    int* bloc = calloc(taille_blocs, sizeof(int));


    /*repartition de l'entree sur le(s) proc(s) du milieu*/
    if (nb_proc%2 == 0) /*si on a un nombre pair de proc, l'entree se partage entre les 2 procs du milieu*/
    {
        if (mon_rang == nb_proc/2-1)
        {
            for (i = 0; i<taille_entree/2; i++)
            {
                bloc[taille_blocs-1-i] = entree[taille_entree/2-i-1];
            }
        }
        if (mon_rang == nb_proc/2)
        {
            for (i = 0; i<taille_entree-taille_entree/2; i++)
            {
                bloc[i] = entree[taille_entree/2+i];
            }
        }
    }

    if (nb_proc%2 == 1) /*l'entree n'est que sur un proc*/
    {
	if (mon_rang == nb_proc/2)
	{
            int debut = (taille_blocs-taille_entree)/2;
            for (i = 0; i<taille_entree; i++)
            {
                bloc[debut+i] = entree[i];
      	    }
	}
    }
/*fprintf(stderr, "proc %d, bloc ini de taille %d : ", mon_rang, taille_blocs);
for (i = 0; i<taille_blocs; i++)
{
fprintf(stderr, "%d, ", bloc[i]);
}
fprintf(stderr, "\n");*/

    /*calcul de la configuration suivante (en boucle)*/

    int* temp = malloc(taille_blocs*sizeof(int));
    int* echange;
    int case_avant;
    int case_apres;
    MPI_Request haut;
    MPI_Request bas;
    for (j = 0; j<nb_steps; j++)
    {
        /*on recupere la case du dessus et celle du dessous pour pouvoir
        appliquer un pas de calcul au block entier*/
        if (mon_rang == 0)
        {
            case_avant = 0;
            MPI_Isend(&bloc[taille_blocs-1], 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &haut);
            MPI_Recv(&case_apres, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Wait(&haut, MPI_STATUS_IGNORE);
        }

        else if (mon_rang == nb_proc-1)
        {
            case_apres = 0;
            MPI_Recv(&case_avant, 1, MPI_INT, mon_rang-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Isend(&bloc[0], 1, MPI_INT, mon_rang-1, 0, MPI_COMM_WORLD, &bas);
            MPI_Wait(&bas, MPI_STATUS_IGNORE);
        }

        else
        {
            MPI_Isend(&bloc[taille_blocs-1], 1, MPI_INT, mon_rang+1, 0, MPI_COMM_WORLD, &haut);
            MPI_Recv(&case_avant, 1, MPI_INT, mon_rang-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Isend(&bloc[0], 1, MPI_INT, mon_rang-1, 0, MPI_COMM_WORLD, &bas);
            MPI_Recv(&case_apres, 1, MPI_INT, mon_rang+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&haut, MPI_STATUS_IGNORE);
            MPI_Wait(&bas, MPI_STATUS_IGNORE);
        }
        /*une etape*/
        total = 4*case_avant+2*bloc[0]+bloc[1];
        temp[0] = tab_loi[total];
        total = 4*bloc[taille_blocs-2]+2*bloc[taille_blocs-1]+case_apres;
        temp[taille_blocs-1] = tab_loi[total];

        for (i = 1; i<taille_blocs-1; i++)
        {
            total = 4*bloc[i-1]+2*bloc[i]+bloc[i+1];
            temp[i] = tab_loi[total];
        }
        echange = temp;
        temp = bloc;
        bloc = echange;

    }
/*fprintf(stderr, "proc %d : ma liste : ", mon_rang);
for (i = 0; i<taille_blocs; i++)
{
fprintf(stderr, "%d, ", bloc[i]);
} 
fprintf(stderr, "\n");*/

    /*calcul le nombre de '1' par proc*/
    int nb_1 = 0;
    for (i = 0; i<taille_blocs; i++)
    {
        if (bloc[i] == 1)
            nb_1++;
    }
    /*nombre total de '1'*/
    int total_1 = 0;
    MPI_Reduce(&nb_1, &total_1, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mon_rang == 0)
    {
        printf("le nombre de 1 apres %d etapes est %d\n", nb_steps, total_1);
    }
    MPI_Finalize();
    return 0;
}
