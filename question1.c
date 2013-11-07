#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<string.h>

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


    /*calcul du nouveau tableau a chaque etape*/
    int* temp = malloc(sizeof(int)*taille_entree);
    int* echange;
    for (j = 0; j<nb_steps; j++)
    {
        for (i = 0; i<taille_entree; i++)
        {
            total = 4*entree[(i-1+taille_entree)%taille_entree]+2*entree[i]+entree[(i+1)%taille_entree];
            temp[i] = tab_loi[total];
        }
        echange = temp;
        temp = entree;
        entree = echange;
    }

    /*impression du resultat*/
    for (i = 0; i<taille_entree-1; i++)
    {
        printf("%d, ", entree[i]);
    }
    printf("%d\n", entree[taille_entree-1]);

    return 0;
}
