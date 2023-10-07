#include <stdio.h>
#include <stdlib.h>

int main() {
    int nx = 100;
    float T = 2.5;
    float h = 1.0/nx;
    float mu = 0.25;
    float dt = mu * h * h;  // Pas de temps
    int nt = T/dt;  // Nombre total d'itérations de temps pour t=2.5

    //creer matrice C,un1 est Un+1, et u et se.
    float C[nx][nx];
    float Un1[nx];
    float u[nx];
    float se[nx];
    int i,j;


    //defini la matrice C

    for(i=0;i<nx;i++){ //initialization
        for(j=0;j<nx;j++){
            C[i][j] = 0;
        }
    }

    for (int m = 0; m < nx; m++) {
        C[m][m] = 1-2 * mu;
        if(m==0)
            C[m][m+1] = mu;
        else if(m==nx-1){
            C[m][m-1] = 2*mu;
        }
        else{
            C[m][m+1] = mu;
            C[m][m-1] = mu;
        }
    }

    // matrice se


    for(i=0;i<nx;i++){
        se[i] = 0;
    }
    se[0] = mu;


    //definir u


    for(int i=0;i<nx;i++){
        u[i] = 0;
    }


    //calcule la matrice n+1s
    for (int i = 0; i < nx; i++) {
        Un1[i] = se[i];
        for (int j = 0; j < nx; j++) {
            Un1[i]= Un1[i] + C[i][j]*u[j];
            u[j] = Un1[i];
        }
    }


    //afficher la matrice C

    printf("la matrice C est: \n");

    for(i=0;i<nx;i++){
        for(j=0;j<nx;j++){
            printf("%.2f ",C[i][j]);
        }
        printf("\n");
    }

    //afficher la matrice Un1

    printf("la matrice Un1 est: \n");

    for(i=0;i<nx;i++){
        printf("%.2f ",Un1[i]);
    }
    printf("\n");

    //afficher la matrice u

    printf("la matrice u est: \n");

    for(i=0;i<nx;i++){
        printf("%.2f ",u[i]);
    }
    printf("\n");





    FILE *fp = fopen("resultat.txt", "w");

    for (int i = 0; i < nx; i++) {
        fprintf(fp, "%lf\n", u[i]);

    }
    fclose(fp);

    return 0;

}
