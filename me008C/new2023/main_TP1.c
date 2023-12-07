#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 700

void remplir_vec_U0(double U0[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        U0[i] = 0;
    }
}
void remplir_mat_C(double C[MAX][MAX], double mu, int nx)
{
    // premiere ligne
    {
        int i = 0;
        C[i][i] = 1 - 2 * mu;
        C[i][i + 1] = mu;
    }

    // l'interieur
    for (int i = 1; i < nx - 1; i++)
    {
        C[i][i - 1] = mu;
        C[i][i] = 1 - 2 * mu;
        C[i][i + 1] = mu;
    }

    // derniere ligne
    {
        int i = nx - 1;
        C[i][i - 1] = 2 * mu;
        C[i][i] = 1 - 2 * mu;
    }
}
void remplir_vec_se(double se[MAX], double mu, int nx)
{
    se[0] = mu;
    for (int i = 1; i < nx; i++)
    {
        se[i] = 0;
    }
}
void remplir_vec_Un_avec_Unplus1(double Un[MAX], double Unplus1[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        Un[i] = Unplus1[i];
    }
}
void iteration_euler_explicite(double C[MAX][MAX], double Un[MAX], double se[MAX], double Unplus1[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        Unplus1[i] = 0;
        for (int j = 0; j < nx; j++)
        {
            Unplus1[i] += C[i][j] * Un[j];
        }
        Unplus1[i] += se[i];
    }
}

void afficher_vect(double x[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
        if (nx < 10)
            printf("%4.4f\t", x[i]);
    if (nx < 10)
        printf("\n");
}
void afficher_mat(double A[MAX][MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            if (nx < 10)
                printf("%.4f\t", A[i][j]);
        }
        if (nx < 10)
            printf("\n");
    }
}
void sauvgarder_vect(double x[MAX], int nx, char nom[MAX])
{
    FILE *fptr;
    fptr = fopen(nom, "w");
    for (int i = 0; i < nx; i++)
    {
        fprintf(fptr, "%1.17e\n", x[i]);
    }

    fclose(fptr);
    printf("Vecteur sauvgarde sous:%s\n", nom);
}
void sauvgarder_vect_et_maillage_et_bords(double x[MAX], int nx, double h, char nom[MAX])
{
    FILE *fptr;
    fptr = fopen(nom, "w");
    fprintf(fptr, "%1.17e\t", 0.0);
    fprintf(fptr, "%1.17e\n", 1.0);
    for (int i = 0; i < nx; i++)
    {
        fprintf(fptr, "%1.17e\t", (i + 1) * h);
        fprintf(fptr, "%1.17e\n", x[i]);
    }

    fclose(fptr);
    printf("Vecteur sauvgarde sous:%s\n", nom);
}
void sauvgarder_valeur_instantane(double x[MAX], int nx, double t, char nom[MAX])
{
    FILE *fptr;
    fptr = fopen(nom, "a");
    fprintf(fptr, "%1.17e\t", t);
    // int indices[] = {10, 20, 25, 50, 75};
    int indices[] = {19, 39, 59, 79, 99};

    for (int i = 0; i < 5; i++)
    {
        fprintf(fptr, "%1.17e\t", x[indices[i]]);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
}

int main()
{
    if (0)
    {
        // Parametres
        int nx = 100;
        double h = 1.0 / nx;
        double dt = 0.25 * h * h;
        double mu = dt / h / h;
        double Tmax = 2.0;
        int nt = Tmax / dt;
        char nom_de_sauvgarde[60] = "u1u.dat";
        // nt = 200;

        double Un[MAX];     // Condition initial
        double se[MAX];     // vecteur se
        double C[MAX][MAX]; // matrice C
        double Unplus1[MAX];

        printf("Parametres:\t nx: %d\t h: %lf\t dt: %lf\t mu: %lf\t Tmax: %lf\t nt: %d\n",
               nx, h, dt, mu, Tmax, nt);

        remplir_vec_U0(Un, nx);
        printf("Vecteur U0:\n");
        // afficher_vect(Un, nx);

        printf("Vecteur se:\n");
        remplir_vec_se(se, mu, nx);
        // afficher_vect(se, nx);

        printf("Matrice C:\n");
        remplir_mat_C(C, mu, nx);
        afficher_mat(C, nx);

        // Integration temporelle
        // nt = 3;
        for (int it = 0; it < nt; it++)
        {
            iteration_euler_explicite(C, Un, se, Unplus1, nx);
            remplir_vec_Un_avec_Unplus1(Un, Unplus1, nx);
            if (nt > 10 && (it + 1) % (nt / 10) == 0)
                printf("EE iter: %d temps: %lf dt: %lf\n", it + 1, (it + 1) * dt, dt);

            sauvgarder_valeur_instantane(Un, nx, (it + 1) * dt, "valeurs_instantanes.dat");
        }

        sauvgarder_vect_et_maillage_et_bords(Un, nx, h, nom_de_sauvgarde);

        printf("Vecteur Un:\n");
        afficher_vect(Un, nx);

        // afficher_vect(Un, nx);
    }

    if (1)
    {
        // Parametres
        int nx = 100;
        double h = 1.0 / nx;
        double dt = 0.25 * h * h;
        double mu = dt / h / h;
        char nom_de_sauvgarde[60] = "u1u.dat";
        // nt = 200;

        double Tmax = 0.001;
        int nt = Tmax / dt;

        double Un[MAX];     // Condition initial
        double se[MAX];     // vecteur se
        double C[MAX][MAX]; // matrice C
        double Unplus1[MAX];

        remplir_vec_U0(Un, nx);
        remplir_vec_se(se, mu, nx);
        remplir_mat_C(C, mu, nx);

        // Integration temporelle
        for (int it = 0; it < nt; it++)
        {

            for (int i = 0; i < nx; i++)
            {
                Unplus1[i] = 0;
                for (int j = 0; j < nx; j++)
                {
                    Unplus1[i] += C[i][j] * Un[j];
                }
                Unplus1[i] += se[i];
            }

            remplir_vec_Un_avec_Unplus1(Un, Unplus1, nx);

            {
                FILE *fptr;
                fptr = fopen("valeurs_instantanees.dat", "a");
                fprintf(fptr, "%1.17e\t%1.17e\n", (it + 1) * dt, Un[nx*2/10-1]);
                fclose(fptr);
            }
        }

        {
            FILE *fptr;
            fptr = fopen("u.dat", "w");
            fprintf(fptr, "%1.17e\t%1.17e\n", 0.0,1.0);
            for (int i = 0; i < nx; i++)
                fprintf(fptr, "%1.17e\t%1.17e\n", (i + 1) * h, Un[i]);
            fclose(fptr);
        }
    }

    return 0;
}