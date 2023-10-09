#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 400

void factoriser_LU(const double A[MAX][MAX], double L[MAX][MAX], double U[MAX][MAX], int n)
{
    // a remplir
}

// descente
void resol_trig_inf(double A[MAX][MAX], double x[MAX], double b[MAX], int n)
{
    // a remplir
}
void resol_trig_inf_tridiag(double l[MAX], double x[MAX], double b[MAX], int n)
{
}

// remontée
void resol_trig_sup(double A[MAX][MAX], double x[MAX], double b[MAX], int n)
{
    // a remplir
}
void resol_trig_sup_tridiag(double u[MAX], double v[MAX], double x[MAX], double b[MAX], int n)
{
    // a remplir
}

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
void remplir_mat_E(double E[MAX][MAX], double mu, int nx)
{
    // premiere ligne
    {
        int i = 0;
        E[i][i] = 1 + 2 * mu;
        E[i][i + 1] = -mu;
    }

    // l'interieur
    for (int i = 1; i < nx - 1; i++)
    {
        E[i][i - 1] = -mu;
        E[i][i] = 1 + 2 * mu;
        E[i][i + 1] = -mu;
    }

    // derniere ligne
    {
        int i = nx - 1;
        E[i][i - 1] = -2 * mu;
        E[i][i] = 1 + 2 * mu;
    }
}
void remplir_mat_A_TD3(double C[MAX][MAX], int nx)
{
    // premiere ligne
    {
        int i = 0;
        C[i][i] = 2.0;
        C[i][i + 1] = -1.0;
    }

    // l'interieur
    for (int i = 1; i < nx - 1; i++)
    {
        C[i][i - 1] = -1.0;
        C[i][i] = 2.0;
        C[i][i + 1] = -1.0;
    }

    // derniere ligne
    {
        int i = nx - 1;
        C[i][i - 1] = -1.0;
        C[i][i] = 2.0;
    }
}

void remplir_vec_b_TD3(double b[MAX], int nx)
{
    b[0] = 1.0;
    for (int i = 1; i < nx; i++)
    {
        b[i] = 0.0;
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
void remplir_vec_si(double se[MAX], double mu, int nx)
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

void remplir_vec_adc(double a[MAX], double d[MAX], double c[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        a[i] = -1.0;
        d[i] = 2.0;
        c[i] = -1.0;
    }
}

void remplir_vec_adc_E(double a[MAX], double d[MAX], double c[MAX], double mu, int nx)
{

    // premiere ligne
    {
        int i = 0;
        d[i] = 1 + 2 * mu;
        c[i] = -mu;
    }

    // l'interieur
    for (int i = 1; i < nx - 1; i++)
    {
        a[i] = -mu;
        d[i] = 1 + 2 * mu;
        c[i] = -mu;
    }

    // derniere ligne
    {
        int i = nx - 1;
        a[i] = -2 * mu;
        d[i] = 1 + 2 * mu;
    }
}

void factoriser_tridiag_luv(double a[MAX], double d[MAX], double c[MAX],
                            double l[MAX], double u[MAX], double v[MAX], int nx)
{
    // a remplir
}
void iteration_euler_explicite(double C[MAX][MAX], double Un[MAX], double se[MAX], double Unplus1[MAX], int nx)
{
    // a remplir
}

void iteration_euler_implicite(double L[MAX][MAX], double U[MAX][MAX], double Un[MAX],
                               double si[MAX], double Unplus1[MAX], int nx)
{
    // a remplir
}
void iteration_euler_implicite_tridiag(double l[MAX], double u[MAX], double v[MAX], double Un[MAX],
                                       double si[MAX], double Unplus1[MAX], int nx)
{
    // a remplir
}

void afficher_vect(double x[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
        printf("%4.2f\t", x[i]);
    printf("\n");
}
void afficher_mat(double A[MAX][MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            printf("%.4f\t", A[i][j]);
        }
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
    int indices[] = {10, 20, 25, 50, 75};

    for (int i = 0; i < 5; i++)
    {
        fprintf(fptr, "%1.17e\t", x[indices[i]]);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
}

int main()
{
    // Parametres
    int nx = 5;
    double h = 1.0 / nx;
    double dt = 2 * h * h;
    double mu = dt / h / h;
    double Tmax = 4;
    int nt = Tmax / dt;
    char nom_de_sauvgarde[60] = "u4i.dat";
    // nt = 200;

    double Un[MAX];                                            // Condition initial
    double si[MAX], b[MAX], x[MAX], y[MAX];                    // vecteurs si,b,x,y
    double A[MAX][MAX], L[MAX][MAX], U[MAX][MAX], E[MAX][MAX]; // matrice A,L,U,E
    double Unplus1[MAX];

    // Exercice 1
    printf("Matrice A:\n");
    remplir_mat_A_TD3(A, nx);
    afficher_mat(A, nx);

    printf("Matrice L:\n");
    factoriser_LU(A, L, U, nx);
    afficher_mat(L, nx);
    printf("Matrice U:\n");
    afficher_mat(U, nx);

    // Exercice 2
    printf("Vecteur b:\n");
    remplir_vec_b_TD3(b, nx);
    afficher_vect(b, nx);

    printf("Vecteur y:\n");
    resol_trig_inf(L, y, b, nx);
    afficher_vect(y, nx);

    printf("Vecteur x:\n");
    resol_trig_sup(U, x, y, nx);
    afficher_vect(x, nx);

    // Exercice 3 - integration temporelle
    nt = 0;
    printf("Parametres:\t nx: %d\t h: %lf\t dt: %lf\t mu: %lf\t Tmax: %lf\t nt: %d\n",
           nx, h, dt, mu, Tmax, nt);

    // factorisation
    printf("Matrice E:\n");
    remplir_mat_E(E, mu, nx);
    remplir_vec_si(si, mu, nx);
    // afficher_mat(E, nx);
    // afficher_vect(si, nx);
    factoriser_LU(E, L, U, nx);

    for (int it = 0; it < nt; it++)
    {
        iteration_euler_implicite(L, U, Un, si, Unplus1, nx);
        remplir_vec_Un_avec_Unplus1(Un, Unplus1, nx);
        if (nt > 10 && (it + 1) % (nt / 10) == 0)
            printf("IE iter: %d temps: %lf dt: %lf\n", it + 1, (it + 1) * dt, dt);
        // afficher_vect(Unplus1, nx);

        // sauvgarder_valeur_instantane(Un, nx, (it + 1) * dt, "valeurs_instantanes.dat");
    }

    sauvgarder_vect_et_maillage_et_bords(Un, nx, h, nom_de_sauvgarde);

    // Exercice 4
    double a[MAX], d[MAX], c[MAX], u[MAX], v[MAX], l[MAX]; // a,b,c; u,v,l

    remplir_vec_adc(a, d, c, nx);
    factoriser_tridiag_luv(a, d, c, l, u, v, nx);
    // printf("Vecteur a:\n");
    // afficher_vect(a, nx);
    // printf("Vecteur d:\n");
    // afficher_vect(d, nx);
    // printf("Vecteur c:\n");
    // afficher_vect(c, nx);
    // printf("Vecteur l:\n");
    // afficher_vect(l, nx);
    // printf("Vecteur u:\n");
    // afficher_vect(u, nx);
    // printf("Vecteur v:\n");
    // afficher_vect(v, nx);

    // printf("Vecteur b:\n");
    // afficher_vect(b, nx);

    // printf("Vecteur y:\n");
    // afficher_vect(y, nx);

    // printf("Vecteur x:\n");
    // afficher_vect(x, nx);

    resol_trig_inf_tridiag(l, y, b, nx);
    resol_trig_sup_tridiag(u, v, x, y, nx);

    // printf("Vecteur y:\n");
    // afficher_vect(y, nx);

    // printf("Vecteur x:\n");
    // afficher_vect(x, nx);

    // Integration temporelle

    remplir_vec_adc_E(a, d, c, mu, nx);
    factoriser_tridiag_luv(a, d, c, l, u, v, nx);

    double U0[MAX];
    remplir_vec_Un_avec_Unplus1(Un, U0, nx);

    for (int it = 0; it < nt; it++)
    {
        // iteration_euler_implicite(L, U, Un, si, Unplus1, nx);
        iteration_euler_implicite_tridiag(l, u, v, Un, si, Unplus1, nx);
        remplir_vec_Un_avec_Unplus1(Un, Unplus1, nx);
        if (nt > 10 && (it + 1) % (nt / 10) == 0)
            printf("IE iter: %d temps: %lf dt: %lf\n", it + 1, (it + 1) * dt, dt);
        // afficher_vect(Unplus1, nx);

        // sauvgarder_valeur_instantane(Un, nx, (it + 1) * dt, "valeurs_instantanes.dat");
    }

    char nom_de_sauvgarde2[60] = "u4_tridiag.dat";
    sauvgarder_vect_et_maillage_et_bords(Un, nx, h, nom_de_sauvgarde2);

    return 0;
}