#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 100

void factoriser_LU(const double A[MAX][MAX], double L[MAX][MAX], double U[MAX][MAX], int n)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < j + 1; i++)
        {
            U[i][j] = A[i][j];
            // printf("%4.2f \t ",U[i][j]);
            // printf("%4.2f \n ",A[i][j]);
            for (int k = 0; k < i; k++)
            {
                U[i][j] -= L[i][k] * U[k][j];
            } // k
        }     // i
        L[j][j] = 1;
        for (int i = j + 1; i < n; i++)
        {
            L[i][j] = A[i][j];
            // printf("%4.2f \n ",L[i][j]);

            for (int k = 0; k < j; k++)
            {
                L[i][j] -= L[i][k] * U[k][j];
            } // k
            L[i][j] = L[i][j] / U[j][j];
        } // i
    }     // j
}

// descente
void resol_trig_inf(double A[MAX][MAX], double x[MAX], double b[MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        x[i] = b[i];
        for (int j = 0; j < i; j++)
        {
            x[i] -= x[j] * A[i][j];
        }
        x[i] = x[i] / A[i][i];
    }
}
void resol_trig_inf_tridiag(double l[MAX], double x[MAX], double b[MAX], int n)
{
    {
        int i = 0;
        x[i] = b[i];
    }
    for (int i = 1; i < n; i++)
        x[i] = b[i] - l[i] * x[i - 1];
}

// remontée
void resol_trig_sup(double A[MAX][MAX], double x[MAX], double b[MAX], int n)
{
    for (int i = n - 1; i > -1; i--)
    {
        x[i] = b[i];
        for (int j = n - 1; j > i; j--)
        {
            x[i] -= x[j] * A[i][j];
        }
        x[i] = x[i] / A[i][i];
    }
}
void resol_trig_sup_tridiag(double u[MAX], double v[MAX], double x[MAX], double b[MAX], int n)
{
    {
        int i = n - 1;
        x[i] = b[i] / u[i];
    }
    for (int i = n - 2; i > -1; i--)
        x[i] = (b[i] - x[i + 1] * v[i]) / u[i];
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

void remplir_Ab_TP3(double C[MAX][MAX], double b[MAX], double h, int nx)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            C[i][j] = 0;
        }
        b[i] = 0;
    }

    // // premiere ligne
    // {
    //     int i = 0;
    //     C[i][i] = 2.0;
    //     C[i][i + 1] = -1.0;
    //     // C[i][nx - 1] = -1; // not the best idea
    // }

    // premiere ligne
    {
        int i = 0;
        C[i][i] = 4.0;
        C[i][i + 1] = 1;
        for (int j = 2; j < nx - 1; j++)
        {
            C[i][j] = 2.0;
        }
        C[i][nx - 1] = 1;

        // for (int j = 0; j < nx; j++)
        // {
        //     C[i][j] = 2.0;
        // }
        // C[i][nx - 1] = 0.0; // not the best idea
        // C[i][nx - 1] = -1; // not the best idea
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
        C[i][0] = -1.0;
        C[i][i - 1] = -1.0;
        C[i][i] = 2.0;
    }
    for (int i = 0; i < nx; i++)
    {
        double x = (i + 1) * h;
        b[i] = h * h * sin(2 * M_PI * x);
    }
}
void remplir_E_si_TP3(double A[MAX][MAX], double b[MAX], double mu, double E[MAX][MAX], double si[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        si[i] = b[i] * mu;
        for (int j = 0; j < nx; j++)
        {
            E[i][j] = A[i][j] * mu;
            if (i == j)
            {
                E[i][j] += 1.0;
            }
        }
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
    {
        int i = 0;
        u[i] = d[i];
        v[i] = c[i];
    }
    for (int i = 1; i < nx - 1; i++)
    {
        l[i] = a[i] / u[i - 1];
        u[i] = d[i] - v[i - 1] * l[i];
        v[i] = c[i];
    }
    {
        int i = nx - 1;
        l[i] = a[i] / u[i - 1];
        u[i] = d[i] - v[i - 1] * l[i];
    }
}

void remplir_Ab_TD3_1(double A[MAX][MAX], double b[MAX])
{
    double C[MAX][MAX] = {{4, -1, 0}, {-1, 4, -1}, {0, -2, 4}};
    double d[MAX] = {3.0 / 16, 4.0 / 16, 6.0 / 16};
    for (int i = 0; i < 3; i++)
    {
        b[i] = d[i];
        for (int j = 0; j < 3; j++)
            A[i][j] = C[i][j];
    }
}
void remplir_Ab_TD3_2(double A[MAX][MAX], double b[MAX])
{
    double d[MAX] = {12, 26, 60};
    for (int i = 0; i < 3; i++)
    {
        b[i] = d[i];
        for (int j = 0; j < 3; j++)
            A[i][j] = pow(j + 1, i + 1);
    }
}
void remplir_MN_Jacobi(double A[MAX][MAX], double M[MAX][MAX], double N[MAX][MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                M[i][j] = A[i][j];
            else
                N[i][j] = -A[i][j];
        }
    }
}
void remplir_MN_Gauss_Seidel(double A[MAX][MAX], double M[MAX][MAX], double N[MAX][MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i + 1; j++)
        {
            M[i][j] = A[i][j];
            // N=
        }
    }
}
void remplir_MN_SOR(double A[MAX][MAX], double M[MAX][MAX], double N[MAX][MAX], double omega, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            M[i][j] = 0.0;
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i + 1; j++)
        {
            if (i == j)

                M[i][j] = A[i][j] / (omega);

            else
                M[i][j] = A[i][j];

            // N=
        }
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

void iteration_euler_implicite(double L[MAX][MAX], double U[MAX][MAX], double Un[MAX],
                               double si[MAX], double Unplus1[MAX], int nx)
{
    double b[MAX], x[MAX], y[MAX];
    for (int i = 0; i < nx; i++)
    {
        b[i] = Un[i] + si[i];
    }

    resol_trig_inf(L, y, b, nx);
    resol_trig_sup(U, Unplus1, y, nx);
}
void iteration_euler_implicite_tridiag(double l[MAX], double u[MAX], double v[MAX], double Un[MAX],
                                       double si[MAX], double Unplus1[MAX], int nx)
{
    double b[MAX], x[MAX], y[MAX];
    for (int i = 0; i < nx; i++)
    {
        b[i] = Un[i] + si[i];
    }

    resol_trig_inf_tridiag(l, y, b, nx);
    resol_trig_sup_tridiag(u, v, Unplus1, y, nx);
}

void residu(double A[MAX][MAX], double b[MAX], double x[MAX], double r[MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        r[i] = b[i];
        for (int j = 0; j < n; j++)
        {
            r[i] -= A[i][j] * x[j];
        }
    }
}
void s_gradient(double A[MAX][MAX], double r[MAX], double s[MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        s[i] = 0;
        for (int j = 0; j < n; j++)
        {
            s[i] += A[i][j] * r[j];
        }
    }
}
double produit_scalaire(double x[MAX], double y[MAX], int n)
{
    double produit = 0;
    for (int i = 0; i < n; i++)
    {
        produit += x[i] * y[i];
    }
    return produit;
}
void ukp1_iter(double u[MAX], double du[MAX], double ukp1[MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        ukp1[i] = u[i] + du[i];
    }
}
void xkp1_gradient(double xk[MAX], double rk[MAX], double alphak, double xkp1[MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        xkp1[i] = xk[i] + rk[i] * alphak;
    }
}
void rkp1_gradient(double rk[MAX], double sk[MAX], double alphak, double rkp1[MAX], int n)
{
    for (int i = 0; i < n; i++)
    {
        rkp1[i] = rk[i] - sk[i] * alphak;
    }
}

double vecteur_L2_norm(double x[MAX], int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += x[i] * x[i];
    }
    return sqrt(norm);
}
void afficher_vect(double x[MAX], int nx)
{
    for (int i = 0; i < nx; i++)
        if (nx < 30)
            printf("%4.4f\t", x[i]);
    if (nx < 30)
        printf("\n");
}
void afficher_mat(double A[MAX][MAX], int nx)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            if (nx < 30)
                printf("%.4f\t", A[i][j]);
        }
        if (nx < 30)
            printf("\n");
    }
}
void initialiser_vecteur(double x[MAX], int nx, double valeur)
{
    for (int i = 0; i < nx; i++)
    {
        x[i] = valeur;
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
void sauvgarder_vect_et_maillage(double x[MAX], int nx, double h, char nom[MAX])
{
    FILE *fptr;
    fptr = fopen(nom, "w");
    // fprintf(fptr, "%1.17e\t", 0.0);
    // fprintf(fptr, "%1.17e\n", 1.0);
    for (int i = 0; i < nx; i++)
    {
        fprintf(fptr, "%1.17e\t", (i + 1) * h);
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
    int indices[] = {6};

    for (int i = 0; i < 1; i++)
    {
        fprintf(fptr, "%1.17e\t", x[indices[i]]);
    }
    fprintf(fptr, "\n");
    fclose(fptr);
}
void TD3_gradient()
{
    // TD3 3.3
    int nx = 3;
    double A[MAX][MAX], b[MAX], rk[MAX], sk[MAX], alphak, xk[MAX];
    printf("Matrice A:\n");
    remplir_mat_A_TD3(A, nx);
    afficher_mat(A, nx);
    b[2] = 4;
    afficher_vect(b, nx);
    residu(A, b, xk, rk, nx);

    for (int i = 0; i < 3; i++)
    {
        s_gradient(A, rk, sk, nx);
        // afficher_vect(sk, nx);
        alphak = produit_scalaire(rk, rk, nx) / produit_scalaire(rk, sk, nx);
        printf("alphak:%lf\n", alphak);
        xkp1_gradient(xk, rk, alphak, xk, nx);
        rkp1_gradient(rk, sk, alphak, rk, nx);
        afficher_vect(xk, nx);
    }
}
void TD3_conj_gradient()
{
    // TD3 3.3
    int nx = 3;
    double A[MAX][MAX], b[MAX], rk[MAX], sk[MAX], alphak, xk[MAX];
    double pk[MAX], qk[MAX], beta, rkp1[MAX];
    printf("Matrice A:\n");
    remplir_mat_A_TD3(A, nx);
    afficher_mat(A, nx);
    b[2] = 4;
    printf("Vecteur b:\n");
    afficher_vect(b, nx);
    residu(A, b, xk, rk, nx);
    residu(A, b, xk, pk, nx);

    for (int i = 0; i < 3; i++)
    {
        printf("======= iter:\t%d\n", i);
        s_gradient(A, pk, qk, nx);
        // afficher_vect(sk, nx);
        alphak = produit_scalaire(rk, rk, nx) / produit_scalaire(pk, qk, nx);
        printf("alphak:%lf\n", alphak);
        xkp1_gradient(xk, pk, alphak, xk, nx);
        rkp1_gradient(rk, qk, alphak, rkp1, nx);
        afficher_vect(xk, nx);
        afficher_vect(rkp1, nx);

        // conj
        beta = produit_scalaire(rkp1, rkp1, nx) / produit_scalaire(rk, rk, nx);
        printf("beta:%lf\n", beta);
        xkp1_gradient(rkp1, pk, beta, pk, nx);
        afficher_vect(pk, nx);
        remplir_vec_Un_avec_Unplus1(rk, rkp1, nx);
    }
}

int main()
{
    // 1
    int a = 4;
    int b = 7;
    int c = 2;

    // Parametres
    int nx = 10 + c;
    double h = 1.0 / nx;
    double dt = 0.3 * h * h;
    double mu = dt / h / h;
    double Tmax = a / 10.0;
    int nt = Tmax / dt;
    char nom_de_sauvgarde[60] = "u_res.dat";
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

    // 2
    double L[MAX][MAX], U[MAX][MAX]; // matrice C
    factoriser_LU(C, L, U, nx);
    printf("L:\n");
    afficher_mat(L, nx);

    // 3
    {
         // Parametres
    int nx = 10+c;
    double h = 1.0 / nx;
    double dt = 1.0 * h * h;
    double mu = dt / h / h;
    double Tmax = a/10.0;
    int nt = Tmax / dt;
    char nom_de_sauvgarde[60] = "u_res_impl.dat";
    // nt = 200;

    double Un[MAX];                                            // Condition initial
    double si[MAX], b[MAX], x[MAX], y[MAX];                    // vecteurs si,b,x,y
    double A[MAX][MAX], L[MAX][MAX], U[MAX][MAX], E[MAX][MAX]; // matrice A,L,U,E
    double Unplus1[MAX];

    // Exercice 1
    printf("Matrice A:\n");
    remplir_mat_A_TD3(A, nx);
    // afficher_mat(A, nx);

    printf("Matrice L:\n");
    factoriser_LU(A, L, U, nx);
    // afficher_mat(L, nx);
    printf("Matrice U:\n");
    // afficher_mat(U, nx);

    // Exercice 2
    printf("Vecteur b:\n");
    remplir_vec_b_TD3(b, nx);
    // afficher_vect(b, nx);

    printf("Vecteur y:\n");
    resol_trig_inf(L, y, b, nx);
    // afficher_vect(y, nx);

    printf("Vecteur x:\n");
    resol_trig_sup(U, x, y, nx);
    // afficher_vect(x, nx);

    // Exercice 3 - integration temporelle
    // nt = 3;
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

        sauvgarder_valeur_instantane(Un, nx, (it + 1) * dt, "valeurs_instantanes_impl.dat");
    }

    sauvgarder_vect_et_maillage_et_bords(Un, nx, h, nom_de_sauvgarde);
    }


    // afficher_vect(Un, nx);
    // if (0)
    // { // TP3
    //     double A[MAX][MAX], b[MAX];
    //     int n = 3;
    //     remplir_Ab_TD3_1(A, b);
    //     // remplir_Ab_TD3_2(A, b); // omega opti approx 1.5
    //     printf("Matrice A:\n");
    //     afficher_mat(A, n);
    //     printf("Vecteur b:\n");
    //     afficher_vect(b, n);

    //     // Jacobi, GS, SOR
    //     double M[MAX][MAX], N[MAX][MAX], r[MAX], u[MAX], du[MAX];
    //     // remplir_MN_Jacobi(A, M, N, n);
    //     // remplir_MN_Gauss_Seidel(A, M, N, n);
    //     double omega_sor = 1.07;
    //     remplir_MN_SOR(A, M, N, omega_sor, n);
    //     printf("Matrice M:\n");
    //     afficher_mat(M, n);
    //     initialiser_vecteur(u, n, 0.0);

    //     int kmax = 10;
    //     double eps = 1e-10;
    //     for (int k = 0; k < kmax; k++)
    //     {
    //         double residu_norm_m1 = vecteur_L2_norm(r, n);
    //         residu(A, b, u, r, n);
    //         // afficher_vect(u, n);
    //         double residu_norm = vecteur_L2_norm(r, n);
    //         // printf("iter: %d\tresidu norm: %4.2e\n", k, residu_norm);
    //         printf("iter: %d\tresidu norm: %4.2e\t rkp1/rk: %4.2e \n", k, residu_norm, residu_norm / residu_norm_m1);
    //         if (residu_norm < eps)
    //             break;
    //         resol_trig_inf(M, du, r, n);
    //         ukp1_iter(u, du, u, n);
    //     }
    // }

    // if (0)
    // { // Application
    //     //  Parametres
    //     int n = 100;
    //     double h = 1.0 / n;
    //     double A[MAX][MAX], b[MAX];
    //     // remplir_Ab_TD3_1(A, b);
    //     remplir_Ab_TP3(A, b, h, n);
    //     printf("Matrice A:\n");
    //     afficher_mat(A, n);
    //     printf("Vecteur b:\n");
    //     afficher_vect(b, n);

    //     // Resolution
    //     //  Jacobi, GS, SOR
    //     double M[MAX][MAX], N[MAX][MAX], r[MAX], u[MAX], du[MAX];
    //     // remplir_MN_Jacobi(A, M, N, n);
    //     // remplir_MN_Gauss_Seidel(A, M, N, n);
    //     double omega_sor = 1.26;
    //     remplir_MN_SOR(A, M, N, omega_sor, n);
    //     printf("Matrice M:\n");
    //     afficher_mat(M, n);
    //     initialiser_vecteur(u, n, 0.0);

    //     int kmax = 8000;
    //     // int kmax = 1;
    //     double eps = 1e-8;
    //     for (int k = 0; k < kmax; k++)
    //     {
    //         double residu_norm_m1 = vecteur_L2_norm(r, n);
    //         residu(A, b, u, r, n);
    //         // afficher_vect(u, n);
    //         double residu_norm = vecteur_L2_norm(r, n);
    //         // printf("iter: %d\tresidu norm: %4.2e\n", k, residu_norm);
    //         printf("iter: %d\tresidu norm: %4.2e\t rkp1/rk: %4.2e \n", k, residu_norm, residu_norm / residu_norm_m1);
    //         if (residu_norm < eps)
    //             break;
    //         resol_trig_inf(M, du, r, n);
    //         ukp1_iter(u, du, u, n);
    //     }
    //     char nom_de_sauvgarde[60] = "u_TP3_stationnaire.dat";
    //     sauvgarder_vect_et_maillage(u, n, h, nom_de_sauvgarde);
    // }

    // // integration temporelle
    // if (1)
    // {
    //     //  Parametres
    //     int n = 100;
    //     double h = 1.0 / n;
    //     double dt = 100 * h * h;
    //     double mu = dt / h / h;
    //     double Tmax = 4;
    //     int nt = Tmax / dt;
    //     nt = 30;
    //     char nom_de_sauvgarde[60] = "u30i.dat";

    //     double A[MAX][MAX], b[MAX], E[MAX][MAX], si[MAX], rhs[MAX];
    //     // remplir_Ab_TD3_1(A, b);
    //     remplir_Ab_TP3(A, b, h, n);
    //     printf("Matrice A:\n");
    //     afficher_mat(A, n);
    //     printf("Vecteur b:\n");
    //     afficher_vect(b, n);
    //     remplir_E_si_TP3(A, b, mu, E, si, n);
    //     printf("Matrice E:\n");
    //     afficher_mat(E, n);
    //     printf("Vecteur si:\n");
    //     afficher_vect(si, n);

    //     // Resolution
    //     //  Jacobi, GS, SOR
    //     double M[MAX][MAX], N[MAX][MAX], r[MAX], u[MAX], du[MAX];
    //     // remplir_MN_Jacobi(E, M, N, n);
    //     // remplir_MN_Gauss_Seidel(E, M, N, n);
    //     double omega_sor = 1.26;
    //     remplir_MN_SOR(E, M, N, omega_sor, n);
    //     printf("Matrice M:\n");
    //     afficher_mat(M, n);

    //     int kmax = 6000;
    //     // int kmax = 1;
    //     double eps = 1e-8;

    //     double Unplus1[MAX], Un[MAX];
    //     initialiser_vecteur(Unplus1, n, 0.0);
    //     initialiser_vecteur(Un, n, 0.0);

    //     // nt = 1;
    //     printf("Parametres:\t nx: %d\t h: %lf\t dt: %lf\t mu: %lf\t Tmax: %lf\t nt: %d\n",
    //            n, h, dt, mu, Tmax, nt);

    //     for (int it = 0; it < nt; it++)
    //     {
    //         xkp1_gradient(Un, si, 1.0, rhs, n);
    //         // afficher_vect(rhs, n);

    //         // resolution iterative
    //         double residu_norm;
    //         for (int k = 0; k < kmax; k++)
    //         {
    //             double residu_norm_m1 = vecteur_L2_norm(r, n);
    //             residu(E, rhs, Unplus1, r, n);
    //             // afficher_vect(u, n);
    //             residu_norm = vecteur_L2_norm(r, n);
    //             // printf("iter: %d\tresidu norm: %4.2e\n", k, residu_norm);
    //             // printf("iter: %d\tresidu norm: %4.2e\t rkp1/rk: %4.2e \n", k, residu_norm, residu_norm / residu_norm_m1);
    //             if (residu_norm < eps)
    //             {
    //                 printf("Conv en iter: %d\tresidu norm: %4.2e\t rkp1/rk: %4.2e \n", k, residu_norm, residu_norm / residu_norm_m1);
    //                 break;
    //             }

    //             resol_trig_inf(M, du, r, n);
    //             ukp1_iter(Unplus1, du, Unplus1, n);
    //         }
    //         if (!(residu_norm <eps))
    //         {
    //             printf("Code converge pas.\n");
    //             break;
    //         }

    //         remplir_vec_Un_avec_Unplus1(Un, Unplus1, n);
    //         if ((nt > 10 && (it + 1) % (nt / 10) == 0) || nt < 1000)
    //             printf("===== IE iter: %d temps: %lf dt: %lf\n", it + 1, (it + 1) * dt, dt);
    //         // afficher_vect(Unplus1, n);

    //         // sauvgarder_valeur_instantane(Un, nx, (it + 1) * dt, "valeurs_instantanes.dat");
    //     }

    //     sauvgarder_vect_et_maillage(Un, n, h, nom_de_sauvgarde);
    // }

    // TD3_conj_gradient();
    // TD3_gradient();

    return 0;
}