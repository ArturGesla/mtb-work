#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 10000
#define NMAX 100

float norme2_vect(float x[MAX], int n)
{
    float var = 0;
    int i;
    for (i = 0; i < n; i++)
        var += (x[i] * x[i]);

    return sqrt(var);
}

void ecritureContour(float Q[MAX], int nx, int ny, float h, float gs, float gw, float ge, float gn)
{
    FILE *fichier1;
    fichier1 = fopen("isoP_cnt.dat", "w");

    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            fprintf(fichier1, "%e\t%e\t%e\n", (ix + 1) * h, (iy + 1) * h, Q[ix + iy * nx]);
        }
        fprintf(fichier1, "\n");
    }

    fclose(fichier1);
}

void ecriture(float Q[MAX], float P[NMAX][NMAX], float x[NMAX], float y[NMAX],
              int nx, int ny, float h, float gs, float gw, float ge, float gn)
{

    int i, j, JJ;

    // remplissage des coordon�es du maillage x, y et de solution P � chaque point du maillage

    FILE *fichier1;
    fichier1 = fopen("isoP_sor.dat", "w");

    for (i = 0; i <= nx + 1; i++)
    {
        x[i] = i * h;
        P[i][0] = gs;
        P[i][ny + 1] = gn;
    }

    for (i = 0; i <= ny + 1; i++)
    {
        y[i] = i * h;
        P[0][i] = gw;
        P[nx + 1][i] = ge;
    }

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            JJ = i + j * nx;
            P[i + 1][j + 1] = Q[JJ];
        }
    }

    if (nx == 5)
    {
        for (i = 0; i <= nx + 1; i++)
            printf("%e ", x[i]);
        printf("\n");

        for (i = 0; i <= ny + 1; i++)
            printf("%e ", y[i]);
        printf("\n");

        for (i = 0; i <= nx + 1; i++)
        {
            for (j = 0; j <= ny + 1; j++)
            {
                printf("%e ", P[i][j]);
            }
            printf("\n");
        }
    }

    // ecriture du maillage et du r�sultat dans un fichier

    // fprintf(fichier1, "%d, %d\n", nx+2, ny+2);
    for (i = 0; i <= nx + 1; i++)
    {
        for (j = 0; j <= ny + 1; j++)
        {
            // fprintf(fichier1, "%e %e %e ", x[i], y[j], P[i][j]);
            fprintf(fichier1, "%e\t", P[i][j]);
        }
        fprintf(fichier1, "\n");
    }

    fclose(fichier1);
}

int main()
{

    // Discretisation
    int nx = 5, ny = 4, nn = nx * ny;
    //  int nx = 30, ny = 25, nn = nx * ny;
    //    int nx = 95, ny = 79, nn = nx * ny;
    float h = 1.5 / (nx + 1);

    //  int nx = 5, ny = 4, nn = nx * ny;
    // //  int nx = 11, ny = 9, nn = nx * ny;
    // //  int nx = 5, ny = 4, nn = nx * ny;
    //     float h = 1.5 / (nx + 1);
    //     float hx = 1.5 / (nx + 1);
    //     float hy = 1.25 / (ny + 1);
    //     if(hx!=hy) return -1;

    // Conditions aux limites
    float gs = 0., gw = 1., ge = 0.2, gn = 0.6, f = 0;
    // float gs = 0.3, gw = 0.5, ge = 1., gn = 0.0;

    // Diagonales non nulles de la matrice, vecteur second membre et r�sidu
    float As[MAX], Aw[MAX], Ap[MAX], Ae[MAX], An[MAX], b[MAX], r[MAX];

    // Pour la r�solution SOR
    int kmax = 1000, kech = 20;
    float eps = 1.e-6, omega = 1.9, nores = 1.;
    float Q[MAX], dQ[MAX];

    // Champ solution 2D
    float P[NMAX][NMAX];

    // Coordonn�es du maillage
    float x[NMAX], y[NMAX];

    // Pour le calcul de flux
    int mx, my;
    float T1 = 293., T2 = 333., lambda = 120;
    float gtv[NMAX], gth[NMAX];
    float IH, IV, QH, QV;

    int i, j, k, JJ;

    // Mise en equation
    // initiation des vecteurs non nuls de la matrice pentadiagonale As, Aw, Ap, Ae, An
    // et du second membre

    for (int i = 0; i < nn; i++)
    {
        Ap[i] = 4.0;
        if ((i + 1) % nx != 0)
            Ae[i] = -1;
        if ((i) % nx != 0)
            Aw[i] = -1;
        if (i < ((nx) * (ny - 1)))
            An[i] = -1;
        if (i > (nx - 1))
            As[i] = -1;

        if (i < nx)
            b[i] += gs;
        if (i > (nx) * (ny - 1) - 1)
            b[i] += gn;
        if (i % nx == 0)
            b[i] += gw;
        if ((i + 1) % nx == 0)
            b[i] += ge;
    }

    FILE *fichier = NULL;
    fichier = fopen("res_sor.dat", "w");

    printf("nx = %d, ny = %d\n", nx, ny);

    // impression de la matrice et du second membre
    if (nx == 5)
    {
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", As[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", Aw[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", Ap[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", Ae[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", An[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", b[i]);
        printf("\n");
    }

    // R�solution par SOR

    // Initialisation

    for (int i = 0; i < nn; i++)
    {
        Q[i] = 0;
        dQ[i] = 0;
    }

    // Boucle conditionnelle
    k = 0;
    nores = 1;
    while (k < kmax && nores > eps)
    {

        for (i = 0; i < nn; i++)
        {

            if (i < 1)
                r[i] = b[i] - Ap[i] * Q[i] - Ae[i] * Q[i + 1] - An[i] * Q[i + nx];
            else if (i < nx)
                r[i] = b[i] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1] - An[i] * Q[i + nx];
            else if (i > nn - nx)
                r[i] = b[i] - As[i] * Q[i - nx] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1];
            else if (i > nn - 1)
                r[i] = b[i] - As[i] * Q[i - nx] - Aw[i] * Q[i - 1] - Ap[i] * Q[i];
            else
                r[i] = b[i] - As[i] * Q[i - nx] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1] - An[i] * Q[i + nx];
        }

        nores = norme2_vect(r, nn);

        if ((k % kech) == 0)
        {
            printf("k = %d, r�sidu = %e\n", k, nores);
            fprintf(fichier, "k = %d, r�sidu = %e\n", k, nores);
            /*k d�signe ici le nombre d'it�rations, nores d�signe ici la norme du r�sidue*/
        }

        for (i = 0; i < nn; i++)
        {

            if (i < 1)
                dQ[i] = omega * r[i] / Ap[i];
            else if (i < nx)
                dQ[i] = omega * (r[i] - Aw[i] * dQ[i - 1]) / Ap[i];
            else
                dQ[i] = omega * (r[i] - As[i] * dQ[i - nx] - Aw[i] * dQ[i - 1]) / Ap[i];

            Q[i] = Q[i] + dQ[i];
        }

        k += 1;
    }

    printf("converegnce en k = %d iterations, r�sidu = %e\n", k, nores);
    fprintf(fichier, "converegnce en k = %d iterations, r�sidu = %e\n", k, nores);

    if (nx == 5)
    {
        for (i = 0; i < nn; i++)
            printf(" %e", Q[i]);
        printf("\n");
    }

    ecriture(Q, P, x, y, nx, ny, h, gs, gw, ge, gn);
    ecritureContour(Q, nx, ny, h, gs, gw, ge, gn);

    // Calcul de flux QH et QV

    printf("IH = %f, IV = %f\n", IH, IV);
    printf("QH = %f, QV = %f\n", QH, QV);

    fclose(fichier);

    return 0;
}
