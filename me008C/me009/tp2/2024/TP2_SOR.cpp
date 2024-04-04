#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 10000
#define NMAX 100
#pragma warning(disable:4996)



float norme2_vect(float x[MAX], int n) {
    float var = 0;
    int i;
    for (i = 0; i < n; i++) var += (x[i] * x[i]);

    return sqrt(var);
}

void ecriture (float Q[MAX], float P[NMAX][NMAX], float x[NMAX], float y[NMAX], int nx, int ny, float h, float gs, float gw, float ge, float gn ) {

    int i, j, JJ;

    // remplissage des coordonées du maillage x, y et de solution P à chaque point du maillage

    FILE* fichier1;
    fichier1 = fopen("isoP_sor.dat", "w");

    for (i = 0; i <= nx+1; i++) {
        x[i] = i * h;
        P[i][0] = gs;
        P[i][ny + 1] = gn;
    }
   
    for (i = 0; i <= ny+1; i++) {
        y[i] = i * h;
        P[0][i] = gw;
        P[nx + 1][i] = ge;
    }

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            JJ = i + j * nx;
            P[i + 1][j + 1] = Q[JJ];
        }
    }

    if (nx == 5) {
        for (i = 0; i <= nx+1; i++) printf("%e ", x[i]);
        printf("\n");

        for (i = 0; i <= ny+1; i++) printf("%e ", y[i]);
        printf("\n");

        for (i = 0; i <= nx + 1; i++) {
            for (j = 0; j <= ny + 1; j++) {
                printf("%e ", P[i][j]);
            }
            printf("\n");
        }
    }

    // ecriture du maillage et du résultat dans un fichier

    fprintf(fichier1, "%d, %d\n", nx+2, ny+2);
    for (i = 0; i <= nx + 1; i++) {
        for (j = 0; j <= ny + 1; j++) {
            fprintf(fichier1, "%e %e %e ", x[i], y[j], P[i][j]);
        }
        fprintf(fichier1, "\n");
    }

    fclose(fichier1);
}

int main() {

    // discretisation 
//	int nx = 5, ny = 4, nn = nx * ny;
//    int nx = 30, ny = 25, nn = nx * ny;
    int nx = 95, ny = 79, nn = nx * ny;
    float h = 1.5 / (nx + 1);

    // condition aux limites
    float gs = 0., gw = 1., ge = 0.2, gn = 0.6, f = 0;
    //float gs = 0.3, gw = 0.5, ge = 1., gn = 0.0;

    // matrice et second membre
    float As[MAX], Aw[MAX], Ap[MAX], Ae[MAX], An[MAX], b[MAX], r[MAX];

    // pour la résolution SOR
    int kmax = 5000, kech = 20; 
    float eps = 1.e-6, omega = 1.9, nores = 1.;
    float Q[MAX], dQ[MAX];

    // solution
    float P[NMAX][NMAX];
    
    // coordonées du maillage
    float x[NMAX], y[NMAX];
 
    // pour le calcul de flux
    int mx, my;
    float T1 = 293., T2 = 333., lambda = 120;
    float gtv[NMAX], gth[NMAX];
    float IH, IV, QH, QV;

    int i, j, k, JJ;


    // Mise en equation
    // initiation de la matrice pentadiagonale As, Aw, Ap, Ae, An et du second membre

    for (i = 0; i < nn; i++) {
        Ap[i] = 4;
        As[i] = Aw[i] = Ae[i] = An[i] = -1;
        b[i] = f*h*h;
    }

    // prise en compte des condistions aux limites

    for (i = 0; i < nx; i++) {
        As[i] = 0;
        b[i] += gs;
        JJ = i + (ny - 1) * nx;
        An[JJ] = 0;
        b[JJ] += gn;
    }
    for (j = 0; j < ny; j++) {
        JJ = j * nx;
        if (JJ >= 0) {
            Aw[JJ] = 0;
            b[JJ] += gw;
        }
        JJ = nx + j * nx - 1;
        if (JJ < nn) {
            Ae[JJ] = 0;
            b[JJ] += ge;
        }
    }

    FILE* fichier = NULL;
    fichier = fopen("res_sor.dat", "w");

    printf("nx = %d, ny = %d\n", nx, ny);

    // impression de la matrice et du second membre
    if (nx == 5) {
        for (int i = 0; i < nn; i++) printf("%.1f, ", As[i]);
        printf("\n");
        for (int i = 0; i < nn; i++) printf("%.1f, ", Aw[i]);
        printf("\n");
        for (int i = 0; i < nn; i++) printf("%.1f, ", Ap[i]);
        printf("\n");
        for (int i = 0; i < nn; i++) printf("%.1f, ", Ae[i]);
        printf("\n");
        for (int i = 0; i < nn; i++) printf("%.1f, ", An[i]);
        printf("\n");
        for (int i = 0; i < nn; i++) printf("%.1f, ", b[i]);
        printf("\n");
    }
    
    // Résolution par SOR

    // initialisation
    
    for (i = 0; i < nn; i++) {
        Q[i] = 0;
        dQ[i] = 0;
    }

    nores = 1.;
    k = 0;

    // boucle conditionnelle

    while (nores >= eps && k < kmax) {

        for (i = 0; i < nn; i++) {
            
            if (i < 1)
                r[i] = b[i] - Ap[i] * Q[i] - Ae[i] * Q[i+1] - An[i] * Q[i+nx];
            else if (i < nx)
                r[i] = b[i]  - Aw[i] * Q[i-1] - Ap[i] * Q[i] - Ae[i] * Q[i+1] - An[i] * Q[i+nx];
            else if (i > nn-nx)
                r[i] = b[i] - As[i] * Q[i-nx] - Aw[i] * Q[i-1] - Ap[i] * Q[i] - Ae[i] * Q[i+1];
            else if (i > nn-1)
                r[i] = b[i] - As[i] * Q[i-nx] - Aw[i] * Q[i-1] - Ap[i] * Q[i];
            else
                r[i] = b[i] - As[i] * Q[i-nx] - Aw[i] * Q[i-1] - Ap[i] * Q[i] - Ae[i] * Q[i+1] - An[i] * Q[i+nx];
            
            /*
            float AsQ=0, AwQ=0, AeQ=0, AnQ=0;
            if ((i - nx) >= 0) AsQ = As[i] * Q[i-nx];
            if ((i - 1) >= 0) AwQ = Aw[i] * Q[i-1];
            if ((i + 1) < nn) AeQ = Ae[i] * Q[i+1];
            if ((i + nx) < nn) AnQ = An[i] * Q[i+nx];
            r[i] = b[i] - AsQ - AwQ - Ap[i] * Q[i] - AeQ - AnQ;
            */
        }

        nores = norme2_vect(r, nn);
    
        if ((k % kech) == 0) {
             printf("k = %d, résidu = %e\n", k, nores);
             fprintf(fichier, "k = %d, résidu = %e\n", k, nores);
             /*k désigne ici le nombre d'itération, nores désigne ici la norme du résidue*/
        }
        
        for (i = 0; i < nn; i++) {
           
            if (i < 1)
                dQ[i] = omega * r[i] / Ap[i];
            else if (i < nx)
                dQ[i] = omega * (r[i] - Aw[i] * dQ[i-1]) / Ap[i];
            else
                dQ[i] = omega * (r[i] - As[i] * dQ[i-nx] - Aw[i] * dQ[i-1]) / Ap[i];
           
            /*
            float AsdQ = 0, AwdQ = 0;
            if ((i - nx) >= 0) AsdQ = As[i] * dQ[i-nx];
            if ((i - 1) >= 0) AwdQ = Aw[i] * dQ[i-1];
            dQ[i] = omega * (r[i] - AsdQ - AwdQ) / Ap[i];
            */

            Q[i] = Q[i] + dQ[i];
        }
 
        k += 1;
     }
 
    printf("converegnce en k = %d iterations, résidu = %e\n", k, nores);
    fprintf(fichier, "converegnce en k = %d iterations, résidu = %e\n", k, nores);
    
    if (nx == 5) {
        for (i = 0; i < nn; i++) printf(" %e", Q[i]);
        printf("\n");
    }

    ecriture(Q, P, x, y, nx, ny, h, gs, gw, ge, gn);

    // Calcul de flux QH et QV

    mx = (nx + 1) / 2;
    my = (ny + 1) / 2;
    printf(" mx = %d, my = %d, h = %f\n", mx, my, h);

    for (i = 0; i <= ny+1; i++) gth[i] = (P[mx+1][i] - P[mx-1][i]) / (2*h);
    for (i = 0; i <= nx+1; i++) gtv[i] = (P[i][my+1] - P[i][my-1]) / (2*h);

    IH = 0;
    for (i = 0; i <= ny; i++) IH = IH - 0.5 * h * (gth[i] + gth[i+1]);
    IV = 0;
    for (i = 0; i <= nx; i++) IV = IV - 0.5 * h * (gtv[i] + gtv[i+1]);

    QH = lambda * (T2 - T1) * IH;
    QV = lambda * (T2 - T1) * IV;

    printf("IH = %f, IV = %f\n", IH, IV);
    printf("QH = %f, QV = %f\n", QH, QV);

    fclose(fichier);

	return 0;

}

