#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 4
#define EPS 1.e-8

double L_oo(double x[])
{
  double m = fabs(x[0]);
  for (int i = 1; i < N; i++)
    if (fabs(x[i]) > m)
      m = fabs(x[i]);
  return m;
}

double prod_scal(double x[], double y[])
{
  double x_scal_y = 0;
  for (int i = 0; i < N; i++)
    x_scal_y += x[i] * y[i];
  return x_scal_y;
}

void prod_Ax(double A[N][N], double x[], double y[])
{
  for (int i = 0; i < N; i++)
  {
    y[i] = 0;
    for (int j = 0; j < N; j++)
    {
      y[i] += A[i][j] * x[j];
    }
  }
}
void prod_cx(double c, double x[], double y[])
{
  for (int i = 0; i < N; i++)
  {
    y[i] = c * x[i];
  }
}
double L_2(double x[])
{
  return sqrt(prod_scal(x, x));
}

void transpose(double A[][N], double At[][N])
{
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      At[i][j] = A[j][i];
  return;
}

void afficher_vect(double x[])
{
  for (int i = 0; i < N; i++)
    printf("%1lf\t", x[i]);
  printf("\n");
  return;
}

void afficher_mat(double A[][N])
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
      printf("%.4lf\t", A[i][j]);
    printf("\n");
  }
  return;
}

// DESCENTE
void resol_trig_inf(double L[][N], double x[], double b[])
{
  for (int i = 0; i < N; i++)
  {
    double l = 0;
    for (int k = 0; k < i; k++)
      l += L[i][k] * x[k];
    x[i] = (b[i] - l) / L[i][i];
  }
  return;
}

// MONTÉE
void resol_trig_sup(double U[][N], double x[], double b[])
{
  for (int i = N - 1; i >= 0; i--)
  {
    double l = 0;
    for (int k = i + 1; k < N; k++)
      l += U[i][k] * x[k];
    x[i] = (b[i] - l) / U[i][i];
  }
  return;
}

void fact_LU(double A[][N], double L[][N], double U[][N])
{
  for (int j = 0; j < N; j++)
  {
    for (int i = 0; i <= j; i++)
    {
      double l = 0;
      for (int k = 0; k < i; k++)
        l += L[i][k] * U[k][j];
      U[i][j] = A[i][j] - l;
    }
    L[j][j] = 1;
    for (int i = j + 1; i < N; i++)
    {
      double l = 0;
      for (int k = 0; k < j; k++)
        l += L[i][k] * U[k][j];
      L[i][j] = (A[i][j] - l) / U[j][j];
    }
  }
  return;
}

double puiss_it(double A[][N], double x[])
{
  double u[N], lambda, lambda_old;
  double y[N];
  lambda = 1;
  // afficher_vect(y);
  // A FAIRE!

  int nit = 40;
  // for (int i = 0; i < nit; i++)
  int i = 0;
  while (fabs(lambda_old - lambda) > 1e-6 && i < nit)
  {
    prod_Ax(A, x, y);
    lambda_old = lambda;
    lambda = L_2(y);
    // printf("it: %d current ev: %4.8f dev: %4.2e \n", i, lambda, fabs(lambda_old - lambda));
    // normalise
    prod_cx(1 / lambda, y, x);
    // printf("current evc: \n");
    // afficher_vect(x);
    i++;
  }
  printf("it: %d current ev: %4.8f dev: %4.2e \n", i, lambda, fabs(lambda_old - lambda));
  printf("Vect:\n");
  afficher_vect(x);

  if (fabs(lambda_old - lambda) > 1e-3)
  {
    printf("========== Error: no convergence. ==========");
  }
  return lambda;
}

double puiss_inv(double A[][N], double x[])
{
  double L[N][N], U[N][N];
  double v[N], y[N];
  double lambda, lambda_old;
  lambda = 1;

  int nit = 40;
  // for (int i = 0; i < nit; i++)
  int i = 0;
  fact_LU(A, L, U);

  while (fabs(lambda_old - lambda) > 1e-15 && i < nit)
  {
    resol_trig_inf(L, v, x);
    resol_trig_sup(U, y, v);
    lambda_old = lambda;
    lambda = L_2(y);
    // printf("it: %d current ev: %4.8f dev: %4.2e \n", i, lambda, fabs(lambda_old - lambda));
    prod_cx(1 / lambda, y, x);
    i++;
  }
  printf("it: %d current ev: %4.8f dev: %4.2e \n", i, 1 / lambda, fabs(lambda_old - lambda));
  printf("Vect:\n");
  afficher_vect(x);

  if (abs(lambda_old - lambda) > 1e-16)
  {
    printf("========== Error: no convergence. ==========");
  }
  return lambda;

  // A FAIRE!
}

void deflation_it(double A[][N], double x[], double lambda)
{
  double At[N][N], xt[N];
  prod_cx(1, x, xt);
  transpose(A, At);
  double l1 = puiss_it(At, xt);
  if (fabs(l1 - lambda) > 1e-14)
    printf("========== Error: diff ev of transpose. Diff: %4.2e==========\n", fabs(l1 - lambda));

  // construct deflated operator
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      A[i][j] = A[i][j] - lambda * x[i] * xt[j] / prod_scal(x, xt);
    }
  }
  // afficher_mat(A);
  // afficher_mat(At);

  // A FAIRE!
}

//===================SOR======================
#include <math.h>
#define MAX 10000
#define NMAX 100

double norme2_vect(double x[MAX], int n)
{
  double var = 0;
  int i;
  for (i = 0; i < n; i++)
    var += (x[i] * x[i]);

  return sqrt(var);
}

void ecritureContour(double Q[MAX], int nx, int ny, double h, double gs, double gw, double ge, double gn)
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

void ecriture(double Q[MAX], double P[NMAX][NMAX], double x[NMAX], double y[NMAX],
              int nx, int ny, double h, double gs, double gw, double ge, double gn)
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

int main(void)
{
  // double x1[N], x2[N], x3[N], x4[N];
  // double A[N][N];
  // double L1, L2, L3, L4;

  // Initialisation de A ici matrice de Vandermonde - A FAIRE

  // // diffusion
  // {
  //   int n = 3;
  //   for (int i = 0; i < n; i++)
  //   {
  //     A[i][i] = 2.0;
  //     if (i > 0)
  //       A[i][i - 1] = -1;
  //     if (i < n - 1)
  //       A[i][i + 1] = -1;
  //   }

  //   afficher_mat(A);
  //   double x1[N] = {0, 0, 1};
  //   prod_cx(1, x1, x2);
  //   prod_cx(1, x1, x3);

  //   // double x0[N] = {1.0 / 2, sqrt(2) / 2, 1.0 / 2};
  //   // afficher_vect(x3);

  //   // afficher_vect(x1);
  //   L1 = puiss_it(A, x1);
  //   deflation_it(A, x1, L1);
  //   L2 = puiss_it(A, x2);
  //   deflation_it(A, x2, L2);
  //   L3 = puiss_it(A, x3);
  // }

  // // Vandermonde
  // {
  //   int n = 4;
  //   for (int i = 0; i < n; i++)
  //   {
  //     for (int j = 0; j < n; j++)
  //     {
  //       A[i][j] = pow((i + 1), j + 1);
  //     }
  //   }
  //   afficher_mat(A);
  //   double x1[N] = {0, 0, 1};
  //   prod_cx(1, x1, x2);
  //   prod_cx(1, x1, x3);
  //   prod_cx(1, x1, x4);

  //   // double x0[N] = {1.0 / 2, sqrt(2) / 2, 1.0 / 2};
  //   // afficher_vect(x3);

  //   // afficher_vect(x1);

  //   // L1 = puiss_it(A, x1);
  //   // deflation_it(A, x1, L1);
  //   // L2 = puiss_it(A, x2);
  //   // deflation_it(A, x2, L2);
  //   // L3 = puiss_it(A, x3);
  //   // deflation_it(A, x3, L3);
  //   // L4 = puiss_it(A, x4);

  //   L1 = puiss_inv(A, x1);
  // }

  // Recherche des valeurs propres par la puissance itérée

  int a = 0;
  int bb = 7;
  int c = 1;
  // int d = 7;

  // ex1
  {
    double x1[N], x2[N], x3[N], x4[N];
    double A[N][N];
    double L1, L2, L3, L4;

    double alpha;
    if (a >= 0 && a <= 3)
      alpha = 0.25;
    if (a >= 4 && a <= 7)
      alpha = 0.5;
    if (a >= 8 && a <= 9)
      alpha = 1.0;

    int n = 4;
    for (int i = 0; i < n; i++)
    {
      A[i][i] = 2.0;
      if (i > 0)
        A[i][i - 1] = -1;
      if (i < n - 1)
        A[i][i + 1] = -1;
    }
    A[0][2] = alpha;

    afficher_mat(A);

    double x0[N] = {0, 0, 1};

    L1 = puiss_it(A, x0);
    L1 = puiss_inv(A, x0);
  }

  // ex2
  {
    // Discretisation
    // int nx = 5, ny = 4, nn = nx * ny;
    int nx = 11, ny = 9, nn = nx * ny;
    // int nx = 23, ny = 19, nn = nx * ny;
    // int nx = 47, ny = 39, nn = nx * ny;
    //  int nx = 30, ny = 25, nn = nx * ny;
    //    int nx = 95, ny = 79, nn = nx * ny;
    double h = 1.5 / (nx + 1);

    //  int nx = 5, ny = 4, nn = nx * ny;
    // //  int nx = 11, ny = 9, nn = nx * ny;
    // //  int nx = 5, ny = 4, nn = nx * ny;
    //     double h = 1.5 / (nx + 1);
    //     double hx = 1.5 / (nx + 1);
    //     double hy = 1.25 / (ny + 1);
    //     if(hx!=hy) return -1;

    // Conditions aux limites
    // double gs = 0., gw = 1., ge = 0.2, gn = 0.6, f = 0;
    double gs = 0.5, gw = 0., ge = 1.0, gn = 0.2, f = 0;
    // double gs = 0.3, gw = 0.5, ge = 1., gn = 0.0;

    // Diagonales non nulles de la matrice, vecteur second membre et r�sidu
    double As[MAX], Aw[MAX], Ap[MAX], Ae[MAX], An[MAX], b[MAX], r[MAX];

    // Pour la r�solution SOR
    int kmax = 10000, kech = 20;
    double eps = 1.e-8, omega = 1.09, nores = 1.;
    double Q[MAX], dQ[MAX];

    if(a<=4)
    omega=1.5+a/10.0;
    if(a>=5)
    omega=1.+a/10.0;
    kmax=10.0 *bb;
    printf("Ex2 ==================================================================================== omega = %4.1f\n",omega);

    // Champ solution 2D
    double P[NMAX][NMAX];

    // Coordonn�es du maillage
    double x[NMAX], y[NMAX];

    // Pour le calcul de flux
    int mx, my;
    double T1 = 293., T2 = 333., lambda = 120;
    double gtv[NMAX], gth[NMAX];
    double IH, IV, QH, QV;

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

    printf("nx = %d, ny = %d\t h=%e\n", nx, ny, h);

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
        // printf("k = %d, residu = %e\n", k, nores);
        fprintf(fichier, "k = %d, residu = %e\n", k, nores);
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

    printf("converegnce en k = %d iterations, residu = %e\n", k, nores);
    fprintf(fichier, "converegnce en k = %d iterations, residu = %e\n", k, nores);

    if (nx == 5)
    {
      for (i = 0; i < nn; i++)
        printf(" %e", Q[i]);
      printf("\n");
    }

    ecriture(Q, P, x, y, nx, ny, h, gs, gw, ge, gn);
    ecritureContour(Q, nx, ny, h, gs, gw, ge, gn);

    // Calcul de flux QH et QV

    {
      QH = 0;
      int ix = (nx - 1) / 2;
      int iy = -1;
      QH += 1 / 4.0 * (Q[ix + iy * nx + 1 + nx] - Q[ix + iy * nx - 1 + nx]);
      for (iy = 0; iy < ny - 1; iy++)
      {
        QH += 1 / 4.0 * (Q[ix + iy * nx + 1] - Q[ix + iy * nx - 1] + Q[ix + iy * nx + 1 + nx] - Q[ix + iy * nx - 1 + nx]);
      }
      iy = ny - 1;
      QH += 1 / 4.0 * (Q[ix + iy * nx + 1] - Q[ix + iy * nx - 1]);
      QH = -QH;
      IH = QH;
      QH = QH * (T2 - T1) * lambda;
    }

    printf("IH = %f, IV = %f\n", IH, IV);
    printf("QH = %f, QV = %f\n", QH, QV);

    fclose(fichier);
  }

   // ex4
  {
    // Discretisation
    // int nx = 5, ny = 4, nn = nx * ny;
    int nx = 11, ny = 9, nn = nx * ny;

    // if(c<=5){
    //   nx=29;
    //   ny=24;
    // }
    // if(c>=6){
    //   nx=23;
    //   ny=19;
    // }
    //  nn = nx * ny;
    // int nx = 23, ny = 19, nn = nx * ny;
    // int nx = 47, ny = 39, nn = nx * ny;
    //  int nx = 30, ny = 25, nn = nx * ny;
    //    int nx = 95, ny = 79, nn = nx * ny;
    double h = 1.5 / (nx + 1);

    //  int nx = 5, ny = 4, nn = nx * ny;
    // //  int nx = 11, ny = 9, nn = nx * ny;
    // //  int nx = 5, ny = 4, nn = nx * ny;
    //     double h = 1.5 / (nx + 1);
    //     double hx = 1.5 / (nx + 1);
    //     double hy = 1.25 / (ny + 1);
    //     if(hx!=hy) return -1;

    // Conditions aux limites
    // double gs = 0., gw = 1., ge = 0.2, gn = 0.6, f = 0;
    double gs = 0.5, gw = 0., ge = 1.0, gn = 0.2, f = 0;
    // double gs = 0.3, gw = 0.5, ge = 1., gn = 0.0;

    // Diagonales non nulles de la matrice, vecteur second membre et r�sidu
    double As[MAX], Aw[MAX], Ap[MAX], Ae[MAX], An[MAX], b[MAX], r[MAX];

    // Pour la r�solution SOR
    int kmax = 10000, kech = 20;
    double eps = 1.e-8, omega = 1.09, nores = 1.;
    double Q[MAX], dQ[MAX];

    if(a<=4)
    omega=1.5+a/10.0;
    if(a>=5)
    omega=1.+a/10.0;
    // kmax=10.0 *bb;
    printf("Ex4 ==================================================================================== omega = %4.1f\n",omega);

    // Champ solution 2D
    double P[NMAX][NMAX];

    // Coordonn�es du maillage
    double x[NMAX], y[NMAX];

    // Pour le calcul de flux
    int mx, my;
    double T1 = 293., T2 = 333., lambda = 120;
    double gtv[NMAX], gth[NMAX];
    double IH, IV, QH, QV;

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

    printf("nx = %d, ny = %d\t h=%e\n", nx, ny, h);

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
        // printf("k = %d, residu = %e\n", k, nores);
        fprintf(fichier, "k = %d, residu = %e\n", k, nores);
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

    printf("converegnce en k = %d iterations, residu = %e\n", k, nores);
    fprintf(fichier, "converegnce en k = %d iterations, residu = %e\n", k, nores);

    if (nx == 5)
    {
      for (i = 0; i < nn; i++)
        printf(" %e", Q[i]);
      printf("\n");
    }

    ecriture(Q, P, x, y, nx, ny, h, gs, gw, ge, gn);
    ecritureContour(Q, nx, ny, h, gs, gw, ge, gn);

    // Calcul de flux QH et QV

    {
      QH = 0;
      int ix = (nx - 1) / 2;
      int iy = -1;
      QH += 1 / 4.0 * (Q[ix + iy * nx + 1 + nx] - Q[ix + iy * nx - 1 + nx]);
      for (iy = 0; iy < ny - 1; iy++)
      {
        QH += 1 / 4.0 * (Q[ix + iy * nx + 1] - Q[ix + iy * nx - 1] + Q[ix + iy * nx + 1 + nx] - Q[ix + iy * nx - 1 + nx]);
      }
      iy = ny - 1;
      QH += 1 / 4.0 * (Q[ix + iy * nx + 1] - Q[ix + iy * nx - 1]);
      QH = -QH;
      IH = QH;
      QH = QH * (T2 - T1) * lambda;
    }

    printf("IH = %f, IV = %f\n", IH, IV);
    printf("QH = %f, QV = %f\n", QH, QV);

    fclose(fichier);
  }

  return 0;
}
