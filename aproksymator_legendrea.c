#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi(double a, double b, int n, int i, double x){
    double p0, p1, p2;
    int index;
    p0 = 1;
    p1 = x;

    if(i == 0)
        return p0;
    if(i == 1)
        return p1;

    for(index = 2; index <= i; index++){
        p2 = (2 - 1/index) * (p1*x) - (1 - 1/index) * p0;
        p0 = p1;
        p1 = p2;
    }

    return p2;
}

/* Pierwsza pochodna fi */
double
dfi(double a, double b, int n, int i, double x){
    double p1;
    double dp0, dp1, dp2;
    int index;
    p1 = x;
    dp0 = 0;
    dp1 = 1;

    if(i == 0)
        return dp0;
    if(i == 1)
        return dp1;
    for(index = 2; index <= i; index++){
        dp2 = (2 - 1/index) * (p1 + dp1*x) - (1 - 1/index) * dp0;
        p1 = fi(a, b, n, index, x);
        dp0 = dp1;
        dp1 = dp2;
    }
    return dp2;
}

/* Druga pochodna fi */
double
d2fi(double a, double b, int n, int i, double x){
    double dp1;
    double d2p0, d2p1, d2p2;
    int index;
    dp1 = 1;
    d2p0 = 0;
    d2p1 = 0;

    if(i == 0)
        return d2p0;
    if(i == 1)
        return d2p1;
    for(index = 2; index <= i; index++){
        d2p2 = (2 - 1/index) * (2*dp1 + d2p1*x) - (1 - 1/index) * d2p0;
        dp1 = dfi(a, b, n, index, x);
        d2p0 = d2p1;
        d2p1 = d2p2;
    }
    return d2p2;
}

/* Trzecia pochodna fi */
double
d3fi(double a, double b, int n, int i, double x){
    double d2p1;
    double d3p0, d3p1, d3p2;
    int index;
    d2p1 = 0;
    d3p0 = 0;
    d3p1 = 0;

    if(i == 0)
        return d3p0;
    if(i == 0)
        return d3p1;
    for(index = 2; index <= i; index++){
        d3p2 = (2 - 1/index) * (3*d2p1 + d3p1*x) - (1 - 1/index) * d3p0;
        d2p1 = d2fi(a, b, n, index, x);
        d3p0 = d3p1;
        d3p1 = d3p2;
    }s
    return d3p2;
}

/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out)
{
    double		h = (b - a) / (n - 1);
    double		h3 = h * h * h;
    int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
    double		hx      [5];
    int		j;

    for (j = 0; j < 5; j++)
        hx[j] = a + h * hi[j];

    fprintf( out, "# nb=%d, i=%d: hi=[", n, i );
    for( j= 0; j < 5; j++ )
        fprintf( out, " %d", hi[j] );
    fprintf( out, "] hx=[" );
    for( j= 0; j < 5; j++ )
        fprintf( out, " %g", hx[j] );
    fprintf( out, "]\n" );
}

void
make_spl(points_t * pts, spline_t * spl)
{

    matrix_t       *eqs= NULL;
    double         *x = pts->x;
    double         *y = pts->y;
    double		a = x[0];
    double		b = x[pts->n - 1];
    int		i, j, k;
    int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

    if( nbEnv != NULL && atoi( nbEnv ) > 0 )
        nb = atoi( nbEnv );

    eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
    {
        FILE           *tst = fopen("debug_base_plot.txt", "w");
        double		dx = (b - a) / (TESTBASE - 1);
        for( j= 0; j < nb; j++ )
            xfi( a, b, nb, j, tst );
        for (i = 0; i < TESTBASE; i++) {
            fprintf(tst, "%g", a + i * dx);
            for (j = 0; j < nb; j++) {
                fprintf(tst, " %g", fi  (a, b, nb, j, a + i * dx));
                fprintf(tst, " %g", dfi (a, b, nb, j, a + i * dx));
                fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
                fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
            }
            fprintf(tst, "\n");
        }
        fclose(tst);
    }
#endif

    for (j = 0; j < nb; j++) {
        for (i = 0; i < nb; i++)
            for (k = 0; k < pts->n; k++)
                add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));

        for (k = 0; k < pts->n; k++)
            add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
    }

#ifdef DEBUG
    write_matrix(eqs, stdout);
#endif

    if (piv_ge_solver(eqs)) {
        spl->n = 0;
        return;
    }
#ifdef DEBUG
    write_matrix(eqs, stdout);
#endif

    if (alloc_spl(spl, nb) == 0) {
        for (i = 0; i < spl->n; i++) {
            double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
            xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
            spl->f[i] = 0;
            spl->f1[i] = 0;
            spl->f2[i] = 0;
            spl->f3[i] = 0;
            for (k = 0; k < nb; k++) {
                double		ck = get_entry_matrix(eqs, k, nb);
                spl->f[i]  += ck * fi  (a, b, nb, k, xx);
                spl->f1[i] += ck * dfi (a, b, nb, k, xx);
                spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
                spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
            }
        }
    }

#ifdef DEBUG
    {
        FILE           *tst = fopen("debug_spline_plot.txt", "w");
        double		dx = (b - a) / (TESTBASE - 1);
        for (i = 0; i < TESTBASE; i++) {
            double yi= 0;
            double dyi= 0;
            double d2yi= 0;
            double d3yi= 0;
            double xi= a + i * dx;
            for( k= 0; k < nb; k++ ) {
                            yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
                            dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
                            d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
                            d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
            }
            fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
        }
        fclose(tst);
    }
#endif

}
