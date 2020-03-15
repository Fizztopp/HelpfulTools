#include "matmul.h"

int idx(int i, int j, int M, int N) {
    return N * i + j;
}

void matmul(double *a, double *b, double *c, int m, int l, int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < l; k++) {
                sum += a[idx(i, k, m, l)] * b[idx(k, j, l, n)];
            }
            c[idx(i, j, m, n)] = sum;
        }
    }
}
