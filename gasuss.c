#include <stdio.h>

#define N  3   //the dimension of equation

void vec_diff(float a[N], float b[N]){
    /* Calcurate the difference of two vectors. Be caution that b[N] changes.*/
    for (int i = 0; i < N; i++){
        b[i] -= a[i];
    }
}

int main(){
    float m[N][N] = {{5,-1,-1}, {2,1,-3},{1,1,1}};    // The matrix
    float b[N] = {0,-5,6};

    printf("The coefficient matrix is : \n");
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            printf("%1.f ", m[i][j]);
            if (j == N-1){
                printf("\n");
            }
        }
    }

    printf("\nUse Gauss method to solve equations : \n");
    for (int i = 0; i < N; i++){
        for (int j = i+1; j < N; j++){
            float coef = m[j][i] / m[i][i];
            float del[N];

            for (int k = 0; k < N; k++){
                del[k] = m[i][k] * coef;
            }
            vec_diff(del, m[j]);
            b[j] -= b[i] * coef;
        }
    }

    for (int i = N -1; i >= 0; i--){
        float x = 1. / m[i][i];
        m[i][i] *= x;
        b[i] *= x;

        for (int j = 0; j < i; j++){
            b[j] -= b[i]*m[j][i];
            m[j][i] = 0;
        }
    }

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            printf("%1.f ", m[i][j]);
            if (j == N - 1){
                printf("\n");
            }
        }
    }

    for (int i = 0; i < N; i++){
        printf("%f ", b[i]);
    }

    return 0;
}
