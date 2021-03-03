#include <iostream>
#include "matrix.h"
#include <math.h>

using namespace std;

void makeMatrix(MATRIX(matrix), float h) 
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (j == i) {
                matrix[i][j] = exp(i * j * h);
            }
            else {
                matrix[i][j] = exp(i * j * h) - i * j * h;
            }
        }
    }
}

void printMatrix(MATRIX(d), int n)
{
    printf("\nMatrix:\n");
    for (int i = 0; i < n; i++)
    {
        if (i == 0) printf("%c", 218);
        else
            if (i == (n - 1)) printf("%c", 192);
            else
                printf("%c", 179);
        for (int j = 0; j < n; j++)
        {
            printf("%11.6f", d[i][j]);
        }
        if (i == 0) printf("   %c", 191);
        else
            if (i == (n - 1)) printf("   %c", 217);
            else
                printf("   %c", 179);
        printf("\n");
    }
}

void makeMatrixE(MATRIX(matrix))
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (j == i) {
                matrix[i][j] = 1;
            }
            else {
                matrix[i][j] = 0;
            }
        }
    }
}

void inverseMatrix(MATRIX(input), MATRIX(output), float* cond) 
{
    float work[4];
    int ipvt[4];
    MATRIX(result);

    decomp(4, input, cond, ipvt, work);
    
    makeMatrixE(result);

    for (int i = 0; i < 4; i++) {
        solve(4, input, &result[i][0], ipvt);
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            output[i][j] = result[j][i];
        }
    }
}

void solveSystemOfLinearEquations(MATRIX(inverse), float* vect, N_TYPE* result)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i] += inverse[i][j] * vect[j];
            //cout << result[i];
            //cout << " result[";
            //cout << i;
            //cout << "]" << endl;
        }
    }
}

void multiMatrix(MATRIX(m1), MATRIX(m2), MATRIX(m3))
{
    for (int k = 0; k < 4; k++)
        for (int j = 0; j < 4; j++)
        {
            m3[k][j] = 0;
            for (int i = 0; i < 4; i++)
                m3[k][j] += m1[k][i] * m2[i][j];
        }
}

void transponeMatrix(MATRIX(input), MATRIX(output))
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            output[i][j] = input[j][i];
        }
    }
}

float norma(float* x)
{
    float sum = 0.0;
    for (int i = 0; i < 4; i++) {
        sum += x[i];
    }
    float result = sqrt(sum);
    return result;
}

void main()
{
    MATRIX(matrix);
    float cond[3];
    MATRIX(inverse);

    MATRIX(copy);
    MATRIX(transpone);
    MATRIX(multiResult);
    float cond2[3];
    MATRIX(inverseMultiResult);
    MATRIX(resultMatrix);

    float h[] = { 0.1, 0.01, 0.001 };
    float vector[4];
    vector[0] = 1.0;
    vector[1] = 2.0;
    vector[2] = 3.0;
    vector[3] = 4.0;

    for (int i = 0; i < 3; i++) {
        //Первый способ
        makeMatrix(matrix, h[i]);
        printMatrix(matrix, 4);

        inverseMatrix(matrix, inverse, &cond[i]);
        //printMatrix(inverse, 4);

        cout << endl;
        cout << "cond = " << cond[i] << endl;

        float x1[] = { 0.0, 0.0, 0.0, 0.0 };

        solveSystemOfLinearEquations(inverse, vector, x1);

        cout << endl;
        cout << "X1:";
        for (int i = 0; i < 4; i++) {
            printf("%12.5f", x1[i]);
        }
        cout << endl;

        //Второй способ
        makeMatrix(copy, h[i]);
        //printMatrix(copy, 4);

        transponeMatrix(copy, transpone);
        //printMatrix(transpone, 4);

        multiMatrix(copy, transpone, multiResult);
        //printMatrix(multiResult, 4);

        inverseMatrix(multiResult, inverseMultiResult, &cond2[i]);
        //printMatrix(inverseMultiResult, 4);

        multiMatrix(inverseMultiResult, transpone, resultMatrix);
        //printMatrix(resultMatrix, 4);

        float x2[] = { 0.0, 0.0, 0.0, 0.0 };

        solveSystemOfLinearEquations(resultMatrix, vector, x2);

        cout << endl;
        cout << "X2:";
        for (int i = 0; i < 4; i++) {
            printf("%12.5f", x2[i]);
        }
        cout << endl;

        //Оценка, вычисляемая в ходе эксперимента
        float x3[4];
        for (int i = 0; i < 4; i++) {
            x3[i] = x1[i] - x2[i];
        }

        float norma1 = norma(x3);
        float norma2 = norma(x1);

        cout << endl;
        cout << "sigma = " << norma1 / norma2 << endl;
    }

}