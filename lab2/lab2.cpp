#include <iostream>
#include "matrix.h"
#include <math.h>
//#include <iomanip>

using namespace std;

void makeMatrix(MATRIX(matrix), float h) 
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (j == i) {
                matrix[i][j] = exp((i + 1.0) * (j + 1.0) * h);
                //cout << endl;
                //cout << setprecision(18) << matrix[i][j] << endl;
            }
            else {
                matrix[i][j] = exp((i + 1.0) * (j + 1.0) * h) - (i + 1.0) * (j + 1.0) * h;
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
            printf("%21.11e", d[i][j]);
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

void inverseMatrix(MATRIX(input), MATRIX(output), double* cond) 
{
    //cout << "--------------------------------" << endl;
    //printMatrix(input, 4);
    double work[4];
    int ipvt[4];
    MATRIX(result);

    decomp(4, input, cond, ipvt, work);

    //printMatrix(input, 4);
    
    makeMatrixE(result);

    //printMatrix(result, 4);

    for (int i = 0; i < 4; i++) {
        solve(4, input, &result[i][0], ipvt);
        //cout << result[i][0] << endl;
        //printMatrix(result, 4);
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            output[i][j] = result[j][i];
        }
    }

    //printMatrix(output, 4);
    //cout << "--------------------------------" << endl;
}

void solveSystemOfLinearEquations(MATRIX(inverse), float* vect, N_TYPE* result)
{
    //printMatrix(inverse, 4);
    //for (int i = 0; i < 4; i++) {
    //    cout << vect[i] << " ";
    //}
    //cout << endl;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i] += inverse[i][j] * vect[j];
        }
        /*cout << " result[";
        cout << i;
        cout << "] = ";
        cout << inverse[i][0];
        cout << " * ";
        cout << vect[0];
        cout << " + ";
        cout << inverse[i][1];
        cout << " * ";
        cout << vect[1];
        cout << " + ";
        cout << inverse[i][2];
        cout << " * ";
        cout << vect[2];
        cout << " + ";
        cout << inverse[i][3];
        cout << " * ";
        cout << vect[3];
        cout << " = ";
        cout << result[i];
        cout << endl;*/
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

double norma(double* x)
{
    double sum = 0.0;
    for (int i = 0; i < 4; i++) {
        //cout << "sum = " << sum << endl;
        sum += x[i] * x[i];
    }
    //cout << "sum = " << sum << endl;
    double result = sqrt(sum);
    //cout << "result = " << result << endl;
    return result;
}

void main()
{
    MATRIX(matrix);
    double cond;
    MATRIX(inverse);

    MATRIX(copy);
    MATRIX(transpone);
    MATRIX(multiResult);
    MATRIX(inverseMultiResult);
    MATRIX(resultMatrix);

    float h[] = { 0.1, 0.01, 0.001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10 };
    float vector[4];
    vector[0] = 1.0;
    vector[1] = 2.0;
    vector[2] = 3.0;
    vector[3] = 4.0;

    for (int i = 0; i < 10; i++) {
        //Первый способ
        makeMatrix(matrix, h[i]);
        printMatrix(matrix, 4);

        inverseMatrix(matrix, inverse, &cond);
        //printMatrix(inverse, 4);

        /*MATRIX(U);
        MATRIX(copy2);
        makeMatrix(copy2, h[i]);
        multiMatrix(copy2, inverse, U);
        printMatrix(U, 4);*/

        cout << endl;
        cout << "cond = " << cond << endl;

        double x1[] = { 0.0, 0.0, 0.0, 0.0 };

        solveSystemOfLinearEquations(inverse, vector, x1);

        cout << endl;
        cout << "X1:";
        for (int i = 0; i < 4; i++) {
            printf("%21.5f", x1[i]);
        }
        cout << endl;

        //Второй способ
        makeMatrix(copy, h[i]);
        //printMatrix(copy, 4);

        transponeMatrix(copy, transpone);
        //printMatrix(transpone, 4);

        multiMatrix(copy, transpone, multiResult);
        //printMatrix(multiResult, 4);

        inverseMatrix(multiResult, inverseMultiResult, &cond);
        //printMatrix(inverseMultiResult, 4);

        multiMatrix(inverseMultiResult, transpone, resultMatrix);
        //printMatrix(resultMatrix, 4);

        double x2[] = { 0.0, 0.0, 0.0, 0.0 };

        solveSystemOfLinearEquations(resultMatrix, vector, x2);

        cout << endl;
        cout << "X2:";
        for (int i = 0; i < 4; i++) {
            printf("%21.5f", x2[i]);
        }
        cout << endl;

        //Оценка, вычисляемая в ходе эксперимента
        double x3[4];
        //cout << "X3:";
        for (int i = 0; i < 4; i++) {
            x3[i] = x1[i] - x2[i];
            //printf("%18.8e", x3[i]);
        }
        //cout << endl;

        double norma1 = norma(x3);
        //cout << norma1 << endl;
        double norma2 = norma(x1);
        //cout << norma2 << endl;

        cout << endl;
        cout << "delta = ";
        printf("%10.4e", norma1 / norma2);
        cout << endl;
    }

}