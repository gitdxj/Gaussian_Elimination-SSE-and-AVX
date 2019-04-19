#include <iostream>
#include <iomanip>
#include <x86intrin.h>
#include <cstdlib>
#include <algorithm>
using namespace std;

void gaussian_elimination(float **Matrix, int N);
void make_matrix(float **Matrix, int N);
void show_matrix(float **Matrix, int N);

int main()
{
//    float *a = new float[4];
//    float *b = new float[4];
//    for(int i = 0; i<4; i++)
//        a[i] = b[i] = i;
//    float c;
//    __m128 t1, t2, sum;
//    t1 = _mm_loadu_ps(a);
//    t2 = _mm_loadu_ps(b);
//    t1 = _mm_mul_ps(t1,t2);
//    sum = _mm_setzero_ps();
//    sum = _mm_add_ps(sum,t1);
//    sum = _mm_hadd_ps(sum,sum);
//    sum = _mm_hadd_ps(sum,sum);
//    _mm_store_ss(&c, sum);
//    cout<<c<<endl;
    int N = 10;
    float **Matrix= new float*[N];
    for(int i = 0; i<N; i++)
        Matrix[i] = new float[N];

    make_matrix(Matrix, N);
    show_matrix(Matrix, N);cout<<endl;
    gaussian_elimination(Matrix, N);
    show_matrix(Matrix, N);

    return 0;
}


void gaussian_elimination_sse(float **Matrix, int N)
/*
para: N为矩阵行数和列数（在这里矩阵均为方阵）

*/
{
    for(int k = 0; k<N; k++)
    {
        __m128 a_k_k = _mm_set_ps1(Matrix[k][k]);
    }
}
void gaussian_elimination(float **Matrix, int N)
/*
1. procedure LU (A)
2. begin
3. 	for k := 1 to n do
4. 		for j := k+1 to n do
5. 			A[k, j] := A[k, j]/A[k, k];
6		A[k, k] := 1.0;
7. 		endfor;
8. 		for i := k + 1 to n do
9. 			for j := k + 1 to n do
10. 				A[i, j] := A[i, j] - A[i, k] × A[k, j ];
11. 		endfor;
12. 		A[i, k] := 0;
13. 	endfor;
14. endfor;
15. end LU
*/
{
    for(int k = 0; k<N; k++)
    {
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            for(int i = k+1; k<N; k++)
                if(0 != Matrix[i][k])
                    for(int j = k; j<N; j++)
                        swap(Matrix[k][j], Matrix[i][j]);
            continue;  // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
        }
        for(int j = k+1; j<N; j++)
            Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
        Matrix[k][k] = 1.0;
        for(int i = k+1; i<N; i++){
            for(int j = k+1; j<N; j++)
                Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
            Matrix[i][k] = 0;
        }
    }

}
void make_matrix(float **Matrix, int N)
{
    srand((unsigned)time(NULL));
    for(int i = 0; i<N; i++)
        for(int j = 0; j<N; j++)
                Matrix[i][j] = rand()%10;

    Matrix[0][0] = 0.0;
}
void show_matrix(float **Matrix, int N)
{
    for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
                cout <<fixed<<setprecision(1) <<setw(6)<<left<< Matrix[i][j];
            cout << endl;
        }
}
