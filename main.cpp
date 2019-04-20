/********************
代码涉及Linux的系统调用
请在Linux下编译运行
*********************/


#include <iostream>
#include <iomanip>  // 用来格式化输出
#include <x86intrin.h>  // 会把本机支持的SSE、AVX库全部导入

#include <cstdlib>  // 用来产生随机数
#include <algorithm>
#include <sys/time.h>
using namespace std;


void make_matrix(float **Matrix, int N);
void show_matrix(float **Matrix, int N);
void gaussian_elimination(float **Matrix, int N);
void gaussian_elimination_sse_u(float **Matrix, int N);  // 未对齐
void gaussian_elimination_sse_a(float **Matrix, int N);  // 对齐
void gaussian_elimination_avx(float **Matrix, int N);
void copy_matrix(float** dst, float** src, int N);
void test();

int main()
{
//    struct timeval start;
//    struct timeval end;
//    unsigned long time_interval;
//
//    const int N = 100;
//    float Matrix[N][N];
//    float **Matrix_A= new float*[N];
//    for(int i = 0; i<N; i++)
//        Matrix_A[i] = new float[N];
//
//    float **Matrix_B= new float*[N];
//    for(int i = 0; i<N; i++)
//        Matrix_B[i] = new float[N];
//
//
//    make_matrix(Matrix_A, N);
////    cout<<"Matrix:"<<endl;show_matrix(Matrix_A, N);cout<<endl;
//
//    copy_matrix(Matrix_B, Matrix_A, N);
//
//    /**regular**/
//    gettimeofday(&start, NULL);
//    gaussian_elimination_sse_u(Matrix_A, N);
//    gettimeofday(&end, NULL);
//    cout << "u algorithm: " << endl;
////    show_matrix(Matrix_A, N);
//
//    time_interval = 1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
//    cout << "time_interval(regular) is " << time_interval << " us"<<endl;
//
//    /**SSE**/
//    gettimeofday(&start, NULL);
//    gaussian_elimination_avx(Matrix_B, N);
//    gettimeofday(&end, NULL);
//    cout << "a algorithm: " << endl;
////    show_matrix(Matrix_B, N);
//
//    time_interval = 1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
//    cout << "time_interval(SSE) is " << time_interval << " us"<<endl;
    test();

    return 0;

    float test[5];
    test[0] = 1;
    test[1] = 2;
    test[2] = 3;
    test[3] = 4;
    test[4] = 5;
//    test[3] = 4;
    __m128 __a = _mm_load_ps(test);
    float arg = 2;
    __m128 __b = _mm_set1_ps(arg);
    __m128 __c = _mm_div_ps(__a, __b);
    _mm_storeu_ps(test, __c);
    for(int i = 0; i<4; i++)
        cout<<test[i]<<endl;
}

void make_matrix(float **Matrix, int N)  // 为矩阵中各位值赋随机值
{
    srand((unsigned)time(NULL));  // 时间作种子
    for(int i = 0; i<N; i++)
        for(int j = 0; j<N; j++)
                Matrix[i][j] = rand()%10;  // 随机值取10以内
}
void show_matrix(float **Matrix, int N)  // 向控制台打印矩阵
{
    for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
                cout << fixed << setprecision(1) << setw(6) << right << Matrix[i][j];
            cout << endl;
        }
    cout << endl;
}
void copy_matrix(float** dst, float** src, int N)  // 把矩阵src的值赋给矩阵dst
{
    for(int i = 0; i<N; i++)
        for(int j = 0; j<N; j++)
            dst[i][j] = src[i][j];
}

void gaussian_elimination(float **Matrix, int N)  // 普通算法的高斯消去法
{
    for(int k = 0; k<N; k++)
    {
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            cout<<"akk = 0"<<endl;
            bool all_zero = false;
            for(int i = k+1; i<N; i++){
                if(0 != Matrix[i][k]){  // 发现第Matrix(i,k)不为0，就让第i行和第k行位置互换
                    for(int j = k; j<N; j++)
                        swap(Matrix[k][j], Matrix[i][j]);
                    cout << "row " << k << "and row" << i << " swapped："<<endl;
//                    show_matrix(Matrix, N); cout<<endl;
                    break;
                }
                else if(N-1 == i)  // k下面的每一行在第k列都是0
                        all_zero = true;
            }
            if(all_zero)
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

/**SSE未对齐的**/
void gaussian_elimination_sse_u(float **Matrix, int N)
/*
:param N为矩阵行数和列数（在这里矩阵均为方阵）
*/
{
    for(int k = 0; k<N; k++)
    {
        // 开始是解决Matrix(k,k)为0的问题
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            cout<<"akk = 0"<<endl;
            bool all_zero = false;
            for(int i = k+1; i<N; i++){
                if(0 != Matrix[i][k]){  // 发现第Matrix(i,k)不为0，就让第i行和第k行位置互换
                    for(int j = k; j<N; j++)
                        swap(Matrix[k][j], Matrix[i][j]);
                    cout << "row " << k << "and row " << i << " swapped："<<endl;
//                    show_matrix(Matrix, N); cout<<endl;
                    break;
                }
                else if(N-1 == i)  // k下面的每一行在第k列都是0
                        all_zero = true;
            }
            if(all_zero)
                continue;  // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
        }

        __m128 A_k_k = _mm_set_ps1(Matrix[k][k]);
//        __m128 A_k_k = _mm_load1_ps(Matrix[k]+k);
        for(int j = N-4; j>k; j-=4)
        {
            __m128 A_k_j = _mm_loadu_ps(Matrix[k]+j);
            A_k_j = _mm_div_ps(A_k_j, A_k_k);
            _mm_storeu_ps(Matrix[k]+j, A_k_j);
        }

        int part_2 = (N-k-1)%4;
        for(int j = k+1; j<k+1+part_2; j++)  // 不能被4整除的部分
            Matrix[k][j] = Matrix[k][j] / Matrix[k][k];

        Matrix[k][k] = 1.0;

        for(int i = k+1; i<N; i++)
        {
            __m128 A_i_k = _mm_set_ps1(Matrix[i][k]);
//             __m128 A_i_k = _mm_load1_ps(Matrix[i]+k);
            for(int j = N-4; j>k; j-=4)
                {
                    __m128 A_k_j = _mm_loadu_ps(Matrix[k]+j);
                    __m128 t = _mm_mul_ps(A_k_j, A_i_k);
                    __m128 A_i_j = _mm_loadu_ps(Matrix[i]+j);
                    A_i_j = _mm_sub_ps(A_i_j, t);
                    _mm_storeu_ps(Matrix[i]+j, A_i_j);
                }
            for(int j = k+1; j<k+1+part_2; j++)  // 不能被4整除的部分
                Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
            Matrix[i][k] = 0.0;
        }
//        if(k == 2)
//            show_matrix(Matrix, N);
    }
}
/*
笔记：
_mm_set_ps1，把寄存器中4个值全部设置为同一浮点数
*/

/**SSE对齐的**/
void gaussian_elimination_sse_a(float** Matrix, int N)
{
    for(int k = 0; k<N; k++)
    {
        // 开始是解决Matrix(k,k)为0的问题
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            cout<<"akk = 0"<<endl;
            bool all_zero = false;
            for(int i = k+1; i<N; i++){
                if(0 != Matrix[i][k]){  // 发现第Matrix(i,k)不为0，就让第i行和第k行位置互换
                    for(int j = k; j<N; j++)
                        swap(Matrix[k][j], Matrix[i][j]);
                    cout << "row " << k << "and row " << i << " swapped："<<endl;
//                    show_matrix(Matrix, N); cout<<endl;
                    break;
                }
                else if(N-1 == i)  // k下面的每一行在第k列都是0
                        all_zero = true;
            }
            if(all_zero)
                continue;  // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
        }

        __m128 A_k_k = _mm_set_ps1(Matrix[k][k]);
//        __m128 A_k_k = _mm_load1_ps(Matrix[k]+k);
        int k1 = k + 4- k%4;
        int n1 = N - N%4;
        for(int j = k1; j < n1; j+=4)
        {
            __m128 A_k_j = _mm_load_ps(Matrix[k]+j);
            A_k_j = _mm_div_ps(A_k_j, A_k_k);
            _mm_store_ps(Matrix[k]+j, A_k_j);
        }

        // 下面两个for循环计算未对齐的部分
        for(int j = k+1; j<k1; j++)
            Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
        for(int j = n1; j<N; j++)
            Matrix[k][j] = Matrix[k][j] / Matrix[k][k];

        Matrix[k][k] = 1.0;

        for(int i = k+1; i<N; i++)
        {
            __m128 A_i_k = _mm_set_ps1(Matrix[i][k]);
//             __m128 A_i_k = _mm_load1_ps(Matrix[i]+k);
            for(int j = k1; j < n1; j+=4)
                {
                    __m128 A_k_j = _mm_load_ps(Matrix[k]+j);
                    __m128 t = _mm_mul_ps(A_k_j, A_i_k);
                    __m128 A_i_j = _mm_load_ps(Matrix[i]+j);
                    A_i_j = _mm_sub_ps(A_i_j, t);
                    _mm_store_ps(Matrix[i]+j, A_i_j);
                }
            for(int j = k+1; j<k1; j++)  // 不能被4整除的部分
                Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
            for(int j = n1; j<N; j++)
                Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
            Matrix[i][k] = 0.0;
        }
//        if(k == 2)
//            show_matrix(Matrix, N);
    }
}

void gaussian_elimination_avx(float **Matrix, int N)
{
    for(int k = 0; k<N; k++)
    {
        // 开始是解决Matrix(k,k)为0的问题
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            cout<<"akk = 0"<<endl;
            bool all_zero = false;
            for(int i = k+1; i<N; i++){
                if(0 != Matrix[i][k]){  // 发现第Matrix(i,k)不为0，就让第i行和第k行位置互换
                    for(int j = k; j<N; j++)
                        swap(Matrix[k][j], Matrix[i][j]);
                    cout << "row " << k << "and row " << i << " swapped："<<endl;
//                    show_matrix(Matrix, N); cout<<endl;
                    break;
                }
                else if(N-1 == i)  // k下面的每一行在第k列都是0
                        all_zero = true;
            }
            if(all_zero)
                continue;  // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
        }

        __m256 A_k_k = _mm256_set1_ps(Matrix[k][k]);
//        __m128 A_k_k = _mm_load1_ps(Matrix[k]+k);
        for(int j = N-8; j>k; j-=8)
        {
            __m256 A_k_j = _mm256_loadu_ps(Matrix[k]+j);
            A_k_j = _mm256_div_ps(A_k_j, A_k_k);
            _mm256_storeu_ps(Matrix[k]+j, A_k_j);
        }

        int part_2 = (N-k-1)%8;
        for(int j = k+1; j<k+1+part_2; j++)  // 不能被4整除的部分
            Matrix[k][j] = Matrix[k][j] / Matrix[k][k];

        Matrix[k][k] = 1.0;

        for(int i = k+1; i<N; i++)
        {
            __m256 A_i_k = _mm256_set1_ps(Matrix[i][k]);
//             __m128 A_i_k = _mm_load1_ps(Matrix[i]+k);
            for(int j = N-8; j>k; j-=8)
                {
                    __m256 A_k_j = _mm256_loadu_ps(Matrix[k]+j);
                    __m256 t = _mm256_mul_ps(A_k_j, A_i_k);
                    __m256 A_i_j = _mm256_loadu_ps(Matrix[i]+j);
                    A_i_j = _mm256_sub_ps(A_i_j, t);
                    _mm256_storeu_ps(Matrix[i]+j, A_i_j);
                }
            for(int j = k+1; j<k+1+part_2; j++)  // 不能被4整除的部分
                Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
            Matrix[i][k] = 0.0;
        }
//        if(k == 2)
//            show_matrix(Matrix, N);
    }
}
void test()
{
    struct timeval start;
    struct timeval end;
    unsigned long time_interval;
    int N = 100000;
    float *data = new float[N];
    for(int i=0; i<N; i++)
        data[i] = i;
    float data2[1024];
    for(int i=0; i<1024; i++)
        data2[i] = i;
    __m128 t_128 = _mm_set_ps1(2.0);

    gettimeofday(&start, NULL);
    for(int i = 0; i<N; i+=4)
        {
            __m128 a_128 = _mm_loadu_ps(data+i);
            a_128 = _mm_div_ps(a_128, t_128);
        }
    gettimeofday(&end, NULL);
    time_interval = 1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    cout << "time_interval(128) is " << time_interval << " us"<<endl;


    __m256 t_256 = _mm256_set1_ps(2.0);

    gettimeofday(&start, NULL);
    for(int i = 0; i<N; i+=8)
        {
            __m256 a_256 = _mm256_loadu_ps(data+i);
            a_256 = _mm256_div_ps(a_256, t_256);
        }
    gettimeofday(&end, NULL);
    time_interval = 1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    cout << "time_interval(256) is " << time_interval << " us"<<endl;
}

