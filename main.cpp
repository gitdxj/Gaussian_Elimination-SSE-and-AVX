/********************
代码涉及Linux的系统调用
请在Linux下编译运行
*********************/
#include <iostream>
#include <iomanip>  // 用来格式化输出
#include <x86intrin.h>  // 会把本机支持的SSE、AVX库全部导入
#include <fstream>
#include <cstdlib>  // 用来产生随机数
#include <algorithm>
#include <sys/time.h>
using namespace std;




void make_matrix(float **Matrix, int N);
void show_matrix(float **Matrix, int N);
void copy_matrix(float** dst, float** src, int N);
bool ls_same(float **a, float**b, int N);
bool swap_rows(float **Matrix, int N, int k);
void gaussian_elimination(float **Matrix, int N);
void gaussian_elimination_sse_unaligned(float **Matrix, int N, bool p45);
void gaussian_elimination_sse_aligned(float **Matrix, int N, bool p45);
void gaussian_elimination_avx_unaligned(float **Matrix, int N, bool p45);
void gaussian_elimination_avx_aligned(float **Matrix, int N, bool p45);

void test();
void test(int N, unsigned long* time_interval);


int main()
{
//    struct timeval start;
//    struct timeval end;
//    unsigned long time_interval;
//
//    const int N = 100;
//    float **Matrix_A= new float*[N];
//    for(int i = 0; i<N; i++)
//        Matrix_A[i] = new float[N];
//
//    float **Matrix_B= new float*[N];
//    for(int i = 0; i<N; i++)
//        Matrix_B[i] = new float[N];
//
//    float **Matrix= new float*[N];
//    for(int i = 0; i<N; i++)
//        Matrix[i] = new float[N];
//
//
//    make_matrix(Matrix_A, N);
////    cout<<"Matrix:"<<endl;show_matrix(Matrix_A, N);cout<<endl;
//
//    copy_matrix(Matrix, Matrix_A, N);
//    copy_matrix(Matrix_B, Matrix_A, N);
//
//    /**regular**/
//    gettimeofday(&start, NULL);
//    gaussian_elimination(Matrix_A, N);
//    gettimeofday(&end, NULL);
//    cout << "u algorithm: " << endl;
////    show_matrix(Matrix_A, N);
//
//    time_interval = 1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
//    cout << "time_interval(regular) is " << time_interval << " us"<<endl;
//
//
//
//    /**SSE**/
//    gettimeofday(&start, NULL);
//    gaussian_elimination_avx_aligned(Matrix_B, N, true);
//    gettimeofday(&end, NULL);
//    cout << "a algorithm: " << endl;
////    show_matrix(Matrix_B, N);
//
//    time_interval = 1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
//    cout << "time_interval(SSE) is " << time_interval << " us"<<endl;


    int C = 9;
    double **data = new double*[C];
    for(int i = 0; i<C; i++)
        data[i] = new double[9];
    for(int i = 0; i<C; i++)
        for(int k = 0; k<9; k++)
            data[i][k] = 0;
    unsigned long *temp = new unsigned long[9];
    int N = 8;
    int iterate_num = 5;
    for(int i = 0; i<C; i++)
    {

        for(int j = 0; j<iterate_num; j++) //运行5次
        {
            test(N, temp);
            for(int k = 0; k<9; k++)
                data[i][k] += temp[k]/(double)iterate_num;
        }

        N = N * 2;
    }

    ofstream outFile;
    outFile.open("data.csv", ios::out);
    outFile << ','<<"regular LU     : "
            << ','<<"SSE unaligned S: "
            << ','<<"SSE unaligned P: "
            << ','<<"SSE   aligned S: "
            << ','<<"SSE   aligned P: "
            << ','<<"AVX unaligned S: "
            << ','<<"AVX unaligned P: "
            << ','<<"AVX   aligned S: "
            << ','<<"AVX   aligned P: " <<','<<endl;
    N = 8;
    for(int i = 0; i<C; i++)
    {
        outFile << N << ',';
        for(int j = 0; j<9; j++)
            outFile<<data[i][j]<<',';
        outFile<<endl;
        N = N * 2;
    }
    outFile.close();
    return 0;


}
/** 为矩阵中各位值赋随机值**/
void make_matrix(float **Matrix, int N)
{
    srand((unsigned)time(NULL));  // 时间作种子
    for(int i = 0; i<N; i++)
        for(int j = 0; j<N; j++)
                Matrix[i][j] = rand()%10;  // 随机值取10以内
}

/**向控制台打印矩阵**/
void show_matrix(float **Matrix, int N)
{
    for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++){
                   cout << fixed/*以小数形式输出（不用科学计数法）*/
                        << setprecision(1)/*保留小数点后一位*/
                        << setw(6)/*指定输出宽度为6，不足用空格补齐*/
                        << right/*向右对齐*/ << Matrix[i][j];
            }
            cout << endl;
        }
    cout << endl;
}

/**判断两个矩阵的值是否相同**/
bool ls_same(float **a, float**b, int N)
{
    for(int i = 0; i<N; i++)
        for(int j = 0; j<N; j++)
            if(a[i][j] != b[i][j])
                return false;
    return true;
}
void copy_matrix(float** dst, float** src, int N)  // 把矩阵src的值赋给矩阵dst
{
    for(int i = 0; i<N; i++)
        for(int j = 0; j<N; j++)
            dst[i][j] = src[i][j];
}


/**
当Matrix(k,k)为0的时候，从k行向下找出一行i使得Matrix(i,k)不为0，将i行和k行互换
因为当Matrix(k,k)为0时，会导致除数为0问题
param N:意味矩阵为大小为N*N
param k:出现Matrix(k,k)==0时的k
若进行了两行间的互换则返回true
未进行互换（意味着从k行向下的每一行在k列的位置都为0）则返回false
**/
bool swap_rows(float **Matrix, int N, int k)
{
//    cout<<"Matrix("<<k<<","<<k<<") = 0"<<endl;
    for(int i = k+1; i<N; i++){
        if(0 != Matrix[i][k]){  // 发现第Matrix(i,k)不为0，就让第i行和第k行位置互换
            for(int j = k; j<N; j++)
                swap(Matrix[k][j], Matrix[i][j]);
//            cout << "row " << k << " and row " << i << " swapped："<<endl;
//            show_matrix(Matrix, N); cout<<endl;
            return true;
        }
        else if(N-1 == i)  // k下面的每一行在第k列都是0
                return false;
    }
}


/**普通的LU高斯消去法**/
void gaussian_elimination(float **Matrix, int N)
{
    for(int k = 0; k<N; k++)
    {
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            bool swapped = swap_rows(Matrix, N, k);

            if(!swapped)   // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
                continue;
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
void gaussian_elimination_sse_unaligned(float **Matrix, int N, bool p45)
{
    for(int k = 0; k<N; k++)
    {
        // 开始是解决Matrix(k,k)为0的问题
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            bool swapped = swap_rows(Matrix, N, k);

            if(!swapped)   // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
                continue;
        }

        __m128 A_k_k = _mm_set_ps1(Matrix[k][k]);
//        __m128 A_k_k = _mm_load1_ps(Matrix[k]+k);
        int part_2 = (N-k-1)%4;
        if(p45){  // p45为true，则把45行的循环向量化
            for(int j = N-4; j>k; j-=4)
            {
                __m128 A_k_j = _mm_loadu_ps(Matrix[k]+j);
                A_k_j = _mm_div_ps(A_k_j, A_k_k);
                _mm_storeu_ps(Matrix[k]+j, A_k_j);
            }
            for(int j = k+1; j<k+1+part_2; j++)  // 不能被4整除的部分
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
            }
        else
            for(int j = k+1; j<N; j++)
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
void gaussian_elimination_sse_aligned(float** Matrix, int N, bool p45)
{
    for(int k = 0; k<N; k++)
    {
        // 开始是解决Matrix(k,k)为0的问题
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            bool swapped = swap_rows(Matrix, N, k);

            if(!swapped)   // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
                continue;
        }

        __m128 A_k_k = _mm_set_ps1(Matrix[k][k]);
//        __m128 A_k_k = _mm_load1_ps(Matrix[k]+k);
//        int k1 = k + 4- k%4;
//        int n1 = N - N%4;
        long matrix_start_addr = (long)(&Matrix[k][0]);
        int offset = (matrix_start_addr%16)/4;
//        int k1 = k - ((k-1)*N+k+offset)%4 +4;
//        int n1 = N - (k*N+offset)%4;
        int k1 = k - (k+offset)%4 +4;
        int n1 = N - (N+offset)%4;
        if(p45){
            for(int j = k1; j < n1; j+=4)
            {
                __m128 A_k_j = _mm_load_ps(Matrix[k]+j);
                A_k_j = _mm_div_ps(A_k_j, A_k_k);
                _mm_store_ps(Matrix[k]+j, A_k_j);
            }
            //下面两个for循环计算不能被4整除的部分
            for(int j = k+1; j<k1; j++)
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
            for(int j = n1; j<N; j++)
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
        }
        else
            for(int j = k+1; j<N; j++)
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
    }
}

/**AVX未对齐的**/
void gaussian_elimination_avx_unaligned(float **Matrix, int N, bool p45)
{
    for(int k = 0; k<N; k++)
    {
        // 开始是解决Matrix(k,k)为0的问题
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            bool swapped = swap_rows(Matrix, N, k);

            if(!swapped)   // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
                continue;
        }

        __m256 A_k_k = _mm256_set1_ps(Matrix[k][k]);

        int part_2 = (N-k-1)%8;

        if(p45){  // 若4到5行的代码向量化
            for(int j = N-8; j>k; j-=8)
            {
                __m256 A_k_j = _mm256_loadu_ps(Matrix[k]+j);
                A_k_j = _mm256_div_ps(A_k_j, A_k_k);
                _mm256_storeu_ps(Matrix[k]+j, A_k_j);
            }
            for(int j = k+1; j<k+1+part_2; j++)  // 不能被8整除的部分
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
        }
        else  // 若4到5行的代码不向量化
            for(int j = k+1; j<N; j++)
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];

        Matrix[k][k] = 1.0;

        for(int i = k+1; i<N; i++)
        {
            __m256 A_i_k = _mm256_set1_ps(Matrix[i][k]);

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

    }
}

/**AVX对齐的**/
void gaussian_elimination_avx_aligned(float **Matrix, int N, bool p45)
{
    for(int k = 0; k<N; k++)
    {
        // 开始是解决Matrix(k,k)为0的问题
        if(0 == Matrix[k][k])  // 如果A(k,k)的位置为0的话，就从后面找一行不为0的互换
        {
            bool swapped = swap_rows(Matrix, N, k);

            if(!swapped)   // 如果下面任何一行的第k列都没有不是0打头的就直接跳下一个k
                continue;
        }

        __m256 A_k_k = _mm256_set1_ps(Matrix[k][k]);

//        int k1 = k + 8- k%8;
//        int n1 = N - N%8;
        long matrix_start_addr = (long)(&Matrix[k][0]);
        int offset = (matrix_start_addr%32)/4;
//        cout<<"offset: "<<offset<<endl;
//        int k1 = k - (k*N+k+offset)%8 + 8;
//        int n1 = N - (k*N+N+offset)%8;
        int k1 = k - (k+offset)%8 + 8;
        int n1 = N - (N+offset)%8;
        if(p45){
            for(int j = k1; j < n1; j+=8)
            {
//                cout<<"("<<k<<","<<j<<")"<<endl;
//                cout<<&Matrix[k][j]<<endl;
//                long addr = (long)&Matrix[k][j];
//                if(0 != addr%32)
//                    cout<<"unaligned"<<endl;

                __m256 A_k_j = _mm256_load_ps(Matrix[k]+j);
                A_k_j = _mm256_div_ps(A_k_j, A_k_k);
                _mm256_store_ps(Matrix[k]+j, A_k_j);
            }
            //下面两个for循环计算不能被8整除的部分
            for(int j = k+1; j<k1; j++)
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
            for(int j = n1; j<N; j++)
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];
        }
        else
            for(int j = k+1; j<N; j++)
                Matrix[k][j] = Matrix[k][j] / Matrix[k][k];


        Matrix[k][k] = 1.0;

        for(int i = k+1; i<N; i++)
        {
            __m256 A_i_k = _mm256_set1_ps(Matrix[i][k]);

            for(int j = k1; j < n1; j+=8)
                {
                    __m256 A_k_j = _mm256_load_ps(Matrix[k]+j);
                    __m256 t = _mm256_mul_ps(A_k_j, A_i_k);
                    __m256 A_i_j = _mm256_load_ps(Matrix[i]+j);
                    A_i_j = _mm256_sub_ps(A_i_j, t);
                    _mm256_store_ps(Matrix[i]+j, A_i_j);
                }
            for(int j = k+1; j<k1; j++)  // 不能被4整除的部分
                Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
            for(int j = n1; j<N; j++)
                Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
            Matrix[i][k] = 0.0;
        }
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

    gettimeofday(&start, NULL);
    float t_32 = 2;
    for(int i = 0; i<N; i++)
    {
        float a_32 = data[i];
        a_32 = a_32/t_32;
    }
    gettimeofday(&end, NULL);
    time_interval = 1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    cout << "time_interval(32) is " << time_interval << " us"<<endl;

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

void test(int N, unsigned long* time_interval)
{

    const int matrix_num = 9;
    float ***Matrix;

    Matrix = new float**[matrix_num]; // 9个矩阵代表对应SSE/AVX，对齐/非对齐，45行是否向量化的八种情况加上平凡的算法
    for(int i = 0; i<matrix_num; i++)
    {
        Matrix[i] = new float*[N];
        for(int j = 0; j<N; j++)
            if(i <= 4)
                Matrix[i][j] = (float*)_mm_malloc(sizeof(float)*N, 16);
            else
                Matrix[i][j] = (float*)_mm_malloc(sizeof(float)*N, 32);
//            Matrix[i][j] = new float[N];
    }
    // 为这9个矩阵赋相同的值
    make_matrix(Matrix[0], N);
    for(int i = 1; i<matrix_num; i++)
        copy_matrix(Matrix[i], Matrix[0], N);




    struct timeval start;
    struct timeval finish;
    unsigned long time;

    gettimeofday(&start, NULL);
    gaussian_elimination(Matrix[0], N);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[0] = time;
    cout << "regular LU     : " << setw(6) << right << time  << " us"<<endl;

    gettimeofday(&start, NULL);
    gaussian_elimination_sse_unaligned(Matrix[1], N, false);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[1] = time;
    cout << "SSE unaligned S: " << setw(6) << right << time  << " us"<<endl;


    gettimeofday(&start, NULL);
    gaussian_elimination_avx_unaligned(Matrix[2], N, true);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[2] = time;
    cout << "SSE unaligned P: " << setw(6) << right << time  << " us"<<endl;


    gettimeofday(&start, NULL);
    gaussian_elimination_sse_aligned(Matrix[3], N, false);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[3] = time;
    cout << "SSE   aligned S: " << setw(6) << right << time  << " us"<<endl;

    gettimeofday(&start, NULL);
    gaussian_elimination_sse_aligned(Matrix[4], N, true);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[4] = time;
    cout << "SSE   aligned P: " << setw(6) << right << time  << " us"<<endl;

    gettimeofday(&start, NULL);
    gaussian_elimination_avx_unaligned(Matrix[5], N, false);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[5] = time;
    cout << "AVX unaligned S: " << setw(6) << right << time  << " us"<<endl;

    gettimeofday(&start, NULL);
    gaussian_elimination_avx_unaligned(Matrix[6], N, true);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[6] = time;
    cout << "AVX unaligned P: " << setw(6) << right << time  << " us"<<endl;

//    cout << &Matrix[7][0][0]<<endl;
    gettimeofday(&start, NULL);
    gaussian_elimination_avx_aligned(Matrix[7], N, false);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[7] = time;
    cout << "AVX   aligned S: " << setw(6) << right << time << " us"<<endl;

    gettimeofday(&start, NULL);
    gaussian_elimination_avx_aligned(Matrix[8], N, true);
    gettimeofday(&finish, NULL);
    time = 1000000*(finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    time_interval[8] = time;
    cout << "AVX   aligned P: " << setw(6) << right << time  << " us"<<endl;

//    for(int i = 1; i<matrix_num; i++)
//        if(!ls_same(Matrix[0], Matrix[i], N))
//            cout << "wrong outcome @ " << "i = "<<i<<endl;

    // 回收内存
    for(int i = 0; i<matrix_num; i++)
        for(int j = 0; j<N; j++)
            delete[] Matrix[i][j];
    for(int i = 0; i<matrix_num; i++)
        delete[] Matrix[i];
    delete[] Matrix;
}

