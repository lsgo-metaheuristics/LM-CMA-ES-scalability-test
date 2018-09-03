#include "math.h"
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include "funsoft.h"

#define  _CRT_SECURE_NO_WARNINGS_GLOBALS 1

using namespace std;

/* 	LM-CMA-ES by Ilya Loshchilov. 2014 (c)
*/

typedef struct
/* random_t
 * sets up a pseudo random number generator instance
 */
{
    /* Variables for Uniform() */
    long int startseed;
    long int aktseed;
    long int aktrand;
    long int *rgrand;

    /* Variables for Gauss() */
    short flgstored;
    double hold;
} random_t;

typedef struct
{
    double value;
    int id;
} sortedvals;

typedef struct
{
    random_t ttime;
    double*	func_tempdata;
    double*	x_tempdata;
    double*	rotmatrix;
    double* func_shiftxi;
    time_t	time_tic_t;
    clock_t	time_tic_c;
    time_t	time_toc_t;
    clock_t	time_toc_c;
} global_t;



typedef Fitness(*FunctionCallback)(int d, double *x);

FunctionCallback functions[] = { Shifted_Sphere, Schwefel_Problem, Shifted_Rosenbrock, Shifted_Rastrigin,
                                 Shifted_Griewank, Shifted_Ackley, f_Schwefel2_22, f_Schwefel1_2, Extended_f_10, f_Bohachevsky, f_Schaffer, f_Hybrid_12,
                                 f_Hybrid_13, f_Hybrid_14, f_Hybrid_15, f_Hybrid_16new, f_Hybrid_17new, f_Hybrid_18new, f_Hybrid_19new, f_FastFractal
                               };

char* functionsNames[] = { "Shifted_Sphere", "Schwefel_Problem", "Shifted_Rosenbrock", "Shifted_Rastrigin",
                           "Shifted_Griewank", "Shifted_Ackley", "Schwefel2_22", "Schwefel1_2", "Extended_f_10", "Bohachevsky", "Schaffer", "Hybrid_12",
                           "Hybrid_13", "Hybrid_14", "Hybrid_15", "Hybrid_16", "Hybrid_17", "Hybrid_18", "Hybrid_19", "FastFractal"
                         };

double functionsDomainBounds[] = { 100.0,  100.0, 100.0,    5.0,  600.0,   32.0, 10.0, 65.536, 100.0, 15.0, 100.0, 100.0, 100.0, 5.0, 10.0, 100.0, 100.0, 5.0, 10.0, 1.0 };
Fitness  functionsOptima[] = { -450.0, -450.0, 390.0, -330.0, -180.0, -140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };



void init_gt(global_t* gt)
{
    gt->func_tempdata = NULL;
    gt->x_tempdata = NULL;
    gt->rotmatrix = NULL;
    gt->func_shiftxi = NULL;
}

void free_gt(global_t* gt)
{
    if (gt->func_tempdata)	{
        delete[] gt->func_tempdata;
        gt->func_tempdata = NULL;
    }
    if (gt->x_tempdata)		{
        delete[] gt->x_tempdata;
        gt->x_tempdata = NULL;
    }
    if (gt->rotmatrix)		{
        delete[] gt->rotmatrix;
        gt->rotmatrix = NULL;
    }
    if (gt->func_shiftxi)	{
        delete[] gt->func_shiftxi;
        gt->func_shiftxi = NULL;
    }
}



/* random_Start(), random_init(), random_exit(), random_Unifo()rm, random_Gauss(), time_tic(), time_tictoc(), time_toc() are adopted
   from Nikolaus Hansen's source code for CMA-ES
*/

void random_exit(random_t *t)
{
    free( t->rgrand);
}

long random_Start( random_t *t, long unsigned inseed)
{
    long tmp;
    int i;

    t->flgstored = 0;
    t->startseed = inseed;
    if (inseed < 1)
        inseed = 1;
    t->aktseed = inseed;
    for (i = 39; i >= 0; --i)
    {
        tmp = t->aktseed/127773;
        t->aktseed = 16807 * (t->aktseed - tmp * 127773)
                     - 2836 * tmp;
        if (t->aktseed < 0) t->aktseed += 2147483647;
        if (i < 32)
            t->rgrand[i] = t->aktseed;
    }
    t->aktrand = t->rgrand[0];
    return inseed;
}

long random_init( random_t *t, long unsigned inseed)
{
    clock_t cloc = clock();

    t->flgstored = 0;
    t->rgrand = new long[32];
    if (inseed < 1) {
        while ((long) (cloc - clock()) == 0)
            ;
        inseed = (long)abs(100.0*time(NULL)+clock());
    }
    return random_Start(t, inseed);
}

unsigned int m_z = 1;
unsigned int m_w = 123;
unsigned int GetUint()
{
    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
    return (m_z << 16) + m_w;
}

double GetUniform()
{
    // 0 <= u < 2^32
    unsigned int u = GetUint();
    // The magic number below is 1/(2^32 + 2).
    // The result is strictly between 0 and 1.
    return (u + 1.0) * 2.328306435454494e-10;
}

int ibit = 31;
int bits[32];
int myrand = 0;
int mask = 1;

void GenerateNewBits()
{
    unsigned int u = GetUint();
    int mask = 1;
    for (int i = 0; i < 31; i++)
    {
        if ((mask&u) >= 1)	bits[i] = 1;
        else				bits[i] = 0;
        mask <<= 1;
        //	printf("%d ", bits[i]);
    }
    ibit = 0;
}

int GetBinary()
{
    unsigned int u = GetUint();
    if (u < 1073741824)	return 0;
    else				return 1;
}

int GetBinaryCheap()
{
    /*	if (ibit == 31)	GenerateNewBits();
    	ibit++;
    	return bits[ibit - 1];*/
    if (ibit == 31)
    {
        mask = 1;
        myrand = GetUint();
    }
    int res = 0;
    if ((mask&myrand) >= 1)	res = 1;
    mask <<= 1;
    return res;
}

int GetRademacherCheap()
{
    //if (ibit == 31)
    {
        mask = 1;
        myrand = GetUint();
    }
    int res = -1;
    if ((mask&myrand) >= 1)	res = 1;
    mask <<= 1;
    return res;
}

double random_Uniform( random_t *t)
{
    return GetUniform();
    long tmp;

    tmp = t->aktseed/127773;
    t->aktseed = 16807 * (t->aktseed - tmp * 127773)
                 - 2836 * tmp;
    if (t->aktseed < 0)
        t->aktseed += 2147483647;
    tmp = t->aktrand / 67108865;
    t->aktrand = t->rgrand[tmp];
    t->rgrand[tmp] = t->aktseed;
    return (double)(t->aktrand)/(2.147483647e9);
}

/* --------------------------------------------------------- */
double random_Gauss(random_t *t)
{
    double x1, x2, rquad, fac;

    if (t->flgstored)
    {
        t->flgstored = 0;
        return t->hold;
    }
    do
    {
        x1 = 2.0 * random_Uniform(t) - 1.0;
        x2 = 2.0 * random_Uniform(t) - 1.0;
        rquad = x1*x1 + x2*x2;
    } while(rquad >= 1 || rquad <= 0);
    fac = sqrt(-2.0*log(rquad)/rquad);
    t->flgstored = 1;
    t->hold = fac * x1;
    return fac * x2;
}

void	time_tic(global_t* t)
{
    t->time_tic_t = time(NULL);	// measure time in seconds
    t->time_tic_c = clock();	// measure time in microseconds up to ~2k seconds
}

double	time_tictoc(global_t* t)
{
    double dt = difftime(t->time_toc_t,t->time_tic_t);
    if (dt < 1000)
        dt = (double)(t->time_toc_c - t->time_tic_c) / CLOCKS_PER_SEC;
    return dt;
}

double	time_toc(global_t* t)
{
    t->time_toc_t = time(NULL);
    t->time_toc_c = clock();
    return time_tictoc(t);
}

void matrix_eye(double* M, int m)
{
    for (int i=0; i<m; i++)
        for (int j=0; j<m; j++)
            if (i == j)	M[i*m + j] = 1.0;
            else		M[i*m + j] = 0.0;
}


// vector res = matrix a X vector b
void matrix_mult_vector(double* res, double* a, double* b,int m)
{
    double val = 0.0;
    for (int i=0; i<m; i++)
    {
        val = 0.0;
        for (int j=0; j<m; j++)
            val += a[i*m + j] * b[j];
        res[i] = val;
    }
}

// maxtrix res = matrix a X matrix b
void matrix_mult_matrix(double* res, double* a, double* b,int m)
{
    double val;
    for (int i=0; i<m; i++)
        for (int j=0; j<m; j++)
        {
            val = 0;
            for (int k=0; k<m; k++)
                val += a[i*m + k] * b[k*m + j];
            res[i*m + j] = val;
        }
}

// matrix res = vector a X vector b
void vector_mult_vector(double* res, double* a, double* b,int m)
{
    for (int i=0; i<m; i++)
        for (int j=0; j<m; j++)
            res[i*m + j] = a[i] * b[j];
}

// vector res = vector a X matrix b
void vector_mult_matrix(double* res, double* a, double* b,int m)
{
    double val;
    for (int i=0; i<m; i++)
    {
        val = 0;
        for (int j=0; j<m; j++)
            val += a[j] * b[j*m + i];
        res[i] = val;
    }
}

double vector_prod(double* a, double* b,int m)
{
    double res = 0.0;
    for (int i=0; i<m; i++)
        res += a[i] * b[i];
    return res;
}

void generateRotationMatrix(double* B, int N, double* tmp1, random_t* rnd)
{
    double* pB;

    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            B[i*N + j] = random_Gauss(rnd);
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<i; j++)
        {
            double ariarj = 0;
            for (int k=0; k<N; k++)
                ariarj = ariarj + B[k*N+i]*B[k*N+j];

            for (int k=0; k<N; k++)
                B[k*N+i] = B[k*N+i] - ariarj * B[k*N+j];
        }
        double normv = 0;
        for (int k=0; k<N; k++)
            normv = normv + B[k*N+i]*B[k*N+i];

        normv = sqrt(normv);
        for(int k=0; k<N; k++)
            B[k*N+i] = B[k*N+i] / normv;
    }
}

double minv(double a, double b)
{
    if (a < b)	return a;
    else		return b;
}

double maxv(double a, double b)
{
    if (a > b)	return a;
    else		return b;
}

double fsphere(double* x, int N)
{
    double Fit = 0;
    for (int i=0; i<N; i++)
        Fit += x[i] * x[i];
    return Fit;
}

double felli(double* x, int N)
{
    double Fit = 0;
    double alpha = pow(10,6.0);
    for (int i=0; i<N; i++)
        Fit += pow(alpha, double(i) / double(N-1) ) * x[i] * x[i];
    return Fit;
}

double felli_fast(double* x, int N, global_t* t)
{
    double Fit = 0;
    if (t->func_tempdata == NULL)
    {
        t->func_tempdata = new double[N];
        double alpha = pow(10,6.0);
        for (int i=0; i<N; i++)
            t->func_tempdata[i] = pow(alpha, double(i) / double(N-1) );
    }

    for (int i=0; i<N; i++)
        Fit += t->func_tempdata[i] * x[i] * x[i];
    return Fit;
}

double fdiscus(double* x, int N)
{
    double Fit = 0;
    Fit = 1e+6 * (x[0] * x[0]);
    for (int i=1; i<N; i++)
        Fit += x[i]*x[i];
    return Fit;
}

double fcigar(double* x, int N)
{
    double Fit = 0;
    for (int i=1; i<N; i++)
        Fit += x[i]*x[i];
    Fit = Fit * 1e+6;
    Fit += x[0] * x[0];
    return Fit;
}

double fdiffpowers_fast(double* x, int N, global_t* t)
{
    double Fit = 0;
    if (t->func_tempdata == NULL)
    {
        t->func_tempdata = new double[N];
        for (int i = 0; i<N; i++)
            t->func_tempdata[i] = double(2.0 + (4.0 * (double(i)) / (N - 1)));
    }
    for (int i = 0; i<N; i++)
        Fit += pow(fabs(x[i]), t->func_tempdata[i]);
    //Fit = sqrt(Fit);
    return Fit;
}

void getRotatedX(double* x, int N, global_t* t)
{
    if (t->x_tempdata == NULL)
        t->x_tempdata = new double[N];
    if (t->rotmatrix == NULL)
    {
        t->rotmatrix = new double[N*N];
        generateRotationMatrix(t->rotmatrix, N, t->x_tempdata, &t->ttime);
    }
    matrix_mult_vector(t->x_tempdata, t->rotmatrix, x, N);
}

double frosen(double* x, int N)
{
    double Fit = 0;
    double tmp1, tmp2;
    double Fit1 = 0;
    double Fit2 = 0;
    //for (int i=0; i<N-1; i++)
    //	Fit += 100 * pow( x[i]*x[i] - x[i+1], 2.0  ) + pow(x[i] - 1.0, 2.0); // function 'pow' is very slow
    for (int i=0; i<N-1; i++)
    {
        tmp1 = x[i]*x[i] - x[i+1];
        tmp2 = x[i] - 1.0;
        Fit1 += tmp1*tmp1;
        Fit2 += tmp2*tmp2;
    }
    Fit = 100*Fit1 + Fit2;
    return Fit;
}





Fitness MyFunc(int FuncId, int N, double* x, global_t* t)
{

    Fitness f = functions[FuncId - 1](N, x);

    return f;
}



int compare(const void * a, const void * b)
{
    if (((sortedvals*)a)->value < ((sortedvals*)b)->value)	return -1;
    if (((sortedvals*)a)->value == ((sortedvals*)b)->value)	return 0;
    if (((sortedvals*)a)->value > ((sortedvals*)b)->value)	return 1;
}

void myqsort(int sz, Fitness* arfitness, int* arindex, sortedvals* arr)
{
    for(int i=0; i<sz; i++)
    {
        arr[i].value = arfitness[i];
        arr[i].id = i;
    }

    qsort( arr, sz, sizeof(sortedvals), compare);
    for(int i=0; i<sz; i++)
    {
        arfitness[i] = arr[i].value;
        arindex[i] = arr[i].id;
    }
}

void invAz(int N, double* Av, int iterator_sz, int* iterator, double* v_arr, double* Lj_arr, double K)
{
    for (int j = 0; j<iterator_sz; j++)
    {
        int jcur = iterator[j];
        double* v_j = &v_arr[jcur * N];
        double v_j_mult_Av = 0;
        for (int p = 0; p<N; p++)
            v_j_mult_Av += v_j[p] * Av[p];
        v_j_mult_Av = Lj_arr[jcur] * v_j_mult_Av;
        for (int p = 0; p<N; p++)
            Av[p] = K * Av[p] - v_j_mult_Av * v_j[p];
    }
}

void UpdateItrs(int itr, int nvectors, int* vec, int* iterator, int* iterator_sz, int* newidx, int* imin, int a1, int a2, double a3, int mupdateperiod)
{
    itr = itr / mupdateperiod;
    *imin = 1;
    if (itr < nvectors)
    {
        iterator[itr] = itr;
    }
    else
    {
        if (nvectors > 1) // otherwise nothing
        {
            int dmin = (vec[iterator[1]] - vec[iterator[0]]) - a1;

            for (int j = 1; j < (nvectors - 1); j++)
            {
                int dcur = (vec[iterator[j + 1]] - vec[iterator[j]]) - (a1 + a2 * pow(double(j) / double(nvectors), a3));

                if (dcur < dmin)
                {
                    dmin = dcur;
                    *imin = j + 1;
                }
            }
            if (dmin >= 0)
                *imin =  0;
            int sav = iterator[*imin];
            for (int j = *imin; j<(nvectors - 1); j++)
                iterator[j] = iterator[j + 1];
            iterator[nvectors - 1] = sav;
        }
    }

    *iterator_sz = itr + 1;
    if (*iterator_sz > nvectors)	*iterator_sz = nvectors;
    *newidx = iterator[*iterator_sz - 1];
    vec[*newidx] = itr*mupdateperiod;
}


void LMCMA(int N, int lambda, int mu, double ccov, double xmin, double xmax, int nvectors,
           double cc, double val_target, double sigma, double c_s, double target_f,
           double maxevals, int FuncId, int inseed, double* output, int printToFile, int sample_symmetry,
           int a1, int a2, double a3, int mupdateperiod,  int sample_type, double dgk_base)
{
    // memory allocation
    // m*n
    double* arx = new double[N*lambda];
    double* v_arr = new double[N*(nvectors)];
    double* pc_arr = new double[N*(nvectors)];
    // n
    double* pc = new double[N];
    double* xmean = new double[N];
    double* xold = new double[N];
    double* z = new double[N];
    double* Az = new double[N];
    double* Av = new double[N];
    // lambda, mu, nvectors
    double* weights = new double[mu];
    int* iterator = new int[nvectors];
    Fitness* arfitness = new Fitness[lambda];
    Fitness* prev_arfitness = new Fitness[lambda];
    int* arindex = new int[lambda];
    Fitness* mixed = new Fitness[2 * lambda];
    int* ranks = new int[2 * lambda];
    int* ranks_tmp = new int[2 * lambda];
    double* Nj_arr = new double[nvectors];
    double* Lj_arr = new double[nvectors];
    sortedvals* arr_tmp = new sortedvals[2 * lambda];
    int* t = new int[nvectors];
    int* vec = new int[nvectors];

    global_t gt;
    init_gt(&gt);
    m_z = inseed+2345;
    m_w = inseed+1234;

    // memory initialization
    random_init(&gt.ttime, inseed);

    double sum_weights = 0;
    for (int i = 0; i<mu; i++)
    {
        weights[i] = log(double(mu + 0.5)) - log(double(1 + i));
        sum_weights = sum_weights + weights[i];
    }
    double mueff = 0;
    for (int i = 0; i<mu; i++)
    {
        weights[i] = weights[i] / sum_weights;
        mueff = mueff + weights[i] * weights[i];
    }
    mueff = 1 / mueff;

    for (int i = 0; i<N; i++)
        pc[i] = 0;

    double K = 1 / sqrt(1 - ccov);
    double M = sqrt(1 - ccov);

    for (int i = 0; i<N; i++)
        xmean[i] = xmin + (xmax - xmin)*random_Uniform(&gt.ttime);

    double counteval = 0;
    int iterator_sz = 0;
    double s = 0;
    int stop = 0;
    int itr = 0;
    int indx = 0;

    Fitness BestF;

    FILE* pFile=NULL;
    if (printToFile == 1)
    {
        char filename[250];
        sprintf_s(filename,"LMCMA%dfunc%d_%d.txt",N,int(FuncId), inseed);
        fopen_s(&pFile, filename, "w");
    }

    if ((printToFile == 1))
    {
        BestF = MyFunc(FuncId, N, &xmean[0], &gt);
        counteval+=1;
        fprintf(pFile,"%g %g\n",counteval, fabs(BestF - target_f));
    }

    while (stop == 0)
    {
        int sign = 1;

        for (int i = 0; i<lambda; i++) // O(lambda*m*n)
        {
            if (sign == 1)
            {
                if (sample_type == 0) // Gaussian
                    for (int k = 0; k < N; k++)	// O(n)
                    {
                        z[k] = random_Gauss(&gt.ttime);
                        Az[k] = z[k];
                    }
                if (0) // Uniform()
                    for (int k = 0; k < N; k++)	// O(n)
                    {
                        if (GetUniform() < 0.5)			z[k] = 1;
                        else							z[k] = -1;
                        Az[k] = z[k];
                    }
                if (sample_type == 1) // RademacherCheap()
                    for (int k = 0; k < N; k++)	// O(n)
                    {
                        z[k] = GetRademacherCheap();
                        Az[k] = z[k];
                    }

                double dgk = dgk_base;
                if (i == 0)	dgk *= 10;
                double gg = dgk*fabs(random_Gauss(&gt.ttime));
                //if (itr == 10000)	printf("%d\n", GetRademacherCheap());
                if (gg > iterator_sz)	gg = iterator_sz;
                int k0 = 0;
                if (iterator_sz > 1)
                    k0 = iterator_sz - gg;

                for (int k = k0; k < iterator_sz; k++)	// O(m*n)
                {
                    int jcur = iterator[k];
                    double* pc_j = &pc_arr[jcur*N];
                    double* v_j = &v_arr[jcur*N];
                    double v_j_mult_z = 0;
                    for (int p = 0; p<N; p++)
                        v_j_mult_z += v_j[p] * z[p];
                    v_j_mult_z = Nj_arr[jcur] * v_j_mult_z;
                    for (int p = 0; p<N; p++)
                        Az[p] = M * Az[p] + v_j_mult_z * pc_j[p];
                }

            }
            for (int k = 0; k<N; k++)	// O(n)
                arx[i*N + k] = xmean[k] + sign*sigma*Az[k];
            if (sample_symmetry) // sample in the opposite direction, seems to work better in most cases AND decreases CPU time by 2.0
                sign = -sign;

            arfitness[i] = MyFunc(FuncId, N, &arx[i*N], &gt);
            counteval = counteval + 1;
            if (counteval == 1)	BestF = arfitness[i];
            if (arfitness[i] < BestF)	BestF = arfitness[i];
			
			if (int(counteval) % 1000 == 0)
				cout << int(counteval) << " " << fabs(BestF - target_f) << endl;


            //	if (int(counteval)%10000000 == 0)
            //		printf("%d %g\n", int(counteval), fabs(BestF- target_f));
        }

        myqsort(lambda, arfitness, arindex, arr_tmp);

        for (int i = 0; i<N; i++)
        {
            xold[i] = xmean[i];
            xmean[i] = 0;
        }

        for (int i = 0; i<mu; i++)
        {
            double* cur_x = &arx[arindex[i] * N];
            for (int j = 0; j<N; j++)
                xmean[j] += weights[i] * cur_x[j];
        }

        double pc_alpha = 1 - cc;
        double pc_beta = sqrt(cc*(2 - cc)*mueff) / sigma;
        for (int i = 0; i < N; i++)
            pc[i] = pc_alpha * pc[i] + pc_beta * (xmean[i] - xold[i]);


        int newidx = -1;
        int imin = -1;
        if (itr%mupdateperiod == 0)
        {

            UpdateItrs(itr, nvectors, vec, iterator, &iterator_sz, &newidx, &imin, a1, a2, a3, mupdateperiod);

            for (int i = 0; i < N; i++)
                pc_arr[newidx*N + i] = pc[i];

            double nv_avr = 0;
            if (imin == 1) imin = 0;
            for (int i = imin; i < iterator_sz; i++)
            {
                int indx = iterator[i];
                for (int j = 0; j < N; j++)
                    Av[j] = pc_arr[indx*N + j];
                invAz(N, Av, i, iterator, v_arr, Lj_arr, K);

                double nv = 0;
                for (int j = 0; j < N; j++)
                {
                    v_arr[indx*N + j] = Av[j];
                    nv += Av[j] * Av[j];
                }
                double sqrtccovnv = sqrt(1 + (ccov / (1 - ccov))*nv);
                nv_avr += nv;
                Nj_arr[indx] = (M / nv)*(sqrtccovnv - 1);
                Lj_arr[indx] = (K / nv)*(1 - (1 / sqrtccovnv));
            }
        }

        if (itr > 0)
        {
            for (int i = 0; i<lambda; i++)
            {
                mixed[i] = arfitness[i];
                mixed[lambda + i] = prev_arfitness[i];
            }

            myqsort(2 * lambda, mixed, ranks, arr_tmp);

            double meanprev = 0;
            double meancur = 0;
            for (int i = 0; i<2 * lambda; i++)
                ranks_tmp[i] = ranks[i];
            for (int i = 0; i<2 * lambda; i++)
                ranks[ranks_tmp[i]] = i;
            for (int i = 0; i<lambda; i++)
            {
                meancur = meancur + ranks[i];
                meanprev = meanprev + ranks[lambda + i];
            }
            meancur = meancur / lambda;
            meanprev = meanprev / lambda;
            double diffv = (meanprev - meancur) / lambda;
            double z1 = diffv - val_target;
            //if (itr % 2 == 0)
            //{
            s = (1 - c_s)*s + c_s*z1;
            double d_s = 1;//2.0*(N-1.0)/N;
            sigma = sigma * exp(s / d_s);
            //}
        }


        for (int i = 0; i<lambda; i++)
            prev_arfitness[i] = arfitness[i];

        if (arfitness[0] < target_f)
            stop = 1;
        if (counteval >= maxevals)
            stop = 1;
        itr = itr + 1;

        if (sigma < 1e-15)
            stop = 1;
		
        if ((printToFile == 1) && (pFile) && ((itr%100 == 0) || (counteval < 1000)))
        {
            fprintf(pFile,"%g %Le\n",counteval, fabs(BestF- target_f));
            fflush (pFile);
        }

    }

    if ((printToFile == 1) && (pFile))
        fprintf(pFile,"%g %Le\n",counteval, fabs(BestF - target_f));

    output[0] = counteval;
    output[1] = fabs(BestF - target_f);

    if (printToFile == 1)
        fclose(pFile);

    random_exit(&gt.ttime);
    delete[] arr_tmp;
    delete[] weights;
    delete[] pc;
    delete[] xmean;
    delete[] xold;
    delete[] z;
    delete[] Az;
    delete[] iterator;
    delete[] v_arr;
    delete[] pc_arr;
    delete[] arx;
    delete[] arfitness;
    delete[] prev_arfitness;
    delete[] arindex;
    delete[] Av;
    delete[] t;
    delete[] vec;
    delete[] mixed;
    delete[] ranks;
    delete[]ranks_tmp;
    delete[] Nj_arr;
    delete[] Lj_arr;
    free_gt(&gt);
}

void sepCMA(int N, int lambda, int mu,double xmin, double xmax, double sigma, Fitness target_f,
            double maxevals, int FuncId, int inseed, double* output, int printToFile)
{
    // memory allocation
    // m*n
    double* arx = new double[N*lambda];
    double* arz = new double[N*lambda];
    // n
    double* diagD = new double[N];
    double* diagC = new double[N];
    double* pc = new double[N];
    double* ps = new double[N];
    double* xmean = new double[N];
    double* xold = new double[N];
    // lambda, mu, nvectors
    double* weights = new double[mu];
    Fitness* arfitness = new Fitness[lambda];
    int* arindex = new int[lambda];
    sortedvals* arr_tmp = new sortedvals[2 * lambda];

    global_t gt;
    init_gt(&gt);

    m_z = inseed+2345;
    m_w = inseed+1234;

    // memory initialization
    random_init(&gt.ttime, inseed);

    for (int i = 0; i<N; i++)
    {
        diagD[i] = 1;
        diagC[i] = 1;
        ps[i] = 0;
        pc[i] = 0;
    }

    double sum_weights = 0;
    for (int i = 0; i<mu; i++)
    {
        weights[i] = log(double(mu + 0.5)) - log(double(1 + i));
        sum_weights = sum_weights + weights[i];
    }
    double mueff = 0;
    for (int i = 0; i<mu; i++)
    {
        weights[i] = weights[i] / sum_weights;
        mueff = mueff + weights[i] * weights[i];
    }
    mueff = 1 / mueff;

    for (int i = 0; i<N; i++)
        xmean[i] = xmin + (xmax - xmin)*random_Uniform(&gt.ttime);

    double c1 = 2.0 / (pow((N + 1.3), 2.0) + mueff);
    double cmu = minv(1.0 - c1, 2.0 * (mueff - 2.0 + 1.0 / mueff) / (pow(N + 2.0, 2.0) + mueff));
    double ccov1_sep = minv(1, c1 * (N + 1.5) / 3.0);
    double ccovmu_sep = minv(1.0 - ccov1_sep, cmu * (N + 1.5) / 3);
    double chiN = sqrt(double(N))*(1.0 - 1.0 / (4.0*N) + 1 / (21.0*N*N));
    double cs = (mueff + 2.0) / (N + mueff + 5.0);
    double damps = 1.0 + 2 * maxv(0, sqrt((mueff - 1.0) / (N + 1.0)) - 1.0) + cs;
    double ccum = (4.0 + mueff/double(N)) / (N+4.0 + 2*mueff/double(N));


    FILE* pFile = NULL;
    if (printToFile == 1)
    {
        char filename[250];
        sprintf_s(filename,"SEPCMA%dfunc%d_%d.txt",N,int(FuncId),inseed);
        fopen_s(&pFile, filename, "w");
    }

    Fitness BestF;
    double counteval = 0;
    int stop = 0;
    int itr = 0;

    if ((printToFile == 1))
    {
        BestF = MyFunc(FuncId, N, &xmean[0], &gt);
        counteval+=1;
        fprintf(pFile,"%g %Le\n",counteval, fabs(BestF - target_f));
    }

    while (stop == 0)
    {
        for (int i = 0; i<lambda; i++) // O(lambda*m*n)
        {
            for (int k = 0; k<N; k++)	// O(n)
            {
                arz[i*N + k] = random_Gauss(&gt.ttime);
                arx[i*N + k] = xmean[k] + sigma * diagD[k] * arz[i*N + k];
            }
            arfitness[i] = MyFunc(FuncId, N, &arx[i*N], &gt);
            counteval = counteval + 1;
            if (counteval == 1)	BestF = arfitness[i];
            if (arfitness[i] < BestF)	BestF = arfitness[i];
        }

        myqsort(lambda, arfitness, arindex, arr_tmp);

        for (int i = 0; i<N; i++)
        {
            xold[i] = xmean[i];
            xmean[i] = 0;
        }

        for (int i = 0; i<mu; i++)
        {
            double* cur_x = &arx[arindex[i] * N];
            for (int j = 0; j<N; j++)
                xmean[j] += weights[i] * cur_x[j];
        }

        double norm_ps = 0;
        for (int i = 0; i<N; i++)
        {
            ps[i] = (1 - cs) * ps[i] + sqrt(cs*(2 - cs)*mueff) * (1. / diagD[i]) * (xmean[i] - xold[i]) / sigma;
            norm_ps += ps[i] * ps[i];
        }
        norm_ps = sqrt(norm_ps);

        for (int i = 0; i<N; i++)
            pc[i] = (1 - ccum)*pc[i] + sqrt(ccum*(2 - ccum)*mueff)*(xmean[i] - xold[i]) / sigma;

        for (int i = 0; i<N; i++)
        {
            double val = 0;
            for (int j = 0; j< mu; j++)
                val += weights[j] * arz[arindex[j] * N + i] * arz[arindex[j] * N + i];
            diagC[i] = (1 - ccov1_sep - ccovmu_sep) * diagC[i] + ccov1_sep * pc[i] * pc[i] + ccovmu_sep * (diagC[i] * val);
        }
        sigma = sigma * exp((cs / damps)*(norm_ps / chiN - 1));

        for (int i = 0; i<N; i++)
            diagD[i] = sqrt(diagC[i]);

        if (arfitness[0] < target_f)
            stop = 1;
        if (counteval >= maxevals)
            stop = 1;
        itr = itr + 1;
        if ((printToFile == 1) && (pFile) && ((itr%100 == 0) || (counteval < 1000)))
        {
            fprintf(pFile,"%d %lf\n",counteval, fabs(BestF - target_f));
            fflush (pFile);
        }

		//if (int(counteval) % 1000 == 0)
			cout << int(counteval) << " " << fabs(BestF - target_f) << endl;
    }


    if ((printToFile == 1) && (pFile))
        fprintf(pFile,"%g %Le\n",counteval, fabs(BestF - target_f));

    output[0] = counteval;
    output[1] = fabs(BestF - target_f);

    if (printToFile == 1)
        fclose(pFile);

    random_exit(&gt.ttime);
    delete[] arr_tmp;
    delete[] weights;
    delete[] pc;
    delete[] xmean;
    delete[] xold;
    delete[] arz;
    delete[] arx;
    delete[] arfitness;
    delete[] arindex;
    delete[] diagD;
    delete[] diagC;
    delete[] ps;
    free_gt(&gt);
}






void CholeskyUpdate(double* dZ, double* Ainv, double* A,
                    double* tmp_vec, double* tmp_vec2,
                    double alpha, double gamma, int N)
{
    //w = Ainv * pc;
    matrix_mult_vector(tmp_vec, Ainv, dZ, N);
    //prodw = norm(w)^2;
    double prodw = vector_prod(tmp_vec, tmp_vec, N);

    // pc * w'
    double k1 = sqrt(alpha);
    double k2 = (sqrt(alpha) / prodw) * (sqrt(1 + (gamma / alpha)*prodw) - 1);
    //A = sqrt(alpha) * A + (sqrt(alpha)/prodw) * (sqrt(1 + (gamma1/alpha)*prodw) -1) * dZ * w';
    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            A[i*N + j] = k1 * A[i*N + j] + k2 * dZ[i] * tmp_vec[j];
    //	A[i*N + j] = k1 * A[i*N + j] + k2 * tmp_mat[i*N + j];

    //Ainv = (1/sqrt(alpha))*Ainv - (1/(sqrt(alpha)*prodw)) * (1 - 1 / (sqrt(1 + (gamma1/alpha)*prodw))) * w * (w' * Ainv );
    k1 = 1 / sqrt(alpha);
    k2 = (1 / (sqrt(alpha)*prodw)) * (1 - 1 / (sqrt(1 + (gamma / alpha)*prodw)));
    vector_mult_matrix(tmp_vec2, tmp_vec, Ainv, N);	// w' * Ainv
    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            Ainv[i*N + j] = k1 * Ainv[i*N + j] - k2 * tmp_vec[i] * tmp_vec2[j];
    //	Ainv[i*N + j] = k1 * Ainv[i*N + j] - k2 * tmp_mat[i*N + j];
}




void runCMAmulambdaCholesky(int N, int lambda, int mu, double sigma, double xmin, double xmax, double maxevals, double target_f, int  FuncId, int verbose, int inseed, double* output, int printToFile)
{
    // memory allocation
    // n*n
    double* A = new double[N*N];
    double* Ainv = new double[N*N];
    // m*n
    double* arx = new double[N*lambda];
    double* arz = new double[N*lambda];
    // n
    double* pc = new double[N];
    double* ps = new double[N];
    double* xmean = new double[N];
    double* zmean = new double[N];
    double* xold = new double[N];
    double* tmp_vec = new double[N];
    double* tmp_vec2 = new double[N];
    // lambda, mu, nvectors
    double* weights = new double[mu];
    Fitness* arfitness = new Fitness[lambda];
    int* arindex = new int[lambda];
    sortedvals* arr_tmp = new sortedvals[2 * lambda];

    global_t gt;
    init_gt(&gt);

    m_z = inseed+2345;
    m_w = inseed+1234;

    // memory initialization
    random_init(&gt.ttime, inseed);

    for (int i = 0; i<N; i++)
    {
        ps[i] = 0;
        pc[i] = 0;
        for (int j = 0; j<N; j++)
            if (i == j)	{
                A[i*N + j] = 1;
                Ainv[i*N + j] = 1;
            }
            else		{
                A[i*N + j] = 0;
                Ainv[i*N + j] = 0;
            }
    }

    double sum_w = 0;
    for (int i = 0; i<mu; i++)
        sum_w += log(i + 1.0);

    double mueff = 0;
    for (int i = 0; i<mu; i++)
    {
        weights[i] = (log(double(mu + 1.0)) - log(double(1.0 + i))) / (mu*log(mu + 1.0) - sum_w);
        mueff = mueff + weights[i] * weights[i];
    }
    mueff = 1 / mueff;

    for (int i = 0; i<N; i++)
        xmean[i] = xmin + (xmax - xmin)*random_Uniform(&gt.ttime);

    double cc = 4.0 / (N + 4.0);
    double ccov = 2.0 / (pow(N + sqrt(2.0), 2.0));
    double chiN = sqrt(double(N))*(1.0 - (1.0 / (4.0*double(N))) + (1.0 / (21.0*double(N*N))));
    double cs = sqrt(mueff) / (sqrt(double(N)) + sqrt(mueff));
    double damps = 1.0 + 2.0*maxv(0, sqrt((mueff - 1.0) / (N + 1.0)) - 1.0) + cs;


    double counteval = 0;
    int stop = 0;
    int itr = 0;

    double BestF;

    FILE* pFile=NULL;
    if (printToFile == 1)
    {
        char filename[250];
        sprintf_s(filename,"CHOLCMA%dfunc%d_%d.txt",N,int(FuncId),inseed);
        fopen_s(&pFile, filename, "w");
    }

    if ((printToFile == 1))
    {
        BestF = MyFunc(FuncId, N, &xmean[0], &gt);
        counteval+=1;
        fprintf(pFile,"%g %le\n",counteval, fabs(BestF - target_f));
    }

    while (stop == 0)
    {
        for (int i = 0; i<lambda; i++) // O(lambda*m*n)
        {
            for (int k = 0; k<N; k++)	// O(n)
            {
                arz[i*N + k] = random_Gauss(&gt.ttime);
            }
            for (int k = 0; k<N; k++)	// O(n^2)
            {
                double Az = 0;
                double* xcur = &arz[i*N];
                double* Acur = &A[k*N];
                for (int p = 0; p<N; p++)
                    Az += Acur[p] * xcur[p];
                arx[i*N + k] = xmean[k] + sigma * Az;
            }
            arfitness[i] = MyFunc(FuncId, N, &arx[i*N], &gt);
            counteval = counteval + 1;
            if (counteval == 1)	BestF = arfitness[i];
            if (arfitness[i] < BestF)	BestF = arfitness[i];
        }

        myqsort(lambda, arfitness, arindex, arr_tmp);

        for (int i = 0; i<N; i++)
        {
            xold[i] = xmean[i];
            xmean[i] = 0;
            zmean[i] = 0;
        }

        for (int i = 0; i<mu; i++)
        {
            double* cur_x = &arx[arindex[i] * N];
            double* cur_z = &arz[arindex[i] * N];
            for (int j = 0; j<N; j++)
            {
                xmean[j] += weights[i] * cur_x[j];
                zmean[j] += weights[i] * cur_z[j];
            }
        }

        double norm_ps = 0;
        for (int i = 0; i<N; i++)
        {
            ps[i] = (1 - cs) * ps[i] + sqrt(cs*(2 - cs)*mueff) * zmean[i];
            norm_ps += ps[i] * ps[i];
        }
        norm_ps = sqrt(norm_ps);

        for (int i = 0; i<N; i++)
            pc[i] = (1 - cc)*pc[i] + sqrt(cc*(2 - cc)*mueff)*(xmean[i] - xold[i]) / sigma;

        double alpha = (1 - ccov);
        double gamma = ccov;
        CholeskyUpdate(pc, Ainv, A, tmp_vec, tmp_vec2, alpha, gamma, N);

        sigma = sigma * exp((cs / damps)*(norm_ps / chiN - 1));

        if (arfitness[0] < target_f)
            stop = 1;
        if (counteval >= maxevals)
            stop = 1;
        itr = itr + 1;

        if ((printToFile == 1) && (pFile) && ((itr%100 == 0) || (counteval < 1000)))
        {
            fprintf(pFile,"%g %le\n",counteval, fabs(BestF - target_f));
            fflush (pFile);
        }
    }

    if ((printToFile == 1) && (pFile))
        fprintf(pFile,"%g %le\n",counteval, fabs(BestF - target_f));

    output[0] = counteval;
    output[1] = fabs(BestF - target_f);

    if (printToFile == 1)
        fclose(pFile);

    output[0] = counteval;
    output[1] = fabs(BestF - target_f);

    random_exit(&gt.ttime);
    delete[] arr_tmp;
    delete[] weights;
    delete[] pc;
    delete[] xmean;
    delete[] xold;
    delete[] arz;
    delete[] arx;
    delete[] arfitness;
    delete[] arindex;
    delete[] tmp_vec;
    delete[] tmp_vec2;
    delete[] ps;
    delete[] A;
    delete[] Ainv;
    delete[] zmean;
    free_gt(&gt);
}


void pveq_qveq(int mu, int N, double* y_vn, double* vn, double* y, double* pvec, double* qvec, double* weights, double norm_v2)
{
    for (int i = 0; i < mu; i++)
    {
        y_vn[i] = 0;
        for (int k = 0; k < N; k++)
            y_vn[i] += vn[k] * y[i*N + k];
    }
    for (int i = 0; i < N; i++)
    {
        pvec[i] = 0;
        for (int k = 0; k < mu; k++)
            pvec[i] += (y[k*N + i] * y[k*N + i] - (norm_v2 / (1.0 + norm_v2))*(vn[i] * y_vn[k] * y[k*N + i]) - 1.0)*weights[k];

        qvec[i] = 0;
        for (int k = 0; k < mu; k++)
            qvec[i] += (y[k*N + i] * y_vn[k] - (vn[i] / 2.0)*(y_vn[k] * y_vn[k] + (1.0 + norm_v2)))*weights[k];
    }
}

void VDCMA(int N, double xmin, double xmax, double sigma, double ftarget,
           double maxeval, int  FuncId, int inseed, double* output, int printToFile)
{

    global_t gt;
    init_gt(&gt);
    m_z = inseed+2345;
    m_w = inseed+1234;

    // memory initialization
    random_init(&gt.ttime, inseed);

    // Constants
    double chiN = sqrt(double(N))*(1.0 - 1.0 / (4.0*N) + 1 / (21.0*N*N));

    // Dynamic parameters
    double* xmean = new double[N];
    double* v = new double[N];
    double* D = new double[N];
    double* pc = new double[N];
    double* ps = new double[N];
    double* vn = new double[N];
    double* vnn = new double[N];
    double* ymean = new double[N];
    double* zmean = new double[N];
    double* invAvnn = new double[N];
    double* pcD = new double[N];

    for (int i = 0; i < N; i++)
        xmean[i] = xmin + (xmax - xmin)*random_Uniform(&gt.ttime);
    //	xmean[i] = 2;// xmin + (xmax - xmin)*random_Uniform(&gt.ttime);

    for (int i = 0; i < N; i++)
    {
        //	v[i] = 0.01*i;//
        v[i] = random_Gauss(&gt.ttime) / (sqrt(double(N)));
        D[i] = 1;
        pc[i] = 0;
        ps[i] = 0;
    }

    // Static parameters
    int lambda = 4 + floor(3 * log(double(N)));
    int mu = lambda / 2;
    double mu_double = double(lambda) / 2.0;

    double* arx = new double[N*lambda];
    double* arz = new double[N];
    double* ary = new double[N*lambda];
    double* arys = new double[N*mu];
    Fitness* arfitness = new Fitness[lambda];
    int* arindex = new int[lambda];
    sortedvals* arr_tmp = new sortedvals[2 * lambda];
    double* y_vn = new double[mu];
    double* pvec = new double[N];
    double* qvec = new double[N];
    double* rvec = new double[N];
    double* svec = new double[N];
    double* A = new double[N];
    double* ngv = new double[N];
    double* ngd = new double[N];
    double* pone = new double[N];
    double* qone = new double[N];

    double* weights = new double[mu];
    double sum_weights = 0;
    for (int i = 0; i<mu; i++)
    {
        weights[i] = log(double(mu_double + 0.5)) - log(double(1 + i));
        sum_weights = sum_weights + weights[i];
    }
    double mueff = 0;
    for (int i = 0; i<mu; i++)
    {
        weights[i] = weights[i] / sum_weights;
        mueff = mueff + weights[i] * weights[i];
    }
    mueff = 1 / mueff;

    double c1_old = 2.0 / (pow(double(N) + 1.3, 2.0) + mueff);
    double cmu_old = minv(1 - c1_old, 2 * (mueff - 2 + 1.0 / mueff) / (pow(N + 2.0,2.0) + mueff));
    double cm = 1;
    double cc = (4.0 + mueff / double(N)) / (N + 4.0 + 2.0 * mueff / double(N));
    double cs = sqrt(mueff) / 2 / (sqrt(double(N)) + sqrt(mueff));
    double ds = 1 + 2 * maxv(0, sqrt((mueff - 1) / (N + 1.0)) - 1) + cs;
    double c1 = (double(N) - 5.0) / 6.0 * c1_old;
    double cmu = minv(1 - c1, (double(N) - 5.0) / 6.0 * cmu_old);

    //Termination criteria
    double minsigma = 1e-15;
    double maxsigma = 1e15;
    //double ftarget = 1e-10;
    //double maxeval = 1e5*N;
    double maxiter = 1e+30;

    double counteval = 0; // number of function evaluations
    double iIter = 0;
    double BestF = 1e+30;

    FILE* pFile=NULL;
    if (printToFile == 1)
    {
        char filename[250];
        sprintf_s(filename,"VDCMA%dfunc%d_%d.txt",N,int(FuncId),inseed);
        fopen_s(&pFile, filename, "w");
    }

    if ((printToFile == 1))
    {
        BestF = MyFunc(FuncId, N, &xmean[0], &gt);
        counteval+=1;
        fprintf(pFile,"%g %g\n",counteval,BestF);
    }

    while (iIter < maxiter)
    {
        iIter = iIter + 1;
        // Compute constants : O(N) multiplication + one sqrt
        //norm_v2 = v'*v;
        double norm_v2 = 0;
        for (int i = 0; i < N; i++)
            norm_v2 += v[i] * v[i];

        double norm_v = sqrt(norm_v2);	// norm_v = sqrt(norm_v2);
        //  vn = v / norm_v;
        for (int i = 0; i < N; i++)
            vn[i] = v[i] / norm_v;

        // vnn = vn .* vn;
        for (int i = 0; i < N; i++)
            vnn[i] = vn[i] * vn[i];


        // Sampling: O(N*lambda) multiplication + one sqrt
        for (int i = 0; i < lambda; i++)
        {
            double vnarz = 0;
            for (int k = 0; k < N; k++)
            {
                //		arz[k] = (-10 + int(iIter)%20)*0.01*(i + 1) / double((k + 1));//
                arz[k] = random_Gauss(&gt.ttime);
                vnarz += vn[k] * arz[k];
            }
            double prodf = (sqrt(1.0 + norm_v2) - 1) * vnarz;
            for (int k = 0; k < N; k++)
                ary[i*N + k] = arz[k] + prodf * vn[k];

            for (int k = 0; k < N; k++)
                arx[i*N + k] = xmean[k] + sigma * D[k] * ary[i*N + k];

            arfitness[i] = MyFunc(FuncId, N, &arx[i*N], &gt);
            counteval = counteval + 1;
            if (counteval == 1)	BestF = arfitness[i];
            if (arfitness[i] < BestF)	BestF = arfitness[i];

            if (int(counteval) % 1000000 == 0)
                printf("%g %g\n", counteval, BestF);
        }

        myqsort(lambda, arfitness, arindex, arr_tmp);
        for (int i = 0; i < mu; i++)
        {
            double* src = &ary[arindex[i] * N];
            double* dst = &arys[i * N];
            for (int k = 0; k < N; k++)
                dst[k] = src[k];
        }

        double ymeanvn = 0;
        for (int k = 0; k < N; k++)
        {
            ymean[k] = 0;
            for (int i = 0; i < mu; i++)
                ymean[k] += weights[i] * arys[i * N + k];
            ymeanvn += ymean[k] * vn[k];
        }

        double prodf = 1.0 / sqrt(1.0 + norm_v2) - 1;
        prodf = prodf * ymeanvn;
        for (int i = 0; i < N; i++)
            zmean[i] = ymean[i] + prodf*vn[i];

        double k1 = 1 - cs;
        double k2 = sqrt(cs * (2 - cs) * mueff);
        double psps = 0;
        for (int i = 0; i < N; i++)
        {
            ps[i] = k1 * ps[i] + k2 * zmean[i];
            psps += ps[i] * ps[i];
        }

        int hsig = 0;
        if ((psps / double(N)) < (2.0 + 4.0 / (double(N) + 1.0))*(1 - pow(1 - cs, 2.0 * iIter)))
            hsig = 1;

        k1 = 1 - cc;
        k2 = hsig * sqrt(cc*(2 - cc)*mueff);
        for (int i = 0; i < N; i++)
            pc[i] = k1*pc[i] + k2*D[i]*ymean[i];

        // Modified Natural Gradients for v and D
        // [alpha A b invAvnn] = alpha_avec_bsca_invavnn(norm_v2, vnn);
        // Compute alpha, A, b and A^(-1)*vnn
        double gamma = 1.0 / sqrt(1.0 + norm_v2);
        double maxvnn = 0;
        for (int i = 0; i < N; i++)
            if (vnn[i] > maxvnn)	maxvnn = vnn[i];
        double alpha = sqrt(norm_v2*norm_v2 + ((1 + norm_v2) / maxvnn)*(2.0 - gamma)) / (2.0 + norm_v2);
        double beta;
        if (alpha < 1)
            beta = (4.0 - (2.0 - gamma) / maxvnn) / pow((1.0 + 2.0 / norm_v2),2.0);
        else				{
            alpha = 1;
            beta = 0;
        }
        double b = 2 * pow(alpha,2.0) - beta;
        for (int i = 0; i < N; i++)
            A[i] = 2.0 - (b + 2.0 * alpha * alpha) * vnn[i];
        for (int i = 0; i < N; i++)
            invAvnn[i] = vnn[i] / A[i];

        // [pvec qvec] = pvec_qvec(norm_v2, vn, arys, weights);
        // Compute pvec and qvec
        pveq_qveq(mu, N, y_vn, vn, arys, pvec, qvec,weights, norm_v2);

        if (hsig == 1)
        {
            for (int i = 0; i < N; i++)
                pcD[i] = pc[i] / D[i];
            double wone = 1;
            pveq_qveq(1, N, y_vn, vn, pcD, pone, qone, &wone, norm_v2);
            for (int i = 0; i < N; i++)
            {
                pvec[i] = cmu * pvec[i] + c1 * pone[i];
                qvec[i] = cmu * qvec[i] + c1 * qone[i];
            }
        }
        else
        {
            for (int i = 0; i < N; i++)
            {
                pvec[i] = cmu * pvec[i];
                qvec[i] = cmu * qvec[i];
            }
        }

        double updatefactor;
        if (cmu + c1*hsig > 0)
        {
            // [ngv ngd] = ngv_ngd(D, norm_v, norm_v2, vn, vnn, alpha, A, b, invAvnn, pvec, qvec);
            // Compute the modified natural gradient for v and D
            double vnqvec = 0;
            for (int i = 0; i < N; i++)
                vnqvec += vn[i] * qvec[i];
            for (int i = 0; i < N; i++)
                rvec[i] = pvec[i] - (alpha / (1.0 + norm_v2)) * ((2.0 + norm_v2)*(vn[i] * qvec[i]) - norm_v2 * vnqvec * vnn[i]);

            double invAvnnrvec = 0;
            double invAvnnvnn = 0;
            for (int i = 0; i < N; i++)
            {
                invAvnnrvec += invAvnn[i] * rvec[i];
                invAvnnvnn += invAvnn[i] * vnn[i];
            }
            for (int i = 0; i < N; i++)
                svec[i] = rvec[i] / A[i] - (b * invAvnnrvec / (1.0 + b*invAvnnvnn))*invAvnn[i];
            double svecvnn = 0;
            for (int i = 0; i < N; i++)
                svecvnn += svec[i] * vnn[i];
            double normngv = 0;
            for (int i = 0; i < N; i++)
            {
                ngv[i] = (qvec[i] - alpha * ((2.0 + norm_v2) * vn[i] * svec[i] - svecvnn * vn[i])) / norm_v;
                ngd[i] = D[i] * svec[i];
                normngv += ngv[i] * ngv[i];
            }
            normngv = sqrt(normngv);


            updatefactor = 1;
            if (0.7 * norm_v / normngv < 1)
                updatefactor = 0.7 * norm_v / normngv;
            double minabsDabsngd = 1e+30;
            for (int i = 0; i < N; i++)
                if (fabs(D[i]) / fabs(ngd[i]) < minabsDabsngd)
                    minabsDabsngd = fabs(D[i]) / fabs(ngd[i]);
            if (0.7 * minabsDabsngd < 1)
                updatefactor = minv(updatefactor, 0.7 *  minabsDabsngd);
        }
        else
        {
            for (int i = 0; i < N; i++)
            {
                ngv[i] = 0;
                ngd[i] = 0;
            }
            updatefactor = 0;
        }

        // Update
        double norm_ps = 0;
        for (int i = 0; i < N; i++)
        {
            xmean[i] += +cm * (sigma * D[i] * ymean[i]);
            norm_ps += ps[i] * ps[i];
        }
        norm_ps = sqrt(norm_ps);

        sigma = sigma * exp(cs / ds * (norm_ps / chiN - 1.0));
        for (int i = 0; i < N; i++)
        {
            v[i] += updatefactor * ngv[i];
            D[i] += updatefactor * ngd[i];
        }
        if (sigma < 1e-15)
        {
            break;
        }
        if (BestF < ftarget)
        {
            break;
        }
        if (counteval > maxeval)
        {
            break;
        }
        if ((printToFile == 1) && (pFile) && ((int(iIter)%100 == 0) || (counteval < 1000)))
        {
            fprintf(pFile,"%g %g\n",counteval,BestF);
            fflush (pFile);
        }
    }

    if ((printToFile == 1) && (pFile))
        fprintf(pFile,"%g %g\n",counteval,BestF);

    output[0] = counteval;
    output[1] = BestF;

    if (printToFile == 1)
        fclose(pFile);

    random_exit(&gt.ttime);
    delete[] xmean;
    delete[] v;
    delete[] D;
    delete[] pc;
    delete[] ps;
    delete[] weights;
    delete[] vn;
    delete[] vnn;
    delete[] arx;
    delete[] arz;
    delete[] ary;
    delete[] arfitness;
    delete[] arindex;
    delete[] arr_tmp;
    delete[] arys;
    delete[] ymean;
    delete[] zmean;
    delete[] invAvnn;
    delete[] y_vn;
    delete[] pvec;
    delete[] qvec;
    delete[] pcD;
    delete[] rvec;
    delete[] svec;
    delete[] A;
    delete[] ngv;
    delete[] ngd;
}


double time_costOfSimpleOperation(int OperationType, double secondsToRun, int N, int verbose, int FuncId)
{
    //global_t tt;
    //global_alloc(&tt,1,N,FuncId,1,0,0,0,NULL,0);

    global_t gt;
    gt.func_tempdata = NULL;
    gt.x_tempdata = NULL;
    gt.rotmatrix = NULL;
    int inseed = 1;
    random_init(&gt.ttime, inseed);

    m_z = inseed + 2345;
    m_w = inseed + 1234;

    double* M = NULL;
    if (OperationType == 1)
        M = new double[N*N];

    bool run = true;
    double totalEvals = 0;
    double delta = 0;
    double N2 = N*N;
    double* V1 = new double[N];
    double* V2 = new double[N];

    for (int j = 0; j<N; j++)
    {
        V1[j] = random_Uniform(&gt.ttime);
        V2[j] = random_Uniform(&gt.ttime);
    }

    if (OperationType == 3) // initialize matrices, etc
    {
        double val = MyFunc((int )FuncId, N, V1, &gt);
    }

    double val = 0;
    time_tic(&gt);
    double factor = N;
    if ((OperationType == 1) || (OperationType == 1))	factor = 1;
    int ij = 0;
    while (run == true)
    {

        double nEvals = (1e+8 / double(N2));
        //	if ((OperationType == 1) || (OperationType == 2))	nEvals = double(1e+8 / (double(N*N)));
        if (nEvals < 1)	nEvals = 1;
        for (double iEval = 0; iEval<nEvals; iEval = iEval + 1)
        {
            if (OperationType == 0)
            {
                val = random_Uniform(&gt.ttime);
                for (int i = 0; i < N; i++)
                    for (int j = 0; j<N; j++)
                        V2[j] = V1[j] * val;
            }
            if (OperationType == 1)
            {
                for (int j = 0; j<N; j++)
                    V1[j] = random_Uniform(&gt.ttime);
                double* pM;
                for (int j = 0; j<N; j++)
                {
                    val = 0.0;
                    pM = &M[j*N];
                    for (int k = 0; k<N; k++)
                        val += pM[k] * V1[k];
                    V2[j] = val;
                }
            }
            if (OperationType == 2)
            {
                for (int k = 0; k<N2; k++)
                    val = random_Gauss(&gt.ttime);
            }
            if (OperationType == 3)
            {
                for (int k = 0; k<N; k++)
                    val = MyFunc((int )FuncId, N, V1, &gt);
            }
            if (OperationType == 4)
            {
                for (int k = 0; k<N2; k++)
                    val = GetRademacherCheap();
            }
        }
        totalEvals = totalEvals + nEvals *factor;
        delta = time_toc(&gt);
        if (delta > secondsToRun)
            run = false;
    }

    double timePerEval = delta / totalEvals;
    if (verbose > 0)
        printf("TotalTime: %e \t NEvals: %e \t PerEval: %e\n", delta, totalEvals, timePerEval);

    if (OperationType == 1)
        delete[] M;
    delete[] V1;
    delete[] V2;
    delete[] gt.func_tempdata;
    delete[] gt.x_tempdata;
    delete[] gt.rotmatrix;
    return timePerEval;
}


int main(int argc, char **argv)
{
    //number of runs
    int num_runs = 25;

    //available number of fitness evaluations

    int FuncId = 1;

    if (argc > 1)
    {
        FuncId = atoi(argv[1]);
    }
    else
    {
        cout << "Please provide the function number as a command line parameter" << endl;
        exit(0);
    }

    Fitness target_f = functionsOptima[FuncId - 1];
    double boundMin = -functionsDomainBounds[FuncId - 1];
    double boundMax = functionsDomainBounds[FuncId - 1];

    //input
    int N = 1000;			//	problem dimension
	double maxevals = N * 5000;		// maximum number of function evaluations allowed, e.g., 1e+6

    int lambda = 4 + floor(3 * log(N));	// 	population size, e.g., 4+floor(3*log(N));
    int mu = floor(lambda / 2);		// 	number of parents, e.g., floor(lambda/2);
    double ccov = 1.0 / (10 * log(N + 1));  // 	learning rate for covariance matrix, e.g., 1/(10*log(N+1))
    double xmin = boundMin;//	x parameters lower bound
    double xmax = boundMax;//	x parameters upper bound
    int nvectors = 4 + floor(3 * log(N));	//	number of stored direction vectors, e.g., nvectors = 4+floor(3*log(N))
    int a1 = N;		// 	parameters to control target distances between stored direction vectors
    int a2 = 0;
    double a3 = 1.0;
    double cc = 1.0 / nvectors;	// learning rate for mean vector's evolution path, e.g., cc = 1/nvectors
    double val_target = 0.25;	// target success rate for new population, e.g., 0.25
    double sigma = 0.1;	// initial step-size, e.g., 0.5
    double c_s = 0.3;	//	decay factor for step-size adaptation, e.g., 0.3 
    int sample_symmetry = 0;	// 1 or 0, to sample symmetrical solutions to save 50% time and sometimes evaluations
    int sample_type = 0;	 	// 0 - Gaussian, 1 - Rademacher
    int mupdateperiod = floor(log(N));
    double dgk_base = 4.0;
    int inseed = 1;		// initial seed for random number generator, e.g., 1
    int algorithmType = 0;
    int printToFile = 1; // 1 or 0
    // output
	cout << "Function f" << FuncId << endl;
	cout << "problem dimension=" << N << endl;
	cout << "maxevals=" << maxevals << endl;

    double output[2];

	cout << scientific << setprecision(8);

	Fitness *bsf_fitness_array = (Fitness*)malloc(sizeof(Fitness) * num_runs);
	Fitness mean_bsf_fitness = 0;
	Fitness std_bsf_fitness = 0;

	for (int j = 0; j < num_runs; j++)
	{
		inseed = j * 10;
		if (algorithmType == 0)
		{
			LMCMA(N, lambda, mu, ccov, xmin, xmax, nvectors, cc, val_target, sigma, c_s, target_f,
				maxevals, FuncId, inseed, output, printToFile, sample_symmetry, a1, a2, a3, mupdateperiod, sample_type, dgk_base);
		}
		if (algorithmType == 1)
		{
			VDCMA(N, xmin, xmax, sigma, target_f, maxevals, FuncId, inseed, output, printToFile);
		}
		if (algorithmType == 2)
		{
			sepCMA(N, lambda, mu, xmin, xmax, sigma, target_f, maxevals, FuncId, inseed, output, printToFile);
		}
		if (algorithmType == 3)
		{
			runCMAmulambdaCholesky(N, lambda, mu, sigma, xmin, xmax, maxevals, target_f, FuncId, 0, inseed, output, printToFile);
		}
		if (algorithmType < 0)
		{
			int OperationType = 0;
			double timetorun = 10;
			if (algorithmType == -1) OperationType = 0;	// check Vector x Scalar Multiplication Cost
			if (algorithmType == -2) OperationType = 1;	// check Matrix x Vector Multiplication Cost
			if (algorithmType == -3) OperationType = 2;	// check Random Gaussian Vector Cost
			if (algorithmType == -4) OperationType = 3;	// check Function Evaluation Cost
			if (algorithmType == -5) OperationType = 4;	// check Rademacher Vector Cost
			output[1] = time_costOfSimpleOperation(OperationType, timetorun, N, 0, FuncId);
			output[0] = maxevals;
		}

		cout << j + 1 << "th run, " << "error value = " << output[1] << endl;
		bsf_fitness_array[j] = output[1];
	}

	for (int j = 0; j < num_runs; j++) mean_bsf_fitness += bsf_fitness_array[j];
	mean_bsf_fitness /= num_runs;

	for (int j = 0; j < num_runs; j++) std_bsf_fitness += pow((mean_bsf_fitness - bsf_fitness_array[j]), 2.0);
	std_bsf_fitness /= num_runs;
	std_bsf_fitness = sqrt(std_bsf_fitness);

	cout << "\nmean = " << mean_bsf_fitness << ", std = " << std_bsf_fitness << endl;
	free(bsf_fitness_array);
  
}

