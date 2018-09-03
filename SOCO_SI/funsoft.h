#ifndef _FUNSOFTCOMPUTING_H

#define _FUNSOFTCOMPUTING_H 1

#ifndef _TFITNESS 
#define _TFITNESS 1
typedef long double Fitness;
#endif 

#define MAX_D 1000

extern Fitness Shifted_Sphere( int dim ,  double* x );
extern Fitness Schwefel_Problem(int dim, double* x);
extern Fitness Shifted_Rosenbrock(int dim, double* x);
extern Fitness Shifted_Rastrigin(int dim, double * x);
extern Fitness Shifted_Griewank(int dim, double * x);
extern Fitness Shifted_Ackley(int dim, double* x);
extern Fitness f_Schwefel2_22(int dim, double *s);
extern Fitness f_Schwefel1_2(int dim, double *s);
extern Fitness Extended_f_10(int dim, double *x);
extern Fitness f_Bohachevsky(int dim, double *s);
extern Fitness f_Schaffer(int dim, double *s);
extern Fitness f_Hybrid_12(int dim, double *s);
extern Fitness f_Hybrid_13(int dim, double *s);
extern Fitness f_Hybrid_14(int dim, double *s);
extern Fitness f_Hybrid_15(int dim, double *s);
extern Fitness f_Hybrid_16new(int dim, double *s);
extern Fitness f_Hybrid_17new(int dim, double *s);
extern Fitness f_Hybrid_18new(int dim, double *s);
extern Fitness f_Hybrid_19new(int dim, double *s);
extern Fitness f_FastFractal(int dim, double *s);

#endif
