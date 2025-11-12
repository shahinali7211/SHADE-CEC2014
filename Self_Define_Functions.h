#ifndef SELF_DEFINE_FUNCTIONS_H_INCLUDED
#define SELF_DEFINE_FUNCTIONS_H_INCLUDED
#include <math.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

struct NewType
{
	double data;
	int id;
};

//const int Population_size = 100;//Population size
//const int timesOfRun = 3;// run time
//const int dim = 30;//the dimension size of the problem to be optimized
//const int MAX_FV = 10000*dim;//the maximum number of fitness evaluations


void cec14_test_func(double *x, double *f, int nx, int mx, int func_num);
bool Compare_NewType(NewType data1, NewType data2);
//void Mutation(double *mutation_vector, double *individual, double *parent1, double *parent2, double *global_best, double MIN, double MAX, double F, int dim);
void Crossover(double *trial_vector, double *individual, double *mutation_vector, double CR, int dim);
void Find_best(double *results, int population_size, double* global_best, double** population, double& global_fitness, int dim);



//Recursive quick sort with index array
template<class VarType>
void sortIndexWithQuickSort(VarType array[], int first, int last, int index[]) {
    VarType x = array[(first + last) / 2];
    int i = first;
    int j = last;
    VarType temp_var = 0;
    int temp_num = 0;

    while (true) {
        while (array[i] < x) i++;
        while (x < array[j]) j--;
        if (i >= j) break;

        temp_var = array[i];
        array[i] = array[j];
        array[j] = temp_var;

        temp_num = index[i];
        index[i] = index[j];
        index[j] = temp_num;

        i++;
        j--;
    }

    if (first < (i - 1)) sortIndexWithQuickSort(array, first, i - 1, index);
    if ((j + 1) < last) sortIndexWithQuickSort(array, j + 1, last, index);
}
#endif // SELF_DEFINE_FUNCTIONS_H_INCLUDED
