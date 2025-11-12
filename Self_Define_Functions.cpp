
#include "Self_Define_Functions.h"

// create new object of class for different functions
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

bool Compare_NewType(NewType data1, NewType data2)
{
	return data1.data < data2.data;
}

using namespace std;
//Fitness Computation
//void Mutation(double *mutation_vector, double *individual, double *parent1, double *parent2, double *global_best, double MIN, double MAX, double F, int dim)
//{
//	int i;
//
//	for (i = 0; i < dim; ++i)
//	{
//		mutation_vector[i] = individual[i] + F * (global_best[i] - individual[i]) + F * (parent1[i] - parent2[i]);
//
//		if (mutation_vector[i] < MIN) {
//			mutation_vector[i] = (MIN + individual[i]) / 2.0;
//		}
//		if (mutation_vector[i] > MAX) {
//			mutation_vector[i] = (MAX + individual[i]) / 2.0;
//		}
//	}
//}


void Crossover(double *trial_vector, double *individual, double *mutation_vector, double CR, int dim)
{
	int i;

	boost::mt19937 generator(time(0)*rand());
	boost::uniform_int<> int_generator(0, dim - 1);
	boost::variate_generator< boost::mt19937&, boost::uniform_int<> > int_number(generator, int_generator);


	boost::uniform_real<> real_generator(0, 1);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > real_number(generator, real_generator);

	int dim_rand = int_number();

	for (i = 0; i < dim; ++i)
	{
		if (i == dim_rand || real_number() <= CR)
			trial_vector[i] = mutation_vector[i];
		else
			trial_vector[i] = individual[i];
	}

}

void Find_best(double *results, int population_size, double* global_best, double** population, double& global_fitness, int dim)
{
	int i;
	NewType *temp = new NewType[population_size];
	for (i = 0; i < population_size; ++i)
	{
		temp[i].data = results[i];
		temp[i].id = i;
	}

	sort(temp, temp + population_size, Compare_NewType);
	global_fitness = temp[0].data;
	memcpy(global_best, population[temp[0].id], sizeof(double)*dim);

	delete[]temp;
}





