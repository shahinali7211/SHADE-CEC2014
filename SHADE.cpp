//#include "./CEC2010/Header.h"
#include "Self_Define_Functions.h"

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

using namespace std;

double* OShift, * M, * y, * z, * x_bound;
int ini_flag = 0, n_flag, func_flag, * SS;

int main(int argc, char* argv[])
{
    int funToRun[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 };
    int funNum = 30;
    int i, j;
    const int Population_size = 100;
    const int timesOfRun = 51;
    const int dim = 30;
    const int MAX_FV = 10000 * dim;

    double CR = 0.9;
    double F = 0.8;
    double* global_best = new double[dim];
    double global_fitness;
    double** trail_vector = new double* [Population_size];
    double** mutation_vector = new double* [Population_size];
    double** population = new double* [Population_size];

    // Archive setup
    const int archive_size = (int)round(Population_size * 2);
    double** archive = new double* [archive_size];
    double* archive_fitness = new double[archive_size];
    int current_archive_size = 0;

    for (i = 0; i < Population_size; ++i) {
        population[i] = new double[dim];
        trail_vector[i] = new double[dim];
        mutation_vector[i] = new double[dim];
    }

    // Initialize archive arrays
    for (i = 0; i < archive_size; ++i) {
        archive[i] = new double[dim];
        archive_fitness[i] = 1e20;
    }

    double* results = new double[Population_size];
    double MAX = 100;
    double MIN = -100;
    double* trail_vector_result = new double[Population_size];
    int FV = 0;

    // SHADE parameters
    int H = dim;
    int k = 0;
    double* memory_CR = new double[H];
    double* memory_F = new double[H];
    vector<double> success_CR;  // Use vectors for dynamic sizing
    vector<double> success_F;
    vector<double> dif_fitness;

    double p = 0.1;
    const double epsilon = 1e-8;

    // new parameters sampling
    double* F_samp = new double[Population_size];
    double* CR_samp = new double[Population_size];

    // For current-to-pbest mutation
    int p_num = max((int)(round(Population_size * p)), 2);
    int* sorted_array = new int[Population_size];
    double* temp_fit = new double[Population_size];

    // --- Open CSV for summary results ---
    ofstream summary("results_50_dim.csv");
    summary << "Function,MeanFitness,StdDev\n";

    // --- Loop over all benchmark functions ---
    for (int function_index = 0; function_index < funNum; function_index++) {
        cout << "\n-------------------------------------------------------" << endl;
        cout << "Function = " << funToRun[function_index]
            << ", Dimension size = " << dim << "\n" << endl;

            boost::mt19937 generator(time(0) * rand());
            boost::uniform_real<> uniform_real_generate_x(MIN, MAX);
            boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);
            boost::uniform_int<> int_generator1(0, Population_size - 1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > int_number1(generator, int_generator1);
            boost::uniform_int<> random_memory_generator(0, H - 1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > random_memory_index(generator, random_memory_generator);

            // store fitness across runs
            double run_results[timesOfRun];
            double sum_fitness = 0.0;

            // --- Run multiple times ---
            for (int run_index = 0; run_index < timesOfRun; run_index++) {
                cout << run_index + 1 << "th run...   ";
                FV = 0;
                current_archive_size = 0;
                k = 0;

                // Clear success vectors for each run
                success_CR.clear();
                success_F.clear();
                dif_fitness.clear();

                for (i = 0; i < H; i++) {
                    memory_CR[i] = 0.5;
                    memory_F[i] = 0.5;
                }

                // initialize population
                for (i = 0; i < Population_size; ++i) {
                    for (j = 0; j < dim; ++j) {
                        population[i][j] = random_real_num_x();
                    }
                }

                for (i = 0; i < Population_size; i++) {
                    cec14_test_func(population[i], &results[i], dim, 1, funToRun[function_index]);
                    results[i] = results[i] - funToRun[function_index] * 100;
                }
                FV += Population_size;
                Find_best(results, Population_size, global_best, population, global_fitness, dim);

                // evolution loop
                while (FV < MAX_FV) {
                    // Clear success vectors for each generation
                    success_CR.clear();
                    success_F.clear();
                    dif_fitness.clear();
                    double sum_dif_fitness = 0.0;

                    // Assign individuals into sorted array
                    for (i = 0; i < Population_size; i++) {
                        sorted_array[i] = i;
                    }
                    // Assign fitness into temp array
                    for (i = 0; i < Population_size; i++) {
                        temp_fit[i] = results[i];
                    }
                    // Sort the temp array and keep the index in sorted_array
                    sortIndexWithQuickSort(temp_fit, 0, Population_size - 1, sorted_array);

                    for (i = 0; i < Population_size; ++i) {
                        int memory_index = random_memory_index();

                        // Sample CR value from normal distribution
                        boost::normal_distribution<> normal_cr_generate(memory_CR[memory_index], 0.1);
                        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > random_normal_cr(generator, normal_cr_generate);
                        if (memory_CR[memory_index] == -1) {
                            CR_samp[i] = 0;
                        }
                        else {
                            CR_samp[i] = random_normal_cr();
                            if (CR_samp[i] > 1) CR_samp[i] = 1;
                            else if (CR_samp[i] < 0) CR_samp[i] = 0;
                        }

                        // Sample F value from Cauchy distribution
                        boost::cauchy_distribution<> cauchy_f_generate(memory_F[memory_index], 0.1);
                        boost::variate_generator<boost::mt19937&, boost::cauchy_distribution<> > random_cauchy_f(generator, cauchy_f_generate);
                        do {
                            F_samp[i] = random_cauchy_f();
                        } while (F_samp[i] <= 0);
                        if (F_samp[i] > 1) F_samp[i] = 1;

                        // Generate pbest index
                        boost::uniform_int<> pbest_generator(0, p_num - 1);
                        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pbest_number(generator, pbest_generator);
                        int pbest_index = sorted_array[pbest_number()];

                        // Generate parent1 index
                        int parent1 = int_number1();
                        while (parent1 == i) parent1 = int_number1();

                        // Generate parent2 from population or archive
                        boost::uniform_int<> int_generator_arc(0, Population_size + current_archive_size - 1);
                        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > int_number_arc(generator, int_generator_arc);
                        int parent2 = int_number_arc();
                        while (parent2 == i || parent2 == parent1) parent2 = int_number_arc();

                        // Current-to-pbest/1 mutation
                        if (parent2 >= Population_size) { // Use archive individual
                            int arc_index = parent2 - Population_size;
                            for (j = 0; j < dim; ++j) {
                                mutation_vector[i][j] = population[i][j] +
                                    F_samp[i] * (population[pbest_index][j] - population[i][j]) +
                                    F_samp[i] * (population[parent1][j] - archive[arc_index][j]);

                                if (mutation_vector[i][j] < MIN) mutation_vector[i][j] = (MIN + population[i][j]) / 2.0;
                                if (mutation_vector[i][j] > MAX) mutation_vector[i][j] = (MAX + population[i][j]) / 2.0;
                            }
                        }
                        else { // Use population individual
                            for (j = 0; j < dim; ++j) {
                                mutation_vector[i][j] = population[i][j] +
                                    F_samp[i] * (population[pbest_index][j] - population[i][j]) +
                                    F_samp[i] * (population[parent1][j] - population[parent2][j]);

                                if (mutation_vector[i][j] < MIN) mutation_vector[i][j] = (MIN + population[i][j]) / 2.0;
                                if (mutation_vector[i][j] > MAX) mutation_vector[i][j] = (MAX + population[i][j]) / 2.0;
                            }
                        }

                        // Crossover - FIXED: use CR_samp[i] instead of CR
                        Crossover(trail_vector[i], population[i], mutation_vector[i], CR_samp[i], dim);

                        cec14_test_func(trail_vector[i], &trail_vector_result[i], dim, 1, funToRun[function_index]);
                        trail_vector_result[i] = trail_vector_result[i] - funToRun[function_index] * 100;
                        FV++;

                        if (trail_vector_result[i] < results[i]) {
                            // Store successful parameters
                            double diff = fabs(results[i] - trail_vector_result[i]);
                            dif_fitness.push_back(diff);
                            success_F.push_back(F_samp[i]);
                            success_CR.push_back(CR_samp[i]);
                            sum_dif_fitness += diff;

                            // Add parent to archive
                            if (current_archive_size < archive_size) {
                                memcpy(archive[current_archive_size], population[i], sizeof(double) * dim);
                                archive_fitness[current_archive_size] = results[i];
                                current_archive_size++;
                            }
                            else {
                                boost::uniform_int<> archive_generator(0, archive_size - 1);
                                int replace_index = archive_generator(generator);
                                memcpy(archive[replace_index], population[i], sizeof(double) * dim);
                                archive_fitness[replace_index] = results[i];
                            }

                            // Update population
                            results[i] = trail_vector_result[i];
                            memcpy(population[i], trail_vector[i], sizeof(double) * dim);
                        }

                        // Update global best
                        if (results[i] < global_fitness) {
                            if (results[i] <= epsilon) {
                                results[i] = 0;
                            }
                            global_fitness = results[i];
                            memcpy(global_best, population[i], sizeof(double) * dim);
                        }
                    }

                    // Update historical memories if successful parameters exist
                    if (!success_CR.empty()) {
                        double temp_sum_sf = 0.0;
                        double temp_sum_cr = 0.0;
                        double weighted_sum_F = 0.0;
                        double weighted_sum_F2 = 0.0;
                        double weighted_sum_CR = 0.0;

                        for (size_t idx = 0; idx < success_CR.size(); idx++) {
                            double weight = dif_fitness[idx] / sum_dif_fitness;
                            weighted_sum_F2 += weight * success_F[idx] * success_F[idx];
                            weighted_sum_F += weight * success_F[idx];
                            weighted_sum_CR += weight * success_CR[idx];
                        }

                        // Update F memory (Lehmer mean)
                        if (weighted_sum_F != 0) {
                            memory_F[k] = weighted_sum_F2 / weighted_sum_F;
                        }

                        // Update CR memory
                        if (weighted_sum_CR == 0) {
                            memory_CR[k] = -1;
                        }
                        else {
                            memory_CR[k] = weighted_sum_CR;
                        }

                        k = (k + 1) % H;
                    }
                }
                cout << "Final error value = " << global_fitness << endl;
                run_results[run_index] = global_fitness;
                sum_fitness += global_fitness;
            }

            // mean & std
            double mean_fitness = sum_fitness / timesOfRun;
            double std_fitness = 0.0;
            for (int r = 0; r < timesOfRun; r++)
                std_fitness += pow(run_results[r] - mean_fitness, 2.0);
            std_fitness = sqrt(std_fitness / timesOfRun);

            cout << "\nMean = " << mean_fitness << ", Std = " << std_fitness << endl;
            summary << "Function " << funToRun[function_index] << "," << mean_fitness << "," << std_fitness << "\n";
    }

    summary.close();

    // release resources
    for (i = 0; i < Population_size; ++i) {
        delete[] population[i];
        delete[] trail_vector[i];
        delete[] mutation_vector[i];
    }
    for (i = 0; i < archive_size; ++i) {
        delete[] archive[i];
    }
    delete[] archive;
    delete[] archive_fitness;
    delete[] population;
    delete[] trail_vector;
    delete[] mutation_vector;
    delete[] results;
    delete[] global_best;
    delete[] memory_CR;
    delete[] memory_F;
    delete[] F_samp;
    delete[] CR_samp;
    delete[] sorted_array;
    delete[] temp_fit;
    delete[] trail_vector_result;

    return 0;
}