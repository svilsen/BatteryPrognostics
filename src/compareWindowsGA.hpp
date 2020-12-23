#ifndef compareWindowsGA_H
#define compareWindowsGA_H

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

struct RandomVariates {
    unsigned int T1, T2, Wmin, Wmax, mutation_window;
    boost::mt19937 rng;
    
    boost::random::uniform_int_distribution<> uniform_T1;
    boost::random::uniform_int_distribution<> uniform_T2;
    boost::random::uniform_int_distribution<> uniform_W;
    boost::random::uniform_int_distribution<> uniform_c;
    boost::random::uniform_int_distribution<> uniform_m;
    
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_T1;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_T2;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_W;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_c;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_m;
    
    RandomVariates();
    RandomVariates(const unsigned int & T1_, const unsigned int & T2_, 
                   const unsigned int & Wmin_, const unsigned int & Wmax_, 
                   const unsigned int & mutation_window_, const unsigned int & seed);
};

class Individual {
public:
    std::vector<unsigned int> Value;
    std::vector<double> Fitness;
    unsigned int M;

    Individual(); 
    Individual(const std::vector<unsigned int> & Value_, const std::vector<double> & Fitness_);
};

class Population {
private:
    unsigned int N;
    
public:
    std::vector<Individual> Individuals;
    
    Population();
    Population(const unsigned int & N_);
    
    std::vector<unsigned int> FindDuplicates(const Individual & NewIndividual);
    int BestDuplicates(const Individual & NewIndividual, const std::vector<unsigned int> & duplicateIndicies);
    void RemoveDuplicates(const std::vector<unsigned int> & duplicateIndicies, const int & maxIndex);
    void Insert(const Individual & NewIndividual);
};

class ReducedWindowsGA {
private:
    std::vector<double> I1, I2;
    int T1, T2;
    unsigned int Wmax;
    double Imin;
    
public:
    std::vector<unsigned int> RI1, RI2;
    int RT1, RT2;
    
    ReducedWindowsGA();
    ReducedWindowsGA(const std::vector<double> & I1_, const std::vector<double> & I2_, 
                     const unsigned int & Wmax_, const double & Imin_);
    
    void Reduce();
};

class CompareWindowsGA {
private:
    RandomVariates RandomNumbers;
    ReducedWindowsGA ReducedWindows;
    std::vector<double> I1, I2, V1, V2, Temp1, Temp2;
    Individual NewIndividual, OldIndividual; //, Child1, Child2, Child3, Child4, Child5, Child6;
    std::vector<Individual> Children;
    
    int T1, T2;
    double restrict_temperature, restrict_voltage;
    unsigned int seed, Wmin, Wmax, N_evolution, N_keep, trace_limit;
    bool trace;

    //
    void CreateNewIndividual();
    std::vector<double> CalculateFitness(const std::vector<unsigned int> & newValue);
    
    void Crossover();
    
    void Mutate(Individual & individual, const unsigned int & m, const unsigned int & column);
    void Mutation();
    
    void Selection();
    
public:
    Population BestIndividuals;
    
    CompareWindowsGA(const std::vector<double> & I1_, const std::vector<double> & I2_,
                     const std::vector<double> & V1_, const std::vector<double> & V2_,
                     const std::vector<double> & Temp1_, const std::vector<double> & Temp2_, 
                     const double & restrict_temperature_, const double & restrict_voltage_,
                     const unsigned int & Wmin_, const unsigned int & Wmax_, 
                     const unsigned int & N_evolution_, const unsigned int & N_keep_, 
                     const bool & trace_, const unsigned int & trace_limit_, 
                     RandomVariates & RandomNumbers_, ReducedWindowsGA & ReducedWindows_);
    
    void EvolveWindows();
};

#endif