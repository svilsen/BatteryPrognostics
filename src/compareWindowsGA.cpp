#include <Rcpp.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "compareWindowsGA.hpp"

// Does individual f dominate g...
bool Dominates(const Individual & f, const Individual & g, unsigned int m)
{
    if (m == f.M)
    {
        return false;
    }
    else
    {
        bool dominates = false;
        if (f.Fitness[m] > g.Fitness[m])
        {
            dominates = true;
        }
        else if (f.Fitness[m] == g.Fitness[m])
        {
            dominates = Dominates(f, g, m + 1);
        }
        
        return dominates;
    }
}

bool IsEqual(const Individual & f, const Individual & g) 
{
    bool is_equal = false;
    double f_0 = static_cast<double>(f.Value[0]);
    double f_1 = static_cast<double>(f.Value[1]);
    
    double g_0 = static_cast<double>(g.Value[0]);
    double g_1 = static_cast<double>(g.Value[1]);
    
    double diff_0 = f_0 - g_0;
    double diff_1 = f_1 - g_1;
    
    if (((std::abs(diff_0) < 1e-8) & (std::abs(diff_1) < 1e-8)) || 
        (f.Value[0] < (g.Value[0] + g.Value[2])) & (f.Value[0] > g.Value[0]) ||
        (g.Value[0] < (f.Value[0] + f.Value[2])) & (g.Value[0] > f.Value[0]) || 
        (f.Value[1] < (g.Value[1] + g.Value[2])) & (f.Value[1] > g.Value[1]) || 
        (g.Value[1] < (f.Value[1] + f.Value[2])) & (g.Value[1] > f.Value[1]))
    {
        is_equal = true;
    }
    
    return is_equal;
}

unsigned int Contained(const Individual & f, const Individual & g) 
{
    int contained = -1;
    if ((f.Value[0] >= g.Value[0]) & ((f.Value[0] + f.Value[2]) <= (g.Value[0] + g.Value[2]))) 
    {
        // f is contained in g
        contained = 0;
    }
    else if ((f.Value[0] <= g.Value[0]) & ((f.Value[0] + f.Value[2]) >= (g.Value[0] + g.Value[2]))) 
    {
        // g is contained in f
        contained = 1;
    }
    
    return contained;
}

// Insertion-sort for populations...
void InsertionSort(std::vector<Individual> & Individuals) 
{ 
    unsigned int N = Individuals.size();
    for (unsigned int n = 1; n < N; n++) 
    {
        Individual key = Individuals[n];
        int m = n - 1;
        bool dominates = Dominates(Individuals[m], key, 0);
        while ((m > -1) & (dominates)) 
        {
            Individuals[m + 1] = Individuals[m];
            m--;
        }
        
        Individuals[m + 1] = key;
    }
} 

// Random variate structure
RandomVariates::RandomVariates() : 
    T1(1), T2(1), Wmin(1), Wmax(2), mutation_window(10), rng(0), 
    uniform_T1(0, T1), generate_uniform_T1(rng, uniform_T1),
    uniform_T2(0, T2), generate_uniform_T2(rng, uniform_T2),
    uniform_W(Wmin, Wmax), generate_uniform_W(rng, uniform_W), 
    uniform_c(1, 2), generate_uniform_c(rng, uniform_c), 
    uniform_m(5, mutation_window), generate_uniform_m(rng, uniform_m) { };
RandomVariates::RandomVariates(const unsigned int & T1_, const unsigned int & T2_, 
                               const unsigned int & Wmin_, const unsigned int & Wmax_, 
                               const unsigned int & mutation_window_, const unsigned int & seed) : 
    T1(T1_), T2(T2_), Wmin(Wmin_), Wmax(Wmax_), mutation_window(mutation_window_), rng(seed), 
    uniform_T1(0, T1 - Wmax - 2), generate_uniform_T1(rng, uniform_T1),
    uniform_T2(0, T2 - Wmax - 2), generate_uniform_T2(rng, uniform_T2),
    uniform_W(Wmin, Wmax), generate_uniform_W(rng, uniform_W), 
    uniform_c(1, 2), generate_uniform_c(rng, uniform_c), 
    uniform_m(std::min(static_cast<int>(10), static_cast<int>(mutation_window / 2)), mutation_window), 
    generate_uniform_m(rng, uniform_m) { }

// Individual class
Individual::Individual() 
{ 
    Value = std::vector<unsigned int>(3, 0);
    Fitness = std::vector<double>(5, -HUGE_VAL);
    M = 5;
}

Individual::Individual(const std::vector<unsigned int> & Value_, const std::vector<double> & Fitness_) :
    Value(Value_), Fitness(Fitness_), M(Fitness.size()) { }

// Population class
Population::Population() 
{
    std::vector<Individual> pop(1);
    Individuals = pop;
}

Population::Population(const unsigned int & N_) : 
    N(N_) 
{ 
    Individuals = std::vector<Individual>(N);
}

std::vector<unsigned int> Population::FindDuplicates(const Individual & NewIndividual)
{
    std::vector<unsigned int> duplicateIndicies;
    for (unsigned int n = 0; n < N; n++) 
    {
        const bool is_equal = IsEqual(NewIndividual, Individuals[n]);
        
        if (is_equal)
        {
            duplicateIndicies.push_back(n);
        }
        else 
        {
            const int contained = Contained(NewIndividual, Individuals[n]);
            if (contained > -1) 
            {
                duplicateIndicies.push_back(n);
            }
        }
    }
    
    return duplicateIndicies;
}

int Population::BestDuplicates(const Individual & NewIndividual, const std::vector<unsigned int> & duplicateIndicies)
{
    const unsigned int M = duplicateIndicies.size();
    Individual MaxIndividual = NewIndividual;
    int MaxIndex = -1;
    
    for (unsigned int m = 0; m < M; m++)
    {
        const Individual & Individual_m = Individuals[duplicateIndicies[m]];
        const bool dominates = Dominates(Individual_m, MaxIndividual, 0);
        if (dominates) 
        {
            MaxIndividual = Individual_m;
        }
    }
    
    return MaxIndex;
}

void Population::RemoveDuplicates(const std::vector<unsigned int> & duplicateIndicies, const int & maxIndex)
{
    const unsigned int M = duplicateIndicies.size();
    const Individual I_empty;
    
    for (unsigned int m = 0; m < M; m++)
    {
        if (m != maxIndex) 
        {
            Individuals[duplicateIndicies[m]] = I_empty;
        }
    }
}


void Population::Insert(const Individual & NewIndividual) 
{
    std::vector<unsigned int> duplicateIndices = FindDuplicates(NewIndividual);
    int maxIndex = BestDuplicates(NewIndividual, duplicateIndices);
    RemoveDuplicates(duplicateIndices, maxIndex);
    InsertionSort(Individuals);
    
    if (maxIndex == -1) 
    {
        bool dominates = Dominates(NewIndividual, Individuals[0], 0);
        if (dominates) 
        {
            bool not_stopped = true;
            unsigned int j = 0;
            while (not_stopped) 
            {
                dominates = Dominates(NewIndividual, Individuals[j + 1], 0);
                if (dominates) 
                {
                    j++;
                    if ((j + 1) == N)
                    {
                        not_stopped = false;
                    } 
                }
                else 
                {
                    not_stopped = false; 
                }
            }
            
            if (j > 0) 
            {   
                for (std::size_t k = 0; k < j; k++)
                {
                    Individuals[k] = Individuals[k + 1];
                }
                
                Individuals[j] = NewIndividual;
            }
            else 
            {
                Individuals[j] = NewIndividual;
            }
        }
    }
}

// Reduce data class 
ReducedWindowsGA::ReducedWindowsGA(const std::vector<double> & I1_, const std::vector<double> & I2_, 
                                   const unsigned int & Wmax_, const double & Imin_) : 
    I1(I1_), I2(I2_), T1(I1.size()), T2(I2.size()), Wmax(Wmax_), Imin(Imin_)
{ 
    RT1 = 0;
    RT2 = 0;
    
    RI1 = std::vector<unsigned int>(T1 - Wmax);
    RI2 = std::vector<unsigned int>(T2 - Wmax);
    
    Reduce();
}

void ReducedWindowsGA::Reduce()
{
    const unsigned int & T = std::max(T1, T2);
    for (unsigned int t = 0; t < (T - Wmax); t++)  
    {
        if (t < (T1 - Wmax)) 
        {
            if (I1[t] > Imin) 
            {
                RI1[RT1] = t;
                RT1++;
            }
        }
        
        if (t < (T2 - Wmax)) 
        {
            if (I2[t] > Imin) 
            {
                RI2[RT2] = t;
                RT2++;
            }
        }
    }
    
    RI1.resize(RT1);
    RI1.shrink_to_fit();
    
    RI2.resize(RT2);
    RI2.shrink_to_fit();
}


// Compare windows GA class
void CompareWindowsGA::CreateNewIndividual()
{
    std::vector<unsigned int> newValue(3, 0);
    
    newValue[0] = ReducedWindows.RI1[RandomNumbers.generate_uniform_T1()];
    newValue[1] = ReducedWindows.RI2[RandomNumbers.generate_uniform_T2()];
    newValue[2] = RandomNumbers.generate_uniform_W();
    
    std::vector<double> newFitness = CalculateFitness(newValue);
    
    Individual newIndividual(newValue, newFitness);
    OldIndividual = NewIndividual;
    NewIndividual = newIndividual;
}

CompareWindowsGA::CompareWindowsGA(const std::vector<double> & I1_, const std::vector<double> & I2_,
                                   const std::vector<double> & V1_, const std::vector<double> & V2_,
                                   const std::vector<double> & Temp1_, const std::vector<double> & Temp2_, 
                                   const double & restrict_temperature_, const double & restrict_voltage_,
                                   const unsigned int & Wmin_, const unsigned int & Wmax_, 
                                   const unsigned int & N_evolution_, const unsigned int & N_keep_,
                                   const bool & trace_, const unsigned int & trace_limit_, 
                                   RandomVariates & RandomNumbers_, ReducedWindowsGA & ReducedWindows_) : 
    I1(I1_), I2(I2_), T1(I1.size()), T2(I2.size()), V1(V1_), V2(V2_), 
    Temp1(Temp1_), Temp2(Temp2_), 
    restrict_temperature(restrict_temperature_), restrict_voltage(restrict_voltage_),
    Wmin(Wmin_), Wmax(Wmax_), 
    N_evolution(N_evolution_), N_keep(N_keep_),  trace(trace_), trace_limit(trace_limit_), 
    RandomNumbers(RandomNumbers_), ReducedWindows(ReducedWindows_)
{ 
    CreateNewIndividual();
    BestIndividuals = Population(N_keep);
    
    Children = std::vector<Individual>(6);
}

void CompareWindowsGA::EvolveWindows()
{
    for (unsigned int n = 0; n < N_evolution; n++)
    {
        if ((n == 0) || (((n + 1) % trace_limit) == 0) || (n == (N_evolution - 1)))
        {
            Rcpp::Rcout << "Iteration: " << n + 1 << "\n";
        }
        
        CreateNewIndividual();
        
        Crossover();
        Mutation();
        
        Selection();
        
        BestIndividuals.Insert(NewIndividual);
        for (unsigned int c = 0; c < Children.size(); c++) 
        {
            BestIndividuals.Insert(Children[c]);
        }
    }
}

std::vector<double> CompareWindowsGA::CalculateFitness(const std::vector<unsigned int> & newValue) 
{
    std::vector<double> newFitness(5, -HUGE_VAL);
    const double abs_difference_voltage = std::abs(V1[newValue[0]] - V2[newValue[1]]) +
        std::abs(V1[newValue[0] + newValue[2]] - V2[newValue[1] + newValue[2]]);
    if (abs_difference_voltage < restrict_voltage) 
    {
        const double abs_difference_temperature = std::abs(Temp1[newValue[0]] - Temp2[newValue[1]]);
        if (abs_difference_temperature < restrict_temperature) 
        {
            newFitness[0] = 0.0;
            newFitness[1] = 0.0;
            newFitness[2] = 0.0;
            newFitness[3] = 0.0;
            newFitness[4] = 0.0;
            
            for (unsigned int w = 0; w < newValue[2]; w++)
            {
                double abs_difference_current = std::abs((I1[newValue[0] + w] - I2[newValue[1] + w]) / I1[newValue[0] + w]);
                if (abs_difference_current > newFitness[0])
                {
                    newFitness[0] = abs_difference_current;
                }
                
                newFitness[4] += std::abs(I1[newValue[0] + w]) / newValue[2];
            }
            
            // 
            newFitness[0] = -newFitness[0];
            newFitness[1] = newValue[2];
            newFitness[2] = -abs_difference_voltage;
            newFitness[3] = -abs_difference_temperature;
        }
    }
    
    return newFitness;
}

void CompareWindowsGA::Crossover() 
{
    unsigned int crossover_point = RandomNumbers.generate_uniform_c();
    std::vector<unsigned int> childValue1(3, 0);
    std::vector<unsigned int> childValue2(3, 0);
    for (unsigned int i = 0; i < 3; i++) 
    {
        if (i < crossover_point) 
        {
            childValue1[i] = OldIndividual.Value[i];
            childValue2[i] = NewIndividual.Value[i];
        }
        else 
        {
            childValue2[i] = OldIndividual.Value[i];
            childValue1[i] = NewIndividual.Value[i];
        }
    }
    
    std::vector<double> childFitness1 = CalculateFitness(childValue1);
    std::vector<double> childFitness2 = CalculateFitness(childValue2);
    
    Individual child1(childValue1, childFitness1);
    Individual child2(childValue2, childFitness2);
    
    Children[0] = child1;
    Children[1] = child2;
}

void CompareWindowsGA::Mutate(Individual & individual, const unsigned int & m, const unsigned int & column) 
{
    std::vector<unsigned int> newValue = individual.Value;
    int min_value = 0;
    int max_value;
    if (column == 0) 
    {
        max_value = ReducedWindows.RI1[ReducedWindows.RT1 - 1];
    }
    else if (column == 1) 
    {
        max_value = ReducedWindows.RI2[ReducedWindows.RT2 - 1];
    }
    else 
    {
        min_value = Wmin;
        max_value = Wmax;
    }
    
    unsigned int newVal = std::min(std::max(min_value, static_cast<int>(newValue[column]) + static_cast<int>(m)), max_value);
    newValue[column] = newVal;
    
    std::vector<double> newFitness = CalculateFitness(newValue);
    Individual newIndividual(newValue, newFitness);
    
    individual = newIndividual;
}

void CompareWindowsGA::Mutation() 
{
    Children[2] = OldIndividual;
    Children[3] = NewIndividual;
    
    Children[4] = Children[0];
    Children[5] = Children[1];
    
    const double mutation_w = RandomNumbers.generate_uniform_m();
    unsigned int column = RandomNumbers.generate_uniform_c();
    for (unsigned int m = -mutation_w; m < mutation_w; m++) 
    {
        for (unsigned int c = 2; c < Children.size(); c++) 
        {
            Mutate(Children[c], m, column);
        }
    }
}

void CompareWindowsGA::Selection() 
{
    bool dominates;
    for (unsigned int c = 0; c < Children.size(); c++) 
    {
        dominates = Dominates(Children[c], NewIndividual, 0);
        if (dominates) 
        {
            NewIndividual = Children[c];
        }
    }
}

//[[Rcpp::export()]]
Rcpp::List compare_windows_ga_cpp(const std::vector<double> & I1, const std::vector<double> & I2,
                                  const std::vector<double> & V1, const std::vector<double> & V2,
                                  const std::vector<double> & Temp1, const std::vector<double> & Temp2,
                                  const unsigned int & mutation_window, const double & restrict_temperature, 
                                  const double & restrict_voltage, 
                                  const unsigned int & Wmin, const unsigned int & Wmax, const unsigned int & Imin,
                                  const unsigned int & N_evolution, const unsigned int & N_keep,
                                  const bool & trace, const unsigned int & trace_limit, 
                                  const unsigned int & seed)
{
    //
    ReducedWindowsGA ReducedWindows(I1, I2, Wmax, Imin);
    
    RandomVariates RandomNumbers(ReducedWindows.RT1 + Wmax + 1, 
                                 ReducedWindows.RT2 + Wmax + 1, 
                                 Wmin, Wmax, mutation_window, seed);
    
    CompareWindowsGA CWGA(I1, I2, V1, V2, Temp1, Temp2, restrict_temperature, restrict_voltage, Wmin, Wmax, 
                          N_evolution, N_keep, trace, trace_limit, RandomNumbers, ReducedWindows);
    
    CWGA.EvolveWindows();
    
    //
    std::vector<Individual> BestIndividuals = CWGA.BestIndividuals.Individuals;
    std::vector<std::vector<unsigned int> > Individuals(N_keep);
    std::vector<std::vector<double> > Fitness(N_keep);
    for (unsigned int n = 0; n < N_keep; n++)
    {
        Individuals[n] = BestIndividuals[n].Value;
        Fitness[n] = BestIndividuals[n].Fitness;
    }
    
    return Rcpp::List::create(Rcpp::Named("Individuals") = Individuals,
                              Rcpp::Named("Fitness") = Fitness);
}

