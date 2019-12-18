//
//
// Description:
//      Class GEP :
//		Gene Expression simulation
//		Encoded Functions:
//			+,-,/,* (these are always used, the rest are optional)
//			E exponential
//			Q square root
//			L less than. L(a,b) = 1 if a<b, 0 otherwise
//			G greater than. G(a,b) = 1 if a>b, 0 otherwise
//			A abs
//			O Log
//			S Sin
//
//Each chromosome is encoded as one or more genes that each have a Head and Tail. 
//Each gene's head contains a set of functions, variables or constants
//Each gene's tail contains terminals or constants
//The length of the Tail depends on the maximum arity of the function set, and
//is computed when the GEP is created with a specified Head length. The Tail
//length = (Head length)*(Arity - 1) + 1
//All functions, variables and constants are encoded as integers
//	Integers >= VARIABLE := a variable
//	Integers = CONSTANT := a constant
//	Integers < CONSTANT := a function
//The chromosome's linking function determines 
//how the gene values are combined to form the
//output value of the chromosome
//Constants are generated in the range 0 to Constant_range
//
// Environment:
//      Software developed at Caltech's Center for Advanced Computing Research
//
// Author List:
//      Julian Bunn                     Original author
//
// Copyright Information:
//      Copyright (C) 2009              California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _GEP_HH
#define _GEP_HH

#include "Chromosome.hh"

#include <string>
#include <iostream>
#include <vector>
#include <utility>

double randomOne();


class GEP 
{
public:
  virtual ~GEP() {}

  GEP(
	 unsigned gene_head_length,
	 unsigned genes_per_chromosome,
	 unsigned constants_per_chromosome,
	 unsigned population_size,
	 double validation_fraction,
	 int Max_Epoch,
	 std::string functions,
	 int seed=0);

  std::string name() const { return "GEP"; }

  bool train(int verbose=0);

  bool reset();

  bool setData(const std::vector<std::vector<double> > *data);

  void setValidationFraction(double frac);

  void print(std::ostream& os) const;
 
  bool printValidation(unsigned cycle);

  void createPopulation();

  void createFunctionSet(const std::string& functions);

  void setRates(double mutation_rate,
		double RIS_rate,
		double IS_rate,
		double One_Point_rate,
		double Two_Point_rate,
		double Whole_Gene_rate,
		double Constant_mutation,
		double Constant_swap,
		double Constant_range);
	  
  // These functions convert between the character 
  // and integer representations of GEP functions
  static char charFromFunction(int);
  static int functionFromChar(char);
  static char charFromInt(int);

  // FUNCTIONS

  std::vector<int> FunctionsInUse;

  std::string FunctionSet;

  enum Function 
    {Plus, Minus, Times, Divide, Sqrt, Exp, LessThan, GreaterThan, Abs, Log, Sin};

  static const int Variable;
  static const int Constant; 
  static const int Max_Depth;

  static const double Very_Large_Double;
  static const double Very_Small_Double;



private:

  // TRAINING AND VALIDATION DATA

  std::vector<std::vector<double> > Dataset;
  std::vector<std::vector<double> > ValidationDataset;

  unsigned gene_head_length_;
  unsigned genes_per_chromosome_;
  unsigned constants_per_chromosome_;
  unsigned population_size_;
  int Max_Epoch_;
  unsigned function_arg_max_;
  std::vector<Chromosome> chromosomes_;

  double validation_fraction_; // fraction of the data to use for validation (not training)

  double mutation_rate_;// how often a chromosome mutates
  double RIS_rate_;	// how often root sequence insertion happens
  double IS_rate_;	// how often sequence insertion happens
  double One_Point_rate_;// sum of the three recombination rates ~ 0.7
  double Two_Point_rate_;// two point recombination (swap a portion between chromosomes)
  double Whole_Gene_rate_;// swap whole genes between chromosomes
  double Constant_mutation_;// How often a chromosome will have all its constants replaced
  double Constant_swap_;// How often two chromosomes will swap a constant
  double Constant_range_;// Range for the generated constants

  Chromosome best_;

  int valPrint_;
};

#endif
