// GeneExpress.cpp : Defines the entry point for the console application.
//

// stdafx is not included in the project files on SourceForge
#include "stdafx.h"

#include <Math.h>
#include "GEP.hh"

int _tmain(int argc, _TCHAR* argv[])
{

  // GEP Parameters
  unsigned gene_head_length = 7;
  unsigned genes_per_chromosome = 1;
  unsigned population_size = 200;
  unsigned constants_per_chromosome = 1;
  int max_epochs = 2000;
  
  std::string functions = "QES";
  
  double mutation_rate = 0.05;
  double RIS_rate = 0.1;
  double IS_rate = 0.1;
  double One_Point_rate = 0.3;
  double Two_Point_rate = 0.3;
  double Whole_Gene_rate = 0.1;
  double Constant_mutation = 0.1;
  double Constant_swap = 0.1;
  double Constant_range = 10.0;

  double validation_fraction = 0.1; // not implemented


  // Generate a training set of multivariate data.
  
  std::vector<std::vector<double> > data;
  std::vector<double> vals;

  for(int i=0;i<50;i++) {
	  double a = randomOne();
	  double b = randomOne();
	  vals.push_back(a);
	  vals.push_back(b);
	  vals.push_back(5*sin(a)+b); // last variable is the function value
	  data.push_back(vals);
	  vals.clear();
  }

  // Create the Gene Expression Program (GEP)

  GEP *gep = new GEP(
				  gene_head_length,
				  genes_per_chromosome,
				  constants_per_chromosome,
				  population_size,
				  validation_fraction,
				  max_epochs,
				  functions);

  // Set the mutation rates etc.

  gep->setRates(
	    mutation_rate,
		RIS_rate,
		IS_rate,
		One_Point_rate,
		Two_Point_rate,
		Whole_Gene_rate,
		Constant_mutation,
		Constant_swap,
		Constant_range);

  // Tell the GEP about the data

  gep->setData(&data);

  // Request the GEP to create a population of chromosomes

  gep->createPopulation();

  int verbose = 1;

  // Request the GEP to train the chromosomes to reproduce the data

  gep->train(verbose);






	return 0;
}

