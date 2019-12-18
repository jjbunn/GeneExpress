
#include "GEP.hh"
#include "Chromosome.hh"

#include <stdio.h>
#include <cassert>
#include <algorithm>
#include <functional>
#include <cmath>
#include <sstream>

using namespace std;

const int GEP::Variable = 500;
const int GEP::Constant = 400; 
const int GEP::Max_Depth = 40;

const double GEP::Very_Large_Double = 1e20;
const double GEP::Very_Small_Double = 1e-20;

double randomOne() {
	// return a random number in the range 0-0.9999999999
	return (double) (rand()-1) / (double) RAND_MAX;
}


GEP::GEP(
	       unsigned gene_head_length,
	       unsigned genes_per_chromosome,
	       unsigned constants_per_chromosome,
	       unsigned population_size,
		   double validation_fraction,
	       int Max_Epoch,
	       string functions,
	       int seed)
  : 
  gene_head_length_(gene_head_length),
  genes_per_chromosome_(genes_per_chromosome),
  constants_per_chromosome_(constants_per_chromosome),
  population_size_(population_size),
  Max_Epoch_(Max_Epoch),
  function_arg_max_(0),
  chromosomes_(population_size,Chromosome()),
  mutation_rate_(0.05), // 0.05
  RIS_rate_(0.1), // 0.1
  IS_rate_(0.1), // 0.1
  One_Point_rate_(0.3), // 0.3
  Two_Point_rate_(0.3), // 0.3
  Whole_Gene_rate_(0.1), // 0.1
  Constant_mutation_(0.1), // 0.1 0.7
  Constant_swap_(0.1), // 0.1
  Constant_range_(10.0) // 10.0
{
  assert (gene_head_length > 1);
  assert (genes_per_chromosome > 0);
  assert (population_size > 1);

  cout << "Begin GEP" << endl 
       << " Genes per Chromosome         = " 
       << genes_per_chromosome_ << endl 
       << " Gene head length             = " <<	gene_head_length_ << endl 
       << " Constants per Chromosome     = " 
       << constants_per_chromosome_ << endl 
       << " Chromosome population size   = " << population_size_ << endl 
       << " Extra functions requested    = " << functions << endl 
       << " Maximum Epochs/Generations   = " << Max_Epoch_ << endl;

  FunctionSet = "+-*/AQELGOS";

  this->createFunctionSet(functions);
}

void GEP::setRates(double mutation_rate,
		      double RIS_rate,
		      double IS_rate,
		      double One_Point_rate,
		      double Two_Point_rate,
		      double Whole_Gene_rate,
		      double Constant_mutation,
		      double Constant_swap,
		      double Constant_range) 
{
  mutation_rate_ = mutation_rate;
  RIS_rate_ = RIS_rate;
  IS_rate_ = IS_rate;
  One_Point_rate_ = One_Point_rate;
  Two_Point_rate_ = Two_Point_rate;
  Whole_Gene_rate_ = Whole_Gene_rate;
  Constant_mutation_ = Constant_mutation;
  Constant_swap_ = Constant_swap;
  Constant_range_ = Constant_range;

  cout <<
    " Mutation Rate                = " << mutation_rate_ << endl <<
    " Root Insertion Rate          = " << RIS_rate_ << endl <<
    " Insertion Rate               = " << IS_rate_ << endl <<
    " One Point Transposition Rate = " << One_Point_rate_ << endl <<
    " Two Point Transposition Rate = " << Two_Point_rate_ << endl <<
    " Whole Gene Transpose Rate    = " << Whole_Gene_rate_ << endl <<
    " Constant Mutation Rate       = " << Constant_mutation_ << endl <<
    " Constant Swap Rate           = " << Constant_swap_ << endl <<
    " Constant Generation Range    = " << Constant_range_ << endl;
  
}


void GEP::createFunctionSet(const std::string& functions) 
{
  // The functions +-*/ are always included
  FunctionsInUse.clear();
  FunctionsInUse.push_back(Plus);
  FunctionsInUse.push_back(Minus);
  FunctionsInUse.push_back(Times);
  FunctionsInUse.push_back(Divide);
  
  cout << "Using functions +-*/";
  for(int i=0;i<functions.size();i++) {
    if( FunctionSet.find(functions[i]) != string::npos ) {
      int f = functionFromChar(functions[i]);
      if( find(FunctionsInUse.begin(),FunctionsInUse.end(),f) 
	  == FunctionsInUse.end()) {
	FunctionsInUse.push_back(f);
	cout << functions[i];
      }
    } 
    else {
      cerr << " Unknown function \"" << functions[i] 
	   << "\" requested." << endl;
    }
  }
  cout << endl;
  function_arg_max_ = 2;
}

bool GEP::setData(const std::vector<std::vector<double> > *_data) {
	Dataset = *_data;
	return true;
}

void GEP::setValidationFraction(double _frac) {
	validation_fraction_ = _frac;
}



int GEP::functionFromChar(char f) 
{
  switch(f) 
    {
    case '+':
      return GEP::Plus;
    case '-':
      return GEP::Minus;
    case '*':
      return GEP::Times;
    case '/':
      return GEP::Divide;
    case 'Q':
      return GEP::Sqrt;
    case 'E':
      return GEP::Exp;
    case 'S':
      return GEP::Sin;
    case 'L':
      return GEP::LessThan;
    case 'G':
      return GEP::GreaterThan;
    case 'A':
      return GEP::Abs;
    case 'O':
      return GEP::Log;

    default:
      return -1;
    }
  return -1;
}


char GEP::charFromFunction(int f) 
{
  switch(f) 
    {
    case GEP::Plus:
      return '+';
    case GEP::Minus:
      return '-';
    case GEP::Times:
      return '*';
    case GEP::Divide:
      return '/';
    case GEP::Sqrt:
      return 'Q';
    case GEP::Exp:
      return 'E';
    case GEP::Sin:
      return 'S';
    case GEP::LessThan:
      return 'L';
    case GEP::GreaterThan:
      return 'G';
    case GEP::Abs:
      return 'A';	
    case GEP::Log:
      return 'O';	
    default:
      return '!';
    }
}


void GEP::createPopulation() 
{
  cout << "Creating Chromosome population of size " 
       << population_size_ << endl;
  chromosomes_.clear();
  chromosomes_.resize(population_size_,Chromosome());

  unsigned size = Dataset.size();
  unsigned dim = Dataset[0].size() - 1; // the last variable is the target chromosome value

  unsigned head_length = gene_head_length_;
  unsigned tail_length = gene_head_length_ * (function_arg_max_ - 1) + 1;
  unsigned ch_length = gene_head_length_ + tail_length;

  cout << "Each gene will have a tail of length " << tail_length 
       << " and a total length of " << ch_length << endl;

  // Create the chromosomes, each of 1 or more genes
  int nchromosome = 0;
  for( vector<Chromosome>::iterator it=chromosomes_.begin();
       it!=chromosomes_.end();it++ ) {
    Chromosome &thisChromosome = *it;
    thisChromosome.setLength(genes_per_chromosome_,
			     constants_per_chromosome_,
			     genes_per_chromosome_*ch_length);

    for( int ch=0;ch<genes_per_chromosome_;ch++ ) {
      vector<int> head;
      vector<int> tail;

      for( int ih=0;ih<head_length;ih++ ) {
	// head can contain functions or terminals (variables or constants)
	// first position in head must be a function
	if( ih == 0 ) {
	  int ipos = rand() % FunctionsInUse.size();
	  head.push_back(FunctionsInUse[ipos]);
	} 
	else {
	  // twice as many head positions will be functions
	  if( randomOne() > 0.5 ) {
	    int ipos = rand() % FunctionsInUse.size();
	    head.push_back(FunctionsInUse[ipos]);
	  } 
	  else {
	    int ipos 
	      = rand() % (dim + constants_per_chromosome_);
	    if( ipos >= dim ) {
	      head.push_back(Constant);
	    } 
	    else {
	      int var = Variable + ipos;
	      head.push_back(var);
	    }
	  }
	}
      }
      
      for(int j=0;j<tail_length;j++) {
	int idim = rand() % (dim + constants_per_chromosome_);
	if( idim >= dim ) {
	  tail.push_back(Constant);
	} 
	else {
	  int var = Variable + idim;
	  tail.push_back(var);
	}
      }
      
      thisChromosome.gene_[ch].setHead(head);
      thisChromosome.gene_[ch].setTail(tail);	
      
    }
    // Generate random constants for this chromosome if required
    if( constants_per_chromosome_ > 0 ) {
      thisChromosome.constants_.clear();
      for(int ic=0;ic<constants_per_chromosome_;ic++) {
	double val = Constant_range_*randomOne();
	thisChromosome.constants_.push_back(val);
      }
      //cout << "  Constants generated: " << thisChromosome.constants_ << endl;
    }
    thisChromosome.print(cout);
    nchromosome++;
  }
}



bool GEP::reset()
{
  return true;
}




bool GEP::printValidation(unsigned cycle)
{
  // no print-out for zero training cycle
 
  if( cycle == 0 ) return true;

  double sigsig = 0;
  double sigbac = 0;
  double bacsig = 0;
  double bacbac = 0;

  double wtot = 0;

  // loop through validation data
  
  int vsize = ValidationDataset.size();
  
  for( int i=0;i<vsize;i++ ) {
    double w = 1.;

    double value = best_.Evaluate(ValidationDataset[i],0);
    wtot += w;
  }

  return true;
}


bool GEP::train(int verbose) {


  unsigned size = Dataset.size();
  assert(size > 0);
  unsigned dim = Dataset[0].size()-1;
  
  unsigned Epoch = 0; // counts the chromosome generations
  
  double best_fitness;
  int lastbest = 0;

 
  // main loop
  while( Epoch < Max_Epoch_ ) {
    
    Epoch++;

    // For each chromosome in the population we calculate its fitness
    // by calculating its specificity and sensitivity over the whole
    // data set
    
    best_fitness = GEP::Very_Large_Double;
    int best_chromosome = 0;
    
    double sumfitness = 0;
    
    for(int ic=0;ic<population_size_;ic++) {
		double sumdiff = 0;
		for(int i=0;i<size;i++) {
			double value = chromosomes_[ic].Evaluate(Dataset[i],verbose);
			double diff = (value-Dataset[i][dim])*(value-Dataset[i][dim]);
			sumdiff += diff;
		}
		sumfitness += sumdiff;
		if(sumdiff < best_fitness) {
			best_fitness = sumdiff;
			best_chromosome = ic;
		}

    }

    double avgfitness = sumfitness / population_size_;

    best_ = chromosomes_[best_chromosome];
	// If there is a perfect match, we break
	if(best_fitness == 0) break;

    if( lastbest != best_chromosome ) {
      lastbest = best_chromosome;	
      if( verbose > 0 ) {
	cout << endl << "Epoch " << Epoch 
	     << " Best chromosome=" << best_chromosome 
	     << " FOM=" << best_fitness << endl;
	best_.print(cout);
      }
    } 
    else if( verbose > 1 ) {
      cout << "Epoch " << Epoch << " Best/Average fitness " << best_fitness 
	   << "/" << avgfitness << endl;
    } 
    else {
      //if(Epoch % 10) cout << "." ;
    }



    if( valPrint_ != 0 && Epoch % valPrint_ == 0 ) {
      if( !this->printValidation(Epoch) ) {
	cerr << "Error printing validation data" << endl;
	return false;
      }
    }

    // Chromosome Mutation
    // Chromosome Transposition
    // a) by Root Insertion Sequence
    // b) by Insertion Sequence (not at the Root)
    // c) by Gene Transposition (not treated yet, and irrelevant for "+" linking function)
    // Chromosome Recombination
    // a) One Point
    // b) Two Point
    // c) Whole Gene
    // Constant mutation
    // Constant swap
		
    for( int ig=0;ig<population_size_;ig++ ) {
       // the best chromosome is always left alone
      if( ig == best_chromosome ) continue;

      Chromosome &chromosome = chromosomes_[ig];

      // Mutation
      if( randomOne() < mutation_rate_ ) {
	if( verbose > 2 ) { 
	  cout << "Before mutation " << endl; 
	  chromosome.print(cout);
	}
	bool mut = chromosome.mutate(FunctionsInUse, 
				     dim, constants_per_chromosome_);
	if(verbose > 2) { 
	  cout << "After mutation " << endl; 
	  chromosome.print(cout);}
      }

      // Root insertion
      if( randomOne() < RIS_rate_ ) {
	if( verbose > 2 ) { 
	  cout << "Before RIS " << endl; 
	  chromosome.print(cout);
	}
	bool ris = chromosome.RIS();
	if( verbose > 2 ) { 
	  cout << "After RIS " << endl; 
	  chromosome.print(cout);
	}
      }
			
      // Insertion
      if( randomOne() < IS_rate_ ) {
	if(verbose > 2) { 
	  cout << "Before IS " << endl; 
	  chromosome.print(cout);
	}
	bool is = chromosome.IS();
	if( verbose > 2 ) { 
	  cout << "After IS " << endl; 
	  chromosome.print(cout);
	}
      }

      // One point recombination
      if( randomOne() < One_Point_rate_) {
	// Pick another chromosome from the population, randomly
	// but not the best chromosome (which must survive unchanged)
	int iother = best_chromosome;
	while( iother == best_chromosome || iother == ig ) 
	  iother = rand() % population_size_;
	Chromosome &otherChromosome = chromosomes_[iother];
	if( verbose > 2 ) { 
	  cout << "Before One Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
	bool io = chromosome.OnePoint(otherChromosome);
	if( verbose > 2 ) { 
	  cout << "After One Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
      }

      // Two point recombination
      if( randomOne() < Two_Point_rate_ ) {
	// Pick another chromosome from the population, randomly
	// but not the best chromosome (which must survive unchanged)
	int iother = best_chromosome;
	while( iother == best_chromosome || iother == ig ) iother = rand() % population_size_;
	Chromosome &otherChromosome = chromosomes_[iother];
	if( verbose > 2 ) { 
	  cout << "Before Two Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
	bool it = chromosome.TwoPoint(otherChromosome);
	if( verbose > 2 ) { 
	  cout << "After Two Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
      }

      // Whole gene recombination
      if( randomOne() < Whole_Gene_rate_ ) {
	// Pick another chromosome from the population, randomly
	// but not the best chromosome (which must survive unchanged)
	int iother = best_chromosome;
	while( iother == best_chromosome || iother == ig ) iother = rand() % population_size_;
	Chromosome &otherChromosome = chromosomes_[iother];
	if( verbose > 2 ) { 
	  cout << "Before Whole Gene " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
	bool wg = chromosome.WholeGene(otherChromosome);
	if(verbose > 2) { 
	  cout << "After Whole Gene " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);}
      }

      // Constant mutation
      if( constants_per_chromosome_ > 0 ) {
	if( randomOne() < Constant_mutation_ ) {
	  if( verbose > 2 ) { 
	    cout << "Before Constant Mutation " << endl; 
	    chromosome.print(cout);
	  }
	  bool cm = chromosome.ConstantMutation(Constant_range_);
	  if( verbose > 2 ) { 
	    cout << "After Constant Mutation " << endl; 
	    chromosome.print(cout);}
	}
	if( randomOne() < Constant_swap_ ) {
	  // Pick another chromosome from the population, randomly
	  // but not the best chromosome (which must survive unchanged)
	  int iother = best_chromosome;
	  while( iother == best_chromosome || iother == ig ) iother = rand() % population_size_;
	  Chromosome &otherChromosome = chromosomes_[iother];
	  if( verbose > 2 ) { 
	    cout << "Before Constant Swap " << endl; 
	    chromosome.print(cout); 
	    otherChromosome.print(cout);
	  }
	  bool cs = chromosome.ConstantSwap(otherChromosome);
	  if( verbose > 2 ) { 
	    cout << "After Constant Swap " << endl; 
	    chromosome.print(cout); 
	    otherChromosome.print(cout);
	  }
	}
      }
    }
  }// while(Epoch < Max_Epoch_)

  cout << endl << "Finished training. Epoch " << Epoch 
       << " Best chromosome FOM " << best_fitness << endl;
  best_.print(cout);
  if( valPrint_ != 0 ) {
    if( !this->printValidation(Epoch) ) {
      cerr << "Error printing validation data" << endl;
      return false;
    }
  }

  cout << " C++ code that implements the trained chromosome: " << endl;
  best_.generateCode();

  return true;
}



char GEP::charFromInt(int i) 
{
  char c = '!';
  if(      i < Constant ) {
    c = charFromFunction(i);
  } 
  else if( i == Constant ) {
    c = '?';
  } 
  else if( i >= Variable ) {
    int ic = (i-Variable) % 26; // we map them all to a-z
    c = 'a' + ic;
  }
  return c;
}


void GEP::print(std::ostream& os) const
{
  os << "Trained GEP " << endl;
  best_.print(os);
}



