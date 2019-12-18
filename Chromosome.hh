//
//
// Description:

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
 
#ifndef _Chromosome_HH
#define _Chromosome_HH

#include <string>
#include <iostream>
#include <vector>


class Gene 
{
public:
  virtual ~Gene() {}

  Gene() 
    : head_(), tail_() {}

  Gene(const Gene& other)
    :
    head_(other.head_),
    tail_(other.tail_)
  {}

  Gene& operator=(const Gene& other) {
    head_ = other.head_;
    tail_ = other.tail_;
    return *this;
  }

  void setHead(const std::vector<int>& head) { head_ = head; }
  void setTail(const std::vector<int>& tail) { tail_ = tail; }

  std::vector<int> getHead() { return head_; }
  std::vector<int> getTail() { return tail_; }

  void generateTree(std::vector<std::vector<int> >& tree) const;
  
  void print(std::ostream& os) const;

private:
  std::vector<int> head_;
  std::vector<int> tail_;
};


class Chromosome 
{
public:
  unsigned num_genes_;
  char link_function_;
  std::vector<Gene> gene_;
  std::vector<double> constants_;
  double fitness_;
  unsigned length_;
  double shift_;
		
  virtual ~Chromosome();

  Chromosome();
  Chromosome(
		unsigned num_genes,
		char link_function);
  
  double Evaluate(const std::vector<double> &data, int verbose) const;
  void generateCode();
  
  bool mutate(const std::vector<int>& funcs, 
	      unsigned dim, 
	      unsigned constants_per_chromosome);
  bool RIS();
  bool IS();
  bool OnePoint(Chromosome &other);
  bool TwoPoint(Chromosome &other);
  bool WholeGene(Chromosome &other);
  bool ConstantMutation(double range);
  bool ConstantSwap(Chromosome &other);

  void setLength(unsigned num_genes, unsigned num_constants, unsigned length) {
    num_genes_ = num_genes;
    gene_.resize(num_genes);
    constants_.resize(num_constants);
    length_ = length;
  }

  void print(std::ostream& os) const;

private:
};

#endif
