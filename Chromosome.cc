//$Id: Chromosome.cc,v 1.1 2008/02/22 23:07:50 narsky Exp $

#include "GEP.hh"
#include "Chromosome.hh"

#include <stdio.h>
#include <cassert>
#include <utility>
#include <algorithm>
#include <functional>
#include <cmath>
#include <sstream>

using namespace std;


Chromosome::~Chromosome() 
{

}


Chromosome::Chromosome() 
  :
  num_genes_(1),
  link_function_('+'),
  gene_(1),
  constants_(),
  fitness_(0),
  length_(0),
  shift_(0)
{}




Chromosome::Chromosome(
			     unsigned num_genes,
			     char link_function)
  :
  num_genes_(num_genes),
  link_function_(link_function),
  gene_(1),
  constants_(),
  fitness_(0),
  length_(0),
  shift_(0)
{
  assert( num_genes >= 1 );
  gene_.resize(num_genes);
}


void Chromosome::generateCode() 
{
  ostringstream outs;
  int iconstant = 0;

  cout << "double Signal(double x[]) {" << endl;

  if( link_function_ == '/' || link_function_ == '*' ) {
    cout << "     double Chromosome_value = 1;" << endl;
  } 
  else {
    cout << "     double Chromosome_value = 0;" << endl;
  }

  for( int ch=0; ch<num_genes_; ch++ ) {

    Gene &gene = gene_[ch];
    
    vector<vector<int> > tree;
    gene.generateTree(tree);
    
    int level = tree.size()-1;
    
    vector<string> vrow, vrowbelow;

    string value;
    double val;

    for( int i=level;i>=0;i-- ) {
      int iposbelow = 0;
      for( int ic=tree[i].size()-1;ic>=0;ic-- ) {
	int c = tree[i][ic];
	if(      c == GEP::Constant ) {
	  // pick the next constant from the chromosome
	  val = constants_.at(iconstant++);
	  if( iconstant >= constants_.size() ) iconstant = 0;
	  outs << val;
	  vrow.push_back(outs.str());
	  outs.str("");
	} 
	else if( c >= GEP::Variable ) {
	  int index = c - GEP::Variable;
	  outs << index;
	  value = "x[" + outs.str() + "]";
	  outs.str("");
	  vrow.push_back(value);
	} 
	else {
	  switch(c) 
	    {
	    case GEP::Sqrt:
	      value = "sqrt(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case GEP::Exp:
	      value = "exp(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case GEP::Sin:
	      value = "sin(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case GEP::Abs:
	      value = "abs(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case GEP::Log:
	      value = "log(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case GEP::Times:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "*" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case GEP::Divide:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "/" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case GEP::Plus:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "+" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case GEP::Minus:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "-" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case GEP::LessThan:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "<" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case GEP::GreaterThan:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ ">" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    default:
	      cerr << "Error parsing value tree" << endl;
	    }
	}
      }
      vrowbelow = vrow;
      vrow.clear();
    }
    assert (vrowbelow.size() == 1);
    
    cout << "     Chromosome_value " << outs.str() << link_function_ 
	 << "= " << vrowbelow[0] << ";" << endl;

  }

  cout << "     Chromosome_value += " << -shift_ << ";" << endl;

  cout << "     return Chromosome_value;" << endl;
  cout << "}" << endl;
}
	

double Chromosome::Evaluate(const std::vector<double>& data, 
			       int verbose) const
{
  // Parse the gene(s) and evaluate each with the given data
  vector<double> gene_value;
  
  int iconstant = 0;
  
  for(int ch=0; ch<num_genes_; ch++) {
    
    const Gene &gene = gene_[ch];
    
    vector<vector<int> > tree;
    gene.generateTree(tree);
    
    int level = tree.size();
    
    if( verbose > 1 ) {
      cout << "Evaluating Gene #" <<  ch << ":";
      gene.print(cout);
      cout << endl;
      //if(constants_.size() > 0) cout << "  Constants " << constants_ << endl;
    }
    
    if( verbose > 1 ) {
      cout << "Finished parsing Gene" << endl;
      for( int i=0;i<level;i++ ) {
	cout << "Level " << i << ": ";
	for( int it=0;it<tree[i].size();it++ ) 
	  cout << GEP::charFromInt(tree[i][it]);
	cout << endl;
      }
    }

    level--;

    vector<double> vrow;
    vector<double> vbelow;
    double val;
    for( int i=level;i>=0;i-- ) {
      int iposbelow = 0;
      for( int ic=tree[i].size()-1;ic>=0;ic-- ) {
	int c = tree[i][ic];
	if(     c == GEP::Constant ) {
	  // pick the next constant from the chromosome
	  val = constants_.at(iconstant++);
	  if( iconstant >= constants_.size() ) iconstant = 0;
	  vrow.push_back(val);	
	} 
	else if( c >= GEP::Variable ) {
	  int index = c - GEP::Variable;
	  val = data.at(index);
	  vrow.push_back(val);
	} 
	else {
	  switch(c) {
	  case GEP::Sqrt:
	    // need to protect against sqrt neg number
	    if( vbelow.at(iposbelow) < 0 ) {
	      val = 0;
	    } 
	    else {
	      val = sqrt(vbelow.at(iposbelow));
	    }
	    iposbelow++;
	    break;
	  case GEP::Abs:
	    // absolute value
	    val = abs(vbelow.at(iposbelow));
	    iposbelow++;
	    break;
	  case GEP::Log:
	    // log base e
	    // need to protect against log of 0 or -ve
		  if(vbelow.at(iposbelow) <= 0) {
			  val = GEP::Very_Small_Double;
		  } else {
			  val = log(vbelow.at(iposbelow));
		  }
	    iposbelow++;
	    break;
	  case GEP::Exp:
	    // exponential
	    val = exp(vbelow.at(iposbelow));
		if(val > GEP::Very_Large_Double) val = GEP::Very_Large_Double;
	    iposbelow++;
	    break;
	  case GEP::Sin:
	    // sin
	    val = sin(vbelow.at(iposbelow));
	    iposbelow++;
	    break;
	  case GEP::LessThan:
	    // Less than
	    val = 0;
	    if(vbelow.at(iposbelow+1) < vbelow.at(iposbelow)) val = 1;
	    iposbelow += 2;
	    break;
	  case GEP::GreaterThan:
	    // Greater than
	    val = 0;
	    if(vbelow.at(iposbelow+1) > vbelow.at(iposbelow)) val = 1;
	    iposbelow += 2;
	    break;
	  case GEP::Times:
	    val = vbelow.at(iposbelow+1) * vbelow.at(iposbelow);
	    iposbelow += 2;
	    break;
	  case GEP::Divide:
	    // need to protect against divide by zero
	    if( vbelow.at(iposbelow) == 0 ) {
	      val = GEP::Very_Large_Double;
	    } 
	    else {
	      val = vbelow.at(iposbelow+1) / vbelow.at(iposbelow);
	    }
		assert(val == val);
	    iposbelow += 2;
	    break;
	  case GEP::Plus:
	    val = vbelow.at(iposbelow+1) + vbelow.at(iposbelow);
	    iposbelow += 2;
	    break;
	  case GEP::Minus:
	    val = vbelow.at(iposbelow+1) - vbelow.at(iposbelow);
	    iposbelow += 2;
	    break;
	  default:
	    cerr << "Error parsing value tree" << endl;
	    return 0;
	  }

	  // general check on value
	  assert(val == val); // for NaN

	  vrow.push_back(val);

	  if(verbose > 2) cout << "Value at level " << i << " is " << val << endl;
	}
      }
      vbelow = vrow;
      vrow.clear();
    }
    assert (vbelow.size() == 1);
    gene_value.push_back(vbelow.at(0));
    if(verbose > 1) cout << "Evaluates to: " << vbelow.at(0) << endl << endl;
  }
  
  // Calculate the return value of the chromosome. 
  // This depends on the linking function and
  // the number of genes. If there is only one gene, then the linking function
  // is not used, and the chromosome has the same value as the gene.
  
  double result = gene_value.at(0);

  if(num_genes_ == 1) return result;
	
  for(int ch=1;ch<num_genes_;ch++) {
    if(      link_function_ == '*' ) {
      result *= gene_value.at(ch);
    } 
    else if( link_function_ == '/' ) {
      result /= gene_value.at(ch);
    } 
    else if( link_function_ == '+' ) {
      result += gene_value.at(ch);
    } 
    else if( link_function_ == '-' ) {
      result -= gene_value.at(ch);
    }
  }

  // shift
  result -= shift_;

  if(verbose > 1) cout << "Chromosome evaluates to: " << result << endl;
  
  return result;
}


void Gene::print(std::ostream& os) const 
{
  vector<int> sequence(head_);
  int lhead = sequence.size();
  sequence.insert(sequence.end(),tail_.begin(),tail_.end());

  for(int is=0;is<sequence.size();is++) {
    if(is == lhead) os << "|";
    int val = sequence.at(is);
    os << GEP::charFromInt(val);
  }
}


void Chromosome::print(std::ostream& os) const 
{
  for(int ch=0;ch<num_genes_;ch++) {
    gene_.at(ch).print(os);
    os << " ";
  }
  os << " L= " << link_function_;
  if( !constants_.empty() ) os << " C=";
  for( int i=0;i<constants_.size();i++ )
    os << " " << constants_[i];
  os << " S= " << shift_;
  os << endl;
}


bool Chromosome::ConstantMutation(double range) 
{
  // pick a constant at random
  int ic = rand() % constants_.size();
  // replace it with a new one
  constants_[ic] = randomOne() * range;
  return true;
}

bool Chromosome::ConstantSwap(Chromosome &other) 
{
  // choose random constants from each chromosome and swap them
  int ic1 = rand() % constants_.size(); 
  int ic2 = rand() % constants_.size(); 
  swap(constants_[ic1],other.constants_[ic2]);
  /*
    double val = constants_[ic1];
    constants_[ic1] = other.constants_[ic2];
    other.constants_[ic2] = val;
  */
  return true;
}

bool Chromosome::WholeGene(Chromosome &other) 
{
  // Choose a random gene from both chromosomes, and swap them
  int ig1 = rand() % num_genes_;
  int ig2 = rand() % num_genes_;
  swap(gene_[ig1],other.gene_[ig2]);
  /*
    vector<int> head1 = gene_[ig1].getHead();
    vector<int> tail1 = gene_[ig1].getTail();
    gene_[ig1].setHead(other.gene_[ig2].getHead());
    gene_[ig1].setTail(other.gene_[ig2].getTail());
    other.gene_[ig2].setHead(head1);
    other.gene_[ig2].setTail(tail1);
  */
  return true;
}


bool Chromosome::OnePoint(Chromosome &other) 
{
  // We choose a random point in the chromosome, and exchange material
  // The assumptions is that chromosomes are the same size as are their genes
  vector<int> chromosome1 = gene_[0].getHead();
  int lhead = chromosome1.size();
  vector<int> tail1 = gene_[0].getTail();
  int ltail = tail1.size();
  chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
  vector<int> chromosome2 = other.gene_[0].getHead();
  vector<int> tail2 = other.gene_[0].getTail();
  chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  
  for(int i=1;i<num_genes_;i++) {
    vector<int> head1 = gene_[i].getHead();
    tail1 = gene_[i].getTail();
    vector<int> head2 = other.gene_[i].getHead();
    tail2 = other.gene_[i].getTail();
    chromosome1.insert(chromosome1.end(),head1.begin(),head1.end());
    chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
    chromosome2.insert(chromosome2.end(),head2.begin(),head2.end());
    chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  }
  assert(chromosome1.size() == chromosome2.size());
  int len = chromosome1.size();
  int pos = rand() % len;
  
  // Swap the contents from position pos to the end of the chromosomes
  swap_ranges(chromosome1.begin()+pos, 
	      chromosome1.end(),
	      chromosome2.begin()+pos);
  
  pos = 0;
  int ig = 0;
  while( pos < len && ig < num_genes_ ) {
    vector<int> new1(chromosome1.begin()+pos,chromosome1.begin()+pos+lhead);
    vector<int> new2(chromosome2.begin()+pos,chromosome2.begin()+pos+lhead);
    gene_[ig].setHead(new1);
    other.gene_[ig].setHead(new2);
    pos += lhead;
    vector<int> new3(chromosome1.begin()+pos,chromosome1.begin()+pos+ltail);
    vector<int> new4(chromosome2.begin()+pos,chromosome2.begin()+pos+ltail);
    gene_[ig].setTail(new3);
    other.gene_[ig].setTail(new4);
    pos += ltail;
    ig++;
  }
  
  return true;
}

bool Chromosome::TwoPoint(Chromosome &other) 
{
  // We choose two random points in the chromosome, and exchange material
  // The assumptions is that chromosomes are the same size as are their genes
  vector<int> chromosome1 = gene_[0].getHead();
  int lhead = chromosome1.size();
  vector<int> tail1 = gene_[0].getTail();
  int ltail = tail1.size();
  chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
  vector<int> chromosome2 = other.gene_[0].getHead();
  vector<int> tail2 = other.gene_[0].getTail();
  chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  
  for(int i=1;i<num_genes_;i++) {
    vector<int> head1 = gene_[i].getHead();
    tail1 = gene_[i].getTail();
    vector<int> head2 = other.gene_[i].getHead();
    tail2 = other.gene_[i].getTail();
    chromosome1.insert(chromosome1.end(),head1.begin(),head1.end());
    chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
    chromosome2.insert(chromosome2.end(),head2.begin(),head2.end());
    chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  }
  assert(chromosome1.size() == chromosome2.size());
  int len = chromosome1.size();
  
  int pos1 = rand() % (len-1);
  int pos2 = pos1 + rand() % (len-pos1);
  
  assert(pos1 <= pos2);
  int nchars = pos2-pos1+1;
  
  // Swap the contents from position pos1 to position pos2
  swap_ranges(chromosome1.begin()+pos1, chromosome1.begin()+pos2,
	      chromosome2.begin()+pos1);
  
  int pos = 0;
  int ig = 0;
  while(pos < len && ig < num_genes_) {
    vector<int> new1(chromosome1.begin()+pos,chromosome1.begin()+pos+lhead);
    vector<int> new2(chromosome2.begin()+pos,chromosome2.begin()+pos+lhead);
    gene_[ig].setHead(new1);
    other.gene_[ig].setHead(new2);
    pos += lhead;
    vector<int> new3(chromosome1.begin()+pos,chromosome1.begin()+pos+ltail);
    vector<int> new4(chromosome2.begin()+pos,chromosome2.begin()+pos+ltail);
    gene_[ig].setTail(new3);
    other.gene_[ig].setTail(new4);
    pos += ltail;
    ig++;
  }
  
  return true;
}

bool Chromosome::mutate(const std::vector<int>& funcs, 
			   unsigned dim, unsigned nconsts) 
{

  // select a random gene
  int ic = rand() % num_genes_;
  vector<int> head = gene_[ic].getHead(); 
  vector<int> tail = gene_[ic].getTail();
  int lhead = head.size();
  int ltail = tail.size();
  // select a random position in the gene
  int posinc = rand() % (lhead+ltail);
  // in the head we can change to a function or a terminal
  // in all cases, we force a real mutation, 
  // i.e. the terminal/function must change
  if( posinc < lhead ) {
    // note that the first position in the head must be a function
    int headval = head.at(posinc);
    if( posinc == 0 || randomOne() < 0.5 ) {
      // mutate to a (different) function
      int ifunc = funcs.at(rand() % funcs.size());
      while(ifunc == headval) ifunc = funcs.at(rand() % funcs.size());
      head[posinc] = ifunc;
    } 
    else {
      int i = rand() % (dim+nconsts);
      int cterm = GEP::Constant;
      if(i < dim) {
	cterm = GEP::Variable + i;	
	while(cterm == headval) 
	  cterm = GEP::Variable + (rand() % dim);
      }
      head[posinc] = cterm;
    }
  } 
  else {
    // in the tail we can only mutate to a terminal or constant
    int tailval = tail[posinc-lhead];
    int i = rand() % (dim+nconsts);
    int cterm = GEP::Constant;
    if(i < dim) {
      cterm = GEP::Variable + (rand() % dim);
      while(cterm == tailval) cterm = GEP::Variable + (rand() % dim);
    }
    tail[posinc-lhead] = cterm;
  }
  gene_[ic].setHead(head);
  gene_[ic].setTail(tail);
  return true;
}


bool Chromosome::RIS() 
{
  // select a random gene
  int ic = rand() % num_genes_;
  vector<int> head = gene_[ic].getHead(); 
  int lhead = head.size();
  // select a random position in the head, but not the start
  int poshead = rand() % lhead;
  while( poshead == 0 ) poshead = rand() % lhead;
  // find first downstream function
  for( int i=poshead;i<head.size();i++ ) {
    if(head[i] < GEP::Constant) { // in the function space
      int len = lhead-i;
      vector<int> ris(head.begin()+i,head.end());
      // move this function headed sequence to the head of the gene
      head.insert(head.begin(),ris.begin(),ris.end());
      // discard the parts of the head that were shifted out
      head.erase(head.begin()+lhead,head.end());
      gene_[ic].setHead(head);
      return true;
    }
  }
  return false;
}

bool Chromosome::IS() 
{
  // select a random gene
  int ic = rand() % num_genes_;
  vector<int> gene(gene_[ic].getHead());
  int lhead = gene.size();
  vector<int> tail(gene_[ic].getTail());
  gene.insert(gene.end(),tail.begin(),tail.end());
  int len = gene.size();
  // select a random position in the gene
  int pos = rand() % (len-1);
  // select a random length for the transposon
  int lseq = 1 + (rand() % (len-pos-1));
  vector<int> transposon(gene.begin()+pos,gene.begin()+pos+lseq);
  // select a random position in the chromosome
  // first select a random gene
  ic = rand() % num_genes_;
  vector<int> head = gene_[ic].getHead();
  pos = rand() % lhead;
  // ensure the position is not at the first position in the head
  while( pos == 0 ) pos = rand() % lhead;
  head.insert(head.begin()+pos,transposon.begin(),transposon.end());
  head.erase(head.begin()+lhead,head.end());
  gene_[ic].setHead(head);
  return true;
}


void Gene::generateTree(std::vector<std::vector<int> >& tree) const 
{
  vector<int> geneVector = head_;
  geneVector.insert(geneVector.end(),tail_.begin(),tail_.end());

  tree.clear();
  tree.resize(GEP::Max_Depth);

  vector<int> nargsreq(GEP::Max_Depth,0);
  int ipos = 0;
  int level = 0;
  nargsreq[level] = 1;
  while( ipos < geneVector.size() ) {
    int c = geneVector[ipos];
    tree[level].push_back(c);
    nargsreq[level]--;
    if(c == GEP::Constant || c >= GEP::Variable) { 
    } else {
      switch(c) {
	case GEP::Sqrt:
	case GEP::Exp:
	case GEP::Abs:
	case GEP::Sin:
	case GEP::Log:
	  nargsreq[level+1] += 1;
	  break;
	case GEP::Plus:
	case GEP::Minus:
	case GEP::Times:
	case GEP::Divide:
	case GEP::LessThan:
	case GEP::GreaterThan:
	  nargsreq[level+1] += 2;
	  break;
	default:
	  cerr << "Unrecognised function in gene \"" << c << "\"" << endl;
	  return;
	}
    }
    if(nargsreq[level] == 0) {
      if(nargsreq[level+1] == 0) break;
      level++;
    }
    ipos++;
    if(level+1 >= GEP::Max_Depth) {
      cerr << "Max expression tree depth " << GEP::Max_Depth 
	   << " reached for C=";
      this->print(cout);
      cout << endl;
      return;
    }
  }
  tree.resize(level+1);

}
