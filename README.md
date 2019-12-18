# GeneExpress
A Machine Learning tool for fitting a functional form to a multivariate dataset.


Release Notes v 0.1
===================

The file "GeneExpress.cpp" contains an example application that uses the
GEP to perform function regression on a set of multivariate data points where
the functional form is known. Specifically, the example uses

f = 5*sin(a)*(1-b))

is used in the example. The data are thus triplets of (a,b,f). Fifty
triplets are generated using random values for (a,b) and the calculated f.

Of course, in general, the analytic function is not known.

Following the calculation of the example data set, the application will 
print the GEP parameters:

Begin GEP
 Genes per Chromosome         = 1
 Gene head length             = 7
 Constants per Chromosome     = 1
 Chromosome population size   = 200
 Extra functions requested    = AOSEL
 Maximum Epochs/Generations   = 2000
Using functions +-*/AOSEL
 Mutation Rate                = 0.05
 Root Insertion Rate          = 0.1
 Insertion Rate               = 0.1
 One Point Transposition Rate = 0.3
 Two Point Transposition Rate = 0.3
 Whole Gene Transpose Rate    = 0.1
 Constant Mutation Rate       = 0.1
 Constant Swap Rate           = 0.1
 Constant Generation Range    = 10
Creating Chromosome population of size 200
Each gene will have a tail of length 8 and a total length of 15

Then it will print each of the 200 randomly generated chromosomes, an example of
which is:

/?aaOA+|?abbb??a  L= + C= 3.94482 S= 0

There is one gene, which uses a single constant C (=3.94482). The head
of the gene (/?aaOA+) prescribes the following formula:

f = 3.94482 / a

i.e. the value of the chromosome will be the value of the first variable
in the data divided into 3.94482. See C.Ferreira's paper(s) for more details
on the formalism.

The linking function used is "+" (and is irrelevant since we use only one gene).

After printing the chromosome population, training will begin. The goal is
to minimise the FOM.







The GEP implements chromosomes/genes and functions as follows:

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


Julian Bunn, Caltech, 2019 
Julian.Bunn@caltech.edu
