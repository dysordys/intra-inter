# intra-inter
Manuscript and all the supporting material for "The effect of intra- and interspecific competition on coexistence in multispecies communities"






## Code
The full code necessary to search for the most/least stabilized community, given a set of matrix elements, is contained in `Optimize_leading_eigenvalue.c`.
This code can be compiled using the command

    gcc Optimize_leading_eigenvalue.c -Wall -O3 -DHAVE_INLINE -lgsl -lgslcblas -o [Out File]

where `[Out File]` is replaced with your desired program name.

The compiled program takes eight inputs:

| Parameter | Description | Value |
|---|---|---|
| `n` | the size of the (square) matrix | integer |
| `filename` | the filename of the matrix to be optimized (stored as a space-delimited text file) | string |
| `seed` | a seed for the random number generator | integer |
| `maximize` | indicates whether to maximize or minimize the leading eigenvalue | 1 or -1 |
| `para1` | first parameter for search algorithm | integer |
| `para2` | second parameter for search algorithm | integer |
| `BvsC` | indicates which matrix elements can be mutated (diagonal or off-diagonal) | 0 or 1 |
| `SearchAlg` | indicates which algorithm to use for the optimization (Genetic algorithm or Hill climber) | 1 or 2 |

*Note: all "or" statements in the Value column are respective to the order indicated by the Description column, and all "integer"s are strictly positive*

For the Genetic Algoritm, the two parameters are the number of individuals in the population and the number of generations to run the algorithm for, respectively, while for the Hill Climber, they are the number of steps to take without finding a better solution before quitting and the number of mutations to try at each step.

Though the program does not return anything to the command line (other than status messages), it will produce a new matrix file in the same directory as the input matrix with ".Bmin", ".Bmax", ".min", or ".max" appended to the input filename depending on the values of `maximize` and `BvsC`:

| `BvsC` | `maximize` | Output Filename |
|---|---|---|
| 1 | 1 | `[filename].Bmax` |
| 1 | -1 | `[filename].Bmin` |
| 0 | 1 | `[filename].max` |
| 0 | -1 | `[filename].min` |

where `[filename]` is the provided input filename.

This new file will have exactly the same matrix elements, but rearranged according to the provided constraints to either minimize or maximize the leading eigenvalue.
The output matrix will not be arranged in any systematic way, but the patterns observed in the manuscript can be elicited by sorting the rows/columns according to the eigenvctor associated with either the leading or trailing eigenvalue.
