# Optimized decision trees
This repository contains an advanced library to induce decision trees. It is capable of constructing (bivariate) oblique decision trees and cross-split decision trees. Furthermore, it handles nominal attributes naturally without treating them as numeric variables.
The Code is written in C++, yet it is intended to be used with Python. SWIG is used to build the extension module.

Note that this code has been developed as part of my PhD project and should merely be regarded as a prototype implementation that I never intended to publish. Not all features have been thoroughly tested and the code certainly contains many bugs that you are likely to encounter. Nonetheless, I still decided to publish it as I am convinced that someone out there might find it helpful. Have fun!

## Building for Linux
### Dependencies
The following dependencies need to be downloaded and installed if necessary:
* g++: Other compilers may work but have not been tested
* SWIG: https://www.swig.org/
* Boost: https://www.boost.org/
* Eigen3 (Version >= 3.4.0): https://eigen.tuxfamily.org/
* liblinear: https://www.csie.ntu.edu.tw/~cjlin/liblinear/
* CLP or Gurobi: https://github.com/coin-or/Clp or https://www.gurobi.com/
* CTPL: https://github.com/vit-vit/CTPL
* EXPRTK: http://www.partow.net/programming/exprtk/
* Eigenmvn: https://github.com/beniz/eigenmvn

In the makefile, the following variables need to be set to the folders containing these dependencies:  

    BOOST_DIR=PATH/TO/BOOST
    EIGEN_DIR=PATH/TO/EIGEN3
    LIBLINEAR_DIR=PATH/TO/LIBLINEAR
    CLP_DIR=PATH/TO/CLP
    CTPL_DIR=PATH/TO/CTPL
    EXPRTK_DIR=PATH/TO/EXPRTK
    EIGENMVN_DIR=PATH/TO/EIGENMVN
    
If you want to use Gurobi, you need to set the following variables (make sure to link to the right libraries depending on your system):

    USEGUROBI=-Duse_gurobi
    GUROBI_INC=-IPATH/TO/GUROBI/include/
    GUROBI_LIB=-LPATH/TO/GUROBI/lib -lgurobi_g++X.X -lgurobiXX

### Installation
Use `make` to compile the library. This will create the 'lib' folder which can be included in your python scripts to use the library.
    
## Usage
* `DecisionTree()`: Intializes the decision tree
* `setParam(param,value)`: Used to set the parameter `param` to the value `v`. `param` and `v` can also be arrays to set multiple parameters simultaneously in one function call.
  * `param`: str or array of str
  * `v`: float/str or array of float/str
* `fit(X,y,nominal_indices)`: Fits a classification tree for the specified trainig data. It is assumed that the labels are integer values starting at zero.
  * `X`: 2d numerical matrix
  * `y`: array of int
  * `nominal_indices` (optional): array of column indces which specify the nominal attributes. Can be omitted
* `fit_reg(X,y,nominal_indices)`: Fits a regression tree for specified trainig data
  * `X`: 2d numerical matrix
  * `y`: array of float
  * `nominal_indices` (optional): array of column indces which specify the nominal attributes.
* `predict(X)`: Used to predict the targets of the input matrix `X` for classification
  * `X`: 2d numerical matrix
* `predict_reg(X)`: Used to predict the targets of the input matrix `X` for regression
  * `X`: 2d numerical matrix
* `prune(X,y)`: Performs reduced error pruning in the classification setting with the pruning set specified by `X` and `y`
  * `X`: 2d numerical matrix
  * `y`: array of numerical labels
* `prune_reg(X,y)`: Performs reduced error pruning in the regression setting with the pruning set specified by `X` and `y`
  * `X`: 2d numerical matrix
  * `y`: array of numerical labels
* `MCCprune(X,y)`: Performs minimum cost-complexity pruning in the classification setting with the pruning set specified by `X` and `y`
  * `X`: 2d numerical matrix
  * `y`: array of numerical labels
* `MCCprune_reg(X,y)`: Performs minimum cost-complexity pruning in the regression setting with the pruning set specified by `X` and `y`
  * `X`: 2d numerical matrix
  * `y`: array of numerical labels
* `export_graphviz(on_edge)`: Exports a fitted decision tree to a str in graphviz format (see, https://graphviz.org/). The parameter `on_edge` specifies whether the rules are written to the edges.
  * `on_edge`: Boolean
* `setAttributeNames(names)` (only for visualization): Sets the names of the attributes 
  * `names`: array of str
* `setYNames(values,categories)` (only for visualization): Assigns str representations to the numerical values in the target vector of the training set
  * `values`: array of int 
  * `categories`: array of str
* `setNominalNames(index,values,categories)` (only for visualization): Assigns str representations to the numerical values in the column of a nominal feature specified by `index`
  * `index`: int
  * `values`: array of int
  * `names`: array of str
* `nodeCount()`: Returns the number of nodes in the tree
* `leafCount()`: Returns the number of leaf nodes
* `getDepth()`: Returns the depth of the tree
* `save(filename)` (experimental): Save a decision tree model to a file
  * `filename`: str
* `load(filename)` (experimental): Load a decision tree model from a file
  * `filename`: str

## Parameters
* **criterion**: {"gini", "entropy", "mse", "mae"}  
Specifies the splitting criterion. Gini and entropy are possible for classification and for regression, either mean squared error or mean absolute error.  
Default: "gini"
* **1d**: {0, 1}  
Specifies whether univariate splits are calculated  
Default: 1
* **2d**: {0, 1}  
Specifies whether the simulated annealing algorithm should be used to determine oblique splits.  
Default: 1
* **sa**: {0, 1}  
Specifies whether the simulated annealing algorithm should be used to determine oblique splits.  
Default: 0
* **ce**: {0, 1}  
Specifies whether the cross-entropy method should be used to determine oblique splits.  
Default: 0
* **xsplit**: {0, 1}  
Specifies whether cross-splits should be used.  
Default: 0
* **lookahead**: {0, 1}  
Specifies whether the lookahead approach should be taken if cross-splits are used.  
Default: 1
* **max_margin**: {0, 1, 2}  
Specifies whether post-processing of oblique splits is carried out. If 0, no post-processing is carried out, otherwise either the l1 or l2 norm is used to solve the underlying generalized svm problem.  
Default: 0
* **max_depth**: int   
Specifies the maximum depth of the decision tree.  
Default: inf
* **min_samples_split**: int  
Specifies the minimum number of observations at a node to perform a split.  
Default: 2
* **min_cases**: int  
Specifies the minimum number of required observations at a node.  
Default: 1
* **max_features**: int  
Maximum number of considered attributes per split. If smaller than the number of features in the dataset, the features are chosen at random.  
Default: inf
* **max_leaf_nodes**: int  
Maximum number of leaf nodes.  
Default: inf
* **min_impurity_decrease**: float  
Minimum required decrease in impurity to perform a split.  
Default: -1
* **min_impurity_split**: float  
Minimum required impurity to perform a split.  
Default: -1
* **random_state**: int  
Sets the seed for the involved randomized processes.  
Default: Current time
* **nominal_timelimit**: int  
Specifies the timelimit in seconds for the branch and bound algorithm to determine nominal splits.  
Default: inf
* **n_threads**: int  
Specifies the number of threads if parallelization is applicable. A value of zero indicates that no parallelization is carried out.  
Default: 0
* **oblique_split_threshold**: int  
Specifies the minimum number of observations required to introduce an oblique split.  
Default: 1
* **random_splits**: {0, 1}  
Splits are calculated completely at random.  
Default: 0
* **normalize**: {0, 1}  
Specifices whether feature scaling should be applied before determining the splits.  
Default: 1
* **nominal_partitions**: int 
Specifices the arity of the nominal partitions.  
Default: 2
* **ce_samples**: int or str
Specifices the number of samples drawn by the cross-entropy method. Can either be an integer or a string specifying a formula that depends on the number of observations "n" at the respective node and the number of numerical features "d".  
Default:  "round(2 * d * log2(n))"
* **cross_entropy_no_improvement**: int or str
Specifices the maximum number of iterations of the cross-entropy method without improvement of the level parameter gamma. Can either be an integer or a string specifying a formula that depends on the number of observations "n" at the respective node and the number of numerical features "d".  
Default: 3
* **cross_entropy_alpha**: [0, 1)
Specifices the update parameter alpha of the cross-entropy method.  
Default: 0.2
* **cross_entropy_rho**: (0, 1)
Specifices the quantil parameter rho of the cross-entropy method.  
Default: 0.1
* **a_start_iterations**: int or str
Specifices the number of start iterations of the simulated annealing method to calculate the initial temperature. Can either be an integer or a string specifying a formula that depends on the number of observations "n" at the respective node and the number of numerical features "d".  
Default: 100
* **a_no_improvement**: int or str
Specifices the maximum number of iterations of the simulated annealing method without improvement of the optimal solution. Can either be an integer or a string specifying a formula that depends on the number of observations "n" at the respective node and the number of numerical features "d".  
Default: "max(d,100)"
* **sa_beta**: (0, 1)
Specifices the beta parameter for the geometric cooling schedule of the simulated annealing method.  
Default: 0.85
* **sa_scale**: {no,geom,equ}
Specifices whether advanced scaling of the input matrix should be carried out for the simulated annealing method. One can choose between no schaling, geometric mean scaling or equilibration
Default: geom
* **sa_max_eta**: int
Maximum number of involved Eta matrices (resulting from the pivot operations) before a new LU factorization is computed.  
Default: 100
