#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>


using namespace std;

struct BasisFunction {
    string AO_type;
    arma::vec exponents;
    arma::vec contracted_coefficients;
    arma::rowvec center;
    arma::vec lmn;
    arma::vec normalization_constants;
};

class AO {

private:
    int num_hydrogen;
    int num_carbon;
    int num_basis;
    int natoms;
    arma::mat coord;
    vector<int> atomic_numbers;

    int count_basis(int num_carbon, int num_hydrogen);
    int count_electrons(int num_carbon, int num_hydrogen);
    BasisFunction createHydrogenBasis(arma::rowvec center); 
    vector<BasisFunction> createCarbonBasis(arma::rowvec center); 
    void updateNormalizationConstants(BasisFunction& basisFunction);
    void updateBasisSet(vector<BasisFunction>& basis_set);



public:
  
    AO(const char* filename);
    ~AO();
    int num_electrons;

    vector<BasisFunction> basis_set;



};

double overlap_integral_1D(double center1, double center2, double alpha1, double alpha2, double l1, double l2);
double overlap_integral_3D(arma::rowvec center1, arma::rowvec center2, double alpha1, double alpha2, arma::vec lmn1, arma::vec lmn2);
double unnormOverlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2, int i, int j);
double overlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2);
arma::mat overlap_matrix(vector<BasisFunction>& basis_set);

