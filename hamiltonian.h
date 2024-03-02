
const double K = 1.75;

arma::mat createhamiltonianEnergy(vector<BasisFunction>& basis_set);
arma::mat fancyH(arma::mat X, arma::mat H);
arma::mat createX(arma::mat overlap_matrix);
arma::mat MO_coefficients(arma::mat X, arma::mat Fancy_H);
double calculateEnergy(arma::mat X, arma::mat Fancy_H, int num_electrons);
double calculateHamiltonianEnergy(AO AO_object, arma::mat Overlap_matrix);

//This is about the Hamiltonian operator calculation process through Armadillo,
//this is also the first time to write the main program and the header file separately, 
//so there may be a slight problem, kindly bear with me!