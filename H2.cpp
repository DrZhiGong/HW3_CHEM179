
#include <iostream>
#include <armadillo>
#include "AO.h"
#include "hamiltonian.h"

using namespace std;


int main() {

    AO H2_ao("H2.txt");

    cout << "Overlap Matrix for H2: " << endl;
    vector<BasisFunction> basis_set = H2_ao.basis_set;

    arma::mat S = overlap_matrix(basis_set);
    S.print();

    cout << "Hamiltonian Matrix for H2: " << endl;
    arma::mat H = createhamiltonianEnergy(basis_set);
    H.print();


    cout << "X matrix for H2: " << endl;
    arma::mat X = createX(overlap_matrix(basis_set));
    X.print();

    cout << "Fancy H matrix for H2: " << endl;
    arma::mat fancy_H = fancyH(X, H);
    fancy_H.print();

    cout << "MO Coefficients C for H2: " << endl;
    MO_coefficients(X, fancy_H).print();

    cout << "Energy for H2: " << endl;
    double energy = calculateHamiltonianEnergy(H2_ao, S);
    cout << energy << endl;

    return 0;
}