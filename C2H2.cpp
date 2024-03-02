
#include <iostream>
#include <armadillo>
#include "AO.h"
#include "hamiltonian.h"

using namespace std;


int main() {

    AO C2H2_ao("C2H2.txt");

    cout << "Overlap Matrix for C2H2: " << endl;
    vector<BasisFunction> basis_set = C2H2_ao.basis_set;

    arma::mat S = overlap_matrix(basis_set);
    S.print();

    cout << "Hamiltonian Matrix for C2H2: " << endl;
    arma::mat H = createhamiltonianEnergy(basis_set);
    H.print();


    cout << "X matrix for C2H2: " << endl;
    arma::mat X = createX(overlap_matrix(basis_set));
    X.print();

    cout << "Fancy H matrix for C2H2: " << endl;
    arma::mat fancy_H = fancyH(X, H);
    fancy_H.print();

    cout << "MO Coefficients C for C2H2: " << endl;
    MO_coefficients(X, fancy_H).print();

    cout << "Energy for C2H2: " << endl;
    double energy = calculateHamiltonianEnergy(C2H2_ao, S);
    cout << energy << endl;




    return 0;
}