#include <iostream>

#include "lml/matrix.hpp"

int main() {
  lml::Matrix<int> A(2, 2);

  A(0, 0) = 1, A(0, 1) = 2;
  A(1, 0) = 3, A(1, 1) = 4;

  std::cout << "A(0,1) = " << A(0, 1) << '\n';

  lml::Matrix<int> B = A + A;

  std::cout << "(A + A)(1,1) = " << B(1, 1) << '\n';

  lml::Matrix<int> C = 2 * A;

  std::cout << "(2 * A)(1,1) = " << C(1, 1) << '\n';

  lml::Matrix<int> D = A * A;

  std::cout << "(A * A)(1,1) = " << D(1, 1) << '\n';

  lml::Matrix<int> E = A.transpose();

  std::cout << "E(0,1) = " << E(0, 1) << " (should be 3) \n";

  lml::Matrix<int> zeros = lml::zeros<int>(2, 3);

  std::cout << "zeros =  " << zeros(0, 0) << zeros(0, 1) << zeros(0, 2)
            << zeros(1, 0) << zeros(1, 1) << zeros(1, 2) << '\n';

  lml::Matrix<int> identity = lml::identity<int>(3);

  std::cout << "identity" << identity(0, 0) << identity(0, 1) << identity(1, 0)
            << identity(1, 1) << identity(2, 2) << '\n';

  std::cout << "det(A) = " << A.det() << '\n';

  lml::Matrix<int> F(3, 3);

  F(0, 0) = 1, F(0, 1) = 2, F(0, 2) = 3;
  F(1, 0) = 4, F(1, 1) = 5, F(1, 2) = 6;
  F(2, 0) = 7, F(2, 1) = 8, F(2, 2) = 9;

  std::cout << "det(F) = " << F.det() << '\n';
  std::cout << "det(F) = " << F.detGaus() << '\n';

  lml::Matrix<int> G = A;
  G.rowswap(0, 1);

  std::cout << "G(0,0) = " << G(0, 0) << " should be 3.\n";

  lml::Matrix<int> H = A;
  H.rowmult(2, 0);

  std::cout << "H(0,0) = " << H(0, 0) << " should be 2.\n";

  std::cout << A(0, 0) << ' ' << A(0, 1) << '\n'
            << A(1, 0) << ' ' << A(1, 1) << '\n';

  lml::Matrix<int> I = A;
  I.rowreduce(0, 1, 0);

  std::cout << "I(0,0) = " << I(1, 0) << ' ' << I(1, 1) << '\n';

  return 0;
}
