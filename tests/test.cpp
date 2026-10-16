#include <iostream>
#include <vector>

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

  lml::Matrix<int> E = A;
  E = E.transpose();

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

  std::cout << "I(1) = " << I(1, 0) << ' ' << I(1, 1) << '\n';

  lml::Matrix<int> J = A;
  J.rowadd(0, 1, 2);

  std::cout << "J(1) = " << J(1, 0) << ' ' << J(1, 1) << " (should be 5 8)\n";

  lml::Matrix<double> K1{{1, 2}, {3, 4}};

  std::cout << "K1 = " << K1(0, 0) << K1(0, 1) << K1(1, 0) << K1(1, 1)
            << " (should be 1234)\n";

  lml::Matrix<double> K(2, 2);
  K(0, 0) = 2, K(0, 1) = 1;
  K(1, 0) = 1, K(1, 1) = 3;

  std::vector<double> b = {5, 6};

  std::vector<double> res = lml::solveAxb(K, b);

  std::cout << "res = {" << res[0] << ", " << res[1]
            << "}(should be {1.8,1.4})\n";

  K = {{1, 1}, {2, 2}};
  b = {2, 4};
  res = lml::solveAxb(K, b);

  std::cout << "res = {" << res[0] << ", " << res[1] << "}(should be {2,0})\n";

  K = {{1, 1}, {2, 2}};
  b = {2, 5};
  res = lml::solveAxb(K, b);

  std::cout << res.empty() << "(should be 1)\n";

  K = {{0, 0}, {0, 0}};
  b = {0, 0};
  res = lml::solveAxb(K, b);

  std::cout << "res = {" << res[0] << ", " << res[1] << "}(should be {0,0})\n";

  K = {{0, 0}, {0, 0}};
  b = {1, 0};
  res = lml::solveAxb(K, b);

  std::cout << res.empty() << "(should be 1)\n";

  return 0;
}
