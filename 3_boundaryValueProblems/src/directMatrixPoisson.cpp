#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <fstream>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void insertCoefficient(int id, int i, int j, double w, std::vector<T> &coeffs,
                       Eigen::VectorXd &b, const Eigen::VectorXd &boundary)
{
    int n = int(boundary.size());
    int id1 = i + j * n;

    if (i == -1 || i == n)
        b(id) -= w * boundary(j); // constrained coefficient
    else if (j == -1 || j == n)
        b(id) -= w * boundary(i); // constrained coefficient
    else
        coeffs.push_back(T(id, id1, w)); // unknown coefficient
}

void buildProblem(std::vector<T> &coefficients, Eigen::VectorXd &b, int n)
{
    int m = n * n;
    b.setZero();
    Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(n, 0, 1);
    Eigen::MatrixXd Y = v.rowwise().replicate(n);
    Eigen::MatrixXd X = v.transpose().colwise().replicate(n);
    Eigen::MatrixXd B = (2 + M_PI * M_PI * (1 - Y.array()) * Y.array()) * sin(M_PI * X.array()) 
                      + (2 + M_PI * M_PI * (1 - X.array()) * X.array()) * sin(M_PI * Y.array());

    Eigen::ArrayXd boundary(n);
    Eigen::ArrayXd z(n);
    boundary.setZero();
    z.setZero();
    for (int i = 0; i < n; i++)
    {
        B(i, 0) = 0.0;      // left
        B(i + m - n) = 0.0; // right
        B(0, i) = 0.0;      // top
        B(n - 1, i) = 0.0;  // bottom
    }
    B = B.reshaped(m, 1);
    b = B;

    // Eigen::ArrayXd boundary = Eigen::ArrayXd::LinSpaced(n, 0,M_PI).sin().pow(2);

    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            int id = i + j * n;
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
            {
                insertCoefficient(id, i, j, 1, coefficients, b, z);
            }
            else
            {
                insertCoefficient(id, i - 1, j, -1, coefficients, b, boundary);
                insertCoefficient(id, i + 1, j, -1, coefficients, b, boundary);
                insertCoefficient(id, i, j - 1, -1, coefficients, b, boundary);
                insertCoefficient(id, i, j + 1, -1, coefficients, b, boundary);
                insertCoefficient(id, i, j,     +4, coefficients, b, boundary);
            }
        }
    }
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Error: expected one and only one argument.\n";
        return -1;
    }

    int n = 101;   // size of the image
    int m = n * n; // number of unknowns (=number of pixels)

    // Assembly:
    std::vector<T> coefficients; // list of non-zeros coefficients
    Eigen::VectorXd b(m);        // the right hand side-vector resulting from the constraints
    buildProblem(coefficients, b, n);

    SpMat A(m, m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    if (n < 5) std::cout << "Here is A:" << std::endl << A << std::endl;

    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(A);
    // Compute the numerical factorization
    solver.factorize(A);
    // Use the factors to solve the linear system
    Eigen::VectorXd x = solver.solve(b);

    std::ofstream saveFile;
    saveFile.open("directMatrixPoisson.txt");
    saveFile.precision(16);
    saveFile << x << std::endl;
    saveFile.close();

    return 0;
}
