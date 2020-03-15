#include <algorithm>
#include <vector>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

/* Wrap function taking Eigen matrices */

Eigen::MatrixXd matmul_eigen(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
    return a * b;
}

/* Wrap object-oriented code */

class CustomMatrix {
public:
    /**
     * Create CustomMatrix and fill with ones.
     */
    CustomMatrix(int m, int n)
    {
        m_ = m;
        n_ = n;

        data_.resize(m * n);
        std::fill(data_.begin(), data_.end(), 1.0);
    }

    /**
     * Perform matrix multiplication
     */
    CustomMatrix dot(const CustomMatrix& other) {
        CustomMatrix result(m_, other.n_);
        for (int i = 0; i < m_; i++) {
            for (int j = 0; j < other.n_; j++) {
                double sum = 0.0;
                for (int k = 0; k < n_; k++) {
                    sum += get(i, k) * other.get(k, j);
                }
                result.set(i, j, sum);
            }
        }
        return result;
    }

private:
    std::vector<double> data_;
    int m_;
    int n_;

    double get(int i, int j) const {
        return data_[n_ * i + j];
    }

    void set(int i, int j, double value) {
        data_[n_ * i + j] = value;
    }
};

namespace py = pybind11;

PYBIND11_MODULE(matmul_pybind, m) {
    // Bind the Eigen function
    m.def("matmul_eigen", &matmul_eigen,
          py::arg("a"),
          py::arg("b"),
          "Perform matrix multiplication between two NumPy arrays");

    // Bind the class
    py::class_<CustomMatrix>(m, "CustomMatrix", "Class representing a matrix.")
        .def(py::init<int, int>(),
             py::arg("m"),
             py::arg("n"),
             "Create a new matrix filled with ones.")
        .def("dot", &CustomMatrix::dot,
             py::arg("other"),
             "Multiply this matrix with other and return the result.");
}
