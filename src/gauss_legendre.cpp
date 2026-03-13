#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <Rcpp.h>
using namespace Rcpp;

// Implementation of Gauss-Legendre quadrature
void gauss_legendre(int n, std::vector<double>& nodes, std::vector<double>& weights) {
    if (n < 2) throw std::runtime_error("n must be >= 2");

    nodes.resize(n);
    weights.resize(n);

    for (int k = 0; k < n; ++k) {
        // Initial guess using Chebyshev nodes
        double x = std::cos(M_PI * (k + 0.75) / (n + 0.5));
        double x_prev;
        int iter = 0;

        // Newton-Raphson iteration
        do {
            x_prev = x;
            double P0 = 1.0, P1 = x;
            for (int i = 2; i <= n; ++i) {
                double Pi = ((2.0*i - 1.0)*x*P1 - (i - 1.0)*P0)/i;
                P0 = P1;
                P1 = Pi;
            }
            double Pn = P1;
            double Pn_der = n * (x*Pn - P0) / (x*x - 1.0);
            x = x_prev - Pn / Pn_der;
            ++iter;
        } while (std::abs(x - x_prev) > 1e-15 && iter < 100);

        nodes[k] = x;

        // Weight
        double P0 = 1.0, P1 = x;
        for (int i = 2; i <= n; ++i) {
            double Pi = ((2.0*i - 1.0)*x*P1 - (i - 1.0)*P0)/i;
            P0 = P1;
            P1 = Pi;
        }
        double Pn_der = n * (x*P1 - P0) / (x*x - 1.0);
        weights[k] = 2.0 / ((1.0 - x*x) * Pn_der * Pn_der);
    }

    // Sort nodes and weights in ascending order
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int i, int j){ return nodes[i] < nodes[j]; });

    std::vector<double> nodes_sorted(n), weights_sorted(n);
    for (int i = 0; i < n; ++i) {
        nodes_sorted[i] = nodes[idx[i]];
        weights_sorted[i] = weights[idx[i]];
    }
    nodes = nodes_sorted;
    weights = weights_sorted;
}

// Rcpp wrapper for calling from R
// [[Rcpp::export]]
List gauss_legendre_cpp(int n, double a = -1.0, double b = 1.0) {
    std::vector<double> nodes, weights;
    gauss_legendre(n, nodes, weights); // nodes in [-1,1]
    // map to [a,b]
    std::vector<double> nodes_mapped(n), weights_mapped(n);
    double c = (a + b) / 2.0;
    double h = (b - a) / 2.0;
    for (int i = 0; i < n; ++i) {
        nodes_mapped[i] = c + h * nodes[i];
        weights_mapped[i] = h * weights[i];
    }
    return List::create(Named("nodes") = nodes_mapped,
                        Named("weights") = weights_mapped);
}