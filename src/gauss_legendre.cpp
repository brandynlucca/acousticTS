#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <Rcpp.h>
using namespace Rcpp;

// Implementation of Gauss-Legendre quadrature
void gauss_legendre(int n, std::vector<double>& nodes, std::vector<double>& weights) {
    if (n < 1) throw std::runtime_error("n must be >= 1");

    if (n == 1) {
        nodes = {0.0};
        weights = {2.0};
        return;
    }

    nodes.resize(n);
    weights.resize(n);

    // Legendre roots and their weights occur in exact symmetric pairs.
    for (int k = 0; k < (n + 1) / 2; ++k) {
        // Initial guess using Chebyshev nodes
        double x = 2 * k + 1 == n ? 0.0 :
            std::cos(M_PI * (k + 0.75) / (n + 0.5));
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

        nodes[k] = -x;
        nodes[n - 1 - k] = x;

        // Weight
        double P0 = 1.0, P1 = x;
        for (int i = 2; i <= n; ++i) {
            double Pi = ((2.0*i - 1.0)*x*P1 - (i - 1.0)*P0)/i;
            P0 = P1;
            P1 = Pi;
        }
        double Pn_der = n * (x*P1 - P0) / (x*x - 1.0);
        weights[k] = 2.0 / ((1.0 - x*x) * Pn_der * Pn_der);
        weights[n - 1 - k] = weights[k];
    }
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
