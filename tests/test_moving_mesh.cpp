#include "LoKI-B/Grid.h"
#include "LoKI-B/GridOps.h"
#include "LoKI-B/LinearAlgebra.h"
#include <Eigen/src/Core/ArithmeticSequence.h>

using namespace loki;

// Sigmoid test function
Vector test_function(const Vector &cells)
{
    return (1. - cells.array().exp()).inverse();
}

Vector monitor_beckett(const Grid &grid, const Vector &eedf, const double beta)
{
    const auto zeta = Grid(grid.nCells(), 1.0);
    auto df_dzeta = cellDerivative(zeta, eedf);

    // Second derivative -> generally works better.
    // df_dzeta = cellDerivative(zeta, df_dzeta);

    const double delta_xi = 1. / grid.nCells();
    const double alpha =
        (delta_xi * (eedf.tail(eedf.size() - 1) - eedf.head(eedf.size() - 1)).cwiseAbs()).cwiseSqrt().sum();

    const Vector w = (1 - beta) * alpha + beta * df_dzeta.cwiseAbs().cwiseSqrt().array();

    return w;
}

Vector moving_mesh(const Grid &grid, const Vector &eedf, const uint depth = 5)
{
    Vector w = Vector::Zero(grid.nCells() + 2);
    w.segment(1, w.size() - 2) = monitor_beckett(grid, eedf, 0.8);

    Matrix monitor_matrix = Matrix::Zero(w.size() - 1, w.size() - 1);

    monitor_matrix.diagonal(-1) = w.segment(1, w.size() - 2);
    monitor_matrix.diagonal() = -(w.tail(w.size() - 1) + w.head(w.size() - 1));
    monitor_matrix.diagonal(1) = w.segment(1, w.size() - 2);

    Vector rhs = Vector::Zero(grid.nCells() + 1);
    rhs(rhs.size() - 1) = 1.;

    monitor_matrix.row(0) = rhs.reverse();
    monitor_matrix.row(monitor_matrix.rows() - 1) = rhs;

    Vector faces = Vector::Zero(grid.nCells() + 1);

    LinAlg::solveTDMA(monitor_matrix, rhs, faces);

    if (depth == 1)
    {
        return faces;
    }

    // TODO: Use an oracle function that is interpolated onto the grid.
    const Grid new_grid(faces, grid.uMax(), false);
    return moving_mesh(new_grid, test_function(new_grid.getCells()), depth - 1);
}

// NOTE: The interpolation functions are written by microsoft copilot.
double interpolate(double x, const Vector &xp, const Vector &fp)
{
    if (x <= xp(0))
        return fp(0);
    if (x >= xp(xp.size() - 1))
        return fp(fp.size() - 1);

    auto it = std::lower_bound(xp.data(), xp.data() + xp.size(), x);
    int idx = it - xp.data() - 1;

    double x0 = xp(idx), x1 = xp(idx + 1);
    double y0 = fp(idx), y1 = fp(idx + 1);

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

Vector interpolate(const Vector &x, const Vector &xp, const Vector &fp)
{
    Vector result(x.size());
    for (int i = 0; i < x.size(); ++i)
    {
        result(i) = interpolate(x(i), xp, fp);
    }
    return result;
}

double mre(const Vector &pred, const Vector &truth)
{
    return (truth - pred).cwiseQuotient(truth).cwiseAbs().mean();
}

int main()
{
    const double u_max = 6.;

    const auto grid = Grid(200, u_max);
    const auto eedf = test_function(grid.getCells());

    const auto faces = moving_mesh(grid, eedf, 3);
    const auto grid_final = Grid(faces, u_max, false);
    const auto eedf_final = test_function(grid_final.getCells());

    const auto grid_oracle = Grid(10000, u_max);
    const auto eedf_oracle = test_function(grid_oracle.getCells());

    const auto uniform_interp = interpolate(grid_oracle.getCells(), grid.getCells(), eedf);
    const auto nonuniform_interp = interpolate(grid_oracle.getCells(), grid_final.getCells(), eedf_final);

    std::cout << mre(uniform_interp, eedf_oracle) << std::endl;
    std::cout << mre(nonuniform_interp, eedf_oracle) << std::endl;
}
