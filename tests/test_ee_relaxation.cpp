/** \file
 *
 *  A demonstration of class ElectronElectronOperator of the LoKI-B project.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2022 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Use class ElectronElectronOperator to solve the time-dependent EEDF for
 *  a system of electrons that only undergo collisions among themselves.
 *  The program shows the resulting Maxwellization of the the initial
 *  (Druyvesteyn) distribution.
 *
 *  This program is run without arguments. See the top of function main() for
 *  the precise simulation settings.
 *  The program produces a file 'eedf.dat' in the current working directory
 *  that consists of blocks that contain f(u). A block is written after every
 *  100 time steps, its header reveals some more information about the state
 *  of affairs.
 *
 *  Note that the EEDF is *not* normalized during the simulation. Since a
 *  conservative formulation is used, the EEDF will stay normalized as time
 *  progresses. This also means that the Boltzmann equation is solved for *all*
 *  grid points (in an iterative method, normalisation can be enforced  by
 *  replacing the equation for one grid point by the normalization condition
 *  before solving the system). The time-discretization uses the backward Euler
 *  approximation.
 *
 *  \author Jan van Dijk
 *  \date   July 2023
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/Operators.h"
#include "LoKI-B/EedfUtilities.h"
#include <iostream>
#include <fstream>

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

void write_eedf(std::ostream& os, const loki::Grid& grid, const loki::Vector& eedf, unsigned iter, double t, double TeV, double norm, double relDiff, double next_dt, int info)
{
    if (iter != 0)
    {
        os << std::endl << std::endl;
    }
    os << "# " << "iteration = " << iter << std::endl;
    os << "# " << "t         = " << t << std::endl;
    os << "# " << "TeV       = " << TeV << std::endl;
    os << "# " << "norm      = " << norm << std::endl;
    os << "# " << "relDiff   = " << relDiff << std::endl;
    os << "# " << "next_dt   = " << next_dt << std::endl;
    os << "# " << "info      = " << info << std::endl;
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        os << grid.getCells()(k) << '\t' << eedf(k) << std::endl;
    }
}

int main(int argc, const char* argv[])
{
    using namespace loki;

    const unsigned nCells = 64;
    const double uMax = 4; // eV
    const double ne = 1e16; // m^-3
    const double n0 = 1e22; // m^-3

    const double shape = 2.0; // 1: Maxwell, 2: Druyvesteyn

    const Grid grid(nCells,uMax);
    ElectronElectronOperator eeOperator(grid);

    /* Use a dense matrix for the system, just like the boltzmannMatrix in LoKI-B.
     * For the present problem, we know it is going to be tridiagonal and we could
     * have used something more optimal, but we want to stay close to the main code.
     * Note that a conversion to a sparse matrix is done before solving the system.
     */
    Matrix M(grid.nCells(),grid.nCells());
    Vector rhs = Vector::Zero(grid.nCells());

    /* Maximum allowed pointwise relative difference in a time step.
     */
    const double tolRelDiff=1e-3;
    const double t_max = 1.0; // s
    const unsigned iter_write = 100;

    double t=0.0; // s (initial value)
    // the initial eedf.
    double TeV = 1.0; // requested value, in eV
    Vector eedf = makePrescribedEDF(grid,shape,TeV);
    /* We know that makePrescribedEDF normalizes the matrix, do not recalculate 'norm'.
     * Note that at this stage 2*<u>/3 will *not* be equal to TeV exactly, because of
     * discretisation effects. The EEDF is normalized to make the cumulative chance 1,
     * the mean energy will not be reproduced exactly. The mismatch is smaller for
     * finer meshes. So recalculate TeV, ensure the actual value is written in the
     * upcoming call to write_eedf.
     */
    double norm=1.0;
    TeV = 2. / 3. * getMeanEnergy(eedf,grid);

    double dt=1e-10; // s (initial value)
    unsigned iter=0;

    std::cout << "Results will be written to file eedf.dat." << std::endl;
    std::ofstream ofs("eedf.dat");
    write_eedf(ofs,grid,eedf,iter,t,TeV,norm, 0, dt, 0);
    while (t<t_max)
    {
        // part 1. Update the target time and the ee-operator data.
        t += dt;
        eeOperator.update_g_ee_AB(grid,eedf,ne,n0);
        TeV = 2. / 3. * getMeanEnergy(eedf,grid);

        /* NOTE: the equation is discretized as in LoKI-B, as:
         *
         *    -(sqrt(u)/(N*gamma)*df/dt - 1/(N*gamma)*dG/du = 0
         *
         * That explains all the minus signs.
         *
         * Since we use the ElectronElectronOperator::discretizeTerm, which
	 * *adds* to the matrix, we clear the matrix first. This needs to
         * be given some thought, but let us first see if this shows up in
         * a profile before trying to optimize again.
         */
        M.setZero();
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            /* part 1: - 1/(N*gamma)*dG/du. This code is the same as that in
             * source/Operators.cpp. This *sets* the elements of matrix M.
             */
            eeOperator.discretizeTerm(M,grid);
            /* part 2: - sqrt(u)/(N*gamma)*df/dt. We use backward Euler, so
             * \f$ df/dt \approx (f-f_old)/dt) \f$. The first terms leads to an
             * additional contribution to M, the second appears as an
             * inhomogemeous term with the same sign on the RHS.
             */
            const double cf=std::sqrt(grid.getCells()(k))/(n0*SI::gamma);
            M.coeffRef(k,k) -= cf/dt;
            rhs(k) = -cf*eedf(k)/dt;
        }

        // Solve the problem. We create a sparse matrix first and solve that.
        SparseMatrix tmp(M.sparseView());
        tmp.makeCompressed();
        Eigen::SparseLU<SparseMatrix,Eigen::AMDOrdering<typename SparseMatrix::StorageIndex>> solver;
        solver.compute(tmp);
        const Vector eedf_old = eedf;
        eedf = solver.solve(rhs);

        /* Calculate the normalization status of the EEDF. Ideally this should stay 1,
         * since the algorithm is conservative and no particles should disappear, up
         * to the order of the machine accuracy.
         */
        norm = eedf.dot(grid.getCells().cwiseSqrt() * grid.du());
        /* Calculate the highest relative change of the eedf in any of the grid points.
         */
        /// \todo Use loki function maxRelDiff. Rename that to somathing less likely to clash.
        double relDiff=0.0;
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            const double localRelDiff = std::abs( (eedf(k)-eedf_old(k))/eedf_old(k) );
            relDiff=std::max(relDiff,localRelDiff);
        }
        /* Use that change to adjust the next timestep. The criterium is that the maximum
         * local relative change in one timestep should be approximately tolRelDiff.
         */
        double dtSuggested = dt*tolRelDiff/relDiff;
        double dtClamped = std::min(2*dt,dtSuggested);
        dt = dtClamped;
        ++iter;

        // write the EEDF every iter_write time steps
        if (iter%iter_write==0)
        {
                write_eedf(ofs,grid,eedf,iter,t,TeV,norm,relDiff,dt,solver.info());
        }
    }
    std::cout << "Simulation finished. Exiting." << std::endl;
    return 0;
}
