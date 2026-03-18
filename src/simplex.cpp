#include "simplex_solver.h"

#include <exception>
#include <iostream>

int main() {
    try {
        const simplex::LinearProgram problem = simplex::readLinearProgram(std::cin, std::cout);
        const simplex::SolveResult result = simplex::solve(problem, true);

        for (const simplex::IterationSnapshot& snapshot : result.snapshots) {
            std::cout << simplex::formatSnapshot(snapshot);
        }

        std::cout << simplex::formatSolution(problem, result);
        return result.status == simplex::SolveStatus::Optimal ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
