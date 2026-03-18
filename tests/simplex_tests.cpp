#include "simplex_solver.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void requireNear(double actual, double expected, double tolerance, const std::string& message) {
    if (std::fabs(actual - expected) > tolerance) {
        throw std::runtime_error(message + " Expected " + std::to_string(expected) +
                                 ", got " + std::to_string(actual) + ".");
    }
}

void testParserHandlesSpacesAndConstants() {
    const simplex::Constraint constraint = simplex::parseConstraint("x1 + 2 <= x2 + 5");
    require(constraint.sense == simplex::ConstraintSense::LessEqual, "Constraint sense parsing failed.");
    require(constraint.coefficients.size() == 2, "Constraint variable count parsing failed.");
    requireNear(constraint.coefficients[0], 1.0, 1e-9, "x1 coefficient mismatch.");
    requireNear(constraint.coefficients[1], -1.0, 1e-9, "x2 coefficient mismatch.");
    requireNear(constraint.rhs, 3.0, 1e-9, "Constraint RHS mismatch.");
}

void testMaximizationProblem() {
    simplex::LinearProgram problem;
    problem.objectiveSense = simplex::ObjectiveSense::Maximize;
    problem.objective = {3.0, 5.0};
    problem.constraints = {
        simplex::parseConstraint("2x1 + 3x2 <= 8"),
        simplex::parseConstraint("2x1 + x2 <= 4")
    };
    problem.variableCount = 2;

    const simplex::SolveResult result = simplex::solve(problem, false);
    require(result.status == simplex::SolveStatus::Optimal, "Expected an optimal solution.");
    requireNear(result.objectiveValue, 40.0 / 3.0, 1e-6, "Objective value mismatch.");
    requireNear(result.variableValues[0], 0.0, 1e-6, "x1 mismatch.");
    requireNear(result.variableValues[1], 8.0 / 3.0, 1e-6, "x2 mismatch.");
}

void testMinimizationWithGreaterEqualConstraints() {
    simplex::LinearProgram problem;
    problem.objectiveSense = simplex::ObjectiveSense::Minimize;
    problem.objective = {3.0, 2.0};
    problem.constraints = {
        simplex::parseConstraint("x1 + x2 >= 4"),
        simplex::parseConstraint("x1 + 2x2 >= 6")
    };
    problem.variableCount = 2;

    const simplex::SolveResult result = simplex::solve(problem, false);
    require(result.status == simplex::SolveStatus::Optimal, "Expected an optimal solution for minimization.");
    requireNear(result.objectiveValue, 8.0, 1e-6, "Minimization objective mismatch.");
    requireNear(result.variableValues[0], 0.0, 1e-6, "Minimization x1 mismatch.");
    requireNear(result.variableValues[1], 4.0, 1e-6, "Minimization x2 mismatch.");
}

void testInfeasibleProblem() {
    simplex::LinearProgram problem;
    problem.objectiveSense = simplex::ObjectiveSense::Maximize;
    problem.objective = {1.0};
    problem.constraints = {
        simplex::parseConstraint("x1 <= 3"),
        simplex::parseConstraint("x1 >= 5")
    };
    problem.variableCount = 1;

    const simplex::SolveResult result = simplex::solve(problem, false);
    require(result.status == simplex::SolveStatus::Infeasible, "Expected infeasible status.");
}

void testUnboundedProblem() {
    simplex::LinearProgram problem;
    problem.objectiveSense = simplex::ObjectiveSense::Maximize;
    problem.objective = {1.0, 1.0};
    problem.constraints = {
        simplex::parseConstraint("x1 - x2 >= 0")
    };
    problem.variableCount = 2;

    const simplex::SolveResult result = simplex::solve(problem, false);
    require(result.status == simplex::SolveStatus::Unbounded, "Expected unbounded status.");
}

}  // namespace

int main() {
    try {
        testParserHandlesSpacesAndConstants();
        testMaximizationProblem();
        testMinimizationWithGreaterEqualConstraints();
        testInfeasibleProblem();
        testUnboundedProblem();
        std::cout << "All simplex tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
}
