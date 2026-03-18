#pragma once

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace simplex {

constexpr double kEpsilon = 1e-9;

enum class ObjectiveSense {
    Maximize,
    Minimize
};

enum class ConstraintSense {
    LessEqual,
    GreaterEqual,
    Equal
};

enum class SolveStatus {
    Optimal,
    Unbounded,
    Infeasible,
    IterationLimit
};

struct AffineExpression {
    std::vector<double> coefficients;
    double constant = 0.0;
};

struct Constraint {
    std::vector<double> coefficients;
    ConstraintSense sense = ConstraintSense::LessEqual;
    double rhs = 0.0;
};

struct LinearProgram {
    ObjectiveSense objectiveSense = ObjectiveSense::Maximize;
    std::vector<double> objective;
    double objectiveConstant = 0.0;
    std::vector<Constraint> constraints;
    std::size_t variableCount = 0;
};

struct IterationSnapshot {
    int phase = 0;
    int iteration = 0;
    std::vector<std::vector<double>> tableau;
    std::vector<std::string> variableNames;
    std::vector<std::string> basicVariableNames;
    std::string enteringVariable;
    std::string leavingVariable;
};

struct SolveResult {
    SolveStatus status = SolveStatus::Optimal;
    double objectiveValue = 0.0;
    std::vector<double> variableValues;
    int iterations = 0;
    std::string message;
    std::vector<IterationSnapshot> snapshots;
};

AffineExpression parseAffineExpression(const std::string& text);
Constraint parseConstraint(const std::string& text);
LinearProgram readLinearProgram(std::istream& input, std::ostream& output);
SolveResult solve(const LinearProgram& problem, bool collectSnapshots = true);
std::string formatSnapshot(const IterationSnapshot& snapshot);
std::string formatSolution(const LinearProgram& problem, const SolveResult& result);
const char* toString(SolveStatus status);

}  // namespace simplex
