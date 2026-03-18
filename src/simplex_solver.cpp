#include "simplex_solver.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace simplex {
namespace {

enum class VariableKind {
    Original,
    Slack,
    Surplus,
    Artificial
};

struct ObjectiveSpec {
    ObjectiveSense sense = ObjectiveSense::Maximize;
    AffineExpression expression;
};

struct SimplexState {
    std::vector<std::vector<double>> tableau;
    std::vector<int> basis;
    std::vector<std::string> variableNames;
    std::vector<VariableKind> variableKinds;
};

struct PhaseOutcome {
    SolveStatus status = SolveStatus::Optimal;
    int iterations = 0;
};

void trimTrailingZeros(std::vector<double>& values);

class ExpressionParser {
public:
    explicit ExpressionParser(std::string_view text) : text_(text) {}

    AffineExpression parse() {
        AffineExpression expression;
        skipSpaces();

        if (position_ >= text_.size()) {
            throw std::invalid_argument("Expression must not be empty.");
        }

        while (position_ < text_.size()) {
            const double sign = parseSign();
            skipSpaces();

            const NumberToken numberToken = parseNumberToken();
            skipSpaces();

            if (matchVariablePrefix()) {
                const std::size_t variableIndex = parseVariableIndex();
                const double coefficient = sign * (numberToken.hasValue ? numberToken.value : 1.0);

                if (expression.coefficients.size() < variableIndex) {
                    expression.coefficients.resize(variableIndex, 0.0);
                }
                expression.coefficients[variableIndex - 1] += coefficient;
            } else if (numberToken.hasValue) {
                expression.constant += sign * numberToken.value;
            } else {
                throw error("Expected a number or variable term.");
            }

            skipSpaces();
            if (position_ >= text_.size()) {
                break;
            }

            const char current = text_[position_];
            if (current != '+' && current != '-') {
                throw error("Expected '+' or '-' between terms.");
            }
        }

        trimTrailingZeros(expression.coefficients);
        return expression;
    }

private:
    struct NumberToken {
        bool hasValue = false;
        double value = 0.0;
    };

    NumberToken parseNumberToken() {
        const std::size_t start = position_;
        bool sawDigit = false;
        bool sawDecimalPoint = false;

        if (position_ < text_.size() && text_[position_] == '.') {
            sawDecimalPoint = true;
            ++position_;
        }

        while (position_ < text_.size()) {
            const char current = text_[position_];
            if (std::isdigit(static_cast<unsigned char>(current)) != 0) {
                sawDigit = true;
                ++position_;
                continue;
            }

            if (current == '.' && !sawDecimalPoint) {
                sawDecimalPoint = true;
                ++position_;
                continue;
            }

            break;
        }

        if (!sawDigit) {
            position_ = start;
            return {};
        }

        return {true, std::stod(std::string(text_.substr(start, position_ - start)))};
    }

    double parseSign() {
        if (position_ >= text_.size()) {
            return 1.0;
        }

        if (text_[position_] == '+') {
            ++position_;
            return 1.0;
        }
        if (text_[position_] == '-') {
            ++position_;
            return -1.0;
        }
        return 1.0;
    }

    bool matchVariablePrefix() const {
        if (position_ >= text_.size()) {
            return false;
        }
        const char current = text_[position_];
        return current == 'x' || current == 'X';
    }

    std::size_t parseVariableIndex() {
        ++position_;
        const std::size_t start = position_;
        while (position_ < text_.size() &&
               std::isdigit(static_cast<unsigned char>(text_[position_])) != 0) {
            ++position_;
        }

        if (start == position_) {
            throw error("Variable names must look like x1, x2, ...");
        }

        const auto value = std::stoul(std::string(text_.substr(start, position_ - start)));
        if (value == 0) {
            throw error("Variable indices start at 1.");
        }
        return value;
    }

    void skipSpaces() {
        while (position_ < text_.size() &&
               std::isspace(static_cast<unsigned char>(text_[position_])) != 0) {
            ++position_;
        }
    }

    std::invalid_argument error(const std::string& message) const {
        std::ostringstream stream;
        stream << message << " Near position " << position_ << " in \"" << text_ << "\".";
        return std::invalid_argument(stream.str());
    }

    std::string_view text_;
    std::size_t position_ = 0;
};

std::string trim(const std::string& text) {
    const auto begin = std::find_if_not(text.begin(), text.end(), [](unsigned char c) {
        return std::isspace(c) != 0;
    });
    const auto end = std::find_if_not(text.rbegin(), text.rend(), [](unsigned char c) {
        return std::isspace(c) != 0;
    }).base();

    if (begin >= end) {
        return {};
    }
    return std::string(begin, end);
}

std::string toLower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return text;
}

void trimTrailingZeros(std::vector<double>& values) {
    while (!values.empty() && std::fabs(values.back()) <= kEpsilon) {
        values.pop_back();
    }
}

void resizeTo(std::vector<double>& values, std::size_t size) {
    if (values.size() < size) {
        values.resize(size, 0.0);
    }
}

ObjectiveSpec parseObjectiveLine(const std::string& line) {
    const std::string cleaned = trim(line);
    const std::string lowered = toLower(cleaned);

    ObjectiveSpec spec;
    std::string expressionText = cleaned;

    const auto stripPrefix = [&](std::string_view prefix) {
        std::string remainder = trim(cleaned.substr(prefix.size()));
        if (!remainder.empty() && remainder.front() == ':') {
            remainder = trim(remainder.substr(1));
        }
        return remainder;
    };

    if (lowered.rfind("maximize", 0) == 0) {
        spec.sense = ObjectiveSense::Maximize;
        expressionText = stripPrefix("maximize");
    } else if (lowered.rfind("max", 0) == 0) {
        spec.sense = ObjectiveSense::Maximize;
        expressionText = stripPrefix("max");
    } else if (lowered.rfind("minimize", 0) == 0) {
        spec.sense = ObjectiveSense::Minimize;
        expressionText = stripPrefix("minimize");
    } else if (lowered.rfind("min", 0) == 0) {
        spec.sense = ObjectiveSense::Minimize;
        expressionText = stripPrefix("min");
    }

    if (expressionText.empty()) {
        throw std::invalid_argument("Objective function must not be empty.");
    }

    spec.expression = parseAffineExpression(expressionText);
    return spec;
}

std::string readRequiredLine(std::istream& input, const std::string& errorMessage) {
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error(errorMessage);
    }
    return line;
}

int readConstraintCount(std::istream& input, std::ostream& output) {
    output << "Enter number of constraints: ";
    const std::string line = trim(readRequiredLine(input, "Could not read the number of constraints."));

    if (line.empty()) {
        throw std::invalid_argument("The number of constraints must not be empty.");
    }

    std::size_t processedCharacters = 0;
    const int value = std::stoi(line, &processedCharacters);
    if (processedCharacters != line.size()) {
        throw std::invalid_argument("The number of constraints must be an integer.");
    }
    if (value < 0) {
        throw std::invalid_argument("The number of constraints cannot be negative.");
    }
    return value;
}

std::size_t computeVariableCount(const LinearProgram& problem) {
    std::size_t variableCount = problem.variableCount;
    variableCount = std::max(variableCount, problem.objective.size());
    for (const Constraint& constraint : problem.constraints) {
        variableCount = std::max(variableCount, constraint.coefficients.size());
    }
    return variableCount;
}

ConstraintSense flipSense(ConstraintSense sense) {
    switch (sense) {
    case ConstraintSense::LessEqual:
        return ConstraintSense::GreaterEqual;
    case ConstraintSense::GreaterEqual:
        return ConstraintSense::LessEqual;
    case ConstraintSense::Equal:
        return ConstraintSense::Equal;
    }

    throw std::logic_error("Unknown constraint sense.");
}

LinearProgram normalizeProblem(const LinearProgram& inputProblem) {
    LinearProgram problem = inputProblem;
    problem.variableCount = computeVariableCount(problem);
    resizeTo(problem.objective, problem.variableCount);

    for (Constraint& constraint : problem.constraints) {
        resizeTo(constraint.coefficients, problem.variableCount);

        if (std::fabs(constraint.rhs) <= kEpsilon) {
            constraint.rhs = 0.0;
        }

        if (constraint.rhs < 0.0) {
            for (double& coefficient : constraint.coefficients) {
                coefficient = -coefficient;
            }
            constraint.rhs = -constraint.rhs;
            constraint.sense = flipSense(constraint.sense);
        }
    }

    return problem;
}

std::vector<double> subtractVectors(std::vector<double> left, const std::vector<double>& right) {
    resizeTo(left, right.size());
    for (std::size_t index = 0; index < right.size(); ++index) {
        left[index] -= right[index];
    }
    trimTrailingZeros(left);
    return left;
}

std::vector<double> computeMaximizationObjective(const LinearProgram& problem) {
    std::vector<double> objective = problem.objective;
    if (problem.objectiveSense == ObjectiveSense::Minimize) {
        for (double& coefficient : objective) {
            coefficient = -coefficient;
        }
    }
    return objective;
}

bool isConstraintSatisfied(const Constraint& constraint, const std::vector<double>& values) {
    double leftHandSide = 0.0;
    for (std::size_t index = 0; index < constraint.coefficients.size(); ++index) {
        leftHandSide += constraint.coefficients[index] * values[index];
    }

    switch (constraint.sense) {
    case ConstraintSense::LessEqual:
        return leftHandSide <= constraint.rhs + kEpsilon;
    case ConstraintSense::GreaterEqual:
        return leftHandSide + kEpsilon >= constraint.rhs;
    case ConstraintSense::Equal:
        return std::fabs(leftHandSide - constraint.rhs) <= kEpsilon;
    }

    return false;
}

double dotProduct(const std::vector<double>& left, const std::vector<double>& right) {
    const std::size_t size = std::min(left.size(), right.size());
    double result = 0.0;
    for (std::size_t index = 0; index < size; ++index) {
        result += left[index] * right[index];
    }
    return result;
}

}  // namespace

AffineExpression parseAffineExpression(const std::string& text) {
    return ExpressionParser(text).parse();
}

Constraint parseConstraint(const std::string& text) {
    const std::string cleaned = trim(text);
    const std::string operators[] = {"<=", ">=", "="};

    std::size_t operatorPosition = std::string::npos;
    std::string relation;
    for (const std::string& candidate : operators) {
        operatorPosition = cleaned.find(candidate);
        if (operatorPosition != std::string::npos) {
            relation = candidate;
            break;
        }
    }

    if (operatorPosition == std::string::npos) {
        throw std::invalid_argument("Constraints must contain one of <=, >= or =.");
    }

    if (cleaned.find(relation, operatorPosition + relation.size()) != std::string::npos) {
        throw std::invalid_argument("Constraints must contain exactly one relation operator.");
    }

    const std::string leftText = trim(cleaned.substr(0, operatorPosition));
    const std::string rightText = trim(cleaned.substr(operatorPosition + relation.size()));

    if (leftText.empty() || rightText.empty()) {
        throw std::invalid_argument("Both sides of a constraint must be provided.");
    }

    const AffineExpression left = parseAffineExpression(leftText);
    const AffineExpression right = parseAffineExpression(rightText);

    Constraint constraint;
    constraint.coefficients = subtractVectors(left.coefficients, right.coefficients);
    constraint.rhs = right.constant - left.constant;

    if (relation == "<=") {
        constraint.sense = ConstraintSense::LessEqual;
    } else if (relation == ">=") {
        constraint.sense = ConstraintSense::GreaterEqual;
    } else {
        constraint.sense = ConstraintSense::Equal;
    }

    return constraint;
}

LinearProgram readLinearProgram(std::istream& input, std::ostream& output) {
    LinearProgram problem;

    output << "Enter objective function (example: max 3x1 + 5x2): ";
    const ObjectiveSpec objective = parseObjectiveLine(
        readRequiredLine(input, "Could not read the objective function."));

    problem.objectiveSense = objective.sense;
    problem.objective = objective.expression.coefficients;
    problem.objectiveConstant = objective.expression.constant;

    const int constraintCount = readConstraintCount(input, output);
    problem.constraints.reserve(static_cast<std::size_t>(constraintCount));

    for (int index = 0; index < constraintCount; ++index) {
        output << "Constraint " << index + 1 << ": ";
        problem.constraints.push_back(
            parseConstraint(readRequiredLine(input, "Could not read a constraint line.")));
    }

    problem.variableCount = computeVariableCount(problem);
    return problem;
}

namespace {

IterationSnapshot makeSnapshot(const SimplexState& state,
                               int phase,
                               int iteration,
                               int enteringColumn,
                               int leavingRow) {
    IterationSnapshot snapshot;
    snapshot.phase = phase;
    snapshot.iteration = iteration;
    snapshot.tableau = state.tableau;
    snapshot.variableNames = state.variableNames;
    snapshot.basicVariableNames.reserve(state.basis.size());

    for (int basisIndex : state.basis) {
        snapshot.basicVariableNames.push_back(state.variableNames[static_cast<std::size_t>(basisIndex)]);
    }

    if (enteringColumn >= 0) {
        snapshot.enteringVariable = state.variableNames[static_cast<std::size_t>(enteringColumn)];
    }
    if (leavingRow >= 0) {
        snapshot.leavingVariable = state.variableNames[static_cast<std::size_t>(state.basis[static_cast<std::size_t>(leavingRow)])];
    }

    return snapshot;
}

int selectEnteringColumn(const std::vector<double>& objectiveRow) {
    const int rhsColumn = static_cast<int>(objectiveRow.size()) - 1;
    for (int column = 0; column < rhsColumn; ++column) {
        if (objectiveRow[static_cast<std::size_t>(column)] < -kEpsilon) {
            return column;
        }
    }
    return -1;
}

int selectLeavingRow(const SimplexState& state, int enteringColumn) {
    const int rowCount = static_cast<int>(state.basis.size());
    const int rhsColumn = static_cast<int>(state.tableau.front().size()) - 1;
    double bestRatio = std::numeric_limits<double>::infinity();
    int bestRow = -1;
    int bestBasisIndex = std::numeric_limits<int>::max();

    for (int row = 0; row < rowCount; ++row) {
        const double pivotCoefficient = state.tableau[static_cast<std::size_t>(row)][static_cast<std::size_t>(enteringColumn)];
        if (pivotCoefficient <= kEpsilon) {
            continue;
        }

        const double ratio = state.tableau[static_cast<std::size_t>(row)][static_cast<std::size_t>(rhsColumn)] / pivotCoefficient;
        const int basisIndex = state.basis[static_cast<std::size_t>(row)];

        if (ratio < bestRatio - kEpsilon ||
            (std::fabs(ratio - bestRatio) <= kEpsilon && basisIndex < bestBasisIndex)) {
            bestRatio = ratio;
            bestRow = row;
            bestBasisIndex = basisIndex;
        }
    }

    return bestRow;
}

void pivot(SimplexState& state, int pivotRow, int pivotColumn) {
    const std::size_t rowIndex = static_cast<std::size_t>(pivotRow);
    const std::size_t columnIndex = static_cast<std::size_t>(pivotColumn);
    const std::size_t columnCount = state.tableau.front().size();
    const double pivotValue = state.tableau[rowIndex][columnIndex];

    for (std::size_t column = 0; column < columnCount; ++column) {
        state.tableau[rowIndex][column] /= pivotValue;
    }

    for (std::size_t row = 0; row < state.tableau.size(); ++row) {
        if (row == rowIndex) {
            continue;
        }

        const double factor = state.tableau[row][columnIndex];
        if (std::fabs(factor) <= kEpsilon) {
            continue;
        }

        for (std::size_t column = 0; column < columnCount; ++column) {
            state.tableau[row][column] -= factor * state.tableau[rowIndex][column];
        }
    }

    state.basis[rowIndex] = pivotColumn;
}

void rebuildObjectiveRow(SimplexState& state, const std::vector<double>& costs) {
    const std::size_t objectiveRow = state.tableau.size() - 1;
    const std::size_t rhsColumn = state.tableau.front().size() - 1;

    std::fill(state.tableau[objectiveRow].begin(), state.tableau[objectiveRow].end(), 0.0);
    for (std::size_t column = 0; column < rhsColumn; ++column) {
        state.tableau[objectiveRow][column] = -costs[column];
    }

    for (std::size_t row = 0; row < state.basis.size(); ++row) {
        const double basisCost = costs[static_cast<std::size_t>(state.basis[row])];
        if (std::fabs(basisCost) <= kEpsilon) {
            continue;
        }

        for (std::size_t column = 0; column <= rhsColumn; ++column) {
            state.tableau[objectiveRow][column] += basisCost * state.tableau[row][column];
        }
    }

    for (double& value : state.tableau[objectiveRow]) {
        if (std::fabs(value) <= kEpsilon) {
            value = 0.0;
        }
    }
}

PhaseOutcome runPhase(SimplexState& state,
                      const std::vector<double>& costs,
                      int phase,
                      bool collectSnapshots,
                      std::vector<IterationSnapshot>& snapshots,
                      int iterationLimit) {
    rebuildObjectiveRow(state, costs);
    PhaseOutcome outcome;

    while (outcome.iterations < iterationLimit) {
        ++outcome.iterations;
        const int enteringColumn = selectEnteringColumn(state.tableau.back());
        const int leavingRow = (enteringColumn >= 0) ? selectLeavingRow(state, enteringColumn) : -1;

        if (collectSnapshots) {
            snapshots.push_back(makeSnapshot(state, phase, outcome.iterations, enteringColumn, leavingRow));
        }

        if (enteringColumn < 0) {
            outcome.status = SolveStatus::Optimal;
            return outcome;
        }
        if (leavingRow < 0) {
            outcome.status = SolveStatus::Unbounded;
            return outcome;
        }

        pivot(state, leavingRow, enteringColumn);
    }

    outcome.status = SolveStatus::IterationLimit;
    return outcome;
}

SimplexState buildInitialState(const LinearProgram& problem) {
    const std::size_t constraintCount = problem.constraints.size();
    const std::size_t originalVariableCount = problem.variableCount;

    std::size_t slackCount = 0;
    std::size_t surplusCount = 0;
    std::size_t artificialCount = 0;

    for (const Constraint& constraint : problem.constraints) {
        switch (constraint.sense) {
        case ConstraintSense::LessEqual:
            ++slackCount;
            break;
        case ConstraintSense::GreaterEqual:
            ++surplusCount;
            ++artificialCount;
            break;
        case ConstraintSense::Equal:
            ++artificialCount;
            break;
        }
    }

    const std::size_t totalVariableCount = originalVariableCount + slackCount + surplusCount + artificialCount;

    SimplexState state;
    state.tableau.assign(constraintCount + 1, std::vector<double>(totalVariableCount + 1, 0.0));
    state.basis.resize(constraintCount, -1);
    state.variableNames.reserve(totalVariableCount);
    state.variableKinds.reserve(totalVariableCount);

    for (std::size_t index = 0; index < originalVariableCount; ++index) {
        state.variableNames.push_back("x" + std::to_string(index + 1));
        state.variableKinds.push_back(VariableKind::Original);
    }

    std::size_t nextVariableIndex = originalVariableCount;
    std::size_t slackId = 1;
    std::size_t surplusId = 1;
    std::size_t artificialId = 1;

    for (std::size_t row = 0; row < constraintCount; ++row) {
        const Constraint& constraint = problem.constraints[row];
        for (std::size_t column = 0; column < originalVariableCount; ++column) {
            state.tableau[row][column] = constraint.coefficients[column];
        }

        switch (constraint.sense) {
        case ConstraintSense::LessEqual:
            state.tableau[row][nextVariableIndex] = 1.0;
            state.basis[row] = static_cast<int>(nextVariableIndex);
            state.variableNames.push_back("s" + std::to_string(slackId++));
            state.variableKinds.push_back(VariableKind::Slack);
            ++nextVariableIndex;
            break;

        case ConstraintSense::GreaterEqual:
            state.tableau[row][nextVariableIndex] = -1.0;
            state.variableNames.push_back("r" + std::to_string(surplusId++));
            state.variableKinds.push_back(VariableKind::Surplus);
            ++nextVariableIndex;

            state.tableau[row][nextVariableIndex] = 1.0;
            state.basis[row] = static_cast<int>(nextVariableIndex);
            state.variableNames.push_back("a" + std::to_string(artificialId++));
            state.variableKinds.push_back(VariableKind::Artificial);
            ++nextVariableIndex;
            break;

        case ConstraintSense::Equal:
            state.tableau[row][nextVariableIndex] = 1.0;
            state.basis[row] = static_cast<int>(nextVariableIndex);
            state.variableNames.push_back("a" + std::to_string(artificialId++));
            state.variableKinds.push_back(VariableKind::Artificial);
            ++nextVariableIndex;
            break;
        }

        state.tableau[row][totalVariableCount] = constraint.rhs;
    }

    if (state.variableNames.size() != totalVariableCount) {
        throw std::logic_error("Failed to build a consistent tableau.");
    }

    return state;
}

std::vector<double> buildPhaseOneCosts(const SimplexState& state) {
    std::vector<double> costs(state.variableNames.size(), 0.0);
    for (std::size_t column = 0; column < state.variableKinds.size(); ++column) {
        if (state.variableKinds[column] == VariableKind::Artificial) {
            costs[column] = -1.0;
        }
    }
    return costs;
}

void removeArtificialBasis(SimplexState& state) {
    std::vector<std::size_t> redundantRows;
    const std::size_t rhsColumn = state.tableau.front().size() - 1;

    for (std::size_t row = 0; row < state.basis.size(); ++row) {
        const int basisIndex = state.basis[row];
        if (state.variableKinds[static_cast<std::size_t>(basisIndex)] != VariableKind::Artificial) {
            continue;
        }

        const auto isCurrentlyBasic = [&](std::size_t column) {
            return std::find(state.basis.begin(), state.basis.end(), static_cast<int>(column)) != state.basis.end();
        };

        std::size_t pivotColumn = rhsColumn;
        for (std::size_t column = 0; column < rhsColumn; ++column) {
            if (state.variableKinds[column] == VariableKind::Artificial || isCurrentlyBasic(column)) {
                continue;
            }
            if (std::fabs(state.tableau[row][column]) > kEpsilon) {
                pivotColumn = column;
                break;
            }
        }

        if (pivotColumn != rhsColumn) {
            pivot(state, static_cast<int>(row), static_cast<int>(pivotColumn));
            continue;
        }

        if (std::fabs(state.tableau[row][rhsColumn]) > kEpsilon) {
            throw std::logic_error("Artificial variable remains basic with non-zero value after phase one.");
        }

        redundantRows.push_back(row);
    }

    for (auto iterator = redundantRows.rbegin(); iterator != redundantRows.rend(); ++iterator) {
        const std::size_t row = *iterator;
        state.tableau.erase(state.tableau.begin() + static_cast<std::ptrdiff_t>(row));
        state.basis.erase(state.basis.begin() + static_cast<std::ptrdiff_t>(row));
    }

    if (state.tableau.empty()) {
        state.tableau.emplace_back();
    }
}

void dropArtificialColumns(SimplexState& state) {
    const std::size_t oldVariableCount = state.variableNames.size();
    std::vector<int> oldToNew(oldVariableCount, -1);
    std::vector<std::string> variableNames;
    std::vector<VariableKind> variableKinds;
    variableNames.reserve(oldVariableCount);
    variableKinds.reserve(oldVariableCount);

    for (std::size_t column = 0; column < oldVariableCount; ++column) {
        if (state.variableKinds[column] == VariableKind::Artificial) {
            continue;
        }
        oldToNew[column] = static_cast<int>(variableNames.size());
        variableNames.push_back(state.variableNames[column]);
        variableKinds.push_back(state.variableKinds[column]);
    }

    std::vector<std::vector<double>> newTableau(state.tableau.size(),
                                                std::vector<double>(variableNames.size() + 1, 0.0));
    const std::size_t oldRhsColumn = state.tableau.front().size() - 1;

    for (std::size_t row = 0; row < state.tableau.size(); ++row) {
        std::size_t newColumn = 0;
        for (std::size_t oldColumn = 0; oldColumn < oldVariableCount; ++oldColumn) {
            if (oldToNew[oldColumn] < 0) {
                continue;
            }
            newTableau[row][newColumn++] = state.tableau[row][oldColumn];
        }
        newTableau[row].back() = state.tableau[row][oldRhsColumn];
    }

    for (int& basisIndex : state.basis) {
        basisIndex = oldToNew[static_cast<std::size_t>(basisIndex)];
    }

    state.tableau = std::move(newTableau);
    state.variableNames = std::move(variableNames);
    state.variableKinds = std::move(variableKinds);
}

std::vector<double> buildPhaseTwoCosts(const SimplexState& state, const std::vector<double>& maximizationObjective) {
    std::vector<double> costs(state.variableNames.size(), 0.0);
    for (std::size_t column = 0; column < maximizationObjective.size() && column < costs.size(); ++column) {
        costs[column] = maximizationObjective[column];
    }
    return costs;
}

SolveResult solveWithoutConstraints(const LinearProgram& problem) {
    SolveResult result;
    result.variableValues.assign(problem.variableCount, 0.0);

    for (const Constraint& constraint : problem.constraints) {
        if (!isConstraintSatisfied(constraint, result.variableValues)) {
            result.status = SolveStatus::Infeasible;
            result.message = "The constraints are infeasible, even without decision variables.";
            return result;
        }
    }

    const std::vector<double> maximizationObjective = computeMaximizationObjective(problem);
    const bool hasImprovingDirection = std::any_of(maximizationObjective.begin(), maximizationObjective.end(), [](double value) {
        return value > kEpsilon;
    });

    if (hasImprovingDirection) {
        result.status = SolveStatus::Unbounded;
        result.message = "The objective can increase without bound because no constraints limit the variables.";
        return result;
    }

    result.status = SolveStatus::Optimal;
    result.objectiveValue = problem.objectiveConstant;
    result.message = "Optimal solution found.";
    return result;
}

void fillVariableValues(const SimplexState& state, SolveResult& result, std::size_t originalVariableCount) {
    result.variableValues.assign(originalVariableCount, 0.0);
    const std::size_t rhsColumn = state.tableau.front().size() - 1;

    for (std::size_t row = 0; row < state.basis.size(); ++row) {
        const int basisIndex = state.basis[row];
        if (basisIndex < 0 || static_cast<std::size_t>(basisIndex) >= originalVariableCount) {
            continue;
        }
        result.variableValues[static_cast<std::size_t>(basisIndex)] = state.tableau[row][rhsColumn];
    }
}

}  // namespace

SolveResult solve(const LinearProgram& inputProblem, bool collectSnapshots) {
    const LinearProgram problem = normalizeProblem(inputProblem);
    SolveResult result;
    result.variableValues.assign(problem.variableCount, 0.0);

    if (problem.variableCount == 0 || problem.constraints.empty()) {
        result = solveWithoutConstraints(problem);
        return result;
    }

    const std::vector<double> maximizationObjective = computeMaximizationObjective(problem);
    SimplexState state = buildInitialState(problem);

    const int iterationLimit =
        std::max<int>(1000, static_cast<int>((problem.variableCount + problem.constraints.size()) * 50));
    const std::vector<double> phaseOneCosts = buildPhaseOneCosts(state);

    if (std::any_of(phaseOneCosts.begin(), phaseOneCosts.end(), [](double value) { return std::fabs(value) > kEpsilon; })) {
        const PhaseOutcome phaseOneOutcome =
            runPhase(state, phaseOneCosts, 1, collectSnapshots, result.snapshots, iterationLimit);
        result.iterations += phaseOneOutcome.iterations;

        if (phaseOneOutcome.status == SolveStatus::IterationLimit) {
            result.status = SolveStatus::IterationLimit;
            result.message = "Phase one hit the iteration limit.";
            return result;
        }
        if (phaseOneOutcome.status == SolveStatus::Unbounded) {
            result.status = SolveStatus::Unbounded;
            result.message = "Phase one became unbounded, which indicates an invalid auxiliary problem.";
            return result;
        }

        if (state.tableau.back().back() < -1e-7) {
            result.status = SolveStatus::Infeasible;
            result.message = "The linear program is infeasible.";
            return result;
        }

        removeArtificialBasis(state);
        dropArtificialColumns(state);
    }

    const std::vector<double> phaseTwoCosts = buildPhaseTwoCosts(state, maximizationObjective);
    const PhaseOutcome phaseTwoOutcome =
        runPhase(state, phaseTwoCosts, 2, collectSnapshots, result.snapshots, iterationLimit);
    result.iterations += phaseTwoOutcome.iterations;
    result.status = phaseTwoOutcome.status;

    if (result.status == SolveStatus::Optimal) {
        fillVariableValues(state, result, problem.variableCount);
        result.objectiveValue = dotProduct(problem.objective, result.variableValues) + problem.objectiveConstant;
        result.message = "Optimal solution found.";
    } else if (result.status == SolveStatus::Unbounded) {
        result.message = "The linear program is unbounded.";
    } else if (result.status == SolveStatus::IterationLimit) {
        result.message = "The solver hit the iteration limit.";
    }

    return result;
}

std::string formatSnapshot(const IterationSnapshot& snapshot) {
    std::ostringstream stream;
    const std::size_t rowCount = snapshot.basicVariableNames.size();
    const std::size_t columnCount = snapshot.variableNames.size();
    const int width = 11;

    stream << "\n===== Phase " << snapshot.phase << " - Iteration " << snapshot.iteration << " =====\n";
    stream << std::setw(width) << "BV";
    for (const std::string& variableName : snapshot.variableNames) {
        stream << std::setw(width) << variableName;
    }
    stream << std::setw(width) << "RHS" << '\n';

    stream << std::setw(width) << "Obj";
    for (std::size_t column = 0; column <= columnCount; ++column) {
        stream << std::setw(width) << std::fixed << std::setprecision(3)
               << snapshot.tableau[rowCount][column];
    }
    stream << '\n';

    for (std::size_t row = 0; row < rowCount; ++row) {
        stream << std::setw(width) << snapshot.basicVariableNames[row];
        for (std::size_t column = 0; column <= columnCount; ++column) {
            stream << std::setw(width) << std::fixed << std::setprecision(3)
                   << snapshot.tableau[row][column];
        }
        stream << '\n';
    }

    if (!snapshot.enteringVariable.empty()) {
        stream << "Entering variable: " << snapshot.enteringVariable << '\n';
    }
    if (!snapshot.leavingVariable.empty()) {
        stream << "Leaving variable: " << snapshot.leavingVariable << '\n';
    }

    return stream.str();
}

std::string formatSolution(const LinearProgram& problem, const SolveResult& result) {
    std::ostringstream stream;
    stream << "\nStatus: " << toString(result.status) << '\n';
    stream << result.message << '\n';
    stream << "Iterations: " << result.iterations << '\n';

    if (result.status == SolveStatus::Optimal) {
        stream << std::fixed << std::setprecision(6);
        stream << "Objective value: " << result.objectiveValue << '\n';
        for (std::size_t index = 0; index < result.variableValues.size(); ++index) {
            stream << "x" << index + 1 << " = " << result.variableValues[index] << '\n';
        }
    }

    stream << "Assumption: x_i >= 0 for all variables.\n";
    stream << "Problem type: "
           << (problem.objectiveSense == ObjectiveSense::Maximize ? "maximization" : "minimization")
           << '\n';
    return stream.str();
}

const char* toString(SolveStatus status) {
    switch (status) {
    case SolveStatus::Optimal:
        return "optimal";
    case SolveStatus::Unbounded:
        return "unbounded";
    case SolveStatus::Infeasible:
        return "infeasible";
    case SolveStatus::IterationLimit:
        return "iteration_limit";
    }

    return "unknown";
}

}  // namespace simplex
