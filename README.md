# Simplex Method Solver

A console-based linear programming solver written in C++17. The project now uses a small reusable solver core, a CLI frontend, CMake, and a lightweight test executable.

## What It Solves

The program assumes non-negative decision variables (`x_i >= 0`) and supports:

- Maximization and minimization objectives
- Constraints with `<=`, `>=`, and `=`
- Expressions with or without spaces, such as `3x1 + 5x2`
- Constants on either side of a constraint, such as `x1 + 2 <= x2 + 5`

Internally, the solver uses a two-phase simplex approach so it can handle cases that need artificial variables, not just the simplest slack-only form.

## Build

```powershell
cmake --preset default
cmake --build --preset default
```

The executable will be created at `build/simplex.exe`.

## Run

```powershell
.\build\simplex.exe
```

Example session:

```text
Enter objective function (example: max 3x1 + 5x2): max 3x1 + 5x2
Enter number of constraints: 2
Constraint 1: 2x1 + 3x2 <= 8
Constraint 2: 2x1 + x2 <= 4
```

## Test

```powershell
ctest --test-dir build --output-on-failure
```

## Project Layout

- `include/simplex_solver.h`: public solver API
- `src/simplex_solver.cpp`: parser and simplex implementation
- `src/simplex.cpp`: command-line entry point
- `tests/simplex_tests.cpp`: smoke and regression tests
