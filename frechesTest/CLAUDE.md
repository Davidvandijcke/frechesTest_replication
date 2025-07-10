# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R package called `frechesTest` that implements an ANOVA-style test for detecting discontinuities (jumps) in conditional Fréchet means of random objects in metric spaces. The package supports multiple metric spaces including probability distributions (Wasserstein distance), covariance/correlation matrices, spherical data, and network Laplacians.

## Common Development Commands

### Build and Check
```bash
# Build the package
R CMD build .

# Check the package (creates frechesTest_0.1.0.tar.gz)
R CMD check frechesTest_0.1.0.tar.gz

# Install locally
R CMD INSTALL .
```

### Using devtools (preferred for development)
```r
# Load all functions for interactive development
devtools::load_all()

# Run all tests
devtools::test()

# Run a specific test file
devtools::test(filter = "covariance")  # runs test-covariance_space.R

# Check the package
devtools::check()

# Build documentation
devtools::document()

# Install the package
devtools::install()
```

### Running Tests
```r
# Run all tests
testthat::test_local()

# Run specific test file
testthat::test_file("tests/testthat/test-general_functionality.R")
```

## Architecture and Key Components

### Core Function Structure
The package revolves around a single main function `frechesTest()` which:
1. Accepts metric-space valued data and a scalar running variable
2. Tests for jumps at a specified cutoff using local polynomial regression
3. Implements bandwidth selection via cross-validation
4. Returns test statistics and p-values

### Key Implementation Details
- **Main function**: `R/frechesTest.R` contains the primary `frechesTest()` function that orchestrates the test
- **Helper functions**: `R/frechesTest_helpers.R` contains internal functions for:
  - One-sided local polynomial weight calculation
  - Fréchet mean and variance estimation
  - Bandwidth selection via K-fold cross-validation
  - OSQP-based projection for Laplacian matrices
  - Metric-specific distance calculations
  
### Metric Space Handling
The package uses the `frechet` package for metric space computations. Different metric spaces require different preprocessing:
- **Density space**: Converts raw data to quantile functions, uses L2-Wasserstein distance
- **Covariance space**: Ensures positive definiteness, supports multiple matrix metrics
- **Sphere space**: Normalizes vectors to unit length, uses geodesic distance
- **Network space**: Handles graph Laplacian matrices with OSQP projection

### Optimization and Caching
- **Density space**: Uses OSQP for monotone quantile function projection with model caching
- **Network space**: Uses OSQP for Laplacian matrix projection with cached models
- **Pre-computation**: Quantiles are pre-calculated for all density objects before analysis

### Testing Structure
Tests are organized by metric space type:
- `test-general_functionality.R`: Core functionality and edge cases
- `test-density_space.R`: Wasserstein distance tests
- `test-covariance_space.R`: Matrix metric tests
- `test-sphere_space.R`: Spherical data tests
- `test-cv_and_undersmoothing.R`: Bandwidth selection tests

Each test file contains multiple test cases verifying both statistical properties (power, size) and implementation correctness.




# Using Gemini CLI for Large Codebase Analysis

When analyzing large codebases or multiple files that might exceed context limits, use the Gemini CLI with its massive
context window. Use `gemini -p` to leverage Google Gemini's large context capacity.

## File and Directory Inclusion Syntax

Use the `@` syntax to include files and directories in your Gemini prompts. The paths should be relative to WHERE you run the
  gemini command:

### Examples:

**Single file analysis:**
gemini -p "@src/main.py Explain this file's purpose and structure"

Multiple files:
gemini -p "@package.json @src/index.js Analyze the dependencies used in the code"

Entire directory:
gemini -p "@src/ Summarize the architecture of this codebase"

Multiple directories:
gemini -p "@src/ @tests/ Analyze test coverage for the source code"

Current directory and subdirectories:
gemini -p "@./ Give me an overview of this entire project"

# Or use --all_files flag:
gemini --all_files -p "Analyze the project structure and dependencies"

Implementation Verification Examples

Check if a feature is implemented:
gemini -p "@src/ @lib/ Has dark mode been implemented in this codebase? Show me the relevant files and functions"

Verify authentication implementation:
gemini -p "@src/ @middleware/ Is JWT authentication implemented? List all auth-related endpoints and middleware"

Check for specific patterns:
gemini -p "@src/ Are there any React hooks that handle WebSocket connections? List them with file paths"

Verify error handling:
gemini -p "@src/ @api/ Is proper error handling implemented for all API endpoints? Show examples of try-catch blocks"

Check for rate limiting:
gemini -p "@backend/ @middleware/ Is rate limiting implemented for the API? Show the implementation details"

Verify caching strategy:
gemini -p "@src/ @lib/ @services/ Is Redis caching implemented? List all cache-related functions and their usage"

Check for specific security measures:
gemini -p "@src/ @api/ Are SQL injection protections implemented? Show how user inputs are sanitized"

Verify test coverage for features:
gemini -p "@src/payment/ @tests/ Is the payment processing module fully tested? List all test cases"

When to Use Gemini CLI

Use gemini -p when:
- Analyzing entire codebases or large directories
- Comparing multiple large files
- Need to understand project-wide patterns or architecture
- Current context window is insufficient for the task
- Working with files totaling more than 100KB
- Verifying if specific features, patterns, or security measures are implemented
- Checking for the presence of certain coding patterns across the entire codebase
