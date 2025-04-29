# Solving a System of Linear Equations Using Jacobi and Gauss-Seidel Methods

def print_matrix(mat):
    for row in mat:
        print('  '.join(map(str, row)))

def create_identity_matrix(size):
    return [[1 if i == j else 0 for j in range(size)] for i in range(size)]

def max_row_sum_norm(matrix):
    return max([sum([abs(cell) for cell in row]) for row in matrix])

def is_dominant_diagonal(matrix):
    for i in range(len(matrix)):
        if abs(matrix[i][i]) < sum(abs(matrix[i][j]) for j in range(len(matrix)) if j != i):
            return False
    return True

def attempt_fix_dominant_diagonal(matrix):
    size = len(matrix)
    assigned_cols = [-1]*size
    new_matrix = [None]*size

    for i in range(size):
        for j in range(size):
            if abs(matrix[i][j]) >= sum(abs(matrix[i][k]) for k in range(size) if k != j):
                if assigned_cols.count(j) == 0:
                    assigned_cols[i] = j
                    break

    if -1 in assigned_cols:
        print("Could not find a dominant diagonal.")
        return matrix

    for i in range(size):
        new_matrix[assigned_cols[i]] = matrix[i]

    return new_matrix

def init_guess(size):
    return [0]*size

def jacobi_solver(coefficients, constants, tol, previous_guess, iteration=1):
    next_guess = []
    for i in range(len(coefficients)):
        sum_other = sum(coefficients[i][j]*previous_guess[j] for j in range(len(coefficients)) if i != j)
        next_val = (constants[i] - sum_other)/coefficients[i][i]
        next_guess.append(next_val)

    print(f"Iteration {iteration}: {next_guess}")

    if all(abs(next_guess[i] - previous_guess[i]) < tol for i in range(len(coefficients))):
        print(f"\nTotal Iterations: {iteration}")
        return

    jacobi_solver(coefficients, constants, tol, next_guess, iteration+1)

def gauss_seidel_solver(coefficients, constants, tol, previous_guess, iteration=1):
    current_guess = previous_guess.copy()

    for i in range(len(coefficients)):
        sum_other = sum(coefficients[i][j]*current_guess[j] for j in range(len(coefficients)) if i != j)
        current_guess[i] = (constants[i] - sum_other)/coefficients[i][i]

    print(f"Iteration {iteration}: {current_guess}")

    if all(abs(current_guess[i] - previous_guess[i]) < tol for i in range(len(coefficients))):
        print(f"\nTotal Iterations: {iteration}")
        return

    gauss_seidel_solver(coefficients, constants, tol, current_guess, iteration+1)

# --- Program Entry Point ---

print("Welcome to the Iterative Methods Solver!\n")

A = [
    [4, 2, 0],
    [2, 10, 4],
    [0, 4, 5]
]
b_vector = [2, 6, 5]
tolerance = 0.00001

method = int(input("\nWhich method would you like to use?\n1. Jacobi\n2. Gauss-Seidel\nPlease enter your choice: "))

if not is_dominant_diagonal(A):
    print("\nNo dominant diagonal detected. Attempting to rearrange...")
    A = attempt_fix_dominant_diagonal(A)
    if not is_dominant_diagonal(A):
        print("\nWarning: Still no dominant diagonal. Convergence is not guaranteed.\n")

if method == 1:
    print("\n--- Solving using Jacobi Method ---\n")
    jacobi_solver(A, b_vector, tolerance, init_guess(len(b_vector)))
elif method == 2:
    print("\n--- Solving using Gauss-Seidel Method ---\n")
    gauss_seidel_solver(A, b_vector, tolerance, init_guess(len(b_vector)))
else:
    print("Invalid choice. Program will exit.")
