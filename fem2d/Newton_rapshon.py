def newton_raphson(f, df, x0, tolerance=1e-7, max_iterations=1000):
    """
    Apply Newton-Raphson method to find the root of a nonlinear equation.

    Parameters:
    f (function): The function for which the root is being found.
    df (function): The derivative of the function f.
    x0 (float): Initial guess for the root.
    tolerance (float): The tolerance for convergence.
    max_iterations (int): The maximum number of iterations to perform.

    Returns:
    float: The root of the function.
    int: The number of iterations performed.
    """
    x_n = x0
    for iteration in range(max_iterations):
        f_x_n = f(x_n)
        df_x_n = df(x_n)
        if df_x_n == 0:
            print("Derivative is zero. No solution found.")
            return None, iteration

        x_next = x_n - f_x_n / df_x_n

        # Check for convergence
        if abs(x_next - x_n) < tolerance:
            return x_next, iteration + 1

        x_n = x_next

    print("Exceeded maximum iterations. No solution found.")
    return x_n, max_iterations

# Example usage:
# Define the function and its derivative


def f(x):
    return x ** 3 - x - 2


def df(x):
    return 3 * x ** 2 - 1


# Initial guess
x0 = 1.5

# Solve using Newton-Raphson
root, iterations = newton_raphson(f, df, x0)
print(f"Root: {root}, found in {iterations} iterations.")
