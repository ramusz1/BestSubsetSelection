
from __future__ import print_function

import sys
sys.path.append("/opt/ibm/ILOG/CPLEX_Studio_Community129/cplex/python/3.7/x86-64_linux")

import cplex
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression


def set_mixed_integer_problem(problem, X, y, k, beta_approximation):
    problem.set_problem_name("integer lasso")
    problem.objective.set_sense(problem.objective.sense.minimize)

    n = X.shape[0]
    p = X.shape[1]

    M = 2 * np.max(np.abs(beta_approximation))
    
    #   min         x^T Q x + c^T x
    #   subject to  Ax <= b
    #               l <= x <= u
    #               some x_i's binary or integral

    I = np.eye(p)
    rvec = np.concatenate( (np.zeros(p), np.ones(p)) )
    # A adds the SOS-1 constraint (b_i, 1 - z_i)
    # cplex has built in sos TODO? (this might be hard because of 1 - z_i subtraction)
    # TODO: add beta norm < M_l to A 
    A = np.vstack((np.hstack((I, -M * I)), 
                   np.hstack((-I, -M * I)),
                   rvec))
    senses = "L" * (2 * p + 1)
    rhs = np.concatenate((np.zeros(2*p), [k]))
    ub = np.concatenate((np.full(p, M), np.ones(p)))
    lb = np.concatenate((np.full(p, -M), np.zeros(p)))
    Q = X.T @ X
    c = - X.T @ y
    # add zeros because of z vector 
    Q = np.hstack((np.vstack((Q, np.zeros((p,p)) )),
                   np.zeros((2*p, p)) ))
    obj = np.concatenate((c, np.zeros(p)))

    z_approximation = (beta_approximation != 0) * 1
    starting_point = np.concatenate((beta_approximation, z_approximation))

    var_types = "C" * p + "B" * p
    var_names = [f'b{i}' for i in range(p)] + [f'z{i}' for i in range(p)]
    problem.variables.add(obj=obj.tolist(), lb=lb, ub=ub, types=var_types, names=var_names)

    all_indices = np.arange(2*p).tolist()
    A_as_sparse_pairs = [cplex.SparsePair(ind=all_indices, val=row.tolist()) for row in A]
    problem.linear_constraints.add(lin_expr=A_as_sparse_pairs, rhs=rhs.tolist(), senses=senses)

    # initial problem solution
    problem.MIP_starts.add(
        cplex.SparsePair(ind=all_indices, val=starting_point.tolist()),
        problem.MIP_starts.effort_level.auto,
        "first_order_solution")

    qmat = [[all_indices, row.tolist()] for row in Q]
    problem.objective.set_quadratic(qmat)


def keep_top_k(x, k):
    ind = np.argsort(np.abs(x))[:-k]
    x[ind] = 0
    return x


def lm(X, y):
    return LinearRegression().fit(X, y).coef_


def largest_eigen_value(X, max_iter=100):
    b = np.random.randn(X.shape[1])
    for _ in range(max_iter):
        v = X @ b
        b = v / np.linalg.norm(v)
    
    return np.linalg.norm(v) / np.linalg.norm(b)


# projected gradients
def get_first_order_solution(problem, X, y, k, n_runs=50, max_iter=1000, tolerance=1e-4):
    n = X.shape[0]
    p = X.shape[1]

    # init output (beta)
    if  p < n :
        beta0 = lm(X,y)
    else:
        beta0 = X.t @ y / np.sum(X * X, axis=1)

    beta0 = keep_top_k(beta0, k) 
    L = largest_eigen_value(X.T @ X)

    # projected gradient runs 
    beta = beta0
    best_crit = np.inf
    for _ in range(n_runs):
        for _ in range(max_iter):
            beta_old = beta
            
            # grad descent step:
            grad = -X.T @ (y - X @ beta)
            beta = beta - grad/L 

            beta = keep_top_k(beta, k)

            if np.linalg.norm(beta - beta_old) / np.max((np.linalg.norm(beta), 1)) < tolerance:
                break

        curr_crit = np.sum( (y - X @ beta)**2 )
        if curr_crit < best_crit:
            print("better beta found in gradient descent!")
            best_crit = curr_crit
            best_beta = beta
        
        beta = beta0 + 2 * np.random.rand(p) * np.max((np.abs(beta0), np.ones(p)), axis=0)

    return best_beta


def best_subset(X, y, k):

    p = cplex.Cplex()

    first_order_solution = get_first_order_solution(p, X, y, k)
    set_mixed_integer_problem(p, X, y, k, first_order_solution)    
    p.solve()
    sol = p.solution
    # solution.get_status() returns an integer code
    print("Solution status = ", sol.get_status(), ":", end=' ')
    # the following line prints the corresponding string
    print(sol.status[sol.get_status()])
    print("Solution value  = ", sol.get_objective_value())

    numrows = p.linear_constraints.get_num()
    for i in range(numrows):
        print("Row %d:  Slack = %10f" % (i, sol.get_linear_slacks(i)))

    miqp_solution = []
    numcols = p.variables.get_num()
    for j in range(numcols):
        print("Column %d:  Value = %10f" % (j, sol.get_values(j)))
        miqp_solution.append(sol.get_values(j))

    p.write("best_subset.lp")

    miqp_solution_beta = np.array(miqp_solution)[:X.shape[1]]

    return first_order_solution, miqp_solution_beta


def standarize(X, y):
    # Paper requires some standarization
    # TODO ?
    pass


def final_error(beta_pred, beta_true):
    return np.sum( (beta_true - beta_pred)**2 )


if __name__ == "__main__":
    n = 1000
    p = 100

    X = np.random.rand(n, p)

    # beta can be whatever
    # beta = np.random.randint(low=0, high=2, size=p)
    beta = np.random.rand(p) * 3
    # TODO intercept handling
    y = X @ beta # + intercept + noise

    # first_order_solution, miqp_solution = best_subset(X, y, k = np.sum(beta) )
    first_order_solution, miqp_solution = best_subset(X, y, k = beta.shape[0] // 2 )
    print('ans: ', beta)
    print('first order solution error: ', final_error(first_order_solution, beta))
    print('miqp solution error: ', final_error(miqp_solution, beta))
