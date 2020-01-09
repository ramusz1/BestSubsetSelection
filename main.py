
from __future__ import print_function

import sys
sys.path.append("/opt/ibm/ILOG/CPLEX_Studio_Community129/cplex/python/3.7/x86-64_linux")

import cplex
import pandas as pd
import numpy as np

def set_mixed_integer_problem(problem, X, y, k, first_order_solution):
    problem.set_problem_name("integer lasso")
    problem.objective.set_sense(problem.objective.sense.minimize)

    n = X.shape[0]
    p = X.shape[1]

    # if k == 0 return zeros ... 

    # M = 2*np.max(np.abs(first_order_solution))
    M = 100
    z = (first_order_solution != 0) * 1
    
    #   min         x^T Q x + c^T x
    #   subject to  Ax <= b
    #               l <= x <= u
    #               some x_i's binary or integral

    Q = X.T @ X
    c = - 2 * X.T @ y

    I = np.eye(p)
    rvec = np.concatenate( (np.zeros(p), np.ones(p)) )
    # A adds the SOS-1 constraint, cplex has built in sos TODO
    # TODO: add beta norm < M_l to A 
    A = np.vstack((np.hstack((I, -M * I)), 
                   np.hstack((-I, -M * I)),
                   rvec))
    senses = "L" * (2 * p + 1)
    rhs = np.concatenate((np.zeros(2*p), [k]))
    ub = np.concatenate((np.full(p, M), np.ones(p)))
    lb = np.concatenate((np.full(p, -M), np.zeros(p)))
    obj = np.concatenate((c, np.zeros(p)))
    var_types = "C" * p + "B" * p
    # TODO start = ...

    # https://stackoverflow.com/questions/38949055/linear-and-quadratic-terms-in-cplex-objective-function
    qmat = [[[x for x in range(p)], row.tolist()] for row in Q]
    problem.objective.set_quadratic(qmat)

    var_names = [f'b{i}' for i in range(p)] + [f'z{i}' for i in range(p)]
    
    # lb = [ -cplex.infinity] * var_num
    problem.variables.add(obj=obj.tolist(), lb=lb, ub=ub, types=var_types, names=var_names)

    indices = np.arange(2*p).tolist()
    A_as_sparse_pairs = [cplex.SparsePair(ind=indices, val=row.tolist()) for row in A]
    problem.linear_constraints.add(lin_expr=A_as_sparse_pairs, rhs=rhs.tolist(), senses=senses)


def keep_top_k(x, k ):
    ind = np.argsort(np.abs(x))[-k:]
    x[ind] = 0
    return x


# projected gradients
def get_first_order_solution(problem, X, y, k):
    n = X.shape[0]
    p = X.shape[1]

    # init output (beta)
    if ( p < n ){
        # beta0 = lsfit .. #TODO
    } else {
        beta0 = X.t @ y / np.sum(X * X, axis=1)
    }
    beta0 = keep_top_k(beta0, k) 
    L = # power method TODO

    # projected gradient runs 
    beta = beta0
    beta_crit = np.inf
    for _ in range(n_runs):
        for _ in range(max_iter):
            beta_old = beta
            
            # grad descent step:
            grad = -X.T @ (y - X@beta)
            beta = beta - grad/L 

            beta = keep_top_k(beta, k)

            if numpy.linalg.norm(beta - beta_old) / np.max((numpy.linalg.norm(), 1)) < tol:
                break

        curr_crit = np.sum( (y - X @ beta)**2 )
        if curr_crit < best_crit:
            best_crit = curr_crit
            best_beta = beta
        
        beta = beta0 + 2 * np.random.rand(p) * np.max((np.abs(beta0), np.ones()))


def best_subset(X, y, k):

    p = cplex.Cplex()

    first_order_solution = get_first_order_solution(p, X, y, k)

    set_mixed_integer_problem(p, X, y, k, first_order_solution)    

    sol = p.solution

    # solution.get_status() returns an integer code
    print("Solution status = ", sol.get_status(), ":", end=' ')
    # the following line prints the corresponding string
    print(sol.status[sol.get_status()])
    print("Solution value  = ", sol.get_objective_value())

    numrows = p.linear_constraints.get_num()

    for i in range(numrows):
        print("Row %d:  Slack = %10f" % (i, sol.get_linear_slacks(i)))

    numcols = p.variables.get_num()

    for j in range(numcols):
        print("Column %d:  Value = %10f" % (j, sol.get_values(j)))

    p.write("miqpex1.lp")


if __name__ == "__main__":
    n = 100
    p = 5

    X = np.random.rand(n, p)
    # beta = np.random.randint(low=0, high=2,size=p)
    beta = np.array([1, 0, 1, 0, 1])
    print(beta)
    # intercept = np.random.rand()
    # noise = np.random.random(n) * 0.01
    y = X @ beta # + intercept + noise
    print(y)
    # check different values of k
    best_subset(X, y, k = 3)
    print(f'ans{beta}')