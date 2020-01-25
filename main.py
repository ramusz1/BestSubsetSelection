import cplex
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression



def set_mixed_integer_problem(problem, X, y, k, xtx, beta_approximation, M):
    problem.set_problem_name("integer lasso")
    problem.objective.set_sense(problem.objective.sense.minimize)

    n = X.shape[0]
    p = X.shape[1]

    A = np.concatenate( (np.zeros(p), np.ones(p)) )
    rhs = [p - k]
    senses = "G"
    ub = np.concatenate((np.full(p, M), np.ones(p)))
    lb = np.concatenate((np.full(p, -M), np.zeros(p)))
    Q = xtx
    c = - X.T @ y
    # add zeros because of z vector 
    Q = np.hstack((np.vstack((Q, np.zeros((p,p)) )),
                   np.zeros((2*p, p)) ))
    obj = np.concatenate((c, np.zeros(p)))

    var_types = "C" * p + "B" * p
    var_names = [f'b{i}' for i in range(p)] + [f'z{i}' for i in range(p)]
    problem.variables.add(obj=obj.tolist(), lb=lb, ub=ub, types=var_types, names=var_names)

    all_indices = np.arange(2*p).tolist()
    problem.linear_constraints.add(lin_expr=[cplex.SparsePair( all_indices, A.tolist())], rhs=rhs, senses=senses)

    if beta_approximation is not None:
        z_inv = (beta_approximation == 0) * 1
        starting_point = np.concatenate((beta_approximation, z_inv))

        # initial problem solution
        problem.MIP_starts.add(
            cplex.SparsePair(ind=all_indices, val=starting_point.tolist()),
            problem.MIP_starts.effort_level.auto,
            "first_order_solution")
        
        sos_weights = []
        for is_zero in z_inv:
            weights = [1,2] if is_zero == 1 else [2,1]
            sos_weights.append(weights)

    else:
        sos_weights = [np.random.permutation([1, 2]).tolist() for _ in range(p)]


    # sos
    for i in range(p):
        problem.SOS.add(type = "1", name = "sos_{i}",
                        SOS = cplex.SparsePair(ind = [i, i + p],
                                               val = sos_weights[i]))


    qmat = [[all_indices, row.tolist()] for row in Q]
    problem.objective.set_quadratic(qmat)


def get_solution(solved_problem):
    sol = solved_problem.solution
    # solution.get_status() returns an integer code
    print("Solution status = ", sol.get_status(), ":", end=' ')
    # the following line prints the corresponding string
    print(sol.status[sol.get_status()])
    print("Solution value  = ", sol.get_objective_value())

    miqp_solution = []
    numcols = solved_problem.variables.get_num()
    for j in range(numcols):
        miqp_solution.append(sol.get_values(j))
    
    beta_length = len(miqp_solution) // 2
    solution_beta = np.array(miqp_solution)[:beta_length]
    solution_z_inv = np.array(miqp_solution)[beta_length:]
    solution_z = 1 - solution_z_inv
    return solution_beta, solution_z


def run_cplex_miqp(X, y, k, xtx, beta0, big_M, timelimit=100):
    y = y.flatten()
    k = int(k)
    if beta0 is not None:
        beta0 = np.array(beta0).flatten()
    p = cplex.Cplex()
    p.parameters.timelimit.set(timelimit)
    set_mixed_integer_problem(p, X, y, k, xtx, beta0, big_M)
    p.solve()
    beta, z = get_solution(p)
    return beta


def tests():
    n = 100
    p = 10
    X = np.random.rand(n, p)
    beta = np.random.rand(p) * 3
    k = 5
    z = np.concatenate((np.ones(k), np.zeros(p-k)))
    beta *= z 
    y = X @ beta 

    beta_solution = run_cplex_miqp(X, y, k, X.T @ X, None, np.inf, 1.5)
    print(beta_solution)
