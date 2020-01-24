import cplex
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression



def set_mixed_integer_problem(problem, X, y, k, beta_approximation, M):
    problem.set_problem_name("integer lasso")
    problem.objective.set_sense(problem.objective.sense.minimize)

    n = X.shape[0]
    p = X.shape[1]

    A = np.concatenate( (np.zeros(p), np.ones(p)) )
    senses = "G"
    # sum_z <= k
    # - sum_z => -k
    # p - sum_z => p - k
    # sum z_inv => p - k
    rhs = [p - k]
    ub = np.concatenate((np.full(p, M), np.ones(p)))
    lb = np.concatenate((np.full(p, -M), np.zeros(p)))
    Q = X.T @ X
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


def keep_top_k(x, k):
    ind = np.argsort(np.abs(x))[:-k]
    x[ind] = 0
    return x


def lm(X, y):
    return LinearRegression().fit(X, y).coef_


def get_fast_beta_approx(X, y, k):
    beta = lm(X,y)
    return keep_top_k(beta, k) 


def largest_eigen_value(X, max_iter=100):
    b = np.random.randn(X.shape[1])
    for _ in range(max_iter):
        v = X @ b
        b = v / np.linalg.norm(v)
    
    return np.linalg.norm(v) / np.linalg.norm(b)


# projected gradients
def get_first_order_solution(X, y, k, n_runs=50, max_iter=1000, tolerance=1e-4):
    n = X.shape[0]
    p = X.shape[1]

    # init output (beta)
    beta0 = get_fast_beta_approx(X, y, k)
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


def abs_corr(x, y):
    return np.corrcoef(x, y)[0][1]


def get_miu(X):
    X_T = X.T
    max_corr = -np.inf
    for i in range(X_T.shape[0]):
        for j in range(i + 1, X_T.shape[0]):
            tmp = abs_corr(X_T[i], X_T[j])
            max_corr = max(tmp, max_corr)

    return max_corr


def theory_driven_big_M(X, y, k):
    miu = get_miu(X)
    if miu * (k - 1) >= 1:
        raise ValueError("miu[k-1] >= 1, switching to mild start")
    gamma_k = 1 - miu * (k - 1) # minimum bound on gamma_k

    ## first expression
    X_T = X.T
    correlations = np.zeros(X_T.shape[0])
    for i in range(X_T.shape[0]):
        correlations[i] = abs_corr(X_T[i], y)
    
    k_highest_correlations = - np.partition( - correlations, k)[:k]
    k_highest_correlations = k_highest_correlations ** 2
    first_expr = 1 / gamma_k * np.sqrt(np.sum(correlations))

    ## second expression
    second_expr = 1 / np.sqrt(gamma_k) * np.linalg.norm(y, ord=2)

    big_M = min(first_expr, second_expr)
    return big_M


def best_subset(X, y, k, start = 'warm'):

    if X.shape[0] < X.shape[1]:
        print("WARNING, case p > n is not properly implemented yet")

    p = cplex.Cplex()
    if 'warm' == start:
        beta_approximation = get_first_order_solution(X, y, k)
        warm_start_tau = 2
        big_M = warm_start_tau * np.max(np.abs(beta_approximation))
    elif 'cold' == start:
        big_M = np.inf
        beta_approximation = None
    elif 'mild' == start:
        beta_approximation = get_fast_beta_approx(X, y, k) 
        mild_start_tau = 5
        big_M = mild_start_tau * np.max(np.abs(beta_approximation))
    elif 'theory' == start:
        big_M = theory_driven_big_M(X, y, k)
        beta_approximation = None
    else :
        raise ValueError("start should be one of: warm, mild, cold")

    print(f"Solving MIQP with {start} start")

    print("[DEBUG] big_M", big_M)

    set_mixed_integer_problem(p, X, y, k, beta_approximation, big_M)

    p.solve()
    sol = p.solution
    # solution.get_status() returns an integer code
    print("Solution status = ", sol.get_status(), ":", end=' ')
    # the following line prints the corresponding string
    print(sol.status[sol.get_status()])
    print("Solution value  = ", sol.get_objective_value())

    miqp_solution = []
    numcols = p.variables.get_num()
    for j in range(numcols):
        miqp_solution.append(sol.get_values(j))
    
    p.write("best_subset.lp")

    miqp_solution_beta = np.array(miqp_solution)[:X.shape[1]]
    miqp_solution_z_inv = np.array(miqp_solution)[X.shape[1]:]
    miqp_solution_z = 1 - miqp_solution_z_inv

    print(f"BETA: {miqp_solution_beta}\nZ: {miqp_solution_z}")

    return beta_approximation, miqp_solution_beta


def run(X, y, k, start):
    y = y.flatten()
    k = int(k)
    beta_approximation, miqp_solution = best_subset(X, y, k, start)
    return beta_approximation, miqp_solution


def tests():
    n = 100
    p = 10
    X = np.random.rand(n, p)
    beta = np.random.rand(p) * 3
    k = 5
    z = np.concatenate((np.ones(k), np.zeros(p-k)))
    beta *= z 
    y = X @ beta 
    for start in ['warm', 'cold', 'mild', 'theory']:
        beta_approximation, miqp_solution = best_subset(X, y, k, start)
        print('ANS beta', beta)
        print('ANS z', z)