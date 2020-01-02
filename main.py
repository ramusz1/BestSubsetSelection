
from __future__ import print_function

import sys
sys.path.append("/opt/ibm/ILOG/CPLEX_Studio_Community129/cplex/python/3.7/x86-64_linux")

import cplex
import pandas as pd
import numpy as np

def setproblemdata(p, X, y, k):
    p.set_problem_name("integer lasso")
    p.objective.set_sense(p.objective.sense.minimize)

    Q = X.T @ X
    c = - X.T @ y

    var_num = X.shape[1]
    variables = [f'b{i}' for i in range(var_num)]
    
    types = 'I' * var_num
    obj = c.tolist()

    lb = [ -cplex.infinity] * var_num
    p.variables.add(obj=obj,types=types, names=variables, lb=lb)

    p.linear_constraints.add(names = ["c0"])
    p.linear_constraints.set_rhs("c0", k)
    ones = np.ones(var_num).tolist()
    p.linear_constraints.set_linear_components("c0", [variables, np.ones(var_num).tolist()])
    p.linear_constraints.set_senses("c0", "L")

    qmat = [[[x for x in range(var_num)], row.tolist()] for row in Q]
    p.objective.set_quadratic(qmat)


def lasso_qpex(X, y, k):

    p = cplex.Cplex()
    setproblemdata(p, X, y, k)

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

    numcols = p.variables.get_num()

    for j in range(numcols):
        print("Column %d:  Value = %10f" % (j, sol.get_values(j)))

    p.write("miqpex1.lp")


if __name__ == "__main__":
    df = pd.read_csv('fitness.txt', delim_whitespace=True)
    X = np.arange(1,100)
    print(X)
    X = np.outer(X, np.array([5,7]))
    print(X)
    y = X @ np.array([2,1])
    print(y)
    # check different values of k
    lasso_qpex(X, y, k = 4)