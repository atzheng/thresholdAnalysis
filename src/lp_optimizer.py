import cplex


# Given a list of thresholds and their objective functions, finds the
# optimal value. num_thresholds is the number of thresholds to return;
# defaults to infinity. M is any value larger than any f1 + f2.
def find_optimal_threshold(thresholds,
                           f1,
                           f2,
                           num_thresholds = float('Inf'),
                           M = None):
    if M is None:
        M = max(f1) + max(f2)
    # Initialize solver and populate data
    problem = populate_solver(init_solver(), f1, f2, M)
    # Get list of thresholds and their lambda ranges
    results = traverse_extreme_points(problem, thresholds, f2)
    return results


# Initializes the solver with some basic parameters
def init_solver():
    problem = cplex.Cplex()
    # Allows running one simplex iteration at a time
    problem.parameters.simplex.limits.iterations.set(1)
    # Specifies use of primal simplex to solve the problem
    problem.parameters.lpmethod.set(
        problem.parameters.lpmethod.values.primal)
    return problem


# Populates constraints and variables in the solver based on some
# input data
def populate_solver(problem, f1, f2, M):
    problem.objective.set_sense(problem.objective.sense.maximize)
    varnames = ['F', 'lambda']
    problem.variables.add(obj = [0.0, M],
                           lb = [0.0, 0.0],
                           names = varnames)

    for f1_val, f2_val in zip(f1, f2):
        problem.linear_constraints.add(
            lin_expr = [cplex.SparsePair(varnames,
                                         [1, f2_val - f1_val])],
            senses = 'G',
            rhs = [f2_val])

    problem.linear_constraints.add(
        lin_expr = [cplex.SparsePair(varnames, [0., 1.])],
        senses = 'R',
        rhs = [1.0],
        range_values = [-1.0])

    return problem


# Uses simplex to iterate through the extreme points
def traverse_extreme_points(problem, thresholds, f2):
    initial_solution = set_initial_basis(problem, f2)
    extreme_points = []
    nonbasic_slacks = [initial_solution]
    while (problem.solution.get_status()
           != problem.solution.status.optimal):
        problem.solve()

        lambda_val = problem.solution.get_values()[1]
        new_nonbasic_slacks = argzero(
            problem.solution.get_linear_slacks())
        threshold_idx = filter(lambda x: x not in nonbasic_slacks,
                               new_nonbasic_slacks)
        if len(threshold_idx) > 0:
            extreme_points.append((lambda_val, threshold_idx[0]))

        nonbasic_slacks = new_nonbasic_slacks

    extreme_points.reverse()
    extreme_points.append((0, initial_solution))
    old_lambda = 1
    results = []
    for lambda_val, threshold_idx in extreme_points[1:]:
        results.append((threshold_idx, old_lambda - lambda_val))
        old_lambda = lambda_val

    return [{'threshold': thresholds[threshold_idx],
             'lambda_range': lambda_range}
            for threshold_idx, lambda_range in results]


# Sets the initial basis to lambda = 0
def set_initial_basis(problem, f2):
    # Find the index of the constratint with the maximum f2 value
    max_f2_idx = argmax(f2)
    # Set basis
    start = problem.start.status
    basic_constraints = [start.basic] * (len(f2) + 1)
    basic_constraints[max_f2_idx] = start.at_lower_bound

    problem.start.set_basis([start.basic, start.at_lower_bound],
                            basic_constraints)
    return max_f2_idx


# returns a list of all indices at which input_list is 0
def argzero(input_list):
    zero_idx = []
    for i, val in enumerate(input_list):
        if val == 0:
            zero_idx.append(i)
    return zero_idx


# Returns the first index corresponding the max value of input_list
def argmax(input_list):
    max_idx = 0
    max_val = 0
    for i, val in enumerate(input_list):
        if val > max_val:
            max_val = val
            max_idx = i
    return max_idx


if __name__ == '__main__':
    import sys
    import csv
    fname = sys.argv[1]
    handle = open(fname)
    dr = csv.DictReader(handle)

    input_data = [row for row in dr]
    thresholds = [float(row['t']) for row in input_data]
    f1 = [float(row['f1']) for row in input_data]
    f2 = [float(row['f2']) for row in input_data]
    solution = find_optimal_threshold(thresholds, f1, f2)

    outfile = open(sys.argv[2], 'wb')
    outfile.write('threshold,lambda_range\n')
    for row in solution:
        outfile.write('%.8f,%.8f\n' %(row['threshold'],
                                      row['lambda_range']))
