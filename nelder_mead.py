import numpy as np

# Python implementation of the Nelder-Mead algorithm

def nelder_mead(f, initial_point, initial_step, tolerance, maxiter, alpha, gamma, rho, sigma):
    
    # Dimension of vector space
    dim = len(initial_point)

    # Create the initial simplex

    # Add initial point and its loss value
    current_simplex = [[initial_point, f(initial_point)]]
    # Add initial step to each dimension of initial point
    for i in range(dim):
        initial_point_copy = np.copy(initial_point)
        initial_point_copy[i] = initial_point_copy[i] + initial_step
        current_simplex.append([initial_point_copy, f(initial_point_copy)])
    
    #print(current_simplex)
    
    iter = 0

    while True:
        
        iter += 1
        print(iter)

        # Order the points by loss value
        ordered_simplex = sorted(current_simplex, key=lambda x: float(x[1]))

        #print(ordered_simplex)

        # Check for termination
        # For now terminate when reached max # of iterations
        if iter == maxiter:
            # Return the simplex
            return ordered_simplex
        
        else:

            # Calculate centroid
            centroid_point = centroid(ordered_simplex[:-1])

            # Calculate reflected point
            worst_point = ordered_simplex[-1][0]
            reflected_point = reflected(centroid_point, alpha, worst_point)
            reflected_point_value = f(reflected_point)

            # Get value of best and 2nd worst points
            best_point_value = ordered_simplex[0][1]
            second_worst_point_value = ordered_simplex[-2][1]
            
            if best_point_value <= reflected_point_value < second_worst_point_value:

                # Replace worst point with reflected
                ordered_simplex.pop()
                ordered_simplex.append([reflected_point, reflected_point_value])
                current_simplex = ordered_simplex
                print('Swapped worst with reflected')
                continue

            elif reflected_point_value < best_point_value:

                # Calculate expanded point
                expanded_point = expanded(centroid_point, gamma, reflected_point)
                expanded_point_value = f(expanded_point)

                if expanded_point_value < reflected_point_value:

                    # Replace worst point with expanded
                    ordered_simplex.pop()
                    ordered_simplex.append([expanded_point, expanded_point_value])
                    current_simplex = ordered_simplex
                    print('Swapped worst with expanded')
                    continue

                else:

                    # Replace worst point with reflected
                    ordered_simplex.pop()
                    ordered_simplex.append([reflected_point, reflected_point_value])
                    current_simplex = ordered_simplex
                    print('Swapped worst with reflected')
                    continue

            else:

                # Get worst point value
                worst_point_value = ordered_simplex[-1][1]

                if reflected_point_value < worst_point_value:

                    # Calculate outside contracted point
                    contracted_point = contraction_outside(centroid_point, rho, reflected_point)
                    contracted_point_value = f(contracted_point)

                    if contracted_point_value < reflected_point_value:

                        # Replace worst point with outside contracted
                        ordered_simplex.pop()
                        ordered_simplex.append([contracted_point, contracted_point_value])
                        current_simplex = ordered_simplex
                        print('Swapped worst with outside contracted')
                        continue

                    else:

                        # Shrink simplex, replace all but the best point with shrunk points
                        current_simplex = shrink(ordered_simplex, sigma, f)
                        print("Shrunk points")
                        continue
                        
                elif reflected_point_value >= worst_point_value:

                    # Calculate inside contracted point
                    contracted_point = contraction_inside(centroid_point, rho, worst_point)
                    contracted_point_value = f(contracted_point)

                    if contracted_point_value < worst_point_value:

                        # Replace worst point with inside contracted
                        ordered_simplex.pop()
                        ordered_simplex.append([contracted_point, contracted_point_value])
                        current_simplex = ordered_simplex
                        print('Swapped worst with inside contracted')
                        continue

                    else:

                        # Shrink simplex, replace all but the best point with shrunk points
                        current_simplex = shrink(ordered_simplex, sigma, f)
                        print('Shrunk points')
                        continue

# Calculate the centroid of the simplex
def centroid(simplex):

    # # of points in simplex
    num_points = len(simplex)

    # Dimension of points
    dim_points = len(simplex[0][0])

    # Sum points in simplex
    point_sum = np.zeros(dim_points)
    for p in simplex:
        point_sum = point_sum + p[0]

    # Return centroid point
    return point_sum / num_points

# Calculate reflected point
def reflected(centroid, alpha, worst_point):
    return centroid + alpha * (centroid - worst_point)

# Calculate expanded point
def expanded(centroid, gamma, reflected):
    return centroid + gamma * (reflected - centroid)

# Calculate outside contracted point
def contraction_outside(centroid, rho, reflected):
    return centroid + rho * (reflected - centroid)

# Calculate inside contracted point
def contraction_inside(centroid, rho, worst):
    return centroid + rho * (worst - centroid)

# Shrink simplex
def shrink(simplex, sigma, f):

    # Get best point
    best_point = simplex[0][0]

    # Replace other points with shrunk points and their loss value
    for p in simplex[1:]:
        p[0] = best_point + sigma * (p[0] - best_point)
        p[1] = f(p[0])

    return simplex

# Termination
def termination(simplex):
    return 1

# Test loss functions
def f(p):
    dim = len(p)
    sum = 0
    for i in range(dim):
        sum += p[i] ** 2
    return sum

def g(p):
    x = p[0]
    y= p[1]
    return (x ** 2 + y - 11) ** 2 + (x + y ** 2 -7) ** 2

def h(p):
    x = p[0]
    y = p[1]
    return 2 * x ** 2 - 1.05 * x ** 4 + x ** 6 / 6 + x * y  + y ** 2

# Test the search
sol = nelder_mead(h, np.array([5. ,5.]), 1., 0.0001, 50, 1, 2, 0.5, 0.5)
print(sol)