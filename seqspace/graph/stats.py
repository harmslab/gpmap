def probability_of_path(paths, transition_matrix):
    """
        Calculate the probability of a path
    """
    if type(paths) == tuple:
        paths = [paths]

    # Get the path length
    path_length = len(paths[0])

    # Build a list of probabilities
    probabilities = list()
    
    # Iterate through all paths in paths-list
    for p in paths:
        
        # Begin by giving this path a probability of 1.
        pi = 1 
        
        # Iterate through edges and multiply by the 
        # transition probability of that edge.
        for i in range(path_length-1):
            pi *= transition_matrix[p[i],p[i+1]]
        
        # Append pi to probabilities
        probabilities.append(pi)
        
    return probabilities