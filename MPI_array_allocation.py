# Script for a function that takes an array and computes its vectors and displacements
import numpy as np
def fallocate(array,size):
    """
    Function that computes the size of each chunk and displacements of the array that is being split across all the cores
    :param array: array to be split, the split occurs on axis=0
    :param size: Number of cores
    :return: split: list of sub arrays split along axis=0
             split_sizes_input: Sizes of the sub arrays
             displacements_input: Displacement vectors
    """
    array = np.ascontiguousarray(array)  # Making sure the array is C-Contiguous
    ndim = len(array.shape) # Number of dimensions of the array
    array_shape = array.shape # Array Shape

    # Splitting array into chunks along axis 0  for memory alocation
    split = np.array_split(array, size, axis=0)
    split_sizes = []

    for i in range(0, len(split), 1):
        split_sizes = np.append(split_sizes, len(split[i]))

    split_sizes_input = split_sizes * np.prod(array.shape[1:])
    displacements_input = np.insert(np.cumsum(split_sizes_input), 0, 0)[0:-1]

    return split, split_sizes_input, displacements_input
