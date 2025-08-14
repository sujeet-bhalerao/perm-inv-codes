from sage.all import Partitions, GelfandTsetlinPatterns
import numpy as np
import scipy.io as sio
import os


def generate_gt_patterns_array(n, d):
    """    Generate Gelfand-Tsetlin patterns for a given n and d.
    Args:
        n (int): The size of the top row of the Gelfand-Tsetlin pattern.
        d (int): The number of rows in the Gelfand-Tsetlin pattern.
    Returns:
        np.ndarray: An array of Gelfand-Tsetlin patterns, each padded to length d.
    """
    patterns = []
    for partition in Partitions(n, max_length=d):
        padded_partition = tuple(partition) + (0,) * (d - len(partition))  
        
        GT = GelfandTsetlinPatterns(top_row=padded_partition)
        for pattern in GT.list():
            padded_GT = []
            for row in pattern:
                padded_row = list(row) + [0] * (d - len(row))
                padded_GT.append(padded_row)
            
            patterns.append(padded_GT)
    
    patterns_array = np.array(patterns)
    print("GT Patterns", patterns_array)
    return patterns_array


def generate_and_save_patterns(n_range, d_range):
    for d in d_range:
        for n in n_range:
            print(f"Generating patterns for n={n}, d={d}")
            gt_patterns_array = generate_gt_patterns_array(n, d)
            
            os.makedirs('gt_patterns', exist_ok=True)
            
            filename = f'gt_patterns/gt_patterns_n{n}_d{d}.mat'
            print(f"Saving GT patterns to {filename}")
            sio.savemat(filename, {'gt_patterns': gt_patterns_array})
            print(f"Saved {filename}")


generate_and_save_patterns(range(35, 41), range(2,4))


