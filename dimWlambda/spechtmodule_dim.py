import scipy.io as sio
from sage.all import Partitions, Partition

def generate_partitions(n, d):
    partitions = []
    for partition in Partitions(n):
        if len(partition) <= d:
            partitions.append(list(partition) + [0]*(d - len(partition)))
    return partitions

def compute_dim_W_lambda(lambda_partition):
    return Partition(lambda_partition).dimension()

n_values = range(35, 41)
d_values = [3]

for n in n_values:
    for d in d_values:
        print(f'Processing n={n}, d={d}...')
        partitions = generate_partitions(n, d)
        dim_W_lambda_list = []
        for partition in partitions:
            lambda_partition = [p for p in partition if p > 0]

            if not lambda_partition:
                dim_W_lambda_list.append(0)
            else:
                dim = compute_dim_W_lambda(lambda_partition)
                dim_W_lambda_list.append(dim)
                
        filename = f'dim_W_lambda_n{n}_d{d}.mat'
        sio.savemat(filename, {'dim_W_lambda': dim_W_lambda_list})
        print(f'Saved: {filename}')
