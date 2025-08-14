function partitions = generate_partitions(n, d)
    partitions = generate_partitions_recursive(n, d, []);
end


function partitions = generate_partitions_recursive(n, d, current_partition)
    partitions = [];
    
    if length(current_partition) == d
        if sum(current_partition) == n
            partitions = [partitions; current_partition];
        end
        return;
    end
    
    max_value = n - sum(current_partition);
    if ~isempty(current_partition)
        max_value = min(max_value, current_partition(end));
    end
    
    for i = max_value:-1:0
        new_partition = [current_partition, i];
        partitions = [partitions; generate_partitions_recursive(n, d, new_partition)];
    end
end
