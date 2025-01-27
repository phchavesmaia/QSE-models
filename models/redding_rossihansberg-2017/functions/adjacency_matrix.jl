function adjacency_matrix(grid_size::Int, diagonal_weight::Float64, straight_weight::Float64)
    "
    Creates a 'grid_size' by 'grid_size' adjacency matrix to construct a 
    network graph where each cell is connected to its adjacent cells (
    i.e. horizontally, vertically, and diagonally)

    "
    n = grid_size^2
    adj_matrix = spzeros(Int64, n, n)  # Use a sparse matrix for efficiency
    weights = spzeros(Float64, n, n) # Use a sparse matrix for efficiency
    for row in 1:grid_size
        for col in 1:grid_size
            node = (row - 1) * grid_size + col  # Linear index of the current cell

            # Iterate over neighbors (row_offset, col_offset)
            for row_offset in -1:1
                for col_offset in -1:1
                    # Skip the cell itself
                    if row_offset == 0 && col_offset == 0
                        continue
                    end

                    neighbor_row = row + row_offset
                    neighbor_col = col + col_offset

                    # Check if the neighbor is within bounds
                    if 1 <= neighbor_row <= grid_size && 1 <= neighbor_col <= grid_size
                        neighbor = (neighbor_row - 1) * grid_size + neighbor_col
                        adj_matrix[node, neighbor] = 1 # create adjacency matrix

                        # Assign weights based on the type of adjacency
                        if abs(row_offset) == 1 && abs(col_offset) == 1  # Diagonal neighbor
                            weights[node, neighbor] = diagonal_weight
                        else  # Straight (vertical/horizontal) neighbor
                            weights[node, neighbor] = straight_weight
                        end
                    end
                end
            end
        end
    end

    return adj_matrix, weights 
end