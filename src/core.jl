function get_indecies(k, n)
    a = collect(1:k-1)
    append!(a, collect(k+1:n))
    return a
end

function update_W!(W, A, B, r, m)
    for k = 1:r
        indecies = get_indecies(k, r)
        ΣWB = sum(i -> W[:, i] * B[i, k], indecies)
        W[:, k] = maximum(hcat(zeros(m, 1), (A[:, k] - ΣWB) / B[k, k]), 2)
    end
end

function update_H!(H, C, D, r, n)
    for k = 1:r
        indecies = get_indecies(k, r)
        ΣDH = sum(i -> D[k, i] * H[i, :], indecies)
        H[k, :] = maximum(hcat(zeros(n, 1), (C[k, :] - ΣDH) / D[k, k]), 2)'
    end
end

function hierarchical_als(X, W, H, r)
    (m, n) = size(X)
    while stopping_criteria
        A = X * H'
        B = H * H'
        update_W!(W, A, B, r, m)
        C = W' * X
        D = W' * W
        update_H!(H, C, D, r, n)
    end
end
