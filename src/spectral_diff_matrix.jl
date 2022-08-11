"""
    spectral_diff_matrix(n, xmin=0, xmax=2π)

Retun the Fourier differentiation matrix for a uniform grid.

It is assumed that the data are given on a uniformly spaced grid from xmin to xmax, with a grid point at xmin but no grid point at xmax.

This routine is based on the matlab code in the DMSuite package by S.C. Reddy and J.A.C. Weideman, available at
http://www.mathworks.com/matlabcentral/fileexchange/29
or here:
http://dip.sun.ac.za/~weideman/research/differ.html  
"""
function spectral_diff_matrix(n, xmin=0, xmax=2π)
    h = 2π / n
    kk = 1 : n - 1
    n1 = Int(floor((n - 1) / 2))
    n2 = Int(ceil((n - 1) / 2))
    D = zeros(n, n)
    if mod(n, 2) == 0
        topc = @. 0.5 / tan((1:n2) * h / 2)
        temp = [topc; -topc[n1 : -1 : 1]]
    else
        topc = @. 0.5 / sin((1:n2) * h / 2)
        temp = [topc; topc[n1 : -1 : 1]]
    end

    col1 = 2π / (xmax - xmin) * [0; ((-1) .^ kk) .* temp]

    # Form a Toeplitz matrix
    for j in 1:n
        D[j, j : end] = -col1[1 : n + 1 - j]
        D[j, 1 : j - 1] = col1[j : -1 : 2]
    end

    return D
end