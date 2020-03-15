#=
Author:
    Damian Hofmann <damian.hofmann@mpsd.mpg.de>

Example code showing how to compute the Chern number (and, along the way,
the Berry flux) for Hamiltonian on a discrete grid in 2d k-space.

This code uses the method described by
    Fukui, Hatsugai, and Suzuki, J. Phys. Soc. Japan 74, 1674 (2005)
    https://doi.org/10.1143/JPSJ.74.1674
=#

# Compute the link variable U_μ[kx,ky]
function link_variable_U(w::Array{Complex128, 4})
    # The eigenvector $\vec w_j(k_x[l],k_y[m])$ is w[l,m,:,j]
    nx, ny, _, nbands = size(w)
    
    # Shift $\vec w_j(\vec k) \to \vec w_j(\vec k + \delta\vec k_\mu)$ for $\mu = x,y$
    w_x = circshift(w, (1, 0, 0, 0))
    w_y = circshift(w, (0, 1, 0, 0))

    U = Array{Complex128}(nx, ny, 2) # Unnormalized link variable

    for ix in 1:nx, iy in 1:ny
        # Compute the overlap matrix $\matr W^{\vec j}_\mu(\vec k)$
        W = Array{Complex128}(nbands, nbands, 2)
        for ib1 in 1:nbands, ib2 in 1:nbands
            W[ib1,ib2,1] = vecdot(w[ix,iy,:,ib1], w_x[ix,iy,:,ib2])
            W[ib1,ib2,2] = vecdot(w[ix,iy,:,ib1], w_y[ix,iy,:,ib2])
        end
        # Compute the unnormalized link variable
        for μ in 1:2
            U[ix,iy,μ] = det(W[:,:,μ])
        end
    end
    
    return U ./ abs.(U) # Return normalized link variable
end

# Compute the Chern field F[kx,ky]
function chern_field_F(U::Array{Complex128, 3})
    U_x = U[:,:,1]
    U_y = U[:,:,2]
    F = U_x .* circshift(U_y, (1, 0)) ./ circshift(U_x, (0, 1)) ./ U_y
    return imag.(log.(F))
end

# Compute the Chern number by summing the Chern field
chern_number_C(F::Array{Float64, 2}) = round(sum(F) / (2π))

# Put everything together to compute the Chern number from
# an array of eigenvectors with dimensions (kx, ky, i, band)
# where i runs over the eigenvector components.
function compute_chern_number(w::Array{Complex128, 4})
    U = link_variable_U(w)
    F = chern_field_F(U)
    return chern_number_C(F)
end

#### EXAMPLE USAGE ####

# Construct the Hamiltonian
ω1(k) = norm(k)
ω2(k) = 1.0
g(k) = exp(im * atan2(k[2], k[1]))

Hamiltonian(k::Array{Float64,1}) =
    [[ω1(k) g(k)]
     [conj(g(k)) ω2(k)]]

# Diagonalize the Hamiltonian
nk = 100
nbands = 2
kvals = linspace(-π, π, nk)
ws = Array{Complex128}(nk, nk, nbands, nbands) # store eigenvectors
for ix in 1:nk, iy in 1:nk
    h = Hamiltonian([kvals[ix], kvals[iy]])
    result = eigfact!(h)
    ws[ix, iy, :, :] = result[:vectors]
end

for band in 1:2
    C = compute_chern_number(ws[:,:,:,[band]])
    println("Chern number for band #$band: $C")
end

C = compute_chern_number(ws)
println("Total Chern number: $C")
