using SparseArrays, MAT, LinearAlgebra

include("FonctionsConstraintPrecond.jl")

# Test sur Stokes-Collection
stokes_path = "stokes-collection"
grids = Dict("Q1-P0" => true,
             "Q1-Q1" => true)

dict = Dict("channel_domain" => true,
            "flow_over_a_backward_facing_step" => true,
            "lid_driven_cavity" =>true,
            "colliding_flow" => true)

for grid in keys(grids)
    if grids[grid]
        for pb in keys(dict)
            if dict[pb]
                print("Reading data... ")
                file = matopen(joinpath(stokes_path,"$pb/$grid/A.mat"))
                A = read(file, "A")
                A = A'
                close(file)
                file = matopen(joinpath(stokes_path,"$pb/$grid/b.mat"))
                b = read(file, "b")[:]
                close(file)
                file = matopen(joinpath(stokes_path,"$pb/$grid/c.mat"))
                c = read(file, "c")[:]
                close(file)
                file = matopen(joinpath(stokes_path,"$pb/$grid/M.mat"))
                M = read(file, "M")
                close(file)
                file = matopen(joinpath(stokes_path,"$pb/$grid/N.mat"))
                N = read(file, "N")
                close(file)
                println("✔")

                println(norm(M-M'))

                m,n = size(A)
                mat = [M A'; A -N]
                D = [b; c]

                t1 = @elapsed x1, stats = solvePrecond(M,A,N,D, "gmres", false, "I", "G", 1e-8, 1e-8, min(m,n));
                println(grid, " -- ", pb, " Sans préconditionneur Nombre d'itérations: ", length(stats.residuals)-1, " Temps: ", t1)

                t2 = @elapsed x2, stats = solvePrecond(M,A,N,D, "gmres", true, "I", "G", 1e-8, 1e-8, min(m,n));
                println(grid, " -- ", pb, " G = I Nombre d'itérations: ", length(stats.residuals)-1, " Temps: ", t2)

                t3 = @elapsed x3, stats = solvePrecond(M,A,N,D, "gmres", true, "Diagonal", "G", 1e-8, 1e-8, min(m,n));
                println(grid, " -- ", pb, " G = Diag(M) Nombre d'itérations: ", length(stats.residuals)-1, " Temps: ", t3)
            end
        end
    end
end