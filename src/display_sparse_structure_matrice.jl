begin
    using Pkg;
    Pkg.activate(mktempdir())
    Pkg.add("PlutoUI")
    Pkg.add("PyPlot")
    Pkg.add("ExtendableSparse")
    Pkg.add("BenchmarkTools")
    using PlutoUI,PyPlot,BenchmarkTools
    using SparseArrays
end;


# A function to handle sizing and return of a pyplot figure

function pyplot(f;width=3,height=3)
        clf()
        f()
        fig=gcf()
        fig.set_size_inches(width,height)
        fig
end

M=sprand(120,120,0.1)
A = sprand(90,120,0.3)
N = spzeros(90,90)
m,n = size(A)
mat = [M A'; A -N]

pyplot(width=3,height=3) do 
    spy(mat)  #,marker=".",markersize=0.5)
end