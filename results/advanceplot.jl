#using CairoMakie
using LaTeXStrings

print(Sys.which("lualatex"))

a = 1
print(a)
using PGFPlotsX

@pgf Axis(
    {
        xlabel = L"x",
        ylabel = L"f(x) = x^2 - x + 4"
    },
    Plot(
        Expression("x^2 - x + 4")
    )
)



