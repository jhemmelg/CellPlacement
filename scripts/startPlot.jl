using DrWatson

@quickactivate("CellPlacement")

using DataFrames
using Plots
using Query
gr(c=:Spectral, size=[1024,768])

df = collect_results(datadir("simulations")) |> DataFrame

length = 300

plot(df[!,"lp"][1][1:length], lc=:blue, legend=:topleft, label="Old Step Distance", title="Startup Comparison")
plot!(df[!,"sp"][1][1:length], lc=:green, label="New Step Distance")
png(plotsdir("CarrStartup.png"))
