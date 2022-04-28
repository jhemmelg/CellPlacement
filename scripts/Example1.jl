using DrWatson

using DataFrames
using Query
using Plots
gr(size=(1024, 768))

df = collect_results(datadir("Analysis1")) |> DataFrame

function plotVar(img; fname="george.png", figTitle="Example")
   lowest = minimum(img)
   highest = maximum(img)
   span = highest - lowest
   normImg = (img .- lowest) ./ span
   size = 100
   plot(Gray.(normImg[1:size,1:size]), title=figTitle)
   png(fname)
end

img1 = df[!,"img1"][1]
plotVar(img1, fname=plotsdir("Example1.png"), figTitle="Example 1")
println("img1 max: ", maximum(img1), ", min: ", minimum(img1))

img2 = df[!, "img2"][1]
plotVar(img2, fname=plotsdir("Example2.png"), figTitle="Example 2")
println("img2 max: ", maximum(img2), ", min: ", minimum(img2))
