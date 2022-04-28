using DrWatson

@quickactivate("CellPlacement")

include("../src/Analyze.jl")
using Statistics
using DataFrames
using Plots
using Images
using FFTW
gr(c=:Spectral, size=[1024,768])

df = collect_results(datadir("simulations"))

index = 1
idealPos = df[!, "ip"][index]
oldPos = df[!, "lp"][index]
newPos = df[!, "sp"][index]
vertNum = df[!, "vertNum"][index]
v = Array(df[!, "v"][index])
numRevPairs = df[!, "numRevPairs"][index]
width = Int32(floor((numRevPairs - 1) * 2))

plot(oldPos[1:width], legend=:topleft)
plot!(idealPos[1:width])


function simFollow(command, tc, start)
   curr = start
   corrFactor = exp(-1 * tc)
   retVal = copy(command)
   for i in 1:size(command)[1]
      err = command[i] - curr
      curr = curr + err * corrFactor
      retVal[i] = curr
   end
   return retVal
end

fol = simFollow(oldPos, 2, 0.0)

plot(oldPos[1:200], legend=:topleft)
plot!(fol[1:200])
plot!(idealPos[1:200] .- fol[1:200])

fol2 = simFollow(newPos, 0.5, 0.0)

plot(newPos[1:200], legend=:topleft)
plot!(fol2[1:200])
plot!(idealPos[1:200] .- fol2[1:200])
