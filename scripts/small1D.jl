using DrWatson

@quickactivate("CellPlacement")

#include("../src/Analyze.jl")
using Statistics
using DataFrames
using Plots
using Images
using FFTW
using Query
gr(c=:Spectral, size=[1024,768])

function maxFreq(df)
   results = zeros(2)
   idealPos = df[!, "ip"][1]
   oldPos = df[!, "lp"][1]
   newPos = df[!, "sp"][1]
   vertNum = df[!, "vertNum"][1]
   v = Array(df[!, "v"][1])
   numRevPairs = df[!, "numRevPairs"][1]

   oldErr = idealPos .- oldPos
   newErr = idealPos .- newPos

   twoRev = Int32(round(vertNum * 2))
   spErr = idealPos[end-twoRev:end] .- newPos[end-twoRev:end]
   lpErr = idealPos[end-twoRev:end] .- oldPos[end-twoRev:end]
   oneRev = Int32((twoRev - 1) / 2)
   relLpErr = lpErr[1:oneRev] .- lpErr[oneRev+1:end-2]
   relSpErr = spErr[1:oneRev] .- spErr[oneRev+1:end-2]

   res1 = abs.(FFTW.rfft(relLpErr[1:1024]))
   res2 = abs.(FFTW.rfft(relSpErr[1:1024]))

   res1, res2
end

function maxFreqs(df)
   circs = sort(unique(df[!,"circ"]))
   freqs = sort(unique(df[!,"engFreq"]))
   results = zeros(Float32, (size(circs)[1], size(freqs)[1], 2))
   for ind1 = 1:size(circs)[1]
      for ind2 = 1:size(freqs)[1]
         df2 = df |> @filter(_.circ == circs[ind1] && _.engFreq == freqs[ind2]) |> DataFrame
         res1, res2 = maxFreq(df2)

         results[ind1, ind2, 1] = maximum(res1)
         results[ind1, ind2, 2] = maximum(res2)
      end
   end
   circs, freqs, results
end

df = collect_results(datadir("simulations")) |> @filter(_.horz == 2020 && _.numRevPairs == 100) |> DataFrame

circs, freqs, results = maxFreqs(df)

plot(freqs, circs, results[:,:,1], st = :surface)
plot(freqs, circs, results[:,:,2], st = :surface)
