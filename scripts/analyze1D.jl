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

df = collect_results(datadir("simulations")) |> @filter(_.horz == 2020) |> DataFrame

index = 4
idealPos = df[!, "ip"][index]
oldPos = df[!, "lp"][index]
newPos = df[!, "sp"][index]
vertNum = df[!, "vertNum"][index]
v = Array(df[!, "v"][index])
numRevPairs = df[!, "numRevPairs"][index]
width = Int32(floor((numRevPairs - 1) * 2))

oldErr = idealPos .- oldPos
newErr = idealPos .- newPos

twoRev = Int32(round(vertNum * 2))
spErr = idealPos[end-twoRev:end] .- newPos[end-twoRev:end]
lpErr = idealPos[end-twoRev:end] .- oldPos[end-twoRev:end]
oneRev = Int32((twoRev - 1) / 2)
relLpErr = lpErr[1:oneRev] .- lpErr[oneRev+1:end-2]
relSpErr = spErr[1:oneRev] .- spErr[oneRev+1:end-2]
plot(relLpErr)
plot(relSpErr)
rfftLp = FFTW.rfft(relLpErr[1:4096])
rfftSp = FFTW.rfft(relSpErr[1:4096])
plot(abs.(rfftLp))
plot(abs.(rfftSp))
