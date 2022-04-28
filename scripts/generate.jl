using DrWatson

@quickactivate("CellPlacement")

include("../src/genPlacement.jl")

function makeExample(d::Dict)
   @unpack circ, horz, vert, engFreq, numRevPairs = d
   vert *= 0.1
   horz *= 0.1

   v, ip = idealPlacement(circ, horz, vert, numRevPairs)
   _, lp = largeStepPlacement(circ, horz, vert, engFreq, numRevPairs)
   _, sp = smallStepPlacement(circ, horz, vert, engFreq, numRevPairs)
   vertNum = actualVert(vert, circ)
   fulld = Dict{String,Any}(copy(d))
   fulld["vertNum"] = vertNum
   fulld["v"] = v
   fulld["ip"] = ip
   fulld["lp"] = lp
   fulld["sp"] = sp
   return fulld
end

# This script will create one data set for the standard
# 700mm circumference cylinder with 70/45 spacing)
allParams = Dict{String,Any}(
   "circ" => [400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500],
   "horz" => Vector(2020:1:2030),
   "vert" => 2020,
   "numRevPairs" => 513,
   "engFreq" => [1700, 3200, 4500, 8100]
)
dicts = dict_list(allParams)

#Threads.@threads for d in dicts
for d in dicts
   println(savename(d, "jld2"))
   println(datadir("simulations", savename(d, "jld2")))
   fulld = makeExample(d)
   @tagsave(datadir("simulations", savename(d, "jld2")), fulld)
end
