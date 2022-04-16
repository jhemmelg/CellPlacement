using DrWatson

@quickactivate("CellPlacement")

include("../src/Placement.jl")

function genPaths(d::Dict)
   @unpack circ, horz, vert, numRevPairs = d

   v, ip = idealPlacement(circ, horz, vert, numRevPairs)
   _, lp = largeStepPlacement(circ, horz, vert, numRevPairs)
   _, sp = smallStepPlacement(circ, horz, vert, numRevPairs)
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
# 700mm circumference cylinder with 70/45 spacing

allParams = Dict{String,Any}(
   "circ" => Vector(1200:25:1500),
   "horz" => 202.0,
   "vert" => 202.0,
   "numRevPairs" => 1000
)
dicts = dict_list(allParams)

for d in dicts
   fulld = genPaths(d)
   @tagsave(datadir("simulations", savename(d, "jld2")), fulld)
end
