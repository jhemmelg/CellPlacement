using DrWatson

@quickactivate("CellPlacement")

using DataFrames
using Query
using Plots
gr(size=(1024,768))

df = collect_results(datadir("Analysis1")) |> DataFrame
minAng1 = minimum(df[!,"ang1"])
maxAng1 = maximum(df[!,"ang1"])
minAng2 = minimum(df[!,"ang2"])
maxAng2 = maximum(df[!,"ang2"])

min1df = df |> @filter(_.ang1 == minAng1) |> DataFrame
max1df = df |> @filter(_.ang1 == maxAng1) |> DataFrame
min2df = df |> @filter(_.ang2 == minAng2) |> DataFrame
max2df = df |> @filter(_.ang2 == maxAng2) |> DataFrame

plot(Gray.(min1df[!,"img1"][1]))
plot(Gray.(max1df[!,"img1"][1]))

maxVar1 = 0.0
minVar1 = 1000000
tot = 0.0
num = 0
for row in eachrow(df)
   record = copy(row)
   variance = maximum(record[:img1]) - minimum(record[:img1])
   global tot += variance
   global num += 1
   if variance > maxVar1
      global maxVar1 = variance
   end
   if variance < minVar1
      global minVar1 = variance
   end
end
avgVar1 = tot / num

maxVar2 = 0.0
minVar2 = 1000000
tot = 0.0
num = 0
for row in eachrow(df)
   record = copy(row)
   variance = maximum(record[:img2]) - minimum(record[:img2])
   global tot += variance
   global num += 1
   if variance > maxVar2
      global maxVar2 = variance
   end
   if variance < minVar2
      global minVar2 = variance
   end
end
avgVar2 = tot / num
