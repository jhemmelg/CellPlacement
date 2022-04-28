using DrWatson
using Statistics

@quickactivate("CellPlacement")

df = collect_results(datadir("Analysis1")) |> DataFrame 

ratios = zeros(Float32, size(df)[1])

i = 1
for row in eachrow(df)
   bigStepVar = maximum(row["img1"]) - minimum(row["img1"])
   smallStepVar = maximum(row["img2"]) - minimum(row["img2"])
   ratio = bigStepVar / smallStepVar
   if ratio > 2.0 && smallStepVar != 0.0
      ratios[i] = ratio
      i += 1
   end
end

ratios = ratios[1:i-1]
println("Max ratio: ", maximum(ratios))
println("Min ratio: ", minimum(ratios))
println("Average  : ", mean(ratios))
println("std      : ", std(ratios))