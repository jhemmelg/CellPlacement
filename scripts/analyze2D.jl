using DrWatson

@quickactivate("CellPlacement")

include("../src/Analyze.jl")
using Statistics
using DataFrames
using Plots
using Images
using FFTW
using Query
gr(c=:Spectral, size=[1024,768])

df = collect_results(datadir("simulations")) |> @filter(_.numRevPairs == 513) |> DataFrame

index = 1
idealPos = df[!, "ip"][index]
oldPos = df[!, "lp"][index]
newPos = df[!, "sp"][index]
vertNum = df[!, "vertNum"][index]
v = Array(df[!, "v"][index])
maxV = maximum(v)
numRevPairs = df[!, "numRevPairs"][index]
width = Int32(floor((numRevPairs - 1) * 2))

function genSpacings(arrV, arrH, vert, meth)
   revs = Int32(size(arrH)[1] / vert)
   advance = Int32(floor(vert))
   height = Int32(vert * 2)
   retVal = zeros(Float32, Int32(height), Int32(revs - 1))
   for col in 1:revs-1
      for row in 1:Int32(vert*2)-1
         i = (col - 1) * height + row
         left = Int32(floor((i + 1) / 2))
         right = Int32(left + advance + rem(i + 1, 2))
         retVal[row,col] = meth(arrH[left], arrV[left], arrH[right], arrV[right])
      end
   end
   retVal
end

function horzOnly(hLeft, vLeft, hRight, vRight)
   hRight - hLeft
end

function trueDist(hLeft, vLeft, hRight, vRight)
   ((hRight - hLeft)^2 + (vRight - vLeft)^2)^0.5
end

function genPeakFreqsImgs(df, index)
   oldPos = df[!, "lp"][index]
   newPos = df[!, "sp"][index]
   vertNum = df[!, "vertNum"][index]
   v = Array(df[!, "v"][index])
   numRevPairs = df[!, "numRevPairs"][index]
   width = Int32(floor((numRevPairs - 1) * 2))
    
   hdiff1 = genSpacings(v, oldPos, vertNum, trueDist)[1:width,1:width]
   hdiff2 = genSpacings(v, newPos, vertNum, trueDist)[1:width,1:width]

   fft1 = FFTW.fft(hdiff1)
   nfft1 = norm.(fft1)
   fft2 = FFTW.fft(hdiff2)
   nfft2 = norm.(fft2)

   nfft1[1,1] = 0
   nfft2[1,1] = 0

   max1, pos1 = findmax(nfft1)
   max2, pos2 = findmax(nfft2)

   #peaks1 = findall(x ->x > max1 ^ 0.5, nfft1)
   #peaks2 = findall(x ->x > max2 ^ 0.5, nfft2)
   peaks1 = findall(x ->x > max1 * 0.1, nfft1)
   peaks2 = findall(x ->x > max2 * 0.1, nfft2)
   println(peaks1)
   println(peaks2)

   new1 = zeros(ComplexF32, size(fft1))
   new2 = zeros(ComplexF32, size(fft2))

   for ind in eachindex(peaks1)
      if (peaks1[ind] != CartesianIndex(513,513))
         new1[peaks1[ind]] = fft1[peaks1[ind]]
      end
   end

   for ind in eachindex(peaks2)
      if (peaks2[ind] != CartesianIndex(513,513))
         new2[peaks2[ind]] = fft2[peaks2[ind]]
      end
   end

   #new1[peaks1] = fft1[peaks1]
   #new2[peaks2] = fft2[peaks2]

   nimg1 = FFTW.ifft(new1)
   nimg2 = FFTW.ifft(new2)

   real.(norm.(nimg1)), real.(norm.(nimg2))
end

#img1, img2 = genPeakFreqsImgs(df, 1)

#println(minimum(img1), ", ", maximum(img1))
#println(minimum(img2), ", ", maximum(img2))

function averageAngle(df, index)
   oldPos = df[!, "lp"][index]
   newPos = df[!, "sp"][index]
   vertNum = df[!, "vertNum"][index]
   v = Array(df[!, "v"][index])
   numRevPairs = df[!, "numRevPairs"][index]
   width = Int32(floor((numRevPairs - 1) * 2))
 
   hdiff1 = genSpacings(v, oldPos, vertNum, horzOnly)[1:width,1:width]
   hdiff2 = genSpacings(v, newPos, vertNum, horzOnly)[1:width,1:width]

   fft1 = FFTW.fft(hdiff1)
   nfft1 = norm.(fft1)
   fft2 = FFTW.fft(hdiff2)
   nfft2 = norm.(fft2)

   max1, _ = findmax(nfft1)
   max2, _ = findmax(nfft2)

   peaks1 = findall(x ->x > max1 ^ 0.5, nfft1)
   peaks2 = findall(x ->x > max2 ^ 0.5, nfft2)

   tan1 = zeros(Float32, size(peaks1))
   for i in 1:size(peaks1)[1]
      tan1[i] = peaks1[i][1] / peaks1[i][2]
   end

   tan2 = zeros(Float32, size(peaks2))
   for i in 1:size(peaks2)[1]
      tan2[i] = peaks2[i][1] / peaks2[i][2]
   end

   tansum1 = 0.0
   for i in 1:size(tan1)[1]
      tansum1 += tan1[i]
   end
   tanavg1 = tansum1 / size(tan1)[1]

   tansum2 = 0.0
   for i in 1:size(tan2)[1]
      tansum2 += tan2[i]
   end
   tanavg2 = tansum2 / size(tan2)[1]

   return atan(tanavg1), atan(tanavg2)
end

#ang1, ang2 = averageAngle(df, 1)

function showVariation(img, ext = 200)
   rimg = real.(img)
   minVal = minimum(rimg)
   maxVal = maximum(rimg)
   plot(Gray.((rimg[1:ext,1:ext] .- minVal) ./ (maxVal - minVal)))
end


circs = [400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]
horzs = Vector(2020:1:2030)
engFreqs = [1700, 3200, 4500, 8100]

for circ in circs
   for horz in horzs
      for engFreq in engFreqs
         d = Dict{String,Any}(["circ" => circ, "horz" => horz, "engFreq" => engFreq])
         df2 = df |> @filter(_.horz == horz && _.circ == circ && _.engFreq == engFreq) |> DataFrame
         img1, img2 = genPeakFreqsImgs(df2, 1)
         ang1, ang2 = averageAngle(df2, 1)
         fulld = Dict{String,Any}(copy(d))
         fulld["img1"] = img1
         fulld["img2"] = img2
         fulld["ang1"] = ang1
         fulld["ang2"] = ang2
         @tagsave(datadir("Analysis1", savename(d, "jld2")), fulld)
      end
   end
end
