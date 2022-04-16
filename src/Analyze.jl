function genSpacings(arrV, arrH, vert, meth)
   revs = size(arrH)[1] / vert
   advance = Int32(floor(vert))
   height = vert * 2
   retVal = zeros(Float32, Int32(revs - 1), Int32(height))
   for col in 1:revs-1
      for row in 1:vert*2
         i = (col - 1) * height + row
         left = floor((i + 1) / 2)
         right = left + advance + rem(i + 1, 2)
         retVal[row,col] = meth(arrH[left], arrV[left], arrH[right], arrV[left])
      end
   end
   retVal
end

function horzOnly(hLeft, vLeft, hRight, vRight)
   hRight - hLeft
end

function distance(hLeft, vLeft, hRight, vRight)
   ((hRight - hLeft)^2 + (vRight - vLeft)^2)^0.5
end