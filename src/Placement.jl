include("constants.jl")

function actualVert(vert, circ)
   cellsPerRev = circ * 1000.0 / vert
   actualVert  = round(cellsPerRev * 2)
   if (rem(actualVert, 2) == 0)
      actualVert = actualVert / 2 + 0.5
   else
      actualVert = actualVert / 2
   end
   actualVert
end

function horzNum(horizontalSpacing)
   exactStepsPerRev = horizontalSpacing / OhioStep
   idealHorzNumber = pulsesPerRev / exactStepsPerRev
   horzNum = round(idealHorzNumber)
end

function revWidth(horzNum)
   (pulsesPerRev / horzNum) * OhioStep
end

function vertSpacing(vertNum, circ)
   circ * 1000 / vertNum
end

function vertPositions(vert, circ, numRevPairs)
   cellsPerTwoRev = Int32(actualVert(vert, circ) * 2)
   actualVerticalSpacing = (circ * 1000.0) / (cellsPerTwoRev / 2)
   vertPositions = Array(range(0.5, numRevPairs * cellsPerTwoRev - 0.5, length = numRevPairs * cellsPerTwoRev)) .* actualVerticalSpacing
   vertPositions = mod.(vertPositions, circ)
end

function idealHorzPositions(revWidth, cellsPerRev, numRevPairs)
   horizontalDistancePerCell = revWidth / cellsPerRev
   idealHorzPositions = Array(range(0.5, 2 * cellsPerRev * numRevPairs - 0.5, length = Int32(2 * cellsPerRev * numRevPairs)) .* horizontalDistancePerCell)
end

function steppedHorzPositions(horzNum, cellsPerRev, numRevPairs)
   stepsPerRev = pulsesPerRev / horzNum
   cellsPerStep = cellsPerRev / stepsPerRev
   steps = floor.((1:2*cellsPerRev*numRevPairs) ./ cellsPerStep) .+ 1
   steppedHorzPositions = steps .* OhioStep
end

function pointOneHorzPositions(revWidth, cellsPerRev, numRevPairs)
   newStep = 0.1
   stepsPerRev = revWidth / newStep
   cellsPerStep = cellsPerRev / stepsPerRev
   cellsPerTwoRev = cellsPerRev * 2
   steps = floor.((1:cellsPerTwoRev*numRevPairs) ./ cellsPerStep) .+ 1
   horzPositions = steps .* newStep
   return horzPositions
end

function idealPlacement(circ, horz, vert, numRevPairs)
   vertPos = vertPositions(vert, circ, numRevPairs)
   cellsPerRev = actualVert(vert, circ)
   horzPos = idealHorzPositions(horz, cellsPerRev, numRevPairs)
   return (vertPos, horzPos)
end

function largeStepPlacement(circ, horz, vert, numRevPairs)
   vertPos = vertPositions(vert, circ, numRevPairs)
   cellsPerRev = actualVert(vert, circ)
   stepHorzNum = horzNum(horz)
   horzPos = steppedHorzPositions(stepHorzNum, cellsPerRev, numRevPairs)
   return (vertPos, horzPos)
end

function smallStepPlacement(circ, horz, vert, numRevPairs)
   vertPos = vertPositions(vert, circ, numRevPairs)
   cellsPerRev = actualVert(vert, circ)
   horzPos = pointOneHorzPositions(horz, cellsPerRev, numRevPairs)
   return (vertPos, horzPos)
end
