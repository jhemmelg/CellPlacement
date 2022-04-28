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
   exactStepsPerRev = horizontalSpacing / 2.0 / OhioStep
   idealHorzNumber = pulsesPerRev / exactStepsPerRev
   horzNum = round(idealHorzNumber)
end

function revWidthByHorzNum(horzNum)
   (pulsesPerRev / horzNum) * OhioStep
end

function vertSpacing(vertNum, circ)
   circ * 1000 / vertNum
end

function vertPositions(vert, circ, numRevPairs)
   cellsPerTwoRev = Int32(actualVert(vert, circ) * 2)
   actualVerticalSpacing = (circ * 1000.0) / (cellsPerTwoRev / 2)
   vertPositions = Array(range(0.5, numRevPairs * cellsPerTwoRev - 0.5, length = numRevPairs * cellsPerTwoRev)) .* actualVerticalSpacing
   vertPositions = mod.(vertPositions, circ * 1000.0)
end


# calculate horizontal position using simulation of carriage horizontal movement
# must generate command array at controllerRate, not head frequency
# run carriage simulation on this data
# then resample at cell centers
#
# Inputs: number of cells around
#         step distance
#         cells per step
#         cell frequency
#         number of two-revs to generate
#         function to run simulation of carriage motion

function posWithCarriageSim(numCellsAround, stepDist, cellsPerStep, engFreq, num2Revs, simFunc)
   println("cellsPerStep = ", cellsPerStep)
   totalCells = Int32(numCellsAround * num2Revs * 2)
   println("totalCells = ", totalCells)
   totalTime = totalCells / engFreq
   totalTicks = Int32(round(totalTime * controllerRate))
   println("totalTicks = ", totalTicks)
   command = zeros(Float32, totalTicks)
   ticksPerStep = cellsPerStep * (controllerRate / engFreq)
   println("ticksPerStep = ", ticksPerStep)
   println("stepDist = ", stepDist)
   for i in 1:totalTicks
      command[i] = stepDist * round(i / ticksPerStep)
   end
   #return command
   simOut = simFunc(command)
   simOut
   pos = zeros(Float32, totalCells)
   for i in 1:totalCells
      time = i / engFreq
      index = Int32(floor(time * controllerRate))
      #println(i, " ", time, " ", index)
      pos[i] = simOut[index]
   end
   pos
end

function steppedHorzPosSim(horzNum, cellsPerRev, numRevPairs, engFreq, simFunc)
   println("horzNum = ", horzNum)
   stepsPerRev = pulsesPerRev / horzNum 
   println("stepsPerRev = ", stepsPerRev)
   cellsPerStep = cellsPerRev / stepsPerRev
   println("cellsPerStep = ", cellsPerStep)
   posWithCarriageSim(cellsPerRev, OhioStep, cellsPerStep, engFreq, numRevPairs, simFunc)
end

#function vertPositions(vert, circ, numRevPairs)
#   cellsPerTwoRev = Int32(actualVert(vert, circ) * 2)
#   actualVerticalSpacing = (circ * 1000.0) / (cellsPerTwoRev / 2)
#   vertPositions = Array(range(0.5, numRevPairs * cellsPerTwoRev - 0.5, length = numRevPairs * cellsPerTwoRev)) .* actualVerticalSpacing
#   vertPositions = mod.(vertPositions, circ)
#end

function idealHorzPositions(horz, cellsPerRev, numRevPairs)
   horizontalDistancePerCell = revWidthByHorzNum(horzNum(horz)) / cellsPerRev
   idealHorzPositions = Array(range(0.5, 2 * cellsPerRev * numRevPairs - 0.5, length = Int32(2 * cellsPerRev * numRevPairs)) .* horizontalDistancePerCell)
end

function idealPlacement(circ, horz, vert, numRevPairs)
   vertPos = vertPositions(vert, circ, numRevPairs)
   cellsPerRev = actualVert(vert, circ)
   horzPos = idealHorzPositions(horz, cellsPerRev, numRevPairs)
   return (vertPos, horzPos)
end

function controlLaw(command, feedback, cState)
   # cState = [lastCommand, lastPosition, currIntegral]
   Kv = 1
   Kpos = 20
   Kp = 200000
   Ki = 1000

   lastOutput, lastCommand, lastPosition, currIntegral = cState

   aVel = feedback - lastPosition
   cVel = command - lastCommand
   velError = Kv * cVel - aVel
   posError = command - feedback
   totalError = (Kpos / 4096) * posError + velError
   newIntegral = currIntegral + totalError
   output = (Kp / 2^23) * totalError + (Ki / 2^23) * newIntegral

   return [output, command, feedback, newIntegral]
end

function nullFilter(sample)
   sample
end

function makeBiquad(num, den)
   mem = zeros(Float64, 2)
   function (sample)
      w = sample - den[1] * mem[1] - den[2] * mem[2]
      y = num[1] * w + num[2] * mem[1] + num[3] * mem[2]
      mem[2] = mem[1]
      mem[1] = w
      y
   end
end

function carrSim(params, carrState, command)
   # params is [f/m, damping]
   # carrState is [position, velocity]
   # command is the applied force
   # simulation runs at 10kS/sec
   sampleTime = 1 / 10000
   accel = command * params[1]
   newVelocity = (1 - params[2]) * carrState[2] + accel * sampleTime
   newPosition = carrState[1] + (1 - params[2]) * carrState[2] * sampleTime + accel * sampleTime^2
   [newPosition, newVelocity]
end

function fullSim(input, startPos, params, contLaw, filt, carrSim)
   cState = [0.0, input[1], startPos, 0.0] # lastCommand, lastOutput, lastPosition, currIntegral
   carrState = [startPos, 0.0]             # position, velocity
   output = zeros(Float32, size(input))
   output[1] = startPos
   for i in 2:size(input)[1]-1
      cState = contLaw(input[i], output[i], cState)
      filtOut = filt(cState[1])
      carrState = carrSim(params, carrState, filtOut)
      output[i+1] = carrState[1]
      cState[3] = output[i]
   end
   output
end

function CarriageSimFunc(command)
   params = [2.225e9, 1.5]
   numerator = [0.019789582118392, 0.039579164236784, 0.019789582118392]
   denominator = [-1.56450402736664, 0.643662333488464]
   filt500 = makeBiquad(numerator, denominator)
   fullSim(command, command[1], params, controlLaw, filt500, carrSim)
end

function largeStepPlacement(circ, horz, vert, engFreq, numRevPairs)
   vertPos = vertPositions(vert, circ, numRevPairs)
   cellsPerRev = actualVert(vert, circ)
   stepHorzNum = horzNum(horz)
   horzPos = steppedHorzPosSim(stepHorzNum, cellsPerRev, numRevPairs, engFreq, CarriageSimFunc)
   return (vertPos, horzPos)
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

function pointOneHorzPosSim(revWidth, cellsPerRev, numRevPairs, engFreq, CarriageSimFunc)
   newStep = 0.1
   stepsPerRev = revWidth / newStep
   countsPerStep = encoderCounts / stepsPerRev
   stepsPerRev = encoderCounts / countsPerStep
   cellsPerStep = cellsPerRev / stepsPerRev
   posWithCarriageSim(cellsPerRev, newStep, cellsPerStep, engFreq, numRevPairs, CarriageSimFunc)
end

function smallStepPlacement(circ, horz, vert, engFreq, numRevPairs)
   vertPos = vertPositions(vert, circ, numRevPairs)
   cellsPerRev = actualVert(vert, circ)
   horzPos = pointOneHorzPosSim(revWidthByHorzNum(horzNum(horz)), cellsPerRev, numRevPairs, engFreq, CarriageSimFunc)
   return (vertPos, horzPos)
end
