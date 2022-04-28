# Use optimization to create accurate model of linear motor response
# overall simulation occurs at 10kS/sec, the controller's sampling rate

using CSV
using DataFrames
using Optim
using Plots
gr(size=[1024,768])

in1 = CSV.read("../linearMotor/posA.csv", DataFrame, header=false)
in1A = Array(in1[!,"Column1"])
out1 = CSV.read("../linearMotor/fdbkA.csv", DataFrame, header=false)
out1A = Array(out1[!, "Column1"])
in2 = CSV.read("../linearMotor/posB.csv", DataFrame, header=false)
in2A = Array(in2[!,"Column1"])
out2 = CSV.read("../linearMotor/fdbkB.csv", DataFrame, header=false)
out2A = Array(out2[!,"Column1"])

# First model - force, mass and damping

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

function controlLaw2(params, command, feedback, cState)
   # cState = [lastCommand, lastPosition, currIntegral]
   Kv, Kpos, Kp, Ki = params[4:end]

   lastOutput, lastCommand, lastPosition, currIntegral = cState

   aVel = feedback - lastPosition
   cVel = command - lastCommand
   velError = Kv * cVel - aVel
   posError = command - feedback
   totalError = (Kpos / 4096) * posError + velError
   newIntegral = currIntegral + totalError / 10000.0
   output = (Kp / 2^23) * totalError + (Ki / 2^23) * newIntegral

   return [output, command, feedback, newIntegral]
end

function nullFilter(sample)
   sample
end

function makeBiquad(num, den, steadyState)
   mem = steadyState * ones(Float64, 2)
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

# state for carriage should just be position and velocity
# acceleration will be instantaneous based on control
# of forcer. Simulation will be run at 10kS/s.
# Position and velocity will be modeled using constant force
# over the sample time, a linear change in velocity and
# parabolic change in position.

# to simulate stiction, simulation could hold controlLaw output
# at zero until actual position shows carriage started moving.

# function for Optim to optimize must run full simulation
# of carriage motion to return the error. It will accept a Vector
# of [force/mass, damping] and return total error between
# simulated carriage path and measured carriage path

function fullSim(input, startPos, params, contLaw, filt, carrSim)
   cState = [0.0, input[1], startPos, 0.0] # lastCommand, lastOutput, lastPosition, currIntegral
   carrState = [startPos, 0.0]             # position, velocity
   output = zeros(Float32, size(input))
   output[1] = startPos
   for i in 2:size(input)[1]-1
      cState = contLaw(input[i], output[i], cState)
      filtOut = filt(cState[1])
      carrState = carrSim(params, carrState, filtOut)
      output[i+1] = round(carrState[1])
      cState[3] = output[i]
   end
   output
end

function fullSim2(input, startPos, params, contLaw2, filt, carrSim2)
   cState = [0.0, input[1], startPos, 0.0] # lastCommand, lastOutput, lastPosition, currIntegral
   carrState = [startPos, 0.0]             # position, velocity
   output = zeros(Float32, size(input))
   output[1] = startPos
   for i in 2:size(input)[1]-1
      cState = contLaw2(params[4:end], input[i], output[i], cState)
      filtOut = filt(cState[1])
      carrState = carrSim2(params, carrState, filtOut)
      output[i+1] = round(carrState[1])
      cState[3] = output[i]
   end
   output
end

function genSim(input, measured, contLaw, filt, carrSim)
   function (params)
      output = fullSim(input, measured[1], params, contLaw, filt, carrSim)

      sum = 0.0
      for i in 1:size(measured)[1]
         err = measured[i] - output[i]
         sum += err^2
      end
      sum^0.5
   end
end

function genSim2(input, measured, contLaw2, filt, carrSim2)
   function (params)
      output = fullSim2(input, measured[1], params, contLaw2, filt, carrSim2)

      sum = 0.0
      for i in 1:size(measured)[1]
         err = measured[i] - output[i]
         sum += err^2
      end
      sum^0.5
   end
end

numerator = [0.019789582118392, 0.039579164236784, 0.019789582118392]
denominator = [-1.56450402736664, 0.643662333488464]

filt500 = makeBiquad(numerator, denominator, in1A[1])

testFunc = genSim(in1A, out1A, controlLaw, filt500, carrSim)

results = optimize(testFunc, [2.0e9, 0.1])

params = Optim.minimizer(results)

filt500 = makeBiquad(numerator, denominator, in1A[1])
simOut = fullSim(in1A, out1A[1], params, controlLaw, filt500, carrSim)

plot(simOut)
plot!(out1A)
plot!(in1A)

filt500 = makeBiquad(numerator, denominator, in2A[1])
testFunc2 = genSim(in2A, out2A, controlLaw, filt500, carrSim)
result2 = optimize(testFunc, [2.0e9, 0.3])
param2 = Optim.minimizer(result2)
filt500 = makeBiquad(numerator, denominator, in2A[1])
sim2Out = fullSim(in2A, out2A[1], param2, controlLaw, nullFilter, carrSim)
plot(sim2Out)
plot!(out2A)

function carrSimStiction(params, carrState, command)
   sampleTime = 1 / 10000
   if abs(carrState[2]) < 0.00001 && abs(command) < abs(params[3])
      return [carrState[1], 0.0]
   end:wait

   accel = command * params[1]
   newVelocity = (1 - params[2]) * carrState[2] + accel * sampleTime
   newPosition = carrState[1] + (1 - params[2]) * carrState[2] * sampleTime + accel * sampleTime^2
   [newPosition, newVelocity]

end

filt500 = makeBiquad(numerator, denominator, in1A[1])
testFunc3 = genSim(in1A, in1A, controlLaw, filt500, carrSimStiction)
result3 = optimize(testFunc3, [1.0e7, 0.1, 1000.0])
param3 = Optim.minimizer(result3)
filt500 = makeBiquad(numerator, denominator, in1A[1])
sim3out = fullSim(in1A, out1A[1], param3, controlLaw, nullFilter, carrSimStiction)
plot(sim3out)
plot!(out1A)
plot!(in1A)

testFunc4 = genSim2(in1A, out1A, controlLaw2, nullFilter, carrSim)
results4 = optimize(testFunc4, [100.0, 1.0, 0.0, 1.0, 20.0, 200000.0, 1000.0])
param4 = Optim.minimizer(results4)
sim4out = fullSim2(in1A, out1A[1], param4, controlLaw2, nullFilter, carrSimStiction)

shift = 110
in1Amod = in1A[1:end-shift]
out1Amod = out1A[shift+1:end]
filt500 = makeBiquad(numerator, denominator, in1Amod[1])
testFunc5 = genSim(in1Amod, out1Amod, controlLaw, filt500, carrSim)
result5 = optimize(testFunc5, [2e9, 1.0])
param5 = Optim.minimizer(result5)
filt500 = makeBiquad(numerator, denominator, in1Amod[1])
sim5out = fullSim(in1Amod, out1Amod[1], param5, controlLaw, filt500, carrSim)
plot(sim5out)
plot!(in1Amod)
plot!(out1Amod)

filt500 = makeBiquad(numerator, denominator, in1A[1])
testFunc6 = genSim(in1A, in1A, controlLaw, filt500, carrSim)
result6 = optimize(testFunc6, [2.0e9, 0.15])
param6 = Optim.minimizer(result6)
filt500 = makeBiquad(numerator, denominator, in1A[1])
sim6Out = fullSim(in1A, in1A[1], param6, controlLaw, filt500, carrSim)
plot(sim6Out)
plot!(in1A)