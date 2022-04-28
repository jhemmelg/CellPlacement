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

# Controller with fixed gains set to machine standard gains
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

# Make a biquad filter function
# A new filter function needs to be made
# for each use as the memory persists
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

# Create a simulation function suitable for optimization
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

# coefficients for a 500Hz low-pass filter
numerator = [0.019789582118392, 0.039579164236784, 0.019789582118392]
denominator = [-1.56450402736664, 0.643662333488464]

filt500 = makeBiquad(numerator, denominator, in1A[1])
testFunc = genSim(in1A, out1A, controlLaw, filt500, carrSim)
results = optimize(testFunc, [2.0e9, 0.1])
params = Optim.minimizer(results)
filt500 = makeBiquad(numerator, denominator, in1A[1])
simOut = fullSim(in1A, out1A[1], params, controlLaw, filt500, carrSim)
plot(in1A, legend=:topleft, label="Input", title="Carriage Simulation")
plot!(simOut, label="Sim")
plot!(out1A, label="Actual")
println(params)

params = [2e9, 0.1]
filt500 = makeBiquad(numerator, denominator, in1A[1])
simOut = fullSim(in1A, out1A[1], params, controlLaw, filt500, carrSim)
plot(in1A, legend=:topleft, label="Input", title="Carriage Simulation")
plot!(simOut, label="Sim")
plot!(out1A, label="Actual")
println(params)
#filt500 = makeBiquad(numerator, denominator, in2A[1])
#testFunc2 = genSim(in2A, out2A, controlLaw, filt500, carrSim)
#result2 = optimize(testFunc, [2.0e9, 0.3])
#param2 = Optim.minimizer(result2)
#filt500 = makeBiquad(numerator, denominator, in2A[1])
#sim2Out = fullSim(in2A, out2A[1], param2, controlLaw, nullFilter, carrSim)
#plot(sim2Out)
#plot!(out2A)
#plot!(in2A)
