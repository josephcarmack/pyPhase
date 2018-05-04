from pyPhase import Simulation as sim

# create a simulation with parameters from an input file
mySim = sim('simparams.dat')

# initialize the simulation
mySim.initialize()

# run the simulation
mySim.execute()
