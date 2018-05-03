'''A very simple phase-field code that numerically solves the Cahn-Hilliard equation.'''

import numpy as np

class InputParams(object):
    '''Input file parser.
    
    Phase-field simulation input parameters:
    
        nx: an integer representing the number of simulation nodes in the x-direction
        ny: an integer representing the number of simulation nodes in the y-direction
        nz: an integer representing the number of simulation nodes in the z-direction
        dt: a float representing the simulation time step
        delta: a float representing the uniform grid spacing in x,y,z-directions
        c0: a float representing the initial average concentration
        M: a float representing the constant Cahn-Hilliard equation mobility
        kap: a float representing the Cahn-Hilliard equation gradient parameter
        a0: a float representing the constant coefficient in the Landau free energy
        a1: a float representing the linear coefficient in the Landau free energy
        a2: a float representing the quadratic coefficient in the Landau free energy
        a3: a float representing the cubic coefficient in the Landau free energy
        a4: a float representing the quartic coefficient in the Landau free energy

    Other Class Attributes:

        outputDir: a string representing the directory in which to write the simulations vtk output files
    '''
    def __init__(self,inputFile):
        '''Return an InputParams object with simulation parameters read from an input file.

        Input Arguments:

            inputFile: a string represting the name of the input file with relative path
                       (i.e. ./relative/path/to/thisIsMyInputFile.txt). Assumes input file
                       has following format: <parameterName> = <value/string>
        '''
        # store simulation parameters in a dictionary
        self.simPar = {
                'nx':int(0),
                'ny':int(0),
                'nz':int(0),
                'dt':float(0.0),
                'delta':float(0.0),
                'c0':float(0.0),
                'M':float(0.0),
                'kap':float(0.0),
                'a0':float(0.0),
                'a1':float(0.0),
                'a2':float(0.0),
                'a3':float(0.0),
                'a4':float(0.0),
                'outputDir':''
                }

        # parse input file
        print('Parsing the input file...')
        with open(inputFile) as inFile:
            # get all non-blank lines from the input file
            lines = (line.rstrip() for line in inFile) # all lines including blank ones
            lines = (line for line in lines if line) # non blank lines
            # assume line has format <parameter> = <value> and parse
            for line in lines:
                param = line.split('=')[0]
                value = line.split('=')[1]
                # remove white space from param and value
                param = param.replace(' ','')
                value = value.replace(' ','')
                # check if param is valid
                if param not in self.simPar.keys():
                    print(str(param)+' is not a valid simulation parameter')
                    raise ValueError
                self.simPar[param] = type(self.simPar[param])(value)

        # print out the values parsed from the input file
        print('Parsed the following simulation parameters:')
        for par in sorted(self.simPar.keys()):
            print('\t'+par+' = '+str(self.simPar[par]))

class ScalarField(object):
    '''An implementation of a simple scalar phase-field order parameter'''
    def __init__(self,simParams):
        '''Returns a ScalarField object initialized with parameters stored in simParams.

        Input Arguments:
            
            simParams: a dictionary containing the simulation parameters.

        Class Attributes:

            data: a numpy array used for storing the ScalarField data
        '''
        dim = (simParams['nx'],simParams['ny'],simParams['nz'])
        self.data = np.zeros(dim)

class Simulation(object):
    '''Run a simple Phase-field simulation.
    
    Simulation attributes:
        inputFileParser: object of class type InputParams that parses the input file
        c: object of type ScalarField that represents a concentration order parameter
    '''
    def __init__(self,inputFile):
        '''Return a Simulation object for running a pyPhase simulation.'''
        self.inputFileParser = InputParams(inputFile)
        self.simPar = self.inputFileParser.simPar
        self.orderParam = ScalarField(self.simPar)

