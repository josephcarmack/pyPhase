'''A very simple phase-field code that numerically solves the Cahn-Hilliard equation.'''



from __future__ import division
import numpy as np
import os



class InputParams(object):
    '''Input file parser.
    
    Phase-field simulation input parameters:
    
        nx: an integer representing the number of simulation nodes in the x-direction
        ny: an integer representing the number of simulation nodes in the y-direction
        nz: an integer representing the number of simulation nodes in the z-direction
        simSteps: an integer representing the number of simulation steps to be taken 
        outputs: an integer representing the number of simulation outputs 
        dt: a float representing the simulation time step
        delta: a float representing the uniform grid spacing in x,y,z-directions
        c0: a float representing the initial average concentration
        M: a float representing the constant Cahn-Hilliard equation mobility
        kap: a float representing the Cahn-Hilliard equation gradient parameter
        a2: a float representing the quadratic coefficient in the Landau free energy
        a4: a float representing the quartic coefficient in the Landau free energy

    Other Class Attributes:

        outputDir: a string representing the directory in which to write the simulations vtk output files
        dataOverwrite: a string that should be "yes" or "no" indicating whether to overwrite data in an
                       already existing output directory.
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
                'simSteps':int(0),
                'outputs':int(0),
                'dt':float(0.0),
                'delta':float(0.0),
                'c0':float(0.0),
                'M':float(0.0),
                'kap':float(0.0),
                'a2':float(0.0),
                'a4':float(0.0),
                'outputDir':'',
                'vtkTagName':'',
                'dataOverwrite':''
                }

        # parse input file
        print('\n\nParsing the input file...')
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
        print('\nParsed the following simulation parameters:\n')
        for par in sorted(self.simPar.keys()):
            print('\t'+par+' = '+str(self.simPar[par]))



class Simulation(object):
    '''Run a simple Phase-field simulation.
    
    Simulation attributes:
        inputFileParser: object of class type InputParams that parses the input file
        c: object of type ScalarField that represents a concentration order parameter
    '''


    def __init__(self,inputFile):
        '''Return a Simulation object for running a pyPhase simulation.'''
        # parse input file
        self.inputFileParser = InputParams(inputFile)
        self.simPar = self.inputFileParser.simPar
        # get simulation parameters
        self.nx = self.simPar['nx']
        self.ny = self.simPar['ny']
        self.nz = self.simPar['nz']
        self.c0 = self.simPar['c0']
        self.a2 = self.simPar['a2']
        self.a4 = self.simPar['a4']
        self.M = self.simPar['M']
        self.dt = self.simPar['dt']
        self.kap = self.simPar['kap']
        self.delta = self.simPar['delta']
        # define order parameter and chemical potential array
        self.dim = (self.nx,self.ny,self.nz)
        self.orderParam = np.zeros(self.dim,np.float64)
        self.chemPot = np.zeros(self.dim,np.float64)
        # define order parameter and checmical potential Fourier space arrays
        self.ordParKspc = np.zeros(self.dim)
        self.chemPotKspc = np.zeros(self.dim)
        self.nxyz = self.nx*self.ny*self.nz


    def _initFields(self):
        '''Initialize the order parameter with small random noise about c0'''
        print('Initializing the order parameter with small random noise about c0...')
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    self.orderParam[i][j][k] = self.c0 + np.random.uniform(-0.1,0.1)


    def _writeVtkOutput(self,arr,timeStamp):
        '''Writes 3-dimensional data to a vtk output file with associated time stamp.
        
        Input Arguments:
            
            arr: a 3-dimensional numpy array containing order parameter data
            timeStamp: an integer representing the output's simulation time stamp
        '''
        fileName = self.simPar['outputDir']+'/'+self.simPar['vtkTagName']+'_'+str(timeStamp)+'.vtk'
        delta = self.simPar['delta']
        with open(fileName,'w') as fout:
            # write vtk header
            fout.write('# vtk DataFile Version 3.1\n')
            fout.write('VTK file containing grid data\n')
            fout.write('ASCII\n')
            fout.write('\n')
            fout.write('DATASET STRUCTURED_POINTS\n')
            fout.write('DIMENSIONS     '+str(arr.shape[0])+'     '+str(arr.shape[1])+'     '+str(arr.shape[2])+'\n')
            fout.write('ORIGIN 1 1 1\n')
            fout.write('SPACING     '+str(delta)+'     '+str(delta)+'     '+str(delta)+'\n')
            fout.write('\n')
            fout.write('POINT_DATA     ' + str(arr.size)+'\n')
            fout.write('SCALARS '+self.simPar['vtkTagName']+' float\n')
            fout.write('LOOKUP_TABLE default\n')

            # write data
            for i in range(self.simPar['nx']):
                for j in range(self.simPar['ny']):
                    for k in range(self.simPar['nz']):
                        fout.write(str(arr[i][j][k])+'\n')


    def initialize(self):
        '''Initialize order parameter and write initial output'''
        print('\nInitializing the simulation:\n')
        self._initFields()
        outDir = self.simPar['outputDir']
        # create output directory
        if os.path.isdir(outDir):
            if self.simPar['dataOverwrite'] == 'yes':
                for f in os.listdir(outDir):
                    filePath = os.path.join(outDir,f)
                    if os.path.isfile(filePath):
                        os.unlink(filePath)
                self._writeVtkOutput(self.orderParam,0)
            else:
                print('data not intended to be overwritten still exists in the "'+outDir+'" directory')
                raise ValueError
        else:
            os.makedirs(outDir)
            self._writeVtkOutput(self.orderParam,0)


    def _updateOrderParam(self,K2):
        '''Compute the time update for the order parameter using the Cahn-Hilliard equation.
        
        Input Arguments:

            K2: numpy array representing the dot product of the discrete frequency vector at each grid node
        '''
        # calculate the chemical potential in real space
        self.chemPot = self.a2*self.orderParam + self.a4*self.orderParam**3
        # transform order parameter and chemical potential to K-space
        self.ordParKspc = np.fft.fftn(self.orderParam)/self.nxyz
        self.chemPotKspc = np.fft.fftn(self.chemPot)/self.nxyz
        # update orderParmeter using semi-implicit Fourier method in K-space
        self.ordParKspc = (self.ordParKspc - self.M*K2*self.chemPotKspc)/(1 + self.M*self.kap*self.dt*K2**2)
        # transform back from K-space to real space
        self.orderParam = np.real(np.fft.ifftn(self.ordParKspc)*self.nxyz)


    def execute(self):
        '''Execute the simulation.'''

        # calculate output interval
        simSteps = self.simPar['simSteps']
        outputs = self.simPar['outputs']
        if simSteps % outputs is not 0:
            print('invalid number of outputs. outputs must evenly divide into simSteps.')
            raise ValueError
        outInterval = simSteps/outputs
        # create fourier space frequency matrix
        kx = np.fft.fftfreq(self.nx)*2.0*np.pi/self.delta
        ky = np.fft.fftfreq(self.ny)*2.0*np.pi/self.delta
        kz = np.fft.fftfreq(self.nz)*2.0*np.pi/self.delta
        KX,KY,KZ = np.meshgrid(kx,ky,kz)
        K2 = KX**2 + KY**2 + KZ**2
        # run simulation for defined number of steps
        for step in range(1,simSteps+1):

            # write output
            if step % outInterval == 0:
                self._writeVtkOutput(self.orderParam,step)

            # update order parameter
            self._updateOrderParam(K2)
