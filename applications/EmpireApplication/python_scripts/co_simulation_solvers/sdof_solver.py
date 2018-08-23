from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
import numpy as np
import json
import os

def CreateSolver(cosim_solver_settings, level):
    return SDoFSolver(cosim_solver_settings, level)

class SDoFSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level, buffer_size=3):
        super(SDoFSolver, self).__init__(cosim_solver_settings, level)

        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        self.mass = parameters["system_parameters"]["mass"]
        self.stiffness = parameters["system_parameters"]["stiffness"]
        self.damping = parameters["system_parameters"]["damping"]

        self.alpha_m = parameters["time_integration_parameters"]["alpha_m"]
        self.alpha_f = parameters["time_integration_parameters"]["alpha_f"]
        self.delta_t = parameters["time_integration_parameters"]["time_step"]

        self.initial_displacement = parameters["initial_values"]["displacement"]
        self.initial_velocity = parameters["initial_values"]["velocity"]

        self.force = parameters["boundary_conditions"]["external_load"]

        #calculate initial acceleration
        factor = self.force - self.stiffness * self.initial_displacement
        self.initial_acceleration = (1/self.mass) * factor

        beta = 0.25 * (1- self.alpha_m + self.alpha_f)**2
        gamma =  0.50 - self.alpha_m + self.alpha_f

        self.LHS = np.array([[1.0, 0.0, -self.delta_t**2 * beta],
                              [0.0, 1.0, -self.delta_t * gamma],
                              [(1-self.alpha_f)*self.stiffness,
                               (1-self.alpha_f) * self.damping,
                               (1-self.alpha_m) * self.mass]])

        self.RHS_matrix = np.array([[1.0, self.delta_t, self.delta_t**2 * (0.5 - beta)],
                                     [0.0, 1.0, self.delta_t*(1-gamma)],
                                     [-self.alpha_f * self.stiffness,
                                      -self.alpha_f * self.damping,
                                      -self.alpha_m * self.mass]])

        self.buffer_size = buffer_size

    def Initialize(self):
        self.x = np.zeros((3, self.buffer_size))

        initial_values = np.array([self.initial_displacement,
                                   self.initial_velocity,
                                   self.initial_acceleration])
        self.dx = initial_values

        #x and dx contain: [displacement, velocity, acceleration]

        if os.path.isfile("results_sdof.txt"):
            os.remove("results_sdof.txt")

        #apply external load as an initial impulse
        self.load_vector = np.array([0,
                                     0,
                                     self.force])

    def OutputSolutionStep(self):
        with open("results_sdof.txt", "a") as results_sdof:
            #outputs displacements
            results_sdof.write(str(self.time) + "\t" + str(self.dx[0]) + "\n")

    def AdvanceInTime(self, current_time):
        # similar to the Kratos CloneTimeStep function
        # advances values along the buffer axis (so rolling columns) using numpy's roll
        self.x = np.roll(self.x,1,axis=1)
        # overwriting at the buffer_idx=0 the newest values
        #buffer_idx = 0
        #self.x[:,buffer_idx] = self.dx
        self.x[:,0] = self.dx

        self.time = current_time + self.delta_t
        return self.time

    def SolveSolutionStep(self):
        ## PMT: ToDo: check with Andreas what this intends to do
        #b = self.RHS_matrix @ self.x[:,0]
        b = np.dot(self.RHS_matrix,self.x[:,0])

        #external load only for testing
        b += self.load_vector
        self.dx = np.linalg.solve(self.LHS, b)

    def GetBufferSize(self):
        return self.buffer_size

    def GetDeltaTime(self):
        return self.delta_t

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            return self.x[:,buffer_idx][0]
        elif identifier == "VELOCITY":
            return self.x[:,buffer_idx][1]
        elif identifier == "ACCELERATION":
            return self.x[:,buffer_idx][2]
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            self.x[:,buffer_idx][0] = value
        elif identifier == "VELOCITY":
            self.x[:,buffer_idx][1] = value
        elif identifier == "ACCELERATION":
            self.x[:,buffer_idx][2] = value
        else:
            raise Exception("Identifier is unknown!")

    def SetData(self, identifier, data):
        if identifier == "LOAD":
            # last index is the external force
            self.load_vector[-1] = data
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            # maybe use buffer index
            self.SetSolutionStepValue("DISPLACEMENT", data,0)
        else:
            raise Exception("Identifier is unknown!")

    def GetData(self, identifier):
        if identifier == "LOAD":
            # last index is the external force
            return self.load_vector[-1]
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            return self.GetSolutionStepValue("DISPLACEMENT",0)
        else:
            raise Exception("Identifier is unknown!")

    def _GetIOName(self):
        return "sdof"

    def _Name(self):
        return self.__class__.__name__
