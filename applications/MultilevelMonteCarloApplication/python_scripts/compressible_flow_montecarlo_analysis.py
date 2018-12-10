from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressFlow

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
# from analysis_stage import AnalysisStage

# Importing derived classes
from potential_flow_analysis import PotentialFlowAnalysis

# Import pycompss
# from pycompss.api.task import task
# from pycompss.api.api import compss_wait_on
# from pycompss.api.parameter import *

# Import exaqute
# from exaqute.ExaquteTaskPyCOMPSs import *   # to exequte with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to exequte with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3
# get_value_from_remote is the equivalent of compss_wait_on
# in the future, when everything is integrated with the it4i team, putting exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC

# Import Monte Carlo library
import mc_utilities as mc

# Import variables class
from cmlmc_utilities_CLASS import StatisticalVariable

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle

class MonteCarloAnalysis(PotentialFlowAnalysis):
    """Main analyis stage for Monte Carlo simulations"""
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(MonteCarloAnalysis,self).__init__(input_model,input_parameters)
        # self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
    
    def _GetSimulationName(self):
        return "Monte Carlo Analysis"

    def ModifyInitialProperties(self):
        '''Introduce here the stochasticity in the Mach number and the angle of attack'''
        Mach = self.sample[0]
        a_infinity = 340 # [m/s] velocity of sound at infinity
        alpha =  self.sample[1]
        v_norm = Mach * a_infinity
        velocity = [v_norm*np.cos(alpha),v_norm*np.sin(alpha),0]
        boundary_processes = self.project_parameters["processes"]["boundary_conditions_process_list"]       
        for i in range(0,boundary_processes.size()):
            python_module = boundary_processes[i]["python_module"].GetString()
            if python_module == "apply_far_field_process":
                self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["velocity_infinity"].SetVector(velocity)


    
##################################################
######## END OF CLASS MONTECARLOANALYSIS #########
##################################################


'''
function generating the random sample
here the sample has a normal distribution
'''
def GenerateSample():
    sample = []
    mean_Mach = 0.7
    std_deviation_Mach = 0.01
    number_samples = 1
    sample.append(np.random.normal(mean_Mach,std_deviation_Mach,number_samples))
    mean_angle_attack = 0.0 # [rad] = 0 [degrees] airfoil already has 5 degrees
    std_deviation_angle_attack = 0.01
    sample.append(np.random.normal(mean_angle_attack,std_deviation_angle_attack,number_samples))
    print("MACH NUMBER = ",sample[0],"ANGLE ATTACK = ",sample[1])
    if sample[0] >= 1.0 or sample[0] <= 0.0 :
        raise Exception ("Mach randomly computed > 1 or < 0")
    return sample



'''
function evaluating the QoI of the problem
'''
def EvaluateQuantityOfInterest(simulation):
    """here we evaluate the QoI of the problem: the lift coefficient"""
    Q = simulation._GetSolver().main_model_part.GetValue(KratosMultiphysics.FRICTION_COEFFICIENT)
    return Q

'''
function executing the problem
input:
        model       : serialization of the model
        parameters  : serialization of the Project Parameters
output:
        QoI         : Quantity of Interest
'''
@ExaquteTask(returns=1)
def ExecuteTask(pickled_model, pickled_parameters):
    '''overwrite the old model serializer with the unpickled one'''
    model_serializer = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    '''overwrite the old parameters serializer with the unpickled one'''
    serialized_parameters = pickle.loads(pickled_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    sample = GenerateSample()
    simulation = MonteCarloAnalysis(current_model,current_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    return QoI


'''
function executing the problem for sample = 1.0
input:
        model       : serialization of the model
        parameters  : serialization of the Project Parameters
output:
        ExactExpectedValueQoI : Quantity of Interest for sample = [0.7,0.0]
'''
@ExaquteTask(returns=1)
def ExecuteExactTask(pickled_model, pickled_parameters):
    '''overwrite the old model serializer with the unpickled one'''
    model_serializer = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    '''overwrite the old parameters serializer with the unpickled one'''
    serialized_parameters = pickle.loads(pickled_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    sample = [0.7,0.0]
    simulation = MonteCarloAnalysis(current_model,current_parameters,sample)
    simulation.Run()
    ExactExpectedValueQoI = EvaluateQuantityOfInterest(simulation)
    return ExactExpectedValueQoI


'''
function serializing the model and the parameters of the problem
the idea is the following:
first from Model/Parameters Kratos object to StreamSerializer Kratos object
second from StreamSerializer Kratos object to pickle string
third from pickle string to StreamSerializer Kratos object
fourth from StreamSerializer Kratos object to Model/Parameters Kratos object
input:
        parameter_file_name   : path of the Project Parameters file
output:
        serialized_model      : model serialized
        serialized_parameters : project parameters serialized
'''
@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def SerializeModelParameters(parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters
    model = KratosMultiphysics.Model()
    # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    fake_sample = [0.7,0.0]
    simulation = MonteCarloAnalysis(model,local_parameters,fake_sample)
    simulation.Initialize()
    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    serialized_parameters = KratosMultiphysics.StreamSerializer()
    serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)
    # pickle dataserialized_data
    pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_parameters = pickle.dumps(serialized_parameters, 2)
    return pickled_model, pickled_parameters


'''
function computing the relative error between the Multilevel Monte Carlo expected value and the exact expected value
input :
        AveragedMeanQoI       : Multilevel Monte Carlo expected value
        ExactExpectedValueQoI : exact expected value
output :
        relative_error        : relative error
'''
@ExaquteTask(returns=1)
def CompareMean(AveragedMeanQoI,ExactExpectedValueQoI):
    relative_error = abs((AveragedMeanQoI - ExactExpectedValueQoI)/ExactExpectedValueQoI)
    return relative_error


if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg = 'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python montecarlo_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python3 montecarlo_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "/home/kratos105b/DataDisk/MultilevelMonteCarloApplication/CompressiblePotentialFlow/ProjectParameters2.json"

    '''create a serialization of the model and of the project parameters'''
    pickled_model,pickled_parameters = SerializeModelParameters(parameter_file_name)
    print("\n############## Serialization completed ##############\n")

    number_samples = 1
    QoI = StatisticalVariable()
    '''to exploit StatisticalVariable UpdateOnePassMeanVariance function we need to initialize a level 0 in values, mean, sample variance and second moment
    and store in this level the informations'''
    QoI.values = [[] for i in range (1)]
    QoI.mean = [[] for i in range (1)]
    QoI.second_moment = [[] for i in range (1)]
    QoI.sample_variance = [[] for i in range (1)]

    for instance in range (0,number_samples):
        QoI.values[0].append(ExecuteTask(pickled_model,pickled_parameters))

    '''Compute mean, second moment and sample variance'''
    for i_sample in range (0,number_samples):
        QoI.UpdateOnepassMeanVariance(0,i_sample)

    QoI = get_value_from_remote(QoI)
    print("\nMC mean lift coefficient = ",QoI.mean)
    
