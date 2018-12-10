from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()
from python_solver import PythonSolver

def CreateSolver(model, custom_settings):
    return PotentialSolver(model, custom_settings["solver_settings"])

class PotentialSolver(PythonSolver):
    def __init__(self, model, custom_settings):
        self.MoveMeshFlag = False

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"        : "model",                   
            "domain_size"            : 2,
            "solver_type": "potential_flow_solver",
            "problem_type"	         : "incompressible",
            "echo_level": 1,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-9,
            "maximum_iterations": 1,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "calculate_solution_norm" : false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[],
            "no_skin_parts"                : [],
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "linear_solver_settings": {
                    "solver_type": "AMGCL",
                    "max_iteration": 400,
                    "gmres_krylov_space_dimension": 500,
                    "smoother_type":"ilu0",
                    "coarsening_type":"ruge_stuben",
                    "coarse_enough" : 5000,
                    "krylov_type": "lgmres",
                    "tolerance": 1e-9,
                    "verbosity": 0,
                    "scaling": false
            }


        }""")

            # "linear_solver_settings"       : {
            #      "solver_type"     : "SkylineLUFactorizationSolver"
            #   }
         
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        model_part_name = self.settings["model_part_name"].GetString()
        super(PotentialSolver,self).__init__(model, self.settings)

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
            self.solver_imports_model_part = False
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)
            self.solver_imports_model_part = True
        
        self.domain_size = custom_settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DENSITY, 1.225)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.WATER_PRESSURE,2.0)#n_parameter
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TEMPERATURE,0.0)#penalty stress
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.INITIAL_PENALTY,0.0)#penalty kutta
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.LAMBDA, 1.4)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SOUND_VELOCITY, 340.0)
        
                    
        #construct the linear solvers
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of LaplacianSolver finished")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.NEGATIVE_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.POSITIVE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.NEGATIVE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.VELOCITY_INFINITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE_DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.LEVEL_SET_DISTANCE)
        
    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.NEGATIVE_POTENTIAL, self.main_model_part)
        
    def Initialize(self):
        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        move_mesh_flag = False #USER SHOULD NOT CHANGE THIS

        if self.settings["problem_type"].GetString() == "incompressible":
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(
                self.main_model_part, 
                time_scheme, 
                self.linear_solver,
                builder_and_solver,
                self.settings["compute_reactions"].GetBool(), 
                self.settings["reform_dofs_at_each_step"].GetBool(), 
                self.settings["calculate_solution_norm"].GetBool(), 
                move_mesh_flag)
        else:
            conv_criteria = KratosMultiphysics.ResidualCriteria(
                self.settings["relative_tolerance"].GetDouble(), 
                self.settings["absolute_tolerance"].GetDouble())
            max_iterations = self.settings["maximum_iterations"].GetInt()
                    
            self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
                self.main_model_part, 
                time_scheme, 
                self.linear_solver,
                conv_criteria,
                max_iterations,
                self.settings["compute_reactions"].GetBool(), 
                self.settings["reform_dofs_at_each_step"].GetBool(), 
                move_mesh_flag)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Check()
    def PrepareModelPart(self):
        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            self.model.AddModelPart(self.main_model_part)
   

    def ImportModelPart(self):
        """This function imports the ModelPart
        """
        if self.solver_imports_model_part:
            # print(self.model.GetModelPart(self.settings["model_part_name"].GetString()))
            self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])
            if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
                #here it would be the place to import restart data if required
                print(self.settings["model_import_settings"]["input_filename"].GetString())
                # IOdir=self.settings["model_import_settings"]["input_filename"].GetString()
            
                # KratosMultiphysics.ModelPartIO(IOdir).ReadModelPart(self.main_model_part)
                throw_errors = False
                # KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part,throw_errors).Execute()
                #here we replace the dummy elements we read with proper elements
                element_replace_settings=KratosMultiphysics.Parameters()
                element_replace_settings.AddEmptyValue("element_replace_settings")
                if (self.settings["problem_type"].GetString() == "incompressible"):
                    if(self.domain_size == 3):
                        element_replace_settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                            {
                            "element_name":"IncompressiblePotentialFlowElement3D4N",
                            "condition_name": "IncompressiblePotentialWallCondition3D3N"
                            }
                            """)
                    elif(self.domain_size == 2):
                        element_replace_settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                            {
                            "element_name":"IncompressiblePotentialFlowElement2D3N",
                            "condition_name": "IncompressiblePotentialWallCondition2D2N"
                            }
                            """)
                    else:
                        raise Exception("Domain size is not 2 or 3!!")
                elif (self.settings["problem_type"].GetString() == "compressible"):
                    if(self.domain_size == 3):
                        element_replace_settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                            {
                            "element_name":"CompressiblePotentialFlowElement3D4N",
                            "condition_name": "CompressiblePotentialWallCondition3D3N"
                            }
                            """)
                    elif(self.domain_size == 2):
                        element_replace_settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                            {
                            "element_name":"CompressiblePotentialFlowElement2D3N",
                            "condition_name": "CompressiblePotentialWallCondition2D2N"
                            }
                            """)
                    else:
                        raise Exception("Domain size is not 2 or 3!!")
                elif (self.settings["problem_type"].GetString() == "incompressible_stresses"):
                    if(self.domain_size == 2):
                        element_replace_settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                            {
                            "element_name":"IncompressibleStressesPotentialFlowElement2D3N",
                            "condition_name": "IncompressiblePotentialWallCondition2D2N"
                            }
                            """)
                    else:
                        raise Exception("Domain size is not 2!!")
                else:
                    raise Exception("Problem type not defined!!")
                
                KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, element_replace_settings["element_replace_settings"]).Execute()
                
            else:
                raise Exception("other input options are not yet implemented")
            print("Solving",self.settings["problem_type"].GetString() ,"case")
            current_buffer_size = self.main_model_part.GetBufferSize()
            if(self.GetMinimumBufferSize() > current_buffer_size):
                self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                    
            print ("model reading finished")
        
    def GetMinimumBufferSize(self):
        return 2;

    def GetComputingModelPart(self):
        return self.main_model_part

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here
    
    def AdvanceInTime(self, current_time):
        dt = 1 #self._ComputeDeltaTime()
        new_time = current_time + dt

        # self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):        
        self.solver.InitializeSolutionStep()


    def SolveSolutionStep(self):
        (self.solver).Solve() 

    def FinalizeSolutionStep(self):        
        self.solver.FinalizeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

