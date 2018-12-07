// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "structural_mechanics_application_variables.h"

/* Utilities */
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"
#include "utilities/intersection_utilities.h"

/* Processes */
#include "custom_processes/impose_axisymmetric_movement_process.h"
#include "custom_processes/compute_center_of_gravity_process.h"

namespace Kratos
{
ImposeAxisymmetricMovementProcess::ImposeAxisymmetricMovementProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"             : "please_specify_model_part_name",
        "new_model_part_name"         : "Axisymmetric_Movement_ModelPart",
        "axisymmetry_axis"            : [0.0,0.0,1.0],
        "max_number_of_searchs"       : 1000,
        "master_variable_name"        : "DISPLACEMENT",
        "slave_variable_name"         : "",
        "master_node_id"              : 0
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeAxisymmetricMovementProcess::Execute()
{
    KRATOS_TRY

    // We execute the different steps
    ExecuteInitialize();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeAxisymmetricMovementProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Getting model parts
    ModelPart& r_root_model_part = mrThisModelPart.GetRootModelPart();
    ModelPart& r_model_part = r_root_model_part.GetSubModelPart(mThisParameters["model_part_name"].GetString());
    const std::string& new_model_part_name = mThisParameters["new_model_part_name"].GetString();
    ModelPart& r_axisymmetric_model_part = new_model_part_name != r_model_part.Name() ? r_model_part.HasSubModelPart(new_model_part_name) ? r_model_part.GetSubModelPart(new_model_part_name) : r_model_part.CreateSubModelPart(new_model_part_name) : r_model_part;

    // Dimension check
    KRATOS_ERROR_IF_NOT(r_root_model_part.GetProcessInfo()[DOMAIN_SIZE] == 3) << "In order to use this process yo need to be working on 3 dimensions" << std::endl;

    // Reorder constrains
    IndexType constraint_id = 1;
    for (auto& r_constrain : r_root_model_part.MasterSlaveConstraints()) {
        r_constrain.SetId(constraint_id);
        ++constraint_id;
    }

    // Getting list of variables
    std::vector<Variable<double>> master_double_list_variables, slave_double_list_variables;
    std::vector<VariableComponent<ComponentType>> master_components_list_variables, slave_components_list_variables;
    std::vector<Variable<array_1d< double, 3>>> master_vector_list_variables, slave_vector_list_variables;
    const std::string& r_master_variable_name = mThisParameters["master_variable_name"].GetString();
    // The master variable
    if(KratosComponents<Variable<double>>::Has(r_master_variable_name)){
        Variable<double> variable = KratosComponents<Variable<double>>::Get(r_master_variable_name);
        master_double_list_variables.push_back(variable);
    } else if (KratosComponents< VariableComponent<ComponentType>>::Has(r_master_variable_name)) {
        VariableComponent<ComponentType> variable = KratosComponents<VariableComponent<ComponentType>>::Get(r_master_variable_name);
        master_components_list_variables.push_back(variable);
    } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(r_master_variable_name)) {
        Variable<array_1d< double, 3>> variable = KratosComponents<Variable<array_1d< double, 3>>>::Get(r_master_variable_name);
        master_vector_list_variables.push_back(variable);
    } else {
        KRATOS_ERROR << "Only double, components and vector variables are allowed in the variables list." ;
    }
    const std::string& r_slave_variable_name = mThisParameters["slave_variable_name"].GetString();
    // We get the slave variable list
    if (r_slave_variable_name != "") {
        if(KratosComponents<Variable<double>>::Has(r_slave_variable_name)){
            Variable<double> variable = KratosComponents<Variable<double>>::Get(r_slave_variable_name);
            slave_double_list_variables.push_back(variable);
        } else if (KratosComponents< VariableComponent<ComponentType>>::Has(r_slave_variable_name)) {
            VariableComponent<ComponentType> variable = KratosComponents<VariableComponent<ComponentType>>::Get(r_slave_variable_name);
            slave_components_list_variables.push_back(variable);
        } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(r_slave_variable_name)) {
            Variable<array_1d< double, 3>> variable = KratosComponents<Variable<array_1d< double, 3>>>::Get(r_slave_variable_name);
            slave_vector_list_variables.push_back(variable);
        } else {
            KRATOS_ERROR << "Only double, components and vector variables are allowed in the variables list." ;
        }
    } else { // Else we consider exactly the same list of variables
        for (auto& var : master_double_list_variables)
            slave_double_list_variables.push_back(var);
        for (auto& var : master_components_list_variables)
            slave_components_list_variables.push_back(var);
    }

    // Getting index of the master node
    int master_node_id = mThisParameters["master_node_id"].GetInt();

    // We iterate over the nodes of the rigid model part
    auto nodes_array = r_axisymmetric_model_part.Nodes();
    const int number_of_nodes = static_cast<int>(nodes_array.size());

    // In case no existing node is assigned we create a new node in the center of gravity of the current model part
    IndexType new_node_id = 1; // Can be used later in case we use a new node created ad hoc
    if (master_node_id == 0) {
        CreateNewNode(new_node_id, r_axisymmetric_model_part);
        mThisParameters["master_node_id"].SetInt(new_node_id);
        master_node_id = mThisParameters["master_node_id"].GetInt();
    }

    // We create an auxiliar model part to generate the axisymmetry (this is a cut in the transversal direction)
    Model& r_model = r_root_model_part.GetModel();
    ModelPart& r_reference_model_part = r_model.CreateModelPart("REFERENCE_AXISYMMETRY_MODEL_PART");
    auto all_local_relations = FillReferenceModelPart(r_reference_model_part, r_axisymmetric_model_part);

    // We create the locator
    auto point_locator = BinBasedFastPointLocator<2>(r_reference_model_part);
    point_locator.UpdateSearchDatabase();

    // List of variables
    const SizeType number_of_double_variables = master_double_list_variables.size();
    const SizeType number_of_components_variables = master_components_list_variables.size();
    const SizeType number_of_vector_variables = master_vector_list_variables.size();

    // Reference constraint
    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    // Get the axis
    const Vector& r_axis_vector = mThisParameters["axisymmetry_axis"].GetVector();

    // Getting tangent plane
    array_1d<double, 3> tangent_xi, tangent_eta;
    MathUtils<double>::OrthonormalBasis(r_axis_vector, tangent_xi, tangent_eta);

    // If we master node ID is zero then we get the center of gravity of the model part
    NodeType::Pointer p_reference_node = r_root_model_part.pGetNode(master_node_id);
    const array_1d<double, 3>& r_reference_node_coordinates = p_reference_node->Coordinates();

    // Creation of the constraints
    #pragma omp parallel
    {
        // Buffer of constraints
        ConstraintContainerType constraints_buffer;
        constraints_buffer.reserve(number_of_nodes);

        // Auxiliar values
        Vector shape_functions;
        Element::Pointer p_element;

        #pragma omp for
        for (int i = 0; i < number_of_nodes; ++i) {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Id() != p_reference_node->Id()) {
                // Compute the radius and position respect the axis
                const array_1d<double, 3>& r_current_node_coordinates = it_node->Coordinates();

                const array_1d<double, 3> vector_points = r_current_node_coordinates - r_reference_node_coordinates;
                const double distance = inner_prod(vector_points, r_axis_vector);
                const array_1d<double, 3> clossest_point = r_reference_node_coordinates + r_axis_vector * distance;
                array_1d<double, 3> axisymmetric_vector = clossest_point - r_current_node_coordinates;
                const double norm = norm_2(axisymmetric_vector);
                // Normalize
                if (norm > std::numeric_limits<double>::epsilon())
                    axisymmetric_vector /= norm;

                const double theta = std::acos(inner_prod(tangent_xi, axisymmetric_vector)/(norm_2(tangent_xi) * norm_2(axisymmetric_vector)));
                const RotationMatrixType rotation_matrix = ComputeRotationMatrix(r_axis_vector, theta);

                const array_1d<double, 3>& coordinates = prod(rotation_matrix, r_current_node_coordinates);
                const bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, mThisParameters["max_number_of_searchs"].GetInt(), 5.0e-2);

                // It is found
                if (is_found) {
                    // The geoemtry of the reference geometry
                    auto& r_geometry = p_element->GetGeometry();
                    const IndexType id_1 = r_geometry[0].Id();
                    const auto& local_relation_1 = all_local_relations[id_1];
                    const IndexType id_2 = r_geometry[1].Id();
                    const auto& local_relation_2 = all_local_relations[id_2];
                    const IndexType id_3 = r_geometry[2].Id();
                    const auto& local_relation_3 = all_local_relations[id_3];

                    std::vector<IndexType> indexes;
                    std::vector<double> coefficients;
                    if (shape_functions[0] > std::numeric_limits<double>::epsilon()) {
                        for (auto& relation : local_relation_1) {
                            indexes.push_back(relation.first);
                            coefficients.push_back(relation.second * shape_functions[0]);
                        }
                    }
                    if (shape_functions[1] > std::numeric_limits<double>::epsilon()) {
                        for (auto& relation : local_relation_2) {
                            indexes.push_back(relation.first);
                            coefficients.push_back(relation.second * shape_functions[1]);
                        }
                    }
                    if (shape_functions[2] > std::numeric_limits<double>::epsilon()) {
                        for (auto& relation : local_relation_3) {
                            indexes.push_back(relation.first);
                            coefficients.push_back(relation.second * shape_functions[2]);
                        }
                    }
                    const SizeType master_size = indexes.size();

                    // We create the constraints
                    if ((number_of_double_variables + number_of_components_variables) > 0) {
                        Matrix relation_matrix(1, master_size);
                        for (IndexType index = 0; index < coefficients.size(); ++index) {
                            relation_matrix(0, index) = coefficients[index];
                        }
                        Vector constant_vector = ZeroVector(1);

                        DofPointerVectorType master_dofs(master_size);
                        DofPointerVectorType slave_dofs(1);
                        for (IndexType i_var = 0; i_var < number_of_double_variables; ++i_var) {
                            for (IndexType index = 0; index < indexes.size(); ++index) {
                                master_dofs[index] = r_root_model_part.pGetNode(indexes[index])->pGetDof(master_double_list_variables[i_var]);
                            }
                            slave_dofs[0] = it_node->pGetDof(slave_double_list_variables[i_var]);
                            auto constraint = r_clone_constraint.Create(constraint_id + (i * number_of_double_variables + i_var) + 1, master_dofs, slave_dofs, relation_matrix, constant_vector);
                            (constraints_buffer).insert((constraints_buffer).begin(), constraint);
                        }
                        for (IndexType i_var = 0; i_var < number_of_components_variables; ++i_var) {
                            for (IndexType index = 0; index < indexes.size(); ++index) {
                                master_dofs[index] = r_root_model_part.pGetNode(indexes[index])->pGetDof(master_components_list_variables[i_var]);
                            }
                            slave_dofs[0] = it_node->pGetDof(slave_components_list_variables[i_var]);
                            auto constraint = r_clone_constraint.Create(constraint_id + (i * number_of_double_variables + i_var) + 1, master_dofs, slave_dofs, relation_matrix, constant_vector);
                            (constraints_buffer).insert((constraints_buffer).begin(), constraint);
                        }
                    }
                    for (IndexType i_var = 0; i_var < number_of_vector_variables; ++i_var) {
//                         auto constraint = r_clone_constraint.Create(constraint_id + (i * number_of_vector_variables + i_var) + 1, *p_reference_node, master_components_list_variables[i_var], *it_node, slave_components_list_variables[i_var], relation, constant);
//                         (constraints_buffer).insert((constraints_buffer).begin(), constraint);
                    }
                }
            }
        }

        // We transfer
        #pragma omp critical
        {
            r_axisymmetric_model_part.AddMasterSlaveConstraints(constraints_buffer.begin(),constraints_buffer.end());
            mrThisModelPart.AddMasterSlaveConstraints(constraints_buffer.begin(),constraints_buffer.end());
        }

        // Remove reference model part
        ClearReferenceModelPart(r_axisymmetric_model_part);
        r_model.DeleteModelPart("REFERENCE_AXISYMMETRY_MODEL_PART");
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeAxisymmetricMovementProcess::CreateNewNode(
    IndexType& rNodeId,
    ModelPart& rAxisymmetricModelPart
    )
{
    // First we compute the center of gravity
    ComputeCenterOfGravityProcess(mrThisModelPart).Execute();

    // We recover the coordinates of the center of gravity
    const array_1d<double, 3>& center_of_gravity_coordinates = mrThisModelPart.GetProcessInfo()[CENTER_OF_GRAVITY];

    // Now we create the new node
    for (auto& r_node : mrThisModelPart.GetRootModelPart().Nodes()) {
        r_node.SetId(rNodeId);
        ++rNodeId;
    }
    rAxisymmetricModelPart.CreateNewNode(rNodeId, center_of_gravity_coordinates[0], center_of_gravity_coordinates[1], center_of_gravity_coordinates[2]);
}

/***********************************************************************************/
/***********************************************************************************/

std::unordered_map<ImposeAxisymmetricMovementProcess::IndexType, ImposeAxisymmetricMovementProcess::LocalRelationMapType> ImposeAxisymmetricMovementProcess::FillReferenceModelPart(
    ModelPart& rReferenceModelPart,
    ModelPart& rAxisymmetricModelPart
    )
{
    // We iterate over the nodes of the rigid model part
    auto nodes_array = rAxisymmetricModelPart.Nodes();
    const int number_of_nodes = static_cast<int>(nodes_array.size());

    // The auxiliar properties
    auto p_prop = rReferenceModelPart.pGetProperties(0);

    // Compute maximum radius and axis coordinate
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<double> max_radius_vector(num_threads, 0.0);
    std::vector<double> max_axis_vector(num_threads, 0.0);
    double radius, axis;

    // If we master node ID is zero then we get the center of gravity of the model part
    ModelPart& r_root_model_part = mrThisModelPart.GetRootModelPart();
    const int master_node_id = mThisParameters["master_node_id"].GetInt();
    NodeType::Pointer p_reference_node = r_root_model_part.pGetNode(master_node_id);
    const array_1d<double, 3>& r_reference_node_coordinates = p_reference_node->Coordinates();

    // Get the axis
    const Vector& r_axis_vector = mThisParameters["axisymmetry_axis"].GetVector();

    #pragma omp parallel for private(radius)
    for (int i = 0; i < number_of_nodes; ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Id() != p_reference_node->Id()) {
            // Compute the radius and position respect the axis
            const array_1d<double, 3>& r_current_node_coordinates = it_node->Coordinates();

            const array_1d<double, 3> vector_points = r_current_node_coordinates - r_reference_node_coordinates;
            radius = inner_prod(vector_points, r_axis_vector);
            const array_1d<double, 3> clossest_point = r_reference_node_coordinates + r_axis_vector * radius;
            array_1d<double, 3> axisymmetric_vector = clossest_point - r_current_node_coordinates;
            axis = std::sqrt(inner_prod(vector_points, vector_points) - std::pow(radius, 2));

            const int id = OpenMPUtils::ThisThread();

            if (radius > max_radius_vector[id])
                max_radius_vector[id] = radius;

            if (axis > max_axis_vector[id])
                max_axis_vector[id] = axis;
        }
    }

    const double max_radius = *std::max_element(max_radius_vector.begin(), max_radius_vector.end());
    const double max_axis = *std::max_element(max_axis_vector.begin(), max_axis_vector.end());

    // Getting tangent plane
    array_1d<double, 3> tangent_xi, tangent_eta;
    MathUtils<double>::OrthonormalBasis(r_axis_vector, tangent_xi, tangent_eta);

    // Creating cutting plane (a two very big triangles)
    PointerVector<Point> points_array_1(3);
    PointerVector<Point> points_array_2(3);
    points_array_1(0) = Kratos::make_shared<Point>(r_reference_node_coordinates + max_axis * r_axis_vector);
    points_array_1(1) = Kratos::make_shared<Point>(r_reference_node_coordinates - max_axis * r_axis_vector);
    points_array_1(2) = Kratos::make_shared<Point>(max_radius * tangent_xi + r_reference_node_coordinates + max_axis * r_axis_vector);
    points_array_2(0) = points_array_1(0);
    points_array_2(1) = points_array_1(1);
    points_array_2(2) = Kratos::make_shared<Point>(max_radius * tangent_xi+ r_reference_node_coordinates - max_axis * r_axis_vector);

    // We define the auxiliar geometry
    auto triangle_1 = Triangle2D3<Point>(points_array_1);
    auto triangle_2 = Triangle2D3<Point>(points_array_2);

    // Auxiliar declarations
    std::unordered_map<IndexType, LocalRelationMapType> all_local_relations;

    array_1d<double, 3> intersection_point;
    int intersection;
    IndexType index_node = 0;
    IndexType index_element = 0;
    // Loop over the elements (something can be done to not cut all the elements, because it is expensive)
    for (auto& r_elem : rAxisymmetricModelPart.Elements()) {

        // The list of created nodes
        std::vector<NodeType::Pointer> created_nodes_list;

        // First we iterate over the edges of the elements
        auto& r_geometry = r_elem.GetGeometry();
        const auto& r_edges = r_geometry.Edges();
        for (auto& r_edge : r_edges) {
            // First intersection try
            intersection = IntersectionUtilities::ComputeTriangleLineIntersection(triangle_1, r_edge[0].Coordinates(), r_edge[0].Coordinates(), intersection_point);

            // Second intersection try
            if (intersection == 0) {
                intersection = IntersectionUtilities::ComputeTriangleLineIntersection(triangle_2, r_edge[0].Coordinates(), r_edge[0].Coordinates(), intersection_point);
            }

            // Create new node
            if (intersection == 1) {
                // Actually creating the new node
                ++index_node;
                auto p_node = rReferenceModelPart.CreateNewNode(index_node, intersection_point[0], intersection_point[1], intersection_point[2]);
                created_nodes_list.push_back(p_node);

                // We get the local coordinates and the shape functions
                Point global_point(intersection_point);
                Point local_point;
                r_edge.PointLocalCoordinates( local_point, global_point);
                Vector N(2);
                N = r_edge.ShapeFunctionsValues(N, local_point);
                LocalRelationMapType this_relation;
                if (N[0] > std::numeric_limits<double>::epsilon()) this_relation.insert({r_edge[0].Id(), N[0]});
                if (N[1] > std::numeric_limits<double>::epsilon()) this_relation.insert({r_edge[1].Id(), N[1]});
                all_local_relations.insert({index_node, this_relation});
            }
        }

        // Now we create a set of triangles with the created nodes
        if (created_nodes_list.size() > 2) { // At least one triangle
            std::vector<IndexType> index_nodes_list(3);
            for (IndexType i_node = 0; i_node < created_nodes_list.size() - 2; ++i_node) {
                ++index_element;
                index_nodes_list[0] = created_nodes_list[i_node]->Id();
                index_nodes_list[1] = created_nodes_list[i_node + 1]->Id();
                index_nodes_list[2] = created_nodes_list[i_node + 2]->Id();
                rReferenceModelPart.CreateNewElement("Element2D3N", index_element, index_nodes_list, p_prop);
            }
        }
    }

    return all_local_relations;
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeAxisymmetricMovementProcess::ClearReferenceModelPart(ModelPart& rReferenceModelPart)
{
    VariableUtils().SetFlag(TO_ERASE, true, rReferenceModelPart.Conditions());
    VariableUtils().SetFlag(TO_ERASE, true, rReferenceModelPart.Elements());
    VariableUtils().SetFlag(TO_ERASE, true, rReferenceModelPart.Nodes());
    rReferenceModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
    rReferenceModelPart.RemoveElementsFromAllLevels(TO_ERASE);
    rReferenceModelPart.RemoveNodesFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/

ImposeAxisymmetricMovementProcess::RotationMatrixType ImposeAxisymmetricMovementProcess::ComputeRotationMatrix(
    const Vector& rAxisVector,
    const double Theta
    )
{
    RotationMatrixType rotation_matrix;

    rotation_matrix(0, 0) = std::cos(Theta) + std::pow(rAxisVector[0], 2) * (1.0 - std::cos(Theta));
    rotation_matrix(0, 1) = rAxisVector[0] * rAxisVector[1] * (1.0 - std::cos(Theta)) - rAxisVector[2] * std::sin(Theta);
    rotation_matrix(0, 2) = rAxisVector[0] * rAxisVector[2] * (1.0 - std::cos(Theta)) - rAxisVector[1] * std::sin(Theta);

    rotation_matrix(1, 0) = rAxisVector[0] * rAxisVector[1] * (1.0 - std::cos(Theta)) - rAxisVector[2] * std::sin(Theta);
    rotation_matrix(1, 1) = std::cos(Theta) + std::pow(rAxisVector[1], 2) * (1.0 - std::cos(Theta));
    rotation_matrix(1, 2) = rAxisVector[1] * rAxisVector[2] * (1.0 - std::cos(Theta)) - rAxisVector[0] * std::sin(Theta);

    rotation_matrix(2, 0) = rAxisVector[2] * rAxisVector[0] * (1.0 - std::cos(Theta)) - rAxisVector[1] * std::sin(Theta);
    rotation_matrix(2, 1) = rAxisVector[1] * rAxisVector[2] * (1.0 - std::cos(Theta)) - rAxisVector[0] * std::sin(Theta);
    rotation_matrix(2, 2) = std::cos(Theta) + std::pow(rAxisVector[2], 2) * (1.0 - std::cos(Theta));

    return rotation_matrix;
}

// class ImposeAxisymmetricMovementProcess
} // namespace Kratos
