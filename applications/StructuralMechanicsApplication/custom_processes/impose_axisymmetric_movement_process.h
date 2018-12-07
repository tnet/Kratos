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

#if !defined(KRATOS_IMPOSE_AXISYMMETRIC_MOVEMENT_PROCESS)
#define KRATOS_IMPOSE_AXISYMMETRIC_MOVEMENT_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
/** 
 * @class ImposeAxisymmetricMovementProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This method assign linear kinematic constrains to a certain submodelpart
 * @details It assigns to the first node of the submodelpart by default or to a certain node if provided
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ImposeAxisymmetricMovementProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImposeAxisymmetricMovementProcess
    KRATOS_CLASS_POINTER_DEFINITION(ImposeAxisymmetricMovementProcess);
    
    /// General type definitions
    typedef Node<3>                                                      NodeType;

    /// General containers type definitions
    typedef ModelPart::MasterSlaveConstraintContainerType ConstraintContainerType;

    /// Component definition
    typedef VectorComponentAdaptor< array_1d< double, 3 > >         ComponentType;

    /// The DoF type definition
    typedef Dof<double>                                                   DofType;

    /// The DoF pointer vector type definition
    typedef std::vector< DofType::Pointer >                  DofPointerVectorType;

    /// Definitions of the integers
    typedef std::size_t                                                 IndexType;
    typedef std::size_t                                                  SizeType;

    /// Definition of the 3x3 bounded rotation matrix
    typedef BoundedMatrix<double, 3, 3>                        RotationMatrixType;

    /// Definition of the local relation map type
    typedef std::unordered_map<IndexType, double>            LocalRelationMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rThisModelPart The model part to compute
     * @param ThisParameters The parameters of configuration
     */
    ImposeAxisymmetricMovementProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~ImposeAxisymmetricMovementProcess() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;
    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ImposeAxisymmetricMovementProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ImposeAxisymmetricMovementProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{
    
    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    
    ModelPart& mrThisModelPart; /// The model part to compute
    Parameters mThisParameters; /// The parameters (can be used for general pourposes)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This creates a new node in the center of gravity as center reference
     * @param rNodeId The reference o the new node created
     * @param rAxisymmetricModelPart The axisymmetric model part
     */
    void CreateNewNode(
        IndexType& rNodeId,
        ModelPart& rAxisymmetricModelPart
        );

    /**
     * @brief This fills the auxiliar model part with the auxiliar elements and nodes
     * @param rReferenceModelPart The reference model part
     * @param rAxisymmetricModelPart The axisymmetric model part
     */
    std::unordered_map<IndexType, LocalRelationMapType> FillReferenceModelPart(
        ModelPart& rReferenceModelPart,
        ModelPart& rAxisymmetricModelPart
        );

    /**
     * @brief This clears the auxiliar model part
     * @param rReferenceModelPart The reference model part
     */
    void ClearReferenceModelPart(ModelPart& rReferenceModelPart);

    /**
     * @brief This method computes the rotation matrix arround an axis
     * @param rAxisVector The vector that define the rotation axis
     * @param Theta The rotation angle
     */
    RotationMatrixType ComputeRotationMatrix(
        const Vector& rAxisVector,
        const double Theta
        );

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ImposeAxisymmetricMovementProcess& operator=(ImposeAxisymmetricMovementProcess const& rOther) = delete;

    /// Copy constructor.
    //ImposeAxisymmetricMovementProcess(ImposeAxisymmetricMovementProcess const& rOther);


    ///@}

}; // Class ImposeAxisymmetricMovementProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ImposeAxisymmetricMovementProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ImposeAxisymmetricMovementProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}
#endif /* KRATOS_IMPOSE_AXISYMMETRIC_MOVEMENT_PROCESS defined */
