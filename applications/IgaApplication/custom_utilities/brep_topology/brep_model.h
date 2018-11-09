#if !defined(KRATOS_BREP_MODEL_H_INCLUDED )
#define  KRATOS_BREP_MODEL_H_INCLUDED

// System includes
#include <vector>

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"

#include "brep_face.h"
#include "brep_edge.h"
#include "brep_vertex.h"

//#include "..\nurbs_brep_modeler.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    /// Holder of all Brep entities.
    /** This class contains the full Brep description needed to describe all CAD-geometries.
    */
    class BrepModel : public IndexedObject, public Flags
    {
    public:
        ///@name Type Definitions
        ///@{
        //KRATOS_CLASS_POINTER_DEFINITION(BrepModel);

        ///@}
        ///@name Life Cycle 
        ///@{

        bool GetIntegrationDomain(
            ModelPart& rModelPart, 
            int& brep_id, 
            const std::string& rType,
            const std::string& rName,
            const int& rPropertiesId,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables);

        std::vector<BrepFace>&   GetFaceVector();
        std::vector<BrepEdge>&   GetEdgeVector();
        std::vector<BrepVertex>& GetVertexVector();

        /// Constructor
        BrepModel::BrepModel(
            int& brep_id,
            double& model_tolerance,
            std::vector<BrepFace>& faces,
            std::vector<BrepEdge>& edges,
            std::vector<BrepVertex>& vertices)
            : m_model_tolerance(model_tolerance),
              m_brep_faces(faces),
              m_brep_edges(edges),
              m_brep_vertices(vertices),
              IndexedObject(brep_id),
              Flags()
        {

        };

        /// Destructor.
        virtual ~BrepModel()
        {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{

        double                  m_model_tolerance;
        std::vector<BrepFace>   m_brep_faces;
        std::vector<BrepEdge>   m_brep_edges;
        std::vector<BrepVertex> m_brep_vertices;

        ///@}
    }; // Class BrepModel

}  // namespace Kratos.

#endif // KRATOS_BREP_MODEL_H_INCLUDED  defined