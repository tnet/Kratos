//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#ifndef KRATOS_COMPUTE_GRADIENT_ADJOINT_PROCESS_H
#define KRATOS_COMPUTE_GRADIENT_ADJOINT_PROCESS_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

#include <string>
#include <iostream>
#include <sstream>

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <unordered_map> //TODO: remove this dependence when Kratos has en internal one
#include <utility>

namespace Kratos
{

class ComputeGradientAdjointProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(ComputeGradientAdjointProcess);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for ComputeGradientAdjointProcess Process
//     ComputeGradientAdjointProcess(ModelPart& rModelPart,
//                      KratosParameters& parameters
//                     ):
//         Process(),
//         mrModelPart(rModelPart),
//         mrOptions(Flags()),
//         mrParameters(parameters)
//     {
//     }
    /// Constructor for ComputeGradientAdjointProcess Process
    ComputeGradientAdjointProcess(ModelPart& rModelPart,
                Matrix& rdRdu,Matrix& rdRdx, Vector& rdFdu                       
                    ):
        Process(),
        mrModelPart(rModelPart),
        mrdRdu(rdRdu),
        mrdRdx(rdRdx),
        mrdFdu(rdFdu)
    {
    }

    /// Destructor.
    ~ComputeGradientAdjointProcess() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;
        std::cout<<"Compute Gradient Adjoint Process Started"<< std::endl;
        
        const double epsilon = 1e-9;

        for(auto it=mrModelPart.ElementsBegin(); it!=mrModelPart.ElementsEnd(); ++it)
        {
            auto geom = it->GetGeometry();
            
            const unsigned int NumNodes=geom.size();
            Matrix LHS_ref = ZeroMatrix(3,3);
            Vector RHS_ref = ZeroVector(3);
            std::vector<std::size_t> id_vector;
            it -> CalculateLocalSystem(LHS_ref,RHS_ref,mrModelPart.GetProcessInfo());
            it -> EquationIdVector(id_vector,mrModelPart.GetProcessInfo());
            const unsigned int NumDof=RHS_ref.size();
            KRATOS_WATCH(NumDof);
            if (it->Is(BOUNDARY)){ 
                Matrix dRdx_elemental = ZeroMatrix(NumDof,NumDof);
                Vector dFdu_elemental = ZeroVector(NumDof);
                Matrix LHS = ZeroMatrix(NumDof,NumDof);
                Vector RHS = ZeroVector(NumDof);
                double vol;
                const array_1d<double, 3> vinfinity = mrModelPart.GetProcessInfo()[VELOCITY_INFINITY];
                const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);
                bounded_matrix<double, 3, 2 > DN_DX;
                array_1d<double, 3 > N;
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, vol);
                Vector normal = ZeroVector(3);
                normal= it -> GetValue(NORMAL);
                for(unsigned int i = 0; i<NumNodes; i++){                    
                    geom[i].GetSolutionStepValue(LEVEL_SET_DISTANCE) += epsilon;
                    it -> CalculateLocalSystem(LHS,RHS,mrModelPart.GetProcessInfo());
                    geom[i].GetSolutionStepValue(LEVEL_SET_DISTANCE) += -epsilon;
                    for (unsigned int j= 0; j<NumDof;j++){
                        dRdx_elemental(i,j) = (RHS_ref(j)-RHS(j))/epsilon;
                        if(NumNodes!=NumDof)
                            dRdx_elemental(i+NumNodes,j) = (RHS_ref(j)-RHS(j))/epsilon;
                    }          
                }

                if(NumDof==NumNodes){
                    array_1d<double,3> phis;  

                    for(unsigned int i = 0; i<NumNodes; i++){   
                        phis(i) = geom[i].GetSolutionStepValue(POSITIVE_POTENTIAL);
                    }
                    dFdu_elemental=-normal(1)*2.0/vinfinity_norm2*(prod(Matrix(prod(DN_DX, trans(DN_DX))),phis)+trans(prod(phis,Matrix(prod(DN_DX, trans(DN_DX))))));
                }else{
                    array_1d<double,3> phis_pos;
                    array_1d<double,3> phis_neg; 
                    Vector dFdu_elemental_pos = ZeroVector(NumNodes); 
                    Vector dFdu_elemental_neg = ZeroVector(NumNodes); 

  
                    for(unsigned int i = 0; i<NumNodes; i++){
                        if (geom[i].GetSolutionStepValue(LEVEL_SET_DISTANCE)>0){
                            phis_pos(i) = geom[i].GetSolutionStepValue(POSITIVE_POTENTIAL);
                            phis_neg(i) = geom[i].GetSolutionStepValue(NEGATIVE_POTENTIAL);
                        }
                        else{
                            phis_pos(i) = geom[i].GetSolutionStepValue(NEGATIVE_POTENTIAL);
                            phis_neg(i) = geom[i].GetSolutionStepValue(POSITIVE_POTENTIAL);
                        }
                    }
                    dFdu_elemental_pos=-normal(1)*2.0/vinfinity_norm2*(prod(Matrix(prod(DN_DX, trans(DN_DX))),phis_pos)+trans(prod(phis_pos,Matrix(prod(DN_DX, trans(DN_DX))))));
                    dFdu_elemental_neg=-normal(1)*2.0/vinfinity_norm2*(prod(Matrix(prod(DN_DX, trans(DN_DX))),phis_neg)+trans(prod(phis_neg,Matrix(prod(DN_DX, trans(DN_DX))))));
                    for(unsigned int i = 0; i<NumNodes; i++){
                        dFdu_elemental(i)=dFdu_elemental_pos(i);
                        dFdu_elemental(i+NumNodes)=dFdu_elemental_neg(i);
                    }
                }

                unsigned int I;
                unsigned int J;
                for (unsigned int i_node = 0; i_node<NumDof;++i_node){
                    for (unsigned int j_node = 0; j_node<NumDof;++j_node){ 
                        I = id_vector[i_node];//geom[i_node].Id();
                        J = id_vector[j_node];//geom[j_node].Id();
                        
                        mrdRdx(I,J) += dRdx_elemental(i_node,j_node);
                        mrdRdu(J,I) += LHS_ref(i_node,j_node); //transposed dRdu
                    }      
                    mrdFdu(I) += dFdu_elemental(i_node);
                }         
            }
        }
        // double number_of_nodes=mrModelPart.NumberOfNodes();
        // Vector lambda_0=ZeroVector(number_of_nodes);
        // Vector lambda=ZeroVector(number_of_nodes);
        // Vector r_0=ZeroVector(number_of_nodes);
        // Vector r=ZeroVector(number_of_nodes);
        // Vector p=ZeroVector(number_of_nodes);
        // r_0=mrdFdu-prod(mrdRdu,lambda_0) ;
        // p=r_0;

        // double norm_diff = 1;
        // double alpha;
        // double beta;
        // int iteration = 0;

        // while (norm_diff > 1e-9){
        //     iteration += 1;
        //     std::cout<<"Iteration: "<<iteration <<std::endl;
        //     alpha=inner_prod(r_0,r_0)/inner_prod(Vector(prod(mrdRdu,p)),p);
        //     lambda=lambda_0+alpha*p;
        //     norm_diff=std::abs(inner_prod(lambda,lambda)-inner_prod(lambda_0,lambda_0));
        //     std::cout<<"norm_diff:" << norm_diff<< std::endl;
        //     r=r_0-alpha*prod(mrdRdu,p);
        //     beta=inner_prod(r,r)/inner_prod(r_0,r_0);
        //     p=r+beta*p;

        //     r_0=r; 
        //     lambda_0=lambda;
        // }

        // for(auto it=mrModelPart.NodesBegin(); it!=mrModelPart.NodesEnd(); ++it)
        // {
        //     if(it->Is(INLET))
        // }
        // for (unsigned int i = 0; i < number_of_nodes;++i){
        //     for (unsigned int j = 0; j < number_of_nodes;++j){
        //         if (!(mrdRdx(i,j)==0.0))
        //             KRATOS_WATCH(mrdRdx(i,j));
        //     }
        // }
        // for (unsigned int i = 0; i < number_of_nodes;++i){
        //     for (unsigned int j = 0; j < number_of_nodes;++j){
        //         if (!(mrdRdu(i,j)==0.0))
        //             KRATOS_WATCH(mrdRdu(i,j));
        //     }
        // }
        // for (unsigned int i = 0; i < number_of_nodes;++i){
        //         if (!(mrdFdu(i)==0.0))
        //             KRATOS_WATCH(mrdFdu(i));
        // }
       
        KRATOS_CATCH("");
    }

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
        return "ComputeGradientAdjointProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeGradientAdjointProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Matrix& mrdRdu;
    Matrix& mrdRdx;
    Vector& mrdFdu;
    Flags mrOptions;


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeGradientAdjointProcess& operator=(ComputeGradientAdjointProcess const& rOther);

    /// Copy constructor.
    ComputeGradientAdjointProcess(ComputeGradientAdjointProcess const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeGradientAdjointProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeGradientAdjointProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_ComputeLift_LEVEL_SET_PROCESS_H
