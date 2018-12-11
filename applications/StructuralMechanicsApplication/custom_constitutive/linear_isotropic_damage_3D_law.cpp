// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//  Collaborators:
//

// Project includes
#include "linear_isotropic_damage_3D_law.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearIsotropicDamage3D::LinearIsotropicDamage3D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

LinearIsotropicDamage3D::LinearIsotropicDamage3D(const LinearIsotropicDamage3D &rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearIsotropicDamage3D::Clone() const
{
    return Kratos::make_shared<LinearIsotropicDamage3D>(LinearIsotropicDamage3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearIsotropicDamage3D::~LinearIsotropicDamage3D()
{
}

//************************************************************************************
//************************************************************************************

bool LinearIsotropicDamage3D::Has(const Variable<bool>& rThisVariable)
{
    if(rThisVariable == INELASTIC_FLAG){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

bool& LinearIsotropicDamage3D::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    mStrainVariable = yield_stress / std::sqrt(young_modulus);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::InitializeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::InitializeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::InitializeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}
//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::InitializeMaterialResponseKirchhoff(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    double strain_variable;
    this->CalculateStressResponse(rValues, strain_variable);
    mStrainVariable = strain_variable;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}
//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponseKirchhoff(
        Kratos::ConstitutiveLaw::Parameters &rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    double strain_variable = 0;
    this->CalculateStressResponse(rValues, strain_variable);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateStressResponse(
        Parameters& rValues,
        double& rStrainVariable)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rValues.GetOptions();
    Vector& r_strain_vector = rValues.GetStrainVector();
    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
        noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    if( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // If we compute the tangent moduli or the stress
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        Vector& r_stress_vector = rValues.GetStressVector();
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

        CalculateConstitutiveTensor(r_constitutive_matrix, r_material_properties);
        noalias(r_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);

        const double strain_norm = std::sqrt(inner_prod(r_stress_vector, r_strain_vector));
        if (strain_norm <= mStrainVariable)
        {
            // ELASTIC
            mInelasticFlag = false;
            rStrainVariable = mStrainVariable;
            const double stress_variable = EvaluateHardeningLaw(rStrainVariable, r_material_properties);
            const double damage_variable = 1. - stress_variable / rStrainVariable;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                r_constitutive_matrix *= (1 - damage_variable);
            }

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
                r_stress_vector *= (1 - damage_variable);
            }
        }
        else
        {
            // INELASTIC
            mInelasticFlag = true;
            rStrainVariable = strain_norm;
            const double stress_variable = EvaluateHardeningLaw(rStrainVariable, r_material_properties);
            const double damage_variable = 1. - stress_variable / rStrainVariable;
            const double hardening_modulus = r_material_properties[ISOTROPIC_HARDENING_MODULUS];
            const double damage_rate = (stress_variable - hardening_modulus * rStrainVariable)
                                       / (rStrainVariable * rStrainVariable * rStrainVariable);

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                r_constitutive_matrix *= (1. - damage_variable);
                r_constitutive_matrix -= damage_rate * outer_prod(r_stress_vector, r_stress_vector);
            }

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
                r_stress_vector *= (1. - damage_variable);
            }
        }
    }
}

//************************************************************************************
//************************************************************************************

double& LinearIsotropicDamage3D::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == STRAIN_ENERGY){
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        if (rParameterValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(r_strain_vector) += rParameterValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        Matrix& r_constitutive_matrix = rParameterValues.GetConstitutiveMatrix();
        CalculateConstitutiveTensor(r_constitutive_matrix, r_material_properties);
        const double stress_like_variable = EvaluateHardeningLaw(mStrainVariable, r_material_properties);
        const double damage_variable = 1. - stress_like_variable / mStrainVariable;

        rValue = 0.5 * ((1. - damage_variable)
                        * inner_prod(r_strain_vector,
                                     prod(r_constitutive_matrix, r_strain_vector)));
    }

    if (rThisVariable == DAMAGE_VARIABLE){
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        const double stress_like_variable = EvaluateHardeningLaw(mStrainVariable, r_material_properties);

        rValue = 1. - stress_like_variable / mStrainVariable;
    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& LinearIsotropicDamage3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
    {
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        if (rParameterValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(r_strain_vector) += rParameterValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        rValue = r_strain_vector;
    }

    if ( rThisVariable == STRAIN) {
    //if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
    //    rThisVariable == ALMANSI_STRAIN_VECTOR) {

        const SizeType space_dimension = this->WorkingSpaceDimension();

        // Compute total deformation gradient
        const Matrix& r_F = rParameterValues.GetDeformationGradientF();
        KRATOS_DEBUG_ERROR_IF(r_F.size1()!= space_dimension || r_F.size2() != space_dimension)
            << "expected size of F " << space_dimension << "x" << space_dimension
            << ", got " << r_F.size1() << "x" << r_F.size2() << std::endl;

        const Matrix C_tensor = prod(trans(r_F),r_F);
        ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, rValue);
    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

double LinearIsotropicDamage3D::EvaluateHardeningLaw(
        double StrainVariable,
        const Properties &rMaterialProperties
)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double inf_yield_stress = rMaterialProperties[INFINITY_YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double strain_variable_init = yield_stress / std::sqrt(young_modulus);
    const double stress_variable_inf = inf_yield_stress / std::sqrt(young_modulus);
    double stress_variable;
    const double tolerance = std::numeric_limits<double>::epsilon();

    if (StrainVariable < strain_variable_init)
        return StrainVariable;
    stress_variable = strain_variable_init + hardening_modulus * (StrainVariable - strain_variable_init);
    if ((hardening_modulus > tolerance && stress_variable > stress_variable_inf) ||
        (hardening_modulus < tolerance && stress_variable < stress_variable_inf))
        stress_variable = stress_variable_inf;
    return stress_variable;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateConstitutiveTensor(
        Matrix &rConstitutiveTensor,
        const Properties &rMaterialProperties
)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double nu = rMaterialProperties[POISSON_RATIO];

    if (rConstitutiveTensor.size1() != 6 || rConstitutiveTensor.size2() != 6)
        rConstitutiveTensor.resize(6, 6, false);
    rConstitutiveTensor.clear();

    rConstitutiveTensor(0, 0) = (E * (1. - nu) / ((1. + nu) * (1. - 2. * nu)));
    rConstitutiveTensor(1, 1) = rConstitutiveTensor(0, 0);
    rConstitutiveTensor(2, 2) = rConstitutiveTensor(0, 0);
    rConstitutiveTensor(3, 3) = rConstitutiveTensor(0, 0) * (1. - 2. * nu) / (2. * (1. - nu));
    rConstitutiveTensor(4, 4) = rConstitutiveTensor(3, 3);
    rConstitutiveTensor(5, 5) = rConstitutiveTensor(3, 3);
    rConstitutiveTensor(0, 1) = rConstitutiveTensor(0, 0) * nu / (1. - nu);
    rConstitutiveTensor(1, 0) = rConstitutiveTensor(0, 1);
    rConstitutiveTensor(0, 2) = rConstitutiveTensor(0, 1);
    rConstitutiveTensor(2, 0) = rConstitutiveTensor(0, 1);
    rConstitutiveTensor(1, 2) = rConstitutiveTensor(0, 1);
    rConstitutiveTensor(2, 1) = rConstitutiveTensor(0, 1);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = this->GetStrainSize();
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension();
}

//************************************************************************************
//************************************************************************************

int LinearIsotropicDamage3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const double tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(INFINITY_YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS));
    KRATOS_CHECK_GREATER(rMaterialProperties[YIELD_STRESS], tolerance);
    KRATOS_CHECK_GREATER(rMaterialProperties[INFINITY_YIELD_STRESS], tolerance);
    KRATOS_CHECK_LESS(rMaterialProperties[ISOTROPIC_HARDENING_MODULUS], 1.);
    KRATOS_CHECK_NOT_EQUAL(rMaterialProperties[ISOTROPIC_HARDENING_MODULUS], tolerance);

    if (rMaterialProperties[ISOTROPIC_HARDENING_MODULUS] > tolerance &&
        rMaterialProperties[INFINITY_YIELD_STRESS] <= rMaterialProperties[YIELD_STRESS])
        KRATOS_ERROR << "If ISOTROPIC_HARDENING_MODULUS is positive, "
            "INFINITY_YIELD_STRESS must be greater than YIELD_STRESS" << std::endl;

    if (rMaterialProperties[ISOTROPIC_HARDENING_MODULUS] < tolerance &&
        rMaterialProperties[INFINITY_YIELD_STRESS] >= rMaterialProperties[YIELD_STRESS])
        KRATOS_ERROR << "If ISOTROPIC_HARDENING_MODULUS is negative, "
                "INFINITY_YIELD_STRESS must be lesser than YIELD_STRESS" << std::endl;

    return 0;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mStrainVariable", mStrainVariable);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mStrainVariable", mStrainVariable);
}

} /* namespace Kratos.*/
