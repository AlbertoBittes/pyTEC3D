/**
 * @file nabla.cpp
 * @brief Nabla: the compiled operator for each case, and the FE data each operator needs.
 */
#include "nabla.h"

#include <libmesh/fe.h>

#include "element_field.h"
#include "logging.h"
#include "nabla_operators.h"

namespace pytec3d
{
    namespace
    {
        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("Nabla: " + message);
        }

        template <typename OutputType, unsigned int Dim, bool Axisymmetric>
        Nabla::Operand identity_function(libMesh::FEAbstract &fe, const Field &field)
        {
            return Identity<OutputType, Dim, Axisymmetric>(static_cast<libMesh::FEGenericBase<OutputType> &>(fe), field);
        }

        template <typename OutputType, unsigned int Dim, bool Axisymmetric>
        Nabla::Operand gradient_function(libMesh::FEAbstract &fe, const Field &field)
        {
            return Gradient<OutputType, Dim, Axisymmetric>(static_cast<libMesh::FEGenericBase<OutputType> &>(fe), field);
        }

        template <typename OutputType, unsigned int Dim, bool Axisymmetric>
        Nabla::Operand curl_function(libMesh::FEAbstract &fe, const Field &field)
        {
            return Curl<OutputType, Dim, Axisymmetric>(static_cast<libMesh::FEGenericBase<OutputType> &>(fe), field);
        }

        template <typename OutputType, unsigned int Dim, bool Axisymmetric>
        Nabla::Operand divergence_function(libMesh::FEAbstract &fe, const Field &field)
        {
            return Divergence<OutputType, Dim, Axisymmetric>(static_cast<libMesh::FEGenericBase<OutputType> &>(fe), field);
        }

        template <typename OutputType, unsigned int Dim, bool Axisymmetric>
        Nabla::Function select(NablaOperator op)
        {
            switch (op)
            {
            case NablaOperator::Identity:
                return &identity_function<OutputType, Dim, Axisymmetric>;
            case NablaOperator::Grad:
                return &gradient_function<OutputType, Dim, Axisymmetric>;
            case NablaOperator::Curl:
                return &curl_function<OutputType, Dim, Axisymmetric>;
            case NablaOperator::Div:
                return &divergence_function<OutputType, Dim, Axisymmetric>;
            }
            fail("unknown operator.");
        }

        template <typename OutputType>
        Nabla::Function select_geometry(NablaOperator op, unsigned int dimension, bool axisymmetric)
        {
            if (dimension == 2)
                return axisymmetric ? select<OutputType, 2, true>(op) : select<OutputType, 2, false>(op);
            if (dimension == 3 && !axisymmetric)
                return select<OutputType, 3, false>(op);
            fail("operators exist for 2D (planar or axisymmetric) and 3D (planar) elements.");
        }
    } // namespace

    Nabla::Nabla(NablaOperator op, const Field &field, unsigned int dimension, bool axisymmetric) : field_(&field), op_(op)
    {
        const std::string name = "'" + std::string(field.name()) + "'";
        if (!field.fe)
            fail("field " + name + " has no finite element.");

        const bool vector_output = field.vector_output;
        switch (op)
        {
        case NablaOperator::Identity:
            break;
        case NablaOperator::Grad:
            if (vector_output && axisymmetric)
                fail("grad of vector field " + name + " is not supported on axisymmetric elements.");
            if (vector_output && field.family != libMesh::LAGRANGE_VEC && field.family != libMesh::L2_LAGRANGE_VEC)
                fail("grad of vector field " + name + " needs a Lagrange vector family.");
            break;
        case NablaOperator::Curl:
            if (vector_output && dimension != 3)
                fail("curl of vector field " + name + " is only supported on 3D elements.");
            if (!vector_output && dimension != 2)
                fail("curl of scalar field " + name + " needs a 2D element.");
            break;
        case NablaOperator::Div:
            if (!vector_output && dimension != 2)
                fail("div of scalar field " + name + " needs a 2D element.");
            break;
        }

        function_ = vector_output ? select_geometry<libMesh::RealGradient>(op, dimension, axisymmetric)
                                  : select_geometry<libMesh::Real>(op, dimension, axisymmetric);

        // Ask the FE for the data this operator reads; it computes them on every reinit().
        if (vector_output)
        {
            auto &vector_fe = static_cast<libMesh::FEVectorBase &>(*field.fe);
            vector_fe.get_phi();
            if (op == NablaOperator::Grad)
                vector_fe.get_dphi();
            else if (op == NablaOperator::Curl)
                vector_fe.get_curl_phi();
            else if (op == NablaOperator::Div)
                vector_fe.get_div_phi();
        }
        else
        {
            auto &scalar_fe = static_cast<libMesh::FEBase &>(*field.fe);
            scalar_fe.get_phi();
            if (op == NablaOperator::Grad || op == NablaOperator::Curl)
                scalar_fe.get_dphi();
        }
        if (axisymmetric)
            field.fe->get_xyz();
    }

    Nabla::Operand Nabla::evaluate() const
    {
        Operand out = function_(*field_->fe, *field_);
        out.scalar = scalar_valued();
        return out;
    }

    bool Nabla::scalar_valued() const
    {
        return op_ == NablaOperator::Div || (op_ == NablaOperator::Identity && !field_->physical_vector);
    }
} // namespace pytec3d
