#ifndef PHARE_MHD_STATE_HPP
#define PHARE_MHD_STATE_HPP

#include "amr/data/field/initializers/field_user_initializer.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"
#include "core/def.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/models/physical_state.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE {
namespace core {
class MHDStateInitializer : public PhysicalStateInitializer {};

template <typename VecFieldT>
class MHDState : public IPhysicalState {
    using field_type = typename VecFieldT::field_type;

   public:
    static constexpr auto dimension = VecFieldT::dimension;

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const {
        return rho.isUsable() and V.isUsable() and B.isUsable() and P.isUsable() and
               M.isUsable() and Etot.isUsable() and J.isUsable();
    }

    NO_DISCARD bool isSettable() const {
        return rho.isSettable() and V.isSettable() and B.isSettable() and P.isSettable() and
               M.isSettable() and Etot.isSettable() and J.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const {
        return std::forward_as_tuple(rho, V, B, P, M, Etot, J);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() {
        return std::forward_as_tuple(rho, V, B, P, M, Etot, J);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    MHDState(PHARE::initializer::PHAREDict const& dict)
        : rho{"rho", MHDQuantity::Scalar::rho},
          V{"V", MHDQuantity::Vector::V},
          B{"B", MHDQuantity::Vector::B},
          P{"P", MHDQuantity::Scalar::P}

          ,
          M{"M", MHDQuantity::Vector::M},
          Etot{"Etot", MHDQuantity::Scalar::Etot}

          ,
          J{"J", MHDQuantity::Vector::J}

          ,
          rhoinit_{
              dict["density"]["initializer"].template to<initializer::InitFunction<dimension>>()},
          Vinit_{dict["velocity"]["initializer"]},
          Binit_{dict["magnetic"]["initializer"]},
          Pinit_{
              dict["pressure"]["initializer"].template to<initializer::InitFunction<dimension>>()} {
    }

    template <typename GridLayout>
    void initialize(GridLayout const& layout) {
        FieldUserFunctionInitializer::initialize(rho, layout, rhoinit_);
        Vinit_.initialize(V, layout);
        Binit_.initialize(B, layout);
        FieldUserFunctionInitializer::initialize(P, layout, Pinit_);
    }

    field_type rho;
    VecFieldT  V;
    VecFieldT  B;
    field_type P;

    VecFieldT  M;
    field_type Etot;

    VecFieldT J;

   private:
    initializer::InitFunction<dimension> rhoinit_;
    VecFieldInitializer<dimension>       Vinit_;
    VecFieldInitializer<dimension>       Binit_;
    initializer::InitFunction<dimension> Pinit_;
};
}  // namespace core
}  // namespace PHARE

#endif  // PHARE_MHD_STATE_HPP
