#ifndef PHARE_MHD_STATE_HPP
#define PHARE_MHD_STATE_HPP

#include "initializer/data_provider.hpp"

#include "core/mhd/mhd_quantities.hpp"
#include "core/models/physical_state.hpp"
#include "core/def.hpp"

namespace PHARE
{
namespace core
{
    class MHDStateInitializer : public PhysicalStateInitializer
    {
    };


    template<typename VecFieldT>
    class MHDState : public IPhysicalState
    {
    public:
        static constexpr auto dimension = VecFieldT::field_type::dimension; 
        
        MHDState(PHARE::initializer::PHAREDict const& dict)
            : rho{dict["rho"]}
            , V{dict["V"]}
            , B{dict["B"]}
            , P{dict["P"]}
            
            , M{"M", MHDQuantity::Vector::M}
            , Etot{"Etot", MHDQuantity::Scalar::Etot}
            
            , J{"J", MHDQuantity::Vector::J}
        {
        }
        
        VecFieldT rho;
        VecFieldT V;
        VecFieldT B;
        VecFieldT P;
        
        VecFieldT M;
        VecFieldT Etot;
        
        VecFieldT J;
        
        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------
        
        NO_DISCARD bool isUsable() const
        {
            return rho.isUsable() and V.isUsable() and B.isUsable() and P.isUsable() and 
                   M.isUsable() and Etot.isUsable() and J.isUsable();
        }
        
        NO_DISCARD bool isSettable() const
        {
            return rho.isSettable() and V.isSettable() and B.isSettable() and P.isSettable() and
                   M.isSettable() and Etot.isSettable() and J.isSettable();
        }
        
        NO_DISCARD auto getCompileTimeResourcesViewList() const
        {
            return std::forward_as_tuple(rho, V, B, P, M, Etot, J);
        }
        
        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            return std::forward_as_tuple(rho, V, B, P, M, Etot, J);
        }

        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

    };
} // namespace core
} // namespace PHARE



#endif // PHARE_MHD_STATE_HPP
