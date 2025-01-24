#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/index/index.hpp"
#include "core/def.hpp"
#include "core/logger.hpp"

#include "initializer/data_provider.hpp"
#include <memory>


namespace PHARE::core
{


template<typename Ions>
class StandardHybridElectronFluxComputer
{
public:
    using VecField   = typename Ions::vecfield_type;
    using Field      = typename Ions::field_type;
    using GridLayout = typename Ions::gridlayout_type;

    StandardHybridElectronFluxComputer(Ions& ions, VecField& J)
        : ions_{ions}
        , J_{J}
        , Ve_{"StandardHybridElectronFluxComputer_Ve", HybridQuantity::Vector::V}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const {
       // auto io = ions_.isUsable();
       // auto j = J_.isUsable();
       // auto ve = Ve_.isUsable();
       // ;
        return ions_.isUsable() && J_.isUsable() && Ve_.isUsable(); }

    NO_DISCARD bool isSettable() const
    {
        return Ve_.isSettable() && ions_.isSettable() && J_.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(Ve_, ions_, J_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(Ve_, ions_, J_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD Field const& density() const
    {
        if (isUsable())
        {
            return ions_.density();
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }

    NO_DISCARD Field& density()
    {
        if (isUsable())
        {
            return ions_.density();
        }
        else
        {
            throw std::runtime_error("Error, cannot return density because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }

    NO_DISCARD VecField& velocity()
    {
        if (isUsable())
        {
            return Ve_;
        }
        else
        {
            throw std::runtime_error("Error, cannot return velocity because "
                                     "StandardHybridElectronFluxComputer is not usable");
        }
    }

    void computeDensity() {}

    void computeBulkVelocity(GridLayout const& layout)
    {
        auto const& Jx  = J_(Component::X);
        auto const& Jy  = J_(Component::Y);
        auto const& Jz  = J_(Component::Z);
        auto const& Vix = ions_.velocity()(Component::X);
        auto const& Viy = ions_.velocity()(Component::Y);
        auto const& Viz = ions_.velocity()(Component::Z);
        auto const& Ni  = ions_.density();

        auto& Vex = Ve_(Component::X);
        auto& Vey = Ve_(Component::Y);
        auto& Vez = Ve_(Component::Z);

        // from Ni because all components defined on primal
        layout.evalOnBox(Ni, [&](auto const&... args) {
            auto arr = std::array{args...};

            auto const JxOnVx = GridLayout::project(Jx, arr, GridLayout::JxToMoments());
            auto const JyOnVy = GridLayout::project(Jy, arr, GridLayout::JyToMoments());
            auto const JzOnVz = GridLayout::project(Jz, arr, GridLayout::JzToMoments());

            Vex(arr) = Vix(arr) - JxOnVx / Ni(arr);
            Vey(arr) = Viy(arr) - JyOnVy / Ni(arr);
            Vez(arr) = Viz(arr) - JzOnVz / Ni(arr);
        });
    }

    auto& getIons() const { return ions_; }

private:
    Ions ions_;
    VecField J_;
    VecField Ve_;
};



template<typename FluxComputer>
class ElectronPressureClosure
{
    using GridLayout = typename FluxComputer::GridLayout;
    using Field      = typename FluxComputer::Field;

public:

    ElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux)
        : flux_{flux}
        , Pe_{"Pe", HybridQuantity::Scalar::P}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const {
        // auto pe = Pe_.isUsable();
        // auto fl = flux_.isUsable();
        // ;
        return Pe_.isUsable() and flux_.isUsable(); }

    NO_DISCARD bool isSettable() const { return Pe_.isSettable(); }

    struct PressureProperty
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(flux_, Pe_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(flux_, Pe_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD Field& pressure()
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - ! isothermal closure pressure not usable");
        return Pe_;
    }

    NO_DISCARD Field const& pressure() const
    {
        if (!Pe_.isUsable())
            throw std::runtime_error("Error - !! isothermal closure pressure not usable");
        return Pe_;
    }

    void virtual computePressure(GridLayout const& /*layout*/) = 0;

protected:
    FluxComputer flux_;
    Field Pe_;
};



template<typename FluxComputer>
class IsothermalElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = typename FluxComputer::GridLayout;
//  using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

public:
    using field_type = Field;

    IsothermalElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux)
        : Super{dict, flux},
          Te_{dict["pressure_closure"]["Te"].template to<double>()}
    {
    }

    void computePressure(GridLayout const& /*layout*/) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - !!! isothermal closure pressure not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [this](auto n) { return n * Te_; });
    }

private:
    double const Te_ = 0;
};



template<typename FluxComputer>
class PolytropicElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = typename FluxComputer::GridLayout;
//  using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

public:
    using field_type = Field;

    PolytropicElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux)
        : Super{dict, flux},
          Te_{dict["pressure_closure"]["Te_"].template to<double>()}
    {
    }

    void computePressure(GridLayout const& /*layout*/) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - !!! isothermal closure pressure not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [this](auto n) { return n * Te_; });
    }

private:
    double const Te_ = 0;
};



template<typename FluxComputer>
std::unique_ptr<ElectronPressureClosure<FluxComputer>> ElectronPressureClosureFactory(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux)
{

    if (dict["pressure_closure"]["name"].template to<std::string>() == "isothermal")
    {
    return std::make_unique<IsothermalElectronPressureClosure<FluxComputer>>(dict, flux);
    }
    else if (dict["pressure_closure"]["name"].template to<std::string>() == "polytropic")
    {
    return std::make_unique<PolytropicElectronPressureClosure<FluxComputer>>(dict, flux);
    }
    return nullptr;

}



template<typename FluxComputer>
class ElectronMomentModel
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

public:
    ElectronMomentModel(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux)
        : dict_{dict}
        , fluxComput_{flux}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict, flux)}
    {
    }

    ElectronMomentModel(ElectronMomentModel const& self)
        : dict_{self.dict_}
        , fluxComput_{self.fluxComput_}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict_, fluxComput_)}
    {
    *pressureClosure_ = *self.pressureClosure_;
    }
   // ElectronMomentModel(ElectronMomentModel const&) = default;
   // ElectronMomentModel& operator=(ElectronMomentModel const&) = default;
  // ~ElectronMomentModel() = default;

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        // auto fc = fluxComput_.isUsable();
        // auto pc = pressureClosure_->isUsable();
        // ;
        return fluxComput_.isUsable() and pressureClosure_->isUsable();
    }

    NO_DISCARD bool isSettable() const { return fluxComput_.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(fluxComput_, *pressureClosure_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(fluxComput_, *pressureClosure_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD Field const& density() const { return fluxComput_.density(); }
    NO_DISCARD VecField const& velocity() const { return fluxComput_.velocity(); }
    NO_DISCARD Field const& pressure() const { return pressureClosure_->pressure(); }

    NO_DISCARD Field& density() { return fluxComput_.density(); }
    NO_DISCARD VecField& velocity() { return fluxComput_.velocity(); }
    NO_DISCARD Field& pressure() { return pressureClosure_->pressure(); }

    void computeDensity() { fluxComput_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout) { pressureClosure_->computePressure(layout); }

    /*auto static deep_copy(ElectronMomentModel & self, initializer::PHAREDict const& dict)  {*/
    /*  assert(self.isUsable());*/
    /*  // auto const& [Ve, ions, J] = self.fluxComput_.getCompileTimeResourcesViewList( );*/
    /*  ElectronMomentModel<FluxComputer> cpy {dict, self.fluxComput_}; // new shared ptr memory*/
    /*  std::get<0>(cpy.getCompileTimeResourcesViewList()) = self.fluxComput_;*/
    /*  std::get<1>(cpy.getCompileTimeResourcesViewList()) = *self.pressureClosure_;*/
    /*  assert(cpy.isUsable());*/
    /*  return cpy;*/
    /*}*/

private:
    initializer::PHAREDict dict_;
    FluxComputer fluxComput_;
    std::unique_ptr<ElectronPressureClosure<FluxComputer>> pressureClosure_;
};



template<typename FluxComputer>
class Electrons : public LayoutHolder<typename FluxComputer::GridLayout>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    /*auto static copy_or_init_model(Electrons & that){*/
    /**/
    /*  if(that.isUsable()) return ElectronMomentModel<FluxComputer>::deep_copy(that.momentModel_, that.dict_);*/
    /*  // else*/
    /*  auto const& [fluxComput, pressureClosure] = that.momentModel_.getCompileTimeResourcesViewList();*/
    /*  return ElectronMomentModel<FluxComputer>{that.dict_, fluxComput};*/
    /*}*/

public:
    Electrons(initializer::PHAREDict const& dict, FluxComputer flux)
        : dict_{dict}
        , momentModel_{dict, flux}
    {
    }

    // Electrons(Electrons const& that) : dict_{that.dict_}, momentModel_{copy_or_init_model(const_cast<Electrons&>(that))}{}

    Electrons(Electrons const& that) = default;


    void update(GridLayout const& layout)
    {
        if (isUsable())
        {
            momentModel_.computeDensity();
            momentModel_.computeBulkVelocity(layout);
            momentModel_.computePressure(layout);
        }
        else
            throw std::runtime_error("Error - Electron  is not usable");
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const {
        // auto mm = momentModel_.isUsable();
        // ;
      return momentModel_.isUsable(); }

    NO_DISCARD bool isSettable() const { return momentModel_.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(momentModel_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(momentModel_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    NO_DISCARD Field const& density() const { return momentModel_.density(); }
    NO_DISCARD VecField const& velocity() const { return momentModel_.velocity(); }
    NO_DISCARD Field const& pressure() const { return momentModel_.pressure(); }

    NO_DISCARD Field& density() { return momentModel_.density(); }
    NO_DISCARD VecField& velocity() { return momentModel_.velocity(); }
    NO_DISCARD Field& pressure() { return momentModel_.pressure(); }

private:
    initializer::PHAREDict dict_;
    ElectronMomentModel<FluxComputer> momentModel_;
};

} // namespace PHARE::core


#endif
