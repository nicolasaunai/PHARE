#ifndef PHARE_ELECTRONS_HPP
#define PHARE_ELECTRONS_HPP

#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
// #include "core/data/grid/gridlayoutdefs.hpp"
// #include "core/utilities/index/index.hpp"
#include "core/def.hpp"
// #include "core/logger.hpp"

#include "initializer/data_provider.hpp"
#include <memory>

#include "core/data/field/initializers/field_user_initializer.hpp"


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

    NO_DISCARD bool isUsable() const { return ions_.isUsable() && J_.isUsable() && Ve_.isUsable(); }

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
        auto const& Ni  = ions_.density(); // gives the particle density, hence the electron density

        auto& Vex = Ve_(Component::X);
        auto& Vey = Ve_(Component::Y);
        auto& Vez = Ve_(Component::Z);

        // from Ni because all components are defined on primal
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
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

public:
    ElectronPressureClosure(PHARE::initializer::PHAREDict const& dict, FluxComputer const& flux, VecField const& B)
        : flux_{flux}
        , B_{B}
        , Pe_{"Pe", HybridQuantity::Scalar::P}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD virtual bool isUsable() const { return Pe_.isUsable() and B_.isUsable() and flux_.isUsable(); }

    NO_DISCARD virtual bool isSettable() const { return Pe_.isSettable() and B_.isUsable() and flux_.isUsable(); }  // TODO flux was not usable before

    struct PressureProperty
    {
        std::string name;
        typename HybridQuantity::Scalar qty;
    };

    using PressureProperties = std::vector<PressureProperty>;

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(flux_, B_, Pe_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(flux_, B_, Pe_); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    void virtual initialize(GridLayout const& layout) = 0;

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
    VecField B_;
    Field Pe_;
};



template<typename FluxComputer>
class IsothermalElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

public:
    using field_type = Field;

    IsothermalElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& flux, VecField const& B)
        : Super{dict, flux, B}
        , Te_{dict["pressure_closure"]["Te"].template to<double>()}
    {
    }

    void initialize(GridLayout const& layout) override {}

    void computePressure(GridLayout const& /*layout*/) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - isothermal closure pressure is not usable");

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
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

    static constexpr std::size_t dim = VecField::dimension;


public:
    using field_type = Field;

    PolytropicElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& flux, VecField const& B)
        : Super{dict, flux, B}
        , gamma_{dict["pressure_closure"]["Gamma"].template to<double>()}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
    {
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    void computePressure(GridLayout const& /*layout*/) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - polytropic pressure closure not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [this](auto n) { return n * 0.1 + 0. * gamma_; });
    }

private:
    double const gamma_ = 5./3.;
    initializer::InitFunction<dim> Pe_init_;
};



template<typename FluxComputer>
class CGLElectronPressureClosure : public ElectronPressureClosure<FluxComputer>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

    using Super = ElectronPressureClosure<FluxComputer>;
    using Super::getCompileTimeResourcesViewList;
    using Super::isSettable;
    using Super::isUsable;

    static constexpr std::size_t dim = VecField::dimension;


public:
    using field_type = Field;

    CGLElectronPressureClosure(PHARE::initializer::PHAREDict const& dict,
                                      FluxComputer const& flux, VecField const& B)
        : Super{dict, flux, B}
        , gamma_{dict["pressure_closure"]["Gamma"].template to<double>()}
        , Pe_init_{dict["pressure_closure"]["Pe"].template to<initializer::InitFunction<dim>>()}
    {
    }

    void initialize(GridLayout const& layout) override
    {
        FieldUserFunctionInitializer::initialize(this->Pe_, layout, Pe_init_);
    }

    void computePressure(GridLayout const& /*layout*/) override
    {
        static_assert(Field::is_contiguous, "Error - assumes Field date is contiguous");

        if (!this->Pe_.isUsable())
            throw std::runtime_error("Error - CGL pressure closure not usable");

        auto const& Ne_ = this->flux_.density();
        std::transform(std::begin(Ne_), std::end(Ne_), std::begin(this->Pe_),
                       [this](auto n) { return n * 0.1; });
    }

private:
    double const gamma_ = 5./3.;
    initializer::InitFunction<dim> Pe_init_;
};



template<typename FluxComputer>
std::unique_ptr<ElectronPressureClosure<FluxComputer>>
ElectronPressureClosureFactory(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux, typename FluxComputer::VecField const& B)
{
    if (dict["pressure_closure"]["name"].template to<std::string>() == "isothermal")
    {
        return std::make_unique<IsothermalElectronPressureClosure<FluxComputer>>(dict, flux, B);
    }
    else if (dict["pressure_closure"]["name"].template to<std::string>() == "polytropic")
    {
        return std::make_unique<PolytropicElectronPressureClosure<FluxComputer>>(dict, flux, B);
    }
    else if (dict["pressure_closure"]["name"].template to<std::string>() == "CGL")
    {
        return std::make_unique<CGLElectronPressureClosure<FluxComputer>>(dict, flux, B);
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
    ElectronMomentModel(PHARE::initializer::PHAREDict const& dict, FluxComputer& flux, VecField const& B)
        : dict_{dict}
        , fluxComput_{flux}
        , B_{B}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict, flux, B)}
    {
    }

    ElectronMomentModel(ElectronMomentModel const& self)
        : dict_{self.dict_}
        , fluxComput_{self.fluxComput_}
        , B_{self.B_}
        , pressureClosure_{ElectronPressureClosureFactory<FluxComputer>(dict_, fluxComput_, B_)}
    {
        *pressureClosure_ = *self.pressureClosure_;
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        // return fluxComput_.isUsable() and pressureClosure_->isUsable();
        return fluxComput_.isUsable() and B_.isUsable() and pressureClosure_->isUsable();
    }

    // NO_DISCARD bool isSettable() const { return fluxComput_.isSettable(); }  // TODO pressureClosure needs also to be settable ?
    NO_DISCARD bool isSettable() const { return fluxComput_.isSettable() and B_.isSettable() and pressureClosure_->isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        // return std::forward_as_tuple(fluxComput_, *pressureClosure_);
        return std::forward_as_tuple(fluxComput_, B_, *pressureClosure_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        // return std::forward_as_tuple(fluxComput_, *pressureClosure_);
        return std::forward_as_tuple(fluxComput_, B_, *pressureClosure_);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    void initialize(GridLayout const& layout) { pressureClosure_->initialize(layout); }

    NO_DISCARD Field const& density() const { return fluxComput_.density(); }
    NO_DISCARD VecField const& velocity() const { return fluxComput_.velocity(); }
    NO_DISCARD Field const& pressure() const { return pressureClosure_->pressure(); }

    NO_DISCARD Field& density() { return fluxComput_.density(); }
    NO_DISCARD VecField& velocity() { return fluxComput_.velocity(); }
    NO_DISCARD Field& pressure() { return pressureClosure_->pressure(); }

    void computeDensity() { fluxComput_.computeDensity(); }
    void computeBulkVelocity(GridLayout const& layout) { fluxComput_.computeBulkVelocity(layout); }
    void computePressure(GridLayout const& layout) { pressureClosure_->computePressure(layout); }

private:
    initializer::PHAREDict dict_;
    FluxComputer fluxComput_;
    VecField B_;
    std::unique_ptr<ElectronPressureClosure<FluxComputer>> pressureClosure_;
};



template<typename FluxComputer>
class Electrons : public LayoutHolder<typename FluxComputer::GridLayout>
{
    using GridLayout = typename FluxComputer::GridLayout;
    using VecField   = typename FluxComputer::VecField;
    using Field      = typename FluxComputer::Field;

public:
    Electrons(initializer::PHAREDict const& dict, FluxComputer flux, VecField const& B)
        : dict_{dict}
        , momentModel_{dict, flux, B}
    {
    }

    Electrons(Electrons const& that) = default;

    void initialize(GridLayout const& layout) { momentModel_.initialize(layout); }

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

    NO_DISCARD bool isUsable() const { return momentModel_.isUsable(); }

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
