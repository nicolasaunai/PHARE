#ifndef PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
#define PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP

#include "core/utilities/parallel.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include <tuple>
#include <string>
#include <vector>
#include <cstddef>
#include <utility>
#include <type_traits>

namespace PHARE::amr
{

// patches of a level with the given resources set on each patch in turn
template<typename GridLayout, typename ResourcesManager, typename Level, typename... Args>
struct LevelAccessor
{
    using GridLayout_t = GridLayout;
    using level_t      = Level;

    LevelAccessor(ResourcesManager& rm, level_t& level, Args&... args_in)
        : rm_{rm}
        , level_{level}
        , args{args_in...}
    {
    }

    auto size() const { return static_cast<std::size_t>(level_.getLocalNumberOfPatches()); }

    // fn(view) per patch, the resources are set on the patch for the duration of the call
    void exec(auto&& fn)
    {
        std::apply(
            [&](auto&... resources) {
                for (auto const& patch : rm_.enumerate(level_, resources...))
                {
                    PatchView view{rm_, *patch, args};
                    fn(view);
                }
            },
            args);
    }

    // fn(view) per patch, patches running concurrently on the threads of core::parallel_for.
    // The model has one set of resource views, which can only point at one patch at a time, so
    // each thread sets its own copy of them: fn must reach the resources through view.args only,
    // and must only write the data of its own patch. Same as exec when one thread is used.
    void concurrentExec(auto&& fn)
    {
        auto const nthreads = core::nbrThreads();
        if (nthreads == 1)
            return exec(fn);

        // made once by the resources manager: per call they would cost this thread, alone at
        // work here, a time that grows with the number of threads
        auto& copies = std::apply(
            [&](auto&... resources) -> auto& { return rm_.viewCopies(nthreads, resources...); },
            args);

        std::vector<patch_t const*> patches; // owned by the level
        patches.reserve(size());
        for (auto const& patch : level_)
            patches.push_back(patch.get());

        core::parallel_for(patches.size(), [&](std::size_t const i) {
            std::apply(
                [&](auto&... resources) {
                    auto const& patch  = *patches[i];
                    auto const onPatch = rm_.setOnPatch(patch, resources...);
                    PatchView view{rm_, patch, std::tuple<Args&...>{resources...}};
                    fn(view);
                },
                copies[core::threadIndex()]);
        });
    }

private:
    using patch_t = std::remove_cvref_t<decltype(*std::declval<level_t&>().getPatch(0))>;

    // only valid during exec, the patch is owned by the level
    struct PatchView
    {
        PatchView(ResourcesManager& rm, patch_t const& patch, std::tuple<Args&...> const& args_in)
            : id{to_string(patch.getGlobalId())}
            , layout{amr::layoutFromPatch<GridLayout_t>(patch)}
            , args{args_in}
            , rm_{rm}
            , patch_{patch}
        {
        }

        auto& patchID() const { return id; }

        void setTime(double const time, auto&... quantities) const
        {
            (rm_.setTime(quantities, patch_, time), ...);
        }

        // cells of this patch owned by another patch of the level, see makeNonOwnedBoxesFor
        auto nonOwnedBoxes(auto const& neighbours) const
        {
            return makeNonOwnedBoxesFor<GridLayout_t::dimension>(patch_, neighbours);
        }

        // this patch domain and its particle ghost box overlaps with the patches of the level, see
        // makeNonLevelGhostBoxFor
        auto nonLevelGhostBoxes(auto const& neighbours) const
        {
            return makeNonLevelGhostBoxFor<GridLayout_t>(patch_, neighbours);
        }

        // patch data not managed as a resource, e.g. the load balancer workload
        template<typename PatchData>
        auto& getPatchDataAs(int const id) const
        {
            return dynamic_cast<PatchData&>(*patch_.getPatchData(id));
        }

        std::string id;
        GridLayout_t layout;
        std::tuple<Args&...> args;

    private:
        ResourcesManager& rm_;
        patch_t const& patch_;
    };

public:
    using PatchView_t = PatchView;

private:
    ResourcesManager& rm_;
    level_t& level_;
    std::tuple<Args&...> args;
};

template<typename GridLayout, typename ResourcesManager, typename Level>
auto make_level_accessor(Level& level, ResourcesManager& rm, auto&... args)
{
    return LevelAccessor<GridLayout, ResourcesManager, Level,
                         std::remove_reference_t<decltype(args)>...>{rm, level, args...};
}

template<typename Model, typename Level>
auto make_model_level_accessor(Level& level, Model& model, auto&... args)
    requires std::is_same_v<std::remove_const_t<Level>, typename Model::amr_types::level_t>
{
    return make_level_accessor<typename Model::gridlayout_type>(level, *model.resourcesManager,
                                                                args...);
}

// fn(view) has run for every patch of the level on return
template<typename Accessor>
void onLevel(Accessor&& accessor, auto&& fn)
{
    accessor.exec(fn);
}

// as onLevel with the patches running concurrently, see LevelAccessor::concurrentExec for what
// fn may do
template<typename Accessor>
void concurrentOnLevel(Accessor&& accessor, auto&& fn)
{
    accessor.concurrentExec(fn);
}

// sum of fn(view) over the patches of the level on this rank
template<typename Accessor>
auto sumOnLevel(Accessor&& accessor, auto&& fn)
{
    using PatchView_t = std::remove_cvref_t<Accessor>::PatchView_t;
    std::decay_t<decltype(fn(std::declval<PatchView_t&>()))> sum{};
    accessor.exec([&](auto& view) { sum += fn(view); });
    return sum;
}

} // namespace PHARE::amr

#endif // PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
