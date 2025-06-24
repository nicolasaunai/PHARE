#ifndef PHARE_DEBUGOD_HPP
#define PHARE_DEBUGOD_HPP

#include "core/def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "amr/data/field/field_data.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/utilities/types.hpp"

#include <SAMRAI/hier/PatchHierarchy.h>

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <source_location>
#include <filesystem>

namespace PHARE::amr
{


template<typename PHARE_TYPES>
class DEBUGOD
{
public:
    static constexpr auto dimension    = PHARE_TYPES::dimension;
    static constexpr auto interp_order = PHARE_TYPES::interp_order;
    using Grid_t                       = PHARE_TYPES::Grid_t;
    using GridLayout_t                 = PHARE_TYPES::GridLayout_t;
    using FieldData_t                  = PHARE::amr::FieldData<GridLayout_t, Grid_t>;

    using Point_t = PHARE::core::Point<double, dimension>;

    struct GodValue
    {
        Point_t coords;
        std::array<std::uint32_t, dimension> loc_index;
        std::array<int, dimension> amr_index;
        double value;
        int rank;
        std::string patchID;
        std::string name;
        int level;
        std::source_location src_loc;
        double time;
        std::string msg;

        // Add other necessary fields and methods as needed
    };
    using GodExtract = std::unordered_map<std::uint32_t, std::vector<GodValue>>;

    void setHierarchy(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hier)
    {
        hierarchy_ = hier;
    }

    bool isActive() const { return hierarchy_ != nullptr; }

    NO_DISCARD auto getTime(std::string name, SAMRAI::hier::Patch const& patch) const
    {
        auto pdata = getPatchData(patch, name);
        if (pdata)
        {
            return pdata->getTime();
        }
        else
        {
            throw std::runtime_error("Patch data not found for " + name);
        }
    }


    NO_DISCARD auto time_is(std::string name, double time) const
    {
        // take first patch of first level
        // it should be enought to get time
        // this limits to looking at times at coarser time steps for now
        auto patch = *(hierarchy_->getPatchLevel(0)->begin());
        auto pdata = getPatchData(*patch, name);
        return core::float_equals(time, pdata->getTime());
    }

    NO_DISCARD auto inspect(std::string name, Point_t const& lower, Point_t const& upper,
                            std::source_location const location = std::source_location::current(),
                            std::string msg                     = "") const
    {
        GodExtract god_values;
        for (auto ilvl = 0u; ilvl < hierarchy_->getNumberOfLevels(); ++ilvl)
        {
            auto level       = hierarchy_->getPatchLevel(ilvl);
            god_values[ilvl] = std::vector<GodValue>{};

            for (auto& patch : *level)
            {
                if (!is_local(*patch))
                {
                    continue;
                }

                auto& field      = getField(*patch, name);
                auto layout      = layoutFromPatch<GridLayout_t>(*patch);
                auto extract_box = PHARE::core::Box<double, dimension>{lower, upper};
                auto patch_ghost_box
                    = phare_box_from<dimension, int>(getPatchData(*patch, name)->getGhostBox());

                auto centering = GridLayout_t::centering(field.physicalQuantity());


                Box<int, dimension> amr_user_box;
                for (auto i = 0u; i < dimension; ++i)
                {
                    amr_user_box.lower[i]
                        = static_cast<int>((extract_box.lower[i]) / layout.meshSize()[i]);
                    amr_user_box.upper[i]
                        = static_cast<int>((extract_box.upper[i]) / layout.meshSize()[i]);
                }


                auto intersected_box = patch_ghost_box * amr_user_box;

                if (!intersected_box)
                {
                    // std::cout << "boxes :\n"
                    //           << "patch_ghost_box: " << patch_ghost_box << "\n"
                    //           << "amr_user_box: " << amr_user_box << "\n"
                    //           << "intersected_box: " << *intersected_box << "\n";
                    continue;
                }

                auto local_box = layout.AMRToLocal(*intersected_box);


                // loop on nodes
                // given the mesh_size_ on root level
                // it is easy to get the level mesh size
                // and given the lower/upper bounds of the coordinates
                // it's easy to iterate over all nodes
                // these if constexpr may be removable
                // with the FieldBox object maybe....

                GodValue gval;

                if constexpr (dimension == 1)
                {
                    //
                }

                else if constexpr (dimension == 2)
                {
                    auto ixStart = local_box.lower[core::dirX];
                    auto ixEnd   = local_box.upper[core::dirX];
                    auto iyStart = local_box.lower[core::dirY];
                    auto iyEnd   = local_box.upper[core::dirY];
                    // std::cout << "ixStart: " << ixStart << " ixEnd: " << ixEnd
                    //           << " iyStart: " << iyStart << " iyEnd: " << iyEnd << "\n";
                    // std::cout << "amr_user_box: " << amr_user_box
                    //           << " intersected_box: " << *intersected_box
                    //           << " local_box: " << local_box << "\n";

                    for (auto ix = ixStart; ix <= ixEnd; ++ix)
                    {
                        for (auto iy = iyStart; iy <= iyEnd; ++iy)
                        {
                            gval.coords = {
                                layout.meshSize()[0]
                                    * (ix + patch_ghost_box.lower[core::dirX]
                                       + (centering[core::dirX] == PHARE::core::QtyCentering::dual
                                              ? 0.5
                                              : 0)),
                                layout.meshSize()[1]
                                    * (iy + patch_ghost_box.lower[core::dirY]
                                       + (centering[core::dirY] == PHARE::core::QtyCentering::dual
                                              ? 0.5
                                              : 0))};
                            gval.value     = field(ix, iy);
                            gval.patchID   = to_string(patch->getGlobalId());
                            gval.rank      = get_rank(*patch);
                            gval.loc_index = {ix, iy};
                            gval.name      = name;
                            gval.level     = ilvl;
                            gval.src_loc   = location;
                            gval.time      = getTime(name, *patch);
                            gval.msg       = msg;
                            gval.amr_index
                                = {amr_user_box.lower[core::dirX], amr_user_box.lower[core::dirY]};
                            // std::cout << "adding value: " << gval.value
                            //           << " at coords: " << gval.coords.str() << " on patch "
                            //           << gval.patchID << " at rank: " << gval.rank << "\n";

                            god_values[ilvl].push_back(gval);
                        }
                    }
                }
                else if constexpr (dimension == 3)
                {
                    // for (auto& node : intersected_box)
                    // {
                    // }
                }
            }
        }

        return god_values;
    }


    NO_DISCARD auto inspect(std::string name, Point_t const& coord)
    {
        return inspect(name, coord, coord);
    }



    void print(GodExtract const& god_values)
    {
        constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
        for (auto& [ilvl, values] : god_values)
        {
            std::cout << "Level " << ilvl << " with nbr values: " << values.size() << "\n";
            for (auto& v : values)
            {
                auto& coords  = v.coords;
                auto& loc_idx = v.loc_index;
                // auto& amr_idx = v.amr_index;
                auto& rank    = v.rank;
                auto& patchID = v.patchID;
                auto& name    = v.name;
                std::cout << name << " at " << coords.str();
                std::cout << std::setprecision(max_precision);
                std::cout << " = " << v.value << " on L" << v.level;
                std::cout << " Rank: " << rank;
                std::cout << " PatchID: " << patchID;
                std::cout << " at " << std::filesystem::path(v.src_loc.file_name()).filename()
                          << ":" << v.src_loc.line();
                std::cout << " at time: " << v.time;
                std::cout << " at loc_index: (" << loc_idx[0] << ", " << loc_idx[1] << ")";
                std::cout << " at amr index: ( " << v.amr_index[0] << ", " << v.amr_index[1] << ")";
                std::cout << " " << v.msg;
                std::cout << "\n";
            }
        }
    }

    static DEBUGOD<PHARE_TYPES>& INSTANCE()
    {
        static DEBUGOD instance;
        return instance;
    }

    // void stop() { god_.release(); }

    // NO_DISCARD auto& god()
    // {
    //     if (!god_)
    //         init();
    //     return *god_;
    // }

private:
    auto get_rank(SAMRAI::hier::Patch const& patch) const
    {
        return patch.getBox().getBoxId().getOwnerRank();
    }


    bool is_local(SAMRAI::hier::Patch const& patch) const
    {
        return get_rank(patch) == PHARE::core::mpi::rank();
    }

    auto getPatchData(SAMRAI::hier::Patch const& patch, std::string const& name) const
    {
        auto db      = SAMRAI::hier::VariableDatabase::getDatabase();
        auto var_id  = db->getVariable(name);
        auto context = db->getContext("default");
        return patch.getPatchData(var_id, context);
    }

    auto& getField(SAMRAI::hier::Patch const& patch, std::string const& name) const
    {
        auto pdata            = getPatchData(patch, name);
        auto const& fielddata = std::dynamic_pointer_cast<FieldData_t>(pdata);
        return fielddata->field;
    }

    DEBUGOD() {}
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
};
}; // namespace PHARE::amr

#endif
