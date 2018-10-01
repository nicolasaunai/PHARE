#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_H
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_H

#include "utilities/meta/meta_utilities.h"

#include <string>
#include <type_traits>
#include <vector>

namespace PHARE
{
template<typename ResourcesUser, typename Attempt = void>
struct has_field : std::false_type
{
};

template<typename ResourcesUser, typename Attempt = void>
struct has_particles : std::false_type
{
};

template<typename ResourcesUser, typename Attempt = void>
struct has_sub_resources : std::false_type
{
};


template<typename ResourcesUser, typename Attempt = void>
struct has_runtime_subresourceuser_list : std::false_type
{
};


template<typename ResourcesUser, typename Attempt = void>
struct has_compiletime_subresourcesuser_list : std::false_type
{
};



/** \brief has_field is a traits that permit to check if a ResourcesUser
 * has field
 *
 */
template<typename ResourcesUser>
struct has_field<ResourcesUser, tryToInstanciate<decltype(
                                    std::declval<ResourcesUser>().getFieldNamesAndQuantities())>>
    : std::true_type
{
};




/** \brief has_particles is a traits that permit to check if a ResourcesUser
 * has particles
 *
 */
template<typename ResourcesUser>
struct has_particles<ResourcesUser, tryToInstanciate<decltype(
                                        std::declval<ResourcesUser>().getParticleArrayNames())>>
    : std::true_type
{
};



/** \brief has_sub_resources is a traits that permit to check if a ResourcesUser
 * have other ResourcesUser
 *
 */




template<typename ResourcesUser>
struct has_runtime_subresourceuser_list<
    ResourcesUser,
    tryToInstanciate<decltype(std::declval<ResourcesUser>().getRunTimeResourcesUserList())>>
    : std::true_type
{
};


template<typename ResourcesUser>
struct has_compiletime_subresourcesuser_list<
    ResourcesUser,
    tryToInstanciate<decltype(std::declval<ResourcesUser>().getCompileTimeResourcesUserList())>>
    : std::true_type
{
};




/** UseResourcePtr is used to select the resources patch data */
struct UseResourcePtr
{
};


/** UseNullPtr is used to select a nullptr with the correct type */
struct UseNullPtr
{
};


template<typename ResourcesUser>
std::vector<std::string> extractNames(ResourcesUser const& user)
{
    std::vector<std::string> names;

    if constexpr (has_field<ResourcesUser>::value)
    {
        auto properties = user.getFieldNamesAndQuantities();

        for (auto const& property : properties)
        {
            names.push_back(property.name);
        }
    }


    else if constexpr (has_particles<ResourcesUser>::value)
    {
        throw std::runtime_error("not implemeted");
    }

    return names;
};




} // namespace PHARE

#endif // PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_H
