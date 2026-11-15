
#include "core/def/phlop.hpp" // IWYU pragma: keep // scope timing
#include "mpi/mpi_utils.hpp"
#include "core/utilities/parallel.hpp"

#include "initializer/data_provider.hpp"

#include "samrai.hpp"

#include <SAMRAI/tbox/SAMRAIManager.h>

#include <stdexcept>

namespace PHARE
{

SamraiLifeCycle::SamraiLifeCycle(int argc, char** argv)
{
    // MPI is told that the process has threads, of which only this one calls it: the threads of
    // core::parallel_for never do
    int provided = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED and core::nbrThreadsFromEnv() > 1)
        throw std::runtime_error(
            "PHARE_THREADS: this MPI cannot be used by a process with threads");
    SAMRAI::tbox::SAMRAI_MPI::init(MPI_COMM_WORLD);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
    // uncomment next line for debugging samrai issues
    // SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort();
    std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
        = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
    SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    PHARE_WITH_PHLOP({
        if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::scope_timer()
                .file_name(".phare/timings/rank." + std::to_string(mpi::rank()) + ".txt")
                .init();
    })
}

SamraiLifeCycle::~SamraiLifeCycle()
{
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
    MPI_Finalize();
}

void SamraiLifeCycle::reset()
{
    PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::startup();
    getRestartManager()->clearRestartItems();
}


SAMRAI::hier::VariableDatabase* SamraiLifeCycle::getDatabase()
{
    return SAMRAI::hier::VariableDatabase::getDatabase();
}

SAMRAI::hier::PatchDataRestartManager* SamraiLifeCycle::getPatchDataRestartManager()
{
    return SAMRAI::hier::PatchDataRestartManager::getManager();
}


SAMRAI::tbox::RestartManager* SamraiLifeCycle::getRestartManager()
{
    return SAMRAI::tbox::RestartManager::getManager();
}


} // namespace PHARE
