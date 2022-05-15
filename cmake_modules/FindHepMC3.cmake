find_package(HepMC3 NO_MODULE)
if(NOT HepMC3_FOUND)
        return()
endif()

add_library(HepMC3 IMPORTED INTERFACE)

set_target_properties(HepMC3
        PROPERTIES
        INTERFACE_LINK_LIBRARIES "${HEPMC3_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${HEPMC3_INCLUDE_DIR}")