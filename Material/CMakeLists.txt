# Adds materials to ErrorEstimationLib
target_sources(Tools PUBLIC
    LCC_LagrangeMultiplier.h
    LCC_LagrangeMultiplier.cpp
    LCC_H1MixPSMaterial.h
    LCC_H1MixPSMaterial.cpp
)

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_LIST_DIR})
