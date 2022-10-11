#find_package(OpenMP REQUIRED)
#include_directories(/usr/local/include)
#link_directories(/usr/local/lib)

# Project to apply error estimation for hybrid H1 approximations
add_pz_target(
  NAME HybridH1
SOURCES main_HybridH1.cpp)
#add_executable(HybridH1 main_HybridH1.cpp)
target_link_libraries(HybridH1 PRIVATE Methods Tools)

add_pz_target(
  NAME Automated
SOURCES main_Automated.cpp)
#add_executable(Automated main_Automated.cpp)
target_link_libraries(Automated PRIVATE Methods Tools)

add_pz_target(
  NAME TimeEfficiency
SOURCES main_TimeEfficiency.cpp)
#add_executable(TimeEfficiency main_TimeEfficiency.cpp)
target_link_libraries(TimeEfficiency PRIVATE Methods Tools)

add_pz_target(
  NAME ErrorMeasurement
SOURCES main_ErrorMeasurement.cpp)
#add_executable(ErrorMeasurement main_ErrorMeasurement.cpp)
target_link_libraries(ErrorMeasurement PRIVATE Methods Tools)

add_pz_target(
        NAME HybridizedMixed
        SOURCES main_HybridizedMixed.cpp)
target_link_libraries(HybridizedMixed PRIVATE Methods Tools)

add_pz_target(
        NAME PragerSynge
        SOURCES main_PragerSynge.cpp)
target_link_libraries(PragerSynge PRIVATE Methods Tools)


# Project to apply error estimation for hybrid H1 approximations
add_pz_target(
  NAME Mixed
SOURCES main_Mixed.cpp)

#add_executable(Mixed main_Mixed.cpp)
target_link_libraries(Mixed PRIVATE Methods Tools)
