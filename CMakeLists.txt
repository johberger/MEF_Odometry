project(minimumenergyfilter)

add_executable(unittests Unit_Tests.cpp MinEnFilter.cpp LieAlgebra.cpp LieGroup.cpp common.cpp)


install(TARGETS unittests RUNTIME DESTINATION bin)
