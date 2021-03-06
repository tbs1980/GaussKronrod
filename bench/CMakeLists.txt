cmake_minimum_required(VERSION 2.8)

project(gauss_konrad_tests CXX)

if(UNIX)
	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -Wfatal-errors -g")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
	endif(CMAKE_COMPILER_IS_GNUCXX)
	#is the compiler clang?
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -g")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++11 -pedantic -w3 -g")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++98 -O3 -march=native ")
	endif()
endif(UNIX)

#Eigen3 headers
find_package(Eigen3 REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MPFRCPP REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${GMP_INCLUDES})
include_directories(${MPFR_INCLUDES})
include_directories(${MPFRCPP_INCLUDES})

# our headers
include_directories(${CMAKE_SOURCE_DIR})

# libs
set(ALL_LIBS ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

add_executable(bench_PH_mpkronrod bench_PH_mpkronrod.cpp)
target_link_libraries(bench_PH_mpkronrod ${ALL_LIBS})
