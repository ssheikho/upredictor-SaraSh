project(ALGLIB)

set(CMAKE_BUILD_TYPE Debug)

set(alg_lib_src
	statistics.cpp
	specialfunctions.cpp
	solvers.cpp
	optimization.cpp
	linalg.cpp
	interpolation.cpp
	integration.cpp
	fasttransforms.cpp
	diffequations.cpp
	dataanalysis.cpp
	ap.cpp
	alglibmisc.cpp
	alglibinternal.cpp
)

include_directories(
	/usr/include
	/usr/include/eigen3
	#/usr/include/nvidia-361
	/usr/include/opencv
	../Linux/Third-Party/libfreenect2/include/libfreenect2
)

add_library(
	alglib SHARED
	${alg_lib_src}
)

add_library(
	alglib STATIC
	${alg_lib_src}
)

set_property(TARGET alglib PROPERTY CXX_STANDARD 11)
