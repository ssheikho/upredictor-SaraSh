cmake_minimum_required(VERSION 2.8)

#cmake_policy(SET CMP0054 OLD)
#cmake_policy(SET CMP0045 OLD)

project(EllipseFit)

include(${CMAKE_CURRENT_LIST_DIR}/../ProjectIncludes.cmake.in)
#set(CMAKE_BUILD_TYPE Debug)

#FIND_PACKAGE(OpenCL)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -std=c++0x")

#/usr/lib/x86_64-linux-gnu/libOpenCL.so

include(FindPkgConfig)
pkg_check_modules(GTK REQUIRED "gtk+-3.0")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GTK_CFLAGS}")

set(ellipse_fit_src
	FitFunctions.cpp
	#LinearAlgebra.cpp
	Plane.cpp		
	EllipseImplicitFit.cpp		
	Ellipse.cpp
	Ellipse3D.cpp
	#RigidTrans2D.cpp
	#RigidTransformation.cpp

	EllipseFit.cpp

	BaysFitFunctions.cpp
	#BayesianFit.cpp
)

include_directories(

	${CMAKE_CURRENT_SOURCE_DIR}/../ALGLIB
	${CMAKE_CURRENT_SOURCE_DIR}/../PathPlanning
	../UBCCamera
	../UBCColorTrackCUDA
	../UBCDART
	../UBCGUI
	../UBCKinect/
	../UBCKinectGUI/
	../UBCSimpleGeometry/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCUtil
	../Linux/Third-Party/dart/src
	../Linux/Third-Party/libfreenect2/
	../Linux/Third-Party/libfreenect2/include/
	../Linux/Third-Party/libfreenect2/include/internal
	../Linux/Third-Party/libfreenect2/include/internal/CL
	../Linux/Third-Party/libfreenect2/include/libfreenect2
	/usr/include
	/usr/include/eigen3
	/usr/include/opencv
	/usr/include/opencv2
	/home/justin/Research/ThirdParty/Ceres/ceres-solver/config
	/home/justin/Research/ThirdParty/Ceres/ceres-solver/include
	#/usr/include/atk-1.0
	#/usr/include/cairo
	#/usr/include/gdk-pixbuf-2.0
	#/usr/include/glib-2.0
	#/usr/include/gtk-3.0
	#/usr/include/pango-1.0
	#/usr/lib/x86_64-linux-gnu/glib-2.0/include/
	${OpenCL_INCLUDE_DIRS}
	${GTK_INCLUDE_DIRS}
	/usr/local/cuda-8.0/include
)

link_directories(

	${CMAKE_CURRENT_SOURCE_DIR}/../ALGLIB/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCCamera/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCColorTrackCUDA/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCDART/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCGUI/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCKinect/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCKinectGUI/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCUtil/
	${CMAKE_CURRENT_SOURCE_DIR}/../Linux/Third-Party/dart/lib
	${CMAKE_CURRENT_SOURCE_DIR}/../Linux/Third-Party/libfreenect2/lib/

	../UBCCamera/
	../UBCColorTrackCUDA/
	../UBCDART/
	../UBCGUI/
	../UBCKinect/
	../UBCKinectGUI/
	${CMAKE_CURRENT_SOURCE_DIR}/../UBCUtil/
	../Linux/Third-Party/dart/lib
	../Linux/Third-Party/libfreenect2/lib/
	#/usr/lib/x86_64-linux-gnu/
	#/usr/local/cuda-8.0/lib64
)

add_executable(
	EllipseFit
	${ellipse_fit_src}
)

target_link_libraries(EllipseFit
	#UBColorTrackCuda
	#ubcdart
	#ubckinectgui
	#ubckinect
	#ubccamera

	#dart
	alglib

	ubcutil
	
	opencv_calib3d
	opencv_contrib
	opencv_core
	opencv_features2d
	opencv_flann
	opencv_gpu
	opencv_highgui
	opencv_imgproc
	opencv_legacy
	opencv_ml
	opencv_objdetect
	opencv_ocl
	opencv_photo
	opencv_stitching
	opencv_superres
	opencv_ts
	opencv_video
	opencv_videostab

	glfw
	GL
	GLEW
	GLU
	glut

	#nvvm
	#cudart
	#OpenCL

	${OpenCL_LIBRARIES}
	#freenect2

	pthread

	${GTK_LDFLAGS}

	glog
	gflags
	ceres

	#ubcgui
)

set_property(TARGET EllipseFit PROPERTY CXX_STANDARD 11)
