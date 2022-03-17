cmake_minimum_required(VERSION 3.20)
project(imgCodecs)

set(CMAKE_CXX_STANDARD 17)
#Adding QT for CMAKE
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
#add_definitions(-DQT_NO_KEYWORDS)

#include_directories(${OpenCV_INCLUDE_DIRS} "/usr/local/include/opencv4")
#include_directories(/home/anirban/miniconda3/include)
#include_directories("/home/anirban/miniconda3/include/python3.7m")
#include_directories(/usr/include/python3.7m)
#link_directories("/usr/include/python3.8")

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

message(STATUS "OpenCV library status:")
message(STATUS " config: ${OPENCV_DIR}")
message(STATUS " version: ${OpenV_VERSION}")
message(STATUS " libraries: ${OpenCV_LIBS}")
message(STATUS " include_path: ${OpenCV_INCLUDE_DIRS}")

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

#target_link_libraries(${PYTHON_Lib} "/usr/include/python3.8")

message("Python_FOUND:${Python_FOUND}")
message("Python_VERSION:${Python_VERSION}")
message("Python_Development_FOUND:${Python_Development_FOUND}")
message("Python_LIBRARIES:${Python_LIBRARIES}")

include_directories(.)
include_directories(cannylib)
include_directories(niftilib)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

qt5_wrap_ui(UI_HEADERS
            chain_coding.ui
            diffresdialog.ui
            eed_4_4d_dialog.ui
            eedinpaintingdialog.ui
            eedsmoothingdialog.ui
            foeedinpaintingdialog.ui
            imgcodecs.ui
            primizedialog.ui
            randommask.ui
            reggridmaskdialog.ui
            slicewiseeeddialog.ui
            timeaxiscombinedialog.ui
            xorresdialog.ui
            )

add_executable(imgCodecs
        ${UI_HEADERS}
        cannylib/CannyEdgeDetector.cpp
        cannylib/CannyEdgeDetector.h
        niftilib/ctpl_stl.h
        niftilib/nifti1.h
        niftilib/nifti1_io.cpp
        niftilib/nifti1_io.h
        niftilib/znzlib.cpp
        niftilib/znzlib.h
        chain_coding.cpp
        chain_coding.h
        derivatives.cpp
        derivatives.h
        diffresdialog.cpp
        diffresdialog.h
        directneighbordilationdialog.cpp
        directneighbordilationdialog.h
        eed_4_4d_dialog.cpp
        eed_4_4d_dialog.h
        eedinpaintingdialog.cpp
        eedinpaintingdialog.h
        eedsmoothingdialog.cpp
        eedsmoothingdialog.h
        fed.cpp
        fed.h
        fed_kappa.h
        foeedinpaintingdialog.cpp
        foeedinpaintingdialog.h
        imgcodecs.cpp
        imgcodecs.h
        inpainting.cpp
        inpainting.h
        main.cpp
        masks.cpp
        masks.h
        primizedialog.cpp
        primizedialog.h
        r_ilh_0_dialog.cpp
        r_ilh_0_dialog.h
        randommask.cpp
        randommask.h
        reggridmaskdialog.cpp
        reggridmaskdialog.h
        slicewiseeeddialog.cpp
        slicewiseeeddialog.h
        smoothing.cpp
        smoothing.h
        supplementary_functions.cpp
        supplementary_functions.h
        timeaxiscombinedialog.cpp
        timeaxiscombinedialog.h
        xorresdialog.cpp
        xorresdialog.h
        )

target_link_libraries(imgCodecs ${OpenCV_LIBS})
target_link_libraries(imgCodecs ${PYTHON_LIBRARIES})

# Use the Widgets module from Qt 5.
target_link_libraries(imgCodecs Qt5::Widgets)


