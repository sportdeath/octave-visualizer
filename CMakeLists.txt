cmake_minimum_required (VERSION 3.4)
project (Visualizer)

#Adding cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/modules/)

# OPENGL
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
set(LIBS ${LIBS} ${OPENGL_LIBRARIES})
if(NOT OPENGL_FOUND)
  message(ERROR “ OPENGL not found!”)
endif(NOT OPENGL_FOUND)

# PortAudio
find_package(PORTAUDIO REQUIRED)
include_directories(${PORTAUDIO_INCLUDE_DIRS})
add_definitions(${PORTAUDIO_DEFINITIONS})
set(LIBS ${LIBS} ${PORTAUDIO_LIBRARIES})
if(NOT PORTAUDIO_FOUND)
  message(ERROR “ PORTAUDIO not found!”)
endif(NOT PORTAUDIO_FOUND)

# FFTW
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
set(LIBS ${LIBS} ${FFTW_LIBRARIES})
if(NOT FFTW_FOUND)
  message(ERROR " FFTW not found!")
endif(NOT FFTW_FOUND)

# SFML
find_package(SFML 2 COMPONENTS system window graphics REQUIRED)
include_directories(${SFML_INCLUDE_DIR})
set(LIBS ${LIBS} ${SFML_LIBRARIES})
if(NOT SFML_FOUND)
  message(ERROR " SFML not found!")
endif(NOT SFML_FOUND)

add_executable(visualizer
	       source/main.cpp
         source/visualizer.cpp
         source/visualizer.hpp
         source/graphics.cpp
         source/graphics.hpp
         source/audio.cpp
         source/audio.hpp
)

target_link_libraries(visualizer
	${LIBS}
)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no c++11 support.")
endif()
