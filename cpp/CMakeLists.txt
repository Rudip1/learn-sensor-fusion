add_library(sensor_fusion STATIC
    src/nmea.cpp
    src/gnss.cpp
    src/geodesy.cpp
    src/gnss_quality.cpp
    src/cloud.cpp
    src/sim.cpp
)
add_library(sensor_fusion::sensor_fusion ALIAS sensor_fusion)
target_include_directories(sensor_fusion PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(sensor_fusion PUBLIC Eigen3::Eigen)
target_compile_features(sensor_fusion PUBLIC cxx_std_17)
if(MSVC)
    target_compile_options(sensor_fusion PRIVATE /W4)
else()
    target_compile_options(sensor_fusion PRIVATE -Wall -Wextra -Wpedantic)
endif()
