FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(Catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

add_executable(sensor_fusion_tests
    test_nmea.cpp
    test_gnss.cpp
    test_geodesy.cpp
    test_gnss_quality.cpp
    test_cloud.cpp
)
target_link_libraries(sensor_fusion_tests PRIVATE sensor_fusion::sensor_fusion Catch2::Catch2WithMain)
target_compile_definitions(sensor_fusion_tests PRIVATE SF_DATA_DIR="${PROJECT_SOURCE_DIR}/data")
catch_discover_tests(sensor_fusion_tests)
