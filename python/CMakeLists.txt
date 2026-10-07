find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE
    sensor_fusion/bindings/module.cpp
    sensor_fusion/bindings/bind_nmea.cpp
    sensor_fusion/bindings/bind_gnss.cpp
    sensor_fusion/bindings/bind_geo.cpp
    sensor_fusion/bindings/bind_cloud.cpp
    sensor_fusion/bindings/bind_sim.cpp
)
target_link_libraries(_core PRIVATE sensor_fusion::sensor_fusion)
install(TARGETS _core LIBRARY DESTINATION sensor_fusion)
