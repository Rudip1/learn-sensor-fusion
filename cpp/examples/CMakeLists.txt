foreach(example 01_nmea_summary 02_convert 03_dop_check)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE sensor_fusion::sensor_fusion)
endforeach()
