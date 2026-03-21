cmake_minimum_required(VERSION 3.14)

find_package(Qt5 REQUIRED COMPONENTS Widgets WebEngineWidgets)

add_executable(prad2qtmon main.cpp)
target_link_libraries(prad2qtmon PRIVATE
    Qt5::Widgets
    Qt5::WebEngineWidgets
)
set_target_properties(prad2qtmon PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    AUTOMOC ON
)
