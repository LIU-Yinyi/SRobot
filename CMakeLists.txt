cmake_minimum_required(VERSION 3.9)
project(SRobot VERSION 2.1.0 DESCRIPTION "Scientific Robot Mathematic Development Library")
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
add_definitions("-DCALCULATE_MODE -DPLATFORM=1 -fPIC -O3")

FILE(GLOB SRCS *.cpp SType/*.cpp)
FILE(GLOB INCS *.h)
FILE(GLOB STYPE_INCS SType/*.h)

add_library(${PROJECT_NAME} SHARED ${SRCS})
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
    
configure_file(SRobot.pc.in SRobot.pc @ONLY)
target_include_directories(${PROJECT_NAME} PUBLIC . ./SType)
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_BINARY_DIR}/SRobot.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
install(FILES ${INCS}
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SRobot)
install(FILES ${STYPE_INCS}
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SRobot/SType)
