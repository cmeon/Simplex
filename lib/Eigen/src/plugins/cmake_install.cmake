# Install script for directory: /home/cmeon/SimplexLP/eigen/Eigen/src/plugins

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/eigen3/Eigen/src/plugins/MatrixCwiseBinaryOps.h;/usr/local/include/eigen3/Eigen/src/plugins/MatrixCwiseUnaryOps.h;/usr/local/include/eigen3/Eigen/src/plugins/ArrayCwiseBinaryOps.h;/usr/local/include/eigen3/Eigen/src/plugins/CommonCwiseBinaryOps.h;/usr/local/include/eigen3/Eigen/src/plugins/CommonCwiseUnaryOps.h;/usr/local/include/eigen3/Eigen/src/plugins/ArrayCwiseUnaryOps.h;/usr/local/include/eigen3/Eigen/src/plugins/BlockMethods.h")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/include/eigen3/Eigen/src/plugins" TYPE FILE FILES
    "/home/cmeon/SimplexLP/eigen/Eigen/src/plugins/MatrixCwiseBinaryOps.h"
    "/home/cmeon/SimplexLP/eigen/Eigen/src/plugins/MatrixCwiseUnaryOps.h"
    "/home/cmeon/SimplexLP/eigen/Eigen/src/plugins/ArrayCwiseBinaryOps.h"
    "/home/cmeon/SimplexLP/eigen/Eigen/src/plugins/CommonCwiseBinaryOps.h"
    "/home/cmeon/SimplexLP/eigen/Eigen/src/plugins/CommonCwiseUnaryOps.h"
    "/home/cmeon/SimplexLP/eigen/Eigen/src/plugins/ArrayCwiseUnaryOps.h"
    "/home/cmeon/SimplexLP/eigen/Eigen/src/plugins/BlockMethods.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

