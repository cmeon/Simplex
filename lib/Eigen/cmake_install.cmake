# Install script for directory: /home/cmeon/SimplexLP/eigen/Eigen

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
   "/usr/local/include/eigen3/Eigen/StdList;/usr/local/include/eigen3/Eigen/UmfPackSupport;/usr/local/include/eigen3/Eigen/QtAlignedMalloc;/usr/local/include/eigen3/Eigen/CholmodSupport;/usr/local/include/eigen3/Eigen/LU;/usr/local/include/eigen3/Eigen/SparseCore;/usr/local/include/eigen3/Eigen/PaStiXSupport;/usr/local/include/eigen3/Eigen/Jacobi;/usr/local/include/eigen3/Eigen/SparseLU;/usr/local/include/eigen3/Eigen/Eigen;/usr/local/include/eigen3/Eigen/Eigenvalues;/usr/local/include/eigen3/Eigen/IterativeLinearSolvers;/usr/local/include/eigen3/Eigen/LeastSquares;/usr/local/include/eigen3/Eigen/Core;/usr/local/include/eigen3/Eigen/StdDeque;/usr/local/include/eigen3/Eigen/OrderingMethods;/usr/local/include/eigen3/Eigen/QR;/usr/local/include/eigen3/Eigen/Array;/usr/local/include/eigen3/Eigen/Geometry;/usr/local/include/eigen3/Eigen/Householder;/usr/local/include/eigen3/Eigen/StdVector;/usr/local/include/eigen3/Eigen/SuperLUSupport;/usr/local/include/eigen3/Eigen/SparseCholesky;/usr/local/include/eigen3/Eigen/SparseQR;/usr/local/include/eigen3/Eigen/MetisSupport;/usr/local/include/eigen3/Eigen/Eigen2Support;/usr/local/include/eigen3/Eigen/SPQRSupport;/usr/local/include/eigen3/Eigen/SVD;/usr/local/include/eigen3/Eigen/Sparse;/usr/local/include/eigen3/Eigen/PardisoSupport;/usr/local/include/eigen3/Eigen/Cholesky;/usr/local/include/eigen3/Eigen/Dense")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/include/eigen3/Eigen" TYPE FILE FILES
    "/home/cmeon/SimplexLP/eigen/Eigen/StdList"
    "/home/cmeon/SimplexLP/eigen/Eigen/UmfPackSupport"
    "/home/cmeon/SimplexLP/eigen/Eigen/QtAlignedMalloc"
    "/home/cmeon/SimplexLP/eigen/Eigen/CholmodSupport"
    "/home/cmeon/SimplexLP/eigen/Eigen/LU"
    "/home/cmeon/SimplexLP/eigen/Eigen/SparseCore"
    "/home/cmeon/SimplexLP/eigen/Eigen/PaStiXSupport"
    "/home/cmeon/SimplexLP/eigen/Eigen/Jacobi"
    "/home/cmeon/SimplexLP/eigen/Eigen/SparseLU"
    "/home/cmeon/SimplexLP/eigen/Eigen/Eigen"
    "/home/cmeon/SimplexLP/eigen/Eigen/Eigenvalues"
    "/home/cmeon/SimplexLP/eigen/Eigen/IterativeLinearSolvers"
    "/home/cmeon/SimplexLP/eigen/Eigen/LeastSquares"
    "/home/cmeon/SimplexLP/eigen/Eigen/Core"
    "/home/cmeon/SimplexLP/eigen/Eigen/StdDeque"
    "/home/cmeon/SimplexLP/eigen/Eigen/OrderingMethods"
    "/home/cmeon/SimplexLP/eigen/Eigen/QR"
    "/home/cmeon/SimplexLP/eigen/Eigen/Array"
    "/home/cmeon/SimplexLP/eigen/Eigen/Geometry"
    "/home/cmeon/SimplexLP/eigen/Eigen/Householder"
    "/home/cmeon/SimplexLP/eigen/Eigen/StdVector"
    "/home/cmeon/SimplexLP/eigen/Eigen/SuperLUSupport"
    "/home/cmeon/SimplexLP/eigen/Eigen/SparseCholesky"
    "/home/cmeon/SimplexLP/eigen/Eigen/SparseQR"
    "/home/cmeon/SimplexLP/eigen/Eigen/MetisSupport"
    "/home/cmeon/SimplexLP/eigen/Eigen/Eigen2Support"
    "/home/cmeon/SimplexLP/eigen/Eigen/SPQRSupport"
    "/home/cmeon/SimplexLP/eigen/Eigen/SVD"
    "/home/cmeon/SimplexLP/eigen/Eigen/Sparse"
    "/home/cmeon/SimplexLP/eigen/Eigen/PardisoSupport"
    "/home/cmeon/SimplexLP/eigen/Eigen/Cholesky"
    "/home/cmeon/SimplexLP/eigen/Eigen/Dense"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/cmeon/SimplexLP/lib/Eigen/src/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

