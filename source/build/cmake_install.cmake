# Install script for directory: /home/qianrui/github/reconstruction/abacus-develop/source

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/qianrui/intelcompile/testabacus")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_base/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_cell/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_symmetry/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_neighbor/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_orbital/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_md/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_deepks/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/module_esolver/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/src_io/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/src_ions/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/src_lcao/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/src_parallel/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/src_pdiag/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/src_pw/cmake_install.cmake")
  include("/home/qianrui/github/reconstruction/abacus-develop/source/build/src_ri/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/qianrui/github/reconstruction/abacus-develop/source/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
