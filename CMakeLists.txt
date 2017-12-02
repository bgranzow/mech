cmake_minimum_required(VERSION 3.0.0)

project(MECH VERSION 1.0.0 LANGUAGES CXX)

include(bob.cmake)
bob_begin_package()
bob_begin_cxx_flags()
bob_cxx11_flags()
bob_end_cxx_flags()

set(MECH_USE_SCOREC_DEFAULT ON)
bob_public_dep(SCOREC)

set(MECH_USE_Trilinos_DEFAULT ON)
bob_public_dep(Trilinos)

list(REVERSE Trilinos_INCLUDE_DIRS)
list(REVERSE Trilinos_TPL_INCLUDE_DIRS)
list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_INCLUDE_DIRS)
list(REVERSE Trilinos_TPL_INCLUDE_DIRS)
list(REVERSE Trilinos_TPL_LIBRARIES)

set(SOURCES
control.cpp
disc.cpp
residual.cpp
)

set(HEADERS
control.hpp
disc.hpp
residual.hpp
)

add_library(mech ${SOURCES})
target_include_directories(mech PUBLIC
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${PETSC_INCLUDE_DIR}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(mech PUBLIC
  SCOREC::core
  teuchoscomm
  teuchoscore
  teuchosparameterlist
  teuchosparser
  teuchosnumerics
  minitensor
  sacado
  ${Trilinos_TPL_LIBRARIES}
  ${PETSC_LIBRARY})
bob_export_target(mech)
install(FILES ${HEADERS} DESTINATION include)

function(add_exe exename srcname)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} PRIVATE mech)
endfunction()

function(serial_test testname exe)
  add_test(
    NAME ${testname}
    COMMAND "./${exe}" ${ARGN})
endfunction()

function(mpi_test testname exe)
  add_test(
    NAME ${testname}
    COMMAND "mpirun" "-np" "4" "./${exe}" ${ARGN})
endfunction()

add_exe(ex1 ex1.cpp)
add_exe(ex2 ex2.cpp)

set(square_mdl "../mesh/square/square.dmg")
set(square_assoc "../mesh/square/square.txt")
set(square_msh_1p "../mesh/square/square-serial.smb")
set(square_msh_4p "../mesh/square/square.smb")

set(cube_mdl "../mesh/cube/cube.dmg")
set(cube_assoc "../mesh/cube/cube.txt")
set(cube_msh_1p "../mesh/cube/cube-serial.smb")
set(cube_msh_4p "../mesh/cube/cube.smb")

serial_test(ex1_square_1p ex1 ${square_mdl} ${square_msh_1p} ${square_assoc})
serial_test(ex1_cube_1p ex1 ${cube_mdl} ${cube_msh_1p} ${cube_assoc})

mpi_test(ex1_square_4p ex1 ${square_mdl} ${square_msh_4p} ${square_assoc})
mpi_test(ex1_cube_4p ex1 ${cube_mdl} ${cube_msh_4p} ${cube_assoc})

bob_end_package()
