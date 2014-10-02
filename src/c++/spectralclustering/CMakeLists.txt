cmake_minimum_required(VERSION 2.4)

if(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0005 OLD)
endif(COMMAND CMAKE_POLICY)


# set paths/ lib directories
set(LIB_DIRS "/usr/lib64")
set(FORTRAN_DIRS "/usr/lib/gcc/x86_64-redhat-linux/4.8.2")
set(ARPACK_LIB ${LIB_DIRS}/libarpack.a)
set(LAPACK_LIB ${LIB_DIRS}/liblapack.a)
set(SUPERLU_LIB ${LIB_DIRS}/libsuperlu_4.3.a)
set(BLAS_LIB ${LIB_DIRS}/libblas.a)
set(FORTRAN_LIB ${FORTRAN_DIRS}/libgfortran.a)
set(ALL_LIBS ${ARPACK_LIB} ${LAPACK_LIB} ${SUPERLU_LIB} ${BLAS_LIB} ${FORTRAN_LIB})

foreach(_lib ${ALL_LIBS})
    message(${_lib})
endforeach()

# include directories
set(ARPACKPP_ROOT "/home/amine/arpack++")
set(ARPACKPP_DIR ${ARPACKPP_ROOT}/include)
set(AREIG_DIR ${ARPACKPP_ROOT}/examples/areig)
set(AREIGSYM_DIR ${ARPACKPP_ROOT}/examples/areig/sym)
set(CSCMAT_DIR ${ARPACKPP_ROOT}/examples/matrices/sym)
set(INCLUDE_DIRS ${ARPACKPP_DIR} ${AREIGSYM_DIR} ${AREIG_DIR} ${CSCMAT_DIR})


include_directories(${INCLUDE_DIRS})
# set compiler flags
set( CMAKE_CXX_FLAGS "-g" )

add_executable(test_eig main.cpp)
target_link_libraries(test_eig ${ALL_LIBS})









