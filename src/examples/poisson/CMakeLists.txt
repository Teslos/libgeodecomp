lgd_generate_sourcelists("./")
set(RELATIVE_PATH "")
include(auto.cmake)

if(WITH_MPI)
  add_executable(libgeodecomp_examples_poisson ${SOURCES})
  set_target_properties(libgeodecomp_examples_poisson PROPERTIES OUTPUT_NAME poisson)
  target_link_libraries(libgeodecomp_examples_poisson ${LOCAL_LIBGEODECOMP_LINK_LIB})
endif()
