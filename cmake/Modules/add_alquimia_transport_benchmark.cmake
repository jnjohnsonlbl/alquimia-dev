# This function adds a reactive transport benchmark test for Alquimia.
function(add_alquimia_transport_benchmark benchmark input_file)
  set(exe ${PROJECT_BINARY_DIR}/drivers/transport)
  add_test(${benchmark} ${exe} ${input_file})
  set_tests_properties(${benchmark} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()
