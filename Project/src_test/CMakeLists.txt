#copiare la riga seguente per aggiungere file .cpp
#list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/<nome file>")

#copiare la riga seguente per aggiungere file .hpp
#list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/<nome file>")

set(src_test_sources ${src_sources} PARENT_SCOPE)
set(src_test_headers ${src_headers} PARENT_SCOPE)
set(src_test_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
