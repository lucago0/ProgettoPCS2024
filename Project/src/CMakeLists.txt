#copiare la riga seguente per aggiungere file .cpp
#list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/<nome file>")
list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/Test.cpp")
list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/GeometryLibrary.cpp")

#copiare la riga seguente per aggiungere file .hpp
#list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/<nome file>")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/GeometryLibrary.hpp")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/Test.hpp")

set(DFN_sources ${DFN_sources} PARENT_SCOPE)
set(DFN_headers ${DFN_headers} PARENT_SCOPE)
set(DFN_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
