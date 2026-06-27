add_library(rgpot_eon STATIC RGPot.cpp RGPotEngine.cpp)
target_include_directories(rgpot_eon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(rgpot_eon PUBLIC WITH_RGPOT)
target_link_libraries(rgpot_eon PUBLIC ${CMAKE_DL_LIBS})
