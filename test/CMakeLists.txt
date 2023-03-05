
set(This Tests)

set(Sources 
    test_spin.cpp
    test_grid.cpp
    )

add_executable(${This} ${Sources})


target_link_libraries(${This} 
                        gtest_main
                
)   
target_link_libraries(${This} 
                        HighFive
                        )

add_test(NAME ${This}
        COMMAND ${This})