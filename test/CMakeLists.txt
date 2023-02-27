
set(This Tests)

set(Sources 
    test_spin.cpp)

add_executable(${This} ${Sources})


target_link_libraries(${This} PUBLIC
                        gtest_main
)   

add_test(NAME ${This}
        COMMAND ${This})