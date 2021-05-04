try_compile(FILESYSTEM_IS_NATIVE
    ${CMAKE_CURRENT_BINARY_DIR}/_filesystem/native
    ${CMAKE_CURRENT_LIST_DIR}/filesystem.cpp
    OUTPUT_VARIABLE FILESYSTEM_NATIVE_ERROR
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

if(NOT FILESYSTEM_IS_NATIVE)
    try_compile(FILESYSTEM_IS_STDCPPFS
        ${CMAKE_CURRENT_BINARY_DIR}/_filesystem/stdcppfs
        ${CMAKE_CURRENT_LIST_DIR}/filesystem.cpp
        LINK_LIBRARIES stdc++fs
        OUTPUT_VARIABLE FILESYSTEM_STDCPPFS_ERROR
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
    )
    if(FILESYSTEM_IS_STDCPPFS)
        set(FILESYSTEM_LIBRARIES stdc++fs)
    else()
        try_compile(FILESYSTEM_IS_CPPFS
            ${CMAKE_CURRENT_BINARY_DIR}/_filesystem/cppfs
            ${CMAKE_CURRENT_LIST_DIR}/filesystem.cpp
            LINK_LIBRARIES c++fs
            OUTPUT_VARIABLE FILESYSTEM_CPPFS_ERROR
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED ON
        )
        if(FILESYSTEM_IS_CPPFS)
            set(FILESYSTEM_LIBRARIES c++fs)
        else()
            message(FATAL_ERROR
                    "No support for C++ filesystem\n"
                    "${FILESYSTEM_NATIVE_ERROR}\n"
                    "${FILESYSTEM_STDCPPFS_ERROR}\n"
                    "${FILESYSTEM_CPPFS_ERROR}\n"
            )
        endif()
    endif()
endif()

