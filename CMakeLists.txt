cmake_minimum_required(VERSION 3.16)
project (itkHomework)

find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

add_executable (itkHomework itkHomework.cxx)

target_link_libraries (itkHomework ${ITK_LIBRARIES})
