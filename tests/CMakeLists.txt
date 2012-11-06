INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

  SET(SOURCES
  testmain.cpp
  testtransportsolver.cpp
  testlist2D.cpp
  smallfunctions.cpp
    )

  ADD_EXECUTABLE(tests ${SOURCES})
  TARGET_LINK_LIBRARIES(tests ${CPPUNIT_LIBRARY})

