PROJECT(dymaxion)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11) #2.8.11 for TARGET_INCLUDE_DIRECTORIES

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()


IF(DEFINED CMAKE_INSTALL_PREFIX)
   SET(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "Choose the install directory.")
ELSE()
   SET(CMAKE_INSTALL_PREFIX /opt/dymaxion/ CACHE STRING "Choose the install directory.")
ENDIF()


FIND_PACKAGE (ITK REQUIRED)
IF ( ITK_FOUND )
      INCLUDE( ${ITK_USE_FILE} ) # https://stackoverflow.com/questions/33583101/cmake-3-4-0-including-multiple-files#33585398
ENDIF( ITK_FOUND )


SET(Libraries
  ${ITK_LIBRARIES}
)


SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags") # use RPATH instead of RUNPATH (currently ignored on Debian-9), -WL tells gcc to pass option to the linker


SET(CurrentExe "equirec2dymaxion")
ADD_EXECUTABLE(${CurrentExe}
  ${CurrentExe}.cxx graypr.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)
