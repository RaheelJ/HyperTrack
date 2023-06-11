#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "GeographicLib::GeographicLib_STATIC" for configuration "Debug"
set_property(TARGET GeographicLib::GeographicLib_STATIC APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::GeographicLib_STATIC PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/Geographic_d.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::GeographicLib_STATIC )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::GeographicLib_STATIC "${_IMPORT_PREFIX}/lib/Geographic_d.lib" )

# Import target "GeographicLib::CartConvert" for configuration "Debug"
set_property(TARGET GeographicLib::CartConvert APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::CartConvert PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/CartConvert_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::CartConvert )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::CartConvert "${_IMPORT_PREFIX}/bin/CartConvert_d.exe" )

# Import target "GeographicLib::ConicProj" for configuration "Debug"
set_property(TARGET GeographicLib::ConicProj APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::ConicProj PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/ConicProj_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::ConicProj )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::ConicProj "${_IMPORT_PREFIX}/bin/ConicProj_d.exe" )

# Import target "GeographicLib::GeodesicProj" for configuration "Debug"
set_property(TARGET GeographicLib::GeodesicProj APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::GeodesicProj PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/GeodesicProj_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::GeodesicProj )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::GeodesicProj "${_IMPORT_PREFIX}/bin/GeodesicProj_d.exe" )

# Import target "GeographicLib::GeoConvert" for configuration "Debug"
set_property(TARGET GeographicLib::GeoConvert APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::GeoConvert PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/GeoConvert_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::GeoConvert )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::GeoConvert "${_IMPORT_PREFIX}/bin/GeoConvert_d.exe" )

# Import target "GeographicLib::GeodSolve" for configuration "Debug"
set_property(TARGET GeographicLib::GeodSolve APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::GeodSolve PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/GeodSolve_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::GeodSolve )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::GeodSolve "${_IMPORT_PREFIX}/bin/GeodSolve_d.exe" )

# Import target "GeographicLib::GeoidEval" for configuration "Debug"
set_property(TARGET GeographicLib::GeoidEval APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::GeoidEval PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/GeoidEval_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::GeoidEval )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::GeoidEval "${_IMPORT_PREFIX}/bin/GeoidEval_d.exe" )

# Import target "GeographicLib::Gravity" for configuration "Debug"
set_property(TARGET GeographicLib::Gravity APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::Gravity PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/Gravity_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::Gravity )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::Gravity "${_IMPORT_PREFIX}/bin/Gravity_d.exe" )

# Import target "GeographicLib::MagneticField" for configuration "Debug"
set_property(TARGET GeographicLib::MagneticField APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::MagneticField PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/MagneticField_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::MagneticField )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::MagneticField "${_IMPORT_PREFIX}/bin/MagneticField_d.exe" )

# Import target "GeographicLib::Planimeter" for configuration "Debug"
set_property(TARGET GeographicLib::Planimeter APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::Planimeter PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/Planimeter_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::Planimeter )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::Planimeter "${_IMPORT_PREFIX}/bin/Planimeter_d.exe" )

# Import target "GeographicLib::RhumbSolve" for configuration "Debug"
set_property(TARGET GeographicLib::RhumbSolve APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::RhumbSolve PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/RhumbSolve_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::RhumbSolve )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::RhumbSolve "${_IMPORT_PREFIX}/bin/RhumbSolve_d.exe" )

# Import target "GeographicLib::TransverseMercatorProj" for configuration "Debug"
set_property(TARGET GeographicLib::TransverseMercatorProj APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(GeographicLib::TransverseMercatorProj PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/TransverseMercatorProj_d.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS GeographicLib::TransverseMercatorProj )
list(APPEND _IMPORT_CHECK_FILES_FOR_GeographicLib::TransverseMercatorProj "${_IMPORT_PREFIX}/bin/TransverseMercatorProj_d.exe" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
