AddTest(
  TARGET base_matrix3
  LIBS base parallel
  SOURCES matrix3_test.cpp
)
AddTest(
  TARGET base_blas_connector
  LIBS base parallel
  SOURCES blas_connector_test.cpp
)
