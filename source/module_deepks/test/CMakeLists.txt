AddTest(
  TARGET test_deepks 
  LIBS deepks
  SOURCES main_deepks.cpp klist_1.cpp LCAO_deepks_test_prep.cpp LCAO_deepks_test.cpp nnr.cpp parallel_orbitals.cpp
)
