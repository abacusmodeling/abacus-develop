AddTest(
  TARGET md_LJ_pot
  SOURCES LJ_pot_test.cpp
)

AddTest(
  TARGET md_func
  SOURCES MD_func_test.cpp
)

AddTest(
  TARGET md_nve
  SOURCES NVE_test.cpp
)

AddTest(
  TARGET md_nvt_ads
  SOURCES NVT_ADS_test.cpp
)

AddTest(
  TARGET md_nvt_nhc
  SOURCES NVT_NHC_test.cpp
)

AddTest(
  TARGET md_msst
  SOURCES MSST_test.cpp
)

AddTest(
  TARGET md_fire
  SOURCES FIRE_test.cpp
)

AddTest(
  TARGET md_lgv
  SOURCES Langevin_test.cpp
)