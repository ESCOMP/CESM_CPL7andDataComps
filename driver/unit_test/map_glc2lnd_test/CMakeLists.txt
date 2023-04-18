set (pfunit_sources
  test_map_glc2lnd.pf
  )

add_pfunit_ctest(map_glc2lnd
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES ${DRV_UNIT_TEST_LIBS})
