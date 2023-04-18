set (pfunit_sources
  test_map_lnd2rof_irrig.pf
  )

add_pfunit_ctest(map_lnd2rof_irrig
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES ${DRV_UNIT_TEST_LIBS})
