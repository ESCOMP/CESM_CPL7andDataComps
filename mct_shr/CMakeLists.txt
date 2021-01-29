# Build a separate list containing the mct wrapper and its dependencies. That
# way, this list can be easily included in unit test builds that link to mct,
# but excluded from builds that do not include mct.
list(APPEND share_mct_sources
  mct_mod.F90
  shr_mct_mod.F90
  shr_pcdf_mod.F90)


sourcelist_to_parent(share_mct_sources)
