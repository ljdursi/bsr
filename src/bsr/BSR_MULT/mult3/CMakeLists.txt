add_executable(mult3
  mod_param.f90
  mod_spinorb.f90
  mod_termexp.f90
  mod_coef.f90
  mod_zoef.f90
  add_res.f90
  conf_loop.f90
  def_ic.f90
  det_mult.f90
  idet_fact.f90
  mult3.f90
  prepare.f90
  rw.f90
  r_conf.f90
  term_loop.f90
  inf_mult.f90)
target_link_libraries(mult3 zconfls shells zcom)
