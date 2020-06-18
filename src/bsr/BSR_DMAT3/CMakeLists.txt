add_executable(bsr_dmat3
  mod_param.f90            mod_cmdata.f90           mod_coeffs.f90           
  mod_dets.f90             mod_dmat.f90             mod_ndef.f90             
  mod_ndet.f90             add_cdata.f90            add_coef.f90             
  bsr_dmat.f90             check_det.f90            dd_out.f90               
  det_fact.f90             dvec_out.f90             d_data.f90               
  d_matr.f90               d_out.f90                d_prep.f90               
  fl_kpol.f90              gen_matrix.f90           gen_zf.f90               
  get_dm.f90               get_dvl.f90              get_dvr.f90              
  idef_type.f90            j_case.f90               ls_case.f90              
  merge_cdata.f90          read_arg.f90             read_bsw.f90             
  read_sol.f90             read_conf.f90            vdet.f90                 
  ifind_pert.f90           inf_bsr_dmat.f90)
target_link_libraries(bsr_dmat3 bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})
