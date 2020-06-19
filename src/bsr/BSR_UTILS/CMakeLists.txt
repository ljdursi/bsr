# extract seperate states from "bound.nnn" file:
add_executable(bound_bsw bound_bsw.f90)
target_link_libraries(bound_bsw bs zconfls zcom ${BLAS_LAPACK_LIBRARIES})

# list of states in the "bound.nnn" files:
add_executable(bound_tab bound_tab.f90)
target_link_libraries(bound_tab zcom)

# hydrogen-like bound and pseudocontinuum states for given B-soline basis:
add_executable(bs_hwf bs_hwf.f90)
target_link_libraries(bs_hwf bs zcom ${BLAS_LAPACK_LIBRARIES})

# calculation of dipole matrix elements between R-matrx states saved in rsol.nnn files
add_executable(bsr_dd bsr_dd.f90)
target_link_libraries(bsr_dd zconfls zcom)

# merging  1.bsw + 2.bsw + 3.bsw + ... --> res.bsw
add_executable(bsw123 bsw123.f90)
target_link_libraries(bsw123 bs zcom ${BLAS_LAPACK_LIBRARIES})

# connverts the bsw-files with orbitals in B-spline repersentation 
# into tab-files suitable for grafic display
add_executable(bsw_tab bsw_tab.f90)
target_link_libraries(bsw_tab bs zcom ${BLAS_LAPACK_LIBRARIES})

# connvert the B-spline results in the MCHF format
add_executable(bsw_w bsw_w.f90)
target_link_libraries(bsw_w bs zcom ${BLAS_LAPACK_LIBRARIES})

# target, rsol.nnn  --> rsol.all
add_executable(collect_rsol collect_rsol.f90)
target_link_libraries(collect_rsol zcom)

# # target, threshoulds -> target.exp
# add_executable(convert_target convert_target.f90)
# target_link_libraries(convert_target zconfls zcom)

# # H.DAT  -->  f-values
# add_executable(f_values f_values.f90)
# target_link_libraries(f_values zconfls zcom)

# different photoionization paramters
add_executable(photo_collect photo_collect.f90)
target_link_libraries(photo_collect zcom)

# # different photoionization paramters
# add_executable(photo_tab_new photo_tab_new.f90)
# target_link_libraries(photo_tab_new zconfls zcom)

# prints B-splines in separate files:
add_executable(print_bs print_bs.f90)
target_link_libraries(print_bs bs zcom ${BLAS_LAPACK_LIBRARIES})

# # H.DAT --> dipole matrix elements
# add_executable(s_values s_values.f90)
# target_link_libraries(s_values zconfls zcom)

# target states information
add_executable(target_states target_states.f90)
target_link_libraries(target_states zconfls zcom)

# Generate B-spline representation for orbitals 
# given in MCHF package format, w-files 
add_executable(w_bsw w_bsw.f90)
target_link_libraries(w_bsw bs zcom ${BLAS_LAPACK_LIBRARIES})

#  f-value calculations between states in bound.nnn files
add_executable(zf_cc_bsr zf_cc_bsr.f90)
target_link_libraries(zf_cc_bsr zconfls zcom)

#  f-value calculations between set of states in the c- and bsw-files
add_executable(zf_bb_bsr zf_bb_bsr.f90)
target_link_libraries(zf_bb_bsr zconfls zcom)

add_executable(target_kpert target_kpert.f90)
target_link_libraries(target_kpert zconfls zcom)
