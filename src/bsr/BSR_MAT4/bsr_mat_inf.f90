!======================================================================
      Subroutine bsr_mat_inf
!======================================================================
!     provide information about bsr_mat program
!----------------------------------------------------------------------
      Character :: A
      Integer :: nu = 99;  Character(80) :: AF = 'bsr_mat_inf'

      iarg = command_argument_count();  if(iarg.eq.0) Return
      Call GET_COMMAND_ARGUMENT(1,A);   if(A.ne.'?') Return
      open(nu,file=AF) 

      write(nu,'(a)') &
'                                                                            ',&
'BSR_MAT - sets up the interaction matrix in B-spline representation         ',&
'                                                                            ',&
'INPUT FILES:   (prepared by bsr_prep, bsr_conf, bsr_breit programs)         ',&
'                                                                            ',&
'knot.dat       -  B-spline grid                                             ',&
'bsr_par        -  input parameters                                          ',&
'target         -  target states and channel information                     ',&
'cfg.nnn        -  configuration list for partial wave  nnn                  ',&
'int_bnk.nnn    -  angular coefficient data bank                             ',&
'target.bsw     -  target w.f. in B-spline basis                             ',&
'target_orb     -  substitution orbitals information                         ',&
'pert_nnn.bsw   -  perturb w.f., if any                                      ',&
'                                                                            ',&
'OUTPUT FILES:                                                               ',&
'                                                                            ',&
'bsr_mat.nnn    -  resulting interaction matrix                              ',&
'mat_log.nnn    -  running information                                       ',&
'int_mat.nnn    -  debug output of integrals (optional)                      ',&
'                                                                            ',&
'MAIN PARAMETERS (with default values)                                       ',&
'                                                                            ',&
'klsp1 = 1  -  first partial wave under consideration                        ',&
'klsp2 = 1  -  last partial wave under consideration                         ',&
'klsp  = 1  -  condider this partial wave only                               ',&
'                                                                            ',&
'iitar = 0  -  target states are supposed to be orthogonomal eigenstates     ',&
'              of target Hamiltonian                                         ',&
'              =1 - target state should be orthogonal, but may not be        ',&
'                   eigenstates                                              ',&
'              =2 - target states may be non-orthpgonal (general case)       ',&
'                                                                            ',&
'mrel =  0  -  relativistic corrections:                                     ',&
'                                                                            ',&
'              mrel=1 - include only one-electron scalar rel. integrals      ',&
'              mrel=2 - plus one-electron spin-orbit interaction             ',&
'              mrel=3 - plus two-electron spin-other-orbit interaction       ',&
'              mrel=4 - plus spin-spin interaction                           ',&
'              mrel=5 - plus orbit-orbit interaction                         ',&
'                                                                            ',&
'              each rel. correction can also be controled by parameters      ',&
'              mso, msoo, mss, moo with values -1,0,+1, meaning              ',&
'              +1 - include anyway, -1 - exculed anyway, 0 - follow mrel     ',&
'                                                                            ',&
'imvc = -1  -  mode for processing the mass-velocity term                    ',&
'              = +1 - include mass-velocity term directly in the Hamiltonian ',&
'              =  0 - exclude anyway                                         ',&
'              = -1 - include mass-velocity term only later,                 ',&
'                     in the BSR_HD program as a first-order correction      ',&
'                                                                            ',&
'izcorr= 0  - if =1, small-radius cut-off will be applied to the spin-orbut  ',&
'             interaction. Recomended for Z > 40.                            ',&
'                                                                            ',&
'zcorr= 1.0 - semiempirical correction to spin-orbit parameters with l=1     ',&
'                                                                            ',&
'EC = ...   - core energy (if needed to be adjusted)                         ',&
'                                                                            ',&
'mk = 7     - maximum multipole index                                        ',&
'                                                                            ',&
'nb = 2000  - number of blocks in module c_data                              ',&
'kb = 5000  - max. number of blocks for given type of integrals              ',&
'mb =  100  - size of block, with total memory needed for c_data             ',&
'             28*mb*nb/(1024d0*1024d0) Mb                                    ',&
'                                                                            ',&
'maxnc = 1000000 - size of buffer for angular coefficients,                  ',&
'                  with memory needed 20*maxnc/(1024d0*1024d0) Mb            ',&
'                                                                            ',&
'eps_ovl = 1d-8  - tolerance for one-elecytron overlaps                      ',&
'eps_c   = 1d-10 - tolerance for integral coefficients                       ',&
'eps_det = 1d-10 - tolerance for total overlap factors                       ',&
'                                                                            ',&
's_ovl   = 0.75  - channel overlap limit for additional orth. conditions     ',&                                          
's_pert  = 0.5   - channel-pertuber limit                                    ',&                                        
'                                                                            ',&
'eps_tar = 1d-6  - tolerance for target energies and overlaps                ',&
'eps_acf = 1D-5  - tolerance for asymptotic coefficients                     ',&
' '

      write(*,'(a)') &
'                                                                            ',&
'BSR_MAT - sets up the interaction matrix in B-spline representation         ',&
'                                                                            ',&
'INPUT FILES are prepared by bsr_prep, bsr_conf programs and by              ',&
'bsr_breit program for angular coefficients (int_bnk.nnn)                    ',&
'                                                                            ',&
'for selected partial waves in range from n1 to n2, call program as:         ',&
'                                                                            ',&
'               bsr_mat  klsp1=n1  klsp2=n2                                  ',&
'                                                                            ',&
'for more detailed description of optional input arguments,                  ',&
'see file "bsr_mat_inf", created after this call                             ',&
'(also check bsr_mat.log and mat_log.nnn files for details and warnings)     ',&                                                              
'                                                                            '
      Call Stop_mpi (0,0,'Read the bsr_mat_inf')                                                               
                                                                            
      End Subroutine bsr_mat_inf                                          
                                                                            








