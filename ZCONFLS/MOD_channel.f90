!===========================================================================
    Module channel
!===========================================================================
!   define scattering channels for one partial wave
!---------------------------------------------------------------------------
    Implicit none

    Integer :: lpar =  0  ! total L for given patial wave
    Integer :: ispar=  0  ! (2*S+1) for given patial wave
    Integer :: ipar =  0  ! parity 
    Integer :: jpar =  0  ! (2*J+1)
    Integer :: nch  =  0  ! number of one-electron channels
    Integer :: mch  =  0  ! max. number of channels
    Integer :: imch =256  ! initial prediction of mch
    Integer :: jmch =256  ! increment for mch   

    Integer, allocatable :: iptar(:)  ! target pointer
    Integer, allocatable :: lch(:)    ! small l for given channel
    Integer, allocatable :: ipch(:)   ! place in the common list of orbitals
    Integer, allocatable :: ipconf(:) ! pointer to the last configuration
    Integer, allocatable :: jkch(:)   ! k-number in case of jK-coupling (2K+1)

    Character(4), allocatable :: ELC(:) ! spectroscopic symbol
	
    Integer :: ncp	       ! number of configurations in perturber
    Integer :: nwp        ! number of orbitals in perturber
    Character(20) :: AFP  ! file-name for perturber

    Integer :: npert=0    ! number of configurations in perturber
    Integer :: mpert=0
    Integer :: ipert=64
    Integer, allocatable :: ippert(:)  ! pointer to the last configuration
                       
    Character(3) :: Tpar

    End Module channel


!=======================================================================    
    Subroutine Allocate_channel(m)
!=======================================================================    
!   allocate arrays in module "channel"
!-----------------------------------------------------------------------    
    Use channel

    Implicit none
    Integer, intent(in) :: m
    Integer, allocatable :: ia(:)
    Character(4), allocatable :: aa(:)

    if(m.le.0) then
     if(Allocated(iptar)) Deallocate(iptar,lch,ipch,ELC,ipconf,jkch)
     mch = 0; nch=0
    elseif(.not.allocated(iptar)) then
     mch = m
     Allocate(iptar(mch),ipconf(mch),lch(mch),ipch(mch), &
              ELC(mch),jkch(mch))
    elseif(m.le.mch) then
     Return
    elseif(nch.gt.0) then
     mch = m
     Allocate(ia(nch))
     ia=iptar(1:nch); Deallocate(iptar)
     Allocate(iptar(mch)); iptar(1:nch)=ia
     ia=ipconf(1:nch); Deallocate(ipconf)
     Allocate(ipconf(mch)); ipconf(1:nch)=ia
     ia=lch(1:nch); Deallocate(lch)
     Allocate(lch(mch)); lch(1:nch)=ia
     ia=ipch(1:nch); Deallocate(ipch)
     Allocate(ipch(mch)); ipch(1:nch)=ia
     ia=jkch(1:nch); Deallocate(jkch)
     Allocate(jkch(mch)); jkch(1:nch)=ia
     Deallocate(ia)
     Allocate(aa(nch))
     aa=ELC(1:nch); Deallocate(ELC)
     Allocate(ELC(mch)); ELC(1:nch)=aa
     Deallocate(aa)
    end if   
    
    End Subroutine Allocate_channel
    

!=======================================================================    
    Subroutine Allocate_pert(m)
!=======================================================================    
!   allocate "pertuber" arrays in module "channel"
!-----------------------------------------------------------------------    
    Use channel

    Implicit none
    Integer, intent(in) :: m
    Integer, allocatable :: ia(:)

    if(m.le.0) then
     if(Allocated(ippert)) Deallocate(ippert)
     mpert = 0; npert=0
    elseif(.not.allocated(ippert)) then
     mpert = m
     Allocate(ippert(0:mpert)); ippert=0
    elseif(m.le.mpert) then
     Return
    else
     mpert = m
     Allocate(ia(0:npert))
     ia=ippert(0:npert); Deallocate(ippert)
     Allocate(ippert(0:mpert)); ippert(0:npert)=ia(0:npert)
     Deallocate(ia)
    end if   
    
    End Subroutine Allocate_pert


!======================================================================
      Subroutine R_channel(nut,klsp)
!======================================================================
!     read from file 'nut' information for channel klsp 
!----------------------------------------------------------------------
      Use channel
      
      Implicit none
      Integer, intent(in) :: nut,klsp
      Character(20) :: AF
      Character(80) :: line
      Integer :: i,j, nlsp, ncfg, nu, Ifind_position

      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'R_channel: nlsp <= 0 '
      if(klsp.gt.nlsp) Stop 'R_channel: klsp > nlsp '
      
      Do i = 1,klsp
       read(nut,'(a)') line
       read(line,*) AF,lpar,ispar,ipar
       AFP=' '; ncp=0; nwp=0 
       j = INDEX(line,'pert')
       if(j.gt.0) read(line(j:),*) AFP,ncp,nwp 
       jpar = -1; if(ispar.eq.0) jpar = lpar + 1
      End do

      i = Ifind_position(nut,'channels:')
      read(nut,*)
      read(nut,*)

      Call Allocate_channel(0)
      Do
       read(nut,'(a)') line; 
       if(line(4:4).ne.'.') Cycle
       read(line(1:3),*) i;  if(i.ne.klsp) Cycle
       read(line,*) AF,AF,AF,AF,nch,AF,AF,ncfg,ncp
       i = nch; Call Allocate_channel(i)
       Do i = 1,nch
        read(nut,*) ELC(i),lch(i),iptar(i),j,ipconf(i),jkch(i)
       End do
       Exit
      End do

      Call Allocate_pert(0)
      if(ncp.le.0) Return

      Call Find_free_unit(nu)
      i = LEN_TRIM(AFP); AF=AFP(1:i)//'.c'
      open(nu,file=AF)
      Call R_pert(nu)
      Close(nu)

      End Subroutine R_channel


!======================================================================
      Subroutine R_pert(nu)
!======================================================================
!     read from unit 'nu' information about perturbers 
!----------------------------------------------------------------------
      Use channel
      
      Implicit none
      Integer, Intent(in) :: nu
      Integer :: i
 
      Call Allocate_pert(0)
      Call Read_ipar(nu,'npert',npert)
      if(npert.gt.0) then      
       Call Allocate_pert(ipert+npert)
       read(nu,*) ippert(1:npert)
      elseif(ncp.gt.0) then
       npert=ncp
       Call Allocate_pert(ipert+npert)
       Do i=0,npert; ippert(i)=i; End do
      end if

      End Subroutine R_pert
