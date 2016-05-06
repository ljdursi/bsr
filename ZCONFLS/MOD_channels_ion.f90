!===========================================================================
      MODULE channels_ion
!===========================================================================
!     define scattering channels 
!---------------------------------------------------------------------------
    
      Implicit none

      Integer :: nlsp =  0                !  number of partial waves               

      Integer, allocatable :: ispar(:)    !  (2*S+1) for partial waves
      Integer, allocatable :: lpar(:)     !  total L 
      Integer, allocatable :: ipar(:)     !  parity (+1 | -1)
      Integer, allocatable :: jpar(:)     !  (2*J+1)
      Integer, allocatable :: nch(:)      !  number of one-electron channels
    
      Integer, allocatable :: iptar(:,:)  !  pointer on the target state
      Integer, allocatable :: lch(:,:)    !  small l for given channel
      Integer, allocatable :: ipch(:,:)   !  pointer on the place in the common list of orbitals

      Integer, allocatable :: ipconf(:,:) !  pointer on the last configuration
      Integer, allocatable :: jkch(:,:)   !  k-number in jK-coupling

      CHARACTER(3), allocatable :: Tpar(:)  ! spectroscopic notation for partial wave
      CHARACTER(4), allocatable :: ELC(:,:) ! spectroscopic symbol for channel orbital

      Integer :: mch  =  0 !  maximum number of one-electron channels

      CHARACTER(20), allocatable :: AFP(:) ! file-name for perturber

      Integer, allocatable :: ncp(:)  !  number of configurations in perturber	
      Integer, allocatable :: nwp(:)  !  number of orbitals in perturber

      END MODULE channels_ion


!=======================================================================    
    Subroutine Allocate_channels_ion
!=======================================================================    
        
    Use channels_ion

    if(nlsp.le.0) Stop ' Allocate_channels: nlsp <= 0 '
    if(mch.le.0) Stop ' Allocate_channels: mch <= 0 '

    Allocate(ispar(nlsp),lpar(nlsp),ipar(nlsp),jpar(nlsp),Tpar(nlsp), &
             nch(nlsp),iptar(nlsp,mch),lch(nlsp,mch),ipch(nlsp,mch),  &
             ELC(nlsp,mch), ipconf(nlsp,mch), jkch(nlsp,mch),         &
             AFP(nlsp),ncp(nlsp),nwp(nlsp))
    lch = 0

    End subroutine Allocate_channels_ion





!======================================================================
      Subroutine R_channels_ion(nut)
!======================================================================
!     read from file 'nut' channels information 
!----------------------------------------------------------------------

      USE channels_ion
      
      Implicit none
      
      Integer, Intent(in) :: nut
      Character(20) :: AF
      Character(80) :: line
      Integer :: i,j,m,Icheck_file, Ifind_position
 
      Call Read_ipar(nut,'max_ch',mch)
      if(mch.le.0) Stop 'R_channel: mch <= 0 '

      Call Read_ipar(nut,'nlsp',nlsp); read(nut,*)
      if(nlsp.le.0) Stop 'R_channel: nlsp <= 0 '

      Call Allocate_channels_ion      


      Do i = 1,nlsp
       read(nut,'(a)') line
       read(line,*) AF,lpar(i),ispar(i),ipar(i)
       AFP(i)=' '; ncp(i)=0; nwp(i)=0 
       j = INDEX(line,'pert')
       if(j.gt.0) read(line(j:),*) AF,ncp(i),nwp(i) 
       j=LEN_TRIM(AF); if(AF(j-1:j).eq.'.c') AF(j-1:j)='  '
       AFP(i)=AF
       jpar(i) = 0; if(ispar(i).eq.0) jpar(i) = lpar(i) + 1
      End do

      i = Ifind_position(nut,'channels:'); read(nut,*)

      Do i = 1,nlsp
       read(nut,*) line
       read(nut,*) AF,AF,AF,AF,nch(i),AF,AF,m,ncp(i)
       Do j = 1,nch(i)
        read(nut,*) ELC(i,j),lch(i,j),iptar(i,j),m,ipconf(i,j),jkch(i,j)
       End do
      End do

      End Subroutine R_channels_ion


!====================================================================
      Integer Function iipar_ion (IL,IS,IP)
!====================================================================
!     partial wave index for give term
!--------------------------------------------------------------------

      USE channels_ion

      Implicit  none

      Integer, intent(in) :: IL,IS,IP
      Integer :: i

      iipar_ion = 0
      Do i = 1,nlsp
       if(lpar(i).ne.IL) Cycle
       if(ispar(i).ne.IS) Cycle
       if(ipar(i).ne.IP) Cycle
       iipar_ion = i
       Exit
      End do

      End Function iipar_ion

