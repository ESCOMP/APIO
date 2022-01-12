module aio_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for APIO (asyncronous PIO interface)
  !----------------------------------------------------------------------------

  use ESMF             , only : ESMF_Mesh, ESMF_GridComp, ESMF_Time, ESMF_TimeInterval
  use ESMF             , only : ESMF_State, ESMF_Clock, ESMF_SUCCESS, ESMF_LOGMSG_INFO
  use ESMF             , only : ESMF_TraceRegionEnter, ESMF_TraceRegionExit
  use ESMF             , only : ESMF_Alarm, ESMF_METHOD_INITIALIZE, ESMF_MethodAdd, ESMF_MethodRemove
  use ESMF             , only : ESMF_TimeGet, ESMF_ClockGet, ESMF_GridCompSetEntryPoint
  use ESMF             , only : ESMF_ClockGetAlarm, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff
  use ESMF             , only : operator(+), ESMF_LogWrite
  use NUOPC            , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model      , only : model_routine_SS        => SetServices
  use NUOPC_Model      , only : model_label_Advance     => label_Advance
  use NUOPC_Model      , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model      , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model      , only : NUOPC_ModelGet, SetVM
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod    , only : SHR_CONST_SPVAL
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_ymd2date
  use shr_mpi_mod      , only : shr_mpi_bcast
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_state_diagnose, chkerr, memcheck
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_advance, shr_strdata_get_stream_domain
  use dshr_strdata_mod , only : shr_strdata_init_from_config
  use dshr_mod         , only : dshr_model_initphase, dshr_init
  use dshr_mod         , only : dshr_state_setscalar, dshr_set_runclock, dshr_check_restart_alarm
  use dshr_mod         , only : dshr_restart_read, dshr_restart_write, dshr_mesh_init
  use dshr_dfield_mod  , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize

  implicit none
  private ! except

  public  :: SetServices
  public  :: SetVM
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: drof_comp_run
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: model_mesh                          ! model mesh
  character(len=CS)            :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: mpicom                              ! mpi communicator
  integer                      :: my_task                             ! my task in mpi communicator mpicom
  logical                      :: masterproc                          ! true of my_task == master_task
  character(len=16)            :: inst_suffix = ""                    ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit                             ! logging unit number
  logical                      :: restart_read
  character(CL)                :: case_name                           ! case name
  character(*) , parameter     :: nullstr = 'null'
                                                                      ! drof_in namelist input
  character(CL)                :: streamfilename = nullstr            ! filename to obtain stream info from
  character(CL)                :: nlfilename = nullstr                ! filename to obtain namelist info from
  character(CL)                :: dataMode = nullstr                  ! flags physics options wrt input data
  character(CL)                :: model_meshfile = nullstr            ! full pathname to model meshfile
  character(CL)                :: model_maskfile = nullstr            ! full pathname to obtain mask from
  character(CL)                :: restfilm = nullstr                  ! model restart file namelist
  integer                      :: nx_global
  integer                      :: ny_global

  logical                      :: diagnose_data = .true.
  integer      , parameter     :: master_task=0                       ! task number of master task
  character(*) , parameter     :: rpfile = 'rpointer.rof'
#ifdef CESMCOUPLED
  character(*) , parameter     :: modName =  "(aio_comp_nuopc)"
#else
  character(*) , parameter     :: modName =  "(cdeps_drof_comp)"
#endif

  ! linked lists
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! model mask and model fraction
  real(r8), pointer            :: model_frac(:) => null()
  integer , pointer            :: model_mask(:) => null()

  ! module pointer arrays
  real(r8), pointer            :: Forr_rofl(:) => null()
  real(r8), pointer            :: Forr_rofi(:) => null()

  character(*) , parameter     :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local varaibles
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=dshr_model_initphase, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, specRoutine=dshr_set_runclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer           :: inst_index ! number of current instance (ie. 1)
    integer           :: nu         ! unit number
    integer           :: ierr       ! error code
    logical           :: exists     ! check for file existence
    type(fldlist_type), pointer :: fldList
    character(len=*),parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(*)    ,parameter :: F00 = "('(" // trim(modName) // ") ',8a)"
    character(*)    ,parameter :: F01 = "('(" // trim(modName) // ") ',a,2x,i8)"
    character(*)    ,parameter :: F02 = "('(" // trim(modName) // ") ',a,l6)"
    !-------------------------------------------------------------------------------

    namelist / drof_nml / datamode, model_meshfile, model_maskfile, &
         restfilm, nx_global, ny_global

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Obtain flds_scalar values, mpi values, multi-instance values and
    ! set logunit and set shr logging to my log file
    call dshr_init(gcomp, 'ROF', sdat, mpicom, my_task, inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
         logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine logical masterproc
    masterproc = (my_task == master_task)

    ! Read drof_nml from nlfilename
    if (masterproc) then
       nlfilename = "drof_in"//trim(inst_suffix)
       open (newunit=nu,file=trim(nlfilename),status="old",action="read")
       read (nu,nml=drof_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(nlfilename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(nlfilename))
       end if

       ! write namelist input to standard out
       write(logunit,F00)' datamode = ',trim(datamode)
       write(logunit,F00)' model_meshfile = ',trim(model_meshfile)
       write(logunit,F00)' model_maskfile = ',trim(model_maskfile)
       write(logunit,F01)' nx_global = ',nx_global
       write(logunit,F01)' ny_global = ',ny_global
       write(logunit,F00)' restfilm = ',trim(restfilm)
    end if

    ! broadcast namelist input
    call shr_mpi_bcast(datamode                  , mpicom, 'datamode')
    call shr_mpi_bcast(model_meshfile            , mpicom, 'model_meshfile')
    call shr_mpi_bcast(model_maskfile            , mpicom, 'model_maskfile')
    call shr_mpi_bcast(nx_global                 , mpicom, 'nx_global')
    call shr_mpi_bcast(ny_global                 , mpicom, 'ny_global')
    call shr_mpi_bcast(restfilm                  , mpicom, 'restfilm')

    ! Validate datamode
    if (trim(datamode) == 'copyall') then
       if (masterproc) write(logunit,*) 'drof datamode = ',trim(datamode)
    else
       call shr_sys_abort(' ERROR illegal drof datamode = '//trim(datamode))
    end if

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldlist_add(fldsExport, "Forr_rofl")
    call dshr_fldlist_add(fldsExport, "Forr_rofi")

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(drof_comp_advertise): Fr_rof '//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine InitializeAdvertise

  !===============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_TIME) :: currTime
    integer         :: current_ymd  ! model date
    integer         :: current_year ! model year
    integer         :: current_mon  ! model month
    integer         :: current_day  ! model day
    integer         :: current_tod  ! model sec into model date
    character(len=*), parameter :: F00   = "('" // trim(modName) // ": ')',8a)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

   end subroutine InitializeRealize

   subroutine ModelFinalize(gcomp, rc)
     type(ESMF_GridComp)  :: gcomp
     integer, intent(out) :: rc
     rc = ESMF_SUCCESS
     if (masterproc) then
        write(logunit,*)
        write(logunit,*) 'drof : end of main integration loop'
        write(logunit,*)
     end if
   end subroutine ModelFinalize
   
end module aio_comp_nuopc
