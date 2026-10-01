subroutine A2TRAVIS

!----------------------------------------------------------------------|
! File input for: beam curvature and width, 
! tor and pol angle
!----------------------------------------------------------------------|

use status, only: NRD
use scalars, only: NA1, RTOR, BTOR, TIME, ROC, SGNIP, SGNBT, ABC, &
   TSTART, TAU
use status, only: TE, NE, FP, XRHO, ZEF, MU, ELON, SHif , IPOL, &
   AMETR, VOLUM, PEECR, CUECR, AREAT, rho_pol, FP_NORM
use read_input, only: AWD, nml_file
use numerical_tools, only: qinterp, integr
use parameters_a2equil, only : equil_now
use pi_const

implicit none

 CHARACTER(Len=256) :: inputfile_name = 'travisinput.data' ! can contain path
 CHARACTER(Len=512) :: output_dir  = 'travisOXscan.inp'
 CHARACTER(Len=256) :: inp_name   ! input file name without path
 CHARACTER(Len=256) :: inp_dir    ! path to input file
 CHARACTER(Len=256) :: CW_dir     ! path to current work directory
 CHARACTER(Len=256) :: exe_dir,exearg
 CHARACTER(Len=3)   :: branch = ' '
 CHARACTER(Len=4)   :: help = ' '
 CHARACTER(Len=5)   :: admin = ' '
 CHARACTER(Len=6)   :: silent = ' '
 CHARACTER(Len=7)   :: verbose = ' '
 CHARACTER(Len=6)   :: no_out = ' '
 CHARACTER(Len=12)  :: no_out_trace = ' '
 CHARACTER(Len=3)   :: printOpt = ' ' 
 CHARACTER(Len=256) :: TEMP_dir = ' '
 CHARACTER(Len=256) :: ipAddress = ' '
 CHARACTER(Len=256) :: astradir = ' '
 CHARACTER(Len=512) :: id_beam(500)
 INTEGER :: nexepth, ios
 INTEGER :: travis_out_inTEMP = 1  ! if 0  -> no output for storing in work directory  
                                   ! if 1  -> save all results (default)
                                   ! if 2  -> save only final results w/o trace data
!                                   ! if 3  -> save all and input data if error met (?)
 INTEGER :: travis_out_Screen = 1
character(len=120) :: fort_name, as_nml, time_str, pecr_file, phi_file, theta_file, &
     pecr_file2, phi_file2, theta_file2

integer n_r, n_beams
double precision r_a(5000),p_abs(5000),j_cd(5000) 

NAMELIST / travis / inputfile_name, output_dir, n_beams, id_beam

write(*,*) 'travis'
!stop
if (TIME-TSTART < tau) return
write(*,*) 'travis'
!stop


as_nml = TRIM(awd) // '/' // TRIM(nml_file)
write(*, *) 'Reading namelist ', TRIM(as_nml)
open(57, FILE=TRIM(as_nml), delim='apostrophe')
read(57, nml=travis, iostat=ios)
close(57)

!=======================================================================
 CALL getParameters(branch,   &
                    inputfile_name, &
                    help,     &
                    admin,    &
                    silent,   &
                    verbose,  &
                    no_out,   &
                    no_out_trace, &
                    printOpt, &
                    TEMP_dir, &
                    ipAddress)
		    
branch = 'ECH'  ! hardwired ECH

!---
     travis_out_inTEMP = 1
 IF (no_out =='no_out') THEN
     travis_out_inTEMP = 0
 ELSEIF (no_out_trace =='no_out_trace') THEN
     travis_out_inTEMP = 2
 ENDIF
!--- 
 IF (silent   =='silent')  travis_out_Screen = 0
 IF (verbose  =='verbose') travis_out_Screen = 1
 IF (printOpt =='yes') CALL exit(0) ! Display options and exit without running TRAVIS IF printOpt
!=======================================================================
 IF (branch == 'NOT') THEN  ! Exit if inputfile or branch name was not found    
   PRINT'(//a,/a)','Cannot open file: ', TRIM(inputfile_name)
   CALL exit(1)
 ELSEIF (branch == 'non') THEN 
   PRINT'(//a,/a)','TRAVIS branch not found,', &
                   'check the first line of input-file: it must contain branch-name.'
   CALL exit(1)
 ENDIF
!=======================================================================
 CALL get_exe_path(exearg,exe_dir,nexepth)
 CALL get_input_path(inputfile_name,inp_dir,inp_name)
 CALL GETCWD(astradir)                     ! get current working directory
 IF(inp_dir /= ' ') CALL CHDIR(inp_dir)  ! set current working directory
 CALL GETCWD(CW_dir)                     ! get current working directory
 print'(a,a)','         Travis executable: ', TRIM(exearg)
 print'(a,a)',' Current working directory: ', TRIM(CW_dir)
 print'(a,a)','         Travis input file: ', TRIM(inp_name)
!=======================================================================
 CALL set_ExpertMode_f77(admin)
 CALL set_TracingOutput_f77(travis_out_inTEMP,travis_out_Screen)
 CALL load_InitDataFile_f77(branch,inp_name)
!---
 CALL set_caller_name_f77(ipAddress, 'TRAVIS-'//branch)
 IF (TEMP_dir .NE. ' ') CALL set_work_dir_f77(TEMP_dir)
!---
 SELECT CASE(branch)
   CASE('OXB')
!     CALL scan_OX_f77(scanfile_name)
   CASE('BXO') 
!     CALL scan_OX_f77(scanfile_name)
   CASE('ECH') 
   ! CALL NeoConductProf_f77()
   CASE('ECE') 
   CASE('RFM') 
 END SELECT
!---
 call set_nTZ_profs_f77(NA1,XRHO(1:NA1),ne(1:na1)*1.e19,Te(1:na1),Zef(1:NA1),'tor_rho','lin')
 CALL run_Beams_f77()
 CALL save_magnetic_config_f77()
 CALL closeErrorReport_f77()
! CALL ExitFromTRAVIS_f77(0)
 CALL CHDIR(astradir)                     ! get current working directory

 call read_output_a2travis(output_dir,n_beams,n_r,r_a,P_abs,j_cd,id_beam(1:n_beams))

 call qinterp(r_a(1:n_r),P_abs(1:n_r),n_r,AMETR(1:NA1)/ABC,PEECR(1:NA1),NA1)
 call qinterp(r_a(1:n_r),j_cd(1:n_r),n_r,AMETR(1:NA1)/ABC,CUECR(1:NA1),NA1)

write(*,*) 'end travis' !,peecr(1:na1),cuecr(1:na1),n_r,r_a(1:n_r),P_abs(1:n_r),j_cd(1:n_r)
return
end subroutine a2travis

subroutine read_output_a2travis(output_dir,n_beams, n_r,r_a,P_abs,j_cd,id_beam)

implicit none

integer n_r, n_beams,j,i,iii
double precision r_a(5000),P_abs(5000),j_cd(5000),lineread(12)
 CHARACTER(Len=512) :: output_dir,id_beam(n_beams) ,dumstring

P_abs=0.
j_cd=0.
 
do j=1,n_beams
 i=0
 iii=0
 open(32,file=trim(output_dir)//'/Pabs_Icd_profiles_'//trim(id_beam(j)))
   read(32,'(A)') dumstring
!   write(*,*) dumstring
   do while (iii==0)
    i=i+1
    read(32,*) lineread
!    write(*,*) lineread
    r_a(i)=lineread(1) 
    P_abs(i)=P_abs(i)+lineread(6) 
    j_cd(i)=j_cd(i)+lineread(8)*1.e-6 
    if (r_a(i)>0.9999999999) iii=1
   enddo
 close(32)
enddo
n_r=i


end subroutine read_output_a2travis

 SUBROUTINE getBranchName(inputfile_name, branch)
!=======================================================================
!* Search the branch name for TRAVIS in the inputfile_name. 
! Only the first line is analized.  
! It must contain the unique name of the branch, one from: 
! ECRH, ECH, ECE, RFM, RfM, OXB, BXO 
!
 CHARACTER(Len=*) :: inputfile_name
 CHARACTER(Len=*) :: branch
 CHARACTER(Len=256) :: header
 INTEGER :: iunit = 11
!=======================================================================
 OPEN (iunit,FILE=TRIM(ADJUSTL(inputfile_name)),ERR=100,STATUS='OLD',ACTION='READ')
 READ (iunit,'(a)') header
 CLOSE(iunit)
 !write(*,*) TRIM(header)
!=======================================================================
! PRESELECTION ONLY FOR DEBUGGER!!!
!---
! branch = 'ECH'; return
!=======================================================================
 branch = 'noname'
 IF (INDEX(header,'ECRH') > 0 .OR. &
     INDEX(header,'ECH')  > 0 .OR. &
     INDEX(header,'ECR')  > 0) THEN 
     branch = 'ECH'
 ELSEIF (INDEX(header,'ECE') > 0) THEN 
     branch = 'ECE'
 ELSEIF (INDEX(header,'OXB') > 0) THEN 
     branch = 'OXB'
 ELSEIF (INDEX(header,'BXO') > 0) THEN 
     branch = 'BXO'
 ELSEIF (INDEX(header,'RFM') > 0 .OR. &
         INDEX(header,'RfM') > 0) THEN 
     branch = 'RFM'
 ENDIF
 RETURN
!---
100 branch = "NOT"
!=======================================================================
 END SUBROUTINE getBranchName

!#######################################################################

 SUBROUTINE getParameters(branch, inputfile_name, help, admin, &
                                  silent, verbose, no_out, no_out_trace, &
                                  printOpt, TEMP_dir, ipAddress)
!=======================================================================
!* This program returns parameters and options from TRAVIS command line: 
! ```
!   travis [options] [inputfile]   
! ```
! See description in the [[DisplayHelp(subroutine)]] above.
!=======================================================================
 CHARACTER(Len=*) :: branch
 CHARACTER(Len=*) :: inputfile_name
 CHARACTER(Len=*) :: help
 CHARACTER(Len=*) :: admin
 CHARACTER(Len=*) :: silent
 CHARACTER(Len=*) :: verbose
 CHARACTER(Len=*) :: no_out
 CHARACTER(Len=*) :: no_out_trace
 CHARACTER(Len=*) :: printOpt
 CHARACTER(Len=*) :: TEMP_dir
 CHARACTER(Len=*) :: ipAddress
 CHARACTER(Len=256) :: cmdLine = ' '
 CHARACTER(Len=256) :: param = ' '
 INTEGER :: length=0, i=0, cnt=0, stat=0 
!=======================================================================
 cnt = command_argument_count()
 DO i=cnt,1,-1
   CALL get_command_argument(i,param,length,stat)
   SELECT CASE(TRIM(param))
     CASE('-h','-help');    help    = 'help'
     CASE('-a','-admin');   admin   = 'admin'
     CASE('-s','-silent');  silent  = 'silent'
     CASE('-v','-verbose'); verbose = 'verbose'
     CASE('-no','-no_out'); no_out  = 'no_out'
     CASE('-nt','-no_out_trace'); no_out_trace = 'no_out_trace'
     CASE('-p','-print-only-options'); printOpt = 'yes'
     CASE default
        IF (INDEX(param,'-IP=') > 0 .OR. INDEX(param,'-I=') > 0) THEN
          ipAddress  = param(1+SCAN(param,'='):LEN(param))
        ELSEIF (INDEX(param,'-TEMP=') > 0 .OR. INDEX(param,'-T=') > 0) THEN
          TEMP_dir  = param(1+SCAN(param,'='):LEN(param))
        ELSE
          inputfile_name = param
        ENDIF   
   END SELECT
 ENDDO
!---
 CALL getBranchName(inputfile_name, branch)
!---
 IF (cnt == 0 .OR. help =='help') THEN
   IF (cnt == 0) print'(/a)',' You can start TRAVIS using parameters.'
!   CALL DisplayHelp()
   IF (cnt == 1) CALL exit(0)  ! exit if it was called with 'help' only
 ENDIF
!---
 IF(verbose=='verbose'.OR.printOpt =='yes') THEN
   CALL get_command(cmdLine,length,stat); 
   print'(a)',' '
   print'(a,a)',' CommandLine:', TRIM(cmdLine)
   print'(a,a)',' inputfile:', TRIM(inputfile_name)
   print'(a,a)',' branch:', branch
   print'(a,a)',' -silent:', silent
   print'(a,a)',' -verbose:', verbose
   print'(a,a)',' -no_out:', no_out
   print'(a,a)',' -no_out_trace:', no_out_trace
   print'(a,a)',' -admin:', admin
   print'(a,a)',' -print-only-options:',printOpt
   print'(a,a)',' -TEMP=', TRIM(TEMP_dir)
   print'(a,a)',' -IP=', TRIM(ipAddress)
   CALL getMagConfName(inputfile_name, param)
   IF (param.NE.' ') print'(a,a)','    ', TRIM(param)
 ENDIF
 IF (silent == 'silent') verbose=' '
 RETURN
!---
 IF (branch == 'NOT') THEN   ! Exit if inputfile or branch name was not found 
    PRINT'(//a,/a)','Cannot open file: ', TRIM(inputfile_name)
    CALL exit(1)
 ELSEIF (branch == 'non') THEN 
    PRINT'(//a,/a)','TRAVIS branch not found,', &
                    'check the first line of an input-file: it must contain correct branch-name.'
    CALL exit(1)
 ENDIF
!=======================================================================
 END SUBROUTINE getParameters
 
!#######################################################################

 SUBROUTINE getMagConfName(inputfile,mname)
!=======================================================================
 CHARACTER(LEN=*) :: inputfile,mname
 mname = ' '
 OPEN(10,file=TRIM(ADJUSTL(inputfile)),err=100,status='old',action='read')
 DO i=1,12
    READ(10,'(a)') mname
    IF(index(mname,'Magnetic_Configuration') > 0 .OR. &
       index(mname,'Magnetic_configuration') > 0) THEN
       mname=adjustl(mname); j=index(mname,' '); 
       mname='Magnetic_configuration:' // mname(j:LEN(mname))
       exit
    ENDIF
 ENDDO
 IF(index(mname,'Magnetic_').LE.0) mname = ' '
 CLOSE(10)
100 RETURN
!=======================================================================
 END SUBROUTINE getMagConfName

!#######################################################################

 SUBROUTINE get_exe_path(exearg,gpth,npth)  ! get exe path from cmd line 
    CHARACTER(len=*), INTENT(out) :: gpth,exearg 
    INTEGER, INTENT(out)          :: npth 
    INTEGER                       :: ilen,istat,npl,npw
    exearg=' ' 
    gpth=' ' 
    npth=0 
    CALL get_command_argument(0,exearg,ilen,istat) !an intrinsic function
    IF (istat == 0) THEN 
        npw  = INDEX(exearg,'\',BACK=.true.) !an intrinsic function, true means search backward
        npl  = INDEX(exearg,'/',BACK=.true.) !...
        npth = MAX(npw,npl)
        gpth = exearg(1:npth) 
    ENDIF 
 END SUBROUTINE

 !########################################################################
 
 SUBROUTINE get_input_path(inp_pathname,inp_dir,inp_name) 
    CHARACTER(len=*), INTENT(in)  :: inp_pathname 
    CHARACTER(len=*), INTENT(out) :: inp_dir,inp_name 
    INTEGER                       :: npth,npl,npw
    inp_dir=' '
    inp_name=inp_pathname
    npw  = INDEX(inp_pathname,'\',BACK=.true.)
    npl  = INDEX(inp_pathname,'/',BACK=.true.)
    npth = MAX(npw,npl)
    inp_dir  = inp_pathname(1:npth) 
    inp_name = inp_pathname(npth+1:LEN(inp_pathname)) 
  END SUBROUTINE

 !#######################################################################
  

