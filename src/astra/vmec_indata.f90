module vmec_indata
!-----------------------------------------------------------------------
! Write VMEC &INDATA namelist for VMEC / LIBSTELL.
!-----------------------------------------------------------------------

implicit none
private
public :: write_vmec_indata, vmec_profile_index, nvmec_prof_max

integer, parameter :: nvmec_prof_max = 100
! spline_cubic / Akima both need at least 4 knots (profile_functions.f
! aborts with 'check s-grid for ...' below that).
integer, parameter :: nvmec_prof_min = 4

! Keys whose value this module supplies; any template line assigning one
! of them is dropped. Matched case-insensitively, blanks squeezed out.
character(len=*), parameter :: OWNED_KEYS = &
    '|PMASS_TYPE|PMASS_FILE|PIOTA_TYPE|PIOTA_FILE|PCURR_TYPE|PCURR_FILE' // &
    '|AM|AI|AC|AM_AUX_S|AM_AUX_F|AI_AUX_S|AI_AUX_F|AC_AUX_S|AC_AUX_F' // &
    '|CURTOR|PHIEDGE|'
character(len=*), parameter :: AXIS_KEYS = &
    '|RAXIS|RAXIS_CC|RAXIS_CS|ZAXIS|ZAXIS_CC|ZAXIS_CS|'

contains

!-----------------------------------------------------------------------
! Pick np <= nwant indices out of ASTRA's 1..na1 radial grid, always
! keeping the first and the last. Uniform stride in *index* space, so the
! radial clustering of the ASTRA grid is carried over
!-----------------------------------------------------------------------
    subroutine vmec_profile_index(na1, nwant, idx, np)
        integer, intent(in)  :: na1, nwant
        integer, intent(out) :: idx(:)
        integer, intent(out) :: np
        integer :: i

! -2 leaves room for the s=0 and s=1 knots write_vmec_indata may have to
! add: ASTRA's grid starts at a half cell, and VMEC evaluates the splines
! on the closed interval [0,1].
        np = min(na1, nwant, nvmec_prof_max - 2, size(idx))
        if (np < nvmec_prof_min) then
            write(*, *) 'vmec_indata: need at least ', nvmec_prof_min, &
                        ' profile points, got ', np
            error stop 'vmec_indata: profile grid too coarse for VMEC splines'
        endif

        do i = 1, np
            idx(i) = 1 + nint(dble(i-1)*dble(na1-1)/dble(np-1))
        enddo
        idx(1)  = 1
        idx(np) = na1
    end subroutine vmec_profile_index

!-----------------------------------------------------------------------
! template  : &INDATA file to start from (boundary, resolution, switches)
! outfile   : namelist actually handed to xvmec2000
! s         : normalised toroidal flux, strictly increasing, s(np) = 1
! pres      : pressure [Pa]                    -> AM_AUX_F
! iota      : rotational transform             -> AI_AUX_F
! curi      : *enclosed* toroidal current [A]  -> AC_AUX_F
! curtor    : total toroidal current [A]
! phiedge   : edge toroidal flux [Wb]
! prof_type : spline flavour, 'cubic_spline' or 'akima_spline'
! nequil    : replaces the $NEQUIL placeholder (VMEC surface count)
! raxis/zaxis/naxis : optional magnetic-axis guess -> RAXIS_CC / ZAXIS_CS
!-----------------------------------------------------------------------
    subroutine write_vmec_indata(template, outfile, s, pres, iota, curi, np, &
                                 curtor, phiedge, prof_type, nequil, &
                                 raxis, zaxis, naxis)
        character(len=*), intent(in) :: template, outfile
        integer,          intent(in) :: np
        double precision, intent(in) :: s(np), pres(np), iota(np), curi(np)
        double precision, intent(in) :: curtor, phiedge
        character(len=*), intent(in) :: prof_type
        integer,          intent(in) :: nequil
        double precision, intent(in), optional :: raxis(:), zaxis(:)
        integer,          intent(in), optional :: naxis

        double precision, parameter :: s_tol = 1.d-8
        double precision :: sv(np+2), pv(np+2), qv(np+2), cv(np+2)
        integer :: uin, uout, ios, i, nax, nv
        logical :: have_axis, saw_slash, drop
        character(len=512) :: line
        character(len=64)  :: key

        if (np < nvmec_prof_min .or. np > nvmec_prof_max - 2) then
            write(*, *) 'vmec_indata: np = ', np, ' outside [', nvmec_prof_min, &
                        ',', nvmec_prof_max - 2, ']'
            error stop 'vmec_indata: bad number of profile points'
        endif
        do i = 2, np
            if (s(i) <= s(i-1)) then
                write(*, *) 'vmec_indata: s not increasing at i = ', i, s(i-1), s(i)
                error stop 'vmec_indata: profile s-grid must be strictly increasing'
            endif
        enddo
        if (s(1) < 0.d0 .or. s(np) > 1.d0 + s_tol) then
            write(*, *) 'vmec_indata: s(1) = ', s(1), ' s(np) = ', s(np)
            error stop 'vmec_indata: profile s-grid outside [0,1]'
        endif

! VMEC's spline_cubic returns iflag=-1 for any abscissa outside [s(1),s(n)]
! and pmass/piota/pcurr are evaluated over the full closed [0,1], so both
! ends must be present exactly. ASTRA's radial grid starts at a half cell
! (XRHO(1) > 0), so the magnetic axis is normally missing and gets
! prepended here by linear extrapolation
        nv = 0
        if (s(1) > s_tol) then
            nv = 1
            sv(1) = 0.d0
            pv(1) = extrap_axis(s(1), s(2), pres(1), pres(2))
            qv(1) = extrap_axis(s(1), s(2), iota(1), iota(2))
            cv(1) = 0.d0        ! no current is enclosed by the magnetic axis
        endif
        do i = 1, np
            sv(nv+i) = s(i)
            pv(nv+i) = pres(i)
            qv(nv+i) = iota(i)
            cv(nv+i) = curi(i)
        enddo
        nv = nv + np
        sv(1) = 0.d0            ! snap off the roundoff at both ends
        if (sv(nv) < 1.d0 - s_tol) then
            nv = nv + 1
            sv(nv) = 1.d0
            pv(nv) = extrap_edge(sv(nv-2), sv(nv-1), pv(nv-2), pv(nv-1))
            qv(nv) = extrap_edge(sv(nv-2), sv(nv-1), qv(nv-2), qv(nv-1))
            cv(nv) = extrap_edge(sv(nv-2), sv(nv-1), cv(nv-2), cv(nv-1))
        else
            sv(nv) = 1.d0
        endif

        have_axis = present(raxis) .and. present(zaxis) .and. present(naxis)
        nax = 0
        if (have_axis) nax = naxis
        have_axis = have_axis .and. nax > 0

        open(newunit=uin, file=trim(template), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*, *) 'vmec_indata: cannot open template ', trim(template)
            error stop 'vmec_indata: missing VMEC template'
        endif
        open(newunit=uout, file=trim(outfile), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*, *) 'vmec_indata: cannot write ', trim(outfile)
            error stop 'vmec_indata: cannot write VMEC input'
        endif

        saw_slash = .false.
        drop      = .false.
        key       = ''
        do
            read(uin, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (is_namelist_end(line)) then
                saw_slash = .true.
                exit                    ! re-emitted after the appended block
            endif
            call subst_nequil(line, nequil)
            if (is_passthrough(line)) then
                write(uout, '(A)') trim(line)
                cycle
            endif
            if (index(line, '=') > 1) then
                key = lhs_key(line)
                drop = index(OWNED_KEYS, '|'//trim(key)//'|') > 0 .or. &
                       (have_axis .and. index(AXIS_KEYS, '|'//trim(key)//'|') > 0)
            else
                if (len_trim(key) == 0) then
                    write(*, *) 'vmec_indata: stray value before any key: ', trim(line)
                    error stop 'vmec_indata: malformed VMEC template'
                endif
            endif
            if (.not. drop) write(uout, '(A)') trim(line)
        enddo
        close(uin)

        if (.not. saw_slash) then
            write(*, *) 'vmec_indata: no terminating "/" in ', trim(template)
            error stop 'vmec_indata: malformed VMEC template'
        endif

        write(uout, '(A)') ''
        write(uout, '(A)') '! ---- profiles written by ASTRA (src/astra/vmec_indata.f90) ----'
        write(uout, '(A)') " PMASS_TYPE = '"//trim(prof_type)//"'"
        write(uout, '(A)') " PIOTA_TYPE = '"//trim(prof_type)//"'"
        write(uout, '(A)') " PCURR_TYPE = '"//trim(prof_type)//"_I'"
        write(uout, '(A,ES23.15)') ' CURTOR  = ', curtor
        write(uout, '(A,ES23.15)') ' PHIEDGE = ', phiedge
        call write_array(uout, 'AM_AUX_S', sv, nv)
        call write_array(uout, 'AM_AUX_F', pv, nv)
        call write_array(uout, 'AI_AUX_S', sv, nv)
        call write_array(uout, 'AI_AUX_F', qv, nv)
        call write_array(uout, 'AC_AUX_S', sv, nv)
        call write_array(uout, 'AC_AUX_F', cv, nv)
        if (have_axis) then
            call write_array(uout, 'RAXIS_CC', raxis, nax)
            call write_array(uout, 'ZAXIS_CS', zaxis, nax)
        endif
        write(uout, '(A)') '/'
        close(uout)
    end subroutine write_vmec_indata

!-----------------------------------------------------------------------
! Replace every occurrence of the $NEQUIL placeholder in a template line
! (vmec_io/vmecinput_template.dat carries it in NS_ARRAY).
    subroutine subst_nequil(line, nequil)
        character(len=*), intent(inout) :: line
        integer,          intent(in)    :: nequil
        character(len=len(line)) :: tmp
        character(len=16) :: num
        integer :: ipos

        write(num, '(I0)') nequil
        do
            ipos = index(line, '$NEQUIL')
            if (ipos == 0) exit
            tmp = line(:ipos-1)//trim(num)//line(ipos+7:)
            line = tmp
        enddo
    end subroutine subst_nequil

!-----------------------------------------------------------------------
! Linear extrapolation from the two innermost knots back to s = 0.
    double precision function extrap_axis(s1, s2, y1, y2)
        double precision, intent(in) :: s1, s2, y1, y2
        extrap_axis = y1 - s1*(y2 - y1)/(s2 - s1)
    end function extrap_axis

!-----------------------------------------------------------------------
! Linear extrapolation from the two outermost knots out to s = 1.
    double precision function extrap_edge(s1, s2, y1, y2)
        double precision, intent(in) :: s1, s2, y1, y2
        extrap_edge = y2 + (1.d0 - s2)*(y2 - y1)/(s2 - s1)
    end function extrap_edge

!-----------------------------------------------------------------------
    subroutine write_array(u, name, a, n)
        integer,          intent(in) :: u, n
        character(len=*), intent(in) :: name
        double precision, intent(in) :: a(:)
        integer :: i, j

        write(u, '(A)') ' '//trim(name)//' ='
        do i = 1, n, 4
            write(u, '(4(1X,ES23.15))') (a(j), j = i, min(i+3, n))
        enddo
    end subroutine write_array

!-----------------------------------------------------------------------
! Blank lines, comments and the '&INDATA' header are copied verbatim.
    logical function is_passthrough(line)
        character(len=*), intent(in) :: line
        character(len=1) :: c

        is_passthrough = .true.
        if (len_trim(line) == 0) return
        c = first_char(line)
        if (c == '!' .or. c == '&' .or. c == '$') return
        is_passthrough = .false.
    end function is_passthrough

!-----------------------------------------------------------------------
    logical function is_namelist_end(line)
        character(len=*), intent(in) :: line
        is_namelist_end = .false.
        if (len_trim(line) == 0) return
        is_namelist_end = first_char(line) == '/'
    end function is_namelist_end

!-----------------------------------------------------------------------
    character(len=1) function first_char(line)
        character(len=*), intent(in) :: line
        character(len=len(line)) :: tmp

        tmp = adjustl(line)
        first_char = tmp(1:1)
    end function first_char

!-----------------------------------------------------------------------
    character(len=64) function lhs_key(line)
        use char_manip, only: clean_string, to_upper

        character(len=*), intent(in) :: line

        integer :: ieq
        character(len=256) :: line_clean

        line_clean = clean_string(line)
        lhs_key = ''
        ieq = index(line_clean, '=')
        if (ieq <= 1) return
        lhs_key = to_upper(line_clean(1: ieq-1))
    end function lhs_key

end module vmec_indata
