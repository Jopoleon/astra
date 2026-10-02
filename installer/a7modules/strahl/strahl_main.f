c     Standalone driver for the STRAHL library from the ASTRA 7 package.
c     ASTRA 8 (sbr/a2strahl.f90) runs "strahl a q" in the instance directory.
      program strahl_main
      implicit none
      integer iarg, narg
      character*1 arg, inp_mode
      character*70 control_file
      logical verbose, only_dens, rad_only, diag_only, rate_out,
     >     temporary, only_tot_dens, quiet, calc_rad, astracall
      control_file = 'strahl.control'
      inp_mode = 'I'
      verbose = .false.
      only_dens = .false.
      only_tot_dens = .false.
      rad_only = .false.
      diag_only = .false.
      calc_rad = .true.
      rate_out = .false.
      temporary = .false.
      quiet = .false.
      astracall = .false.
      narg = command_argument_count()
      do iarg = 1, narg
         call get_command_argument(iarg, arg)
         if (arg.eq.'a'.or.arg.eq.'A') inp_mode = 'A'
         if (arg.eq.'q'.or.arg.eq.'Q') quiet = .true.
         if (arg.eq.'n'.or.arg.eq.'N') only_dens = .true.
         if (arg.eq.'g'.or.arg.eq.'G') only_tot_dens = .true.
         if (arg.eq.'r'.or.arg.eq.'R') calc_rad = .false.
         if (arg.eq.'v'.or.arg.eq.'V') verbose = .true.
      enddo
      if (only_tot_dens) only_dens = .true.
      call strahl(verbose, only_dens, only_tot_dens, inp_mode,
     >     rad_only, calc_rad, diag_only, rate_out, control_file,
     >     quiet, temporary, astracall)
      end
c     The ASTRA 7 in-process coupling (sbr/as7_strahl_4imp.f) is not
c     used in standalone mode (astracall = .false.).
      subroutine astraout
      stop 'astraout: ASTRA 7 coupling not available'
      end
      subroutine astraset
      stop 'astraset: ASTRA 7 coupling not available'
      end
      subroutine astra_tend
      stop 'astra_tend: ASTRA 7 coupling not available'
      end
      subroutine initstrahl
      stop 'initstrahl: ASTRA 7 coupling not available'
      end
      subroutine initgrid2
      stop 'initgrid2: ASTRA 7 coupling not available'
      end
