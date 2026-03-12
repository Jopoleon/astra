module machine_config

use json_module, only : json_file
use io_mod, only: machine
  
implicit none

type(json_file) :: config
character(len=120) :: json_cfg
logical :: cfg_exists

contains

!---------------------------------------------------------------------
    subroutine config_read

    json_cfg = 'exp/cnf/' // trim(machine) // '_description_in.json'
    call config%initialize(compact_reals=.true.)
    call config%load(filename=json_cfg)
    INQUIRE(FILE=TRIM(json_cfg), EXIST=cfg_exists)

    end subroutine config_read

!---------------------------------------------------------------------
    subroutine config_close

    call config%destroy()

    end subroutine config_close

end module machine_config
