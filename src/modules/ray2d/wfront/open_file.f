c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

      subroutine open_file( f_timeout, error, fname_timeout )

      implicit none

      integer f_timeout, error

      error = 0
      f_timeout = 0
      fname_timeout = "time.out"
      open( f_timeout, file=fname_timeout, err=900 )

      goto 999
 900  write(*,*) "Could not open file ", fname_timeout
      error = 1

 999  end

