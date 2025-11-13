program test_init_driver
  use iso_c_binding
  implicit none
  interface
    subroutine c_input() bind(C, name="input")
    end subroutine c_input
    subroutine c_initialize() bind(C, name="initialize")
    end subroutine c_initialize
  end interface

  call c_input()
  write(*,'(I0,4X)',advance='no') 1
  print '(A)', 'Input Finished. Beginning initialisation.'
  call c_initialize()
  print '(A)', 'Initialisation complete.'
end program test_init_driver


