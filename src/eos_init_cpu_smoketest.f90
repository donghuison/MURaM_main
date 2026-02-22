program eos_init_cpu_smoketest
   use iso_fortran_env, only: real64
   use muram_eos_cpu, only: eos_init_cpu, eos_finalize_cpu, n_eps, n_lr, n_lp3, n_s3, eps_off, ss_off, &
      t_interp, p_interp, s_interp, ne_interp, rhoi_interp, amb_interp, d3_interp, eps3_interp
   implicit none

   character(len=512) :: eos_file
   integer :: argc
   logical :: found
   real(real64) :: ee, dd, tt, pp, ss, ne, rhoi, amb
   real(real64) :: dd_back, ee_back

   argc = command_argument_count()
   if (argc >= 1) then
      call get_command_argument(1, eos_file)
      inquire(file=trim(eos_file), exist=found)
      if (.not. found) then
         print *, "ERROR: EOS file not found:", trim(eos_file)
         print *, "Hint: run with an absolute path or from project root."
         error stop 1
      end if
   else
      call resolve_default_eos_path(eos_file, found)
      if (.not. found) then
         print *, "ERROR: Could not find default EOS file."
         print *, "Tried:"
         print *, "  RUNDIR/Uppsala_mergedeos_PI_A.dat"
         print *, "  ../RUNDIR/Uppsala_mergedeos_PI_A.dat"
         print *, "  ../../RUNDIR/Uppsala_mergedeos_PI_A.dat"
         print *, "Please pass EOS path as first argument."
         error stop 1
      end if
   end if

   print *, "=== eos_init_cpu_smoketest ==="
   print *, "EOS file:", trim(eos_file)

   call eos_init_cpu(trim(eos_file))

   print *, "eos_init_cpu: initialized OK"
   print *, "n_eps n_lr n_lp3 n_s3 =", n_eps, n_lr, n_lp3, n_s3
   print *, "eps_off ss_off =", eps_off, ss_off

   ee = 1.0d13
   dd = 1.0d-7

   tt = t_interp(ee, dd)
   pp = p_interp(ee, dd)
   ss = s_interp(ee, dd)
   ne = ne_interp(ee, dd)
   rhoi = rhoi_interp(ee, dd)
   amb = amb_interp(ee, dd)

   dd_back = d3_interp(pp, ss)
   ee_back = eps3_interp(pp, ss)

   print *, "sample forward EOS at ee=1e13, rho=1e-7:"
   print *, "T =", tt
   print *, "P =", pp
   print *, "S =", ss
   print *, "ne =", ne
   print *, "rhoi =", rhoi
   print *, "amb =", amb

   print *, "sample inverse EOS from (P,S):"
   print *, "rho_back =", dd_back
   print *, "eps_back =", ee_back

   call eos_finalize_cpu()
   print *, "eos_finalize_cpu: done"

contains

   subroutine resolve_default_eos_path(path, found_path)
      character(len=*), intent(out) :: path
      logical, intent(out) :: found_path

      character(len=512) :: c1, c2, c3

      c1 = "RUNDIR/Uppsala_mergedeos_PI_A.dat"
      c2 = "../RUNDIR/Uppsala_mergedeos_PI_A.dat"
      c3 = "../../RUNDIR/Uppsala_mergedeos_PI_A.dat"

      inquire(file=trim(c1), exist=found_path)
      if (found_path) then
         path = c1
         return
      end if

      inquire(file=trim(c2), exist=found_path)
      if (found_path) then
         path = c2
         return
      end if

      inquire(file=trim(c3), exist=found_path)
      if (found_path) then
         path = c3
         return
      end if

      path = c1
   end subroutine resolve_default_eos_path

end program eos_init_cpu_smoketest
