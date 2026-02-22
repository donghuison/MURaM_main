module muram_eos_cpu
   use iso_fortran_env, only: int32, real32, real64
   implicit none
   private

   public :: eos_init_cpu, eos_finalize_cpu
   public :: t_interp, p_interp, s_interp, ne_interp, rhoi_interp, amb_interp
   public :: d3_interp, eps3_interp
   public :: n_eps, n_lr, n_lp3, n_s3, eps_off, ss_off

   integer, parameter :: ik = int32
   integer, parameter :: rk = real64

   logical :: eos_ready = .false.

   integer(ik) :: n_eps = 0_ik
   integer(ik) :: n_lr = 0_ik
   integer(ik) :: n_lp3 = 0_ik
   integer(ik) :: n_s3 = 0_ik

   real(rk) :: eps_off = 0.0_rk
   real(rk) :: ss_off = 0.0_rk

   real(rk) :: del_eps = 0.0_rk
   real(rk) :: del_lr = 0.0_rk
   real(rk) :: del_lp3 = 0.0_rk
   real(rk) :: del_s3 = 0.0_rk

   real(rk) :: inv_del_eps = 0.0_rk
   real(rk) :: inv_del_lr = 0.0_rk
   real(rk) :: inv_del_lp3 = 0.0_rk
   real(rk) :: inv_del_s3 = 0.0_rk

   real(rk), allocatable :: xeps(:), xlr(:), xlp3(:), xs3(:)
   real(rk), allocatable :: p_eostab(:,:), t_eostab(:,:), s_eostab(:,:)
   real(rk), allocatable :: ne_eostab(:,:), rhoi_eostab(:,:), amb_eostab(:,:)
   real(rk), allocatable :: eps3_eostab(:,:), d3_eostab(:,:)

contains

   subroutine eos_init_cpu(eos_file)
      character(len=*), intent(in) :: eos_file

      integer :: unit_id
      integer :: ios
      integer(ik) :: i, j
      integer(ik) :: n_fwd, n_inv
      integer(ik) :: ind1, ind2, ind3, ind4, ind5, ind6
      real(real32) :: header32(14)
      real(real32), allocatable :: buf32(:)
      real(rk) :: eps0, eps1, lr0, lr1, lp0, lp1, s0, s1
      character(len=512) :: iomsg

      if (eos_ready) call eos_finalize_cpu()

      open(newunit=unit_id, file=trim(eos_file), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios, iomsg=iomsg)
      if (ios /= 0) then
         write(*,*) "ERROR: eos_init_cpu failed to open EOS file:", trim(eos_file)
         write(*,*) "I/O message:", trim(iomsg)
         error stop 1
      end if

      read(unit_id) header32

      n_eps = int(header32(1), kind=ik)
      n_lr = int(header32(2), kind=ik)
      n_lp3 = int(header32(3), kind=ik)
      n_s3 = int(header32(4), kind=ik)

      eps0 = real(header32(5), rk)
      eps1 = real(header32(6), rk)
      lr0 = real(header32(7), rk)
      lr1 = real(header32(8), rk)
      lp0 = real(header32(9), rk)
      lp1 = real(header32(10), rk)
      s0 = real(header32(11), rk)
      s1 = real(header32(12), rk)

      eps_off = real(header32(13), rk)
      ss_off = real(header32(14), rk)

      del_eps = (eps1 - eps0) / real(n_eps - 1_ik, rk)
      del_lr = (lr1 - lr0) / real(n_lr - 1_ik, rk)
      del_lp3 = (lp1 - lp0) / real(n_lp3 - 1_ik, rk)
      del_s3 = (s1 - s0) / real(n_s3 - 1_ik, rk)

      inv_del_eps = 1.0_rk / del_eps
      inv_del_lr = 1.0_rk / del_lr
      inv_del_lp3 = 1.0_rk / del_lp3
      inv_del_s3 = 1.0_rk / del_s3

      allocate(xeps(n_eps), xlr(n_lr), xlp3(n_lp3), xs3(n_s3))
      allocate(p_eostab(n_eps, n_lr), t_eostab(n_eps, n_lr), s_eostab(n_eps, n_lr))
      allocate(ne_eostab(n_eps, n_lr), rhoi_eostab(n_eps, n_lr), amb_eostab(n_eps, n_lr))
      allocate(eps3_eostab(n_lp3, n_s3), d3_eostab(n_lp3, n_s3))

      do i = 1_ik, n_eps
         xeps(i) = eps0 + real(i - 1_ik, rk) * del_eps
      end do
      do i = 1_ik, n_lr
         xlr(i) = lr0 + real(i - 1_ik, rk) * del_lr
      end do
      do i = 1_ik, n_lp3
         xlp3(i) = lp0 + real(i - 1_ik, rk) * del_lp3
      end do
      do i = 1_ik, n_s3
         xs3(i) = s0 + real(i - 1_ik, rk) * del_s3
      end do

      n_fwd = n_eps * n_lr
      allocate(buf32(6 * n_fwd))
      read(unit_id) buf32

      do i = 1_ik, n_eps
         do j = 1_ik, n_lr
            ind1 = (j - 1_ik) * n_eps + i
            ind2 = ind1 + n_fwd
            ind3 = ind2 + n_fwd
            ind4 = ind3 + n_fwd
            ind5 = ind4 + n_fwd
            ind6 = ind5 + n_fwd
            p_eostab(i, j) = real(buf32(ind1), rk)
            t_eostab(i, j) = real(buf32(ind2), rk)
            s_eostab(i, j) = real(buf32(ind3), rk)
            ne_eostab(i, j) = real(buf32(ind4), rk)
            rhoi_eostab(i, j) = real(buf32(ind5), rk)
            amb_eostab(i, j) = real(buf32(ind6), rk)
         end do
      end do
      deallocate(buf32)

      n_inv = n_lp3 * n_s3
      allocate(buf32(2 * n_inv))
      read(unit_id) buf32

      do i = 1_ik, n_lp3
         do j = 1_ik, n_s3
            ind1 = (j - 1_ik) * n_lp3 + i
            ind2 = ind1 + n_inv
            eps3_eostab(i, j) = real(buf32(ind1), rk)
            d3_eostab(i, j) = real(buf32(ind2), rk)
         end do
      end do
      deallocate(buf32)

      close(unit_id)
      eos_ready = .true.
   end subroutine eos_init_cpu

   subroutine eos_finalize_cpu()
      if (allocated(xeps)) deallocate(xeps)
      if (allocated(xlr)) deallocate(xlr)
      if (allocated(xlp3)) deallocate(xlp3)
      if (allocated(xs3)) deallocate(xs3)

      if (allocated(p_eostab)) deallocate(p_eostab)
      if (allocated(t_eostab)) deallocate(t_eostab)
      if (allocated(s_eostab)) deallocate(s_eostab)
      if (allocated(ne_eostab)) deallocate(ne_eostab)
      if (allocated(rhoi_eostab)) deallocate(rhoi_eostab)
      if (allocated(amb_eostab)) deallocate(amb_eostab)

      if (allocated(eps3_eostab)) deallocate(eps3_eostab)
      if (allocated(d3_eostab)) deallocate(d3_eostab)

      eos_ready = .false.
   end subroutine eos_finalize_cpu

   function t_interp(ee, dd) result(temp)
      real(rk), intent(in) :: ee, dd
      real(rk) :: temp
      integer(ik) :: i_e, i_d
      real(rk) :: f_e, f_d, loge, logr, logt

      call ensure_ready()
      loge = log(ee + eps_off)
      logr = log(dd)
      call uniform_cell(loge, xeps(1), inv_del_eps, n_eps, i_e, f_e)
      call uniform_cell(logr, xlr(1), inv_del_lr, n_lr, i_d, f_d)
      logt = bilinear_value(t_eostab, i_e, f_e, i_d, f_d)
      temp = exp(logt)
   end function t_interp

   function p_interp(ee, dd) result(pres)
      real(rk), intent(in) :: ee, dd
      real(rk) :: pres
      integer(ik) :: i_e, i_d
      real(rk) :: f_e, f_d, loge, logr, logp

      call ensure_ready()
      loge = log(ee + eps_off)
      logr = log(dd)
      call uniform_cell(loge, xeps(1), inv_del_eps, n_eps, i_e, f_e)
      call uniform_cell(logr, xlr(1), inv_del_lr, n_lr, i_d, f_d)
      logp = bilinear_value(p_eostab, i_e, f_e, i_d, f_d)
      pres = exp(logp)
   end function p_interp

   function s_interp(ee, dd) result(ss)
      real(rk), intent(in) :: ee, dd
      real(rk) :: ss
      integer(ik) :: i_e, i_d
      real(rk) :: f_e, f_d, loge, logr

      call ensure_ready()
      loge = log(ee + eps_off)
      logr = log(dd)
      call uniform_cell(loge, xeps(1), inv_del_eps, n_eps, i_e, f_e)
      call uniform_cell(logr, xlr(1), inv_del_lr, n_lr, i_d, f_d)
      ss = bilinear_value(s_eostab, i_e, f_e, i_d, f_d) - ss_off
   end function s_interp

   function ne_interp(ee, dd) result(ne)
      real(rk), intent(in) :: ee, dd
      real(rk) :: ne
      integer(ik) :: i_e, i_d
      real(rk) :: f_e, f_d, loge, logr

      call ensure_ready()
      loge = log(ee + eps_off)
      logr = log(dd)
      call uniform_cell(loge, xeps(1), inv_del_eps, n_eps, i_e, f_e)
      call uniform_cell(logr, xlr(1), inv_del_lr, n_lr, i_d, f_d)
      ne = exp(bilinear_value(ne_eostab, i_e, f_e, i_d, f_d))
   end function ne_interp

   function rhoi_interp(ee, dd) result(rhoi)
      real(rk), intent(in) :: ee, dd
      real(rk) :: rhoi
      integer(ik) :: i_e, i_d
      real(rk) :: f_e, f_d, loge, logr

      call ensure_ready()
      loge = log(ee + eps_off)
      logr = log(dd)
      call uniform_cell(loge, xeps(1), inv_del_eps, n_eps, i_e, f_e)
      call uniform_cell(logr, xlr(1), inv_del_lr, n_lr, i_d, f_d)
      rhoi = exp(bilinear_value(rhoi_eostab, i_e, f_e, i_d, f_d))
   end function rhoi_interp

   function amb_interp(ee, dd) result(amb)
      real(rk), intent(in) :: ee, dd
      real(rk) :: amb
      integer(ik) :: i_e, i_d
      real(rk) :: f_e, f_d, loge, logr

      call ensure_ready()
      loge = log(ee + eps_off)
      logr = log(dd)
      call uniform_cell(loge, xeps(1), inv_del_eps, n_eps, i_e, f_e)
      call uniform_cell(logr, xlr(1), inv_del_lr, n_lr, i_d, f_d)
      amb = exp(bilinear_value(amb_eostab, i_e, f_e, i_d, f_d))
   end function amb_interp

   function d3_interp(pp, ss) result(dd)
      real(rk), intent(in) :: pp, ss
      real(rk) :: dd
      integer(ik) :: i_p, i_s
      real(rk) :: f_p, f_s, logp, ss1

      call ensure_ready()
      logp = log(pp)
      ss1 = ss + ss_off
      call uniform_cell(logp, xlp3(1), inv_del_lp3, n_lp3, i_p, f_p)
      call uniform_cell(ss1, xs3(1), inv_del_s3, n_s3, i_s, f_s)
      dd = exp(bilinear_value(d3_eostab, i_p, f_p, i_s, f_s))
   end function d3_interp

   function eps3_interp(pp, ss) result(ee)
      real(rk), intent(in) :: pp, ss
      real(rk) :: ee
      integer(ik) :: i_p, i_s
      real(rk) :: f_p, f_s, logp, ss1

      call ensure_ready()
      logp = log(pp)
      ss1 = ss + ss_off
      call uniform_cell(logp, xlp3(1), inv_del_lp3, n_lp3, i_p, f_p)
      call uniform_cell(ss1, xs3(1), inv_del_s3, n_s3, i_s, f_s)
      ee = exp(bilinear_value(eps3_eostab, i_p, f_p, i_s, f_s)) - eps_off
   end function eps3_interp

   subroutine ensure_ready()
      if (.not. eos_ready) then
         error stop "EOS not initialized. Call eos_init_cpu first."
      end if
   end subroutine ensure_ready

   pure subroutine uniform_cell(x, x0, inv_dx, n, idx, frac)
      real(rk), intent(in) :: x, x0, inv_dx
      integer(ik), intent(in) :: n
      integer(ik), intent(out) :: idx
      real(rk), intent(out) :: frac

      real(rk) :: u

      u = (x - x0) * inv_dx
      idx = int(floor(u), kind=ik) + 1_ik
      if (idx < 1_ik) idx = 1_ik
      if (idx > n - 1_ik) idx = n - 1_ik

      frac = u - real(idx - 1_ik, rk)
      if (frac < 0.0_rk) frac = 0.0_rk
      if (frac > 1.0_rk) frac = 1.0_rk
   end subroutine uniform_cell

   pure function bilinear_value(tbl, ix, fx, iy, fy) result(val)
      real(rk), intent(in) :: tbl(:,:)
      integer(ik), intent(in) :: ix, iy
      real(rk), intent(in) :: fx, fy
      real(rk) :: val
      real(rk) :: w00, w10, w01, w11

      w00 = (1.0_rk - fx) * (1.0_rk - fy)
      w10 = fx * (1.0_rk - fy)
      w01 = (1.0_rk - fx) * fy
      w11 = fx * fy

      val = w00 * tbl(ix, iy) + w10 * tbl(ix + 1_ik, iy) + &
         w01 * tbl(ix, iy + 1_ik) + w11 * tbl(ix + 1_ik, iy + 1_ik)
   end function bilinear_value

end module muram_eos_cpu
