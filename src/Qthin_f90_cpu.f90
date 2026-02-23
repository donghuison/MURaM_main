module qthin_chianti_cpu_mod
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  real(dp), parameter :: X_H_default    = 0.7_dp
  real(dp), parameter :: ne_ref_default = 1.0e9_dp

  logical, save :: initialized = .false.
  integer, save :: ntab = 0
  real(dp), allocatable, save :: logT_tab(:), Q_tab(:)
  real(dp), save :: lgT0, del0

contains

  integer function free_unit()
    integer :: u
    logical :: opened
    do u = 10, 999
      inquire(unit=u, opened=opened)
      if (.not. opened) then
        free_unit = u
        return
      end if
    end do
    stop 'free_unit: no free unit available'
  end function free_unit


  subroutine init_radloss_chianti(filename, ne_ref)
    character(len=*), intent(in) :: filename
    real(dp), intent(in), optional :: ne_ref

    integer :: unit, i, ios
    real(dp) :: Tval, Qval, ne_ref_use

    if (present(ne_ref)) then
      ne_ref_use = ne_ref
    else
      ne_ref_use = ne_ref_default
    end if

    if (allocated(logT_tab)) deallocate(logT_tab)
    if (allocated(Q_tab))    deallocate(Q_tab)

    unit = free_unit()
    open(unit=unit, file=filename, status='old', form='formatted', iostat=ios)
    if (ios /= 0) stop 'init_radloss_chianti: cannot open Radloss_Chianti.dat'

    read(unit, *, iostat=ios) ntab
    if (ios /= 0 .or. ntab < 2) stop 'init_radloss_chianti: bad ntab'

    allocate(logT_tab(ntab), Q_tab(ntab))
    do i = 1, ntab
      read(unit, *, iostat=ios) Tval, Qval
      if (ios /= 0) stop 'init_radloss_chianti: bad table row'
      logT_tab(i) = log(Tval)
      Q_tab(i)    = Qval * ne_ref_use * ne_ref_use
    end do
    close(unit)

    lgT0 = logT_tab(1)
    del0 = 1.0_dp / (logT_tab(2) - logT_tab(1))

    initialized = .true.
  end subroutine init_radloss_chianti


  subroutine get_radloss_qthin_cpu(temp, rho, pres, ghosts, p_cut, qthin, radloss_file, X_H, ne_ref)
    real(dp), intent(in)  :: temp(:,:,:)
    real(dp), intent(in)  :: rho(:,:,:)
    real(dp), intent(in)  :: pres(:,:,:)
    integer, intent(in)   :: ghosts(3)
    real(dp), intent(in)  :: p_cut
    real(dp), intent(out) :: qthin(:,:,:)
    character(len=*), intent(in), optional :: radloss_file
    real(dp), intent(in), optional :: X_H, ne_ref

    integer :: i, j, k, n, n1, n2
    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: ibeg, iend, jbeg, jend, kbeg, kend
    real(dp) :: X_H_use, ne_ref_use, inv_pmax, ne_par
    real(dp) :: t1, t2, r1, r2, tmin, tmax
    real(dp) :: t_a, t_b, ff, ts, pr, pt, rr, qq, qloss
    real(dp) :: weight
    character(len=256) :: file_use

    if (present(X_H)) then
      X_H_use = X_H
    else
      X_H_use = X_H_default
    end if
    if (present(ne_ref)) then
      ne_ref_use = ne_ref
    else
      ne_ref_use = ne_ref_default
    end if

    if (.not. initialized) then
      if (present(radloss_file)) then
        file_use = radloss_file
      else
        file_use = 'Radloss_Chianti.dat'
      end if
      call init_radloss_chianti(trim(file_use), ne_ref_use)
    end if

    ilo = lbound(temp, 1); ihi = ubound(temp, 1)
    jlo = lbound(temp, 2); jhi = ubound(temp, 2)
    klo = lbound(temp, 3); khi = ubound(temp, 3)

    if (ghosts(1) < 1) stop 'get_radloss_qthin_cpu: ghosts(1) must be >= 1'
    if (any(ghosts < 0)) stop 'get_radloss_qthin_cpu: ghosts must be >= 0'

    ibeg = ilo + ghosts(1)
    iend = ihi - ghosts(1)
    jbeg = jlo + ghosts(2)
    jend = jhi - ghosts(2)
    kbeg = klo + ghosts(3)
    kend = khi - ghosts(3)

    if (ibeg > iend) stop 'get_radloss_qthin_cpu: empty i physical range'
    if (jbeg > jend) stop 'get_radloss_qthin_cpu: empty j physical range'
    if (kbeg > kend) stop 'get_radloss_qthin_cpu: empty k physical range'
    if (ibeg-1 < ilo .or. iend+1 > ihi) stop 'get_radloss_qthin_cpu: need at least 1 ghost in i'

    inv_pmax = 1.0_dp / p_cut
    ne_par = sqrt(0.5_dp * (1.0_dp + X_H_use) * X_H_use) * 6.0e23_dp / ne_ref_use

    qthin = 0.0_dp

    do k = kbeg, kend
      do j = jbeg, jend

        do i = ibeg, iend + 1
          t1 = log(temp(i-1, j, k))
          t2 = log(temp(i,   j, k))
          r1 = log(rho(i-1, j, k) * ne_par)
          r2 = log(rho(i,   j, k) * ne_par)

          tmin = min(t1, t2)
          tmax = max(t1, t2)

          n1 = int((tmin - lgT0) * del0) + 1
          n2 = int((tmax - lgT0) * del0) + 1

          if (n1 < 1) n1 = 1
          if (n2 > ntab-1) n2 = ntab-1

          do n = n1, n2
            t_a = max(tmin, logT_tab(n))
            t_b = min(tmax, logT_tab(n+1))

            if (t_b - t_a > 1.0e-6_dp) then
              ff = (t_b - t_a) / (tmax - tmin)
              ts = 0.5_dp * (t_a + t_b)
              pr = (t2 - ts) / (t2 - t1)
              pt = (logT_tab(n+1) - ts) / (logT_tab(n+1) - logT_tab(n))
            else if (t_b - t_a >= 0.0_dp) then
              ff = 1.0_dp
              ts = 0.5_dp * (t_a + t_b)
              pr = 0.5_dp
              pt = (logT_tab(n+1) - ts) / (logT_tab(n+1) - logT_tab(n))
            else
              ff = 0.0_dp
              pr = 0.5_dp
              pt = 0.5_dp
            end if

            rr = exp(pr*r1 + (1.0_dp - pr)*r2)
            qq = pt*Q_tab(n) + (1.0_dp - pt)*Q_tab(n+1)
            qloss = -rr*rr * qq * ff

            qthin(i-1, j, k) = qthin(i-1, j, k) + qloss * pr
            qthin(i,   j, k) = qthin(i,   j, k) + qloss * (1.0_dp - pr)
          end do
        end do

        do i = ibeg, iend
          weight = max(0.0_dp, 1.0_dp - (pres(i, j, k) * inv_pmax)**2)
          qthin(i, j, k) = qthin(i, j, k) * weight
        end do

      end do
    end do

  end subroutine get_radloss_qthin_cpu

end module qthin_chianti_cpu_mod

