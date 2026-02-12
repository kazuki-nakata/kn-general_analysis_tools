! ===== stats_utils_1d.f90 =====
module stats_utils_1d
  implicit none
  private
  public :: calculate_mad, calculate_var

contains

  !--- クイックセレクトで 0-origin k番目を返す（破壊的） ---
  subroutine kth_smallest_1d(a, n, k0, val)
    implicit none
    real(4), intent(inout) :: a(*)   ! 長さ >= n
    integer(4), intent(in) :: n, k0   ! 0 <= k < n
    real(4), intent(out)   :: val
    integer(4) :: left, right, i, j, k, mid
    real(4) :: pivot, tmp

    k=k0
    if (n <= 0) then
      val = 0.0; return
    end if

    left = 1; right = n
    do
      if (right - left <= 16) then
        do i = left+1, right
          tmp = a(i); j = i-1
          do while (j >= left .and. a(j) > tmp)
            a(j+1) = a(j); j = j-1
          end do
          a(j+1) = tmp
        end do
        val = a(left + k)
        return
      end if

      mid   = (left + right) / 2
      pivot = a(mid)
      i = left; j = right
      do
        do while (a(i) < pivot); i = i + 1; end do
        do while (a(j) > pivot); j = j - 1; end do
        if (i <= j) then
          tmp = a(i); a(i) = a(j); a(j) = tmp
          i = i + 1; j = j - 1
        end if
        if (i > j) exit
      end do

      if (left + k <= j) then
        right = j
      else if (left + k >= i) then
        k = k - (i - left)
        left = i
      else
        val = a(left + k); return
      end if
    end do
  end subroutine kth_smallest_1d

  !--- 1D中央値 ---
  subroutine median_1d(v, n, med)
    real(4), intent(in)  :: v(*)
    integer(4), intent(in) :: n
    real(4), intent(out) :: med
    real(4), allocatable :: w(:)
    real(4) :: a, b

    if (n <= 0) then
      med = 0.0; return
    end if

    allocate(w(n)); w = v(1:n)
    if (mod(n,2) == 1) then
      call kth_smallest_1d(w, n, n/2, med)
    else
      call kth_smallest_1d(w, n, n/2-1, a)
      call kth_smallest_1d(w, n, n/2,   b)
      med = 0.5*(a + b)
    end if
    deallocate(w)
  end subroutine median_1d

  !--- 1D MAD（未スケール）: median(|x - median(x)|) ---
  subroutine calculate_mad(v, n, mad_out)
    real(4), intent(in)  :: v(*)
    integer(4), intent(in) :: n
    real(4), intent(out) :: mad_out
    real(4), allocatable :: w(:), d(:)
    real(4) :: med
    integer(4) :: i

    if (n <= 0) then
      mad_out = 0.0; return
    end if

    allocate(w(n)); w = v(1:n)
    call median_1d(w, n, med)

    allocate(d(n))
    do i = 1, n
      d(i) = abs(v(i) - med)
    end do
    call median_1d(d, n, mad_out)

    deallocate(w); deallocate(d)
  end subroutine calculate_mad

  !--- 1D 分散（unbiased: /(n-1), biased: /(n)）---
  subroutine calculate_var(v, n, var_out, unbiased)
    real(4), intent(in)  :: v(*)
    integer(4), intent(in) :: n
    real(4), intent(out) :: var_out
    logical,  intent(in), optional :: unbiased
    logical :: ub
    integer(4) :: i
    real(8) :: mean8, s2

    if (n <= 1) then
      var_out = 0.0; return
    end if
    ub = .true.; if (present(unbiased)) ub = unbiased

    mean8 = 0.0d0
    do i = 1, n
      mean8 = mean8 + dble(v(i))
    end do
    mean8 = mean8 / dble(n)

    s2 = 0.0d0
    do i = 1, n
      s2 = s2 + (dble(v(i)) - mean8)**2
    end do

    if (ub) then
      var_out = real(s2 / dble(n-1), kind=4)
    else
      var_out = real(s2 / dble(n),   kind=4)
    end if
  end subroutine calculate_var



end module stats_utils_1d
