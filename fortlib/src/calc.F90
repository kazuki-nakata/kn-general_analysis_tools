MODULE Calc
  IMPLICIT NONE
  INTEGER,PARAMETER :: null = 999
  REAL(4),PARAMETER :: Undef = 999
CONTAINS
  
subroutine estimate_sigma_2d(img, nx, ny, ix, iy, nz, stat_type, ndw, dst, stat)
  use stats_utils_1d, only: calculate_mad, calculate_var
   ! 推定位置 (ix,iy) ごとに、勾配画像 gradx/grady から局所ノイズσを推定
  ! stat_type: 1=MAD, 2=Var
  ! ndw: 窓サイズ（奇数を推奨）  dst: サンプリング間引き（未使用なら1）
  implicit none
  integer(4), intent(in)  :: nx, ny, nz, ndw, dst, stat_type
  integer(4), intent(in)  :: ix(nz), iy(nz)
  real(4),    intent(in)  :: img(nx,ny)
  real(4),    intent(out) :: stat(nz)
  real(4), allocatable :: buf(:)
  integer(4) :: i,j,k, l,nsel
  integer(4) :: i0, i1, j0, j1, nwi, nwj,cap, r
  real(4)    :: stat0, gain
   print *,"nx=",nx,"ny=",ny,"nz=",nz

  ! 外部の統計関数（あなたの実装が既にある前提）
  ! 2次元配列を受けてスカラーを返す想定です。1次元想定なら reshape してください。

  ! 安全のため初期化
  stat = 0.0
  r    = (ndw-1)/2

  do k = 1, nz
    ! --- 窓の境界（1-origin）---
    i0 = max(1, ix(k)-r);  i1 = min(nx, ix(k)+r)
    j0 = max(1, iy(k)-r);  j1 = min(ny, iy(k)+r)

    ! 最大要素数（間引き反映）を見積もってバッファ確保
    cap = ((i1-i0)/dst + 1) * ((j1-j0)/dst + 1)
    if (cap < 1) cap = 1
    allocate(buf(cap))

   nsel = 0
   do i = i0, i1, dst
      do j = j0, j1, dst
         nsel = nsel + 1
         buf(nsel) = img(i,j)
      end do
   end do

      select case (stat_type)
      case (1)   ! MADベース（標準正規前提の 0.67448975 でスケール）
        call calculate_mad( buf, nsel, stat0 )
        gain = 0.67448975
      case (2)   ! 分散ベース（不偏なら ddof=1 を内部で処理）
        call calculate_var( buf, nsel, stat0 )
        gain = 1.0
      case default
        ! 未知の指定は MAD 扱い
        call calculate_mad( buf, nsel, stat0 )
        gain = 0.67448975
      end select

      stat(k) = stat0/gain
      deallocate(buf)
  end do

end subroutine estimate_sigma_2d

subroutine odr_quad_map_2d(src, tgt, nx, ny, ix, iy, nz, ndw, dst,          &
                           sigma_src, sigma_tgt, maxit, tol,                &
                           pred, r2, verr, coef, info)
  ! 各位置 (ix,iy) の周囲 ndw×ndw 窓から、src→tgt の局所写像 y ≈ b0 + b1 x + b2 x^2 を
  ! ODR（誤差つき回帰）で推定。dst で間引きサンプリング。
  ! ノイズが弱く/非線形が弱い場合は ODR の近似（重み付き最小二乗）で高速に収束。
  ! 収束しない／データが少ない等は Deming（線形EIV）にフォールバックします。
  use Numerical_Module, only:deming_linear_from_buf, solve3x3_sym
  implicit none
  integer(4), intent(in)    :: nx, ny, nz, ndw, dst, maxit
  integer(4), intent(in)    :: ix(nz), iy(nz)
  real(4),    intent(in)    :: src(nx,ny), tgt(nx,ny)
  real(4),    intent(in)    :: sigma_src(nx,ny), sigma_tgt(nx,ny)  ! 観測ノイズσマップ
  real(4),    intent(in)    :: tol                                  ! 収束閾値（例:1e-5〜1e-3）
  real(4),    intent(out)   :: pred(nz), r2(nz), verr(nz)
  real(4),    intent(out)   :: coef(3,nz)   ! (b0,b1,b2)
  integer(4), intent(out)   :: info(nz)     ! 0:OK, >0:フォールバック/問題コード

  integer(4) :: k, i0, i1, j0, j1, i, j, nsel, cap, r, step, it
  real(4), allocatable :: x(:), y(:), sx(:), sy(:), w(:)
  real(8) :: b0, b1, b2, b0_old, b1_old, b2_old, gprime, wi
  real(8) :: Ey, Vy, ss_res, ss_tot, var_res, sig2_loc
  real(8) :: H(3,3), rhs(3), db(3)
  real(8), parameter :: EPS = 1.0d-12

  pred = 0.0; r2 = 0.0; verr = 0.0; coef = 0.0; info = 0
  r    = (ndw-1)/2
  step = max(1, dst)

  do k = 1, nz
    ! --- 窓境界 ---
    i0 = max(1, ix(k)-r);  i1 = min(nx, ix(k)+r)
    j0 = max(1, iy(k)-r);  j1 = min(ny, iy(k)+r)

    ! 最大長を見積もって一時バッファ確保
    cap = ((i1-i0)/step + 1) * ((j1-j0)/step + 1)
    if (cap < 6) then
      info(k) = 10   ! データ不足
      cycle
    end if
    allocate(x(cap), y(cap), sx(cap), sy(cap), w(cap))

    ! --- 窓→1次元化（dst間引き） ---
    nsel = 0
    do i = i0, i1, step
      do j = j0, j1, step
        nsel = nsel + 1
        x(nsel)  = src(i,j)
        y(nsel)  = tgt(i,j)
        sx(nsel) = sigma_src(i,j)
        sy(nsel) = sigma_tgt(i,j)
      end do
    end do
    if (nsel < 6) then
      info(k) = 11
      deallocate(x,y,sx,sy,w); cycle
    end if

    ! --- 初期値：Deming線形（b2=0） ---
    call deming_linear_from_buf(x, y, sx, sy, nsel, b0, b1)
    b2 = 0.0d0

    ! --- 反復：ODRの近似（重み付き二乗和最小化）
    ! 目的関数 ~ sum_i ( y_i - f(x_i; b) )^2 / ( sy_i^2 + (f'(x_i))^2 sx_i^2 )
    ! => 重み w_i = 1 / ( sy_i^2 + (b1 + 2*b2*x_i)^2 * sx_i^2 )
    do it = 1, maxit
      b0_old = b0; b1_old = b1; b2_old = b2

      ! 重み計算 & 正規方程式の蓄積（3×3）
      H = 0.0d0; rhs = 0.0d0
      do i = 1, nsel
        gprime = dble(b1) + 2.0d0*dble(b2)*dble(x(i))
        wi = 1.0d0 / ( dble(sy(i))**2 + (gprime*gprime)*(dble(sx(i))**2) + EPS )
        w(i) = real(wi,kind=4)

        ! 基底 h = [1, x, x^2]
        H(1,1) = H(1,1) + wi
        H(1,2) = H(1,2) + wi*dble(x(i))
        H(1,3) = H(1,3) + wi*dble(x(i))*dble(x(i))
        H(2,2) = H(2,2) + wi*dble(x(i))*dble(x(i))
        H(2,3) = H(2,3) + wi*dble(x(i))**3
        H(3,3) = H(3,3) + wi*dble(x(i))**4

        rhs(1) = rhs(1) + wi*dble(y(i))
        rhs(2) = rhs(2) + wi*dble(x(i))*dble(y(i))
        rhs(3) = rhs(3) + wi*(dble(x(i))**2)*dble(y(i))
      end do
      ! 対称成分を埋める
      H(2,1) = H(1,2); H(3,1) = H(1,3); H(3,2) = H(2,3)

      ! 解く：H * [b0 b1 b2]^T = rhs
      db(1) = 0.0d0; db(2) = 0.0d0; db(3) = 0.0d0
      call solve3x3_sym(H, rhs, db)   ! db を解ベクトルとして使用
      b0 = db(1); b1 = db(2); b2 = db(3)

      ! 収束判定
      if ( abs(b0-b0_old) < tol .and. abs(b1-b1_old) < tol .and. abs(b2-b2_old) < tol ) exit
    end do

    ! --- 失敗時は Deming 線形へフォールバック ---
    if (it > maxit) then
      call deming_linear_from_buf(x, y, sx, sy, nsel, b0, b1)
      b2 = 0.0d0
      info(k) = 1     ! フォールバック
    end if

    ! --- 係数・予測・指標 ---
    coef(1,k) = real(b0,kind=4)
    coef(2,k) = real(b1,kind=4)
    coef(3,k) = real(b2,kind=4)

    pred(k) = coef(1,k) + coef(2,k)*src(ix(k),iy(k)) + coef(3,k)*(src(ix(k),iy(k))**2)

    ! R^2 と残差分散（ターゲットノイズを差し引き）
    Ey = 0.0d0; Vy = 0.0d0
    do i = 1, nsel
      Ey = Ey + dble(y(i))
    end do
    Ey = Ey / dble(nsel)
    do i = 1, nsel
      Vy = Vy + (dble(y(i)) - Ey)**2
    end do
    Vy = Vy / dble(max(1,nsel-1))

    ss_res = 0.0d0
    do i = 1, nsel
      ss_res = ss_res + ( dble(y(i)) - (b0 + b1*dble(x(i)) + b2*(dble(x(i))**2)) )**2
    end do
    ss_tot = Vy * dble(max(1,nsel-1))
    if (ss_tot > 0.0d0) then
      r2(k) = real( max(0.0d0, min(1.0d0, 1.0d0 - ss_res/(ss_tot+EPS))), kind=4 )
    else
      r2(k) = 0.0
    end if

    var_res = ss_res / dble(max(1,nsel-1))

    ! 窓内のターゲットノイズ分散の平均を引く
    sig2_loc = 0.0d0
    do i = i0, i1, step
      do j = j0, j1, step
        sig2_loc = sig2_loc + dble(sigma_tgt(i,j))**2
      end do
    end do
    sig2_loc = sig2_loc / dble(nsel)

    verr(k) = real( max(0.0d0, var_res - sig2_loc), kind=4 )

    deallocate(x,y,sx,sy,w)
  end do
end subroutine odr_quad_map_2d


END MODULE Calc

