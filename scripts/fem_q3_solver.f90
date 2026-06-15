program fem_q3_solver
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nel = 80
  integer, parameter :: nn = nel + 1
  integer, parameter :: ndof = 2 * nn
  integer, parameter :: nmodes_out = 24
  integer, parameter :: nt = 2001
  integer, parameter :: nsnap = 401
  real(dp), parameter :: length_bar = 4.0_dp
  real(dp), parameter :: alpha = 0.50_dp
  real(dp), parameter :: kappa = 20.0_dp
  real(dp), parameter :: v0_bar = 1.0_dp
  real(dp), parameter :: tau_max = 16.0_dp

  integer :: i, j, nred, b_node, c_node, d_node, b_dof, c_dof, d_dof
  integer, allocatable :: map_full_to_red(:), red_to_full(:)
  real(dp), allocatable :: k_full(:, :), m_full(:, :)
  real(dp), allocatable :: k_red(:, :), m_red(:, :)
  real(dp), allocatable :: lchol(:, :), inv_l(:, :), astd(:, :)
  real(dp), allocatable :: eigvec_std(:, :), eigval(:), modes(:, :)
  real(dp), allocatable :: freq(:), xi(:), modal_amp(:), y_nodes(:, :)
  real(dp) :: h, modal_mass, sgn

  call execute_command_line('cmd /c if not exist outputs\fem mkdir outputs\fem')

  h = length_bar / real(nel, dp)
  b_node = nel / 4
  c_node = nel / 2
  d_node = 3 * nel / 4
  b_dof = 2 * b_node + 1
  c_dof = 2 * c_node + 1
  d_dof = 2 * d_node + 1

  allocate(k_full(ndof, ndof), m_full(ndof, ndof))
  k_full = 0.0_dp
  m_full = 0.0_dp

  call assemble_beam(k_full, m_full, h)
  m_full(c_dof, c_dof) = m_full(c_dof, c_dof) + alpha
  k_full(d_dof, d_dof) = k_full(d_dof, d_dof) + kappa

  allocate(map_full_to_red(ndof))
  map_full_to_red = 0
  nred = 0
  do i = 1, ndof
    if (i /= b_dof) then
      nred = nred + 1
      map_full_to_red(i) = nred
    end if
  end do
  allocate(red_to_full(nred))
  do i = 1, ndof
    if (map_full_to_red(i) > 0) red_to_full(map_full_to_red(i)) = i
  end do

  allocate(k_red(nred, nred), m_red(nred, nred))
  do i = 1, nred
    do j = 1, nred
      k_red(i, j) = k_full(red_to_full(i), red_to_full(j))
      m_red(i, j) = m_full(red_to_full(i), red_to_full(j))
    end do
  end do

  allocate(lchol(nred, nred), inv_l(nred, nred), astd(nred, nred))
  call cholesky_lower(m_red, lchol, nred)
  call invert_lower(lchol, inv_l, nred)
  astd = matmul(inv_l, matmul(k_red, transpose(inv_l)))
  astd = 0.5_dp * (astd + transpose(astd))

  allocate(eigvec_std(nred, nred), eigval(nred))
  call jacobi_symmetric(astd, eigval, eigvec_std, nred)
  call sort_eigensystem(eigval, eigvec_std, nred)

  allocate(modes(nred, nmodes_out), freq(nmodes_out), modal_amp(nmodes_out))
  do i = 1, nmodes_out
    modes(:, i) = matmul(transpose(inv_l), eigvec_std(:, i))
    modal_mass = dot_product(modes(:, i), matmul(m_red, modes(:, i)))
    modes(:, i) = modes(:, i) / sqrt(modal_mass)
    sgn = mode_orientation(modes(:, i), nred)
    modes(:, i) = sgn * modes(:, i)
    freq(i) = sqrt(max(eigval(i), 0.0_dp))
    modal_amp(i) = alpha * v0_bar * modes(map_full_to_red(c_dof), i) / freq(i)
  end do

  allocate(xi(nn), y_nodes(nn, nmodes_out))
  do i = 1, nn
    xi(i) = h * real(i - 1, dp)
    do j = 1, nmodes_out
      y_nodes(i, j) = full_mode_value(modes(:, j), red_to_full, nred, 2 * (i - 1) + 1)
    end do
  end do

  call write_parameters(nred, h, b_node, c_node, d_node)
  call write_frequencies(freq, modes, m_red, map_full_to_red, c_dof, d_dof, nmodes_out, nred)
  call write_modes(xi, y_nodes, nn, nmodes_out)
  call write_c_response(freq, modal_amp, modes, map_full_to_red(c_dof), nmodes_out)
  call write_snapshots(freq, modal_amp, xi, y_nodes, nn, nmodes_out)

contains

  subroutine assemble_beam(k, m, h)
    real(dp), intent(inout) :: k(:, :), m(:, :)
    real(dp), intent(in) :: h
    real(dp) :: ke(4, 4), me(4, 4)
    integer :: e, a, b, edof(4)

    ke = reshape([ &
      12.0_dp / h**3,  6.0_dp / h**2, -12.0_dp / h**3,  6.0_dp / h**2, &
       6.0_dp / h**2,  4.0_dp / h,    -6.0_dp / h**2,  2.0_dp / h, &
     -12.0_dp / h**3, -6.0_dp / h**2,  12.0_dp / h**3, -6.0_dp / h**2, &
       6.0_dp / h**2,  2.0_dp / h,    -6.0_dp / h**2,  4.0_dp / h ], [4, 4])

    me = h / 420.0_dp * reshape([ &
      156.0_dp,  22.0_dp * h,   54.0_dp, -13.0_dp * h, &
       22.0_dp * h, 4.0_dp * h**2, 13.0_dp * h, -3.0_dp * h**2, &
       54.0_dp,  13.0_dp * h,  156.0_dp, -22.0_dp * h, &
      -13.0_dp * h, -3.0_dp * h**2, -22.0_dp * h, 4.0_dp * h**2 ], [4, 4])

    do e = 0, nel - 1
      edof = [2 * e + 1, 2 * e + 2, 2 * (e + 1) + 1, 2 * (e + 1) + 2]
      do a = 1, 4
        do b = 1, 4
          k(edof(a), edof(b)) = k(edof(a), edof(b)) + ke(a, b)
          m(edof(a), edof(b)) = m(edof(a), edof(b)) + me(a, b)
        end do
      end do
    end do
  end subroutine assemble_beam

  subroutine cholesky_lower(a, l, n)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n, n)
    real(dp), intent(out) :: l(n, n)
    integer :: i, j, k
    real(dp) :: sumv
    l = 0.0_dp
    do i = 1, n
      do j = 1, i
        sumv = a(i, j)
        do k = 1, j - 1
          sumv = sumv - l(i, k) * l(j, k)
        end do
        if (i == j) then
          if (sumv <= 0.0_dp) stop 'Mass matrix is not positive definite.'
          l(i, j) = sqrt(sumv)
        else
          l(i, j) = sumv / l(j, j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine invert_lower(l, inv_l, n)
    integer, intent(in) :: n
    real(dp), intent(in) :: l(n, n)
    real(dp), intent(out) :: inv_l(n, n)
    integer :: i, j, k
    real(dp) :: sumv
    inv_l = 0.0_dp
    do j = 1, n
      do i = 1, n
        sumv = merge(1.0_dp, 0.0_dp, i == j)
        do k = 1, i - 1
          sumv = sumv - l(i, k) * inv_l(k, j)
        end do
        inv_l(i, j) = sumv / l(i, i)
      end do
    end do
  end subroutine invert_lower

  subroutine jacobi_symmetric(a, d, v, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(n, n)
    real(dp), intent(out) :: d(n), v(n, n)
    integer :: i, j, p, q, iter, max_iter
    real(dp) :: apq, app, aqq, tau, t, c, s, temp, max_off

    v = 0.0_dp
    do i = 1, n
      v(i, i) = 1.0_dp
    end do
    max_iter = 80 * n * n
    do iter = 1, max_iter
      max_off = 0.0_dp
      p = 1
      q = 2
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(a(i, j)) > max_off) then
            max_off = abs(a(i, j))
            p = i
            q = j
          end if
        end do
      end do
      if (max_off < 1.0e-9_dp) exit

      app = a(p, p)
      aqq = a(q, q)
      apq = a(p, q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t * t)
      s = t * c

      do j = 1, n
        if (j /= p .and. j /= q) then
          temp = a(j, p)
          a(j, p) = c * temp - s * a(j, q)
          a(p, j) = a(j, p)
          a(j, q) = s * temp + c * a(j, q)
          a(q, j) = a(j, q)
        end if
      end do
      a(p, p) = c * c * app - 2.0_dp * s * c * apq + s * s * aqq
      a(q, q) = s * s * app + 2.0_dp * s * c * apq + c * c * aqq
      a(p, q) = 0.0_dp
      a(q, p) = 0.0_dp

      do j = 1, n
        temp = v(j, p)
        v(j, p) = c * temp - s * v(j, q)
        v(j, q) = s * temp + c * v(j, q)
      end do
    end do
    do i = 1, n
      d(i) = a(i, i)
    end do
  end subroutine jacobi_symmetric

  subroutine sort_eigensystem(d, v, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: d(n), v(n, n)
    integer :: i, j, p
    real(dp) :: tmp, col(n)
    do i = 1, n - 1
      p = i
      do j = i + 1, n
        if (d(j) < d(p)) p = j
      end do
      if (p /= i) then
        tmp = d(i)
        d(i) = d(p)
        d(p) = tmp
        col = v(:, i)
        v(:, i) = v(:, p)
        v(:, p) = col
      end if
    end do
  end subroutine sort_eigensystem

  real(dp) function mode_orientation(mode, n)
    integer, intent(in) :: n
    real(dp), intent(in) :: mode(n)
    integer :: i, imax
    real(dp) :: vmax
    imax = 1
    vmax = abs(mode(1))
    do i = 2, n
      if (abs(mode(i)) > vmax) then
        vmax = abs(mode(i))
        imax = i
      end if
    end do
    if (mode(imax) >= 0.0_dp) then
      mode_orientation = 1.0_dp
    else
      mode_orientation = -1.0_dp
    end if
  end function mode_orientation

  real(dp) function full_mode_value(mode, red_to_full, n, full_dof)
    integer, intent(in) :: n, red_to_full(n), full_dof
    real(dp), intent(in) :: mode(n)
    integer :: i
    full_mode_value = 0.0_dp
    do i = 1, n
      if (red_to_full(i) == full_dof) then
        full_mode_value = mode(i)
        return
      end if
    end do
  end function full_mode_value

  subroutine write_parameters(nred, h, b_node, c_node, d_node)
    integer, intent(in) :: nred, b_node, c_node, d_node
    real(dp), intent(in) :: h
    integer :: u
    open(newunit=u, file='outputs/fem/fem_parameters.csv', status='replace', action='write')
    write(u, '(a)') 'parameter,value'
    write(u, '(a,i0)') 'nel,', nel
    write(u, '(a,i0)') 'nodes,', nn
    write(u, '(a,i0)') 'reduced_dofs,', nred
    write(u, '(a,es24.16)') 'h,', h
    write(u, '(a,es24.16)') 'alpha,', alpha
    write(u, '(a,es24.16)') 'kappa,', kappa
    write(u, '(a,es24.16)') 'v0_bar,', v0_bar
    write(u, '(a,i0)') 'b_node,', b_node
    write(u, '(a,i0)') 'c_node,', c_node
    write(u, '(a,i0)') 'd_node,', d_node
    close(u)
  end subroutine write_parameters

  subroutine write_frequencies(freq, modes, m_red, map, c_dof, d_dof, nm, n)
    integer, intent(in) :: nm, n, map(:), c_dof, d_dof
    real(dp), intent(in) :: freq(nm), modes(n, nm), m_red(n, n)
    integer :: u, i
    real(dp) :: mm, yc, yd
    open(newunit=u, file='outputs/fem/fem_frequencies.csv', status='replace', action='write')
    write(u, '(a)') 'mode,omega_bar,beta,mass_norm,Y_C,Y_D'
    do i = 1, nm
      mm = dot_product(modes(:, i), matmul(m_red, modes(:, i)))
      yc = modes(map(c_dof), i)
      yd = modes(map(d_dof), i)
      write(u, '(i0,5(",",es24.16))') i, freq(i), sqrt(freq(i)), mm, yc, yd
    end do
    close(u)
  end subroutine write_frequencies

  subroutine write_modes(xi, y_nodes, nnode, nm)
    integer, intent(in) :: nnode, nm
    real(dp), intent(in) :: xi(nnode), y_nodes(nnode, nm)
    integer :: u, i, j
    real(dp) :: maxv(nm)
    do j = 1, nm
      maxv(j) = maxval(abs(y_nodes(:, j)))
    end do
    open(newunit=u, file='outputs/fem/fem_modes.csv', status='replace', action='write')
    write(u, '(a)', advance='no') 'xi'
    do j = 1, nm
      write(u, '(a,i0)', advance='no') ',mode_', j
    end do
    write(u, *)
    do i = 1, nnode
      write(u, '(es24.16)', advance='no') xi(i)
      do j = 1, nm
        write(u, '(",",es24.16)', advance='no') y_nodes(i, j) / maxv(j)
      end do
      write(u, *)
    end do
    close(u)
  end subroutine write_modes

  subroutine write_c_response(freq, amp, modes, c_red, nm)
    integer, intent(in) :: nm, c_red
    real(dp), intent(in) :: freq(nm), amp(nm), modes(:, :)
    integer :: u, i, j
    real(dp) :: tau, disp, vel
    open(newunit=u, file='outputs/fem/fem_c_response.csv', status='replace', action='write')
    write(u, '(a)') 'tau,w_c,v_c'
    do i = 1, nt
      tau = tau_max * real(i - 1, dp) / real(nt - 1, dp)
      disp = 0.0_dp
      vel = 0.0_dp
      do j = 1, nm
        disp = disp + modes(c_red, j) * amp(j) * sin(freq(j) * tau)
        vel = vel + modes(c_red, j) * amp(j) * freq(j) * cos(freq(j) * tau)
      end do
      write(u, '(3(es24.16,:,","))') tau, disp, vel
    end do
    close(u)
  end subroutine write_c_response

  subroutine write_snapshots(freq, amp, xi, y_nodes, nnode, nm)
    integer, intent(in) :: nnode, nm
    real(dp), intent(in) :: freq(nm), amp(nm), xi(nnode), y_nodes(nnode, nm)
    integer :: u, it, i, j
    real(dp) :: tau, disp
    open(newunit=u, file='outputs/fem/fem_snapshots.csv', status='replace', action='write')
    write(u, '(a)') 'tau,xi,w'
    do it = 1, nsnap
      tau = tau_max * real(it - 1, dp) / real(nsnap - 1, dp)
      do i = 1, nnode
        disp = 0.0_dp
        do j = 1, nm
          disp = disp + y_nodes(i, j) * amp(j) * sin(freq(j) * tau)
        end do
        write(u, '(3(es24.16,:,","))') tau, xi(i), disp
      end do
    end do
    close(u)
  end subroutine write_snapshots

end program fem_q3_solver
