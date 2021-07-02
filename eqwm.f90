!------------------------------------------------!
! EQWM                                           !
!------------------------------------------------!
Module EQWM
  
  ! Modules
  Use iso_fortran_env, Only : Int32, Int64
  !Use interpolation

  ! prevent implicit typing
  Implicit None

  ! private and public variables
  Private
  Public :: solve_u_eqwm
  
  ! Declarations
Contains 
  
  Subroutine grid(ncv, stretch_factor, h_wm, y_cv, y_fa)
    
    Real   (Int64), Intent(In) :: stretch_factor, h_wm
    Integer(Int32), Intent(In) :: ncv        
    Real   (Int64), Dimension(:), Intent(Out) :: y_cv, y_fa

    Integer(Int32) :: nfa, ifa, icv
    Real   (Int64) :: y_max, tmp_inc
    
    nfa = ncv + 1
    
    y_fa(1) = 0.0d0
    tmp_inc = 1.0d0
    
    Do ifa = 2,nfa
       y_fa(ifa) = y_fa(ifa-1) + tmp_inc
       y_cv(ifa-1) = 0.5d0* ( y_fa(ifa) + y_fa(ifa-1) )
       tmp_inc = tmp_inc*stretch_factor
    End Do
    
    ! rescale grid so that y_cv(max) = h_wm
    ! LES condition is enforced at the centroid of the last cell. 
    y_max = Maxval(y_cv)
    Do ifa = 1,nfa
       y_fa(ifa) = y_fa(ifa)/y_max*h_wm
    Enddo
    Do icv = 1,ncv
       y_cv(icv) = y_cv(icv)/y_max*h_wm
    Enddo
        
  End Subroutine grid
  
  Subroutine thomas_algorithm(n_ta,a_ta,b_ta,c_ta,d_ta,x_ta)
    !        a - sub-diagonal (means it is the diagonal below the main diagonal)
    !        b - the main diagonal
    !        c - sup-diagonal (means it is the diagonal above the main diagonal)
    !        d - right part
    !        x - the answer
    !        n - number of equations
    
    Integer(Int32), Intent(In) :: n_ta
    Real   (Int64), Dimension(n_ta), Intent(In) :: a_ta, b_ta, c_ta, d_ta
    Real   (Int64), Dimension(n_ta), Intent(Out) :: x_ta
    Real   (Int64), Dimension(n_ta) :: cp_ta, dp_ta
    Real   (Int64) :: m_ta
    Integer(Int32) :: i_ta, j_ta
    
    ! initialize c-prime and d-prime
    cp_ta(1) = c_ta(1)/b_ta(1)
    dp_ta(1) = d_ta(1)/b_ta(1)
    
    ! solve for vectors c-prime and d-prime
    Do i_ta = 2,n_ta
       m_ta = b_ta(i_ta)-cp_ta(i_ta-1)*a_ta(i_ta)
       cp_ta(i_ta) = c_ta(i_ta)/m_ta
       dp_ta(i_ta) = (d_ta(i_ta)-dp_ta(i_ta-1)*a_ta(i_ta))/m_ta
       !if (a_ta(i_ta) /= a_ta(i_ta)) then
       ! write(*,*) m_ta
       !end if
    Enddo
    ! initialize x
    x_ta(n_ta) = dp_ta(n_ta)
    ! solve for x from the vectors c-prime and d-prime
    Do j_ta = n_ta-1, 1, -1
       x_ta(j_ta) = dp_ta(j_ta)-cp_ta(j_ta)*x_ta(j_ta+1)
       !write(*,*) j_ta
    End do
    
  End subroutine thomas_algorithm
  
  Subroutine solve_u_eqwm(umag_les, mu_lam, rho_eq, h_wm, ncv, stretch_factor, &
       tau_wall, iter, tau_wall_in, tt)
    
    Real   (Int64), Intent(In)  :: umag_les, mu_lam, rho_eq, h_wm
    Real   (Int64), Intent(In)  :: stretch_factor, tau_wall_in, tt
    Real   (Int64), Intent(Out) :: tau_wall
    Integer(Int32), Intent(In)  :: ncv
    Integer(Int32), Intent(Out) :: iter

    Real   (Int64) :: kappa, Aplus, tol, max_iter, tau_wall_lam
    Real   (Int64) :: tau_wall_prev, nu_eq, utau_eq, D, mut, diag, superdiag, subdiag
    Integer(Int32) :: done, nfa, i, ifa, icv, j

    Real(Int64), Dimension(:),   Allocatable :: y_cv, y_fa, ones_vect, mu_fa
    Real(Int64), Dimension(:,:), Allocatable :: A
    Real(Int64), Dimension(:),   Allocatable :: b, aa, bb, cc, u_wm
    
    kappa    = 0.41d0
    Aplus    = 17.0d0
    tol      = 0.0000001d0
    max_iter = 100
    
    ! reference laminar solution
    tau_wall_lam = mu_lam*umag_les/h_wm
    
    ! build implicit WM grid for ODE solve 
    Allocate(y_cv(ncv+1), y_fa(ncv+1))

    Call grid( ncv, stretch_factor, h_wm, y_cv, y_fa )
    nfa = Size(y_fa)
    !if (tt>=0.001293302871109800) then
    !            write(*,*) y_cv(2), y_fa(2)
    !end if
    iter = 0
    done = 0
    tau_wall = tau_wall_in
    !write(*,*) tau_wall
    Allocate(ones_vect(nfa), mu_fa(nfa))
    Allocate(A(ncv,ncv))
    Allocate(b(ncv))
    Allocate(aa(ncv), bb(ncv), cc(ncv), u_wm(ncv))
    Do while (done == 0)
       tau_wall_prev = tau_wall
       !write(*,*) iter, tau_wall_prev
       !
       !                ! populate total viscosity at face
       Do i = 1, nfa
          ones_vect(i) = 1.0d0
       Enddo
       mu_fa = mu_lam * ones_vect
       Do ifa=2,nfa
          nu_eq         = mu_lam / rho_eq
          utau_eq       = (tau_wall / rho_eq)**0.5d0
          D          = 1.0d0-dexp(-1d0*y_fa(ifa)*utau_eq/nu_eq/Aplus)
          D          = D**2.0d0
          mut        = rho_eq * kappa * y_fa(ifa) * utau_eq * D
          mu_fa(ifa) = mu_fa(ifa) + mut
       Enddo

       ! construct momentum matrix system
       ! solution vector is u_wm at y_cv. 
       ! top bc .. 
       A(ncv, ncv) = 1.0d0              
       b(ncv)      = umag_les

       ! internal cvs 
       Do icv = 2,ncv-1
          superdiag = mu_fa(icv+1)/( y_cv(icv+1) - y_cv(icv  ) )
          subdiag   = mu_fa(icv  )/( y_cv(icv  ) - y_cv(icv-1) )
          diag      = -1.0d0* (superdiag+subdiag)
          A(icv,icv) = diag
          A(icv, icv+1) = superdiag
          A(icv, icv-1) = subdiag
          b(icv)        = 0.0d0
       Enddo

       ! wall bc (enforced in eqn for icv=1, wall CV). 
       b(1)      = 0.0d0
       superdiag = mu_fa(2)/( y_cv(2)-y_cv(1) )
       diag      = -1.0d0* ( superdiag + mu_fa(1)/( y_cv(1)-y_fa(1) ) )
       A(1,1)    = diag
       A(1,2)    = superdiag
       A(ncv, ncv-1) = 0.0d0

       aa(1)   = 0.0d0
       cc(ncv) = 0.0d0
       
       Do i = 1, ncv
          Do j = 1, ncv
             If (i==j) then
                bb(i) = A(i,j)
             Else if (i-1==j) then
                aa(i) = A(i,j)
                !if (abs(A(i,j))>100) then
                !write(*,*) A(i,j), i, j
                !end if
                !tt = tt+1
             Else if (i+1==j) then
                cc(i) = A(i,j)
                !write(*,*) j
                !ttt = ttt+1
             End if
          End do
       End do

       Call thomas_algorithm(ncv,aa,bb,cc,b,u_wm)
       !                u_wm = A\b;
       !                % update tau_wall
       tau_wall = mu_fa(1) * ( u_wm(1) - 0.0d0)/(y_cv(1)-y_fa(1))
       !if (tau_wall /= tau_wall .and. iter==0) then
       ! write(*,*) tau_wall
       !end if      
       !write(*,*) iter, tau_wall
 
       If ( Abs( (tau_wall - tau_wall_prev)/tau_wall_lam) < tol) then
          done = 1
       End if

       If (done == 0 ) then
          iter = iter + 1
       End if
       If (iter > max_iter ) then
          done = 1       
       End if

    End do

    Deallocate(y_cv, y_fa)
    Deallocate(ones_vect, mu_fa)
    Deallocate(A)
    Deallocate(b)
    Deallocate(aa, bb, cc, u_wm)
    
  End Subroutine solve_u_eqwm
  
End Module EQWM
