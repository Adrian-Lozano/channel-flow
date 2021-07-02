!------------------------------------------!
!        Module for interpolation          !
!------------------------------------------!
Module interpolation

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global,          Only : weight_y_0, weight_y_1

  ! prevent implicit typing
  Implicit None

Contains

  !-------------------------------------------------------!
  !            Linear interpolation of u in x             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  ! Parameter: di-> 1=faces to center  (normal   average) !
  !                 2=centers to faces (weighted average) !
  !-------------------------------------------------------!
  ! uniform mesh assumed
  Subroutine interpolate_x(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer(Int32) :: n(3), n1, n2, n3

    n  = Shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ui = 0d0
    ui( 1:n1-1, 1:n2, 1:n3) = 0.5d0*( u(1:n1-1,:,:) + u(2:n1,:,:) )

  End Subroutine interpolate_x

  !-------------------------------------------------------!
  !            Linear interpolation of u in y             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  ! Parameter: di-> 1=faces to center  (normal   average) !
  !                 2=centers to faces (weighted average) !
  !-------------------------------------------------------!
  Subroutine interpolate_y(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer (Int32) :: n(3), n1, n2, n3, i1, i3

    n  = Shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ui = 0d0
    If ( di==1 ) Then
       ! faces to centers (normal average)
       ui(1:n1, 1:n2-1, 1:n3) = 0.5d0*( u(:,1:n2-1,:) + u(:,2:n2,:) ) 
    Elseif ( di==2 ) Then
       ! centers to faces (weighted average)
       Do i1 = 1, n1
          Do i3 = 1, n3
             ui(i1, 1:n2-1, i3) = weight_y_0*u(i1,1:n2-1,i3) + weight_y_1*u(i1,2:n2,i3)
          End Do
       End Do
    Else
       Stop 'Error: invalid interpolation'
    End If
   
  End Subroutine interpolate_y

  !-------------------------------------------------------!
  !            Linear interpolation of u in z             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  ! Parameter: di-> 1=faces to center  (normal   average) !
  !                 2=centers to faces (weighted average) !
  !-------------------------------------------------------!
  ! uniform mesh assumed
  Subroutine interpolate_z(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer(Int32) :: n(3), n1, n2, n3

    n  = Shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ui = 0d0
    ui(1:n1, 1:n2, 1:n3-1) = 0.5d0*( u(:,:,1:n3-1) + u(:,:,2:n3) )
    
  End Subroutine interpolate_z

  !-------------------------------------------------------!
  !        General second order interpolation in x        !
  !                                                       !
  ! Input : y,u,yi                                        !
  ! Output: ui                                            !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine interpolate_x_2nd(y,u,yi,ui)

    Real    (Int64), Intent(In)  :: u (:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Real    (Int64), Intent(In)  :: y(:), yi(:)

    Integer(Int32) :: n(3), n1, n2, n3, j, jj, nn(1), ny, j2
    Real   (Int64) :: w0, w1, w2, yref

    n  = shape(u)
    n2 = n(2)
    n3 = n(3)

    nn = shape(yi)
    n1 = nn(1)

    ny = size(y)

    If ( ny/=n(1) ) Stop 'Error wrong size!'

    ui = 0d0    
    ! middle points
    Do j = 2, n1-1
       yref = yi(j)
       Do j2 = 3, ny
          If ( y(j2)>=yref ) Then
             jj = j2-1
             Exit
          End If
       End Do
       w0   = ( yref - y(jj)   )*( yref - y(jj+1) )/( ( y(jj-1) - y(jj)   )*( y(jj-1) - y(jj+1) ) )
       w1   = ( yref - y(jj-1) )*( yref - y(jj+1) )/( ( y(jj)   - y(jj-1) )*( y(jj)   - y(jj+1) ) )
       w2   = ( yref - y(jj-1) )*( yref - y(jj)   )/( ( y(jj+1) - y(jj-1) )*( y(jj+1) - y(jj)   ) )
       ui(j, 1:n2, 1:n3) = w0*u(jj-1,1:n2,1:n3) + w1*u(jj,1:n2,1:n3) +  w2*u(jj+1,1:n2,1:n3)
    End Do

    ! first point
    yref = yi(1)
    jj   = 2
    Do j2 = 3, ny
       If ( y(j2)>=yref ) Then
          jj = j2-1
          Exit
       End If
    End Do
    If ( jj==0 ) jj = ny
    w0   = ( yref - y(jj)   )*( yref - y(jj+1) )/( ( y(jj-1) - y(jj)   )*( y(jj-1) - y(jj+1) ) )
    w1   = ( yref - y(jj-1) )*( yref - y(jj+1) )/( ( y(jj  ) - y(jj-1) )*( y(jj)   - y(jj+1) ) )
    w2   = ( yref - y(jj-1) )*( yref - y(jj  ) )/( ( y(jj+1) - y(jj-1) )*( y(jj+1) - y(jj)   ) )
    ui(1, 1:n2, 1:n3) = w0*u(jj-1,1:n2,1:n3) + w1*u(jj,1:n2,1:n3) +  w2*u(jj+1,1:n2,1:n3)

    ! last point
    yref = yi(n1)
    jj   = ny-1
    Do j2 = 3, ny
       If ( y(j2)>=yref ) Then
          jj = j2-1
          Exit
       End If
    End Do
    w0   = ( yref - y(jj)   )*( yref - y(jj+1) )/( ( y(jj-1) - y(jj)   )*( y(jj-1) - y(jj+1) ) )
    w1   = ( yref - y(jj-1) )*( yref - y(jj+1) )/( ( y(jj  ) - y(jj-1) )*( y(jj)   - y(jj+1) ) )
    w2   = ( yref - y(jj-1) )*( yref - y(jj  ) )/( ( y(jj+1) - y(jj-1) )*( y(jj+1) - y(jj)   ) )
    ui(n1, 1:n2, 1:n3) = w0*u(jj-1,1:n2,1:n3) + w1*u(jj,1:n2,1:n3) +  w2*u(jj+1,1:n2,1:n3)

  End Subroutine interpolate_x_2nd

  !-------------------------------------------------------!
  !        General second order interpolation in z        !
  !                                                       !
  ! Input : y,u,yi                                        !
  ! Output: ui                                            !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine interpolate_z_2nd(y,u,yi,ui)

    Real    (Int64), Intent(In)  :: u (:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Real    (Int64), Intent(In)  :: y(:), yi(:)

    Integer(Int32) :: n(3), n1, n2, n3, j, jj, nn(1), ny, j2
    Real   (Int64) :: w0, w1, w2, yref

    n  = shape(u)
    n1 = n(1)
    n2 = n(2)

    nn = shape(yi)
    n3 = nn(1)

    ny = size(y)

    If ( ny/=n(3) ) Stop 'Error wrong size!'

    ui = 0d0    
    ! middle points
    Do j = 2, n3-1
       yref = yi(j)
       Do j2 = 3, ny
          If ( y(j2)>yref ) Then
             jj = j2-1
             Exit
          End If
       End Do
       w0   = ( yref - y(jj)   )*( yref - y(jj+1) )/( ( y(jj-1) - y(jj)   )*( y(jj-1) - y(jj+1) ) )
       w1   = ( yref - y(jj-1) )*( yref - y(jj+1) )/( ( y(jj)   - y(jj-1) )*( y(jj)   - y(jj+1) ) )
       w2   = ( yref - y(jj-1) )*( yref - y(jj)   )/( ( y(jj+1) - y(jj-1) )*( y(jj+1) - y(jj)   ) )
       ui(1:n1, 1:n2, j) = w0*u(1:n1,1:n2,jj-1) + w1*u(1:n1,1:n2,jj) + w2*u(1:n1,1:n2,jj+1)
    End Do

    ! first point
    yref = yi(1)
    jj   = 2
    Do j2 = 3, ny
       If ( y(j2)>yref ) Then
          jj = j2-1
          Exit
       End If
    End Do
    w0   = ( yref - y(jj)   )*( yref - y(jj+1) )/( ( y(jj-1) - y(jj)   )*( y(jj-1) - y(jj+1) ) )
    w1   = ( yref - y(jj-1) )*( yref - y(jj+1) )/( ( y(jj  ) - y(jj-1) )*( y(jj)   - y(jj+1) ) )
    w2   = ( yref - y(jj-1) )*( yref - y(jj  ) )/( ( y(jj+1) - y(jj-1) )*( y(jj+1) - y(jj)   ) )
    ui(1:n1, 1:n2, 1) = w0*u(1:n1,1:n2,jj-1) + w1*u(1:n1,1:n2,jj) + w2*u(1:n1,1:n2,jj+1)

    ! last point
    yref = yi(n3)
    jj   = ny-1
    Do j2 = 3, ny
       If ( y(j2)>yref ) Then
          jj = j2-1
          Exit
       End If
    End Do
    w0   = ( yref - y(jj)   )*( yref - y(jj+1) )/( ( y(jj-1) - y(jj)   )*( y(jj-1) - y(jj+1) ) )
    w1   = ( yref - y(jj-1) )*( yref - y(jj+1) )/( ( y(jj  ) - y(jj-1) )*( y(jj)   - y(jj+1) ) )
    w2   = ( yref - y(jj-1) )*( yref - y(jj  ) )/( ( y(jj+1) - y(jj-1) )*( y(jj+1) - y(jj)   ) )
    ui( 1:n1, 1:n2, n3) = w0*u(1:n1,1:n2,jj-1) + w1*u(1:n1,1:n2,jj) + w2*u(1:n1,1:n2,jj+1)

  End Subroutine interpolate_z_2nd

End Module interpolation
