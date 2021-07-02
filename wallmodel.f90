!------------------------------------!
!     Module for LES wall-models     !
!------------------------------------!
Module wallmodel

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use interpolation
  Use subgrid
  Use boundary_conditions
  Use EQWM

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------!
  !                                              !
  !              Select wall model               !
  !                                              !
  !----------------------------------------------!
  Subroutine compute_wall_model(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    ! select model
    If     ( iwall_model==1 ) Then 
       ! set constant
       Call compute_constant_alpha
    Elseif ( iwall_model==2 ) Then
       ! Exact wall model with Robin boundary condition
       Call compute_exact_Robin(U_,V_,W_,nu_t_)
    Elseif ( iwall_model==3 ) Then
       ! Equilibrium wall model (Larsson)
       Call compute_equilibrium_wall_model(U_,V_,W_,nu_t_)
    Elseif ( iwall_model==4 ) Then 
       ! Exact wall model with Neumann boundary condition
       Call compute_exact_Neumann(U_,V_,W_,nu_t_) 
    Elseif ( iwall_model==5 ) Then
       ! Boses's wall model
       Call compute_bose_model(U_,V_,W_,nu_t_)       
    Elseif ( iwall_model==6 ) Then 
       ! generic model 
       Call compute_generic_dynamic_model(U_,V_,W_,nu_t_,gwm_cu,gwm_c_test,gwm_diag,gwm_ai)     
    Elseif ( iwall_model==7 ) Then 
       ! generic model 
       Call compute_generic_dynamic_model_2(U_,V_,W_,nu_t_,gwm_diag,gwm_cu,gwm_ai)
    Elseif ( iwall_model==8 ) Then 
       ! generic model 
       Call compute_generic_dynamic_model_3(U_,V_,W_,nu_t_,gwm_diag,gwm_cu,gwm_ai)
    End If
    
    ! penetration vs. no penetration at the wall
    ! 1 -> penetration    (V/=0) 
    ! 0 -> no penetration (V =0)
    If ( penetration == 0 ) Then
       alpha_y = 0d0
    End If

    ! alpha limiter
    If ( .False. ) Then
       alpha_x = Max( alpha_x(2,1,2), 0d0 )
       alpha_x = Min( alpha_x(2,1,2), 1d0 )
       
       alpha_y = Max( alpha_y(2,1,2), 0d0 )
       alpha_y = Min( alpha_y(2,1,2), 1d0 )
       
       alpha_z = Max( alpha_z(2,1,2), 0d0 )
       alpha_z = Min( alpha_z(2,1,2), 1d0 )
    End If

    ! compute boundary conditions for pseudo-pressure
    If ( iwall_model>0 ) Then
       Call compute_pseudo_pressure_bc_for_robin_bc
    End If
    
  End Subroutine compute_wall_model

  !--------------------------------------------------------------------!
  !                                                                    !
  !                        Boses' wall-model                           !
  !                                                                    !
  !          Compute slip lengths alpha_x, alpha_y and alpha_z         !
  !                  for ui = alpha_i dui/dy at the wall               !
  !                                                                    !
  ! Equation:                                                          !
  !                                                                    !
  !  alpha^2*( fil_size_2^2*dhat(ui)/dn*dhat(uj)/dn - dui/dn*duj/dn ) +!
  !            Tij - hat(tau_ij) = hat(ui*uj) - ui*uj                  !
  !                                                                    !
  !  n->normal to the wall                                             !
  !                                                                    !
  ! Definitions:                                                       !
  !                                                                    !
  !  LWij    = hat(ui*uj) - ui*uj - Tij + hat(tau_ij)                  !
  !  MWij    = fil_size_2^2*hat(dui/dn)*hat(duj/dn) - dui/dn*duj/dn    !
  !  tau_ij = -2*nu_t*SWij                                             !
  !  Tij    = -2*nu_t_hat*hat(SWij)                                    !
  !                                                                    !
  ! alpha = alpha_x = alpha_y = alpha_z                                !
  ! alpha computed at V positions                                      !
  !                                                                    !
  ! Input:  U_,V_,W_,nu_t (flow velocities)                            !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)                   !
  !                                                                    !
  !--------------------------------------------------------------------!
  Subroutine compute_bose_model(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    ! local variables
    Integer(Int32) :: i, j, k

    !--------------------------------------------------------------------------!
    ! Part 1: Compute LWij = hat(ui*uj) - ui*uj
    
    ! interpolate velocity to U and W to V location (at wall)
    call interpolate_x(U_,     term_1(2:nx,:,:),1) 
    call interpolate_y(term_1, term  (:,1:ny,:),2)
    
    call interpolate_z(W_,     term_1(:,:,2:nz),1) 
    call interpolate_y(term_1, term_2(:,1:ny,:),2) 

    ! fill in missing values (periodicity)
    Call apply_periodic_bc_x(term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! SWij = ui*uj at v location (at wall)
    ! U at V location -> term  (2:nx-1,1:ny,:)
    ! W at V location -> term_2(:,1:ny,2:nzg-1)
    ! bottom wall 
    SWij(:,2:3,:,1) = term  (2:nxg,1:2,2:nzg) * term  (2:nxg,1:2,2:nzg)  ! u^2
    SWij(:,2:3,:,2) = V_    (2:nxg,1:2,2:nzg) * V_    (2:nxg,1:2,2:nzg)  ! v^2
    SWij(:,2:3,:,3) = term_2(2:nxg,1:2,2:nzg) * term_2(2:nxg,1:2,2:nzg)  ! w^2
    SWij(:,2:3,:,4) = term  (2:nxg,1:2,2:nzg) * V_    (2:nxg,1:2,2:nzg)  ! uv
    SWij(:,2:3,:,5) = term  (2:nxg,1:2,2:nzg) * term_2(2:nxg,1:2,2:nzg)  ! uw
    SWij(:,2:3,:,6) = V_    (2:nxg,1:2,2:nzg) * term_2(2:nxg,1:2,2:nzg)  ! vw
    ! top wall
    SWij(:,nyg-2:nyg-1,:,1) = term  (2:nxg,ny-1:ny,2:nzg) * term  (2:nxg,ny-1:ny,2:nzg) ! u^2
    SWij(:,nyg-2:nyg-1,:,2) = V_    (2:nxg,ny-1:ny,2:nzg) * V_    (2:nxg,ny-1:ny,2:nzg) ! v^2
    SWij(:,nyg-2:nyg-1,:,3) = term_2(2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! w^2
    SWij(:,nyg-2:nyg-1,:,4) = term  (2:nxg,ny-1:ny,2:nzg) * V_    (2:nxg,ny-1:ny,2:nzg) ! uv
    SWij(:,nyg-2:nyg-1,:,5) = term  (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! uw
    SWij(:,nyg-2:nyg-1,:,6) = V_    (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! vw

    ten_buf(2:nxg,2:nyg,2:nzg,1:6) = SWij;
    Do i = 1,6 
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do

    ! Compute LWij = hat(ui*uj) -> stored in LWij(3:nxg-1,2:nyg-1,3:nzg-1,:)
    Call filter_tensor_xzy(ten_buf(1:nxg,        2:5,1:nzg,:),LWij(2:nxg-1,        2:5,2:nzg-1,:)) 
    Call filter_tensor_xzy(ten_buf(1:nxg,nyg-4:nyg-1,1:nzg,:),LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:)) 

    ! Compute LWij = hat(ui*uj) - ui*uj at the wall 
    ! Valid for LWij(2:nxg-1,(/2,nyg-1/),2:nzg-1)
    LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) - SWij(2:nxg-1,    2,2:nzg-1,:)
    LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-1,2:nzg-1,:)

    !-----------------------------------------------------------------------------!
    ! Part 2: Compute -Tij + hat(tau_ij) at the wall and add it to LWij
    If ( Dirichlet_nu_t==0 ) Then
       
       If (est_Tij == 0) Then
          
          !----------compute hat(tau_ij)       
          ! compute SWij at V locations in the first two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
          Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
          
          ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
          Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)
          
          ! tau_ij = -2*nu_t*SWij -> SWij
          Do i = 2, nxg
             Do k = 2, nzg
                Do j = 1, 6
                   SWij(i,      2:3,k,j) = -2d0*term_2(i,    1:2,k)*SWij(i,      2:3,k,j) ! bottom 
                   SWij(i,nyg-1:nyg,k,j) = -2d0*term_2(i,ny-1:ny,k)*SWij(i,nyg-1:nyg,k,j) ! top
                End Do
             End Do
          End Do
          
          ten_buf(2:nxg,2:nyg,2:nzg,1:6) = SWij;
          Do i = 1,6
             Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
             Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
             Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
          End Do
          
          ! hat(tau_ij) -> filter tau_ij 
          !     tau_ij  -> SWij(:,2:4,:) (bottom) and SWij(:,nyg-2:nyg  ,:)
          ! hat(tau_ij) -> MWij(:,2:4,:) (bottom) and MWij(:,nyg-3:nyg-1,:)
          Call filter_tensor_xzy(ten_buf(1:nxg,      2:4,1:nzg,:),MWij(2:nxg-1,        2:4,2:nzg-1,:)) 
          Call filter_tensor_xzy(ten_buf(1:nxg,nyg-2:nyg,1:nzg,:),MWij(2:nxg-1,nyg-3:nyg-1,2:nzg-1,:)) 
          
          ! add it
          LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + MWij(2:nxg-1,    2,2:nzg-1,:)
          LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + MWij(2:nxg-1,nyg-1,2:nzg-1,:)
          
          !----------compute T_ij
          ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
          Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
          Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
          Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
          
          ! apply periodicity in x and z
          Call apply_periodic_bc_x(Uff,1)
          Call apply_periodic_bc_z(Uff,1)
          Call update_ghost_interior_planes(Uff,1)
          Call apply_periodic_bc_x(Vff,2)
          Call apply_periodic_bc_z(Vff,2)
          Call update_ghost_interior_planes(Vff,2)
          Call apply_periodic_bc_x(Wff,2)
          Call apply_periodic_bc_z(Wff,3)
          Call update_ghost_interior_planes(Wff,3)
          
          ! compute eddy viscosity for filtered velocities
          Call compute_eddy_viscosity(Uff,Vff,Wff,avg_nu_t_hat,nu_t_hat)
          
          ! interpolate nu_t_hat from cell centers to cell faces (== cell edges because averaged in xz)
          Call interpolate_y(nu_t_hat,term_2(:,1:ny,:),2)    
          
          ! compute hat(SWij) at V locations in two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
          Call compute_Sij_at_V_location_wall(Uff,Vff,Wff,SWij)
          
          ! Tij = -2*nu_t_hat*hat(SWij) -> SWij
          Do i = 2, nxg
             Do k = 2, nzg
                Do j = 1, 6
                   SWij(i,  2,k,j) = -2d0*term_2(i, 1,k)*SWij(i,  2,k,j) ! bottom 
                   SWij(i,nyg,k,j) = -2d0*term_2(i,ny,k)*SWij(i,nyg,k,j) ! top
                End Do
             End Do
          End Do
          
          ! add it
          LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) - SWij(2:nxg-1,  2,2:nzg-1,:)
          LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg,2:nzg-1,:)
          
       Elseif (est_Tij == 1) Then
          
          ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
          Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
          Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
          Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
          
          ! apply periodicity in x and z
          Call apply_periodic_bc_x(Uff,1)
          Call apply_periodic_bc_z(Uff,1)
          Call update_ghost_interior_planes(Uff,1)
          Call apply_periodic_bc_x(Vff,2)
          Call apply_periodic_bc_z(Vff,2)
          Call update_ghost_interior_planes(Vff,2)
          Call apply_periodic_bc_x(Wff,2)
          Call apply_periodic_bc_z(Wff,3)
          Call update_ghost_interior_planes(Wff,3)
          
          Call compute_Tij_estimate(Uff,Vff,Wff,SWij)
          
          LWij(2:nxg-1,    2,2:nzg-1,1) = LWij(2:nxg-1,    2,2:nzg-1,1) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,  2,2:nzg-1,1) 
          LWij(2:nxg-1,    2,2:nzg-1,2) = LWij(2:nxg-1,    2,2:nzg-1,2) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,  2,2:nzg-1,2) 
          LWij(2:nxg-1,    2,2:nzg-1,3) = LWij(2:nxg-1,    2,2:nzg-1,3) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,  2,2:nzg-1,3) 
          
          LWij(2:nxg-1,nyg-1,2:nzg-1,1) = LWij(2:nxg-1,nyg-1,2:nzg-1,1) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,nyg,2:nzg-1,1) 
          LWij(2:nxg-1,nyg-1,2:nzg-1,2) = LWij(2:nxg-1,nyg-1,2:nzg-1,2) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,nyg,2:nzg-1,2) 
          LWij(2:nxg-1,nyg-1,2:nzg-1,3) = LWij(2:nxg-1,nyg-1,2:nzg-1,3) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,nyg,2:nzg-1,3) 
          
          LWij(2:nxg-1,    2,2:nzg-1,4) = LWij(2:nxg-1,    2,2:nzg-1,4) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,  2,2:nzg-1,4)
          LWij(2:nxg-1,    2,2:nzg-1,5) = LWij(2:nxg-1,    2,2:nzg-1,5) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,  2,2:nzg-1,5)
          LWij(2:nxg-1,    2,2:nzg-1,6) = LWij(2:nxg-1,    2,2:nzg-1,6) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,  2,2:nzg-1,6) 
          
          LWij(2:nxg-1,nyg-1,2:nzg-1,4) = LWij(2:nxg-1,nyg-1,2:nzg-1,4) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,nyg,2:nzg-1,4)
          LWij(2:nxg-1,nyg-1,2:nzg-1,5) = LWij(2:nxg-1,nyg-1,2:nzg-1,5) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,nyg,2:nzg-1,5)
          LWij(2:nxg-1,nyg-1,2:nzg-1,6) = LWij(2:nxg-1,nyg-1,2:nzg-1,6) - (fil_size_2*Delta)**2d0 * SWij(2:nxg-1,nyg,2:nzg-1,6)
          
       Else ! est_Tij == 2 Tij - \hat{tij} = \hat{uv} - \hat{u}\hat{v}
          
          Stop 'Error: unknown est_Tij!'
          
       End If
       
    End If
    
    !-----------------------------------------------------------------------------!
    ! Part 3: Compute MWij = fil_size_2^2*hat(dui/dy)*hat(duj/dy) - dui/dy*duj/dy

    ! interpolate velocity to cell centers
    call interpolate_x(U_,term  (2:nx,:,:),1) 
    call interpolate_z(W_,term_2(:,:,2:nz),1) 

    Call apply_periodic_bc_x(term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! Compute dui/dn at V location (at wall)
    ! SWij(:,:,:,1) -> dU/dn
    ! SWij(:,:,:,2) -> dV/dn
    ! SWij(:,:,:,3) -> dW/dn 
    Do i = 2, nxg
      Do k = 2, nzg

        SWij(i,2,k,1)     = (term  (i,2,k)     - term  (i,1,k)  ) / (yg(2)-yg(1))             ! dU/dy at lower wall
        SWij(i,2,k,2)     = (V_    (i,2,k)     - V_    (i,1,k)  ) / (y (2)-y (1))             ! dV/dy at lower wall (1st order)
        SWij(i,2,k,3)     = (term_2(i,2,k)     - term_2(i,1,k)  ) / (yg(2)-yg(1))             ! dW/dy at lower wall

        SWij(i,nyg-1,k,1) = -(term  (i,nyg-1,k) - term  (i,nyg,k)) / (yg(nyg)-yg(nyg-1))      ! dU/dy at upper wall
        SWij(i,nyg-1,k,2) = -(V_    (i,ny -1,k) - V_    (i,ny ,k)) / (y ( ny)-y ( ny-1))      ! dV/dy at upper wall (1st order)
        SWij(i,nyg-1,k,3) = -(term_2(i,nyg-1,k) - term_2(i,nyg,k)) / (yg(nyg)-yg(nyg-1))      ! dW/dy at upper wall

        SWij(i,3,k,1)     = (term  (i,3,k)     - term  (i,2,k)  ) / (yg(3)-yg(2))             ! dU/dy at lower wall+1
        SWij(i,3,k,2)     = (V_    (i,3,k)     - V_    (i,2,k)  ) / (y (3)-y (2))             ! dV/dy at lower wall+1 (1st order)
        SWij(i,3,k,3)     = (term_2(i,3,k)     - term_2(i,2,k)  ) / (yg(3)-yg(2))             ! dW/dy at lower wall+1

        SWij(i,nyg-2,k,1) = -(term  (i,nyg-2,k) - term  (i,nyg-1,k)) / (yg(nyg-1)-yg(nyg-2))  ! dU/dy at upper wall-1
        SWij(i,nyg-2,k,2) = -(V_    (i, ny-2,k) - V_    (i, ny-1,k)) / (y ( ny-1) -y ( ny-2)) ! dV/dy at upper wall-1 (1st order)
        SWij(i,nyg-2,k,3) = -(term_2(i,nyg-2,k) - term_2(i,nyg-1,k)) / (yg(nyg-1)-yg(nyg-2))  ! dW/dy at upper wall-1

      End Do
    End Do

    ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
    Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
    Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
    Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
       
    ! apply periodicity in x and z
    Call apply_periodic_bc_x(Uff,1)
    Call apply_periodic_bc_z(Uff,1)
    Call update_ghost_interior_planes(Uff,1)
    Call apply_periodic_bc_x(Vff,2)
    Call apply_periodic_bc_z(Vff,2)
    Call update_ghost_interior_planes(Vff,2)
    Call apply_periodic_bc_x(Wff,2)
    Call apply_periodic_bc_z(Wff,3)
    Call update_ghost_interior_planes(Wff,3)

    call interpolate_x(Uff,term  (2:nx,:,:),1) 
    call interpolate_z(Wff,term_2(:,:,2:nz),1) 

    Call apply_periodic_bc_x(term  ,2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    Do i = 2, nxg
      Do k = 2, nzg

        SWij(i,2,k,4)     = (term  (i,2,k)     - term  (i,1,k)  ) / (yg(2)-yg(1))            ! dU/dy at lower wall
        SWij(i,2,k,5)     = (Vff   (i,2,k)     - Vff   (i,1,k)  ) / (y (2)-y (1))            ! dV/dy at lower wall (1st order)
        SWij(i,2,k,6)     = (term_2(i,2,k)     - term_2(i,1,k)  ) / (yg(2)-yg(1))            ! dW/dy at lower wall

        SWij(i,nyg-1,k,4) = (term  (i,nyg,k)   - term  (i,nyg-1,k)) / (yg(nyg)-yg(nyg-1))    ! dU/dy at upper wall
        SWij(i,nyg-1,k,5) = (Vff   (i,ny ,k )  - Vff   (i,ny -1,k)) / (y ( ny)-y ( ny-1))    ! dV/dy at upper wall (1st order)
        SWij(i,nyg-1,k,6) = (term_2(i,nyg,k)   - term_2(i,nyg-1,k)) / (yg(nyg)-yg(nyg-1))    ! dW/dy at upper wall

      End Do
    End Do

    ! MWij = fil_size_2^2*hat(dui/dn)*hat(duj/dn) - dui/dn*duj/dn at wall
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! Diagonal elements
          Do j = 1, 3
             MWij(i,    2,k,j) = fil_size_2**2d0*SWij(i,    2,k,j+3)**2d0 - SWij(i,    2,k,j)**2d0 ! bottom
             MWij(i,nyg-1,k,j) = fil_size_2**2d0*SWij(i,nyg-1,k,j+3)**2d0 - SWij(i,nyg-1,k,j)**2d0 ! top
          End Do
          ! Off-diagonal elements
          MWij(i,    2,k,4) = fil_size_2**2d0*SWij(i,    2,k,4)*SWij(i,    2,k,5) - SWij(i,    2,k,1)*SWij(i,    2,k,2)
          MWij(i,nyg-1,k,4) = fil_size_2**2d0*SWij(i,nyg-1,k,4)*SWij(i,nyg-1,k,5) - SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,2)

          MWij(i,    2,k,5) = fil_size_2**2d0*SWij(i,    2,k,4)*SWij(i,    2,k,6) - SWij(i,    2,k,1)*SWij(i,    2,k,3)
          MWij(i,nyg-1,k,5) = fil_size_2**2d0*SWij(i,nyg-1,k,4)*SWij(i,nyg-1,k,6) - SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,3)

          MWij(i,    2,k,6) = fil_size_2**2d0*SWij(i,    2,k,5)*SWij(i,    2,k,6) - SWij(i,    2,k,2)*SWij(i,    2,k,3)
          MWij(i,nyg-1,k,6) = fil_size_2**2d0*SWij(i,nyg-1,k,5)*SWij(i,nyg-1,k,6) - SWij(i,nyg-1,k,2)*SWij(i,nyg-1,k,3)
       End Do
    End Do
    
    !-----------------------------------------------------------------------------!
    ! Part 4: Least squares solution term = LWij*MWij, term_1 = MWij*MWij

    term   = 0d0
    term_1 = 0d0
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! Diagonal elements
          Do j = 1, 3
             ! LWij*MWij
             term  (i, 1,k) = term  (i, 1,k) + LWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term  (i,ny,k) = term  (i,ny,k) + LWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
             ! MWij*MWij
             term_1(i, 1,k) = term_1(i, 1,k) + MWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term_1(i,ny,k) = term_1(i,ny,k) + MWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
          End Do
          ! Off-diagonal
          Do j = 4, 6
             ! LWij*MWij
             term  (i, 1,k) = term  (i, 1,k) + 2d0*LWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term  (i,ny,k) = term  (i,ny,k) + 2d0*LWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
             ! MWij*MWij
             term_1(i, 1,k) = term_1(i, 1,k) + 2d0*MWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term_1(i,ny,k) = term_1(i,ny,k) + 2d0*MWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
          End Do
       End Do
    End Do
    
    !-----------------------------------------------------------------------------!
    ! Part 5: Compute alpha
 
    ! Average over wall direction. Alpha^2 = <LWij*MWij> / <MWij*MWij>
    alpha_x = Sum( term  (2:nxg-1,1,2:nzg-1) + term  (2:nxg-1,ny,2:nzg-1) )
    alpha_y = Sum( term_1(2:nxg-1,1,2:nzg-1) + term_1(2:nxg-1,ny,2:nzg-1) )            

    If (myid == 0) Then
       Call MPI_Reduce(MPI_IN_PLACE,alpha_x(1,1,1),1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(MPI_IN_PLACE,alpha_y(1,1,1),1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    Else
       Call MPI_Reduce(alpha_x(1,1,1),0,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(alpha_y(1,1,1),0,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    End If

    Call Mpi_bcast (  alpha_x(1,1,1),1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  alpha_y(1,1,1),1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    alpha_x = alpha_x(1,1,1) / alpha_y(1,1,1)  

    ! Clipping if necessary
    alpha_x = Max( alpha_x,0d0 )
    alpha_x = alpha_x**0.5d0
    
    ! alpha_y
    alpha_y = alpha_x(2,1,2)
    
    ! alpha_z
    alpha_z = alpha_x(2,1,2) 
    
  End Subroutine compute_bose_model

  !--------------------------------------------------------------!
  !                                                              !
  !              Momentum balance at the wall                    !
  !                                                              !
  !     Compute slip lengths alpha_x, alpha_y and alpha_z        !
  !           for ui = alpha_i dui/dy at the wall                !
  !                                                              !
  ! alpha_x = alpha_y = alpha_z                                  !
  !                                                              !
  ! From x-momentum balance at the walls:                        !
  !       -2<UV+tau_UV> = dPdx*Ly - 2*nu*<dU/dy>                 !
  !    -> -2<UV> = dPdx*Ly - 2*<(nu+nu_t)*dU/dy>                 !
  !                                                              !
  ! Robin condition at the wall:                                 !
  !        U = alpha dU/dy, V = alpha dV/dy                      !
  !                                                              !
  ! Condition (everything is at the wall):                       !
  !      alpha^2 = (<nu*dU/dy>-dPdx*Ly/2-<tau_UV>)/<dU/dy*dV/dy> !
  !                                                              !
  ! Input:  U_,V_,W_,nu_t (flow velocities)                      !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)             !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_exact_Robin(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    ! local variables
    Real   (Int64) :: dU, dUV, mtau
    Real   (Int64) :: dU_local, dUV_local, mtau_local
    Integer(Int32) :: i, j, k

    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at cell edges in xy plane)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)

    ! Compute term_4 = tau_UV at cell edges
    ! term_2 -> nu_t
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1:nx,j,1:nzg)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    mtau_local = Sum( term_4(2:nx-1,1,2:nzg-1) - term_4(2:nx-1,ny,2:nzg-1) )
    Call MPI_Allreduce(mtau_local,mtau,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    mtau = mtau/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dU/dy*dV/dy> at the wall    

    ! interpolate V to cell edges
    Call interpolate_x(V_,term_1,1)

    ! product: dU/dy*dV/dy (first order)
    ! term_1 -> V
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k))/(yg(nyg) - yg(nyg-1))*(term_1(i,ny,k) - term_1(i,ny-1,k))/(y(ny) - y(ny-1))
       End Do
    End Do
    
    ! average <dU/dy*dV/dy>
    dUV_local = Sum( term(2:nx-1,1:2,2:nzg-1) )
    Call MPI_Allreduce(dUV_local,dUV,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUV = dUV/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall    

    ! product
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
       End Do
    End Do
    
    ! average <dU/dy>
    dU_local = Sum( term(2:nx-1,1:2,2:nzg-1) )
    Call MPI_Allreduce(dU_local,dU,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dU = dU/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 5: compute alpha_i

    alpha_x(:,:,:) = (nu*dU - dPdx_ref*Ly/2d0 - mtau)/dUV

    ! clipping 
    alpha_x = Max( alpha_x(2,1,2), 0d0 )
    alpha_x = ( alpha_x(2,1,2) )**0.5d0

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2)
    
  End Subroutine compute_exact_Robin

  !---------------------------------------------------!
  !       Set constant alpha_x=alpha_y=alpha_z        !
  !---------------------------------------------------!
  Subroutine compute_constant_alpha

    Real(Int64) :: omega_alpha
    
    omega_alpha = 2d0*pi*freq_mult/(y(2)-y(1))*dPdx**0.5d0

    alpha_x = alpha_mean_x*( 1d0 + alpha_std*dsin(omega_alpha*t) ) 
    alpha_y = alpha_mean_y*( 1d0 + alpha_std*dsin(omega_alpha*t) ) 
    alpha_z = alpha_mean_z*( 1d0 + alpha_std*dsin(omega_alpha*t) ) 

  end Subroutine compute_constant_alpha

  !--------------------------------------------------------------!
  !                                                              !
  !               Equilibrium wall model (Larsson)               !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_equilibrium_wall_model(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    ! local variables
    Integer(Int32) :: i, k, iter
    Real   (Int64) :: tau_eqwm, u_les, tau_eqwm_x, tau_eqwm_z, nu_t_ik

    ! interpolate W to U position 
    ! bottom
    temp_3(:,1,:) = W_(:,jref_eqwm,:)
    Call interpolate_z_2nd(z,  temp_3, zg, temp_1)
    Call interpolate_x_2nd(xg, temp_1,  x, temp_2)
    ! top
    temp_3(:,1,:) = W_(:,nyg-jref_eqwm+1,:)
    Call interpolate_z_2nd(z,  temp_3, zg, temp_1)
    Call interpolate_x_2nd(xg, temp_1,  x, temp_4)

    ! loop
    n_nan_eqwm  = 0
    n_iter_eqwm = 0
    Do i = 1,nx
       Do k = 1, nzg
          !--------------------Bottom wall
          ! input velocity
          u_les = ( U_(i,jref_eqwm,k)**2d0 + temp_2(i,1,k)**2d0 )**0.5d0 
          ! solve ODE
          Call solve_u_eqwm(u_les, nu, 1.0d0, yg(jref_eqwm), 100, 1.10d0, tau_eqwm, iter, tau_eqwm_initial(i,1,k), t)
          ! sanity check 
          If ( Isnan(tau_eqwm) ) Then
             n_nan_eqwm = n_nan_eqwm + 1 
             tau_eqwm   = tau_mean_eqwm
          end If
          If ( iter==100 ) Then
             tau_eqwm    = tau_mean_eqwm
             n_iter_eqwm = n_iter_eqwm + 1
          end If
          tau_eqwm_initial(i,1,k) = tau_eqwm
          ! distribute stress in x and z
          tau_eqwm_x = tau_eqwm*U_    (i,jref_eqwm,k)/u_les
          tau_eqwm_z = tau_eqwm*temp_2(i,        1,k)/u_les
          ! eddy viscosity
          nu_t_ik = 0.5d0*( 0.5d0*( nu_t_(i,1,k) + nu_t_(i+1,1,k) ) +  0.5d0*( nu_t_(i,2,k) + nu_t_(i+1,2,k) ) )
          ! constants for Neumann boundary conditions
          alpha_x     (i,1,k) = tau_eqwm_x/( nu + nu_t_ik )
          alpha_z_temp(i,1,k) = tau_eqwm_z/( nu + nu_t_ik )
          !--------------------top wall
          ! input velocity
          u_les = ( U_(i,nyg-jref_eqwm+1,k)**2d0 + temp_4(i,1,k)**2d0 )**0.5d0 
          ! solve ODE
          Call solve_u_eqwm(u_les, nu, 1.0d0, yg(jref_eqwm), 100, 1.10d0, tau_eqwm, iter, tau_eqwm_initial(i,2,k), t)
          ! sanity check 
          If ( Isnan(tau_eqwm) ) Then
             n_nan_eqwm = n_nan_eqwm + 1 
             tau_eqwm   = tau_mean_eqwm
          end If
          If ( iter==100 ) Then
             tau_eqwm    = tau_mean_eqwm
             n_iter_eqwm = n_iter_eqwm + 1
          end If
          tau_eqwm_initial(i,2,k) = tau_eqwm
          ! distribute stress in x and z
          tau_eqwm_x = tau_eqwm*U_    (i,nyg-jref_eqwm+1,k)/u_les
          tau_eqwm_z = tau_eqwm*temp_4(i,              1,k)/u_les
          ! eddy viscosity
          nu_t_ik = 0.5d0*( 0.5d0*( nu_t_(i,nyg,k) + nu_t_(i+1,nyg,k) ) +  0.5d0*( nu_t_(i,nyg-1,k) + nu_t_(i+1,nyg-1,k) ) )
          ! constants for Neumann boundary conditions
          alpha_x     (i,2,k) = tau_eqwm_x/( nu + nu_t_ik )
          alpha_z_temp(i,2,k) = tau_eqwm_z/( nu + nu_t_ik )
       End do
    End do
    
    ! interpolate boundary condition for W
    Do i = 1, nxg-2
       Do k = 1, nz-1
          alpha_z(i,:,k) = 0.25d0*( alpha_z_temp(i,:,k) + alpha_z_temp(i+1,:,k) + alpha_z_temp(i,:,k+1) + alpha_z_temp(i+1,:,k+1) )
       End Do
    End Do
    alpha_z(nxg  ,:, :) = alpha_z(3,:,:)
    alpha_z(nxg-1,:, :) = alpha_z(2,:,:)
    alpha_z( :   ,:,nz) = alpha_z(:,:,2)

    ! check NaNs
    If ( Mod(istep,nmonitor)==0 ) Then
       If (myid == 0) Then
          Call MPI_Reduce(MPI_IN_PLACE,n_nan_eqwm,1,MPI_integer,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Else
          Call MPI_Reduce(n_nan_eqwm,0,1,MPI_integer,MPI_sum,0,MPI_COMM_WORLD,ierr)
       End If
    End If
    
  End Subroutine compute_equilibrium_wall_model

  !--------------------------------------------------------------!
  !                                                              !
  !          Equilibrium wall model (Larsson) with slip          !
  !                                                              !
  !  Compute wall stress from Eq. WM and apply slip boundary     !
  !  condition                                                   !
  !                                                              !
  !--------------------------------------------------------------!
  Subroutine compute_exact_Neumann(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    ! local variables
    Integer(Int32) :: i, k, iter
    Real   (Int64) :: tau_eqwm, u_les, tau_eqwm_x, tau_eqwm_z, nu_t_ik

    ! interpolate W to U position 
    ! bottom
    temp_3(:,1,:) = W_(:,jref_eqwm,:)
    Call interpolate_z_2nd(z,  temp_3, zg, temp_1)
    Call interpolate_x_2nd(xg, temp_1,  x, temp_2)
    ! top
    temp_3(:,1,:) = W_(:,nyg-jref_eqwm+1,:)
    Call interpolate_z_2nd(z,  temp_3, zg, temp_1)
    Call interpolate_x_2nd(xg, temp_1,  x, temp_4)
    
    nu_t_ik = Sum(nu_t_(2:nxg-1,1:2,2:nzg-1))
    If (myid == 0) Then
      Call MPI_Reduce(MPI_IN_PLACE,nu_t_ik,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    Else
      Call MPI_Reduce(nu_t_ik,0,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    End If
    Call Mpi_bcast ( nu_t_ik,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    nu_t_ik = nu_t_ik / Real((nxg-2)*(nzg_global-2)*2,8)

    alpha_x = dPdx_ref / (nu+nu_t_ik)
   
    ! interpolate boundary condition for W
    !Do i = 1, nxg-2
    !   Do k = 1, nz-1
    !      alpha_z(i,:,k) = 0.25d0*( alpha_z_temp(i,:,k) + alpha_z_temp(i+1,:,k) + alpha_z_temp(i,:,k+1) + alpha_z_temp(i+1,:,k+1) )
    !   End Do
    !End Do
    !alpha_z(nxg  ,:, :) = alpha_z(3,:,:)
    !alpha_z(nxg-1,:, :) = alpha_z(2,:,:)
    !alpha_z( :   ,:,nz) = alpha_z(:,:,2)
    alpha_z = 0d0

  End Subroutine compute_exact_Neumann

  !----------------------------------------------------------------------!
  !                                                                      !
  !              Dynamic model using momentum equations                  !
  !                                                                      !
  !          Compute slip lengths alpha_x, alpha_y and alpha_z           !
  !                  for ui = alpha_i dui/dy at the wall                 !
  !                                                                      !
  ! Equation:                                                            !
  !                                                                      !
  !  alpha^2 = <LiWj*MWij>/<MWij*MWij>                                   !
  !  In case nodiag == 1 : Do not include diagonal terms                 !
  !                                                                      !
  !    MWij = dui/dn*duj/dn                                              !
  !    LWij = ui*uj + a1 * ( ui*uj                     + tauij         ) ! 
  !                 + a2 * ( hat(ui*uj)                + hat(tauij)    ) !
  !                 + a3 * ( hat(hat(ui*uj))           + hat(hat(tauij)) !
  !                 + a4 * ( hat(ui)*hat(uj)           + Tij           ) !
  !                 + a5 * ( hat(hat(ui)*hat(uj))      + hat(Tij)      ) !
  !                 + a6 * ( hat(hat(ui))*hat(hat(uj)) + TTij          ) !
  !                                                                      !
  !  n->normal to the wall                                               !
  !                                                                      !
  ! Definitions:                                                         !
  !                                                                      !
  ! alpha = alpha_x = alpha_y = alpha_z                                  !
  ! alpha computed at V positions                                        !
  !                                                                      !
  ! Input:  U_,V_,W_,nu_t (flow velocities)                              !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)                     !
  !                                                                      !
  !----------------------------------------------------------------------!
  Subroutine compute_generic_dynamic_model(U_,V_,W_,nu_t_,cu,c_test,c_diag,ai)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    Integer(Int32), Intent(In) :: cu, c_test, c_diag
    Integer(Int32), Dimension(6), Intent(In) :: ai

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: alpha_x_1, alpha_x_2, alpha_y_1, alpha_y_2

    !--------------------------------------------------------------------------!
    ! Part 1: Compute LWij and MWij
    Call compute_LWij_MWij(U_,V_,W_,nu_t_,cu,c_test,ai)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   If (myid==0) Then
!      Write(*,*) 'Luu',LWij(2:4,nyg-1,2,1)
!      Write(*,*) 'Lvv',LWij(2:4,nyg-1,2,2)
!      Write(*,*) 'Lww',LWij(2:4,nyg-1,2,3)
!      Write(*,*) 'Luv',LWij(2:4,nyg-1,2,4)
!      Write(*,*) 'Luw',LWij(2:4,nyg-1,2,5)
!      Write(*,*) 'Lvw',LWij(2:4,nyg-1,2,6)
!      Write(*,*) '------------'
!      Write(*,*) 'Muu',MWij(2:4,nyg-1,2,1)
!      Write(*,*) 'Mvv',MWij(2:4,nyg-1,2,2)
!      Write(*,*) 'Mww',MWij(2:4,nyg-1,2,3)
!      Write(*,*) 'Muv',MWij(2:4,nyg-1,2,4)
!      Write(*,*) 'Muw',MWij(2:4,nyg-1,2,5)
!      Write(*,*) 'Mvw',MWij(2:4,nyg-1,2,6)
!   end If
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    !-----------------------------------------------------------------------------!
    ! Part 2: Least squares solution term = LWij*MWij, term_1 = MWij*MWij
    term   = 0d0
    term_1 = 0d0
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! Diagonal elements
          If ( c_diag == 1 ) Then 
             Do j = 1, 3
                ! LWij*MWij
                term  (i, 1,k) = term  (i, 1,k) + 1d0*LWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
                term  (i,ny,k) = term  (i,ny,k) + 1d0*LWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
                ! MWij*MWij
                term_1(i, 1,k) = term_1(i, 1,k) + 1d0*MWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
                term_1(i,ny,k) = term_1(i,ny,k) + 1d0*MWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
             End Do
          End If
          ! Off-diagonal x2
          Do j = 4, 6
             ! LWij*MWij
             term  (i, 1,k) = term  (i, 1,k) + 2d0*LWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term  (i,ny,k) = term  (i,ny,k) + 2d0*LWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
             ! MWij*MWij
             term_1(i, 1,k) = term_1(i, 1,k) + 2d0*MWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term_1(i,ny,k) = term_1(i,ny,k) + 2d0*MWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
          End Do
       End Do
    End Do
    
    !-----------------------------------------------------------------------------!
    ! Part 3: Compute alpha
    alpha_x_1 = Sum( term  (2:nxg-1, 1,2:nzg-1) ) 
    alpha_x_2 = Sum( term  (2:nxg-1,ny,2:nzg-1) )
    alpha_y_1 = Sum( term_1(2:nxg-1, 1,2:nzg-1) ) 
    alpha_y_2 = Sum( term_1(2:nxg-1,ny,2:nzg-1) ) 
    
    alpha_x   = alpha_x_1 + alpha_x_2
    alpha_y   = alpha_y_1 + alpha_y_2 

    If (myid == 0) Then
      Call MPI_Reduce(MPI_IN_PLACE,alpha_x(1,1,1),1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
      Call MPI_Reduce(MPI_IN_PLACE,alpha_y(1,1,1),1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    Else
      Call MPI_Reduce(alpha_x(1,1,1),0,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
      Call MPI_Reduce(alpha_y(1,1,1),0,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
    End If

    Call Mpi_bcast (  alpha_x(1,1,1),1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  alpha_y(1,1,1),1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    alpha_x = alpha_x(1,1,1) / alpha_y(1,1,1)  

    alpha_x = Max( alpha_x,0d0 )
    alpha_x = alpha_x**0.5d0

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2) 

  End Subroutine compute_generic_dynamic_model
 
  !----------------------------------------------------------------------!
  !                                                                      !
  !              Dynamic model using momentum equations                  !
  !                                                                      !
  !          Compute slip lengths alpha_x, alpha_y and alpha_z           !
  !                  for ui = alpha_i dui/dy at the wall                 !
  !                                                                      !
  ! Equation:                                                            !
  !                                                                      !
  ! Input:  U_,V_,W_,nu_t (flow velocities)                              !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)                     !
  !                                                                      !
  !----------------------------------------------------------------------!
  Subroutine compute_generic_dynamic_model_2(U_,V_,W_,nu_t_,c_diag,c_case,ai)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    Integer(Int32), Intent(In) :: c_case, c_diag
    Integer(Int32), Dimension(6), Intent(In) :: ai

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dU, dUV, mtau
    Real   (Int64) :: dU_local, dUV_local, mtau_local
    Real   (Int64) :: LWijMWij_local, MWijMWij_local, LWijMWij, MWijMWij
    Real   (Int64) :: tauw_local, tauw
    Real   (Int64), Dimension(5) :: Cw

    !--------------------------------------------------------------------------!
    ! Part 1: Compute LWij and MWij
    Call compute_LWij_MWij_CW(U_,V_,W_,nu_t_,c_case,ai)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   If (myid==0) Then
!      Write(*,*) 'Luu',LWij(2:4,2,2,1)
!      Write(*,*) 'Lvv',LWij(2:4,2,2,2)
!      Write(*,*) 'Lww',LWij(2:4,2,2,3)
!      Write(*,*) 'Luv',LWij(2:4,2,2,4)
!      Write(*,*) 'Luw',LWij(2:4,2,2,5)
!      Write(*,*) 'Lvw',LWij(2:4,2,2,6)
!      Write(*,*) '------------'
!      Write(*,*) 'Muu',MWij(2:4,2,2,1)
!      Write(*,*) 'Mvv',MWij(2:4,2,2,2)
!      Write(*,*) 'Mww',MWij(2:4,2,2,3)
!      Write(*,*) 'Muv',MWij(2:4,2,2,4)
!      Write(*,*) 'Muw',MWij(2:4,2,2,5)
!      Write(*,*) 'Mvw',MWij(2:4,2,2,6)
!   end If
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    !-----------------------------------------------------------------------------!
    ! Part 2: Least squares solution term = LWij*MWij, term_1 = MWij*MWij
    term   = 0d0
    term_1 = 0d0
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          ! Diagonal elements
          If ( c_diag == 1 ) Then 
             Do j = 1, 3
                ! LWij*MWij
                term  (i, 1,k) = term  (i, 1,k) + 1d0*LWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
                term  (i,ny,k) = term  (i,ny,k) + 1d0*LWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
                ! MWij*MWij
                term_1(i, 1,k) = term_1(i, 1,k) + 1d0*MWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
                term_1(i,ny,k) = term_1(i,ny,k) + 1d0*MWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
             End Do
          End If
          ! Off-diagonal x2
          Do j = 4, 6
             ! LWij*MWij
             term  (i, 1,k) = term  (i, 1,k) + 2d0*LWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term  (i,ny,k) = term  (i,ny,k) + 2d0*LWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
             ! MWij*MWij
             term_1(i, 1,k) = term_1(i, 1,k) + 2d0*MWij(i,    2,k,j)*MWij(i,    2,k,j) ! bottom
             term_1(i,ny,k) = term_1(i,ny,k) + 2d0*MWij(i,nyg-1,k,j)*MWij(i,nyg-1,k,j) ! top
          End Do
       End Do
    End Do
    
    !-----------------------------------------------------------------------------!
    ! Part 3: Compute alpha
    
    LWijMWij_local = Sum( term  (2:nxg-1, 1,2:nzg-1) ) !+ Sum( term  (2:nxg-1,ny,2:nzg-1) )
    MWijMWij_local = Sum( term_1(2:nxg-1, 1,2:nzg-1) ) !+ Sum( term_1(2:nxg-1,ny,2:nzg-1) )

    Call MPI_Allreduce(LWijMWij_local,LWijMWij,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_Allreduce(MWijMWij_local,MWijMWij,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)

    Cw = LWijMWij / MWijMWij

    Call compute_delta_tau(U_,V_,W_,nu_t_,Cw,ai)

    tauw_local = Sum( LWij(2:nxg-1, 2,2:nzg-1, 4) ) !+ Sum( LWij(2:nxg-1,nyg-1,2:nzg-1, 4) )
    Call MPI_Allreduce(tauw_local,tauw,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    tauw = tauw / Real((nxg_global-2)*(nzg_global-2),8)


    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at cell edges in xy plane)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)

    ! Compute term_4 = tau_UV at cell edges
    ! term_2 -> nu_t
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1:nx,j,1:nzg)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    mtau_local = Sum( term_4(2:nx-1,1,2:nzg-1) - term_4(2:nx-1,ny,2:nzg-1) )
    Call MPI_Allreduce(mtau_local,mtau,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    mtau = mtau/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dU/dy*dV/dy> at the wall    

    ! interpolate V to cell edges
    Call interpolate_x(V_,term_1,1)

    ! product: dU/dy*dV/dy (first order)
    ! term_1 -> V
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k))/(yg(nyg) - yg(nyg-1))*(term_1(i,ny,k) - term_1(i,ny-1,k))/(y(ny) - y(ny-1))
       End Do
    End Do
    
    ! average <dU/dy*dV/dy>
    dUV_local = Sum( term(2:nx-1,1:2,2:nzg-1) )
    Call MPI_Allreduce(dUV_local,dUV,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUV = dUV/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall    

    ! product
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
       End Do
    End Do
    
    ! average <dU/dy>
    dU_local = Sum( term(2:nx-1,1:2,2:nzg-1) )
    Call MPI_Allreduce(dU_local,dU,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dU = dU/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 5: compute alpha_i

    alpha_x(:,:,:) = (nu*dU - tauw*Ly/2d0 - mtau)/dUV

    ! clipping 
    alpha_x = Max( alpha_x(2,1,2), 0d0 )
    alpha_x = ( alpha_x(2,1,2) )**0.5d0

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2)

  End Subroutine compute_generic_dynamic_model_2

  !----------------------------------------------------------------------!
  !                                                                      !
  !              Dynamic model using momentum equations                  !
  !                                                                      !
  !          Compute slip lengths alpha_x, alpha_y and alpha_z           !
  !                  for ui = alpha_i dui/dy at the wall                 !
  !                                                                      !
  ! Equation:                                                            !
  !                                                                      !
  ! Input:  U_,V_,W_,nu_t (flow velocities)                              !
  ! Output: alpha_x, alpha_y, alpha_z (slip lengths)                     !
  !                                                                      !
  !----------------------------------------------------------------------!
  Subroutine compute_generic_dynamic_model_3(U_,V_,W_,nu_t_,c_diag,c_case,ai)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_

    Integer(Int32), Intent(In) :: c_case, c_diag
    Integer(Int32), Dimension(6), Intent(In) :: ai

    ! local variables
    Integer(Int32) :: ii, jj, i, j, k
    Real   (Int64) :: dU, dUV, mtau
    Real   (Int64) :: dU_local, dUV_local, mtau_local
    Real   (Int64) :: LWijMWij_local, MWijMWij_local, LWijMWij, MWijMWij 
    Real   (Int64) :: tauw_local, tauw
    Integer(Int32) :: dim_temp

    Real   (Int64), Dimension(5)   :: Cw, LM
    Real   (Int64), Dimension(5,5) :: MM
    Integer(Int32), Dimension(6)   :: ai_temp

    Real   (Int64), Dimension(2:nxg-1,2:nyg-1,2:nzm+1, 6, 5) :: MWij_temp

    ! local variables for LAPACK solver
    Integer(Int32) :: INFO
    Integer(Int32), Dimension(5) :: IPIV 

    ai_temp(1) = 0
    ai_temp(2) = 0
    ai_temp(3) = 0
    ai_temp(4) = 0
    ai_temp(5) = 0

    dim_temp   = 0
    !--------------------------------------------------------------------------!
    Do ii = 1,5
       If ( ai(ii)/=0 ) Then
          dim_temp   = dim_temp+1
          ai_temp(ii) = ai(ii)
          ! Part 1: Compute LWij and MWij
          Call compute_LWij_MWij_CW(U_,V_,W_,nu_t_,c_case,ai_temp)
          MWij_temp(:,:,:,:,dim_temp) = MWij
          ai_temp(ii) = 0
       End If
    End Do

    Do ii = 1,dim_temp
          !-----------------------------------------------------------------------------!
          ! Part 2: Least squares solution term = LWij*MWij, term_1 = MWij*MWij
          term   = 0d0
          term_1 = 0d0
          Do i = 2, nxg-1
             Do k = 2, nzg-1
                ! Diagonal elements
                If ( c_diag == 1 ) Then 
                   Do j = 1, 3
                      ! LWij*MWij
                      term  (i, 1,k) = term  (i, 1,k) + 1d0*LWij(i,    2,k,j)*MWij_temp(i,    2,k,j,ii) ! bottom
                      term  (i,ny,k) = term  (i,ny,k) + 1d0*LWij(i,nyg-1,k,j)*MWij_temp(i,nyg-1,k,j,ii) ! top
                   End Do
                End If
                ! Off-diagonal x2
                Do j = 4, 6
                   ! LWij*MWij
                   term  (i, 1,k) = term  (i, 1,k) + 2d0*LWij(i,    2,k,j)*MWij_temp(i,    2,k,j,ii) ! bottom
                   term  (i,ny,k) = term  (i,ny,k) + 2d0*LWij(i,nyg-1,k,j)*MWij_temp(i,nyg-1,k,j,ii) ! top
                End Do
             End Do
          End Do
          LWijMWij_local = Sum( term  (2:nxg-1, 1,2:nzg-1) ) !+ Sum( term  (2:nxg-1,ny,2:nzg-1) )
          Call MPI_Allreduce(LWijMWij_local,LWijMWij,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)

          LM(ii) = LWijMWij

       Do jj = 1,dim_temp
          !-----------------------------------------------------------------------------!
          ! Part 2: Least squares solution term = LWij*MWij, term_1 = MWij*MWij
          term   = 0d0
          term_1 = 0d0
          Do i = 2, nxg-1
             Do k = 2, nzg-1
                ! Diagonal elements
                If ( c_diag == 1 ) Then 
                   Do j = 1, 3
                      ! MWij*MWij
                      term_1(i, 1,k) = term_1(i, 1,k) + 1d0*MWij_temp(i,    2,k,j,ii)*MWij_temp(i,    2,k,j,jj) ! bottom
                      term_1(i,ny,k) = term_1(i,ny,k) + 1d0*MWij_temp(i,nyg-1,k,j,ii)*MWij_temp(i,nyg-1,k,j,jj) ! top
                   End Do
                End If
                ! Off-diagonal x2
                Do j = 4, 6
                   ! MWij*MWij
                   term_1(i, 1,k) = term_1(i, 1,k) + 2d0*MWij_temp(i,    2,k,j,ii)*MWij_temp(i,    2,k,j,jj) ! bottom
                   term_1(i,ny,k) = term_1(i,ny,k) + 2d0*MWij_temp(i,nyg-1,k,j,ii)*MWij_temp(i,nyg-1,k,j,jj) ! top
                End Do
             End Do
          End Do
    
          MWijMWij_local = Sum( term_1(2:nxg-1, 1,2:nzg-1) ) !+ Sum( term_1(2:nxg-1,ny,2:nzg-1) )
          Call MPI_Allreduce(MWijMWij_local,MWijMWij,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
          MM(ii,jj) = MWijMWij
       End Do
    End Do

!If (myid == 0) write(*,*) 'MM', MM(1:dim_temp,1:dim_temp)
!If (myid == 0) write(*,*) 'LM', LM(1:dim_temp)


    If (dim_temp > 1) Then
       Call DGESV(dim_temp, 1, MM(1:dim_temp,1:dim_temp), dim_temp, IPIV, LM(1:dim_temp), dim_temp, INFO) 
    Else
       LM(1) = LM(1) / MM(1,1)
    End If

    dim_temp = 0
    Do ii = 1,5
       If ( ai(ii)/=0 ) Then
          dim_temp   = dim_temp+1
          Cw(ii) = LM(dim_temp) 
       End If
    End Do

!If (myid == 1) Write(*,*) Cw

    Call compute_delta_tau(U_,V_,W_,nu_t_,Cw,ai)

    tauw_local = Sum( LWij(2:nxg-1, 2,2:nzg-1, 4) ) !+ Sum( LWij(2:nxg-1,nyg-1,2:nzg-1, 4) )
    Call MPI_Allreduce(tauw_local,tauw,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    tauw = tauw / Real((nxg_global-2)*(nzg_global-2),8)

!If (myid == 1) Write(*,*) tauw 

    !---------------------------------------------------------!
    ! Part 1: compute tau_UV -> term_4 (at cell edges in xy plane)
   
    ! compute Suv at cell edges -> term_3(1:nx,1:ny,1:nzg)
    Do i = 1, nx
       Do j = 1, ny
          Do k = 1, nzg
             ! 1/2*(dU/dy + dV/dx)
             term_3(i,j,k) = 0.5d0*( (U_(i,  j+1,k) - U_(i,j,k))/( yg(j+1)-yg(j) ) + &
                                     (V_(i+1,j  ,k) - V_(i,j,k))/( xg(i+1)-xg(i) ) ) 
          End Do
       End Do
    End Do
    
    ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)

    ! Compute term_4 = tau_UV at cell edges
    ! term_2 -> nu_t
    ! term_3 -> Suv
    Do j = 1, ny
       term_4(1:nx,j,1:nzg) = -2d0*term_2(1:nx,j,1:nzg)*term_3(1:nx,j,1:nzg)
    End Do
    
    ! average tau_UV
    mtau_local = Sum( term_4(2:nx-1,1,2:nzg-1) - term_4(2:nx-1,ny,2:nzg-1) )
    Call MPI_Allreduce(mtau_local,mtau,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    mtau = mtau/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 2: compute <dU/dy*dV/dy> at the wall    

    ! interpolate V to cell edges
    Call interpolate_x(V_,term_1,1)

    ! product: dU/dy*dV/dy (first order)
    ! term_1 -> V
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k))/(yg(2) - yg(1))*(term_1(i,2,k) - term_1(i,1,k))/(y(2) - y(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k))/(yg(nyg) - yg(nyg-1))*(term_1(i,ny,k) - term_1(i,ny-1,k))/(y(ny) - y(ny-1))
       End Do
    End Do
    
    ! average <dU/dy*dV/dy>
    dUV_local = Sum( term(2:nx-1,1:2,2:nzg-1) )
    Call MPI_Allreduce(dUV_local,dUV,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dUV = dUV/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 3: compute <dU/dy> at the wall    

    ! product
    term(:,1:2,:) = 0d0
    Do i = 2, nx-1
       Do k = 2, nzg-1
          ! bottom wall
          term(i,1,k) =  (U_(i,2,k) - U_(i,1,k)) / (yg(2) - yg(1))
          ! top wall (sign changed for average)
          term(i,2,k) = -(U_(i,nyg,k) - U_(i,nyg-1,k)) / (yg(nyg) - yg(nyg-1))
       End Do
    End Do
    
    ! average <dU/dy>
    dU_local = Sum( term(2:nx-1,1:2,2:nzg-1) )
    Call MPI_Allreduce(dU_local,dU,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    dU = dU/Real( 2*(nx_global-2)*(nzg_global-2), 8)

    !---------------------------------------------------------!
    ! Part 5: compute alpha_i

    alpha_x(:,:,:) = (nu*dU - tauw*Ly/2d0 - mtau)/dUV

    ! clipping 
    alpha_x = Max( alpha_x(2,1,2), 0d0 )
    alpha_x = ( alpha_x(2,1,2) )**0.5d0

    ! alpha_y
    alpha_y = alpha_x(2,1,2)

    ! alpha_z
    alpha_z = alpha_x(2,1,2)

  End Subroutine compute_generic_dynamic_model_3

  !-------------------------------------------------------------!
  !                                                             !
  !            Compute Neumann boundary conditions              !
  !    for pseudo-pressure when slip-wall model is active       !
  !                                                             !
  ! This has to be called every sub-step                        !
  !                                                             !
  ! Conditions bottom wall:                                     !
  !                                                             !
  !     V1 = V1*  -  (p2-p1)/(yg(2)-yg(1))                      !
  !     V2 = V2*  -  (p3-p2)/(yg(3)-yg(2))                      !
  !     V1 = alpha_y*(V2-V1)/(y (2)-y (1))                      !
  !                                                             !
  !     => p1 = p_b2*p2 + p_bc3*p3                              !
  !        p_bc2   = 1 + beta*Delta_r                           !
  !        p_bc3   =   - beta*Delta_r                           !
  !        Delta_r = ( yg(2)-yg(1) )/( yg(3)-yg(2) )            !
  !        beta    = alpha_y/(alpha_y + y(2)-y(1) )             !
  !                                                             !
  !     V* -> velocity without pressure                         !
  !                                                             !
  ! Equation for first interior points:                         !
  !                                                             !
  !    (a + c*p_bc3)*p3 + (b + c*p_bc2)*p2 = rhs_p2             !
  !                                                             !
  !                                                             !
  ! Conditions top wall: (n->ny, ng->nyg)                       !
  !                                                             !
  !     V(n)   = V(  n)* -(p(ng)  -p(ng-1))/(yg(ng)  -yg(ng-1)) !
  !     V(n-1) = V(n-1)* -(p(ng-1)-p(ng-2))/(yg(ng-1)-yg(ng-2)) !
  !     V(n)   = alpha_y*(V(n)-V(n-1))/(y(n)-y(n-1))            !
  !                                                             !
  !     => p(ng)    = p_bn1*p(ng-2) + p_bcn*p(ng-1)             !
  !        p_bcn    = 1d0 + beta*Delta_r                        !
  !        p_bcn1   =     - beta*Delta_r                        !
  !        alpha_y  = - alpha_y -> same alpha_y each wall!      !
  !        Delta_r  = (yg(ng)-yg(ng-1))/(yg(ng-1)-yg(ng-2))     !
  !        beta     = alpha_y/( alpha_y - (y(n)-y(n-1)) )       !
  !                                                             !
  ! Equation for last interior points:                          !
  !                                                             !
  !   (b+a*p_bcn)*p(ng-1) + (c+a*p_bcn1)*p(ng-2) = rhs_p2(ng-1) !
  !                                                             !
  !                                                             !
  ! Assumed:                                                    !
  !      - alpha_y is not function of (x,z)                     !
  !      - V1* = alpha_y*(V2*-V1*)/(y(2)-y(1))                  !
  !        Otherwise the rhs for P must be modified             !
  !                                                             !
  !-------------------------------------------------------------!
  Subroutine compute_pseudo_pressure_bc_for_robin_bc

    ! local variables
    Real   (Int64) :: a, b, c
    Real   (Int64) :: beta, Delta_r, alphad
    Real   (Int64) :: p_bc2, p_bc3, p_bcn, p_bcn1
    Integer(Int64) :: j    

    ! bottom wall
    j        = 2 
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    Delta_r  = ( yg(2)-yg(1) )/( yg(3)-yg(2) )
    alphad   = alpha_y(2,1,2)
    beta     = alphad/(alphad + y(2)-y(1) )
    p_bc2    = 1d0 + beta*Delta_r
    p_bc3    =     - beta*Delta_r

    Dyy(2,2) = b + c*p_bc2
    Dyy(2,3) = a + c*p_bc3
    
    ! top wall
    j        = nyg-1
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    alphad   = -alpha_y(2,2,2)
    Delta_r  = ( yg(nyg) - yg(nyg-1) )/( yg(nyg-1) - yg(nyg-2) )
    beta     = alphad/( alphad - (y(ny)-y(ny-1)) )
    p_bcn    = 1d0 + beta*Delta_r
    p_bcn1   =     - beta*Delta_r
    
    Dyy(nyg-1,nyg-1) = b + a*p_bcn
    Dyy(nyg-1,nyg-2) = c + a*p_bcn1
    
  End Subroutine compute_pseudo_pressure_bc_for_robin_bc
      
  !------------------------------------------------------------------------------!
  !                                                                              !
  !                     Compute estimate of \Delta \tau_ij                       !
  !                     for the first two cell at the wall                       !
  !                                                                              !
  ! Assumed uniform mesh in y                                                    !
  !                                                                              !
  ! Input:  U, V, W, (flow velocities)                                           !
  ! Output: SWij(2:nxg,2:3,2:nzg) (bottom)                                       !
  !         SWij(2:nxg,nyg-1:nyg,2:nzg) (top)                                    !
  !                                                                              !
  !------------------------------------------------------------------------------!
  Subroutine compute_Tij_estimate(U_,V_,W_,SWij_)
    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_    

    Real(Int64), Dimension(2:nxg,2:nyg,2:nzg,6), Intent(Out) :: SWij_ 

    ! local variables
    Integer(Int32) :: i, j, k, jv
    Integer(Int32) :: jindex(4), jj
       ! indices for diagonal part
       jindex(1) = 2
       jindex(2) = 3
       jindex(3) = nyg-1
       jindex(4) = nyg

       ! Compute diagonal SWij
       Do k = 2, nzg-1
          Do jj = 1, 4
             j  = jindex(jj)
             Do i = 2, nxg-1
                ! These are at cell center
                ! dU/dx
                ten_buf_2(i,j,k,1) = 0.5d0*( (U_(i,j,k) - U_(i-1,j,k))/(x(i)-x(i-1)) + (U_(i,j-1,k) - U_(i-1,j-1,k))/(x(i)-x(i-1)) )
                ! dV/dy first order
                jv = j
                If ( j>nyg/2 ) jv = j-1
                ten_buf_2(i,j,k,2) = (V_(i,jv,k) - V_(i,jv-1,k))/(y(jv)-y(jv-1)) 
                ! dW/dz
                ten_buf_2(i,j,k,3) = 0.5d0*( (W_(i,j,k) - W_(i,j,k-1))/(z(k)-z(k-1)) + (W_(i,j-1,k) - W_(i,j-1,k-1))/(z(k)-z(k-1)) )
             End Do
         End Do
       End Do

       ! Compute off-diagonal SWij
       Do k = 2, nzg-1
          Do jj = 1, 4
             j  = jindex(jj)
             jv = j
             If ( j>nyg/2 ) jv = j-1
             Do i = 2, nxg-1
                ! These are at cell edges
                ! dU/dy 
                                     !(i,j)
                ten_buf_2(i,j,k,4) = 0.5d0*( (U_(i-1,j,k)-U_(i-1,j-1,k))/(yg(j)-yg(j-1)) + & 
                                     !(i+1,j)
                                             (U_(i-1+1,j,k)-U_(i-1+1,j-1,k))/(yg(j)-yg(j-1))  )   
                ! dU/dz 
                ten_buf_2(i,j,k,5) = 0.25d0*( &
                     !(i,k)
                     (U_(i-1,j,k)-U_(i-1,j,k-1))/(zg(k)-zg(k-1)) + &  
                     !(i+1,k)
                     (U_(i-1+1,j,k)-U_(i-1+1,j,k-1))/(zg(k)-zg(k-1)) + &  
                     !(i,k+1)
                     (U_(i-1,j,k+1)-U_(i-1,j,k-1+1))/(zg(k+1)-zg(k-1+1)) + &  
                     !(i+1,k+1)
                     (U_(i-1+1,j,k+1)-U_(i-1+1,j,k-1+1))/(zg(k+1)-zg(k-1+1))  )    

                ! dV/dz 
                                        !(j,k)
                ten_buf_2(i,j,k,6) = 0.5d0*( (V_(i,jv-1,k)-V_(i,jv-1,k-1))/(zg(k)-zg(k-1)) + & 
                                        !(j,k+1)
                                             (V_(i,jv-1,k+1)-V_(i,jv-1,k-1+1))/(zg(k+1)-zg(k-1+1)) )

                ! dV/dx 
                ten_buf_2(i,j,k,7) = 0.5d0*( (V_(i,jv-1,k)-V_(i-1,jv-1,k))/(xg(i)-xg(i-1)) + & 
                                     !(i+1,j)
                                             (V_(i+1,jv-1,k)-V_(i-1+1,jv-1,k))/(xg(i+1)-xg(i-1+1)) )   
                ! dW/dx
                ten_buf_2(i,j,k,8) = 0.25d0*( &
                     !(i,k)
                     (W_(i,j,k-1)-W_(i-1,j,k-1))/(xg(i)-xg(i-1)) + &  
                     !(i+1,k)
                     (W_(i+1,j,k-1)-W_(i-1+1,j,k-1))/(xg(i+1)-xg(i-1+1)) + &  
                     !(i,k+1)
                     (W_(i,j,k-1+1)-W_(i-1,j,k-1+1))/(xg(i)-xg(i-1)) + &  
                     !(i+1,k+1)
                     (W_(i+1,j,k-1+1)-W_(i-1+1,j,k-1+1))/(xg(i+1)-xg(i-1+1)) )    

                ! dW/dy 
                                        !(j,k)
                ten_buf_2(i,j,k,9) = 0.5d0*( (W_(i,j,k-1)-W_(i,j-1,k-1))/(yg(j)-yg(j-1)) + & 
                                        !(j,k+1)
                                             (W_(i,j,k-1+1)-W_(i,j-1,k-1+1))/(yg(j)-yg(j-1)) ) 
             End Do
          End Do
       End Do
    Do i = 1,9
       Call apply_periodic_bc_x(ten_buf_2(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf_2(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf_2(:,:,:,i),4)
    End Do

    ! dU/dx*dU/dx + dU/dy*dU/dy + dU/dz*dU/dz
    SWij_(2:nxg,2:3      ,2:nzg,1) = ten_buf_2(2:nxg,2:3      ,2:nzg,1)**2d0 + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,4)**2d0 + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,5)**2d0
    SWij_(2:nxg,nyg-1:nyg,2:nzg,1) = ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,1)**2d0 + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,4)**2d0 + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,5)**2d0
    ! dV/dx*dV/dx + dV/dy*dV/dy + dV/dz*dV/dz
    SWij_(2:nxg,2:3      ,2:nzg,2) = ten_buf_2(2:nxg,2:3      ,2:nzg,2)**2d0 + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,6)**2d0 + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,7)**2d0
    SWij_(2:nxg,nyg-1:nyg,2:nzg,2) = ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,2)**2d0 + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,6)**2d0 + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,7)**2d0
    ! dW/dx*dW/dx + dW/dy*dW/dy + dW/dz*dW/dz
    SWij_(2:nxg,2:3      ,2:nzg,3) = ten_buf_2(2:nxg,2:3      ,2:nzg,3)**2d0 + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,8)**2d0 + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,9)**2d0
    SWij_(2:nxg,nyg-1:nyg,2:nzg,3) = ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,3)**2d0 + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,8)**2d0 + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,9)**2d0

    ! dU/dx*dV/dx + dU/dy*dV/dy + dU/dz*dV/dz
    SWij_(2:nxg,2:3      ,2:nzg,4) = ten_buf_2(2:nxg,2:3      ,2:nzg,1) * ten_buf_2(2:nxg,2:3      ,2:nzg,7) + & 
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,4) * ten_buf_2(2:nxg,2:3      ,2:nzg,2) + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,5) * ten_buf_2(2:nxg,2:3      ,2:nzg,6) 
    SWij_(2:nxg,nyg-1:nyg,2:nzg,4) = ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,1) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,7) + & 
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,4) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,2) + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,5) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,6) 
    ! dU/dx*dW/dx + dU/dy*dW/dy + dU/dz*dW/dz
    SWij_(2:nxg,2:3      ,2:nzg,5) = ten_buf_2(2:nxg,2:3      ,2:nzg,1) * ten_buf_2(2:nxg,2:3      ,2:nzg,8) + & 
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,4) * ten_buf_2(2:nxg,2:3      ,2:nzg,9) + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,5) * ten_buf_2(2:nxg,2:3      ,2:nzg,3) 
    SWij_(2:nxg,nyg-1:nyg,2:nzg,5) = ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,1) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,8) + & 
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,4) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,9) + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,5) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,3) 
    ! dW/dx*dV/dx + dW/dy*dV/dy + dW/dz*dV/dz
    SWij_(2:nxg,2:3      ,2:nzg,6) = ten_buf_2(2:nxg,2:3      ,2:nzg,7) * ten_buf_2(2:nxg,2:3      ,2:nzg,8) + & 
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,2) * ten_buf_2(2:nxg,2:3      ,2:nzg,9) + &
                                    ten_buf_2(2:nxg,2:3      ,2:nzg,6) * ten_buf_2(2:nxg,2:3      ,2:nzg,3) 
    SWij_(2:nxg,nyg-1:nyg,2:nzg,6) = ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,7) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,8) + & 
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,2) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,9) + &
                                    ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,6) * ten_buf_2(2:nxg,nyg-1:nyg,2:nzg,3) 

  End Subroutine compute_Tij_estimate

  !----------------------------------------------------------------------!
  !                    Compute LWij and MWij                             !
  !                                                                      !
  !  MWij = dui/dn*duj/dn - c_test*fil_size_2^2*dhat(ui)/dn*dhat(uj)/dn  !
  !                                                                      !
  !  LWij = cu*ui*uj - c_test*hat(ui)*hat(uj)                            !
  !         + a1 * ( ui*uj                     + tauij          )        !  
  !         + a2 * ( hat(ui*uj)                + hat(tauij)     )        !
  !         + a3 * ( hat(hat(ui*uj))           + hat(hat(tauij) )        !
  !         + a4 * ( hat(ui)*hat(uj)           + Tij            )        !
  !         + a5 * ( hat(hat(ui)*hat(uj))      + hat(Tij)       )        !
  !         + a6 * ( hat(hat(ui))*hat(hat(uj)) + TTij           )        !
  !                                                                      !
  ! Input:  U_,V_,W_,nu_t_                                               !
  ! Output: LWij, MWij                                                   !
  !                                                                      !
  !----------------------------------------------------------------------!
  Subroutine compute_LWij_MWij(U_,V_,W_,nu_t_,cu,c_test,ai)

    Real   (Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real   (Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real   (Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real   (Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_
    Integer(Int32), Dimension(6),           Intent(In) :: ai
    Integer(Int32),                         Intent(In) :: cu,c_test

    ! local variables
    Integer(Int32) :: i, j, k, a1, a2, a3, a4, a5, a6
    Real   (Int64) :: alpha_x_1_local, alpha_x_1
    Real   (Int64) :: alpha_y_1_local, alpha_y_1

    a1 = ai(1)
    a2 = ai(2)
    a3 = ai(3)
    a4 = ai(4)
    a5 = ai(5)
    a6 = ai(6)


    !----------------------------------------------------------------------------!
    ! Compute LWij = (a1+cu)*ui*uj + a2*hat(ui*uj)      + a3*hat(hat(ui*uj))     ! 
    !                              + a4*hat(ui)*hat(uj) + a5*hat(hat(ui)*hat(uj))!
    !                              + a6*hat(hat(ui))*hat(hat(uj)                 !
    !----------------------------------------------------------------------------!    

    !----------------------------------------------------------------------------!
    ! Part 1: ui*uj, hat(ui*uj) and hat(hat(ui*uj))
    LWij = 0d0
    If ( (a1+cu)/=0 .or. a2/=0 .or. a3/=0 ) Then

      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(U_,     term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      
      Call interpolate_z(W_,     term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2) 

      ! fill in missing values (periodicity)
      Call apply_periodic_bc_x(term,  2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)

      ! SWij = ui*uj at v location (at wall)
      ! U at V location -> term  (2:nx-1,1:ny,:)
      ! W at V location -> term_2(:,1:ny,2:nzg-1)
      ! bottom wall 
      SWij(:,2:5,:,1)         = term  (2:nxg,1:4,2:nzg) * term  (2:nxg,1:4,2:nzg)  ! u^2
      SWij(:,2:5,:,2)         = V_    (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! v^2
      SWij(:,2:5,:,3)         = term_2(2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! w^2
      SWij(:,2:5,:,4)         = term  (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! uv
      SWij(:,2:5,:,5)         = term  (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! uw
      SWij(:,2:5,:,6)         = V_    (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! vw
      ! top wall
      SWij(:,nyg-4:nyg-1,:,1) = term  (2:nxg,ny-3:ny,2:nzg) * term  (2:nxg,ny-3:ny,2:nzg) ! u^2
      SWij(:,nyg-4:nyg-1,:,2) = V_    (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! v^2
      SWij(:,nyg-4:nyg-1,:,3) = term_2(2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! w^2
      SWij(:,nyg-4:nyg-1,:,4) = term  (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! uv
      SWij(:,nyg-4:nyg-1,:,5) = term  (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! uw
      SWij(:,nyg-4:nyg-1,:,6) = V_    (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! vw
      
      ! ui*uj
      LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a1+cu,8)*SWij(2:nxg-1,    2,2:nzg-1,:)
      LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a1+cu,8)*SWij(2:nxg-1,nyg-1,2:nzg-1,:)
      
      If ( a2/=0 .or. a3/=0 ) Then
         
         ten_buf(2:nxg,2:nyg,2:nzg,:) = SWij(:,:,:,:)
         Do i = 1,6
            Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
            Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
            Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
         End Do
         
         ! Compute SWij = hat(ui*uj) -> stored in SWij(2:nxg-1,2:nyg-1,2:nzg-1,:)
         Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),SWij(2:nxg-1,2:nyg-1,2:nzg-1,:)) 

         ! add it
         LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a2,8)*SWij(2:nxg-1,    2,2:nzg-1,:)
         LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a2,8)*SWij(2:nxg-1,nyg-1,2:nzg-1,:)
         
         If ( a3/= 0 ) Then
            
            ! SWij = hat(hat(ui*uj))
            ten_buf(2:nxg,2:nyg,2:nzg,:) = SWij(:,:,:,:)
            Do i = 1,6
               Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
               Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
               Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
            End Do
            Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),SWij(2:nxg-1,2:nyg-1,2:nzg-1,:))
            
            ! add it
            LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a3,8)*SWij(2:nxg-1,    2,2:nzg-1,:)
            LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a3,8)*SWij(2:nxg-1,nyg-1,2:nzg-1,:)

         End If

      End If
   End If

   !--------------------------------------------------------------------------!
   ! Part 2: hat(ui)*hat(uj), hat(hat(ui)*hat(uj)) and hathat(ui)*hathat(uj)   
   If ( (a4-c_test)/=0 .or. a5/= 0 .or. a6/=0 ) Then

      ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
      Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
      Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
      Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
         
      ! apply periodicity in x and z
      Call apply_periodic_bc_x(Uff,1)
      Call apply_periodic_bc_z(Uff,1)
      Call update_ghost_interior_planes(Uff,1)
      Call apply_periodic_bc_x(Vff,2)
      Call apply_periodic_bc_z(Vff,2)
      Call update_ghost_interior_planes(Vff,2)
      Call apply_periodic_bc_x(Wff,2)
      Call apply_periodic_bc_z(Wff,3)
      Call update_ghost_interior_planes(Wff,3)
      
      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(Uff,    term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      Call interpolate_z(Wff,    term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2)

      Call apply_periodic_bc_x(term  ,2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)
  
      ! hat(ui)*hat(uj)
      SWij(:,2:3        ,:,1) = term  (2:nxg,1:2   ,2:nzg) * term  (2:nxg,1:2     ,2:nzg) ! hat(u)^2
      SWij(:,2:3        ,:,2) = Vff   (2:nxg,1:2   ,2:nzg) * Vff   (2:nxg,1:2     ,2:nzg) ! hat(v)^2
      SWij(:,2:3        ,:,3) = term_2(2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(w)^2
      SWij(:,2:3        ,:,4) = term  (2:nxg,1:2   ,2:nzg) * Vff   (2:nxg,1:2     ,2:nzg) ! hat(u)hat(v)
      SWij(:,2:3        ,:,5) = term  (2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(u)hat(w)
      SWij(:,2:3        ,:,6) = Vff   (2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(v)hat(w)

      SWij(:,nyg-2:nyg-1,:,1) = term  (2:nxg,ny-1:ny,2:nzg) * term  (2:nxg,ny-1:ny,2:nzg) ! hat(u)^2
      SWij(:,nyg-2:nyg-1,:,2) = Vff   (2:nxg,ny-1:ny,2:nzg) * Vff   (2:nxg,ny-1:ny,2:nzg) ! hat(v)^2
      SWij(:,nyg-2:nyg-1,:,3) = term_2(2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(w)^2
      SWij(:,nyg-2:nyg-1,:,4) = term  (2:nxg,ny-1:ny,2:nzg) * Vff   (2:nxg,ny-1:ny,2:nzg) ! hat(u)hat(v)
      SWij(:,nyg-2:nyg-1,:,5) = term  (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(u)hat(w)
      SWij(:,nyg-2:nyg-1,:,6) = Vff   (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(v)hat(w)

      ! add it
      LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a4-c_test,8)*SWij(2:nxg-1,    2,2:nzg-1,:)
      LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a4-c_test,8)*SWij(2:nxg-1,nyg-1,2:nzg-1,:)

      If ( a5/=0 ) Then
         
         ! Compute SWij = hat(hat(ui)*hat(uj)) -> stored in SWij(3:nxg-1,2:nyg-1,3:nzg-1,:)         
         ten_buf(2:nxg,2:nyg,2:nzg,:) = SWij(:,:,:,:)
         Do i = 1,6
            Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
            Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
            Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
         End Do
         Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),SWij(2:nxg-1,2:nyg-1,2:nzg-1,:)) 
         
         ! add it
         LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a5,8)*SWij(2:nxg-1,    2,2:nzg-1,:)
         LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a5,8)*SWij(2:nxg-1,nyg-1,2:nzg-1,:)

      End If

      If ( a6/=0 ) Then
         
         ! hat(hat(ui))
         Call filter_xzy( Uff(1:nx, 1:nyg,1:nzg), term(2:nx-1 ,1:nyg,2:nzg-1) )

         Call apply_periodic_bc_x(term(1:nx ,1:nyg,1:nzg),1)
         Call apply_periodic_bc_z(term(1:nx ,1:nyg,1:nzg),1)
         Call update_ghost_interior_planes(term(1:nx ,1:nyg,1:nzg),1)
         
         Call interpolate_x(term(1:nx ,1:nyg,1:nzg), term_1(2:nx,:,:),1) 
         Call interpolate_y(term_1                 , term  (:,1:ny,:),2)
         
         Call filter_xzy( Wff(1:nxg,1:nyg,1:nz ), term_2(2:nxg-1,1:nyg,2:nz-1 ) )
         Call apply_periodic_bc_x(term_2(1:nxg,1:nyg,1:nz ),2)
         Call apply_periodic_bc_z(term_2(1:nxg,1:nyg,1:nz ),3)
         Call update_ghost_interior_planes(term_2(1:nxg,1:nyg,1:nz ),3)
         
         Call interpolate_z(term_2(1:nxg,1:nyg,1:nz ), term_1(:,:,2:nz),1) 
         Call interpolate_y(term_1                   , term_2(:,1:ny,:),2)

         Call filter_xzy( Vff(1:nxg,1:ny, 1:nzg), term_1(2:nxg-1,1:ny ,2:nzg-1) )
         
         Call apply_periodic_bc_x(term_1(1:nxg,1:ny, 1:nzg),2)
         Call apply_periodic_bc_z(term_1(1:nxg,1:ny, 1:nzg),2)
         Call update_ghost_interior_planes(term_1(1:nxg,1:ny, 1:nzg),2)
         
         Call apply_periodic_bc_x(term  ,2)
         Call apply_periodic_bc_z(term_2,4)
         Call update_ghost_interior_planes(term_2,4)
         
         ! hat(hat(ui))*hat(hat(uj))
         SWij(:,2    ,:,1) = term  (2:nxg,1 ,2:nzg) * term  (2:nxg,1 ,2:nzg) ! hathat(u)^2
         SWij(:,2    ,:,2) = term_1(2:nxg,1 ,2:nzg) * term_1(2:nxg,1 ,2:nzg) ! hathat(v)^2
         SWij(:,2    ,:,3) = term_2(2:nxg,1 ,2:nzg) * term_2(2:nxg,1 ,2:nzg) ! hathat(w)^2
         SWij(:,2    ,:,4) = term  (2:nxg,1 ,2:nzg) * term_1(2:nxg,1 ,2:nzg) ! hathat(u)hathat(v)
         SWij(:,2    ,:,5) = term  (2:nxg,1 ,2:nzg) * term_2(2:nxg,1 ,2:nzg) ! hathat(u)hathat(w)
         SWij(:,2    ,:,6) = term_1(2:nxg,1 ,2:nzg) * term_2(2:nxg,1 ,2:nzg) ! hathat(v)hathat(w)

         SWij(:,nyg-1,:,1) = term  (2:nxg,ny,2:nzg) * term  (2:nxg,ny,2:nzg) ! hathat(u)^2
         SWij(:,nyg-1,:,2) = term_1(2:nxg,ny,2:nzg) * term_1(2:nxg,ny,2:nzg) ! hathat(v)^2
         SWij(:,nyg-1,:,3) = term_2(2:nxg,ny,2:nzg) * term_2(2:nxg,ny,2:nzg) ! hathat(w)^2
         SWij(:,nyg-1,:,4) = term  (2:nxg,ny,2:nzg) * term_1(2:nxg,ny,2:nzg) ! hathat(u)hathat(v)
         SWij(:,nyg-1,:,5) = term  (2:nxg,ny,2:nzg) * term_2(2:nxg,ny,2:nzg) ! hathat(u)hathat(w)
         SWij(:,nyg-1,:,6) = term_1(2:nxg,ny,2:nzg) * term_2(2:nxg,ny,2:nzg) ! hathat(v)hathat(w)
         
         ! add it
         LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a6,8)*SWij(2:nxg-1,    2,2:nzg-1,:)
         LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a6,8)*SWij(2:nxg-1,nyg-1,2:nzg-1,:)
         
      End If
   End If
   
   !-----------------------------------------------------------------------------!
   ! Part 3: Compute SGS tensor contributions to LWij
   If ( Dirichlet_nu_t==0 ) Then
      
      If (est_Tij == 0) then

        If ( a1/=0 .or. a2/=0 .or. a3/=0 ) Then
           
           !----------compute tau_ij, hat(tau_ij) and hat(hat(tau_ij))   
           ! compute SWij at V locations in the first two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
           Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
           
           ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
           Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)
            
           ! tau_ij = -2*nu_t*SWij -> SWij
           Do i = 2, nxg
              Do k = 2, nzg
                 Do j = 1, 6
                    SWij(i,2:nyg,k,j) = -2d0*term_2(i,1:ny,k)*SWij(i,2:nyg,k,j)
                 End Do
              End Do
           End Do
           
           ! add it
           LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a1,8) * SWij(2:nxg-1,  2,2:nzg-1,:)
           LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a1,8) * SWij(2:nxg-1,nyg,2:nzg-1,:)
           
           If ( a2/=0 .or. a3/=0 ) Then
            
              ! hat(tau_ij) -> filter tau_ij 
              ten_buf(2:nxg,2:nyg,2:nzg,1:6) = SWij;
              Do i = 1,6
                 Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
                 Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
                 Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
              End Do
              Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg,1:nzg,:),SWij(2:nxg-1,2:nyg,2:nzg-1,:)) 
              
              ! add it
              LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a2,8) * SWij(2:nxg-1,  2,2:nzg-1,:)
              LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a2,8) * SWij(2:nxg-1,nyg,2:nzg-1,:)
              
              If ( a3/=0 )  Then

                 ! hathat(tau_ij) -> filter hat(tau_ij)
                 ten_buf(2:nxg,2:nyg,2:nzg,1:6) = SWij;
                 Do i = 1,6
                    Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
                    Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
                    Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
                 End Do
                 Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg,1:nzg,:),SWij(2:nxg-1,2:nyg,2:nzg-1,:)) 
                 
                 ! add it
                 LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a3,8) * SWij(2:nxg-1,  2,2:nzg-1,:)
                 LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a3,8) * SWij(2:nxg-1,nyg,2:nzg-1,:)
              End If

           End If
        End If
        
        If ( a4/=0 .or. a5/=0 .or. a6/=0 ) Then

           !----------compute T_ij
           ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
           Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
           Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
           Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
           
           ! apply periodicity in x and z
           Call apply_periodic_bc_x(Uff,1)
           Call apply_periodic_bc_z(Uff,1)
           Call update_ghost_interior_planes(Uff,1)
           Call apply_periodic_bc_x(Vff,2)
           Call apply_periodic_bc_z(Vff,2)
           Call update_ghost_interior_planes(Vff,2)
           Call apply_periodic_bc_x(Wff,2)
           Call apply_periodic_bc_z(Wff,3)
           Call update_ghost_interior_planes(Wff,3)
           
           ! compute eddy viscosity for filtered velocities
           Call compute_eddy_viscosity(Uff,Vff,Wff,avg_nu_t_hat,nu_t_hat)
           
           ! interpolate nu_t_hat from cell centers to cell faces (== cell edges because averaged in xz)
           Call interpolate_y(nu_t_hat,term_2(:,1:ny,:),2)    
           
           ! compute hat(SWij) at V locations in two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
           Call compute_Sij_at_V_location_wall(Uff,Vff,Wff,SWij)
           
           ! Tij = -2*nu_t_hat*hat(SWij) -> SWij
           Do i = 2, nxg
              Do k = 2, nzg
                 Do j = 1, 6
                    SWij(i,2:nyg,k,j) = -2d0*term_2(i,1:ny,k)*SWij(i,2:nyg,k,j) 
                 End Do
              End Do
           End Do
           
           ! add it
           LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a4,8) * SWij(2:nxg-1,  2,2:nzg-1,:)
           LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a4,8) * SWij(2:nxg-1,nyg,2:nzg-1,:)

           If ( a5/=0 ) Then
            
              ! hat(T_ij) 
              ten_buf(2:nxg,2:nyg,2:nzg,1:6) = SWij;
              Do i = 1,6
                 Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
                 Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
                 Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
              End Do
              Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg,1:nzg,:),SWij(2:nxg-1,2:nyg,2:nzg-1,:)) 
              
              ! add it
              LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a5,8) * SWij(2:nxg-1,  2,2:nzg-1,:)
              LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a5,8) * SWij(2:nxg-1,nyg,2:nzg-1,:)

           End If

           If ( a6/=0 )  Then

              ! hat(hat(ui))
              Call filter_xzy( Uff(1:nx, 1:nyg,1:nzg), Ufff(2:nx-1 ,1:nyg,2:nzg-1) )
              Call filter_xzy( Vff(1:nxg,1:ny, 1:nzg), Vfff(2:nxg-1,1:ny ,2:nzg-1) )
              Call filter_xzy( Wff(1:nxg,1:nyg,1:nz ), Wfff(2:nxg-1,1:nyg,2:nz-1 ) )

              Call apply_periodic_bc_x(Ufff,1)
              Call apply_periodic_bc_z(Ufff,1)
              Call update_ghost_interior_planes(Ufff,1)
              
              Call apply_periodic_bc_x(Vfff,2)
              Call apply_periodic_bc_z(Vfff,2)
              Call update_ghost_interior_planes(Vfff,2)
              
              Call apply_periodic_bc_x(Wfff,2)
              Call apply_periodic_bc_z(Wfff,3)
              Call update_ghost_interior_planes(Wfff,3)

              
              ! compute hat(SWij) at V locations in two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
              Call compute_Sij_at_V_location_wall(Ufff,Vfff,Wfff,SWij)

              ! compute eddy viscosity for filtered velocities
              Call compute_eddy_viscosity(Ufff,Vfff,Wfff,avg_nu_t_hat,nu_t_hat)
              
              ! interpolate nu_t_hat from cell centers to cell faces (== cell edges because averaged in xz)
              Call interpolate_y(nu_t_hat,term_2(:,1:ny,:),2)    
              
              ! TTij = -2*nu_t_hat*hat(SWij) -> SWij
              Do i = 2, nxg
                 Do k = 2, nzg
                    Do j = 1, 6
                       SWij(i,2:nyg,k,j) = -2d0*term_2(i,1:ny,k)*SWij(i,2:nyg,k,j) 
                    End Do
                 End Do
              End Do
              
              ! add it
              LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + Real(a6,8) * SWij(2:nxg-1,  2,2:nzg-1,:)
              LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + Real(a6,8) * SWij(2:nxg-1,nyg,2:nzg-1,:)
            
           End If

        End If
        
     Else

        Stop "Cannot use this model for Tij here!"

     End If

  End If ! Dirichlet_nu_t==0
  
  !-----------------------------------------------------------------------------!
  ! Part 3: Compute MWij = dui/dy*duj/dy - c_test*fil_size_2^2*dhat(ui)/dy*dhat(uj)/dy

  !-----------------------------------------------------------------------!
  ! compute dui/dy*duj/dy

  ! interpolate velocity to cell centers
  Call interpolate_x(U_,term  (2:nx,:,:),1) 
  Call interpolate_z(W_,term_2(:,:,2:nz),1) 
  
  Call apply_periodic_bc_x(term  ,2)
  Call apply_periodic_bc_z(term_2,4)
  Call update_ghost_interior_planes(term_2,4)

  ! Compute dui/dn at V location (at wall)
  ! SWij(:,:,:,1) -> dU/dn
  ! SWij(:,:,:,2) -> dV/dn
  ! SWij(:,:,:,3) -> dW/dn 
  Do i = 2, nxg
     Do k = 2, nzg
        SWij(i,2,k,1)     = (term  (i,2,k)     - term  (i,1,k)  ) / (yg(2)-yg(1))            ! dU/dy at lower wall
        SWij(i,2,k,2)     = (V_    (i,2,k)     - V_    (i,1,k)  ) / (y (2)-y (1))            ! dV/dy at lower wall (1st order)
        SWij(i,2,k,3)     = (term_2(i,2,k)     - term_2(i,1,k)  ) / (yg(2)-yg(1))            ! dW/dy at lower wall
        
        SWij(i,nyg-1,k,1) = (term  (i,nyg,k) - term  (i,nyg-1,k)) / (yg(nyg)-yg(nyg-1))      ! dU/dy at upper wall
        SWij(i,nyg-1,k,2) = (V_    (i,ny ,k) - V_    (i,ny -1,k)) / (y ( ny)-y ( ny-1))      ! dV/dy at upper wall (1st order)
        SWij(i,nyg-1,k,3) = (term_2(i,nyg,k) - term_2(i,nyg-1,k)) / (yg(nyg)-yg(nyg-1))      ! dU/dy at upper wall
     End Do
  End Do
  
  Do i = 2, nxg-1
     Do k = 2, nzg-1
        MWij(i,    2,k,1) = SWij(i,    2,k,1)*SWij(i,    2,k,1)
        MWij(i,    2,k,2) = SWij(i,    2,k,2)*SWij(i,    2,k,2)
        MWij(i,    2,k,3) = SWij(i,    2,k,3)*SWij(i,    2,k,3)
        MWij(i,    2,k,4) = SWij(i,    2,k,1)*SWij(i,    2,k,2)
        MWij(i,    2,k,5) = SWij(i,    2,k,1)*SWij(i,    2,k,3)
        MWij(i,    2,k,6) = SWij(i,    2,k,2)*SWij(i,    2,k,3)
        
        MWij(i,nyg-1,k,1) = SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,1)
        MWij(i,nyg-1,k,2) = SWij(i,nyg-1,k,2)*SWij(i,nyg-1,k,2)
        MWij(i,nyg-1,k,3) = SWij(i,nyg-1,k,3)*SWij(i,nyg-1,k,3)
        MWij(i,nyg-1,k,4) = SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,2)
        MWij(i,nyg-1,k,5) = SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,3)
        MWij(i,nyg-1,k,6) = SWij(i,nyg-1,k,2)*SWij(i,nyg-1,k,3)
     End Do
  End Do

  !-----------------------------------------------------------------------!
  ! compute dhat(ui)/dy*dhat(uj)/dy
  If ( Abs(c_test)>0 ) Then

     ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
     Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
     Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
     Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
     
     ! apply periodicity in x and z
     Call apply_periodic_bc_x(Uff,1)
     Call apply_periodic_bc_z(Uff,1)
     Call update_ghost_interior_planes(Uff,1)
     Call apply_periodic_bc_x(Vff,2)
     Call apply_periodic_bc_z(Vff,2)
     Call update_ghost_interior_planes(Vff,2)
     Call apply_periodic_bc_x(Wff,2)
     Call apply_periodic_bc_z(Wff,3)
     Call update_ghost_interior_planes(Wff,3)
     
     ! interpolate velocity to cell centers
     Call interpolate_x(Uff,term  (2:nx,:,:),1) 
     Call interpolate_z(Wff,term_2(:,:,2:nz),1) 
     
     Call apply_periodic_bc_x(term  ,2)
     Call apply_periodic_bc_z(term_2,4)
     Call update_ghost_interior_planes(term_2,4)
     
     ! Compute dhat(ui)/dn at V location (at wall)
     ! SWij(:,:,:,1) -> dhat(U)/dn
     ! SWij(:,:,:,2) -> dhat(V)/dn
     ! SWij(:,:,:,3) -> dhat(W)/dn 
     Do i = 2, nxg
        Do k = 2, nzg
           SWij(i,2,k,1)     = (term  (i,2,k)     - term  (i,1,k)  ) / (yg(2)-yg(1))            ! dhat(U)/dy at lower wall
           SWij(i,2,k,2)     = (Vff   (i,2,k)     - Vff   (i,1,k)  ) / (y (2)-y (1))            ! dhat(V)/dy at lower wall (1st order)
           SWij(i,2,k,3)     = (term_2(i,2,k)     - term_2(i,1,k)  ) / (yg(2)-yg(1))            ! dhat(W)/dy at lower wall
           
           SWij(i,nyg-1,k,1) = (term  (i,nyg,k) - term  (i,nyg-1,k)) / (yg(nyg)-yg(nyg-1))      ! dhat(U)/dy at upper wall
           SWij(i,nyg-1,k,2) = (Vff   (i,ny ,k) - Vff   (i,ny -1,k)) / (y ( ny)-y ( ny-1))      ! dhat(V)/dy at upper wall (1st order)
           SWij(i,nyg-1,k,3) = (term_2(i,nyg,k) - term_2(i,nyg-1,k)) / (yg(nyg)-yg(nyg-1))      ! dhat(W)/dy at upper wall
        End Do
     End Do
     
     Do i = 2, nxg-1
        Do k = 2, nzg-1
           MWij(i,    2,k,1) = MWij(i,    2,k,1) - c_test*fil_size_2**2d0*SWij(i,    2,k,1)*SWij(i,    2,k,1)
           MWij(i,    2,k,2) = MWij(i,    2,k,2) - c_test*fil_size_2**2d0*SWij(i,    2,k,2)*SWij(i,    2,k,2)
           MWij(i,    2,k,3) = MWij(i,    2,k,3) - c_test*fil_size_2**2d0*SWij(i,    2,k,3)*SWij(i,    2,k,3)
           MWij(i,    2,k,4) = MWij(i,    2,k,4) - c_test*fil_size_2**2d0*SWij(i,    2,k,1)*SWij(i,    2,k,2)
           MWij(i,    2,k,5) = MWij(i,    2,k,5) - c_test*fil_size_2**2d0*SWij(i,    2,k,1)*SWij(i,    2,k,3)
           MWij(i,    2,k,6) = MWij(i,    2,k,6) - c_test*fil_size_2**2d0*SWij(i,    2,k,2)*SWij(i,    2,k,3)           
           MWij(i,nyg-1,k,1) = MWij(i,nyg-1,k,1) - c_test*fil_size_2**2d0*SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,1)
           MWij(i,nyg-1,k,2) = MWij(i,nyg-1,k,2) - c_test*fil_size_2**2d0*SWij(i,nyg-1,k,2)*SWij(i,nyg-1,k,2)
           MWij(i,nyg-1,k,3) = MWij(i,nyg-1,k,3) - c_test*fil_size_2**2d0*SWij(i,nyg-1,k,3)*SWij(i,nyg-1,k,3)
           MWij(i,nyg-1,k,4) = MWij(i,nyg-1,k,4) - c_test*fil_size_2**2d0*SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,2)
           MWij(i,nyg-1,k,5) = MWij(i,nyg-1,k,5) - c_test*fil_size_2**2d0*SWij(i,nyg-1,k,1)*SWij(i,nyg-1,k,3)
           MWij(i,nyg-1,k,6) = MWij(i,nyg-1,k,6) - c_test*fil_size_2**2d0*SWij(i,nyg-1,k,2)*SWij(i,nyg-1,k,3)
        End Do
     End Do

  End If

End Subroutine compute_LWij_MWij


  !----------------------------------------------------------------------!
  !                    Compute LWij and MWij                             !
  !                                                                      !
  !                                                                      !
  !  icase = 1-->tau_ij      (lv1), hat(tau_ij) (lv2)                    !
  !          2-->tau_ij      (lv1), T_ij        (lv2)                    !
  !          3-->hat(tau_ij) (lv1), T_ij        (lv2)                    !
  !                                                                      !
  !                                                                      !
  !  !!!!!!!!!! only a1, a2, a3 implemented at this time !!!!!!!!!!!     !
  !                                                                      !
  !                                                                      !
  !  MWij = -Delta^2 * (a1*|S|*Sij + a2*Sik*Skj + a3*Rik*Rkj +           ! 
  !                     a4*(Sik*Rkj - Rik*Skj) + a5*1/|S|*(SSR-RSS))     !
  !                                                            at lv1    !
  !         +Delta^2 * (a1*|S|*Sij + a2*Sik*Skj + a3*Rik*Rkj +           ! 
  !                     a4*(Sik*Rkj - Rik*Skj) + a5*1/|S|*(SSR-RSS))     !
  !                                                            at lv2    !
  !                                                                      !
  !  LWij =  (-ui*uj - tauSGS_ij + 2*nu*Sij)                             ! 
  !                                                            at lv1    !
  !         -(-ui*uj - tauSGS_ij + 2*nu*Sij)                             ! 
  !                                                            at lv2    !
  !                                                                      !
  !                                                                      !
  ! Input:  U_,V_,W_,nu_t_                                               !
  ! Output: LWij, MWij                                                   !
  !                                                                      !
  !----------------------------------------------------------------------!
  Subroutine compute_LWij_MWij_CW(U_,V_,W_,nu_t_,icase,ai)

    Real   (Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real   (Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real   (Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real   (Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_
    Integer(Int32), Dimension(6),           Intent(In) :: ai
    Integer(Int32),                         Intent(In) :: icase 

    ! local variables
    Integer(Int32) :: i, j, k, l, a1, a2, a3, a4, a5, ind(3,3) 
    Real   (Int64) :: r_sign(3,3), Delta

    Delta = ((y(2)-y(1)) * (x(2)-x(1)) * (z(2)-z(1)))**(1d0/3d0)

    a1 = ai(1)
    a2 = ai(2)
    a3 = ai(3)
    a4 = ai(4)
    a5 = ai(5)

    ! index for tensor 
    ind(1,1) = 1;
    ind(2,2) = 2;
    ind(3,3) = 3;
    ind(1,2) = 4;
    ind(2,1) = 4;
    ind(1,3) = 5;
    ind(3,1) = 5;
    ind(2,3) = 6;
    ind(3,2) = 6;

    ! sign for antisymmetric tensor Rij
    r_sign(1,1) =  1d0;
    r_sign(2,2) =  1d0;
    r_sign(3,3) =  1d0;
    r_sign(1,2) =  1d0;
    r_sign(2,1) = -1d0;
    r_sign(1,3) =  1d0;
    r_sign(3,1) = -1d0;
    r_sign(2,3) =  1d0;
    r_sign(3,2) = -1d0;


    ! Compute LWij 
    LWij = 0d0
    If     (icase == 1) Then
      !--------------------Part 1------------------------------!
      ! Compute the stress for hat(tau_w)

      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(U_,     term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      
      Call interpolate_z(W_,     term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2) 

      ! fill in missing values (periodicity)
      Call apply_periodic_bc_x(term,  2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)

      ! SWij = ui*uj at v location (at wall)
      ! U at V location -> term  (2:nx-1,1:ny,:)
      ! W at V location -> term_2(:,1:ny,2:nzg-1)
      ! bottom wall 
      SWij(:,2:5,:,1)         = term  (2:nxg,1:4,2:nzg) * term  (2:nxg,1:4,2:nzg)  ! u^2
      SWij(:,2:5,:,2)         = V_    (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! v^2
      SWij(:,2:5,:,3)         = term_2(2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! w^2
      SWij(:,2:5,:,4)         = term  (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! uv
      SWij(:,2:5,:,5)         = term  (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! uw
      SWij(:,2:5,:,6)         = V_    (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! vw
      ! top wall
      SWij(:,nyg-4:nyg-1,:,1) = term  (2:nxg,ny-3:ny,2:nzg) * term  (2:nxg,ny-3:ny,2:nzg) ! u^2
      SWij(:,nyg-4:nyg-1,:,2) = V_    (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! v^2
      SWij(:,nyg-4:nyg-1,:,3) = term_2(2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! w^2
      SWij(:,nyg-4:nyg-1,:,4) = term  (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! uv
      SWij(:,nyg-4:nyg-1,:,5) = term  (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! uw
      SWij(:,nyg-4:nyg-1,:,6) = V_    (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! vw
      
      ! ui*uj
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) + SWij(2:nxg-1,        2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) + SWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:)

      ! compute tau_ij, hat(tau_ij) and hat(hat(tau_ij))   
      ! compute SWij at V locations in the first two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
           
      ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
      Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)
            
      ! tau_ij = -2*nu_t*SWij -> SWij
      Do i = 2, nxg
        Do k = 2, nzg
           Do j = 1, 6
              SWij(i,2:nyg,k,j) = -2d0*(nu+term_2(i,1:ny,k))*SWij(i,2:nyg,k,j)
           End Do
        End Do
      End Do
           
      ! add it
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) + SWij(2:nxg-1,      2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,:)

      ! filter the -LWij
      ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,1:6) = LWij;
      Do i = 1,6
        Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
        Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
        Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
      End Do
      Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),LWij(2:nxg-1,2:nyg-1,2:nzg-1,:)) 


      !--------------------Part 2------------------------------!
      ! Compute the stress for tau_w

      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(U_,     term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      
      Call interpolate_z(W_,     term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2) 

      ! fill in missing values (periodicity)
      Call apply_periodic_bc_x(term,  2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)

      ! SWij = ui*uj at v location (at wall)
      ! U at V location -> term  (2:nx-1,1:ny,:)
      ! W at V location -> term_2(:,1:ny,2:nzg-1)
      ! bottom wall 
      SWij(:,2:5,:,1)         = term  (2:nxg,1:4,2:nzg) * term  (2:nxg,1:4,2:nzg)  ! u^2
      SWij(:,2:5,:,2)         = V_    (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! v^2
      SWij(:,2:5,:,3)         = term_2(2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! w^2
      SWij(:,2:5,:,4)         = term  (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! uv
      SWij(:,2:5,:,5)         = term  (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! uw
      SWij(:,2:5,:,6)         = V_    (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! vw
      ! top wall
      SWij(:,nyg-4:nyg-1,:,1) = term  (2:nxg,ny-3:ny,2:nzg) * term  (2:nxg,ny-3:ny,2:nzg) ! u^2
      SWij(:,nyg-4:nyg-1,:,2) = V_    (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! v^2
      SWij(:,nyg-4:nyg-1,:,3) = term_2(2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! w^2
      SWij(:,nyg-4:nyg-1,:,4) = term  (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! uv
      SWij(:,nyg-4:nyg-1,:,5) = term  (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! uw
      SWij(:,nyg-4:nyg-1,:,6) = V_    (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! vw
      
      ! ui*uj
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,        2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:)

      ! compute tau_ij, hat(tau_ij) and hat(hat(tau_ij))   
      ! compute SWij at V locations in the first two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
           
      ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
      Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)
            
      ! tau_ij = -2*nu_t*SWij -> SWij
      Do i = 2, nxg
        Do k = 2, nzg
           Do j = 1, 6
              SWij(i,2:nyg,k,j) = -2d0*(nu+term_2(i,1:ny,k))*SWij(i,2:nyg,k,j)
           End Do
        End Do
      End Do
           
      ! add it
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,      2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,:)
    Elseif (icase == 2) Then
      !--------------------Part 1------------------------------!
      ! Compute the stress for tau_w

      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(U_,     term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      
      Call interpolate_z(W_,     term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2) 

      ! fill in missing values (periodicity)
      Call apply_periodic_bc_x(term,  2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)

      ! SWij = ui*uj at v location (at wall)
      ! U at V location -> term  (2:nx-1,1:ny,:)
      ! W at V location -> term_2(:,1:ny,2:nzg-1)
      ! bottom wall 
      SWij(:,2:5,:,1)         = term  (2:nxg,1:4,2:nzg) * term  (2:nxg,1:4,2:nzg)  ! u^2
      SWij(:,2:5,:,2)         = V_    (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! v^2
      SWij(:,2:5,:,3)         = term_2(2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! w^2
      SWij(:,2:5,:,4)         = term  (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! uv
      SWij(:,2:5,:,5)         = term  (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! uw
      SWij(:,2:5,:,6)         = V_    (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! vw
      ! top wall
      SWij(:,nyg-4:nyg-1,:,1) = term  (2:nxg,ny-3:ny,2:nzg) * term  (2:nxg,ny-3:ny,2:nzg) ! u^2
      SWij(:,nyg-4:nyg-1,:,2) = V_    (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! v^2
      SWij(:,nyg-4:nyg-1,:,3) = term_2(2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! w^2
      SWij(:,nyg-4:nyg-1,:,4) = term  (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! uv
      SWij(:,nyg-4:nyg-1,:,5) = term  (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! uw
      SWij(:,nyg-4:nyg-1,:,6) = V_    (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! vw
      
      ! ui*uj
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,        2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:)

      ! compute tau_ij, hat(tau_ij) and hat(hat(tau_ij))   
      ! compute SWij at V locations in the first two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
           
      ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
      Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)
            
      ! tau_ij = -2*nu_t*SWij -> SWij
      Do i = 2, nxg
        Do k = 2, nzg
           Do j = 1, 6
              SWij(i,2:nyg,k,j) = -2d0*(nu+term_2(i,1:ny,k))*SWij(i,2:nyg,k,j)
           End Do
        End Do
      End Do
           
      ! add it
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,      2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,:)

      !--------------------Part 2------------------------------!
      ! Compute the stress for T_w

      ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
      Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
      Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
      Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
         
      ! apply periodicity in x and z
      Call apply_periodic_bc_x(Uff,1)
      Call apply_periodic_bc_z(Uff,1)
      Call update_ghost_interior_planes(Uff,1)
      Call apply_periodic_bc_x(Vff,2)
      Call apply_periodic_bc_z(Vff,2)
      Call update_ghost_interior_planes(Vff,2)
      Call apply_periodic_bc_x(Wff,2)
      Call apply_periodic_bc_z(Wff,3)
      Call update_ghost_interior_planes(Wff,3)
      
      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(Uff,    term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      Call interpolate_z(Wff,    term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2)

      Call apply_periodic_bc_x(term  ,2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)
  
      ! hat(ui)*hat(uj)
      SWij(:,2:3        ,:,1) = term  (2:nxg,1:2   ,2:nzg) * term  (2:nxg,1:2     ,2:nzg) ! hat(u)^2
      SWij(:,2:3        ,:,2) = Vff   (2:nxg,1:2   ,2:nzg) * Vff   (2:nxg,1:2     ,2:nzg) ! hat(v)^2
      SWij(:,2:3        ,:,3) = term_2(2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(w)^2
      SWij(:,2:3        ,:,4) = term  (2:nxg,1:2   ,2:nzg) * Vff   (2:nxg,1:2     ,2:nzg) ! hat(u)hat(v)
      SWij(:,2:3        ,:,5) = term  (2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(u)hat(w)
      SWij(:,2:3        ,:,6) = Vff   (2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(v)hat(w)

      SWij(:,nyg-2:nyg-1,:,1) = term  (2:nxg,ny-1:ny,2:nzg) * term  (2:nxg,ny-1:ny,2:nzg) ! hat(u)^2
      SWij(:,nyg-2:nyg-1,:,2) = Vff   (2:nxg,ny-1:ny,2:nzg) * Vff   (2:nxg,ny-1:ny,2:nzg) ! hat(v)^2
      SWij(:,nyg-2:nyg-1,:,3) = term_2(2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(w)^2
      SWij(:,nyg-2:nyg-1,:,4) = term  (2:nxg,ny-1:ny,2:nzg) * Vff   (2:nxg,ny-1:ny,2:nzg) ! hat(u)hat(v)
      SWij(:,nyg-2:nyg-1,:,5) = term  (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(u)hat(w)
      SWij(:,nyg-2:nyg-1,:,6) = Vff   (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(v)hat(w)

      ! add it
      LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + SWij(2:nxg-1,    2,2:nzg-1,:)
      LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + SWij(2:nxg-1,nyg-1,2:nzg-1,:)


      ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
      Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
      Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
      Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
           
      ! apply periodicity in x and z
      Call apply_periodic_bc_x(Uff,1)
      Call apply_periodic_bc_z(Uff,1)
      Call update_ghost_interior_planes(Uff,1)
      Call apply_periodic_bc_x(Vff,2)
      Call apply_periodic_bc_z(Vff,2)
      Call update_ghost_interior_planes(Vff,2)
      Call apply_periodic_bc_x(Wff,2)
      Call apply_periodic_bc_z(Wff,3)
      Call update_ghost_interior_planes(Wff,3)
           
      ! compute eddy viscosity for filtered velocities
      Call compute_eddy_viscosity(Uff,Vff,Wff,avg_nu_t_hat,nu_t_hat)
           
      ! interpolate nu_t_hat from cell centers to cell faces (== cell edges because averaged in xz)
      Call interpolate_y(nu_t_hat,term_2(:,1:ny,:),2)    
           
      ! compute hat(SWij) at V locations in two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
      Call compute_Sij_at_V_location_wall(Uff,Vff,Wff,SWij)
           
      ! Tij = -2*nu_t_hat*hat(SWij) -> SWij
      Do i = 2, nxg
        Do k = 2, nzg
          Do j = 1, 6
            SWij(i,2:nyg,k,j) = -2d0*(nu+term_2(i,1:ny,k))*SWij(i,2:nyg,k,j) 
          End Do
        End Do
      End Do
           
      ! add it
      LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + SWij(2:nxg-1,  2,2:nzg-1,:)
      LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + SWij(2:nxg-1,nyg,2:nzg-1,:)
    Elseif (icase == 3) Then
      !--------------------Part 1------------------------------!
      ! Compute the stress for hat(tau_w)

      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(U_,     term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      
      Call interpolate_z(W_,     term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2) 

      ! fill in missing values (periodicity)
      Call apply_periodic_bc_x(term,  2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)

      ! SWij = ui*uj at v location (at wall)
      ! U at V location -> term  (2:nx-1,1:ny,:)
      ! W at V location -> term_2(:,1:ny,2:nzg-1)
      ! bottom wall 
      SWij(:,2:5,:,1)         = term  (2:nxg,1:4,2:nzg) * term  (2:nxg,1:4,2:nzg)  ! u^2
      SWij(:,2:5,:,2)         = V_    (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! v^2
      SWij(:,2:5,:,3)         = term_2(2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! w^2
      SWij(:,2:5,:,4)         = term  (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! uv
      SWij(:,2:5,:,5)         = term  (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! uw
      SWij(:,2:5,:,6)         = V_    (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! vw
      ! top wall
      SWij(:,nyg-4:nyg-1,:,1) = term  (2:nxg,ny-3:ny,2:nzg) * term  (2:nxg,ny-3:ny,2:nzg) ! u^2
      SWij(:,nyg-4:nyg-1,:,2) = V_    (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! v^2
      SWij(:,nyg-4:nyg-1,:,3) = term_2(2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! w^2
      SWij(:,nyg-4:nyg-1,:,4) = term  (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! uv
      SWij(:,nyg-4:nyg-1,:,5) = term  (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! uw
      SWij(:,nyg-4:nyg-1,:,6) = V_    (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! vw
      
      ! ui*uj
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,        2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:)

      ! compute tau_ij, hat(tau_ij) and hat(hat(tau_ij))   
      ! compute SWij at V locations in the first two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
           
      ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
      Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)
            
      ! tau_ij = -2*nu_t*SWij -> SWij
      Do i = 2, nxg
        Do k = 2, nzg
           Do j = 1, 6
              SWij(i,2:nyg,k,j) = -2d0*(nu+term_2(i,1:ny,k))*SWij(i,2:nyg,k,j)
           End Do
        End Do
      End Do
           
      ! add it
      LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,      2:5,2:nzg-1,:)
      LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,:)

      ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,1:6) = LWij;
      Do i = 1,6
        Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
        Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
        Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
      End Do
      Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),LWij(2:nxg-1,2:nyg-1,2:nzg-1,:)) 

      !--------------------Part 2------------------------------!
      ! Compute the stress for T_w 

      ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
      Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
      Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
      Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
         
      ! apply periodicity in x and z
      Call apply_periodic_bc_x(Uff,1)
      Call apply_periodic_bc_z(Uff,1)
      Call update_ghost_interior_planes(Uff,1)
      Call apply_periodic_bc_x(Vff,2)
      Call apply_periodic_bc_z(Vff,2)
      Call update_ghost_interior_planes(Vff,2)
      Call apply_periodic_bc_x(Wff,2)
      Call apply_periodic_bc_z(Wff,3)
      Call update_ghost_interior_planes(Wff,3)
      
      ! interpolate velocity to U and W to V location (at wall)
      Call interpolate_x(Uff,    term_1(2:nx,:,:),1) 
      Call interpolate_y(term_1, term  (:,1:ny,:),2)
      Call interpolate_z(Wff,    term_1(:,:,2:nz),1) 
      Call interpolate_y(term_1, term_2(:,1:ny,:),2)

      Call apply_periodic_bc_x(term  ,2)
      Call apply_periodic_bc_z(term_2,4)
      Call update_ghost_interior_planes(term_2,4)
  
      ! hat(ui)*hat(uj)
      SWij(:,2:3        ,:,1) = term  (2:nxg,1:2   ,2:nzg) * term  (2:nxg,1:2     ,2:nzg) ! hat(u)^2
      SWij(:,2:3        ,:,2) = Vff   (2:nxg,1:2   ,2:nzg) * Vff   (2:nxg,1:2     ,2:nzg) ! hat(v)^2
      SWij(:,2:3        ,:,3) = term_2(2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(w)^2
      SWij(:,2:3        ,:,4) = term  (2:nxg,1:2   ,2:nzg) * Vff   (2:nxg,1:2     ,2:nzg) ! hat(u)hat(v)
      SWij(:,2:3        ,:,5) = term  (2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(u)hat(w)
      SWij(:,2:3        ,:,6) = Vff   (2:nxg,1:2   ,2:nzg) * term_2(2:nxg,1:2     ,2:nzg) ! hat(v)hat(w)

      SWij(:,nyg-2:nyg-1,:,1) = term  (2:nxg,ny-1:ny,2:nzg) * term  (2:nxg,ny-1:ny,2:nzg) ! hat(u)^2
      SWij(:,nyg-2:nyg-1,:,2) = Vff   (2:nxg,ny-1:ny,2:nzg) * Vff   (2:nxg,ny-1:ny,2:nzg) ! hat(v)^2
      SWij(:,nyg-2:nyg-1,:,3) = term_2(2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(w)^2
      SWij(:,nyg-2:nyg-1,:,4) = term  (2:nxg,ny-1:ny,2:nzg) * Vff   (2:nxg,ny-1:ny,2:nzg) ! hat(u)hat(v)
      SWij(:,nyg-2:nyg-1,:,5) = term  (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(u)hat(w)
      SWij(:,nyg-2:nyg-1,:,6) = Vff   (2:nxg,ny-1:ny,2:nzg) * term_2(2:nxg,ny-1:ny,2:nzg) ! hat(v)hat(w)

      ! subtract it
      LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + SWij(2:nxg-1,    2,2:nzg-1,:)
      LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + SWij(2:nxg-1,nyg-1,2:nzg-1,:)


      ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
      Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
      Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
      Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
           
      ! apply periodicity in x and z
      Call apply_periodic_bc_x(Uff,1)
      Call apply_periodic_bc_z(Uff,1)
      Call update_ghost_interior_planes(Uff,1)
      Call apply_periodic_bc_x(Vff,2)
      Call apply_periodic_bc_z(Vff,2)
      Call update_ghost_interior_planes(Vff,2)
      Call apply_periodic_bc_x(Wff,2)
      Call apply_periodic_bc_z(Wff,3)
      Call update_ghost_interior_planes(Wff,3)
           
      ! compute eddy viscosity for filtered velocities
      Call compute_eddy_viscosity(Uff,Vff,Wff,avg_nu_t_hat,nu_t_hat)
           
      ! interpolate nu_t_hat from cell centers to cell faces (== cell edges because averaged in xz)
      Call interpolate_y(nu_t_hat,term_2(:,1:ny,:),2)    
           
      ! compute hat(SWij) at V locations in two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
      Call compute_Sij_at_V_location_wall(Uff,Vff,Wff,SWij)
           
      ! Tij = -2*nu_t_hat*hat(SWij) -> SWij
      Do i = 2, nxg
        Do k = 2, nzg
          Do j = 1, 6
            SWij(i,2:nyg,k,j) = -2d0*(nu+term_2(i,1:ny,k))*SWij(i,2:nyg,k,j) 
          End Do
        End Do
      End Do
           
      ! add it
      LWij(2:nxg-1,    2,2:nzg-1,:) = LWij(2:nxg-1,    2,2:nzg-1,:) + SWij(2:nxg-1,  2,2:nzg-1,:)
      LWij(2:nxg-1,nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-1,2:nzg-1,:) + SWij(2:nxg-1,nyg,2:nzg-1,:)

    Else
      Stop('Error: icase value not defined')
    End If
       
    ! Compute MWij
    MWij = 0d0
    If     (icase == 1) Then
      !--------------------Part 3------------------------------!
      ! Compute the stress for Delta(hat(tau_w)) 
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
      Call compute_Rij_at_V_location_wall(U_,V_,W_,Rij )
      S = 0d0
      Do i = 1,6
        If (i .le. 3) Then
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        Else 
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + 2d0 * SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + 2d0 * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        End If
      End Do
      S = S**0.5d0

      If ( a1/=0 ) Then
         Do i = 1,6
            MWij(2:nxg-1,2:4,2:nzg-1,i) = MWij(2:nxg-1,2:4,2:nzg-1,i) - Real(a1,8) * S(2:nxg-1,2:4,2:nzg-1) * SWij(2:nxg-1,2:4,2:nzg-1,i)
            MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) - Real(a1,8) * S(2:nxg-1,nyg-4:nyg-1,2:nzg-1)* SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
         End Do
      End If

      If ( a2/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a2,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If
  
      If ( a3/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a3,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a4/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a4,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a5/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               Do l = 1,3
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*SWij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
               End Do
             End Do
           End Do
         End Do
         Do i = 1,6
            MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) - Real(a5,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,i)/S(2:nxg-1,2:nyg-1,2:nzg-1)
         End Do
      End If

      ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,1:6) = MWij;
      Do i = 1,6
        Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
        Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
        Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
      End Do
      Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),MWij(2:nxg-1,2:nyg-1,2:nzg-1,:)) 

      !--------------------Part 4------------------------------!
      ! Compute the stress for Delta((tau_w)) 
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
      Call compute_Rij_at_V_location_wall(U_,V_,W_,Rij )
      S = 0d0
      Do i = 1,6
        If (i .le. 3) Then
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        Else 
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + 2d0 * SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + 2d0 * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        End If
      End Do
      S = S**0.5d0

      If ( a1/=0 ) Then
         Do i = 1,6
            MWij(2:nxg-1,2:4,2:nzg-1,i) = MWij(2:nxg-1,2:4,2:nzg-1,i) + Real(a1,8) * S(2:nxg-1,2:4,2:nzg-1) * SWij(2:nxg-1,2:4,2:nzg-1,i)
            MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) + Real(a1,8) * S(2:nxg-1,nyg-4:nyg-1,2:nzg-1)* SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
         End Do
      End If

      If ( a2/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a2,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If
  
      If ( a3/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a3,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a4/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a4,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a5/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               Do l = 1,3
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*SWij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
               End Do
             End Do
           End Do
         End Do
         Do i = 1,6
            MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) + Real(a5,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,i)/S(2:nxg-1,2:nyg-1,2:nzg-1)
         End Do
      End If

    Elseif (icase == 2) Then
      !--------------------Part 3------------------------------!
      ! Compute the stress for Delta((tau_w)) 
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
      Call compute_Rij_at_V_location_wall(U_,V_,W_,Rij )
      S = 0d0
      Do i = 1,6
        If (i .le. 3) Then
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        Else 
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + 2d0 * SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + 2d0 * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        End If
      End Do
      S = S**0.5d0

      If ( a1/=0 ) Then
         Do i = 1,6
            MWij(2:nxg-1,2:4,2:nzg-1,i) = MWij(2:nxg-1,2:4,2:nzg-1,i) + Real(a1,8) * S(2:nxg-1,2:4,2:nzg-1) * SWij(2:nxg-1,2:4,2:nzg-1,i)
            MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) + Real(a1,8) * S(2:nxg-1,nyg-4:nyg-1,2:nzg-1)* SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
         End Do
      End If

      If ( a2/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a2,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If
  
      If ( a3/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a3,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a4/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a4,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a5/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               Do l = 1,3
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*SWij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
               End Do
             End Do
           End Do
         End Do
         Do i = 1,6
            MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) + Real(a5,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,i)/S(2:nxg-1,2:nyg-1,2:nzg-1)
         End Do
      End If
      !--------------------Part 4------------------------------!
      ! Compute the stress for Delta(T_w) 

      ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
      Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
      Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
      Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
           
      ! apply periodicity in x and z
      Call apply_periodic_bc_x(Uff,1)
      Call apply_periodic_bc_z(Uff,1)
      Call update_ghost_interior_planes(Uff,1)
      Call apply_periodic_bc_x(Vff,2)
      Call apply_periodic_bc_z(Vff,2)
      Call update_ghost_interior_planes(Vff,2)
      Call apply_periodic_bc_x(Wff,2)
      Call apply_periodic_bc_z(Wff,3)
      Call update_ghost_interior_planes(Wff,3)
           
      Call compute_Sij_at_V_location_wall(Uff,Vff,Wff,SWij)
      Call compute_Rij_at_V_location_wall(Uff,Vff,Wff,Rij)     

      S = 0d0
      Do i = 1,6
        If (i .le. 3) Then
          S(2:nxg-1,2:5,2:nzg-1) = S(2:nxg-1,2:5,2:nzg-1) + SWij(2:nxg-1,2:5,2:nzg-1,i) * SWij(2:nxg-1,2:5,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        Else 
          S(2:nxg-1,2:5,2:nzg-1) = S(2:nxg-1,2:5,2:nzg-1) + 2d0 * SWij(2:nxg-1,2:5,2:nzg-1,i) * SWij(2:nxg-1,2:5,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + 2d0 * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        End If
      End Do
      S = S**0.5d0

      If ( a1/=0 ) Then
         Do i = 1,6
            MWij(2:nxg-1,2:5,2:nzg-1,i) = MWij(2:nxg-1,2:5,2:nzg-1,i) - Real(a1,8) * S(2:nxg-1,2:5,2:nzg-1) * SWij(2:nxg-1,2:5,2:nzg-1,i)
            MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) - Real(a1,8) * S(2:nxg-1,nyg-4:nyg-1,2:nzg-1)* SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
         End Do
      End If

      If ( a2/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a2,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If
  
      If ( a3/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a3,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a4/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a4,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a5/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               Do l = 1,3
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*SWij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
               End Do
             End Do
           End Do
         End Do
         Do i = 1,6
            MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) - Real(a5,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,i)/S(2:nxg-1,2:nyg-1,2:nzg-1)
         End Do
      End If

    Elseif (icase == 3) Then
      !--------------------Part 3------------------------------!
      ! Compute the stress for Delta(hat(tau_w)) 
      Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
      Call compute_Rij_at_V_location_wall(U_,V_,W_,Rij )
      S = 0d0
      Do i = 1,6
        If (i .le. 3) Then
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        Else 
          S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + 2d0 * SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + 2d0 * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        End If
      End Do
      S = S**0.5d0

      If ( a1/=0 ) Then
         Do i = 1,6
            MWij(2:nxg-1,2:4,2:nzg-1,i) = MWij(2:nxg-1,2:4,2:nzg-1,i) + Real(a1,8) * S(2:nxg-1,2:4,2:nzg-1) * SWij(2:nxg-1,2:4,2:nzg-1,i)
            MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) + Real(a1,8) * S(2:nxg-1,nyg-4:nyg-1,2:nzg-1)* SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
         End Do
      End If

      If ( a2/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a2,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If
  
      If ( a3/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a3,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a4/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Real(a4,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a5/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               Do l = 1,3
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*SWij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
               End Do
             End Do
           End Do
         End Do
         Do i = 1,6
            MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) + Real(a5,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,i)/S(2:nxg-1,2:nyg-1,2:nzg-1)
         End Do
      End If

      ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,1:6) = MWij;
      Do i = 1,6
        Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
        Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
        Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
      End Do
      Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),MWij(2:nxg-1,2:nyg-1,2:nzg-1,:)) 

      !--------------------Part 4------------------------------!
      ! Compute the stress for Delta(T_w) 

      ! filter velocities: use Uff because Uf is used in compute_eddy_viscosity
      Call filter_xzy( U_(1:nx, 1:nyg,1:nzg), Uff(2:nx-1 ,1:nyg,2:nzg-1) )
      Call filter_xzy( V_(1:nxg,1:ny, 1:nzg), Vff(2:nxg-1,1:ny ,2:nzg-1) )
      Call filter_xzy( W_(1:nxg,1:nyg,1:nz ), Wff(2:nxg-1,1:nyg,2:nz-1 ) )
           
      ! apply periodicity in x and z
      Call apply_periodic_bc_x(Uff,1)
      Call apply_periodic_bc_z(Uff,1)
      Call update_ghost_interior_planes(Uff,1)
      Call apply_periodic_bc_x(Vff,2)
      Call apply_periodic_bc_z(Vff,2)
      Call update_ghost_interior_planes(Vff,2)
      Call apply_periodic_bc_x(Wff,2)
      Call apply_periodic_bc_z(Wff,3)
      Call update_ghost_interior_planes(Wff,3)
           
      Call compute_Sij_at_V_location_wall(Uff,Vff,Wff,SWij)
      Call compute_Rij_at_V_location_wall(Uff,Vff,Wff,Rij)     

      S = 0d0
      Do i = 1,6
        If (i .le. 3) Then
          S(2:nxg-1,2:5,2:nzg-1) = S(2:nxg-1,2:5,2:nzg-1) + SWij(2:nxg-1,2:5,2:nzg-1,i) * SWij(2:nxg-1,2:5,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        Else 
          S(2:nxg-1,2:5,2:nzg-1) = S(2:nxg-1,2:5,2:nzg-1) + 2d0 * SWij(2:nxg-1,2:5,2:nzg-1,i) * SWij(2:nxg-1,2:5,2:nzg-1,i)
          S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + 2d0 * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
        End If
      End Do
      S = S**0.5d0

      If ( a1/=0 ) Then
         Do i = 1,6
            MWij(2:nxg-1,2:5,2:nzg-1,i) = MWij(2:nxg-1,2:5,2:nzg-1,i) - Real(a1,8) * S(2:nxg-1,2:5,2:nzg-1) * SWij(2:nxg-1,2:5,2:nzg-1,i)
            MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) - Real(a1,8) * S(2:nxg-1,nyg-4:nyg-1,2:nzg-1)* SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
         End Do
      End If

      If ( a2/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a2,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If
  
      If ( a3/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a3,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a4/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
               ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
             End Do
           End Do
         End Do
         MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) - Real(a4,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
      End If

      If ( a5/=0 ) Then
         ten_buf = 0d0
         Do i = 1,3
           Do j = 1,3
             Do k = 1,3
               Do l = 1,3
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*SWij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
                 ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
               End Do
             End Do
           End Do
         End Do
         Do i = 1,6
            MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) - Real(a5,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,i)/S(2:nxg-1,2:nyg-1,2:nzg-1)
         End Do
      End If
    Else
      Stop('Error: icase value not defined')
    End If

    ! Multiply MWij by -Delta^2
    MWij = - Delta**2d0 * MWij 

  End Subroutine compute_LWij_MWij_CW


  Subroutine compute_delta_tau(U_,V_,W_,nu_t_,Cw,ai)


    Real   (Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real   (Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real   (Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real   (Int64), Dimension(nxg,nyg,nzg), Intent(In) :: nu_t_
    Integer(Int32), Dimension(6),           Intent(In) :: ai
    Real   (Int64), Dimension(5),           Intent(In) :: Cw 

    ! local variables
    Integer(Int32) :: i, j, k, l, a1, a2, a3, a4, a5, ind(3,3) 
    Real   (Int64) :: r_sign(3,3), Delta

    Delta = ((y(2)-y(1)) * (x(2)-x(1)) * (z(2)-z(1)))**(1d0/3d0)

    a1 = ai(1)
    a2 = ai(2)
    a3 = ai(3)
    a4 = ai(4)
    a5 = ai(5)

    ind(1,1) = 1;
    ind(2,2) = 2;
    ind(3,3) = 3;
    ind(1,2) = 4;
    ind(2,1) = 4;
    ind(1,3) = 5;
    ind(3,1) = 5;
    ind(2,3) = 6;
    ind(3,2) = 6;

    r_sign(1,1) =  1d0;
    r_sign(2,2) =  1d0;
    r_sign(3,3) =  1d0;
    r_sign(1,2) =  1d0;
    r_sign(2,1) = -1d0;
    r_sign(1,3) =  1d0;
    r_sign(3,1) = -1d0;
    r_sign(2,3) =  1d0;
    r_sign(3,2) = -1d0;


    ! Compute MWij <-- tau 
    MWij = 0d0

    ! Compute the stress for Delta(hat(tau_w)) 
    Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
    Call compute_Rij_at_V_location_wall(U_,V_,W_,Rij )
    S = 0d0
    Do i = 1,6
      If (i .le. 3) Then
        S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
        S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
      Else 
        S(2:nxg-1,2:4,2:nzg-1) = S(2:nxg-1,2:4,2:nzg-1) + 2d0 * SWij(2:nxg-1,2:4,2:nzg-1,i) * SWij(2:nxg-1,2:4,2:nzg-1,i)
        S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) = S(2:nxg-1,nyg-4:nyg-1,2:nzg-1) + 2d0 * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i) * SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
      End If
    End Do
    S = S**0.5d0

    If ( a1/=0 ) Then
       Do i = 1,6
         MWij(2:nxg-1,2:4,2:nzg-1,i) = MWij(2:nxg-1,2:4,2:nzg-1,i) + Cw(1) * Real(a1,8) * S(2:nxg-1,2:4,2:nzg-1) * SWij(2:nxg-1,2:4,2:nzg-1,i)
         MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,i) + Cw(1) * Real(a1,8) * S(2:nxg-1,nyg-4:nyg-1,2:nzg-1)* SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,i)
       End Do
    End If

    If ( a2/=0 ) Then
      ten_buf = 0d0
      Do i = 1,3
        Do j = 1,3
          Do k = 1,3
            ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
            ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
          End Do
        End Do
      End Do
      MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Cw(2) * Real(a2,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
   End If
  
   If ( a3/=0 ) Then
      ten_buf = 0d0
      Do i = 1,3
        Do j = 1,3
          Do k = 1,3
            ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
            ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
          End Do
        End Do
      End Do
      MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Cw(3) * Real(a3,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
    End If

    If ( a4/=0 ) Then
      ten_buf = 0d0
      Do i = 1,3
        Do j = 1,3
          Do k = 1,3
            ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
            ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,j))
            ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*r_sign(k,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
            ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,j))
          End Do
        End Do
      End Do
      MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,:) + Cw(4) * Real(a4,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,:)
    End If

    If ( a5/=0 ) Then
      ten_buf = 0d0
      Do i = 1,3
        Do j = 1,3
          Do k = 1,3
            Do l = 1,3
              ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
              ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,2:5,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,2:5,2:nzg-1,ind(i,k))*SWij(2:nxg-1,2:5,2:nzg-1,ind(k,l))*SWij(2:nxg-1,2:5,2:nzg-1,ind(l,j))
              ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) + &
                    SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*r_sign(l,j)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
              ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) = ten_buf(2:nxg-1,nyg-4:nyg-1,2:nzg-1,ind(i,j)) - &
                    r_sign(i,k)*Rij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(i,k))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(k,l))*SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,ind(l,j))
            End Do
          End Do
        End Do
      End Do
      Do i = 1,6
        MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) = MWij(2:nxg-1,2:nyg-1,2:nzg-1,i) + Cw(5) * Real(a5,8) * ten_buf(2:nxg-1,2:nyg-1,2:nzg-1,i)/S(2:nxg-1,2:nyg-1,2:nzg-1)
      End Do
    End If


    LWij = Delta**2d0 * MWij 

    ! interpolate velocity to U and W to V location (at wall)
    Call interpolate_x(U_,     term_1(2:nx,:,:),1) 
    Call interpolate_y(term_1, term  (:,1:ny,:),2)
      
    Call interpolate_z(W_,     term_1(:,:,2:nz),1) 
    Call interpolate_y(term_1, term_2(:,1:ny,:),2) 

    ! fill in missing values (periodicity)
    Call apply_periodic_bc_x(term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! SWij = ui*uj at v location (at wall)
    ! U at V location -> term  (2:nx-1,1:ny,:)
    ! W at V location -> term_2(:,1:ny,2:nzg-1)
    ! bottom wall 
    SWij(:,2:5,:,1)         = term  (2:nxg,1:4,2:nzg) * term  (2:nxg,1:4,2:nzg)  ! u^2
    SWij(:,2:5,:,2)         = V_    (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! v^2
    SWij(:,2:5,:,3)         = term_2(2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! w^2
    SWij(:,2:5,:,4)         = term  (2:nxg,1:4,2:nzg) * V_    (2:nxg,1:4,2:nzg)  ! uv
    SWij(:,2:5,:,5)         = term  (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! uw
    SWij(:,2:5,:,6)         = V_    (2:nxg,1:4,2:nzg) * term_2(2:nxg,1:4,2:nzg)  ! vw
    ! top wall
    SWij(:,nyg-4:nyg-1,:,1) = term  (2:nxg,ny-3:ny,2:nzg) * term  (2:nxg,ny-3:ny,2:nzg) ! u^2
    SWij(:,nyg-4:nyg-1,:,2) = V_    (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! v^2
    SWij(:,nyg-4:nyg-1,:,3) = term_2(2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! w^2
    SWij(:,nyg-4:nyg-1,:,4) = term  (2:nxg,ny-3:ny,2:nzg) * V_    (2:nxg,ny-3:ny,2:nzg) ! uv
    SWij(:,nyg-4:nyg-1,:,5) = term  (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! uw
    SWij(:,nyg-4:nyg-1,:,6) = V_    (2:nxg,ny-3:ny,2:nzg) * term_2(2:nxg,ny-3:ny,2:nzg) ! vw
      
    ! ui*uj
    LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,        2:5,2:nzg-1,:)
    LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:)

    !----------compute tau_ij, hat(tau_ij) and hat(hat(tau_ij))   
    ! compute SWij at V locations in the first two cells-> SWij(2:nxg,2:3,2:nzg) (bottom) and SWij(2:nxg,nyg-1:nyg,2:nzg) (top)
    Call compute_Sij_at_V_location_wall(U_,V_,W_,SWij)
           
    ! interpolate nu_t from cell centers to cell faces (== cell edges because averaged in xz)
    Call interpolate_y(nu_t_,term_2(:,1:ny,:),2)
            
    ! tau_ij = -2*nu_t*SWij -> SWij
    Do i = 2, nxg
      Do k = 2, nzg
        Do j = 1, 6
          SWij(i,2:nyg,k,j) = -2d0*(nu+term_2(i,1:ny,k))*SWij(i,2:nyg,k,j)
        End Do
      End Do
    End Do
           
    ! add it
    LWij(2:nxg-1,        2:5,2:nzg-1,:) = LWij(2:nxg-1,        2:5,2:nzg-1,:) - SWij(2:nxg-1,      2:5,2:nzg-1,:)
    LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) = LWij(2:nxg-1,nyg-4:nyg-1,2:nzg-1,:) - SWij(2:nxg-1,nyg-3:nyg,2:nzg-1,:)
  End Subroutine compute_delta_tau


End Module wallmodel

