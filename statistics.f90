!--------------------------------------------!
! Module for computing some basic statistics !
!--------------------------------------------!
Module statistics

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use interpolation
  Use input_output
  Use mass_flow

  ! prevent implicit typing
  Implicit None

Contains

  !------------------------------------------!
  ! Compute some basic statistics on the fly !
  !------------------------------------------!
  Subroutine compute_statistics    

    Integer(Int32) :: jj
    Real   (Int64) :: dUdy_wall_b, dUdy_wall_t
    Real   (Int64) ::   UV_wall_b,   UV_wall_t
    Real   (Int64) :: U1, U2, Un1, Un2
    Real   (Int64) :: U1_local, U2_local, Un1_local, Un2_local
    Real   (Int64) :: nu_t_avg_b_local, nu_t_avg_t_local
    Real   (Int64) :: nu_t_avg_b, nu_t_avg_t

    ! if pressure not computed 
    pressure_computed = .False.

    ! statistics computed at grid y -> U and W interpolated    
    if ( Mod(istep,nstats)==0 .Or. istep==1 ) Then

       ! Compute actual pressure (should be called first, uses term_1,...)
       !Call compute_pressure       
       ! now computed in projection.f90
       pressure_computed = .True.
       
       ! interpolate U in y -> term_1
       Call interpolate_y(U,term_1,2)

       ! interpolate W in y -> term_2
       Call interpolate_y(W,term_2,2)

       ! interpolate U in x (already in y) -> term
       Call interpolate_x(term_1,term,1)

       ! interpolate P in y -> term_3
      ! Call interpolate_y(P,term_3,2)

       ! compute local statistics
       Do jj=1,ny

          Umean  (jj) = Sum(term_1 (1:nx-2, jj,1:nzg-2) )
          Wmean  (jj) = Sum(term_2 (1:nxg-2,jj,1:nz-2 ) )
          Vmean  (jj) = Sum(     V (2:nxg-1,jj,2:nzg-1) )

          U2mean (jj) = Sum(term_1 (1:nx-2, jj,1:nzg-2)**2d0 )
          W2mean (jj) = Sum(term_2 (1:nxg-2,jj,1:nz-2 )**2d0 )
          V2mean (jj) = Sum(     V (2:nxg-1,jj,2:nzg-1)**2d0 )

          UVmean (jj) = Sum( V(2:nxg-1,jj,2:nzg-1)*term(1:nxg-2,jj,2:nzg-1) )

          Pmean  (jj) = Sum( term_3 (1:nxg-2, jj, 1:nzg-2)      )
          P2mean (jj) = Sum( term_3 (1:nxg-2, jj, 1:nzg-2)**2d0 )

       End Do

       ! reduce statatistics between processors      
       IF ( myid==0 ) Then

          Call MPI_Reduce(MPI_IN_PLACE,Umean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,Vmean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,Wmean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(MPI_IN_PLACE,U2mean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,V2mean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,W2mean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(MPI_IN_PLACE,UVmean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

          Call MPI_Reduce(MPI_IN_PLACE,Pmean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,P2mean,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
                    
       Else

          Call MPI_Reduce(Umean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(Vmean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(Wmean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(U2mean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(V2mean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(W2mean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(UVmean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

          Call MPI_Reduce(Pmean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(P2mean,0,ny,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
       End IF

       ! These statistics are only good for processor 0
       Umean  = Umean/Real( ( nx_global-2)*(nzg_global-2), 8)
       Vmean  = Vmean/Real( (nxg_global-2)*(nzg_global-2), 8)
       Wmean  = Wmean/Real( (nxg_global-2)*( nz_global-2), 8)

       U2mean = U2mean/Real( ( nx_global-2)*(nzg_global-2), 8)
       V2mean = V2mean/Real( (nxg_global-2)*(nzg_global-2), 8)
       W2mean = W2mean/Real( (nxg_global-2)*( nz_global-2), 8)

       UVmean = UVmean/Real( (nxg_global-2)*(nzg_global-2) ,8)

       Pmean  = UVmean/Real( (nxg_global-2)*(nzg_global-2) ,8)
       P2mean = UVmean/Real( (nxg_global-2)*(nzg_global-2) ,8)

       ! Mean derivative at the walls (CHECK THIS PLEASE)
       U1_local  = Sum( U(2:nx-1,1,2:nzg-1)     )
       U2_local  = Sum( U(2:nx-1,2,2:nzg-1)     )
       Un1_local = Sum( U(2:nx-1,nyg-1,2:nzg-1) )
       Un2_local = Sum( U(2:nx-1,nyg  ,2:nzg-1) )
      
       Call MPI_Reduce( U1_local, U1,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce( U2_local, U2,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(Un1_local,Un1,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(Un2_local,Un2,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       U1  = U1 /Real( (nx-2)*(nzg_global-2), 8 )
       U2  = U2 /Real( (nx-2)*(nzg_global-2), 8 )
       Un1 = Un1/Real( (nx-2)*(nzg_global-2), 8 )
       Un2 = Un2/Real( (nx-2)*(nzg_global-2), 8 )
       
       dUdy_wall_b = ( U2 -  U1)/( yg(  2) - yg(    1)) 
       dUdy_wall_t = (Un2 - Un1)/( yg(nyg) - yg(nyg-1))

       ! Mean eddy viscosity (this is approx, rigorously <2*nu_t*Sij> )
       nu_t_avg_b_local = Sum( nu_t(2:nxg-1,    1:2  ,2:nzg-1) )
       nu_t_avg_t_local = Sum( nu_t(2:nxg-1,nyg-1:nyg,2:nzg-1) )
       Call MPI_Reduce(nu_t_avg_b_local,nu_t_avg_b,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Call MPI_Reduce(nu_t_avg_t_local,nu_t_avg_t,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       nu_t_avg_b = nu_t_avg_b/Real( (nxg-2)*(nzg_global-2), 8 )/2d0
       nu_t_avg_t = nu_t_avg_t/Real( (nxg-2)*(nzg_global-2), 8 )/2d0

       ! Mean Reynolds stress at the walls
       UV_wall_b = UVmean( 1) 
       UV_wall_t = UVmean(ny)

       ! friction velocity
       utau = ( ( UV_wall_t - UV_wall_b - (nu+nu_t_avg_t)*dUdy_wall_t + (nu+nu_t_avg_b)*dUdy_wall_b )/Ly )**0.5d0

       ! friction Reynolds number
       Retau = utau*(y(ny)-y(1))/2d0/nu

       ! mean mass flow in x
       Call compute_mean_mass_flow_U(U,Qflow_x)
       Call compute_mean_mass_flow_V(V,Qflow_y)

       ! write statistics
       Call output_statistics

       ! Sanity check 
       If ( Any( Isnan(U) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(V) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(W) ) ) Stop 'Error: NaNs!'
       
    End If
    
  End Subroutine compute_statistics

End Module statistics
