module cabaret_solver_class

	use kind_parameters
	use global_data
	use data_manager_class
	use computational_domain_class
	use computational_mesh_class
	use boundary_conditions_class
	use field_pointers
	use table_approximated_real_gas_class
	use thermophysical_properties_class	
	use chemical_kinetics_solver_class
	use chemical_properties_class
	use viscosity_solver_class
	use large_particles_method
	use fourier_heat_transfer_solver_class
	use fickean_diffusion_solver_class

	implicit none

	include "omp_lib.h"

	private
	public	:: cabaret_solver, cabaret_solver_c

	type(field_scalar_flow)	,target	::	rho_f_new, p_f_new, e_i_f_new, v_s_f_new, E_f_f_new, T_f_new
	type(field_vector_flow)	,target	::	Y_f_new, v_f_new
	
	type cabaret_solver
		logical			:: diffusion_flag, viscosity_flag, heat_trans_flag, reactive_flag, sources_flag, gas_dynamics_flag
		real(dkind)		:: courant_fraction
		real(dkind)		:: time, time_step
		
		type(chemical_kinetics_solver)		:: chem_kin_solver
		type(diffusion_solver)				:: diff_solver
		type(heat_transfer_solver)			:: heat_trans_solver
		type(viscosity_solver)				:: viscosity_solver
		type(table_approximated_real_gas)	:: state_eq

		type(computational_domain)				:: domain
		type(chemical_properties_pointer)		:: chem
		type(thermophysical_properties_pointer)	:: thermo
		type(computational_mesh_pointer)		:: mesh
		type(boundary_conditions_pointer)		:: boundary

		type(field_scalar_cons_pointer)	:: rho	, T	, p	, v_s, gamma, E_f	, e_i ,mol_mix_conc
		type(field_scalar_flow_pointer)	:: gamma_f_new, rho_f_new, p_f_new, e_i_f_new, v_s_f_new, E_f_f_new, T_f_new
		
		type(field_scalar_cons_pointer)	:: E_f_prod_chem, E_f_prod_heat, E_f_prod_visc

		type(field_vector_cons_pointer)	:: v, Y	, v_prod_visc
		type(field_vector_flow_pointer)	:: v_f_new, Y_f_new
		
		type(field_vector_cons_pointer)	:: Y_prod_chem, Y_prod_diff
		
		! Conservative variables
		real(dkind) ,dimension(:,:,:)	,allocatable    :: rho_old, p_old, E_f_old, e_i_old, E_f_prod, rho_prod, v_s_old, gamma_old
		real(dkind)	,dimension(:,:,:,:)	,allocatable	:: v_old, Y_old, v_prod, Y_prod

		! Flow variables
		real(dkind) ,dimension(:,:,:,:,:)	,allocatable    :: v_f, Y_f
		real(dkind) ,dimension(:,:,:,:)		,allocatable    :: rho_f, p_f, e_i_f, E_f_f, v_s_f
		! Quasi invariants

	contains
		procedure	,private	:: apply_boundary_conditions_main
		procedure				:: solve_problem
		procedure				:: calculate_time_step
		procedure				:: get_time_step
		procedure				:: get_time
		procedure				:: check_symmetry
		procedure				:: check_symmetry_2
	end type

	interface	cabaret_solver_c
		module procedure	constructor
	end interface

contains

	type(cabaret_solver)	function constructor(manager,calculation_time, diffusion_flag, viscosity_flag, heat_trans_flag, reactive_flag, sources_flag, gas_dynamics_flag, courant_fraction)
		type(data_manager)						,intent(inout)	:: manager
		real(dkind)								,intent(in)		:: calculation_time
		logical									,intent(in)		:: diffusion_flag
		logical									,intent(in)		:: viscosity_flag
		logical									,intent(in)		:: heat_trans_flag
		logical									,intent(in)		:: reactive_flag
		logical									,intent(in)		:: sources_flag
		logical									,intent(in)		:: gas_dynamics_flag
		real(dkind)								,intent(in)		:: courant_fraction

		type(field_scalar_cons_pointer)	:: scal_ptr
		type(field_vector_cons_pointer)	:: vect_ptr
		type(field_tensor_cons_pointer)	:: tens_ptr		

		type(field_scalar_flow_pointer)	:: scal_f_ptr		
		type(field_vector_flow_pointer)	:: vect_f_ptr
		
		
		real(dkind)				:: spec_summ
		integer	,dimension(3,2)	:: loop_bound
		integer					:: dimensions, species_number
		integer					:: i, j, k, dim, spec

		loop_bound		= manager%domain%long_cycle
		dimensions		= manager%domain%dimensions
		species_number	= manager%chemistry%chem_ptr%species_number
		
		constructor%diffusion_flag		= diffusion_flag
		constructor%viscosity_flag		= viscosity_flag
		constructor%heat_trans_flag		= heat_trans_flag
		constructor%reactive_flag		= reactive_flag
		constructor%sources_flag		= sources_flag
		constructor%gas_dynamics_flag	= gas_dynamics_flag
		constructor%courant_fraction	= courant_fraction

		constructor%domain				= manager%domain
		constructor%chem%chem_ptr		=> manager%chemistry%chem_ptr
		constructor%boundary%bc_ptr		=> manager%boundary_conditions_pointer%bc_ptr
		constructor%mesh%mesh_ptr		=> manager%computational_mesh_pointer%mesh_ptr
		constructor%thermo%thermo_ptr	=> manager%thermophysics%thermo_ptr

		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'density')
		constructor%rho%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'temperature')
		constructor%T%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'pressure')
		constructor%p%s_ptr					=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'full_energy')
		constructor%E_f%s_ptr				=> scal_ptr%s_ptr
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'internal_energy')
		constructor%e_i%s_ptr				=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'mixture_molar_concentration')
		constructor%mol_mix_conc%s_ptr		=> scal_ptr%s_ptr		
		
		call manager%create_scalar_field(rho_f_new	,'density_flow'				,'rho_f_new')
		constructor%rho_f_new%s_ptr => rho_f_new		
		call manager%create_scalar_field(p_f_new	,'pressure_flow'			,'p_f_new')
		constructor%p_f_new%s_ptr => p_f_new
		call manager%create_scalar_field(e_i_f_new	,'internal_energy_flow'		,'e_i_f_new')
		constructor%e_i_f_new%s_ptr => e_i_f_new
		call manager%create_scalar_field(E_f_f_new	,'full_energy_flow'			,'E_f_f_new')
		constructor%E_f_f_new%s_ptr => E_f_f_new
		call manager%create_scalar_field(v_s_f_new	,'velocity_of_sound_flow'	,'v_s_f_new')
		constructor%v_s_f_new%s_ptr => v_s_f_new
		call manager%create_scalar_field(T_f_new	,'temperature_flow'			,'T_f_new')
		constructor%T_f_new%s_ptr => T_f_new
		
		call manager%create_vector_field(Y_f_new,'specie_molar_concentration_flow'	,'Y_f_new',	'chemical')
		constructor%Y_f_new%v_ptr => Y_f_new		
		call manager%create_vector_field(v_f_new,'velocity_flow'					,'v_f_new',	'spatial')
		constructor%v_f_new%v_ptr => v_f_new	
		
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'velocity')
		constructor%v%v_ptr				=> vect_ptr%v_ptr		
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'specie_molar_concentration')
		constructor%Y%v_ptr				=> vect_ptr%v_ptr
		
		constructor%state_eq	=	table_approximated_real_gas_c(manager)

		call manager%get_flow_field_pointer(scal_f_ptr,vect_f_ptr,'adiabatic_index_flow')
		constructor%gamma_f_new%s_ptr	=> scal_f_ptr%s_ptr				
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'adiabatic_index')
		constructor%gamma%s_ptr			=> scal_ptr%s_ptr		
		call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'velocity_of_sound')
		constructor%v_s%s_ptr			=> scal_ptr%s_ptr
		
		if (reactive_flag) then
			constructor%chem_kin_solver		= chemical_kinetics_solver_c(manager)
			call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'energy_production_chemistry')
			constructor%E_f_prod_chem%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'specie_production_chemistry')
			constructor%Y_prod_chem%v_ptr			=> vect_ptr%v_ptr
		end if
		
		if (diffusion_flag) then
			constructor%diff_solver			= diffusion_solver_c(manager)
			call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'specie_production_diffusion')
			constructor%Y_prod_diff%v_ptr			=> vect_ptr%v_ptr
		end if
		
		if (heat_trans_flag) then
			constructor%heat_trans_solver	= heat_transfer_solver_c(manager)
			call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'energy_production_heat_transfer')
			constructor%E_f_prod_heat%s_ptr			=> scal_ptr%s_ptr
		end if		

		if(viscosity_flag) then
			constructor%viscosity_solver			= viscosity_solver_c(manager)
			call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'energy_production_viscosity')
			constructor%E_f_prod_visc%s_ptr			=> scal_ptr%s_ptr
			call manager%get_cons_field_pointer(scal_ptr,vect_ptr,tens_ptr,'velocity_production_viscosity')
			constructor%v_prod_visc%v_ptr			=> vect_ptr%v_ptr
		end if		
		
		allocate(constructor%rho_old(	loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))
										
		allocate(constructor%p_old(		loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))
										
		allocate(constructor%E_f_old(	loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))	
										
		allocate(constructor%e_i_old(	loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))
										
		allocate(constructor%E_f_prod(	loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))		
										
		allocate(constructor%rho_prod(	loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))
										
		allocate(constructor%v_s_old(	loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))										
										
		allocate(constructor%gamma_old(	loop_bound(1,1):loop_bound(1,2), &
										loop_bound(2,1):loop_bound(2,2), &
										loop_bound(3,1):loop_bound(3,2)))
	
		allocate(constructor%v_old(		dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))										

		allocate(constructor%v_prod(	dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))										
										
										
		allocate(constructor%Y_old(		species_number					, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))	
										
		allocate(constructor%Y_prod(	species_number					, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))
			
										
		
		allocate(constructor%rho_f(		dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))
										
		allocate(constructor%p_f(		dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))
										
		allocate(constructor%E_f_f(		dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))	
										
		allocate(constructor%e_i_f(		dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))
										
		allocate(constructor%v_s_f(		dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))	
										
		allocate(constructor%v_f(		dimensions						, &
										dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))	
										
										
		allocate(constructor%Y_f(		species_number					, &
										dimensions						, &
										loop_bound(1,1):loop_bound(1,2)	, &
										loop_bound(2,1):loop_bound(2,2)	, &
										loop_bound(3,1):loop_bound(3,2)))		
							
										
		loop_bound(:,1) = manager%domain%short_cycle(:,1) 
		loop_bound(:,2) = manager%domain%long_cycle(:,2)
										
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)
			do dim = 1, dimensions
				if (( (i*I_m(dim,1) +j*I_m(dim,2) +k*I_m(dim,3)) > loop_bound(dim,1)) .and. ( (i*I_m(dim,1) +j*I_m(dim,2) +k*I_m(dim,3)) < loop_bound(dim,2))) then
					constructor%p_f(dim,i,j,k)		=	0.5_dkind * (constructor%p%s_ptr%cells(i,j,k)	+ constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
					constructor%rho_f(dim,i,j,k)	=	0.5_dkind * (constructor%rho%s_ptr%cells(i,j,k)	+ constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
					spec_summ = 0.0_dkind
					do spec = 1, species_number
						constructor%Y_f(spec,dim,i,j,k) = 0.5_dkind * (constructor%Y%v_ptr%pr(spec)%cells(i,j,k) + constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
						spec_summ = spec_summ + max(constructor%Y_f(spec,dim,i,j,k), 0.0_dkind)
					end do

					do spec = 1,species_number
						constructor%Y_f(spec,dim,i,j,k) = max(constructor%Y_f(spec,dim,i,j,k), 0.0_dkind) / spec_summ
					end do
				end if
				
				if ((i*I_m(dim,1) +j*I_m(dim,2) +k*I_m(dim,3)) == loop_bound(dim,1)) then
					constructor%p_f(dim,i,j,k) = constructor%p%s_ptr%cells(i,j,k)
					constructor%rho_f(dim,i,j,k) = constructor%rho%s_ptr%cells(i,j,k)
					do spec = 1, species_number
						constructor%Y_f(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i,j,k)
					end do					
				end if
				if ((i*I_m(dim,1) +j*I_m(dim,2) +k*I_m(dim,3)) == loop_bound(dim,2)) then
					constructor%p_f(dim,i,j,k) = constructor%p%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
					constructor%rho_f(dim,i,j,k) = constructor%rho%s_ptr%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
					do spec = 1, species_number
						constructor%Y_f(spec,dim,i,j,k) = constructor%Y%v_ptr%pr(spec)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))
					end do
				end if
			end do
        end do										
		end do	
		end do	
		
		constructor%E_f_f	= 0.0_dkind
		constructor%e_i_f	= 0.0_dkind
		constructor%v_s_f	= 0.0_dkind
		constructor%v_f		= 0.0_dkind
		
		constructor%time		= calculation_time
		
		call constructor%apply_boundary_conditions_main()
		
	end function

	subroutine solve_problem(this)
		class(cabaret_solver)	,intent(inout)	:: this
    
		real(dkind)	,dimension(2)	:: r, q, r_corrected, q_corrected, r_new, q_new
		real(dkind)	,dimension(this%domain%dimensions,2)			:: v_inv, v_inv_corrected, v_inv_new
		real(dkind)	,dimension(this%domain%dimensions)				:: v_inv_half, v_inv_old
		real(dkind)	,dimension(this%chem%chem_ptr%species_number,2)	:: Y_inv, Y_inv_corrected, y_inv_new
		real(dkind)	,dimension(this%chem%chem_ptr%species_number)	:: Y_inv_half, Y_inv_old
		real(dkind)					:: r_half, q_half, R_old, Q_old
		real(dkind)					:: G_half, G_half_lower, G_half_higher
		
		real(dkind)	,dimension(3)	:: characteristic_speed
		real(dkind)	:: v_f_approx, v_s_f_approx
		real(dkind)	:: v_f_approx_lower, v_f_approx_higher

		real(dkind)	:: g_inv
		real(dkind)	:: f, corr, alpha = 0.0
		real(dkind) :: max_inv, min_inv
		real(dkind) :: mean_higher, mean_lower
		real(dkind)	:: sources
		real(dkind)	:: mean_sources
		real(dkind)	:: summ_frac
		real(dkind)	:: energy_source = 180000000.0_dkind
		real(dkind)	,save	:: energy_output = 0.0_dkind
		
		!real(dkind)	:: r_inf, q_inf, s_inf
		real(dkind)	,parameter	:: u_inf		= 0.0
		real(dkind)	,parameter	:: p_inf		= 101325.000000000_dkind
		real(dkind)	,parameter	:: rho_inf		= 0.487471044003493_dkind
		real(dkind)	,parameter	:: c_inf		= 539.709011600784_dkind
		real(dkind)	,parameter	:: gamma_inf	= 1.40137_dkind
		real(dkind)			:: g_inf		= 1.0_dkind / rho_inf / c_inf
		real(dkind)	,save	:: q_inf		= u_inf - p_inf / rho_inf / c_inf
		
		
		real(dkind)	:: spec_summ, rho_Y
		
		integer	,dimension(3,2)	:: loop_bound
		integer :: i,j,k,plus,dim,dim1,dim2,spec,iter		
		logical	:: flag
		
		call this%calculate_time_step()

		this%time = this%time + this%time_step		

		associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					e_i			=> this%e_i%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
					gamma		=> this%gamma%s_ptr		, &
					
					v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &
					
					v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					E_f_f		=> this%E_f_f	, &
					e_i_f		=> this%e_i_f	, &
					p_f			=> this%p_f		, &
					v_s_f		=> this%v_s_f	, & 
					Y_f			=> this%Y_f		, &
					
					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					E_f_old		=> this%E_f_old	,&
					Y_old		=> this%Y_old	,&
					p_old		=> this%p_old	,&
					v_s_old		=> this%v_s_old	,&
					gamma_old	=> this%gamma_old	, &
					
					p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					v_s_f_new	=> this%v_s_f_new%s_ptr		, &
					E_f_f_new	=> this%E_f_f_new%s_ptr		, &
					e_i_f_new	=> this%e_i_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr		, &
					gamma_f_new	=> this%gamma_f_new%s_ptr	, &
					
					v_prod_visc		=> this%v_prod_visc%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					E_f_prod_chem 	=> this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat	=> this%E_f_prod_heat%s_ptr	, &
					E_f_prod_visc	=> this%E_f_prod_visc%s_ptr	, &
					E_f_prod		=> this%E_f_prod			, &
					v_prod			=> this%v_prod				, &
					Y_prod			=> this%Y_prod				, &
					rho_prod		=> this%rho_prod			, &

					mesh		=> this%mesh%mesh_ptr)
										
		call this%apply_boundary_conditions_main()						

		!g_inf		= 1.0_dkind / rho_inf / c_inf 
		
		if (this%heat_trans_flag)	call this%heat_trans_solver%solve_heat_transfer(this%time_step)
		if (this%diffusion_flag)	call this%diff_solver%solve_diffusion(this%time_step)
		if (this%viscosity_flag)	call this%viscosity_solver%solve_viscosity(this%time_step)

		loop_bound		= this%domain%short_cycle
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)
		
			E_f_prod(i,j,k) = 0.0_dkind
			rho_prod(i,j,k) = 0.0_dkind
			Y_prod(:,i,j,k) = 0.0_dkind
			v_prod(:,i,j,k)	= 0.0_dkind	
			
			if (this%heat_trans_flag)	E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_heat%cells(i,j,k) 

			do spec = 1,this%chem%chem_ptr%species_number
				if (this%diffusion_flag)	Y_prod(spec,i,j,k)	= Y_prod(spec,i,j,k)	+ Y_prod_diff%pr(spec)%cells(i,j,k)
			end do
			
			if (this%viscosity_flag)	E_f_prod(i,j,k) = E_f_prod(i,j,k) + E_f_prod_visc%cells(i,j,k) 
			if (this%viscosity_flag)	then
				do dim = 1, this%domain%dimensions
					v_prod(dim,i,j,k)	= v_prod(dim,i,j,k) + v_prod_visc%pr(dim)%cells(i,j,k) 
				end do
			end if
		end do	
		end do
		end do
	
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)
			rho_old(i,j,k)			=	rho%cells(i,j,k)
			E_f_old(i,j,k)			=	E_f%cells(i,j,k)	
			p_old(i,j,k)			=	p%cells(i,j,k)
			v_s_old(i,j,k)			=	v_s%cells(i,j,k)
			gamma_old(i,j,k)		=	gamma%cells(i,j,k)

			rho%cells(i,j,k)		=	0.0_dkind 
			do dim = 1,this%domain%dimensions
				rho%cells(i,j,k)	=	rho%cells(i,j,k)	- (	rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) -	rho_f(dim,i,j,k)*v_f(dim,dim,i,j,k))/mesh%cell_size(1)
			end do
			rho%cells(i,j,k)		=	rho_old(i,j,k)	+  0.5_dkind*rho_prod(i,j,k) + 0.5_dkind*this%time_step * rho%cells(i,j,k)
			
			spec_summ = 0.0_dkind
			do spec = 1,this%chem%chem_ptr%species_number
				Y_old(spec,i,j,k)			=	Y%pr(spec)%cells(i,j,k)
				Y%pr(spec)%cells(i,j,k)		=	0.0_dkind !		
				do	dim = 1,this%domain%dimensions
					Y%pr(spec)%cells(i,j,k)	=  Y%pr(spec)%cells(i,j,k)	-	(		rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*Y_f(spec,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	&
																				-	rho_f(dim,i,j,k)*Y_f(spec,dim,i,j,k)*v_f(dim,dim,i,j,k))/mesh%cell_size(1)
				end do	
				Y%pr(spec)%cells(i,j,k)		=   rho_old(i,j,k)  *	Y_old(spec,i,j,k) + 0.5_dkind*Y_prod(spec,i,j,k) + 0.5_dkind*this%time_step * Y%pr(spec)%cells(i,j,k)
				Y%pr(spec)%cells(i,j,k)		=  Y%pr(spec)%cells(i,j,k)/rho%cells(i,j,k)
	
				spec_summ = spec_summ + max(Y%pr(spec)%cells(i,j,k), 0.0_dkind)
			end do

			do spec = 1,this%chem%chem_ptr%species_number
				y%pr(spec)%cells(i,j,k) = max(y%pr(spec)%cells(i,j,k), 0.0_dkind) / spec_summ 
			end do

			do dim = 1,this%domain%dimensions
				v_old(dim,i,j,k)			=	v%pr(dim)%cells(i,j,k)
				v%pr(dim)%cells(i,j,k)		=	0.0_dkind 
				do dim1 = 1,this%domain%dimensions
					v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(		rho_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))*v_f(dim,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))*v_f(dim1,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))	&
																				-	rho_f(dim1,i,j,k)*v_f(dim,dim1,i,j,k)*v_f(dim1,dim1,i,j,k)) /mesh%cell_size(1)
				end do
				v%pr(dim)%cells(i,j,k)	=	v%pr(dim)%cells(i,j,k)	-	 ( p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - p_f(dim,i,j,k)) /mesh%cell_size(1)
				v%pr(dim)%cells(i,j,k)	=	rho_old(i,j,k)*v_old(dim,i,j,k) + 0.5_dkind * v_prod(dim,i,j,k) + 0.5_dkind*this%time_step*v%pr(dim)%cells(i,j,k)
				v%pr(dim)%cells(i,j,k)	=	v%pr(dim)%cells(i,j,k) /rho%cells(i,j,k)
			end do	

			E_f%cells(i,j,k)		=	0.0_dkind
			do dim = 1,this%domain%dimensions
				E_f%cells(i,j,k)			= 	E_f%cells(i,j,k)		-	((rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*E_f_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+	p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))*v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))															&
																		-	(rho_f(dim,i,j,k)*E_f_f(dim,i,j,k)																		+	p_f(dim,i,j,k))									*v_f(dim,dim,i,j,k))/mesh%cell_size(1)
			end do	
			E_f%cells(i,j,k)		=	rho_old(i,j,k) * E_f_old(i,j,k)  + 0.5_dkind * E_f_prod(i,j,k)  + 0.5_dkind*this%time_step* E_f%cells(i,j,k)
			E_f%cells(i,j,k)		=	E_f%cells(i,j,k) / rho%cells(i,j,k) 
			
		end do
		end do
		end do
		
		! ********** Conservative variables state eq *******************	
		call this%state_eq%apply_state_equation() 
		! **************************************************************
		
		loop_bound		= this%domain%short_cycle
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)		
		
			do dim = 1,this%domain%dimensions
			
				G_half			= 1.0_dkind / (v_s%cells(i,j,k)*rho%cells(i,j,k))
				
				if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,1) ) then
					G_half_lower	= 1.0_dkind / (v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))*rho%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))
				end if
				
				if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,2) ) then
					G_half_higher	= 1.0_dkind / (v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))*rho%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
				end if
			
				! *********** Riemann quasi invariants *************************
				r(1) = v_f(dim,dim,i,j,k)									+ G_half*p_f(dim,i,j,k)
				r(2) = v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+ G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
				q(1) = v_f(dim,dim,i,j,k)									- G_half*p_f(dim,i,j,k)
				q(2) = v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))  - G_half*p_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
              
				R_half   = v%pr(dim)%cells(i,j,k)	+ G_half*(p%cells(i,j,k))
				Q_half   = v%pr(dim)%cells(i,j,k)	- G_half*(p%cells(i,j,k))
			
				do dim1 = 1,this%domain%dimensions
					if (dim1 == dim) then
						v_inv(dim1,1)		= p_f(dim1,i,j,k)										- v_s%cells(i,j,k)**2*rho_f(dim1,i,j,k)
						v_inv(dim1,2)		= p_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3))	- v_s%cells(i,j,k)**2*rho_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) 
						v_inv_half(dim1)	= p%cells(i,j,k)										- v_s%cells(i,j,k)**2*rho%cells(i,j,k)
					else
						v_inv(dim1,1)		= v_f(dim1,dim,i,j,k)
						v_inv(dim1,2)		= v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
						v_inv_half(dim1)	= v%pr(dim1)%cells(i,j,k)
					end if
				end do
				
				R_old   = v_old(dim,i,j,k)	+ G_half*p_old(i,j,k)
				Q_old   = v_old(dim,i,j,k)	- G_half*p_old(i,j,k)
			
				do dim1 = 1,this%domain%dimensions
					v_inv_old(dim)	=  p_old(i,j,k)	- v_s_old(i,j,k)**2*rho_old(i,j,k) 
				end do

				! ******************* Linear interpolation *********************
				
				r_new(1) = 2.0_dkind*R_half - r(2)
				q_new(1) = 2.0_dkind*Q_half - q(2)
				
				r_new(2) = 2.0_dkind*R_half - r(1)
				q_new(2) = 2.0_dkind*Q_half - q(1)				
			
				do dim1 = 1,this%domain%dimensions
					v_inv_new(dim1,1)	= 2.0_dkind*v_inv_half(dim1) - v_inv(dim1,2)
					v_inv_new(dim1,2)	= 2.0_dkind*v_inv_half(dim1) - v_inv(dim1,1)
				end do

				 ! **************** Non-linear flow correction ******************
	
				f = 0.0_dkind
				do dim1 = 1,this%domain%dimensions
					f = f +	v%pr(dim1)%cells(i,j,k) * v%pr(dim1)%cells(i,j,k)
				end do
				f = (0.25_dkind * rho_prod(i,j,k) * (E_f%cells(i,j,k) - f + p%cells(i,j,k) / rho%cells(i,j,k)) - 0.25_dkind * E_f_prod(i,j,k))
				do dim1 = 1,this%domain%dimensions
					f = f +	v%pr(dim1)%cells(i,j,k) * v_prod(dim1,i,j,k)
				end do				
				f = f * (gamma%cells(i,j,k) - 1.0_dkind)

				!f = (0.25_dkind * E_f_prod(i,j,k) - 0.25_dkind * v%pr(dim)%cells(i,j,k) * v_prod(dim,i,j,k)) * (gamma%cells(i,j,k) - 1.0_dkind)
				!g_inv =  f /  rho%cells(i,j,k) / v_s%cells(i,j,k) + 0.25_dkind * v_prod(dim,i,j,k) /  rho%cells(i,j,k)
				
				max_inv = 0.0_dkind
				min_inv = 0.0_dkind
			
				g_inv = - f /  rho%cells(i,j,k) / v_s%cells(i,j,k) - 0.25_dkind * rho_prod(i,j,k) * (v%pr(dim)%cells(i,j,k) - v_s%cells(i,j,k)) / rho%cells(i,j,k) + v_prod(dim,i,j,k) /  rho%cells(i,j,k)

				do dim1 = 1,this%domain%dimensions
					if ( dim1 /= dim ) then
						g_inv = g_inv - 0.5_dkind * this%time_step*v_s%cells(i,j,k) * (v_f(dim1,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - v_f(dim1,dim1,i,j,k))/mesh%cell_size(1)
						g_inv = g_inv - 0.5_dkind * this%time_step*v%pr(dim1)%cells(i,j,k)*(v_f(dim,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - v_f(dim,dim1,i,j,k))/mesh%cell_size(1)
						g_inv = g_inv - 0.5_dkind * this%time_step*v%pr(dim1)%cells(i,j,k)/ rho%cells(i,j,k) / v_s%cells(i,j,k) *(p_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - p_f(dim1,i,j,k))/mesh%cell_size(1) 					
					end if
				end do
				
				

				max_inv = max(r(1),R_half,r(2)) + g_inv
				min_inv = min(r(1),R_half,r(2)) + g_inv

				if ((min_inv <= r_new(1)).and.(r_new(1) <= max_inv))    r_corrected(1) = r_new(1)
				if (r_new(1) < min_inv)                                 r_corrected(1) = min_inv
				if (max_inv < r_new(1))                                 r_corrected(1) = max_inv
                
				if ((min_inv <= r_new(2)).and.(r_new(2) <= max_inv))    r_corrected(2) = r_new(2)
				if (r_new(2) < min_inv)                                 r_corrected(2) = min_inv
				if (max_inv < r_new(2))                                 r_corrected(2) = max_inv               

				g_inv =  f /  rho%cells(i,j,k) / v_s%cells(i,j,k) - 0.25_dkind * rho_prod(i,j,k) * ( v%pr(dim)%cells(i,j,k) + v_s%cells(i,j,k)) / rho%cells(i,j,k) + v_prod(dim,i,j,k) /  rho%cells(i,j,k)
				
				do dim1 = 1,this%domain%dimensions
					if ( dim1 /= dim ) then
						g_inv = g_inv + 0.5_dkind * this%time_step*v_s%cells(i,j,k) * (v_f(dim1,dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3)) - v_f(dim1,dim1,i,j,k))/mesh%cell_size(1)
						g_inv = g_inv - 0.5_dkind * this%time_step*v%pr(dim1)%cells(i,j,k)*(v_f(dim,dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3)) - v_f(dim,dim1,i,j,k))/mesh%cell_size(1)
						g_inv = g_inv + 0.5_dkind * this%time_step*v%pr(dim1)%cells(i,j,k)/ rho%cells(i,j,k) / v_s%cells(i,j,k) *(p_f(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3)) - p_f(dim1,i,j,k))/mesh%cell_size(1) 
					end if
				end do

				max_inv = 0.0_dkind
				min_inv = 0.0_dkind				
				
				!g_inv = - f /  rho%cells(i,j,k) / v_s%cells(i,j,k) + 0.25_dkind * v_prod(dim,i,j,k) /  rho%cells(i,j,k)
				
				max_inv = max(q(1),Q_half,q(2)) + g_inv
				min_inv = min(q(1),Q_half,q(2)) + g_inv
			
				if ((min_inv <= q_new(1)).and.(q_new(1) <= max_inv))    q_corrected(1) = q_new(1)
				if (q_new(1) < min_inv)                                 q_corrected(1) = min_inv
				if (max_inv < q_new(1))                                 q_corrected(1) = max_inv
                    
				if ((min_inv <= q_new(2)).and.(q_new(2) <= max_inv))	q_corrected(2) = q_new(2)
				if (q_new(2) < min_inv)									q_corrected(2) = min_inv
				if (max_inv < q_new(2))									q_corrected(2) = max_inv

				do dim2 = 1,this%domain%dimensions
					if(dim2 == dim) then
						g_inv =  - f 
						do dim1 = 1,this%domain%dimensions
							if ( dim1 /= dim ) then
								g_inv = g_inv + 0.5_dkind * this%time_step*v_s%cells(i,j,k)*  v_s%cells(i,j,k) * v%pr(dim1)%cells(i,j,k) * (rho_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - rho_f(dim1,i,j,k))/mesh%cell_size(1)
								g_inv = g_inv - 0.5_dkind * this%time_step*v%pr(dim1)%cells(i,j,k) * (p_f(dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - p_f(dim1,i,j,k))/mesh%cell_size(1) 
							end if
						end do
					else
						g_inv = 0.0_dkind
						do dim1 = 1,this%domain%dimensions
							if ( dim1 /= dim ) then
								g_inv = g_inv + 0.5_dkind * this%time_step * rho%cells(i,j,k) * v%pr(dim1)%cells(i,j,k) * (v_f(dim2,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - v_f(dim2,dim1,i,j,k))/mesh%cell_size(1)
							end if
						end do
						
						g_inv = g_inv + v%pr(dim2)%cells(i,j,k) * 0.25_dkind * rho_prod(i,j,k) + 0.5_dkind * this%time_step*(p_f(dim2,i+i_m(dim2,1),j+i_m(dim2,2),k+i_m(dim2,3)) - p_f(dim2,i,j,k))/mesh%cell_size(1)
						
						g_inv =  (v_prod(dim2,i,j,k) - g_inv) / rho%cells(i,j,k) 
					end if
    
					max_inv = max(v_inv(dim2,1),v_inv_half(dim2),v_inv(dim2,2)) + g_inv
					min_inv = min(v_inv(dim2,1),v_inv_half(dim2),v_inv(dim2,2)) + g_inv
						
					corr = 0.0_dkind
					if((max_inv - min_inv) /= 0.0 ) corr	= (alpha)/(max_inv - min_inv)
					max_inv = max_inv + corr
					min_inv = min_inv - corr						
					
					if ((min_inv <= v_inv_new(dim2,1)).and.(v_inv_new(dim2,1) <= max_inv))		v_inv_corrected(dim2,1) = v_inv_new(dim2,1)
					if (v_inv_new(dim2,1) < min_inv)											v_inv_corrected(dim2,1) = min_inv
					if (max_inv < v_inv_new(dim2,1))											v_inv_corrected(dim2,1) = max_inv
                
					if ((min_inv <= v_inv_new(dim2,2)).and.(v_inv_new(dim2,2) <= max_inv))		v_inv_corrected(dim2,2) = v_inv_new(dim2,2)
					if (v_inv_new(dim2,2) < min_inv)											v_inv_corrected(dim2,2) = min_inv
					if (max_inv < v_inv_new(dim2,2))											v_inv_corrected(dim2,2) = max_inv   											
    
				end do

				!max_inv = 0.0_dkind
				!min_inv = 0.0_dkind				
				!
				!g_inv =  f
				!
				!max_inv = max(v_inv(dim,1),v_inv_half(dim),v_inv(dim,2)) + g_inv
				!min_inv = min(v_inv(dim,1),v_inv_half(dim),v_inv(dim,2)) + g_inv
				!	
				!if ((min_inv <= v_inv_new(dim,1)).and.(v_inv_new(dim,1) <= max_inv))	v_inv_corrected(dim,1) = v_inv_new(dim,1)
				!if (v_inv_new(dim,1) < min_inv)											v_inv_corrected(dim,1) = min_inv
				!if (max_inv < v_inv_new(dim,1))											v_inv_corrected(dim,1) = max_inv
    !            
				!if ((min_inv <= v_inv_new(dim,2)).and.(v_inv_new(dim,2) <= max_inv))	v_inv_corrected(dim,2) = v_inv_new(dim,2)
				!if (v_inv_new(dim,2) < min_inv)											v_inv_corrected(dim,2) = min_inv
				!if (max_inv < v_inv_new(dim,2))											v_inv_corrected(dim,2) = max_inv   
				
				! ************* Lower edge *************
				! ************* Approximated velocity and speed of sound *******

				if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,1) ) then
				
					if	((.not.((abs(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	<	abs(	v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))))	.and.(	  v%pr(dim)%cells(i,j,k)	>		v_s%cells(i,j,k))))	&
					.and.(.not.((	 v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	<			-v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	.and.(abs(v%pr(dim)%cells(i,j,k))	< abs(	v_s%cells(i,j,k)))))) then 
					
						v_f_approx		= 0.5_dkind*(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ v%pr(dim)%cells(i,j,k))
						v_s_f_approx	= 0.5_dkind*(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ v_s%cells(i,j,k))

						characteristic_speed(1) = v_f_approx + v_s_f_approx
						characteristic_speed(2) = v_f_approx - v_s_f_approx
						characteristic_speed(3) = v_f_approx
				
						if (( characteristic_speed(1) >= 0.0_dkind )	.and.&
							( characteristic_speed(2) < 0.0_dkind )		.and.&
							( characteristic_speed(3) >= 0.0_dkind )) then				
				
							p_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k)	-	q_corrected(1)				)/ (G_half_lower + G_half)
							rho_f_new%cells(dim,i,j,k)			=	(rho_f_new%cells(dim,i,j,k)	+	p_f_new%cells(dim,i,j,k)	)/ (v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2)
					
							do dim1 = 1,this%domain%dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	(v_f_new%pr(dim1)%cells(dim,i,j,k)	+ G_half_lower	*	q_corrected(1))	/ (G_half_lower + G_half)
								end if
							end do

						end if		
		
						if (( characteristic_speed(1) >= 0.0_dkind ).and.&
							( characteristic_speed(2) < 0.0_dkind ).and.&
							( characteristic_speed(3) < 0.0_dkind )) then
					
							p_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k)	-	q_corrected(1)		)	/ (G_half_lower + G_half)
							rho_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k)	-	v_inv_corrected(dim,1))	/ (v_s%cells(i,j,k)**2)
					
							do dim1 = 1,this%domain%dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	(v_f_new%pr(dim1)%cells(dim,i,j,k)	+  G_half_lower	*	q_corrected(1))	/ (G_half_lower + G_half)
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corrected(dim1,1)
								end if
							end do
							continue
						end if
				

						if (( characteristic_speed(1) >= 0.0_dkind ).and.&
							( characteristic_speed(2) >= 0.0_dkind ).and.&
							( characteristic_speed(3) >= 0.0_dkind )) then
					
						end if
				
						if (( characteristic_speed(1) < 0.0_dkind ).and.&
							( characteristic_speed(2) < 0.0_dkind ).and.&
							( characteristic_speed(3) < 0.0_dkind )) then
					
							p_f_new%cells(dim,i,j,k)			= 0.5_dkind * (r_corrected(1) - q_corrected(1)) / G_half
					
							rho_f_new%cells(dim,i,j,k)			= (p_f_new%cells(dim,i,j,k) - v_inv_corrected(dim,1)) / (v_s%cells(i,j,k)**2)
					
							do dim1 = 1,this%domain%dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	0.5_dkind * (r_corrected(1) + q_corrected(1))
								else
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corrected(dim1,1)
								end if
							end do
						end if
						
						if (abs(characteristic_speed(3)) < 1.0E-10_dkind) then
							rho_f_new%cells(dim,i,j,k) = rho_f(dim,i,j,k)
							do dim1 = 1,this%domain%dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_f(dim1,dim,i,j,k)
								end if
							end do
						end if
						
					end if
				end if
					
				! ************* Higher edge *************
				! ************* Approximated velocity and speed of sound *******
				
				
				if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,2) ) then
				
					if ((	.not.((abs(	v%pr(dim)%cells(i,j,k))	<	abs(v_s%cells(i,j,k)))	.and.(		v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))		>		v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))))		&
					.and.(	.not.((		v%pr(dim)%cells(i,j,k)	<		-v_s%cells(i,j,k))	.and.(abs(	v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))	< abs(	v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))))))then
				
						v_f_approx		= 0.5_dkind*(v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+ v%pr(dim)%cells(i,j,k))
						v_s_f_approx	= 0.5_dkind*(v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			+ v_s%cells(i,j,k))
				
						characteristic_speed(1) = v_f_approx + v_s_f_approx
						characteristic_speed(2) = v_f_approx - v_s_f_approx
						characteristic_speed(3) = v_f_approx
				
						if (( characteristic_speed(1) >= 0.0_dkind )	.and.&
							( characteristic_speed(2) < 0.0_dkind )		.and.&
							( characteristic_speed(3) >= 0.0_dkind )) then				
				
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=		r_corrected(2)			!/ (G_half_higher + G_half)
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	-	v_inv_corrected(dim,2)	!/ (v_s%cells(i,j,k)**2)
					
							do dim1 = 1,this%domain%dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	(G_half_higher	*	r_corrected(2))	!/ (G_half_higher + G_half)
								else
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	v_inv_corrected(dim1,2)
									if ( characteristic_speed(3) == 0 ) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
									end if
								end if
							end do
						end if		
		
						if (( characteristic_speed(1) >= 0.0_dkind ).and.&
							( characteristic_speed(2) < 0.0_dkind ).and.&
							( characteristic_speed(3) < 0.0_dkind )) then
					
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	r_corrected(2)			!/ (G_half_higher + G_half)
						
							do dim1 = 1,this%domain%dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	G_half_higher	*	r_corrected(2)	!/ (G_half_higher + G_half)
								end if
							end do
						end if
				

						if (( characteristic_speed(1) >= 0.0_dkind ).and.&
							( characteristic_speed(2) >= 0.0_dkind ).and.&
							( characteristic_speed(3) >= 0.0_dkind )) then
					
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= 0.5_dkind * (r_corrected(2) - q_corrected(2)) / G_half
					
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = (p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(dim,2)) / (v_s%cells(i,j,k)**2)
					
							do dim1 = 1,this%domain%dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	0.5_dkind * (r_corrected(2) + q_corrected(2))
								else
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	v_inv_corrected(dim1,2)
									if ( characteristic_speed(3) == 0 ) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
									end if
								end if
							end do
						end if
				
						if (( characteristic_speed(1) < 0.0_dkind ).and.&
							( characteristic_speed(2) < 0.0_dkind ).and.&
							( characteristic_speed(3) < 0.0_dkind )) then
						end if
					end if
				end if
					
				!**************************** Sound points *****************************
				!**************************** Lower edge *******************************
				
				if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,1) ) then
				
					if (((abs(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	<	abs(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))) &
					.and.(	  v%pr(dim)%cells(i,j,k)	>	v_s%cells(i,j,k))) )then 
					
						if (characteristic_speed(3) < 0.0_dkind) then
							rho_f_new%cells(dim,i,j,k) = (p_f_new%cells(dim,i,j,k) - v_inv_corrected(dim,1)) / (v_s%cells(i,j,k)**2)
						end if
    
					end if
					if (((	 v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	 <	  -v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)))	&
					.and.	(abs(v%pr(dim)%cells(i,j,k)) < abs(v_s%cells(i,j,k)))))then 
				
						v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.5_dkind*(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))/v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v%pr(dim)%cells(i,j,k)/v_s%cells(i,j,k)) &
																*0.5_dkind*(v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3)) + v_s%cells(i,j,k))
    
						p_f_new%cells(dim,i,j,k)			= (v_f_new%pr(dim)%cells(dim,i,j,k) - q_corrected(1))/G_half
					
						if( characteristic_speed(3) >= 0.0_dkind ) then
							rho_f_new%cells(dim,i,j,k) = (rho_f_new%cells(dim,i,j,k) + p_f_new%cells(dim,i,j,k)) / (v_s%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))**2)
						else
							rho_f_new%cells(dim,i,j,k) = (p_f_new%cells(dim,i,j,k) - v_inv_corrected(dim,1)) / (v_s%cells(i,j,k)**2)
						end if
					end if	
				end if
				
				!**************************** Higher edge *******************************
				
				if ((I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,2) ) then
				
					if (((abs(v%pr(dim)%cells(i,j,k))	<	abs(v_s%cells(i,j,k))) &
					.and.(	  v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	>	v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))) )then 
					
						v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	0.5_dkind*(v%pr(dim)%cells(i,j,k)/v_s%cells(i,j,k) + v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))/v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) &
																								*0.5_dkind*(v_s%cells(i,j,k) + v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))
    
						p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= (r_corrected(2) - v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/G_half
				
						if (characteristic_speed(3) >= 0.0_dkind) then
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))		= (p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(dim,2)) / (v_s%cells(i,j,k)**2)
						end if
    
					end if
					if (((	 v%pr(dim)%cells(i,j,k)	 <	  -v_s%cells(i,j,k))	&
					.and.	(abs(v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))) < abs(v_s%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))))))then 
	   
						if( characteristic_speed(3) >= 0.0_dkind ) then
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = - v_inv_corrected(dim,2)
						end if
					end if	
	   
				end if
				
				
				!**************************** Boundary conditions *****************************
				
				if ((I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) == loop_bound(dim,1)) then
    
						v_f_approx		= v%pr(dim)%cells(i,j,k)
						v_s_f_approx	= v_s%cells(i,j,k)						
    
						characteristic_speed(1) = v_f_approx + v_s_f_approx
						characteristic_speed(2) = v_f_approx - v_s_f_approx
						characteristic_speed(3) = v_f_approx				
    
						if (( characteristic_speed(1) >= 0.0_dkind )	.and.&
							( characteristic_speed(2) < 0.0_dkind )		.and.&
							( characteristic_speed(3) >= 0.0_dkind )) then				
    
							v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.0_dkind 
							
							do dim1 = 1,this%domain%dimensions
								if( dim1 /= dim) then							
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v%pr(dim1)%cells(i,j,k)
								end if
							end do		
    
							p_f_new%cells(dim,i,j,k)			=	-q_corrected(1)	/ G_half
							
							rho_f_new%cells(dim,i,j,k)			=	rho%cells(i,j,k)
							
						end if
						

						if (( characteristic_speed(1) >= 0.0_dkind )	.and.&
							( characteristic_speed(2) < 0.0_dkind )		.and.&
							( characteristic_speed(3) < 0.0_dkind )) then	
							
							p_f_new%cells(dim,i,j,k)			=	-q_corrected(1)/ G_half
							
							rho_f_new%cells(dim,i,j,k)			=	(p_f_new%cells(dim,i,j,k)	-	v_inv_corrected(dim,1))	/ (v_s%cells(i,j,k)**2)
					
							v_f_new%pr(dim)%cells(dim,i,j,k)	=	0.0_dkind
							
							do dim1 = 1,this%domain%dimensions
								if( dim1 /= dim) then
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_inv_corrected(dim1,1)
								end if
							end do
						
						end if
					
						if (abs(characteristic_speed(3)) < 1.0E-10_dkind) then
							rho_f_new%cells(dim,i,j,k) = rho_f(dim,i,j,k)
							do dim1 = 1,this%domain%dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i,j,k)	=	v_f(dim1,dim,i,j,k)
								end if
							end do
						end if
						
				end if 				
    
				if ((I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) == loop_bound(dim,2)) then
				
						v_f_approx		= v%pr(dim)%cells(i,j,k)
						v_s_f_approx	= v_s%cells(i,j,k)						
    
						characteristic_speed(1) = v_f_approx + v_s_f_approx
						characteristic_speed(2) = v_f_approx - v_s_f_approx
						characteristic_speed(3) = v_f_approx				
				  
						if (( characteristic_speed(1) >= 0.0_dkind )	.and.&
							( characteristic_speed(2) < 0.0_dkind )		.and.&
							( characteristic_speed(3) >= 0.0_dkind )) then				
							
							v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dkind
							
							do dim1 = 1, this%domain%dimensions
								if (dim1 /= dim) then
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v_inv_corrected(dim1,2)
								end if
							end do
								
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= r_corrected(2)/G_half
							
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= (p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(dim,2)) / (v_s%cells(i,j,k)**2)
										
						end if 

						if (( characteristic_speed(1) >= 0.0_dkind )	.and.&
							( characteristic_speed(2) < 0.0_dkind )		.and.&
							( characteristic_speed(3) < 0.0_dkind )) then
      
							v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = 0.0_dkind
							
							do dim1 = 1, this%domain%dimensions
								if (dim1 /= dim) then
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v%pr(dim1)%cells(i,j,k)
								end if
							end do
						
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= rho%cells(i,j,k)
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= r_corrected(2)/G_half
      
						end if

						! ********************************************************************
						
						if (abs(characteristic_speed(3)) < 1.0E-10_dkind) then
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
							do dim1 = 1,this%domain%dimensions
								if ( dim /= dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
								end if
							end do
						end if	
				end if 					
				
				if (I_m(dim,1)*i == loop_bound(dim,2)) then
				
						v_f_approx		= v%pr(dim)%cells(i,j,k)
						v_s_f_approx	= v_s%cells(i,j,k)						
    
						characteristic_speed(1) = v_f_approx + v_s_f_approx
						characteristic_speed(2) = v_f_approx - v_s_f_approx
						characteristic_speed(3) = v_f_approx
						
						! ******************* Acoustic outlet ***********************************						
						
						if (( characteristic_speed(1) >= 0.0_dkind )	.and.&
							( characteristic_speed(2) < 0.0_dkind )) then		!.and.&
							!( characteristic_speed(3) >= 0.0_dkind )) then		
						
							v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	v%pr(dim)%cells(i,j,k)	!r_corrected(2) - p%cells(i,j,k)*G_half	

							!(g_inf * r_corrected(2) - G_half * g_inf * p_inf) / (G_half + g_inf) 
							!v%pr(dim)%cells(i,j,k)/v_s%cells(i,j,k) * ( 2.0_dkind * v_s%cells(i,j,k) - v_s_f(dim,i,j,k)) 
							!0.5_dkind * (r_corrected(2) + q_corrected(2)) !sqrt(sqrt(((p%cells(i,j,k)-p_inf)*(rho%cells(i,j,k)-rho_inf)/(rho%cells(i,j,k)*rho_inf))**2))
							
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	 p%cells(i,j,k) !(r_corrected(2) - v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/G_half
							
							!(r_corrected(2) + g_inf * p_inf)/( G_half + g_inf )
							!(r_corrected(2) - v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)))/G_half
							
							if ( characteristic_speed(3) >= 0.0_dkind ) then
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	rho%cells(i,j,k)!(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(dim,2)) / (v_s%cells(i,j,k)**2)	!rho%cells(i,j,k)	!(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(dim,2)) / (v_s%cells(i,j,k)**2)
							else
								rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	rho%cells(i,j,k)!(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - (p_inf) + c_inf**2 * rho_inf) / (c_inf**2)			!rho%cells(i,j,k)	!(p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - p_inf + c_inf**2 * rho_inf) / (c_inf**2)
							end if

							
							if ( v_f_new%pr(dim)%cells(dim,i-1,j,k) > 0.0_dkind) then
								continue
							end if
														
							!	E_f_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			=	E_f%cells(i,j,k)

						end if
						
						! ******************* Shock outlet ***********************************
	   
						if (( characteristic_speed(1) >= 0.0_dkind ).and.&
							( characteristic_speed(2) >= 0.0_dkind ).and.&
							( characteristic_speed(3) >= 0.0_dkind )) then
					
							p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= 0.5_dkind * (r_corrected(2) - q_corrected(2)) / G_half
					
							rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = (p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) - v_inv_corrected(dim,2)) / (v_s%cells(i,j,k)**2)
					
							do dim1 = 1,this%domain%dimensions
								if ( dim == dim1 ) then 
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	0.5_dkind * (r_corrected(2) + q_corrected(2))
								else
									v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	=	v_inv_corrected(dim1,2)
									if ( characteristic_speed(3) == 0 ) then
										v_f_new%pr(dim1)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = v_f(dim1,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
									end if
								end if
							end do
						end if
				end if

			end do	
			
		end do
		end do
		end do
        ! ************************************************  		

		loop_bound		= this%domain%short_cycle
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)		
		
			do dim = 1,this%domain%dimensions
	
				v_f_approx_lower		= 0.5_dkind*(v%pr(dim)%cells(i-I_m(dim,1),j-I_m(dim,2),k-I_m(dim,3))	+ v%pr(dim)%cells(i,j,k)) !- v_f(dim,dim,i,j,k)
				v_f_approx_higher		= 0.5_dkind*(v%pr(dim)%cells(i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	+ v%pr(dim)%cells(i,j,k)) !- v_f(dim,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))

				do spec = 1,this%chem%chem_ptr%species_number
					y_inv(spec,1)		= Y_f(spec,dim,i,j,k)									* rho_f(dim,i,j,k) !rho_f_new%cells(dim,i,j,k)											
					y_inv(spec,2)		= Y_f(spec,dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	* rho_f(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) !rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
					Y_inv_half(spec)	= Y%pr(spec)%cells(i,j,k)								* rho%cells(i,j,k)	
				end do	
				
				do spec = 1,this%chem%chem_ptr%species_number
					y_inv_new(spec,1)	= 2.0_dkind*Y_inv_half(spec) - y_inv(spec,2)
					y_inv_new(spec,2)	= 2.0_dkind*Y_inv_half(spec) - y_inv(spec,1)
				end do	
				
				do spec = 1,this%chem%chem_ptr%species_number
					y_inv_old(spec)		= Y_old(spec,i,j,k)	
				end do
				
				do spec = 1,this%chem%chem_ptr%species_number
					g_inv = 0.25_dkind * Y_prod(spec,i,j,k) 
					
					do dim1 = 1,this%domain%dimensions
						if ( dim1 /= dim ) then
							g_inv = g_inv - this%time_step * v%pr(dim1)%cells(i,j,k) * rho%cells(i,j,k) * (Y_f(spec,dim1,i+I_m(dim1,1),j+I_m(dim1,2),k+I_m(dim1,3)) - Y_f(spec,dim1,i,j,k))/mesh%cell_size(1)
						end if
					end do

					max_inv	= max(y_inv(spec,1),Y_inv_half(spec),y_inv(spec,2)) +  g_inv 
					min_inv	= min(y_inv(spec,1),Y_inv_half(spec),y_inv(spec,2)) +  g_inv 
					
					if ((min_inv <= y_inv_new(spec,1)).and.(y_inv_new(spec,1) <= max_inv))		y_inv_corrected(spec,1) = y_inv_new(spec,1)
					if (y_inv_new(spec,1) < min_inv)											y_inv_corrected(spec,1) = min_inv
					if (max_inv < y_inv_new(spec,1))											y_inv_corrected(spec,1) = max_inv
                
					if ((min_inv <= y_inv_new(spec,2)).and.(y_inv_new(spec,2) <= max_inv))		y_inv_corrected(spec,2) = y_inv_new(spec,2)
					if (y_inv_new(spec,2) < min_inv)											y_inv_corrected(spec,2) = min_inv
					if (max_inv < y_inv_new(spec,2))											y_inv_corrected(spec,2) = max_inv
				
				end do					
				
				spec_summ = 0.0_dkind
				do spec = 1,this%chem%chem_ptr%species_number
					if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,2) ) then
						!if ( v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) >= 0.0_dkind) then
						if (v_f_approx_higher >= 0.0_dkind) then
							Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = y_inv_corrected(spec,2) / rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
						end if
					end if
					
					if ( (I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) /= loop_bound(dim,1) ) then
						!if ( v_f_new%pr(dim)%cells(dim,i,j,k) < 0.0_dkind) then
						if (v_f_approx_lower < 0.0_dkind) then
							Y_f_new%pr(spec)%cells(dim,i,j,k) = y_inv_corrected(spec,1) / rho_f_new%cells(dim,i,j,k) 
						end if					
					end if

					if ((I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) == loop_bound(dim,1)) then
						if ( v_f_new%pr(dim)%cells(dim,i,j,k) < 0.0_dkind) then
							Y_f_new%pr(spec)%cells(dim,i,j,k) = y_inv_corrected(spec,1) / rho_f_new%cells(dim,i,j,k)
						else
							Y_f_new%pr(spec)%cells(dim,i,j,k) = Y%pr(spec)%cells(i,j,k)  
						end if
					end if
					
					if ((I_m(dim,1)*i + I_m(dim,2)*j + I_m(dim,3)*k) == loop_bound(dim,2)) then
						!if ( v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) > 0.0_dkind) then
						if ( v_f_approx_higher >= 0.0_dkind ) then
							Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = Y%pr(spec)%cells(i,j,k) !y_inv_corrected(spec,2) / rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))
						else
							Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3)) = Y%pr(spec)%cells(i,j,k)
						end if
					end if
					
					if (abs(v_f_approx_lower) < 1.0E-10_dkind) then
						Y_f_new%pr(spec)%cells(dim,i,j,k) = Y_f(spec,dim,i,j,k)
					end if
				
					spec_summ = spec_summ + max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dkind)
				end do

				do spec = 1,this%chem%chem_ptr%species_number
					Y_f_new%pr(spec)%cells(dim,i,j,k) = max(Y_f_new%pr(spec)%cells(dim,i,j,k), 0.0_dkind) / spec_summ
				end do

			end do	
			
		end do
		end do
		end do			
		
		! ******************* Eqn of state ***************
		call this%state_eq%apply_state_equation_flow_variables() 
		
        ! ************************************************  
		loop_bound		= this%domain%short_cycle
        ! *********** Conservative variables calculation ***************
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)
  
			
			rho%cells(i,j,k)		= 0.0_dkind
			do dim = 1,this%domain%dimensions
				mean_higher	= 0.5_dkind*(rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) *v_f(dim,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) + rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))  *v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))) 
				mean_lower	= 0.5_dkind*(rho_f(dim,i,j,k)   *v_f(dim,dim,i,j,k) +   rho_f_new%cells(dim,i,j,k)    *v_f_new%pr(dim)%cells(dim,i,j,k))
			
				rho%cells(i,j,k)	=	rho%cells(i,j,k)	- (	mean_higher - mean_lower )	/mesh%cell_size(1)
			end do
			rho%cells(i,j,k)		=	rho_old(i,j,k) + rho_prod(i,j,k) + this%time_step*rho%cells(i,j,k)
  
			spec_summ = 0.0_dkind
			do spec = 1,this%chem%chem_ptr%species_number
				y%pr(spec)%cells(i,j,k)		= 0.0_dkind
				do	dim = 1,this%domain%dimensions
					mean_higher	= 0.5_dkind*( rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) * y_f(spec,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) *v_f(dim,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))  &
											+ rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))  * y_f_new%pr(spec)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) *v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))
					
					mean_lower	= 0.5_dkind*( rho_f(dim,i,j,k)   *y_f(spec,dim,i,j,k)   *v_f(dim,dim,i,j,k)   + rho_f_new%cells(dim,i,j,k)    * y_f_new%pr(spec)%cells(dim,i,j,k)   *v_f_new%pr(dim)%cells(dim,i,j,k))
  
					y%pr(spec)%cells(i,j,k)	=  y%pr(spec)%cells(i,j,k)	-	(mean_higher - mean_lower )/mesh%cell_size(1)
				end do
				
				y%pr(spec)%cells(i,j,k)		=	rho_old(i,j,k) * y_old(spec,i,j,k) + y_prod(spec,i,j,k) + this%time_step * y%pr(spec)%cells(i,j,k)
				y%pr(spec)%cells(i,j,k)		=	y%pr(spec)%cells(i,j,k)	/ rho%cells(i,j,k)
				
				spec_summ = spec_summ + max(Y%pr(spec)%cells(i,j,k), 0.0_dkind)
			end do
			
			do spec = 1,this%chem%chem_ptr%species_number
				y%pr(spec)%cells(i,j,k) = max(y%pr(spec)%cells(i,j,k), 0.0_dkind) / spec_summ 
			end do


			do dim = 1,this%domain%dimensions
				v%pr(dim)%cells(i,j,k)		=	0.0_dkind 
				do dim1 = 1,this%domain%dimensions
					mean_higher	= 0.5_dkind*(		rho_f(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3)) *v_f(dim,dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3))*v_f(dim1,dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3))  &
												+	rho_f_new%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3))  *v_f_new%pr(dim)%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3)) *v_f_new%pr(dim1)%cells(dim1,i+i_m(dim1,1),j+i_m(dim1,2),k+i_m(dim1,3))  )
					mean_lower	= 0.5_dkind*(		rho_f(dim1,i,j,k) *v_f(dim,dim1,i,j,k)*v_f(dim1,dim1,i,j,k)  &
												+	rho_f_new%cells(dim1,i,j,k)  *v_f_new%pr(dim)%cells(dim1,i,j,k) *v_f_new%pr(dim1)%cells(dim1,i,j,k)  )
												
					v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(	mean_higher - mean_lower)	/mesh%cell_size(1)
				end do
  
				mean_higher	= 0.5_dkind*(p_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))	+	p_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)) )
				mean_lower	= 0.5_dkind*(p_f(dim,i,j,k) 									+	p_f_new%cells(dim,i,j,k) )
												
				v%pr(dim)%cells(i,j,k)	=  v%pr(dim)%cells(i,j,k)	-	(	mean_higher - mean_lower)	/mesh%cell_size(1)
				
				v%pr(dim)%cells(i,j,k)	= rho_old(i,j,k)*v_old(dim,i,j,k) + v_prod(dim,i,j,k)  +  this%time_step * v%pr(dim)%cells(i,j,k) 
				v%pr(dim)%cells(i,j,k)	= v%pr(dim)%cells(i,j,k) / rho%cells(i,j,k) 
			end do	
  
			e_f%cells(i,j,k)		=	0.0_dkind  
			do dim = 1,this%domain%dimensions
				mean_higher	= 0.5_dkind*(		(rho_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))*E_f_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))	+	p_f(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))*v_f(dim,dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))	&
											+	(rho_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))*E_f_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3))	+	p_f_new%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))*v_f_new%pr(dim)%cells(dim,i+i_m(dim,1),j+i_m(dim,2),k+i_m(dim,3)))
											
				mean_lower	= 0.5_dkind*(		(rho_f(dim,i,j,k)*E_f_f(dim,i,j,k)						+	p_f(dim,i,j,k))*v_f(dim,dim,i,j,k)	&
											+	(rho_f_new%cells(dim,i,j,k)*E_f_f_new%cells(dim,i,j,k)	+	p_f_new%cells(dim,i,j,k))*v_f_new%pr(dim)%cells(dim,i,j,k))	
			
				E_f%cells(i,j,k)			= 	E_f%cells(i,j,k)		-	(	mean_higher - mean_lower)	/mesh%cell_size(1)
			end do	
			E_f%cells(i,j,k) = rho_old(i,j,k) * E_f_old(i,j,k) +  E_f_prod(i,j,k) + this%time_step * e_f%cells(i,j,k)
			E_f%cells(i,j,k) = E_f%cells(i,j,k) /rho%cells(i,j,k)
  
        end do
		end do
		end do
        ! **************************************************************

		if (this%reactive_flag)	then
			call this%chem_kin_solver%solve_chemical_kinetics(this%time_step)
		
			do k = loop_bound(3,1),loop_bound(3,2)
			do j = loop_bound(2,1),loop_bound(2,2)
			do i = loop_bound(1,1),loop_bound(1,2)
				E_f%cells(i,j,k) = E_f%cells(i,j,k) + E_f_prod_chem%cells(i,j,k) 
				do spec = 1,this%chem%chem_ptr%species_number
					Y%pr(spec)%cells(i,j,k)	 = Y%pr(spec)%cells(i,j,k)	+ Y_prod_chem%pr(spec)%cells(i,j,k)
				end do
			end do
			end do
			end do			
		end if
		
		call this%state_eq%apply_state_equation() 		
		
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)			
			do dim = 1,this%domain%dimensions
				if (I_m(dim,1)*i == loop_bound(dim,2)) then
					v_f_new%pr(dim)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= v%pr(dim)%cells(i,j,k)	
					p_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= p%cells(i,j,k)
					rho_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= rho%cells(i,j,k)
					E_f_f_new%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))			= E_f%cells(i,j,k)
					do spec = 1,this%chem%chem_ptr%species_number
						Y_f_new%pr(spec)%cells(dim,i+I_m(dim,1),j+I_m(dim,2),k+I_m(dim,3))	= y%pr(spec)%cells(i,j,k)
					end do
				end if
			end do
		end do
		end do
		end do	

		do spec = 1,this%chem%chem_ptr%species_number
			Y_f(spec,:,:,:,:)	= Y_f_new%pr(spec)%cells(:,:,:,:)		
		end do
		
		do dim = 1,this%domain%dimensions
			v_f(dim,:,:,:,:)     = v_f_new%pr(dim)%cells(:,:,:,:)
		end do
		
		v_s_f	= v_s_f_new%cells(:,:,:,:)
        p_f     = p_f_new%cells(:,:,:,:)
        rho_f   = rho_f_new%cells(:,:,:,:)
        e_i_f   = e_i_f_new%cells(:,:,:,:)
		E_f_f	= E_f_f_new%cells(:,:,:,:)		
		
		continue
		
		end associate
		
	
	end subroutine

	subroutine apply_boundary_conditions_main(this)

		class(cabaret_solver)		,intent(inout)		:: this

		integer	,dimension(3,2)	:: loop_bound

		integer	:: sign, bound_number
		integer :: i,j,k,plus,dim,dim1,dim2,specie_number

		associate(  T				=> this%T%s_ptr					, &
					mol_mix_conc	=> this%mol_mix_conc%s_ptr		, &
					p				=> this%p%s_ptr					, &
					rho				=> this%rho%s_ptr				, &
					v				=> this%v%v_ptr					, &
					v_f_new			=> this%v_f_new%v_ptr			, &
					Y				=> this%Y%v_ptr					, &
					bc				=> this%boundary%bc_ptr			, &
					mesh			=> this%mesh%mesh_ptr)

		!$omp parallel default(none)  private(i,j,k,plus,dim,dim1,sign,loop_bound,bound_number) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(T,p,rho,v,Y,mesh,bc)
			loop_bound		= this%domain%short_cycle
		!$omp do collapse(3) schedule(static)

			do k = loop_bound(3,1),loop_bound(3,2)
			do j = loop_bound(2,1),loop_bound(2,2)
			do i = loop_bound(1,1),loop_bound(1,2)
				if(bc%bc_markers(i,j,k) == 0) then
					do dim = 1,this%domain%dimensions
						do plus = 1,2
							sign			= (-1)**plus
							bound_number	= bc%bc_markers(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))
							if( bound_number /= 0 ) then

								p%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= p%cells(i,j,k)
								rho%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	= rho%cells(i,j,k)
								T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= T%cells(i,j,k)
								mol_mix_conc%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))		= mol_mix_conc%cells(i,j,k)
								
								do dim1 = 1, this%domain%dimensions
									if(dim1 == dim) then
										v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = - v%pr(dim1)%cells(i,j,k)
									else
										v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i,j,k)
									end if
								end do

								do specie_number = 1, this%chem%chem_ptr%species_number
									Y%pr(specie_number)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3))	=	Y%pr(specie_number)%cells(i,j,k)
								end do

								select case(bc%boundary_types(bound_number)%type_name)
									case('wall')
										if(bc%boundary_types(bound_number)%conductive) T%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = bc%boundary_types(bound_number)%wall_temperature
										if(.not.bc%boundary_types(bound_number)%slip) then
											do dim1 = 1, this%domain%dimensions
												if(dim1 /= dim) then
													v%pr(dim1)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = v%pr(dim1)%cells(i-sign*I_m(dim,1),j-sign*I_m(dim,2),k-sign*I_m(dim,3)) - 6.0_dkind * v%pr(dim1)%cells(i,j,k)
												end if
											end do
										end if
									case ('outlet')
										v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) =  sign*sqrt(abs((p%cells(i,j,k) - bc%boundary_types(bound_number)%farfield_pressure)*(rho%cells(i,j,k) - bc%boundary_types(bound_number)%farfield_density)/bc%boundary_types(bound_number)%farfield_density/rho%cells(i,j,k)))
									case ('inlet')
										v%pr(dim)%cells(i+sign*I_m(dim,1),j+sign*I_m(dim,2),k+sign*I_m(dim,3)) = -sign*sqrt(abs((p%cells(i,j,k) - bc%boundary_types(bound_number)%farfield_pressure)*(rho%cells(i,j,k) - bc%boundary_types(bound_number)%farfield_density)/bc%boundary_types(bound_number)%farfield_density/rho%cells(i,j,k)))
								end select

							end if
						end do
					end do
				end if
			end do
			end do
			end do

		!$omp end do nowait
		!$omp end parallel

		end associate

	end subroutine
	
	subroutine calculate_time_step(this)
		class(cabaret_solver)	,intent(inout)	:: this
		
		real(dkind)	:: delta_t_interm, time_step, velocity_value

		integer	,dimension(3,2)	:: loop_bound
		integer	:: sign
		integer :: i,j,k,dim

		time_step			= 1.0e-06_dkind

		associate(  v				=> this%v%v_ptr		, &
					v_s				=> this%v_s%s_ptr		, &
					bc				=> this%boundary%bc_ptr	, &
					mesh			=> this%mesh%mesh_ptr)
		
		!$omp parallel default(none)  private(i,j,k,dim,delta_t_interm,velocity_value,loop_bound) , &
		!$omp& firstprivate(this)	,&
		!$omp& shared(v,v_s,mesh,bc,time_step)
			loop_bound		= this%domain%short_cycle
		!$omp do collapse(3) schedule(static) reduction(min:time_step)
					
		do k = loop_bound(3,1),loop_bound(3,2)
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1),loop_bound(1,2)
			if(this%boundary%bc_ptr%bc_markers(i,j,k) == 0) then
				velocity_value		= 0.0_dkind
				do dim = 1,this%domain%dimensions
					velocity_value = velocity_value + v%pr(dim)%cells(i,j,k)*v%pr(dim)%cells(i,j,k)
				end do
				delta_t_interm = minval(mesh%cell_size,mesh%cell_size > 0.0_dkind) / (sqrt(velocity_value) + v_s%cells(i,j,k))
				if (delta_t_interm < time_step) then
					time_step = delta_t_interm
				end if
			end if
		end do
		end do
		end do
	
		!$omp end do nowait
		!$omp end parallel
		
		!if (this%courant_fraction * time_step < 5e-09_dkind) then
		!	this%time_step = 0.1_dkind * this%courant_fraction * time_step
		!else
			this%time_step = this%courant_fraction * time_step
		!end if
		!print *, this%courant_fraction * time_step
			
		end associate
			
	end subroutine
	
	
	pure function get_time_step(this)
		real(dkind)						:: get_time_step
		class(cabaret_solver)	,intent(in)		:: this

		get_time_step = this%time_step
	end function

	pure function get_time(this)
		real(dkind)						:: get_time
		class(cabaret_solver)	,intent(in)		:: this

		get_time = this%time
	end function

	
	subroutine check_symmetry(this)
		class(cabaret_solver)	,intent(in)		:: this
	
		real(dkind)	:: max_error, error
		integer		:: max_error_i,max_error_j,max_error_k
		
		integer	,dimension(3,2)	:: loop_bound
		integer :: i,j,k,dim
	
		associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					e_i			=> this%e_i%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
					gamma		=> this%gamma%s_ptr		, &
					
					v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &
					
					v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					E_f_f		=> this%E_f_f	, &
					e_i_f		=> this%e_i_f	, &
					p_f			=> this%p_f		, &
					v_s_f		=> this%v_s_f	, & 
					Y_f			=> this%Y_f		, &

					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					E_f_old		=> this%E_f_old	,&
					Y_old		=> this%Y_old	,&
					p_old		=> this%p_old	,&
					v_s_old		=> this%v_s_old	,&
					gamma_old	=> this%gamma_old	, &
					
					p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					v_s_f_new	=> this%v_s_f_new%s_ptr		, &
					E_f_f_new	=> this%E_f_f_new%s_ptr		, &
					e_i_f_new	=> this%e_i_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr		, &
					
					v_prod_visc		=> this%v_prod_visc%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					E_f_prod_chem 	=> this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat	=> this%E_f_prod_heat%s_ptr	, &
					E_f_prod_visc	=> this%E_f_prod_visc%s_ptr	, &
					E_f_prod		=> this%E_f_prod			, &
					Y_prod			=> this%Y_prod				, &
					rho_prod		=> this%rho_prod			, &

					mesh		=> this%mesh%mesh_ptr)
					
				
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
					
		loop_bound		= this%domain%short_cycle

		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (rho%cells(i,j,1) /= rho%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1)) then
				error = abs(rho%cells(i,j,1)  - rho%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'density error', max_error , 	max_error_i,max_error_j,  loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (p%cells(i,j,1) /= p%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1)) then
				error = abs(p%cells(i,j,1)  - p%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'pressure error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (E_f%cells(i,j,1) /= E_f%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1)) then
				error = abs(E_f%cells(i,j,1)  - E_f%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'energy error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j	

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(v%pr(1)%cells(i,j,1)) /= abs(v%pr(1)%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v%pr(1)%cells(i,j,1))  - abs(v%pr(1)%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'velocity_x error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j	
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(v%pr(2)%cells(i,j,1)) /= abs(v%pr(2)%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v%pr(2)%cells(i,j,1))  - abs(v%pr(2)%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'velocity_y error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j				

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(Y%pr(1)%cells(i,j,1)) /= abs(Y%pr(1)%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(Y%pr(1)%cells(i,j,1))  - abs(Y%pr(1)%cells(loop_bound(1,2) + 1 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'conc_1 error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(v_f_new%pr(1)%cells(1,i,j,1)) /= abs(v_f_new%pr(1)%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v_f_new%pr(1)%cells(1,i,j,1))  - abs(v_f_new%pr(1)%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow velocity_x error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1)+1,loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(v_f_new%pr(1)%cells(2,i,j,1)) /= abs(v_f_new%pr(1)%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2- j , 1))) then
				error = abs(abs(v_f_new%pr(1)%cells(2,i,j,1))  - abs(v_f_new%pr(1)%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow velocity_x error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(v_f_new%pr(2)%cells(1,i,j,1)) /= abs(v_f_new%pr(2)%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v_f_new%pr(2)%cells(1,i,j,1))  - abs(v_f_new%pr(2)%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow velocity_y error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j			
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1)+1,loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(v_f_new%pr(2)%cells(2,i,j,1)) /= abs(v_f_new%pr(2)%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2 - j , 1))) then
				error = abs(abs(v_f_new%pr(2)%cells(2,i,j,1))  - abs(v_f_new%pr(2)%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow velocity_y error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(rho_f_new%cells(1,i,j,1)) /= abs(rho_f_new%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(rho_f_new%cells(1,i,j,1))  - abs(rho_f_new%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow density error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1)+1,loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(rho_f_new%cells(2,i,j,1)) /= abs(rho_f_new%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2 - j , 1))) then
				error = abs(abs(rho_f_new%cells(2,i,j,1))  - abs(rho_f_new%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow density error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j			
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(p_f_new%cells(1,i,j,1)) /= abs(p_f_new%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(p_f_new%cells(1,i,j,1))  - abs(p_f_new%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow pressure error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j			
		

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1)+1,loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(p_f_new%cells(2,i,j,1)) /= abs(p_f_new%cells(2,loop_bound(1,2) +1 - i , loop_bound(2,2) + 2 - j , 1))) then
				error = abs(abs(p_f_new%cells(2,i,j,1))  - abs(p_f_new%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow pressure error' , max_error , 		max_error_i,max_error_j, 	loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0		
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(Y_f_new%pr(1)%cells(1,i,j,1)) /= abs(Y_f_new%pr(1)%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(Y_f_new%pr(1)%cells(1,i,j,1))  - abs(Y_f_new%pr(1)%cells(1,loop_bound(1,2) + 2 - i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow conc_1 error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1)+1,loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(Y_f_new%pr(1)%cells(2,i,j,1)) /= abs(Y_f_new%pr(1)%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2- j , 1))) then
				error = abs(abs(Y_f_new%pr(1)%cells(2,i,j,1))  - abs(Y_f_new%pr(1)%cells(2,loop_bound(1,2) + 1 - i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow conc_1 error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j			

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
					
		loop_bound		= this%domain%short_cycle

		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (rho%cells(i,j,1) /= rho%cells(j,i,1)) then
				error = abs(rho%cells(i,j,1)  - rho%cells(j,i,1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y density error', max_error , 	max_error_i,max_error_j

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (p%cells(i,j,1) /= p%cells(j,i,1)) then
				error = abs(p%cells(i,j,1)  - p%cells(j,i,1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y pressure error' , max_error , 		max_error_i,max_error_j	
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (E_f%cells(i,j,1) /= E_f%cells(j,i,1)) then
				error = abs(E_f%cells(i,j,1)  - E_f%cells(j,i,1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y energy error' , max_error , 		max_error_i,max_error_j

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(v%pr(1)%cells(i,j,1)) /= abs(v%pr(2)%cells(j,i,1))) then
				error = abs(abs(v%pr(1)%cells(i,j,1))  - abs(v%pr(2)%cells(j,i,1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y velocity error' , max_error , 		max_error_i,max_error_j

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)/2

			if (abs(Y%pr(1)%cells(i,j,1)) /= abs(Y%pr(1)%cells(j,i,1))) then
				error = abs(abs(Y%pr(1)%cells(i,j,1))  - abs(Y%pr(1)%cells(j,i,1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y conc_1 error' , max_error , 		max_error_i,max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(v_f_new%pr(1)%cells(1,i,j,1)) /= abs(v_f_new%pr(2)%cells(2,j,i,1))) then
				error = abs(abs(v_f_new%pr(1)%cells(1,i,j,1))  - abs(v_f_new%pr(2)%cells(2,j,i,1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y flow velocity_x error' , max_error , max_error_i,max_error_j	

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(v_f_new%pr(2)%cells(1,i,j,1)) /= abs(v_f_new%pr(1)%cells(2,j,i,1))) then
				error = abs(abs(v_f_new%pr(2)%cells(1,i,j,1))  - abs(v_f_new%pr(1)%cells(2,j,i,1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y flow velocity_y error' , max_error , 		max_error_i,max_error_j	

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(rho_f_new%cells(1,i,j,1)) /= abs(rho_f_new%cells(2,j,i,1))) then
				error = abs(abs(rho_f_new%cells(1,i,j,1))  - abs(rho_f_new%cells(2,j,i,1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y flow density error' , max_error , 		max_error_i,max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(p_f_new%cells(1,i,j,1)) /= abs(p_f_new%cells(2,j,i,1))) then
				error = abs(abs(p_f_new%cells(1,i,j,1))  - abs(p_f_new%cells(2,j,i,1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y flow pressure error' , max_error , 		max_error_i,max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0		
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)/2+1

			if (abs(Y_f_new%pr(1)%cells(1,i,j,1)) /= abs(Y_f_new%pr(1)%cells(2,j,i,1))) then
				error = abs(abs(Y_f_new%pr(1)%cells(1,i,j,1))  - abs(Y_f_new%pr(1)%cells(2,j,i,1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x_y flow conc_1 error' , max_error , 		max_error_i,max_error_j
		
		end associate		
		
		
	end subroutine
	
	subroutine check_symmetry_2(this)
		class(cabaret_solver)	,intent(in)		:: this
	
		real(dkind)	:: max_error, error
		integer		:: max_error_i,max_error_j,max_error_k
		
		integer	,dimension(3,2)	:: loop_bound
		integer :: i,j,k,dim
	
		associate(	rho			=> this%rho%s_ptr		, &
					p			=> this%p%s_ptr			, &
					E_f			=> this%E_f%s_ptr		, &
					e_i			=> this%e_i%s_ptr		, &
					v_s			=> this%v_s%s_ptr		, &
					gamma		=> this%gamma%s_ptr		, &
					
					v			=> this%v%v_ptr	, &
					Y			=> this%Y%v_ptr	, &
					
					v_f			=> this%v_f		, &
					rho_f		=> this%rho_f	, &
					E_f_f		=> this%E_f_f	, &
					e_i_f		=> this%e_i_f	, &
					p_f			=> this%p_f		, &
					v_s_f		=> this%v_s_f	, & 
					Y_f			=> this%Y_f		, &

					rho_old		=> this%rho_old	, &
					v_old		=> this%v_old	, &
					E_f_old		=> this%E_f_old	,&
					Y_old		=> this%Y_old	,&
					p_old		=> this%p_old	,&
					v_s_old		=> this%v_s_old	,&
					gamma_old	=> this%gamma_old	, &
					
					p_f_new		=> this%p_f_new%s_ptr		, &	
					rho_f_new	=> this%rho_f_new%s_ptr		, &
					v_s_f_new	=> this%v_s_f_new%s_ptr		, &
					E_f_f_new	=> this%E_f_f_new%s_ptr		, &
					e_i_f_new	=> this%e_i_f_new%s_ptr		, &
					Y_f_new		=> this%Y_f_new%v_ptr		, &
					v_f_new		=> this%v_f_new%v_ptr		, &
					
					v_prod_visc		=> this%v_prod_visc%v_ptr	, &
					Y_prod_chem		=> this%Y_prod_chem%v_ptr	, &
					Y_prod_diff		=> this%Y_prod_diff%v_ptr	, &
					E_f_prod_chem 	=> this%E_f_prod_chem%s_ptr	, &
					E_f_prod_heat	=> this%E_f_prod_heat%s_ptr	, &
					E_f_prod_visc	=> this%E_f_prod_visc%s_ptr	, &
					E_f_prod		=> this%E_f_prod			, &
					Y_prod			=> this%Y_prod				, &
					rho_prod		=> this%rho_prod			, &

					mesh		=> this%mesh%mesh_ptr)
					
				
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
					
		loop_bound		= this%domain%short_cycle

		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)

			if (rho%cells(i,j,1) /= rho%cells(i , loop_bound(2,2) + 1 - j , 1)) then
				error = abs(rho%cells(i,j,1)  - rho%cells(i , loop_bound(2,2) + 1 - j , 1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'density error', max_error , 	max_error_i,max_error_j,  loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)

			if (p%cells(i,j,1) /= p%cells(i , loop_bound(2,2) + 1 - j , 1)) then
				error = abs(p%cells(i,j,1)  - p%cells(i , loop_bound(2,2) + 1 - j , 1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'pressure error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)

			if (E_f%cells(i,j,1) /= E_f%cells(i , loop_bound(2,2) + 1 - j , 1)) then
				error = abs(E_f%cells(i,j,1)  - E_f%cells(i , loop_bound(2,2) + 1 - j , 1))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'energy error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j	

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(v%pr(1)%cells(i,j,1)) /= abs(v%pr(1)%cells(i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v%pr(1)%cells(i,j,1))  - abs(v%pr(1)%cells(i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'velocity_x error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j	
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(v%pr(2)%cells(i,j,1)) /= abs(v%pr(2)%cells(i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v%pr(2)%cells(i,j,1))  - abs(v%pr(2)%cells(i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'velocity_y error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j				

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(Y%pr(1)%cells(i,j,1)) /= abs(Y%pr(1)%cells(i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(Y%pr(1)%cells(i,j,1))  - abs(Y%pr(1)%cells(i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'conc_1 error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)

			if (abs(v_f_new%pr(1)%cells(1,i,j,1)) /= abs(v_f_new%pr(1)%cells(1,i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v_f_new%pr(1)%cells(1,i,j,1))  - abs(v_f_new%pr(1)%cells(1,i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow velocity_x error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1)+1,loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(v_f_new%pr(1)%cells(2,i,j,1)) /= abs(v_f_new%pr(1)%cells(2,i , loop_bound(2,2) + 2- j , 1))) then
				error = abs(abs(v_f_new%pr(1)%cells(2,i,j,1))  - abs(v_f_new%pr(1)%cells(2, i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow velocity_x error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)
		do i = loop_bound(1,1)+1,loop_bound(1,2)

			if (abs(v_f_new%pr(2)%cells(1,i,j,1)) /= abs(v_f_new%pr(2)%cells(1,i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(v_f_new%pr(2)%cells(1,i,j,1))  - abs(v_f_new%pr(2)%cells(1,i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow velocity_y error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j			
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1)+1,loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(v_f_new%pr(2)%cells(2,i,j,1)) /= abs(v_f_new%pr(2)%cells(2,i , loop_bound(2,2) + 2 - j , 1))) then
				error = abs(abs(v_f_new%pr(2)%cells(2,i,j,1))  - abs(v_f_new%pr(2)%cells(2,i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow velocity_y error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)

			if (abs(rho_f_new%cells(1,i,j,1)) /= abs(rho_f_new%cells(1, i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(rho_f_new%cells(1,i,j,1))  - abs(rho_f_new%cells(1, i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow density error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(rho_f_new%cells(2,i,j,1)) /= abs(rho_f_new%cells(2,i , loop_bound(2,2) + 2 - j , 1))) then
				error = abs(abs(rho_f_new%cells(2,i,j,1))  - abs(rho_f_new%cells(2,i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow density error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j			
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)

			if (abs(p_f_new%cells(1,i,j,1)) /= abs(p_f_new%cells(1, i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(p_f_new%cells(1,i,j,1))  - abs(p_f_new%cells(1, i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow pressure error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j			
		

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(p_f_new%cells(2,i,j,1)) /= abs(p_f_new%cells(2, i , loop_bound(2,2) + 2 - j , 1))) then
				error = abs(abs(p_f_new%cells(2,i,j,1))  - abs(p_f_new%cells(2, i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow pressure error' , max_error , 		max_error_i,max_error_j, 	loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j

		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0		
		
		do j = loop_bound(2,1),loop_bound(2,2)/2
		do i = loop_bound(1,1)+1,loop_bound(1,2)

			if (abs(Y_f_new%pr(1)%cells(1,i,j,1)) /= abs(Y_f_new%pr(1)%cells(1, i , loop_bound(2,2) + 1 - j , 1))) then
				error = abs(abs(Y_f_new%pr(1)%cells(1,i,j,1))  - abs(Y_f_new%pr(1)%cells(1, i , loop_bound(2,2) + 1 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'x flow conc_1 error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 2 - max_error_i, loop_bound(2,2) + 1 - max_error_j		
		
		max_error = 0.0_dkind
		max_error_i = 0
		max_error_j	= 0
		
		do j = loop_bound(2,1),loop_bound(2,2)/2+1
		do i = loop_bound(1,1),loop_bound(1,2)

			if (abs(Y_f_new%pr(1)%cells(2,i,j,1)) /= abs(Y_f_new%pr(1)%cells(2, i , loop_bound(2,2) + 2- j , 1))) then
				error = abs(abs(Y_f_new%pr(1)%cells(2,i,j,1))  - abs(Y_f_new%pr(1)%cells(2, i , loop_bound(2,2) + 2 - j , 1)))
				if ( error > max_error ) then
					max_error = error
					max_error_i = i
					max_error_j = j
				end if
			end if
		end do
		end do					
		
		write(*,'(A20,E14.7,4I5)') 'y flow conc_1 error' , max_error , 		max_error_i,max_error_j, 	 loop_bound(1,2) + 1 - max_error_i, loop_bound(2,2) + 2 - max_error_j			

			
		end associate		
		
		
	end subroutine	
end module
