! ## File: dynamics.f90
! ## - main program: running the SIS dynamics, based on SIS II algorithm.
! ## See README.md for more information and use
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article "Simulation of Markovian epidemic models on large networks"
! Copyright (C) 2016 Wesley Cota, Silvio C. Ferreira
! 
! Please cite the above cited paper as reference to our code.
! 
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------
! Author    : Wesley Cota
! Email     : wesley.cota@ufv.br
! Date      : July 2016
!-----------------------------------------------------------------------------
! See README.md for more details
! This code is available at <https://github.com/wcota/dynSIS>

program dynSIS
use mod_netdata
use mod_random
use mod_read_tools
implicit none
    
    ! File read/write arguments and parameter
    character*1024                :: f_output, f_temp
    integer, parameter            :: und_output = 10
    
    
    ! Dynamics samples variables
    integer                       :: dynp_i, dynp_sam
    real*8                        :: dynp_lb
    integer                       :: dynp_tmax
    real*8                        :: dynp_pINI ! fraction of network first infected (random)
    
    ! SIS-II - Dynamics Variables
    real*8                        :: dyn_p
    real*8                        :: dyn_t, dyn_dt
    integer, allocatable          :: dyn_ocp(:), dyn_sig(:)
    integer                       :: dyn_voc, dyn_sk ! # of infected vertex and # of infected edges
    
    ! SIS-II - Network structure variables
    integer                       :: net_kmax
    
    ! Output variables and average measures
    real*8, allocatable           :: avg_rho(:), avg_t(:)
    integer, allocatable          :: avg_sam(:)
    integer                       :: dyn_dt_pos, dyn_dt_pos_max
    
    call print_info('###############################################################################')
    call print_info('#### Simulation of Markovian epidemic models on networks: SIS-II algorithm ####')
    call print_info('##============ Copyright(C) 2016 Wesley Cota, Silvio C. Ferreira ============##')
    call print_info('##======= This code is available at <https://github.com/wcota/dynSIS> =======##')
    call print_info('##======= Please cite the above cited paper as reference to our code. =======##')
    call print_info('##=== This code is under GNU General Public License. Please see README.md ===##')
    call print_info('###############################################################################')
    
    ! initial value of input counter to read command arguments, used by mod_read_tools
    inp_pos = 1 
    
    ! Files arguments read
    call read_arg(f_input)
    call read_arg(f_output)
    
    ! Read network to memory
    call readEdges()
    
    ! To be used in the SIS-II algorithm.
    net_kmax = maxval(net_k)
    
    ! Initate the random generator
    call print_info('')
    call print_progress('Generating random seed')
    call random_ini()
    call print_done()
    
    ! We are ready! All Network data is here, now we need to read the dynamical parameters.
    call print_info('')
    call print_info('Now we need to read the dynamical parameters.')
    call read_dyn_parameters()
    
    ! Let's run the dynamics. But first, we allocate the average matrices
    allocate(avg_rho(dynp_tmax), avg_t(dynp_tmax), avg_sam(dynp_tmax))
    avg_rho = 0d0
    avg_t = 0d0
    avg_sam = 0
    dyn_dt_pos_max = 1
    call print_info('')
    call print_info('Running dynamics...')
    do dynp_i=1,dynp_sam
        write(f_temp,*) dynp_i
        call print_progress('Sample # '//trim(adjustl(f_temp)))
    
        ! Run dynamics
        call read_dyn_run()
        
        ! Open file and write info
        open(und_output,file=f_output)
        write(und_output,'(a)') "#@ Network file: "//trim(adjustl(f_input))
        write(und_output,'(a,i7)') "#@ Number of nodes: ", net_N
        write(und_output,'(a,i7)') "#@ Number of edges: ", net_skk
        write(und_output,'(a,i7)') "#! Samples: ", dynp_i
        write(und_output,'(a,f11.5)') "#! Infection rate lambda: ", dynp_lb
        write(und_output,'(a,i7)') "#! Maximum time steps: ", dynp_tmax
        write(und_output,'(a,f11.5)') "#! Fraction of infected vertices (initial condition): ", dynp_pINI
        do dyn_dt_pos = 1, dyn_dt_pos_max
            write(und_output,*) 1d0*avg_t(dyn_dt_pos)/avg_sam(dyn_dt_pos) , 1d0*avg_rho(dyn_dt_pos)/dynp_i
        enddo
        ! Close file
        close(und_output)
        call print_done()
    enddo
    
    call print_info('')
    call print_info('Everything ok!')
    call print_info('Input file: '//trim(adjustl(f_input)))
    call print_info('Output file: '//trim(adjustl(f_output)))
    
contains

    subroutine read_dyn_parameters()
        call read_i(dynp_sam,"How much dynamics samples?")
        call read_f(dynp_lb,"Value of infection rate lambda (mu is defined as equal to 1)")
        call read_i(dynp_tmax,"Maximum time steps (it stops if the absorbing state is reached)")
        call read_f(dynp_pINI,"Fraction of infected vertices on the network as initial condition (is random to each sample)")
        
        allocate(dyn_ocp(net_N),dyn_sig(net_N))
    end subroutine

    subroutine read_dyn_run()
    
        integer :: pos_ocp
        integer :: ver, pos_nei
        real*8  :: rnd
    
        call random_initial_condition()
        
        dyn_t = 0d0
        dyn_dt_pos = 1
        
        dyn_time_loop : do while (dyn_t <= dynp_tmax)
            
            dyn_p = 1d0*dyn_voc/ (dyn_voc + 1d0*dynp_lb * dyn_sk)
            rnd = random_d()
            
            if (rnd < dyn_p) then
                pos_ocp = random_int(1,dyn_voc)
                ver = dyn_ocp(pos_ocp)
                
                ! healed
                dyn_sig(ver) = 0
                dyn_sk = dyn_sk - net_k(ver)
                dyn_ocp(pos_ocp) = dyn_ocp(dyn_voc) ! change positions
                dyn_voc = dyn_voc - 1
            else
                select_neig : do
                    pos_ocp = random_int(1,dyn_voc)
                    ver = dyn_ocp(pos_ocp)
                    if (random_d() < 1d0*net_k(ver)/(1d0*net_kmax)) exit select_neig
                enddo select_neig
                
                pos_nei = random_int(net_ini(ver) , net_ini(ver) + net_k(ver) - 1)
                ver = net_con(pos_nei)
                
                if (dyn_sig(ver) == 0) then
                    dyn_sig(ver) = 1
                    dyn_sk = dyn_sk + net_k(ver)
                    dyn_voc = dyn_voc + 1
                    dyn_ocp(dyn_voc) = ver
                endif
            endif
            
            if (dyn_voc == 0) exit dyn_time_loop
            
            dyn_dt = 1d0/(dyn_voc + 1d0*dynp_lb * dyn_sk)
            
            dyn_t = dyn_t + dyn_dt
            
            do while (dyn_t >= dyn_dt_pos)
                dyn_dt_pos_max = max(dyn_dt_pos,dyn_dt_pos_max)
                avg_rho(dyn_dt_pos) = avg_rho(dyn_dt_pos) + 1d0*dyn_voc/net_N
                avg_t(dyn_dt_pos) = avg_t(dyn_dt_pos) + dyn_t
                avg_sam(dyn_dt_pos) = avg_sam(dyn_dt_pos) + 1
                dyn_dt_pos = dyn_dt_pos + 1
            enddo
            
        enddo dyn_time_loop
        
    end subroutine
    
    subroutine random_initial_condition()
        integer :: ver, vti
        
        dyn_sig = 0
        dyn_ocp = 0
        dyn_voc = 0
        dyn_sk = 0
        do vti = 1, int(net_N*dynp_pINI)
            vti_ver : do
                ver = random_int(1,net_N)
                if (dyn_sig(ver) == 0) then
                    dyn_voc = dyn_voc + 1
                    dyn_ocp(dyn_voc) = ver
                    dyn_sig(ver) = 1
                    dyn_sk = dyn_sk + net_k(ver)
                    exit vti_ver
                endif
            enddo vti_ver
        enddo
    end subroutine
    
end program
