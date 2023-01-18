!=======================================================================
! This is part of the 2DECOMP&FFT library
!
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil)
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the FFTW (version 3.x) implementation of the FFT library

module decomp_2d_fft

   use decomp_2d_constants
   use decomp_2d  ! 2D decomposition module
   use iso_c_binding

   implicit none

   include "fftw3.f"

   private        ! Make everything private unless declared public

   ! engine-specific global variables
   integer, save :: plan_type = FFTW_MEASURE

   ! FFTW plans
   ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
   ! For c2c transforms:
   !     use plan(-1,j) for  forward transform;
   !     use plan( 1,j) for backward transform;
   ! For r2c/c2r transforms:
   !     use plan(0,j) for r2c transforms;
   !     use plan(2,j) for c2r transforms;
   type(C_PTR), save :: plan(-1:2, 3)

   ! This is defined in fftw3.f03 but not in fftw3.f
   interface
      subroutine fftw_cleanup() bind(C, name='fftw_cleanup')
         import
      end subroutine fftw_cleanup
   end interface

   integer, parameter, public :: D2D_FFT_BACKEND = D2D_FFT_BACKEND_FFTW3

   ! common code used for all engines, including global variables,
   ! generic interface definitions and several subroutines
#include "fft_common.f90"

   ! Return a FFTW3 plan for multiple 1D c2c FFTs in X direction
   subroutine c2c_1m_x_plan(plan1, decomp, isign)

      implicit none

      type(C_PTR), intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: isign

      complex(mytype), allocatable, dimension(:, :, :) :: a1

      allocate (a1(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)))

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft(plan1, 1, decomp%xsz(1), &
                               decomp%xsz(2)*decomp%xsz(3), a1, decomp%xsz(1), 1, &
                               decomp%xsz(1), a1, decomp%xsz(1), 1, decomp%xsz(1), &
                               isign, plan_type)
#else
      call sfftw_plan_many_dft(plan1, 1, decomp%xsz(1), &
                               decomp%xsz(2)*decomp%xsz(3), a1, decomp%xsz(1), 1, &
                               decomp%xsz(1), a1, decomp%xsz(1), 1, decomp%xsz(1), &
                               isign, plan_type)
#endif

      deallocate (a1)

      return
   end subroutine c2c_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D c2c FFTs in Y direction
   subroutine c2c_1m_y_plan(plan1, decomp, isign)

      implicit none

      type(C_PTR), intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: isign

      complex(mytype), allocatable, dimension(:, :) :: a1

      ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
      ! done one Z-plane at a time. So plan for 2D data sets here.

      allocate (a1(decomp%ysz(1), decomp%ysz(2)))

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft(plan1, 1, decomp%ysz(2), decomp%ysz(1), &
                               a1, decomp%ysz(2), decomp%ysz(1), 1, a1, decomp%ysz(2), &
                               decomp%ysz(1), 1, isign, plan_type)
#else
      call sfftw_plan_many_dft(plan1, 1, decomp%ysz(2), decomp%ysz(1), &
                               a1, decomp%ysz(2), decomp%ysz(1), 1, a1, decomp%ysz(2), &
                               decomp%ysz(1), 1, isign, plan_type)
#endif

      deallocate (a1)

      return
   end subroutine c2c_1m_y_plan

   ! Return a FFTW3 plan for multiple 1D c2c FFTs in Z direction
   subroutine c2c_1m_z_plan(plan1, decomp, isign)

      implicit none

      type(C_PTR), intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: isign

      complex(mytype), allocatable, dimension(:, :, :) :: a1

      allocate (a1(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)))

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft(plan1, 1, decomp%zsz(3), &
                               decomp%zsz(1)*decomp%zsz(2), a1, decomp%zsz(3), &
                               decomp%zsz(1)*decomp%zsz(2), 1, a1, decomp%zsz(3), &
                               decomp%zsz(1)*decomp%zsz(2), 1, isign, plan_type)
#else
      call sfftw_plan_many_dft(plan1, 1, decomp%zsz(3), &
                               decomp%zsz(1)*decomp%zsz(2), a1, decomp%zsz(3), &
                               decomp%zsz(1)*decomp%zsz(2), 1, a1, decomp%zsz(3), &
                               decomp%zsz(1)*decomp%zsz(2), 1, isign, plan_type)
#endif

      deallocate (a1)

      return
   end subroutine c2c_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D r2c FFTs in X direction
   subroutine r2c_1m_x_plan(plan1, decomp_ph, decomp_sp)

      implicit none

      type(C_PTR), intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

      real(mytype), allocatable, dimension(:, :, :) :: a1
      complex(mytype), allocatable, dimension(:, :, :) :: a2

      allocate (a1(decomp_ph%xsz(1), decomp_ph%xsz(2), decomp_ph%xsz(3)))
      allocate (a2(decomp_sp%xsz(1), decomp_sp%xsz(2), decomp_sp%xsz(3)))
#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%xsz(1), &
                                   decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
                                   decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
                                   plan_type)
#else
      call sfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%xsz(1), &
                                   decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
                                   decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
                                   plan_type)
#endif
      deallocate (a1, a2)

      return
   end subroutine r2c_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
   subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

      implicit none

      type(C_PTR), intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

      complex(mytype), allocatable, dimension(:, :, :) :: a1
      real(mytype), allocatable, dimension(:, :, :) :: a2

      allocate (a1(decomp_sp%xsz(1), decomp_sp%xsz(2), decomp_sp%xsz(3)))
      allocate (a2(decomp_ph%xsz(1), decomp_ph%xsz(2), decomp_ph%xsz(3)))
#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%xsz(1), &
                                   decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
                                   decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
                                   plan_type)
#else
      call sfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%xsz(1), &
                                   decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
                                   decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
                                   plan_type)
#endif
      deallocate (a1, a2)

      return
   end subroutine c2r_1m_x_plan

   ! Return a FFTW3 plan for multiple 1D r2c FFTs in Z direction
   subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

      implicit none

      type(C_PTR), intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

      real(mytype), allocatable, dimension(:, :, :) :: a1
      complex(mytype), allocatable, dimension(:, :, :) :: a2

      allocate (a1(decomp_ph%zsz(1), decomp_ph%zsz(2), decomp_ph%zsz(3)))
      allocate (a2(decomp_sp%zsz(1), decomp_sp%zsz(2), decomp_sp%zsz(3)))
#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
                                   decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
#else
      call sfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
                                   decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
#endif
      deallocate (a1, a2)

      return
   end subroutine r2c_1m_z_plan

   ! Return a FFTW3 plan for multiple 1D c2r FFTs in Z direction
   subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph)

      implicit none

      type(C_PTR), intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

      complex(mytype), allocatable, dimension(:, :, :) :: a1
      real(mytype), allocatable, dimension(:, :, :) :: a2

      allocate (a1(decomp_sp%zsz(1), decomp_sp%zsz(2), decomp_sp%zsz(3)))
      allocate (a2(decomp_ph%zsz(1), decomp_ph%zsz(2), decomp_ph%zsz(3)))

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
                                   decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, plan_type)
#else
      call sfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
                                   decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
                                   decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, plan_type)
#endif
      deallocate (a1, a2)

      return
   end subroutine c2r_1m_z_plan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time initialisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_fft_engine

      implicit none

      call decomp_2d_fft_log("FFTW (version 3.x)")

      if (format == PHYSICAL_IN_X) then

         ! For C2C transforms
         call c2c_1m_x_plan(plan(-1, 1), ph, FFTW_FORWARD)
         call c2c_1m_y_plan(plan(-1, 2), ph, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(-1, 3), ph, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(1, 3), ph, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(1, 2), ph, FFTW_BACKWARD)
         call c2c_1m_x_plan(plan(1, 1), ph, FFTW_BACKWARD)

         ! For R2C/C2R tranforms
         call r2c_1m_x_plan(plan(0, 1), ph, sp)
         call c2c_1m_y_plan(plan(0, 2), sp, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(0, 3), sp, FFTW_FORWARD)
         call c2c_1m_z_plan(plan(2, 3), sp, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(2, 2), sp, FFTW_BACKWARD)
         call c2r_1m_x_plan(plan(2, 1), sp, ph)

      else if (format == PHYSICAL_IN_Z) then

         ! For C2C transforms
         call c2c_1m_z_plan(plan(-1, 3), ph, FFTW_FORWARD)
         call c2c_1m_y_plan(plan(-1, 2), ph, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(-1, 1), ph, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(1, 1), ph, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(1, 2), ph, FFTW_BACKWARD)
         call c2c_1m_z_plan(plan(1, 3), ph, FFTW_BACKWARD)

         ! For R2C/C2R tranforms
         call r2c_1m_z_plan(plan(0, 3), ph, sp)
         call c2c_1m_y_plan(plan(0, 2), sp, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(0, 1), sp, FFTW_FORWARD)
         call c2c_1m_x_plan(plan(2, 1), sp, FFTW_BACKWARD)
         call c2c_1m_y_plan(plan(2, 2), sp, FFTW_BACKWARD)
         call c2r_1m_z_plan(plan(2, 3), sp, ph)

      end if

      return
   end subroutine init_fft_engine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time finalisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine finalize_fft_engine

      implicit none

      integer :: i, j

      do j = 1, 3
         do i = -1, 2
#ifdef DOUBLE_PREC
            call dfftw_destroy_plan(plan(i, j))
#else
            call sfftw_destroy_plan(plan(i, j))
#endif
         end do
      end do

      call fftw_cleanup()

      return
   end subroutine finalize_fft_engine

   ! Following routines calculate multiple one-dimensional FFTs to form
   ! the basis of three-dimensional FFTs.

   ! c2c transform, multiple 1D FFTs in x direction
   subroutine c2c_1m_x(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR), intent(IN) :: plan1

#ifdef DOUBLE_PREC
      call dfftw_execute_dft(plan1, inout, inout)
#else
      call sfftw_execute_dft(plan1, inout, inout)
#endif

      return
   end subroutine c2c_1m_x

   ! c2c transform, multiple 1D FFTs in y direction
   subroutine c2c_1m_y(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR), intent(IN) :: plan1

      integer :: k, s3

      ! transform on one Z-plane at a time
      s3 = size(inout, 3)
      do k = 1, s3
#ifdef DOUBLE_PREC
         call dfftw_execute_dft(plan1, inout(:, :, k), inout(:, :, k))
#else
         call sfftw_execute_dft(plan1, inout(:, :, k), inout(:, :, k))
#endif
      end do

      return
   end subroutine c2c_1m_y

   ! c2c transform, multiple 1D FFTs in z direction
   subroutine c2c_1m_z(inout, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR), intent(IN) :: plan1

#ifdef DOUBLE_PREC
      call dfftw_execute_dft(plan1, inout, inout)
#else
      call sfftw_execute_dft(plan1, inout, inout)
#endif

      return
   end subroutine c2c_1m_z

   ! r2c transform, multiple 1D FFTs in x direction
   subroutine r2c_1m_x(input, output)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT)  ::  input
      complex(mytype), dimension(:, :, :), intent(OUT) :: output

#ifdef DOUBLE_PREC
      call dfftw_execute_dft_r2c(plan(0, 1), input, output)
#else
      call sfftw_execute_dft_r2c(plan(0, 1), input, output)
#endif

      return

   end subroutine r2c_1m_x

   ! r2c transform, multiple 1D FFTs in z direction
   subroutine r2c_1m_z(input, output)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT)  ::  input
      complex(mytype), dimension(:, :, :), intent(OUT) :: output

#ifdef DOUBLE_PREC
      call dfftw_execute_dft_r2c(plan(0, 3), input, output)
#else
      call sfftw_execute_dft_r2c(plan(0, 3), input, output)
#endif

      return

   end subroutine r2c_1m_z

   ! c2r transform, multiple 1D FFTs in x direction
   subroutine c2r_1m_x(input, output)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT)  ::  input
      real(mytype), dimension(:, :, :), intent(OUT) :: output

#ifdef DOUBLE_PREC
      call dfftw_execute_dft_c2r(plan(2, 1), input, output)
#else
      call sfftw_execute_dft_c2r(plan(2, 1), input, output)
#endif

      return

   end subroutine c2r_1m_x

   ! c2r transform, multiple 1D FFTs in z direction
   subroutine c2r_1m_z(input, output)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: input
      real(mytype), dimension(:, :, :), intent(OUT) :: output

#ifdef DOUBLE_PREC
      call dfftw_execute_dft_c2r(plan(2, 3), input, output)
#else
      call sfftw_execute_dft_c2r(plan(2, 3), input, output)
#endif

      return

   end subroutine c2r_1m_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D FFT - complex to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2c(in, out, isign)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: in
      complex(mytype), dimension(:, :, :), intent(OUT) :: out
      integer, intent(IN) :: isign

#ifndef OVERWRITE
      complex(mytype), allocatable, dimension(:, :, :) :: wk1
#endif

      if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_FORWARD .OR. &
          format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_BACKWARD) then

         ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
         call c2c_1m_x(in, plan(isign, 1))
#else
         allocate (wk1(ph%xsz(1), ph%xsz(2), ph%xsz(3)))
         wk1 = in
         call c2c_1m_x(wk1, plan(isign, 1))
#endif

         ! ===== Swap X --> Y; 1D FFTs in Y =====

         if (dims(1) > 1) then
#ifdef OVERWRITE
            call transpose_x_to_y(in, wk2_c2c, ph)
#else
            call transpose_x_to_y(wk1, wk2_c2c, ph)
#endif
            call c2c_1m_y(wk2_c2c, plan(isign, 2))
         else
#ifdef OVERWRITE
            call c2c_1m_y(in, plan(isign, 2))
#else
            call c2c_1m_y(wk1, plan(isign, 2))
#endif
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_c2c, out, ph)
         else
#ifdef OVERWRITE
            call transpose_y_to_z(in, out, ph)
#else
            call transpose_y_to_z(wk1, out, ph)
#endif
         end if
         call c2c_1m_z(out, plan(isign, 3))

      else if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_BACKWARD &
               .OR. &
               format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_FORWARD) then

         ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
         call c2c_1m_z(in, plan(isign, 3))
#else
         allocate (wk1(ph%zsz(1), ph%zsz(2), ph%zsz(3)))
         wk1 = in
         call c2c_1m_z(wk1, plan(isign, 3))
#endif

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
#ifdef OVERWRITE
            call transpose_z_to_y(in, wk2_c2c, ph)
#else
            call transpose_z_to_y(wk1, wk2_c2c, ph)
#endif
            call c2c_1m_y(wk2_c2c, plan(isign, 2))
         else  ! out==wk2_c2c if 1D decomposition
#ifdef OVERWRITE
            call transpose_z_to_y(in, out, ph)
#else
            call transpose_z_to_y(wk1, out, ph)
#endif
            call c2c_1m_y(out, plan(isign, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_c2c, out, ph)
         end if
         call c2c_1m_x(out, plan(isign, 1))

      end if

#ifndef OVERWRITE
      deallocate (wk1)
#endif

      return
   end subroutine fft_3d_c2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D forward FFT - real to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_r2c(in_r, out_c, kind0)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT) :: in_r
      complex(mytype), dimension(:, :, :), intent(OUT) :: out_c
      integer, dimension(:), intent(IN), optional :: kind0

      integer, dimension(3) :: kind


      ! ===== transformation type =====
      ! dimension 1, 2, 3 = direction x, y z
      ! 0 = normal fft
      ! 1 = cosine transform
      ! 2 = sine transform

      if(PRESENT(kind0)) then
         kind(1:3) = kind0(1:3)
      else
         kind(1:3) = 0 
      end if

      if( kind(1) * kind(2) * kind(3) /= 0) then
         errorcode = 11
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, 'Invalid input for fft_3d_r2c.')
      end if

      
      if (format == PHYSICAL_IN_X) then

         if (kind(1) /= 0 ) then

            ! ===== 1D DFTs in X; Swap X --> Y  =====
            call r2r_dft_1m_x(in_r, wk1_r2r, kind)
            call transpose_x_to_y(wk1_r2r, wk2_r, sp)

            if (kind(2) /= 0 ) then 

               ! xyz = 110
               ! ===== 1D DFTs in Y =====
               call r2r_dft_1m_y(wk2_r, wk2_r2r, kind)
               ! ===== Swap Y --> Z; D FFTs in Z  =====
               call transpose_y_to_z(wk2_r2r, wk3_r, sp)
               call r2c_1m_z(wk3_r, out_c) 

            else if (kind(2) == 0 )

               if(kind(3) == 0) then  

                  ! xyz = 100
                  ! ===== 1D DFTs in Y =====
                  call r2c_1m_y(wk2_r, wk2_r2c) 
                  ! ===== Swap Y --> Z; 1D FFTs in Z =====
                  call transpose_y_to_z(wk2_r2c, out_c, sp)
                  call c2c_1m_z(out_c, plan(0, 3))

               else if(kind(3) /= 0) then 

                  ! xyz = 101
                  ! ===== Swap Y --> Z; 1D DFTs in Z =====
                  call transpose_y_to_z(wk2_r, wk3_r, sp)
                  call r2r_dft_1m_z(wk3_r, wk3_r2r, kind)
                  ! ===== Swap Z --> Y; 1D FFTs in Y; Swap Y --> Z =====
                  call transpose_z_to_y(wk3_r2r, wk2_r, sp)
                  call r2c_1m_y(wk2_r, wk2_r2c) 
                  call transpose_y_to_z(wk2_r2c, out_c, sp)

               end if

            end if

         else if (kind(1) == 0 ) then

            if(kind(2) /= 0) then

               ! ===== Swap X --> Y; 1D DFTs in Y; Swap Y --> Z ===== 
               call transpose_x_to_y(in_r, in2_r, ph)
               call r2r_dft_1m_y(in2_r, wk2_r2r, kind)
               call transpose_y_to_z(wk2_r2r, wk3_r, sp)

               if(kind(3) == 0) then  

                  ! xyz = 010
                  ! ===== 1D FFTs in Z =====
                  call r2c_1m_z(wk3_r, wk3_r2c) 
                  ! ===== Swap Z --> Y --> X; 1D FFTs in X; Swap X --> Y --> Z ===== 
                  call transpose_z_to_y(wk3_r2c, wk2_r2c, sp)
                  call transpose_y_to_x(wk2_r2c, wk1_r2c, sp)
                  call c2c_1m_x(wk1_r2c, plan(0, 1))
                  call transpose_x_to_y(wk1_r2c, wk2_r2c, sp)
                  call transpose_y_to_z(wk2_r2c, out_c,   sp)

               else if(kind(3) /= 0) then

                  ! xyz = 011
                  ! ===== 1D DFTs in Z =====
                  call r2r_dft_1m_z(wk3_r, wk3_r2r, kind)
                  ! ===== Swap Z --> Y --> X; 1D FFTs in X; Swap X --> Y --> Z ===== 
                  call transpose_z_to_y(wk3_r2r, wk2_r2r, sp)
                  call transpose_y_to_x(wk2_r2r, wk1_r2r, sp)
                  call r2c_1m_x(wk1_r2r, wk1_r2c) 
                  call transpose_x_to_y(wk1_r2c, wk2_r2c, sp)
                  call transpose_y_to_z(wk2_r2c, out_c,   sp)

               end if


            else if(kind(2) == 0) then

               if(kind(3) == 0) then

                  ! xyz = 000, restore to original
                  ! ===== 1D FFTs in X =====
                  call r2c_1m_x(in_r, wk13)
                  ! ===== Swap X --> Y; 1D FFTs in Y =====
                  call transpose_x_to_y(wk13, wk2_r2c, sp)
                  call c2c_1m_y(wk13, plan(0, 2))
                  ! ===== Swap Y --> Z; 1D FFTs in Z =====
                  call transpose_y_to_z(wk13, out_c, sp)
                  call c2c_1m_z(out_c, plan(0, 3))

               else if(kind(3) /= 0) then

                  ! xyz = 001
                  ! ===== Swap X --> Y --> Z; 1D DFTs in Z ===== 
                  call transpose_x_to_y(in_r,  in2_r, ph)
                  call transpose_y_to_z(in2_r, in3_r, ph)
                  call r2r_dft_1m_z(in3_r, wk3_r2r, kind)
                  ! ===== Swap Z --> Y; 1D FFTs in Y =====
                  call transpose_z_to_y(wk3_r2r, wk2_r2r, sp)
                  call r2c_1m_y(wk2_r2r, wk2_r2c)
                  ! ===== Swap Y --> X; 1D FFTs in X; Swap X --> Y --> Z =====
                  call transpose_y_to_x(wk2_r2c, wk1_r2c, sp)
                  call c2c_1m_x(wk1_r2c, plan(0, 1))
                  call transpose_x_to_y(wk1_r2c, wk2_r2c, sp)
                  call transpose_y_to_z(wk2_r2c, out_c,   sp)

               end if
            end if


         end if

      else if (format == PHYSICAL_IN_Z) then

         if (kind(3) /= 0 ) then

            ! ===== 1D DFTs in Z; Swap Z --> Y  =====
            call r2r_dft_1m_z(in_r, wk3_r2r, kind)
            call transpose_z_to_y(wk3_r2r, wk2_r, sp)

            if (kind(2) /= 0 ) then 

               ! zyx = 110
               ! ===== 1D DFTs in Y =====
               call r2r_dft_1m_y(wk2_r, wk2_r2r, kind)
               ! ===== Swap Y --> X; 1D FFTs in X  =====
               call transpose_y_to_x(wk2_r2r, wk1_r, sp)
               call r2c_1m_x(wk1_r, out_c) 

            else if (kind(2) == 0 )

               if(kind(1) == 0) then  

                  ! zyx = 100
                  ! ===== 1D DFTs in Y =====
                  call r2c_1m_y(wk2_r, wk2_r2c) 
                  ! ===== Swap Y --> X; 1D FFTs in X =====
                  call transpose_y_to_x(wk2_r2c, out_c, sp)
                  call c2c_1m_x(out_c, plan(0, 1))

               else if(kind(1) /= 0) then 

                  ! zyx = 101
                  ! ===== Swap Y --> X; 1D DFTs in X =====
                  call transpose_y_to_x(wk2_r, wk1_r, sp)
                  call r2r_dft_1m_x(wk1_r, wk1_r2r, kind)
                  ! ===== Swap X --> Y; 1D FFTs in Y; Swap Y --> X =====
                  call transpose_x_to_y(wk1_r2r, wk2_r, sp)
                  call r2c_1m_y(wk2_r, wk2_r2c) 
                  call transpose_y_to_x(wk2_r2c, out_c, sp)

               end if

            end if

         else if (kind(3) == 0 ) then

            if(kind(2) /= 0) then

               ! ===== Swap Z --> Y; 1D DFTs in Y; Swap Y --> X ===== 
               call transpose_z_to_y(in_r, in2_r, ph)
               call r2r_dft_1m_y(in2_r, wk2_r2r, kind)
               call transpose_y_to_x(wk2_r2r, wk1_r, sp)

               if(kind(1) == 0) then  

                  ! zyx = 010
                  ! ===== 1D FFTs in X =====
                  call r2c_1m_x(wk1_r, wk1_r2c) 
                  ! ===== Swap X --> Y --> Z; 1D FFTs in Z; Swap Z --> Y --> X ===== 
                  call transpose_x_to_y(wk1_r2c, wk2_r2c, sp)
                  call transpose_y_to_z(wk2_r2c, wk3_r2c, sp)
                  call c2c_1m_z(wk3_r2c, plan(0, 3))
                  call transpose_z_to_y(wk3_r2c, wk2_r2c, sp)
                  call transpose_y_to_x(wk2_r2c, out_c,   sp)

               else if(kind(1) /= 0) then

                  ! zxy = 011
                  ! ===== 1D DFTs in X =====
                  call r2r_dft_1m_x(wk1_r, wk1_r2r, kind)
                  ! ===== Swap X --> Y --> Z; 1D FFTs in Z; Swap Z--> Y --> X ===== 
                  call transpose_x_to_y(wk1_r2r, wk2_r2r, sp)
                  call transpose_y_to_x(wk2_r2r, wk3_r2r, sp)
                  call r2c_1m_z(wk3_r2r, wk3_r2c) 
                  call transpose_z_to_y(wk3_r2c, wk2_r2c, sp)
                  call transpose_y_to_x(wk2_r2c, out_c,   sp)

               end if


            else if(kind(2) == 0) then

               if(kind(1) == 0) then

                  ! zyx = 000, restore to original
                  ! ===== 1D FFTs in Z =====
                  call r2c_1m_z(in_r, wk13)
                  ! ===== Swap Z --> Y; 1D FFTs in Y =====
                  call transpose_z_to_y(wk13, wk2_r2c, sp)
                  call c2c_1m_y(wk2_r2c, plan(0, 2))
                  ! ===== Swap Y --> X; 1D FFTs in X =====
                  call transpose_y_to_x(wk2_r2c, out_c, sp)
                  call c2c_1m_x(out_c, plan(0, 1))

               else if(kind(1) /= 0) then

                  ! zyx = 001
                  ! ===== Swap Z --> Y --> X; 1D DFTs in X ===== 
                  call transpose_z_to_y(in_r,  in2_r, ph)
                  call transpose_y_to_x(in2_r, in1_r, ph)
                  call r2r_dft_1m_x(in1_r, wk1_r2r, kind)
                  ! ===== Swap X --> Y; 1D FFTs in Y =====
                  call transpose_x_to_y(wk1_r2r, wk2_r2r, sp)
                  call r2c_1m_y(wk2_r2r, wk2_r2c)
                  ! ===== Swap Y --> Z; 1D FFTs in Z; Swap Z --> Y --> X =====
                  call transpose_y_to_z(wk2_r2c, wk3_r2c, sp)
                  call c2c_1m_z(wk3_r2c, plan(0, 3))
                  call transpose_z_to_y(wk3_r2c, wk2_r2c, sp)
                  call transpose_y_to_x(wk2_r2c, out_c,   sp)

               end if

            end if

         end if
      
      end if

      return
   end subroutine fft_3d_r2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2r(in_c, out_r, kind0)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: in_c
      real(mytype), dimension(:, :, :), intent(OUT) :: out_r
      integer, dimension(:), intent(IN), optional :: kind0

      integer, dimension(3) :: kind


      ! ===== transformation type =====
      ! dimension 1, 2, 3 = direction x, y z
      ! 0 = normal fft
      ! 1 = cosine transform
      ! 2 = sine transform

      if(PRESENT(kind0)) then
         kind(1:3) = kind0(1:3)
      else
         kind(1:3) = 0 
      end if

      if( kind(1) * kind(2) * kind(3) /= 0) then
         errorcode = 11
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, 'Invalid input for fft_3d_r2c.')
      end if

      if (format == PHYSICAL_IN_X) then

         if (kind(3) == 0 ) then

            if (kind(2) == 0 ) then 

               ! ===== 1D FFTs in Z; Swap Z --> Y  =====
               call c2c_1m_z(in_c, plan(2, 3))
               call transpose_z_to_y(in_c, wk2_c2c, sp)

               if(kind(1) == 0) then 

                  ! zyx = 000, restore to the original
                  ! ===== 1D FFTs in Y  =====
                  call c2c_1m_y(wk2_c2c, plan(2, 2))
                  ! =====Swap Y --> X; 1D FFTs in X  =====
                  call transpose_y_to_x(wk2_c2c, wk1_c2c, sp)
                  call c2r_1m_x(wk1_c2c, out_r)

               else if(kind(1) /= 0) then

                  ! zyx = 001
                  ! ===== 1D DFTs in Y  =====
                  call c2r_1m_y(wk2_c2c, wk2_c2r)
                  ! =====Swap Y --> X; 1D FFTs in X  =====
                  call transpose_y_to_x(wk2_c2r, wk1_c2r, sp)
                  call r2r_dft_1m_x(wk1_c2r, out_r)

               end if

            else if (kind(2) /= 0 )

               if(kind(1) == 0) then  

                  ! zyx = 010
                  ! ===== 1D FFTs in Z; Swap Z --> Y --> X =====
                  call c2c_1m_z(in_c, plan(2, 3))
                  call transpose_z_to_y(in_c,    wk2_c2c, sp)
                  call transpose_y_to_x(wk2_c2c, wk1_c2c, sp)
                  ! ===== 1D FFTs in X  =====
                  call c2r_1m_x(wk1_c2c, wk1_c2c)



                  ! ===== Swap Y --> X; 1D FFTs in X =====
                  call transpose_y_to_x(wk2_r2c, out_c, sp)
                  call c2c_1m_x(out_c, plan(0, 1))

               else if(kind(1) /= 0) then 

                  ! zyx = 101
                  ! ===== Swap Y --> X; 1D DFTs in X =====
                  call transpose_y_to_x(wk2_r, wk1_r, sp)
                  call r2r_dft_1m_x(wk1_r, wk1_r2r, kind)
                  ! ===== Swap X --> Y; 1D FFTs in Y; Swap Y --> X =====
                  call transpose_x_to_y(wk1_r2r, wk2_r, sp)
                  call r2c_1m_y(wk2_r, wk2_r2c) 
                  call transpose_y_to_x(wk2_r2c, out_c, sp)

               end if

            end if

         else if (kind(3) /= 0 ) then

            if(kind(2) /= 0) then

               ! ===== Swap Z --> Y; 1D DFTs in Y; Swap Y --> X ===== 
               call transpose_z_to_y(in_r, in2_r, ph)
               call r2r_dft_1m_y(in2_r, wk2_r2r, kind)
               call transpose_y_to_x(wk2_r2r, wk1_r, sp)

               if(kind(1) == 0) then  

                  ! zyx = 010
                  ! ===== 1D FFTs in X =====
                  call r2c_1m_x(wk1_r, wk1_r2c) 
                  ! ===== Swap X --> Y --> Z; 1D FFTs in Z; Swap Z --> Y --> X ===== 
                  call transpose_x_to_y(wk1_r2c, wk2_r2c, sp)
                  call transpose_y_to_z(wk2_r2c, wk3_r2c, sp)
                  call c2c_1m_z(wk3_r2c, plan(0, 3))
                  call transpose_z_to_y(wk3_r2c, wk2_r2c, sp)
                  call transpose_y_to_x(wk2_r2c, out_c,   sp)

               else if(kind(1) /= 0) then

                  ! zxy = 011
                  ! ===== 1D DFTs in X =====
                  call r2r_dft_1m_x(wk1_r, wk1_r2r, kind)
                  ! ===== Swap X --> Y --> Z; 1D FFTs in Z; Swap Z--> Y --> X ===== 
                  call transpose_x_to_y(wk1_r2r, wk2_r2r, sp)
                  call transpose_y_to_x(wk2_r2r, wk3_r2r, sp)
                  call r2c_1m_z(wk3_r2r, wk3_r2c) 
                  call transpose_z_to_y(wk3_r2c, wk2_r2c, sp)
                  call transpose_y_to_x(wk2_r2c, out_c,   sp)

               end if


            else if(kind(2) == 0) then

               if(kind(1) == 0) then

                  ! zyx = 000, restore to original
                  ! ===== 1D FFTs in Z =====
                  call r2c_1m_z(in_r, wk13)
                  ! ===== Swap Z --> Y; 1D FFTs in Y =====
                  call transpose_z_to_y(wk13, wk2_r2c, sp)
                  call c2c_1m_y(wk2_r2c, plan(0, 2))
                  ! ===== Swap Y --> X; 1D FFTs in X =====
                  call transpose_y_to_x(wk2_r2c, out_c, sp)
                  call c2c_1m_x(out_c, plan(0, 1))

               else if(kind(1) /= 0) then

                  ! zyx = 001
                  ! ===== Swap Z --> Y --> X; 1D DFTs in X ===== 
                  call transpose_z_to_y(in_r,  in2_r, ph)
                  call transpose_y_to_x(in2_r, in1_r, ph)
                  call r2r_dft_1m_x(in1_r, wk1_r2r, kind)
                  ! ===== Swap X --> Y; 1D FFTs in Y =====
                  call transpose_x_to_y(wk1_r2r, wk2_r2r, sp)
                  call r2c_1m_y(wk2_r2r, wk2_r2c)
                  ! ===== Swap Y --> Z; 1D FFTs in Z; Swap Z --> Y --> X =====
                  call transpose_y_to_z(wk2_r2c, wk3_r2c, sp)
                  call c2c_1m_z(wk3_r2c, plan(0, 3))
                  call transpose_z_to_y(wk3_r2c, wk2_r2c, sp)
                  call transpose_y_to_x(wk2_r2c, out_c,   sp)

               end if

            end if

         end if

      else if (format == PHYSICAL_IN_Z) then


      end if

      



      

   end subroutine

      































#ifndef OVERWRITE
      complex(mytype), allocatable, dimension(:, :, :) :: wk1
#endif

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
         call c2c_1m_z(in_c, plan(2, 3))
#else
         allocate (wk1(sp%zsz(1), sp%zsz(2), sp%zsz(3)))
         wk1 = in_c
         call c2c_1m_z(wk1, plan(2, 3))
#endif

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
         call transpose_z_to_y(in_c, wk2_r2c, sp)
#else
         call transpose_z_to_y(wk1, wk2_r2c, sp)
#endif
         call c2c_1m_y(wk2_r2c, plan(2, 2))

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_r2c, wk13, sp)
            call c2r_1m_x(wk13, out_r)
         else
            call c2r_1m_x(wk2_r2c, out_r)
         end if

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
         call c2c_1m_x(in_c, plan(2, 1))
#else
         allocate (wk1(sp%xsz(1), sp%xsz(2), sp%xsz(3)))
         wk1 = in_c
         call c2c_1m_x(wk1, plan(2, 1))
#endif

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
#ifdef OVERWRITE
            call transpose_x_to_y(in_c, wk2_r2c, sp)
#else
            call transpose_x_to_y(wk1, wk2_r2c, sp)
#endif
            call c2c_1m_y(wk2_r2c, plan(2, 2))
         else  ! in_c==wk2_r2c if 1D decomposition
#ifdef OVERWRITE
            call c2c_1m_y(in_c, plan(2, 2))
#else
            call c2c_1m_y(wk1, plan(2, 2))
#endif
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_r2c, wk13, sp)
         else
#ifdef OVERWRITE
            call transpose_y_to_z(in_c, wk13, sp)
#else
            call transpose_y_to_z(wk1, wk13, sp)
#endif
         end if
         call c2r_1m_z(wk13, out_r)

      end if

#ifndef OVERWRITE
      deallocate (wk1)
#endif

      return
   end subroutine fft_3d_c2r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D DFT - real to real
   ! condition: all 3-D with non-periodic b.c.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_r2r(in_r, out_c, iforward, kind0)

      implicit none

      real(mytype), dimension(:, :, :), intent(INOUT) :: in_r
      complex(mytype), dimension(:, :, :), intent(OUT) :: out_c
      integer, intent(IN) :: iforward
      integer, dimension(:), intent(IN), optional :: kind0

      integer, dimension(3) :: kind


      ! ===== transformation type =====
      ! dimension 1, 2, 3 = direction x, y z
      ! 0 = normal fft
      ! 1 = cosine transform
      ! 2 = sine transform
      if(PRESENT(kind0)) then
         kind(1:3) = kind0(1:3)
      else
         kind(1:3) = 0 
      end if

      if( kind(1) * kind(2) * kind(3) == 0) then
         errorcode = 11
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, 'Invalid input for fft_3d_r2r.'
      end if

      
      if ( (format == PHYSICAL_IN_X .and. iforward == +1 ) .or. &
           (format == PHYSICAL_IN_Z .and. iforward == -1)) then
      
         ! ===== 1D DFTs in X =====
         call r2r_dft_1m_x(in_r, wk1_r, kind)

         ! ===== Swap X --> Y; 1D DFTs in Y =====
         call transpose_x_to_y(wk1_r, wk2_r, sp)
         call r2r_dft_1m_y(wk2_r, wk22_r, kind)

         !===== Swap Y --> Z, 1D DFTs in Z   =====
         call transpose_y_to_z(wk22_r, wk3_r, sp)
         call r2r_dft_1m_z(wk3_r, wk33_r, kind) 

      else if ( (format == PHYSICAL_IN_Z .and. iforward == +1 ) .or. &
                (format == PHYSICAL_IN_X .and. iforward == -1)) then

         ! ===== 1D DFTs in Z =====
         call r2r_dft_1m_z(in_r, wk3_r, kind)

         ! ===== Swap Z --> Y; 1D DFTs in Y =====
         call transpose_z_to_y(wk3_r, wk2_r, sp)
         call r2r_dft_1m_y(wk2_r, wk22_r, kind)

         !===== Swap Y --> X, 1D DFTs in X   =====
         call transpose_y_to_x(wk22_r, wk1_r, sp)
         call r2r_dft_1m_x(wk1_r, wk11_r, kind) 
      else
         errorcode = 12
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, 'Invalid input for fft_3d_r2r.'
      end if


      return
   end subroutine fft_3d_r2r


   

end module decomp_2d_fft
