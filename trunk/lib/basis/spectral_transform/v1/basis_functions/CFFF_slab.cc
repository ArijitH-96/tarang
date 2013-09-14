/* Tarang-2
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-2 .
 *
 * Tarang-2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-2; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */
/*! \file fourier.cc 
 * 
 * @sa fourier.h
 * 
 * @author  M. K. Verma
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		ArrayIFFT
 */ 


#include "spectral_transform.h"

//*********************************************************************************************
void SpectralTransform::Init_CFFF_SLAB()
{
    if (Ny > 1) {
        local_Nx = Nx/2+1;
        local_Ny = Ny/numprocs;
        local_Nz = Nz;
        
        local_Nx_start = 0;
        local_Ny_start = my_id * local_Ny;
        local_Nz_start = 0;
    }
    else if (Ny == 1) {
        local_Nx = Nx/2+1;
        local_Ny = 1;
        local_Nz = Nz/numprocs;
        
        local_Nx_start = 0;
        local_Ny_start = 0;
        local_Nz_start = my_id * local_Nz;
    }
    
    if (Ny>1) {
        X_3d.resize(local_Ny,Nz,Nx);
        Xr_3d.resize(local_Ny,Nz,2*Nx);
    }
    else if (Ny==1) {
        X_2d.resize(local_Nz,Nx);
        Xr_2d.resize(local_Nz,2*Nx);
    } 
    
	
    //Initialize plans
    
    if (Ny > 1) {
        c2c_xyz_forward_plan = FFTW_MPI_PLAN_DFT_3D_DP(Ny, Nz, Nx,reinterpret_cast<FFTW_COMPLEX_DP*>(Xr_3d.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(X_3d.data()), MPI_COMM_WORLD, FFTW_FORWARD, FFTW_PLAN_FLAG);
        
        c2c_xyz_inverse_plan = FFTW_MPI_PLAN_DFT_3D_DP(Ny, Nz, Nx,reinterpret_cast<FFTW_COMPLEX_DP*>(X_3d.data()), reinterpret_cast<FFTW_COMPLEX_DP*>(Xr_3d.data()), MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_PLAN_FLAG);

    }
    
    else if (Ny == 1) {
        /*	r2c_xz_plan = FFTW_MPI_PLAN_DFT_R2C_2D_DP(Nx, Nz, reinterpret_cast<DP*>(plane_xz.data()),reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG);
         
         c2r_xz_plan = FFTW_MPI_PLAN_DFT_C2R_2D_DP(Nx, Nz,reinterpret_cast<FFTW_COMPLEX_DP*>(plane_xz.data()), reinterpret_cast<DP*>(plane_xz.data()),MPI_COMM_WORLD, FFTW_PLAN_FLAG); */
    }
};
//**********************************************************************************************

void SpectralTransform::Zero_pad_last_plane_CFFF_SLAB(Array<complx,3> Ar)
{
}

//*********************************************************************************************

void SpectralTransform::Norm_CFFF_SLAB(Array<complx,3> A)
{
	A = A/(DP(Nx) * DP(Ny) * DP(Nz));
}


//*********************************************************************************************

// Ar: yxz, A: xyz
void SpectralTransform::Forward_transform_CFFF_SLAB(Array<complx,3> Ar, Array<complx,3> A)
{
    FTc2c_xyz(Ar, A);
	Norm_CFFF_SLAB(A);
}

void SpectralTransform::Forward_transform_CFFF_SLAB(Array<complx,2> Ar, Array<complx,2> A)
{
    //FTc2c_xz(Ar, A);
	//Norm_CFFF_SLAB(A);
}


//*********************************************************************************************


void SpectralTransform::Inverse_transform_CFFF_SLAB(Array<complx,3> A, Array<complx,3> Ar)
{
    IFTc2c_xyz(A, Ar);
}

void SpectralTransform::Inverse_transform_CFFF_SLAB(Array<complx,2> A, Array<complx,2> Ar)
{
    //IFTc2c_xz(A, Ar);
}


//********************************	End of four_tr.cc *****************************************


