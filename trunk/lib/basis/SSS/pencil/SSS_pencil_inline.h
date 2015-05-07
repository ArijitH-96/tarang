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

/*! \file sincosfour_inline.h 
 * 
 * @brief Inline functions to compute array indices \f$ \vec{i} \f$ given wavenumber  
 *		\f$ \vec{k} \f$ and viceversa.
 *		Also contains other useful functions like Kmagnitude, Max radius etc.
 *
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the grid index 
 *		\f$ \vec{i} \f$. Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber 
 *		\f$ \vec{k} \f$ using \f$ K_i = k_i * f_i \f$  where \f$ f_i \f$ is the kfactor[i]. 
 *
 * lx, ly, lz = local array indices. <BR>
 * i1, i2, i3 = global or total array indices.
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 *		or gridwavenumber depending on the switch WAVENOACTUAL or WAVENOGRID.
 * 
 * @author  M. K. Verma
 * @version 4.0 Parallel
 * @date	August 2008
 * @bug		No known bugs
 */


#include "SSS_pencil.h"

//*********************************************************************************************

/*! @brief	Get grid waveno kx given first local array index lx.
 * 
 *  if i1 <=Ni/2, ki=i1; else kx=i1-Ni.  <BR>
 *	i1= local_N1_start + lx.
 * 
 * \param lx  first local index of an array

 */
inline int SSS_PENCIL::Get_kx(int lx) {return  (lx <= Nx/2) ? lx : (lx-Nx); } 
 

/*! @brief	Get local array index lx given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 *	i1= local_N1_start + lx.
 * 
 * \param	kx  grid wavenumber along x
 * \return	lx  local array index along x
 */
inline int SSS_PENCIL::Get_lx(int kx) {return (kx >= 0) ? kx : (kx + Nx); } 


/*! @brief	Get array index ix given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 * 
 * \param	kx  grid wavenumber along x
 * \return	i1  array index along x
 */
inline int SSS_PENCIL::Get_ix(int kx) {return  (kx >= 0) ? kx : (kx + Nx);  }		


/*! @brief	Get grid waveno ky given first local array index ly.
 * 
 *  If ly<=N2/2, k2=ly; else ky=ly-N2.
 * 
 * \param ly  second local index of an array
 * \return kx corresponding to lx
 */
inline int SSS_PENCIL::Get_ky(int ly) {return  ((ly_start + ly) <= Ny/2) ? (ly_start + ly) : (ly_start + ly-Ny); }

/*! @brief	Get local array index ly given grid waveno ky.
 * 
 *  If ky>=0, ly=ky; else ly=ky+N2.
 * 
 * \param	ky  grid wavenumber along y
 * \return	ly  local array index along y
 */
inline int SSS_PENCIL::Get_ly(int ky) { return  (ky >= 0) ? (ky-ly_start) : (ky + Ny-ly_start);  } 

inline int SSS_PENCIL::Get_iy(int ky) { return  (ky >= 0) ? ky : (ky + Ny); }

inline int SSS_PENCIL::Get_kz(int lz)  {return lz + lz_start;} 

inline int SSS_PENCIL::Get_lz(int kz)  {return kz - lz_start;}

inline int SSS_PENCIL::Get_iz(int kz)  {return kz;}

inline bool SSS_PENCIL::Probe_in_me(int kx, int ky, int kz) 
{  
	int ly = Get_ly(ky);
	int lz = Get_lz(kz);
	
	return ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) );
}



inline Complex SSS_PENCIL::Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> A)
{
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) )
		return A_real(lx, ly, lz);

	return 0;
}

inline TinyVector<Complex,3> SSS_PENCIL::Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	Array<Real,3> Ax_real(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ay_real(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Az_real(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) )
		return TinyVector<Complex,3>(Ax_real(lx, ly, lz), Ay_real(lx, ly, lz), Az_real(lx, ly, lz));

	return TinyVector<Complex,3>(0,0,0);
}


// ASSIGN...
inline void SSS_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Complex field)
{
	cout << "MYERROR: SSS_PENCIL::Assign_spectral_field(); Use real data type " << endl;
}

inline void SSS_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
	cout << "MYERROR: SSS_PENCIL::Assign_spectral_field(); Use real data type " << endl;
}

inline void SSS_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Real field)
{
	
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) )
		A_real(lx, ly, lz) = field;
}

inline void SSS_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	
	Array<Real,3> Ax_real(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ay_real(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Az_real(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) ) {
		Ax_real(lx, ly, lz) = V(0);
        Ay_real(lx, ly, lz) = V(1);
        Az_real(lx, ly, lz) = V(2);
	}
}


inline void SSS_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Complex field)
{
	cout << "MYERROR: SSS_PENCIL::Add_spectral_field(); Use real data type " << endl;
}

inline void SSS_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
	cout << "MYERROR: SSS_PENCIL::Add_spectral_field(); Use real data type " << endl;
}

inline void SSS_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A,Real field)
{
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) )
		A_real(lx, ly, lz) += field;
}

inline void SSS_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	
	Array<Real,3> Ax_real(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ay_real(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Az_real(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) ) {
		Ax_real(lx, ly, lz) += V(0);
        Ay_real(lx, ly, lz) += V(1);
        Az_real(lx, ly, lz) += V(2);
	}
}


// LOCAL...
inline Complex SSS_PENCIL::Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A)
{
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	return A_real(lx, ly, lz);
}

inline TinyVector<Complex,3> SSS_PENCIL::Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	Array<Real,3> Ax_real(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ay_real(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Az_real(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	return TinyVector<Complex,3>(Ax_real(lx, ly, lz), Ay_real(lx, ly, lz), Az_real(lx, ly, lz));
}


// ASSIGN...
inline void SSS_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Complex field)
{
	cout << "MYERROR: SSS_PENCIL::Assign_local_spectral_field(); Use real data type " << endl;
}

inline void SSS_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
	cout << "MYERROR: SSS_PENCIL::Assign_local_spectral_field(); Use real data type " << endl;
}

inline void SSS_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Real field)
{
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) )
		A_real(lx, ly, lz) = field;
}

inline void SSS_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	Array<Real,3> Ax_real(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ay_real(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Az_real(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) ) {
		Ax_real(lx, ly, lz) = V(0);
        Ay_real(lx, ly, lz) = V(1);
        Az_real(lx, ly, lz) = V(2);
	}
}


inline void SSS_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Complex field)
{
	cout << "MYERROR: SSS_PENCIL::Add_local_spectral_field(); Use real data type " << endl;
	
}

inline void SSS_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
	
	cout << "MYERROR: SSS_PENCIL::Add_local_spectral_field(); Use real data type " << endl;
}

inline void SSS_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Real field)
{
	Array<Real,3> A_real(reinterpret_cast<Real*>(A.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) )
		A_real(lx, ly, lz) += field;
}

inline void SSS_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	
	Array<Real,3> Ax_real(reinterpret_cast<Real*>(Ax.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Ay_real(reinterpret_cast<Real*>(Ay.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	Array<Real,3> Az_real(reinterpret_cast<Real*>(Az.data()), shape_complex_array*shape(1,1,2), neverDeleteData);
	
	if ( ((ly >= 0) && (ly < maxly)) && ((lz >= 0) && (lz < maxlz)) ) {
		Ax_real(lx, ly, lz) += V(0);
        Ay_real(lx, ly, lz) += V(1);
        Az_real(lx, ly, lz) += V(2);
	}
}


// REAL - SPACE FNS
inline int SSS_PENCIL::Get_lx_real_space(int rx)  {return  rx - my_x_pcoord_real*maxrx;}

inline int SSS_PENCIL::Get_ly_real_space(int ry)  {return  ry - my_y_pcoord_real*maxry;}

inline int SSS_PENCIL::Get_lz_real_space(int rz) {return rz;}


inline int SSS_PENCIL::Get_rx_real_space(int lx)  {return  lx + my_x_pcoord_real*maxrx;}

inline int SSS_PENCIL::Get_ry_real_space(int ly)  {return  ly + my_y_pcoord_real*maxry; }		

inline int SSS_PENCIL::Get_rz_real_space(int lz) {return lz;}



inline bool SSS_PENCIL::Probe_in_me_real_space(int rx, int ry, int rz) 
{
	int lx = Get_lx_real_space(rx);
	int ly = Get_ly_real_space(ry);
	
	return ( ((lx >= 0) && (lx < maxrx)) && ((ly >= 0) && (ly < maxry)) );
}


// lz is the coord of the complex array
inline Real SSS_PENCIL::Get_real_field(int rx, int ry, int rz, Array<Real,3> A)
{	
	int lx = Get_lx_real_space(rx);
	int ly = Get_ly_real_space(ry);
	int lz = Get_lz_real_space(rz);
	
	if ( ((lx >= 0) && (lx < maxrx)) && ((ly >= 0) && (ly < maxry)) )
		return (A(lx, ly, lz));

	return 0;
}

inline TinyVector<Real,3> SSS_PENCIL::Get_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az)
{
	
	int lx = Get_lx_real_space(rx);
	int ly = Get_ly_real_space(ry);
	int lz = Get_lz_real_space(rz);
	
	if ( ((lx >= 0) && (lx < maxrx)) && ((ly >= 0) && (ly < maxry)) )
		return TinyVector<Real,3>(Ax(lx, ly, lz), Ay(lx, ly, lz), Az(lx, ly, lz));
	
	return TinyVector<Real,3>(0,0,0);
}


// lz is the coord of the complex array
inline void SSS_PENCIL::Assign_real_field(int rx, int ry, int rz, Array<Real,3> A, Real field)
{	
	int lx = Get_lx_real_space(rx);
	int ly = Get_ly_real_space(ry);
	int lz = Get_lz_real_space(rz);
	
	if ( ((lx >= 0) && (lx < maxrx)) && ((ly >= 0) && (ly < maxry)) )
		A(lx, ly, lz) = field;
}

inline void SSS_PENCIL::Assign_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az, TinyVector<Real,3> V)
{
	
	int lx = Get_lx_real_space(rx);
	int ly = Get_ly_real_space(ry);
	int lz = Get_lz_real_space(rz);
	
	if ( ((lx >= 0) && (lx < maxrx)) && ((ly >= 0) && (ly < maxry)) ) {
		Ax(lx, ly, lz) = V(0);
		Ay(lx, ly, lz) = V(1);
		Az(lx, ly, lz) = V(2);
	}    
	
}


/**********************************************************************************************
 
 Compute Wavenumber
 
 ***********************************************************************************************/


inline void SSS_PENCIL::Wavenumber(int lx, int ly, int lz, TinyVector<Real,3> &K)
{
	K = Get_kx(lx)*kfactor[1],  Get_ky(ly)*kfactor[2], Get_kz(lz)*kfactor[3];
}


	// Complex K; The imaginary part is zero.  Written to use cross function of blitz.
	// Omega = cross(V,K).
inline void SSS_PENCIL::Wavenumber(int lx, int ly, int lz, TinyVector<Complex,3> &K)
{
	K = Complex(Get_kx(lx)*kfactor[1], 0.0), Complex(Get_ky(ly)*kfactor[2], 0.0), Complex(Get_kz(lz)*kfactor[3], 0.0);
}


/**********************************************************************************************

		  If wavenos computed using actual wavenumber:  Ki = Kfactor[i]*grid[i]

***********************************************************************************************/


///  WAVENOACTUAL: \f$ K = \sqrt{K_x^2 + K_y^2 + K_z^2} \f$
inline Real SSS_PENCIL::Kmagnitude(int lx, int ly, int lz)
{ 
	if	(global.field.waveno_switch)
		return sqrt( pow2(Get_kx(lx)*kfactor[1]) + pow2(Get_ky(ly)*kfactor[2]) + pow2(Get_kz(lz)*kfactor[3]) ); 
	
	else 
		return sqrt( pow2(Get_kx(lx)) + pow2(Get_ky(ly)) + pow2(Get_kz(lz)) );
}


/// WAVENOACTUAL -- Radius of the smallest sphere that contains the wavenumber K box Ni's. <BR>
/// The range of kx=[0:N1-1];  ki = [-Ni/2+1 : Ni/2] along perp directions.
inline int SSS_PENCIL::Min_radius_outside() 
{
	if	(global.field.waveno_switch)
		return (int) ceil(sqrt(  pow2((Nx-1)*kfactor[1]) + pow2((Ny-1)*kfactor[2]) + pow2((Nz-1)*kfactor[3]) ));
	
	else 
		return (int) ceil(sqrt(  pow2(Nx-1) + pow2((Ny-1)) + pow2((Nz-1)) ));
}


/// WAVENOACTUAL -- Radius of the largest sphere that fits inside the wavenumber K box Ni's. <BR>
/// The range of kx=[0:N1-1];  ki = [-Ni/2+1 : Ni/2] along perp directions.
inline int SSS_PENCIL::Max_radius_inside() 
{
	int ans = 1;
	Real Kmag;
	
	if	(global.field.waveno_switch) {
        Kmag = min( (Nx-1)*kfactor[1], (Ny-1)*kfactor[2]);
		Kmag = min(Kmag, ((Ny-1))*kfactor[3]); 

		ans = ((int) Kmag);
	}
	
	else {
        ans = min(Nx-1, (Ny-1));
		ans = min(ans, (Nz-1));
	}
	
	return ans;
}		


/*! \brief WAVENOACTUAL -- Returns the approximate number of modes in a wavenumber 
 *			K shell of radius "radius".
 * 
 *  We divide the area in Fourier space by volume of unit lattice \f$ \Pi_i f_i \f$, 
 *			where \f$ f_i \f$ is the factor[i]. In 1D, The area in NOT divided by kfactor.
 *
 * \param  radius
 * \return The number of modes in a shell of radius. In 2D, it is quarter circle (kx, ky>= 0). 
 *			In 3D, it is (1/8) sphere with (kx,kz>=0).
 */
inline Real SSS_PENCIL::Approx_number_modes_in_shell(int radius)
{
	if (global.field.waveno_switch)
        return (4*M_PI*radius*radius)/(kfactor[1]*kfactor[2]*kfactor[3]);	
    
    else 
        return (4*M_PI*radius*radius);
}


//*********************************************************************************************



/*! \brief Returns multiplication factor for computing enregy spectrum etc. 
 *  
 * Modal energy  = \f$ E(k) =  |A(k_x, \vec{k}_{\perp})|^{2} \f$  if  \f$ (k_x > 0) \f$.  <BR>
 * Modal energy  = \f$ E(k) =  |A(0, \vec{k}_{\perp})|^{2}/2 \f$  if  \f$ (k_x = 0) \f$.  <BR>
 * In simulation we double the energy of most of the modes because complex conjugates modes with 
 * -ky are not stored in the simulation.  
 * The modes on the xy-plane are not doubled because their c.c. are already counted.
 * 
 * \param  lx, ly, lz
 * \return Multiplication factor for computing enregy spectrum etc.
 */
inline Real SSS_PENCIL::Multiplicity_factor(int lx, int ly, int lz)
{

	int kx = Get_kx(lx);
	
    // All zeros
    if ((kx==0) && (ly==0) && (lz==0))
        return 0.5;
        
    // All nonzero
    else if ((kx>0) && (ly>0) && (lz>0))
        return 4.0;
    
    // Any two are zero
    else if ((kx*ly==0) && (ly*lz==0) && (lz*kx==0))
        return 1.0;
    
    // Only one is zero
    else 
        return 2.0;
}



/**********************************************************************************************

		Modal energy

***********************************************************************************************/
/// Modal energy  = \f$ E(k) =  |A(k_x, \vec{k}_{\perp})|^{2} \f$  if  \f$ (k_x > 0) \f$.  <BR>
/// Modal energy  = \f$ E(k) =  |A(0, \vec{k}_{\perp})|^{2}/2 \f$  if  \f$ (k_x = 0) \f$.  <BR>
inline Real SSS_PENCIL::Modal_energy(int lx, int ly, int lz, Array<Complex,3> A)
{
    
    if (lz%2 == 1)
        return pow2(imag(A(lx, ly, lz/2)));
    
    else
        return pow2(real(A(lx, ly, lz/2)));
}

/**********************************************************************************************

	Get Modal helicity for (lx,ly,lz).
 Modal helicity = 0 for sincos basis.
 Equialent Fourier modes either purely real or purely imag.

***********************************************************************************************/

inline Real SSS_PENCIL::Get_Modal_helicity
(
	int lx, int ly, int lz, 
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az
)
{
	return 0;	
}

/**********************************************************************************************

	Compute Modal Vorticity
	Helicity is zero, so we are turning of the vorticity function at present.

***********************************************************************************************/


inline void SSS_PENCIL::Compute_Modal_vorticity
(
	int lx, int ly, int lz, 
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
	TinyVector<Complex,3> &vorticity
)
{
	/*
	TinyVector<Real,3> K;
	TinyVector<Real,3> V;
    
	
	V = Get_local_spectral_field(Ax, lx, ly, lz), Get_local_spectral_field(Ay, lx, ly, lz), Get_local_spectral_field(Az, lx, ly, lz); 
    
	// Actually it is vect(A)/I, but it cancels with I K x vect(A)/I
	
	Wavenumber(lx, ly, lz, K);
	
	vorticity(0) = complex<Real>(1,0)* (K(1)*V(2) - K(2)*V(1));
	vorticity(1) = complex<Real>(1,0)* (K(2)*V(0) - K(0)*V(2));
	vorticity(2) = complex<Real>(1,0)* (K(0)*V(1) - K(1)*V(0));
	 */
	
	vorticity = Complex(0,0) , Complex(0,0), Complex(0,0);
}



inline void SSS_PENCIL::Compute_Modal_vorticity_y_component
(
    int lx, int ly, int lz, 
    Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
    Complex &vort_y
)
{
	/*
	TinyVector<Real,3> K;
	TinyVector<Real,3> V;
	
	V = Get_local_spectral_field(Ax, lx, ly, lz), Get_local_spectral_field(Ay, lx, ly, lz), Get_local_spectral_field(Az, lx, ly, lz); 
	// Actually it is vect(A)/I, but it cancels with I K x vect(A)/I
	
	Wavenumber(lx, ly, lz, K);
	
	vort_y = complex<Real>(1,0)* (K(2)*V(0) - K(0)*V(2));
	 */
	
	vort_y = 0;
}



//*********************************************************************************************

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  \f$ K_{||} = K_1 \f$.		
inline Real SSS_PENCIL::AnisKpll(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_kx(lx)*kfactor[1]); 
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_ky(ly) * kfactor[2]);
	
	else if (global.field.anisotropy_dirn == 3)
		return (Get_kz(lz) * kfactor[3]);
		
	else
		return 0;		// for -Wall
}

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  \f$ K_\perp =\sqrt{K_2^2 + K_3^2} \f$.			
inline Real SSS_PENCIL::AnisKperp(int lx, int ly, int lz)
{
	if (global.field.anisotropy_dirn == 1)
		return sqrt( pow2(Get_ky(ly) * kfactor[2]) + pow2(Get_kz(lz)*kfactor[3]) ); 
	
	else if (global.field.anisotropy_dirn == 2)
		return sqrt( pow2(Get_kx(lx)*kfactor[1]) + pow2(Get_kz(lz)*kfactor[3]) );
		
	else if (global.field.anisotropy_dirn == 3)
		return sqrt( pow2(Get_kx(lx)*kfactor[1]) + pow2(Get_ky(ly)*kfactor[2]) );
			
	else
		return 0;		// for -Wall
}

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  horizontal direction 1, \f$ K_{h1} = K_2 \f$.										
inline Real SSS_PENCIL::AnisKh1(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_ky(ly) * kfactor[2]); 
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_kz(lz) * kfactor[3]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (Get_kx(lx) * kfactor[1]);
			
	else
		return 0;		// for -Wall
}

/// 3D == Anisotropic axis along x1: for anisotropic energy spectrum and 
///			energy transfer calculations,  horizontal direction 2, \f$ K_{h2} = K_3 \f$.				
inline Real SSS_PENCIL::AnisKh2(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_kz(lz) * kfactor[3]);  
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_kx(lx) * kfactor[1]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (Get_ky(ly) * kfactor[2]);
			
	else
		return 0;		// for -Wall
}
			
/// Cylindrical: Anis_min_Kpll
inline Real SSS_PENCIL::Anis_min_Kpll() 
{ 
	return 0.0;
}
				
/// Cylindrical: Anis_max_Kpll
inline Real SSS_PENCIL::Anis_max_Kpll() 
{ 
	
	Real maxKpll = 0.0;
	
    if (global.field.anisotropy_dirn == 1)
        maxKpll = ((Nx-1) * kfactor[1]); 
    
    else if (global.field.anisotropy_dirn == 2)
        maxKpll = ((Ny-1) * kfactor[2]); 
        
    else if (global.field.anisotropy_dirn == 3)
        maxKpll = ((Nz-1) * kfactor[3]); 
	
	return maxKpll;
}

	
/// 3D Cylindrical: Anis_max_Krho_radius_inside the wavenumber box.
inline int SSS_PENCIL::Anis_max_Krho_radius_inside() 			
{
	Real Kmag = 0.0;
	
	if (global.program.alias_option == "ALIAS")
	{
		if (global.field.anisotropy_dirn == 1)
			Kmag = min( (Ny-1)*kfactor[2], (Nz-1)*kfactor[3] ); 
		
		else if (global.field.anisotropy_dirn == 2)
			Kmag = min( (Nx-1)*kfactor[1], (Nz-1)*kfactor[3] );
			
		else if (global.field.anisotropy_dirn == 3)
			Kmag = min( (Nx-1)*kfactor[1], (Ny-1)*kfactor[2] ); 
	}
		
	else
	{
		if (global.field.anisotropy_dirn == 1)
			Kmag = min( (2*Ny/3)*kfactor[2], (2*Nz/3)*kfactor[3] );
		
		else if (global.field.anisotropy_dirn == 2)
			Kmag = min( (2*Nx/3)*kfactor[1], (2*Nz/3)*kfactor[3] );
			
		else if (global.field.anisotropy_dirn == 3)
			Kmag = min( (2*Nx/3)*kfactor[1], (2*Ny/3)*kfactor[2] ); 
	}
	
	return ((int) Kmag);	
}

// Max polar angle
inline Real SSS_PENCIL::Get_max_polar_angle() 
{	
	
	return M_PI/2;
}			
	
//*********************************************************************************************

/*! \brief Returns the angle K vector makes with the anisotropic axis 
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  lx, ly, lz (3D)
 * \return \f$ \tan^{-1}(K_{\perp}/K_{||}) \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline Real SSS_PENCIL::AnisKvect_polar_angle(int lx, int ly, int lz)
{
	Real kkpll, kkperp;
	
	kkpll = AnisKpll(lx, ly, lz);
	kkperp = AnisKperp(lx, ly, lz);
	
	return Get_polar_angle(kkperp, kkpll);
}


/*! \brief 3D: Returns the azimutal angle.
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  lx, ly, lz (3D)
 * \return \f$ \tan^{-1}(Ky}/Kx \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline Real SSS_PENCIL::AnisKvect_azimuthal_angle(int lx, int ly, int lz)
{
	
	Real kkh1 = AnisKh1(lx, ly, lz);
	Real kkh2 = AnisKh2(lx, ly, lz);
	
	return Get_azimuthal_angle(kkh1, kkh2);
}			
						

//=================================== End of inline functions ==========================//



