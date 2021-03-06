/* 
 * heliangliang, USTB, 2012.03.11, CopyRight Reserved
 * see "Reconstructing surfaces of particle-based fluids using anisotropic kernels (2010)"
*/

#include"surface_field.h"

void surface_field::surface_setup(vec3d min, vec3d max, double r_p, double r_h, double cell_size)
{
	//volume_p = 4.0/3.0*M_PI*r_p*r_p*r_p;
	volume_p = r_p*r_p*r_p;
	radius_h = r_h;

	// field grid setup
	field_grid_cell_size = cell_size;
	field_grid_min = min - r_h * 0.75;
	field_grid_cell_num = ((max-min+r_h*1.5)/cell_size).to<INT_T>() + 1;
	field_grid_divisor.set(1.0/cell_size, 1.0/cell_size, 1.0/cell_size);
	field_grid.resize(field_grid_cell_num.x*field_grid_cell_num.y*field_grid_cell_num.z);
	field_grid_num = field_grid.size();

	// partcile grid setup
	h_field = 2*r_h; // Equation (11), ri=k*hi
	grid_min = min - h_field * 1.25;
	grid_cell_num = ( (max-min+h_field*2.5)/h_field ).to<INT_T>() + 1;
	grid_divisor.set( 1.0/h_field, 1.0/h_field, 1.0/h_field );
	grid_first.resize( grid_cell_num.x*grid_cell_num.y*grid_cell_num.z );
	grid_first.assign(grid_first.size(), -1);

	p_next.resize(0);
	particle_num = p_next.size();

	// Equation (15)
	k_r = 4.0;
	k_s = 1.0/(0.16*h_field*h_field) * 1.6; // k_s=1.0/(0.16*h_field*h_field) such that ||k_s*Ci||=~=1, we set ||k_s*Ci||=~=1-2
	k_n = 0.95; // set k_n=~=1 such that the isolated particle don't disappear
	N_e = 25;
	// Equation (6)
	lambda = 1.0; // the closer lambda is to 1, the smoother the surface, Otherwise, the surface is faithful to the particles

} // void surface_field::surface_setup(vec3d min, vec3d max, double h, double cell_size)

const vector<RE_TYPE>& surface_field::surface_isotro(const vector<vec3d>& position)
{
	vec3<INT_T> p_n_start, p_n_end, p_n;
	vec3d p_p_start, p_p;
	INT_T i, k, num=position.size();
	double r, h_square=radius_h*radius_h;
	INT_T h_cell_num=INT_T(radius_h/field_grid_cell_size);
	
	field_grid.assign(field_grid.size(),RE_TYPE(0.0));
	for( i=0; i<num; ++i ) {
		p_n_end = p_n_start = ((position[i]-field_grid_min)*field_grid_divisor).to<INT_T>();
		p_n_start -= h_cell_num; p_n_end += h_cell_num+1;
		for( k=0; k<3; ++k ) {
			if(p_n_start[(int)k]<0) p_n_start[(int)k]=0;
			if(p_n_end[(int)k]>field_grid_cell_num[(int)k]-1) p_n_end[(int)k]=field_grid_cell_num[(int)k]-1;
		}
		p_p_start = field_grid_min + (p_n_start.to<double>()*field_grid_cell_size);
		for( p_n.z=p_n_start.z,p_p.z=p_p_start.z; p_n.z<=p_n_end.z; ++p_n.z,p_p.z+=field_grid_cell_size ) {
			for( p_n.y=p_n_start.y,p_p.y=p_p_start.y; p_n.y<=p_n_end.y; ++p_n.y,p_p.y+=field_grid_cell_size ) {
				k = (p_n.z * field_grid_cell_num.y + p_n.y) * field_grid_cell_num.x + p_n_start.x;
				for( p_n.x=p_n_start.x,p_p.x=p_p_start.x; p_n.x<=p_n_end.x; ++p_n.x,p_p.x+=field_grid_cell_size ) {
					r = (position[i]-p_p).len_square();
					if( r<=h_square ) {
						r = sqrt(r);
						field_grid[k] += RE_TYPE(ker_spline(r));
					}
					++ k;
				}
			}
		}
	} // for( i=0; i<num; ++i )
	return field_grid;
} // vector<RE_TYPE>& surface_field::surface_isotro(const vector<vec3d>& position)

#define M_COPY_Array2D(a, b) \
	a[0][0]=b[0][0], a[0][1]=b[0][1], a[0][2]=b[0][2], \
	a[1][0]=b[1][0], a[1][1]=b[1][1], a[1][2]=b[1][2], \
	a[2][0]=b[2][0], a[2][1]=b[2][1], a[2][2]=b[2][2]
#define MAX_NEIGHBOR 1024

const vector<RE_TYPE>& surface_field::surface_aniso(const vector<vec3d>& position)
{
	INT_T i, j, k, cell_neighbor[27], num=position.size();
	vec3d sum_x, xw, x_bar; double sum_w;
	INT_T neighbor_num, neighbor[MAX_NEIGHBOR]; double neighbor_w[MAX_NEIGHBOR];
	double weight, h_field_square=h_field*h_field;
	mat33d G_i, C_i, U, S, V;
	Array2D<double> tnt_a(3,3), tnt_u(3,3), tnt_v(3,3);
	Array1D<double> tnt_s(3);
	vec3<INT_T> p_n_start, p_n_end, p_n;
	vec3d p_p_start, p_p;
	double r, h_square=radius_h*radius_h;
	INT_T h_cell_num=INT_T(radius_h/field_grid_cell_size);

	INT_T out_factor = num/100;

	grid_insert(position);
	field_grid.assign(field_grid.size(),RE_TYPE(0.0));
	for( i=0; i<num; ++i ) {
		if((i+1)%out_factor==0){
			std::cout << "--  " << (i+1)/out_factor << "%\n";
		}
		grid_neighbors(position[i], cell_neighbor);
		sum_x.set(0.0,0.0,0.0); sum_w = 0.0; neighbor_num = 0;
		for( k=0; k<27; ++k ) { // Equation(10)
			for( j=grid_first[cell_neighbor[k]]; j>=0; j=p_next[j] ) {
				if( i==j ) {
					continue;
				}
				weight = (position[i]-position[j]).len_square();
				if( weight<=h_field_square ) {
					weight = weight_wij(sqrt(weight));
					sum_x += position[j] * weight;
					sum_w += weight;
					if( neighbor_num<MAX_NEIGHBOR ) {
						neighbor[neighbor_num] = j;
						neighbor_w[neighbor_num++] = weight;
					}
				}
			}
		} // for( k=0; k<27; ++k )
		if( neighbor_num>=N_e ) {
			xw = sum_x / sum_w; // Equation(10)
			//xw = position[i];
			x_bar = position[i]*(1-lambda) + xw*lambda; // Equation (6)
			C_i = mat33d::O;
			for( k=0; k<neighbor_num; ++k ) {
				j = neighbor[k]; weight = neighbor_w[k];
				//C_i_h += (position[j]-xw).pro_trans_self() * weight;
				C_i += (position[j]-x_bar).pro_trans_self() * weight; // change from Equation (6)
			}
			C_i /= sum_w; // Equation (9)
			M_COPY_Array2D(tnt_a, C_i);
			SVD<double> svd(tnt_a); // Equation (12)
			svd.getSingularValues(tnt_s);
			tnt_s[1] = std::max(tnt_s[1], tnt_s[0]/k_r);
			tnt_s[2] = std::max(tnt_s[2], tnt_s[0]/k_r); // Equation (15)
			tnt_s[0] *= k_s; tnt_s[1] *= k_s; tnt_s[2] *= k_s;
			//double a=tnt_s[0], b=tnt_s[1], c=tnt_s[2];
			svd.getU(tnt_u); svd.getV(tnt_v);
			M_COPY_Array2D(U, tnt_u); M_COPY_Array2D(V, tnt_v);
			S = mat33d::O; S[0][0] = 1.0/tnt_s[0]; S[1][1] = 1.0/tnt_s[1]; S[2][2] = 1.0/tnt_s[2];
			//S /= radius_h; // Equation (16)
			G_i = U.product(S).product(V.transpose());
		} else { // if( neighbor_num>=N_e )
			x_bar = position[i];
			G_i = mat33d::E*( k_n*(3-2*neighbor_num/double(N_e)) );
			tnt_s[0] = tnt_s[1] = tnt_s[2] = 1.0/k_n;
		} // if( neighbor_num>=N_e )

		p_n_end = p_n_start = ((x_bar-field_grid_min)*field_grid_divisor).to<INT_T>();
		p_n_start -= INT_T(h_cell_num*tnt_s[0]); p_n_end += INT_T(h_cell_num*tnt_s[0])+1;
		for( k=0; k<3; ++k ) {
			if(p_n_start[(int)k]<0) p_n_start[(int)k]=0;
			if(p_n_end[(int)k]>field_grid_cell_num[(int)k]-1) p_n_end[(int)k]=field_grid_cell_num[(int)k]-1;
		}
		p_p_start = field_grid_min + p_n_start.to<double>()*field_grid_cell_size;
		for( p_n.z=p_n_start.z,p_p.z=p_p_start.z; p_n.z<=p_n_end.z; ++p_n.z,p_p.z+=field_grid_cell_size ) {
			for( p_n.y=p_n_start.y,p_p.y=p_p_start.y; p_n.y<=p_n_end.y; ++p_n.y,p_p.y+=field_grid_cell_size ) {
				k = (p_n.z * field_grid_cell_num.y + p_n.y) * field_grid_cell_num.x + p_n_start.x;
				for( p_n.x=p_n_start.x,p_p.x=p_p_start.x; p_n.x<=p_n_end.x; ++p_n.x,p_p.x+=field_grid_cell_size ) {
					r = (p_p-x_bar).pro_matrix_left(G_i).len_square();
					if( r<=h_square ) {
						r = sqrt(r);
						field_grid[k] += RE_TYPE(G_i.determinant()*ker_spline(r));
						//field_grid[k] += RE_TYPE(ker_spline(r));
					}
					++ k;
				}
			}
		}
	} // for( i=0; i<num; ++i )
	return field_grid;
} // vector<RE_TYPE>& surface_field::surface_aniso(const vector<vec3d>& position)

void surface_field::grid_insert(const vector<vec3d>& position)
{
	if( p_next.size()!=position.size() ) {
		p_next.resize(position.size());
		particle_num = p_next.size();
	}

	grid_first.assign(grid_first.size(),-1);
	for( INT_T i=0,ngrid; i<particle_num; ++i ) {
		ngrid = ( INT_T((position[i].z-grid_min.z)*grid_divisor.z)*grid_cell_num.y
				+ INT_T((position[i].y-grid_min.y)*grid_divisor.y) ) * grid_cell_num.x
				+ INT_T((position[i].x-grid_min.x)*grid_divisor.x);
		p_next[i] = grid_first[ngrid];
		grid_first[ngrid] = i;
	}
}

void surface_field::grid_neighbors(const vec3d& p, INT_T* r)
{
	r[0] = ( (INT_T((p.z-grid_min.z)*grid_divisor.z)-1)*grid_cell_num.y
		   + (INT_T((p.y-grid_min.y)*grid_divisor.y)-1) ) * grid_cell_num.x
		   + (INT_T((p.x-grid_min.x)*grid_divisor.x)-1);	// 0 0 0
	r[1] = r[0] + 1;				// 0 0 1
	r[2] = r[0] + 2;				// 0 0 2
	r[3] = r[0] + grid_cell_num.x;	// 0 1 0
	r[4] = r[3] + 1;				// 0 1 1
	r[5] = r[3] + 2;				// 0 1 2
	r[6] = r[3] + grid_cell_num.x;	// 0 2 0
	r[7] = r[6] + 1;				// 0 2 1
	r[8] = r[6] + 2;				// 0 2 2
	r[9] = r[0] + grid_cell_num.x*grid_cell_num.y;	// 1 0 0
	r[10] = r[9] + 1;				// 1 0 1
	r[11] = r[9] + 2;				// 1 0 2
	r[12] = r[9] + grid_cell_num.x;	// 1 1 0
	r[13] = r[12] + 1;				// 1 1 1
	r[14] = r[12] + 2;				// 1 1 2
	r[15] = r[12] + grid_cell_num.x;// 1 2 0
	r[16] = r[15] + 1;				// 1 2 1
	r[17] = r[15] + 2;				// 1 2 2
	r[18] = r[9] + grid_cell_num.x*grid_cell_num.y;	// 2 0 0
	r[19] = r[18] + 1;				// 1 0 1
	r[20] = r[18] + 2;				// 1 0 2
	r[21] = r[18] + grid_cell_num.x;// 1 1 0
	r[22] = r[21] + 1;				// 1 1 1
	r[23] = r[21] + 2;				// 1 1 2
	r[24] = r[21] + grid_cell_num.x;// 1 2 0
	r[25] = r[24] + 1;				// 1 2 1
	r[26] = r[24] + 2;				// 1 2 2
}

/* B-cubic spline kernel W_spline
 * "Smoothed Particles: A new paradigm for animating highly deformable bodies"
 *
 *            / 1 - (6/h^2) r^2 + (6/h^3) r^3  0 - h/2
 * 8/(pi*h^3) | (2/h^3) (h-r)^3                h/2 - h
 *            \ 0                              r<0 or r>h
*/
inline double surface_field::ker_spline(double r)
{
	static double param1 = 8.0 / ( M_PI*pow(radius_h,3.0) ) * volume_p;
	static double param2 = -48.0 / ( M_PI*pow(radius_h,5.0) ) * volume_p;
	static double param3 = 48.0 / ( M_PI*pow(radius_h,6.0) ) * volume_p;
	static double param4 = 16.0 / ( M_PI*pow(radius_h,6.0) ) * volume_p;
	static double h_half = radius_h*0.5;

	/*if(r<0.0 || r>h) {
		return 0.0;
	}*/

	return (r<h_half) ? ( param1 + param2* r*r + param3* r*r*r ) :
						( param4 * (radius_h-r)*(radius_h-r)*(radius_h-r) );
} // double surface_field::ker_spline(double r)

// Equation (11)
inline double surface_field::weight_wij(double r)
{
	static double param = 1.0 / (h_field*h_field*h_field);
	//static double param = 1.0 / (radius_h*2.0);

	return ( 1-r*r*r*param );
	//return (1-r*param)*(1-r*param)*(1-r*param);
}

