// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file container.hh
 * \brief Header file for the container_base template and related classes. */

#ifndef VOROPP_CONTAINER_HH
#define VOROPP_CONTAINER_HH

#include "config.hh"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class voropp_loop;
class radius_poly;
class wall;

/** \brief A class representing the whole simulation region.
 *
 * The container class represents the whole simulation region. The
 * container constructor sets up the geometry and periodicity, and divides
 * the geometry into rectangular grid of blocks, each of which handles the
 * particles in a particular area. Routines exist for putting in particles,
 * importing particles from standard input, and carrying out Voronoi
 * calculations. */
template<class r_option>
class container_base {
	public:
		container_base(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi);
		~container_base();
		void draw_particles(const char *filename);
		void draw_particles_rgl(const char *filename,std::vector<int>*res_id);
		void draw_particles();
		void draw_particles(ostream &os);
		void draw_particles_rgl(ostream &os,std::vector<int>*res_id);
		void draw_particles_pov(const char *filename);
		void draw_particles_pov();
		void draw_particles_pov(ostream &os);
		void import(istream &is);
		inline void import();
		inline void import(const char *filename);
		void region_count();
		void clear();
		void draw_cells_gnuplot(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void draw_cells_gnuplot(const char *filename);
		void draw_cells_pov(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		void draw_cells_rgl(const char *filename,std::vector<int>*res_id,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void draw_cells_pov(const char *filename);
		inline void draw_cells_rgl(const char *filename,std::vector<int>*res_id);
		void store_cell_volumes(fpoint *bb);
		fpoint packing_fraction(fpoint *bb,fpoint cx,fpoint cy,fpoint cz,fpoint r);
		fpoint packing_fraction(fpoint *bb,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		fpoint sum_cell_volumes();
		void compute_all_cells();
		void print_all(ostream &os);
		void print_all();
		void print_all(const char *filename);
		void print_all_neighbor(ostream &os);
		void print_all_neighbor();
		void print_all_neighbor(const char *filename);
		void print_all_custom(const char *format,ostream &os);
		void print_all_custom(const char *format);
		void get_res_volume(
			vector<int>* res_types,
			vector<int>* res_id,
			vector<string>* at_type,
			
			vector<fpoint> 	*volume,
			vector<int> 	*unpacked,
			vector<fpoint> 	*surface,
			vector<fpoint>	*surfA,
			vector<fpoint>	*surfG,
			vector<fpoint>	*surfC,
			vector<fpoint>	*surfU,
			vector<int> 	*orderA,
			vector<int> 	*orderG,
			vector<int> 	*orderC,
			vector<int> 	*orderU,

			vector<fpoint> 	*bvolume,
			vector<int> 	*bunpacked,
			vector<fpoint> 	*bsurface,
			vector<fpoint>	*bsurfA,
			vector<fpoint>	*bsurfG,
			vector<fpoint>	*bsurfC,
			vector<fpoint>	*bsurfU,
			vector<int> 	*borderA,
			vector<int> 	*borderG,
			vector<int> 	*borderC,
			vector<int> 	*borderU,
	
vector<fpoint> *vatom0,
	vector<fpoint> *vatom1,
	vector<fpoint> *vatom2,
	vector<fpoint> *vatom3,
	vector<fpoint> *vatom4,
	vector<fpoint> *vatom5,
	vector<fpoint> *vatom6,
	vector<fpoint> *vatom7,
	vector<fpoint> *vatom8,
	vector<fpoint> *vatom9,
	vector<fpoint> *vatom10,
	vector<fpoint> *vatom11,
	vector<fpoint> *vatom12,
	vector<fpoint> *vatom13,
	vector<fpoint> *vatom14,
	vector<fpoint> *vatom15,
	vector<fpoint> *vatom16,
	vector<fpoint> *vatom17,
	vector<fpoint> *vatom18,
	vector<fpoint> *vatom19,
	vector<fpoint> *vatom20,
	vector<fpoint> *vatom21,
	vector<fpoint> *vatom22,
	vector<fpoint> *vatom23,
	vector<fpoint> *vatom24,
	vector<fpoint> *vatom25,
	vector<fpoint> *vatom26,

	vector<int> *atom0,
	vector<int> *atom1,
	vector<int> *atom2,
	vector<int> *atom3,
	vector<int> *atom4,
	vector<int> *atom5,
	vector<int> *atom6,
	vector<int> *atom7,
	vector<int> *atom8,
	vector<int> *atom9,
	vector<int> *atom10,
	vector<int> *atom11,
	vector<int> *atom12,
	vector<int> *atom13,
	vector<int> *atom14,
	vector<int> *atom15,
	vector<int> *atom16,
	vector<int> *atom17,
	vector<int> *atom18,
	vector<int> *atom19,
	vector<int> *atom20,
	vector<int> *atom21,
	vector<int> *atom22,
	vector<int> *atom23,
	vector<int> *atom24,
	vector<int> *atom25,
	vector<int> *atom26,



			vector<int>	*parent,
			vector<int>	*child,
			vector<int>	*order
			);
		void print_all_custom(const char *format,const char *filename);
		template<class n_option>
		inline bool compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s);
		template<class n_option>
		bool compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z);
		template<class n_option>
		inline bool compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s);
		template<class n_option>
		bool compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z,fpoint r);
		void add_wall(wall &w);
		bool point_inside(fpoint x,fpoint y,fpoint z);
		bool point_inside_walls(fpoint x,fpoint y,fpoint z);
	protected:
		/** The minimum x coordinate of the container. */
		const fpoint ax;
		/** The maximum x coordinate of the container. */
		const fpoint bx;
		/** The minimum y coordinate of the container. */
		const fpoint ay;
		/** The maximum y coordinate of the container. */
		const fpoint by;
		/** The minimum z coordinate of the container. */
		const fpoint az;
		/** The maximum z coordinate of the container. */
		const fpoint bz;
		/** The inverse box length in the x direction, set to
		 * nx/(bx-ax). */
		const fpoint xsp;
		/** The inverse box length in the y direction, set to
		 * ny/(by-ay). */
		const fpoint ysp;
		/** The inverse box length in the z direction, set to
		 * nz/(bz-az). */
		const fpoint zsp;
		/** The number of boxes in the x direction. */
		const int nx;
		/** The number of boxes in the y direction. */
		const int ny;
		/** The number of boxes in the z direction. */
		const int nz;
		/** A constant, set to the value of nx multiplied by ny, which
		 * is used in the routines which step through boxes in
		 * sequence. */
		const int nxy;
		/** A constant, set to the value of nx*ny*nz, which is used in
		 * the routines which step through boxes in sequence. */
		const int nxyz;
		/** The number of boxes in the x direction for the searching mask. */
		const int hx;
		/** The number of boxes in the y direction for the searching mask. */
		const int hy;
		/** The number of boxes in the z direction for the searching mask. */
		const int hz;
		/** A constant, set to the value of hx multiplied by hy, which
		 * is used in the routines which step through mask boxes in
		 * sequence. */
		const int hxy;
		/** A constant, set to the value of hx*hy*hz, which is used in
		 * the routines which step through mask boxes in sequence. */
		const int hxyz;
		/** A boolean value that determines if the x coordinate in
		 * periodic or not. */
		const bool xperiodic;
		/** A boolean value that determines if the y coordinate in
		 * periodic or not. */
		const bool yperiodic;
		/** A boolean value that determines if the z coordinate in
		 * periodic or not. */
		const bool zperiodic;
		/** This sets the current value being used to mark tested blocks
		 * in the mask. */
		unsigned int mv;
		/** The current number of wall objects, initially set to zero. */
		int wall_number;
		/** The current amount of memory allocated for walls. */
		int current_wall_size;
		/** This object contains all the functions for handling how
		 * the particle radii should be treated. If the template is
		 * instantiated with the radius_mono class, then this object
		 * contains mostly blank routines that do nothing to the
		 * cell computation, to compute the basic Voronoi diagram.
		 * If the template is instantiated with the radius_poly calls,
		 * then this object provides routines for modifying the
		 * Voronoi cell computation in order to create the radical
		 * Voronoi tessellation. */
		r_option radius;
		/** The amount of memory in the array structure for each
		 * particle. This is set to 3 when the basic class is
		 * initialized, so that the array holds (x,y,z) positions. If
		 * the container class is initialized as part of the derived
		 * class container_poly, then this is set to 4, to also hold
		 * the particle radii. */
		int sz;
		/** The position of the first element on the search list to be
		 * considered. */
		int s_start;
		/** The position of the last element on the search list to be
		 * considered. */
		int s_end;
		/** The current size of the search list. */
		int s_size;		
		/** This array holds the number of particles within each
		 * computational box of the container. */
		int *co;
		/** This array holds the maximum amount of particle memory for
		 * each computational box of the container. If the number of
		 * particles in a particular box ever approaches this limit,
		 * more is allocated using the add_particle_memory() function.
		 */
		int *mem;
		/** This array is used during the cell computation to determine
		 * which blocks have been considered. */
		unsigned int *mask;
		/** This array is used to store the list of blocks to test during
		 * the Voronoi cell computation. */
		int *sl;		
		/** An array to hold the minimum distances associated with the
		 * worklists. This array is initialized during container
		 * construction, by the initialize_radii() routine. */
		fpoint *mrad;
		/** This array holds pointers to any wall objects that have
		 * been added to the container. */
		wall **walls;
		/** This array holds the numerical IDs of each particle in each
		 * computational box. */
		int **id;
		/** A two dimensional array holding particle positions. For the
		 * derived container_poly class, this also holds particle
		 * radii. */
		fpoint **p;

		template<class n_option>
		inline void print_all_internal(voronoicell_base<n_option> &c,ostream &os);
		template<class n_option>
		void print_all_custom_internal(voronoicell_base<n_option> &c,const char *format,ostream &os);
		template<class n_option>
		inline bool initialize_voronoicell(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
		void add_particle_memory(int i);
		void add_list_memory();
	private:
#include "worklist.hh"
		template<class n_option>
		inline bool corner_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint zl,fpoint xh,fpoint yh,fpoint zh);
		template<class n_option>
		inline bool edge_x_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint zl,fpoint x1,fpoint yh,fpoint zh);
		template<class n_option>
		inline bool edge_y_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint zl,fpoint xh,fpoint y1,fpoint zh);
		template<class n_option>
		inline bool edge_z_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint z0,fpoint xh,fpoint yh,fpoint z1);
		template<class n_option>
		inline bool face_x_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint z0,fpoint y1,fpoint z1);
		template<class n_option>
		inline bool face_y_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint z0,fpoint x1,fpoint z1);
		template<class n_option>
		inline bool face_z_test(voronoicell_base<n_option> &c,fpoint x0,fpoint y0,fpoint zl,fpoint x1,fpoint y1);
		inline void initialize_radii();
		inline void compute_minimum(fpoint &minr,fpoint &xlo,fpoint &xhi,fpoint &ylo,fpoint &yhi,fpoint &zlo,fpoint &zhi,int ti,int tj,int tk);
		inline bool compute_min_max_radius(int di,int dj,int dk,fpoint fx,fpoint fy,fpoint fz,fpoint gx,fpoint gy,fpoint gz,fpoint& crs,fpoint mrs);
		friend class voropp_loop;
		friend class radius_poly;
};

/** \brief A class encapsulating all routines specifically needed in
 * the standard Voronoi tessellation.
 *
 * This class encapsulates all the routines that are required for carrying out
 * a standard Voronoi tessellation that would be appropriate for a monodisperse
 * system. When the container class is instantiated using this class, all
 * information about particle radii is switched off. Since all these functions
 * are declared inline, there should be no loss of speed. */
class radius_mono {
	public:
		/** The number of floating point numbers allocated for each
		 * particle in the container, set to 3 for this case for the x,
		 * y, and z positions. */
		const int mem_size;
		/** This constructor sets a pointer back to the container class
		 * that created it, and initializes the mem_size constant to 3.
		 * \param[in] icc a pointer the container class that created
		 *                this class. */
		radius_mono(container_base<radius_mono> *icc) : mem_size(3), cc(icc) {};
		inline void import(istream &is);
		/** This is a blank placeholder function that does nothing. */
		inline void store_radius(int i,int j,fpoint r) {};
		/** This is a blank placeholder function that does nothing. */
		inline void clear_max() {};
		/** This is a blank placeholder function that does nothing. */
		inline void init(int s,int i) {};
		inline fpoint volume(int ijk,int s);
		inline fpoint cutoff(fpoint lrs);
		inline fpoint scale(fpoint rs,int t,int q);
		/** This is a blank placeholder function that does nothing. */
		inline void print(ostream &os,int ijk,int q,bool later=true) {};
		inline void rad(ostream &os,int l,int c);
	private:
		container_base<radius_mono> *cc;
};

/** \brief A class encapsulating all routines specifically needed in the
 * Voronoi radical tessellation.
 *
 * This class encapsulates all the routines that are required for carrying out
 * the radical Voronoi tessellation that is appropriate for polydisperse sphere.
 * When the container class is instantiated with this class, information about particle
 * radii is switched on. */
class radius_poly {
	public:
		/** The number of floating point numbers allocated for each
		 * particle in the container, set to 4 for this case for the x,
		 * y, and z positions, plus the radius. */
		const int mem_size;
		/** This constructor sets a pointer back to the container class
		 * that created it, and initializes the mem_size constant to 4.
		 */
		radius_poly(container_base<radius_poly> *icc) :
			mem_size(4), cc(icc), max_radius(0) {};
		inline void import(istream &is);
		inline void store_radius(int i,int j,fpoint r);
		inline void clear_max();
		inline void init(int ijk,int s);
		inline fpoint volume(int ijk,int s);
		inline fpoint cutoff(fpoint lrs);
		inline fpoint scale(fpoint rs,int t,int q);
		inline void print(ostream &os,int ijk,int q,bool later=true);
		inline void rad(ostream &os,int l,int c);
	private:
		container_base<radius_poly> *cc;
		fpoint max_radius,crad,mul;
};

/** \brief A class to handle loops on regions of the container handling
 * non-periodic and periodic boundary conditions.
 *
 * Many of the container routines require scanning over a rectangular sub-grid
 * of blocks, and the routines for handling this are stored in the voropp_loop
 * class. A voropp_loop class can first be initialized to either calculate the
 * subgrid which is within a distance r of a vector (vx,vy,vz), or a subgrid
 * corresponding to a rectangular box. The routine inc() can then be
 * successively called to step through all the blocks within this subgrid.
 */
class voropp_loop {
	public:
		template<class r_option>
		voropp_loop(container_base<r_option> *q);
		inline int init(fpoint vx,fpoint vy,fpoint vz,fpoint r,fpoint &px,fpoint &py,fpoint &pz);
		inline int init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax,fpoint &px,fpoint &py,fpoint &pz);
		inline int inc(fpoint &px,fpoint &py,fpoint &pz);
		/** The current block index in the x direction, referencing a
		 * real cell in the range 0 to nx-1. */
		int ip;
		/** The current block index in the y direction, referencing a
		 * real cell in the range 0 to ny-1. */
		int jp;
		/** The current block index in the z direction, referencing a
		 * real cell in the range 0 to nz-1. */
		int kp;
	private:
		int i,j,k,ai,bi,aj,bj,ak,bk,s;
		int aip,ajp,akp,inc1,inc2;
		inline int step_mod(int a,int b);
		inline int step_div(int a,int b);
		inline int step_int(fpoint a);
		fpoint apx,apy,apz;
		const fpoint sx,sy,sz,xsp,ysp,zsp,ax,ay,az;
		const int nx,ny,nz,nxy,nxyz;
		const bool xperiodic,yperiodic,zperiodic;
};

/** \brief Pure virtual class from which wall objects are derived.
 *
 * This is a pure virtual class for a generic wall object. A wall object
 * can be specified by deriving a new class from this and specifying the
 * functions.*/
class wall {
	public:
		virtual ~wall() {};
		/** A pure virtual function for testing whether a point is
		 * inside the wall object. */
		virtual bool point_inside(fpoint x,fpoint y,fpoint z) = 0;
		/** A pure virtual function for cutting a cell without
		 * neighbor-tracking with a wall. */
		virtual bool cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) = 0;
		/** A pure virtual function for cutting a cell with
		 * neighbor-tracking enabled with a wall. */
		virtual bool cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) = 0;
};

/** The basic container class. */
typedef container_base<radius_mono> container;

/** The polydisperse container class. */
typedef container_base<radius_poly> container_poly;
#endif
