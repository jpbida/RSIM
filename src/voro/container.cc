// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file container.cc
 * \brief Function implementations for the container_base template and related
 * classes. */

#include "cell.hh"
#include "container.hh"

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (xa,xb) the minimum and maximum x coordinates.
 * \param[in] (ya,yb) the minimum and maximum y coordinates.
 * \param[in] (za,zb) the minimum and maximum z coordinates.
 * \param[in] (xn,yn,zn) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] (xper,yper,zper) flags setting whether the container is periodic
 *                             in each coordinate direction.
 * \param[in] memi the initial memory allocation for each block. */
template<class r_option>
container_base<r_option>::container_base(fpoint xa,fpoint xb,fpoint ya,
		fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,
		bool xper,bool yper,bool zper,int memi)
	: ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
	xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),zsp(zn/(zb-za)),nx(xn),ny(yn),nz(zn),
	nxy(xn*yn),nxyz(xn*yn*zn),hx(xper?2*xn+1:xn),hy(yper?2*yn+1:yn),
	hz(zper?2*zn+1:zn),hxy(hx*hy),hxyz(hx*hy*hz),
	xperiodic(xper),yperiodic(yper),zperiodic(zper),
	mv(0),wall_number(0),current_wall_size(init_wall_size),radius(this),
	sz(radius.mem_size),s_size(3*(3+hxy+hz*(hx+hy))),
	co(new int[nxyz]),mem(new int[nxyz]),mask(new unsigned int[hxyz]),
	sl(new int[s_size]),mrad(new fpoint[hgridsq*seq_length]),
	walls(new wall*[init_wall_size]),id(new int*[nxyz]),p(new fpoint*[nxyz]) {
	int l;
	for(l=0;l<nxyz;l++) co[l]=0;
	for(l=0;l<nxyz;l++) mem[l]=memi;
	for(l=0;l<hxyz;l++) mask[l]=0;
	for(l=0;l<nxyz;l++) id[l]=new int[memi];
	for(l=0;l<nxyz;l++) p[l]=new fpoint[sz*memi];

	// Precompute the radius table used in the cell construction
	initialize_radii();
}

/** The container destructor frees the dynamically allocated memory. */
template<class r_option>
container_base<r_option>::~container_base() {
	int l;
	for(l=0;l<nxyz;l++) delete [] p[l];
	for(l=0;l<nxyz;l++) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] walls;
	delete [] mrad;
	delete [] sl;
	delete [] mask;
	delete [] mem;
	delete [] co;
}

/** Dumps all the particle positions and identifies to a file.
 * \param[in] os an output stream to write to. */
template<class r_option>
void container_base<r_option>::draw_particles(ostream &os) {
	int c,l,i;
	for(l=0;l<nxyz;l++) for(c=0;c<co[l];c++) {
		os << id[l][c];
		for(i=sz*c;i<sz*(c+1);i++) os << " " << p[l][i];
		os << "\n";
	}
}

template<class r_option>
void container_base<r_option>::draw_particles_rgl(ostream &os,std::vector<int> * res_id) {
	int c,l,i;
	for(l=0;l<nxyz;l++) for(c=0;c<co[l];c++) {
		os << res_id->at(id[l][c]);
		for(i=sz*c;i<sz*(c+1);i++) os << " " << p[l][i];
		os << "\n";
	}
}


/** An overloaded version of the draw_particles() routine, that just prints
 * to standard output. */
template<class r_option>
void container_base<r_option>::draw_particles() {
	draw_particles(cout);
}

/** An overloaded version of the draw_particles() routine, that outputs
 * the particle positions to a file.
 * \param[in] filename the file to write to. */
template<class r_option>
void container_base<r_option>::draw_particles(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_particles(os);
	os.close();
}

template<class r_option>
void container_base<r_option>::draw_particles_rgl(const char *filename,std::vector<int> *res_id) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_particles_rgl(os,res_id);
	os.close();
}

/** Dumps all the particle positions in POV-Ray format.
 * \param[in] os an output stream to write to. */
template<class r_option>
void container_base<r_option>::draw_particles_pov(ostream &os) {
	int l,c;
	for(l=0;l<nxyz;l++) for(c=0;c<co[l];c++) {
		os << "// id " << id[l][c] << "\n";
		os << "sphere{<" << p[l][sz*c] << "," << p[l][sz*c+1] << ","
		   << p[l][sz*c+2] << ">,";
		radius.rad(os,l,c);
		os << "}\n";
	}
}

/** An overloaded version of the draw_particles_pov() routine, that just prints
 * to standard output. */
template<class r_option>
void container_base<r_option>::draw_particles_pov() {
	draw_particles_pov(cout);
}

/** An overloaded version of the draw_particles_pov() routine, that outputs
 * the particle positions to a file.
 * \param[in] filename the file to write to. */
template<class r_option>
void container_base<r_option>::draw_particles_pov(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_particles_pov(os);
	os.close();
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
template<class r_option>
void container_base<r_option>::put(int n,fpoint x,fpoint y,fpoint z) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
			if(co[i]==mem[i]) add_particle_memory(i);
			p[i][sz*co[i]]=x;p[i][sz*co[i]+1]=y;p[i][sz*co[i]+2]=z;
			radius.store_radius(i,co[i],0.5);
			id[i][co[i]++]=n;
		}
	}
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle.*/
template<class r_option>
void container_base<r_option>::put(int n,fpoint x,fpoint y,fpoint z,fpoint r) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
			if(co[i]==mem[i]) add_particle_memory(i);
			p[i][sz*co[i]]=x;p[i][sz*co[i]+1]=y;p[i][sz*co[i]+2]=z;
			radius.store_radius(i,co[i],r);
			id[i][co[i]++]=n;
		}
	}
}

/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
template<class r_option>
void container_base<r_option>::add_particle_memory(int i) {
	int *idp;fpoint *pp;
	int l,nmem=2*mem[i];
#if VOROPP_VERBOSE >=3
	cerr << "Particle memory in region " << i << " scaled up to " << nmem << endl;
#endif
	if(nmem>max_particle_memory)
		voropp_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new fpoint[sz*nmem];
	for(l=0;l<sz*co[i];l++) pp[l]=p[i][l];
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
}

/** Add list memory. */
template<class r_option>
inline void container_base<r_option>::add_list_memory() {
	int i,j=0,*ps;
	ps=new int[s_size*2];
#if VOROPP_VERBOSE >=2
	cerr << "List memory scaled up to " << s_size*2 << endl;
#endif
	if(s_start<=s_end) {
		for(i=s_start;i<s_end;i++) ps[j++]=sl[i];
	} else {
		for(i=s_start;i<s_size;i++) ps[j++]=sl[i];
		for(i=0;i<s_end;i++) ps[j++]=sl[i];
	}
	s_size*=2;
	s_start=0;s_end=j;
	delete [] sl;sl=ps;
}

/** Import a list of particles from standard input.
 * \param[in] is a standard input stream to read from. */
template<class r_option>
void container_base<r_option>::import(istream &is) {
	radius.import(is);
}

/** An overloaded version of the import routine, that reads the standard input.
 */
template<class r_option>
inline void container_base<r_option>::import() {
	import(cin);
}

/** An overloaded version of the import routine, that reads in particles from
 * a particular file.
 * \param[in] filename the name of the file to read from. */
template<class r_option>
inline void container_base<r_option>::import(const char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	if(is.fail()) voropp_fatal_error("Unable to open file for import",VOROPP_FILE_ERROR);
	import(is);
	is.close();
}

/** Outputs the number of particles within each region. */
template<class r_option>
void container_base<r_option>::region_count() {
	int i,j,k,ijk=0;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++)
		cout << "Region (" << i << "," << j << "," << k << "): " << co[ijk++] << " particles" << endl;
}

/** Clears a container of particles. */
template<class r_option>
void container_base<r_option>::clear() {
	for(int ijk=0;ijk<nxyz;ijk++) co[ijk]=0;
	radius.clear_max();
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in gnuplot format.
 * \param[in] filename the name of the file to write to.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box. */
template<class r_option>
void container_base<r_option>::draw_cells_gnuplot(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px,py,pz;
	voropp_loop l1(this);
	int q,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px;y=p[s][sz*q+1]+py;z=p[s][sz*q+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				if(compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)) c.draw_gnuplot(os,x,y,z);
			}
		}
	} while((s=l1.inc(px,py,pz))!=-1);
	os.close();
}

/** An overloaded version of draw_cells_gnuplot() that computes the Voronoi
 * cells for the entire simulation region and saves the output in gnuplot
 * format.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
void container_base<r_option>::draw_cells_gnuplot(const char *filename) {
	draw_cells_gnuplot(filename,ax,bx,ay,by,az,bz);
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in POV-Ray format.
 * \param[in] filename the name of the file to write to.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box. */
template<class r_option>
void container_base<r_option>::draw_cells_pov(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px,py,pz;
	voropp_loop l1(this);
	int q,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(q=0;q<co[s];q++) {
			os << "// cell " << id[s][q] << "\n";
			x=p[s][sz*q]+px;y=p[s][sz*q+1]+py;z=p[s][sz*q+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				if(compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)) c.draw_pov(os,x,y,z);
			}
		}
	} while((s=l1.inc(px,py,pz))!=-1);
	os.close();
}

template<class r_option>
void container_base<r_option>::draw_cells_rgl(const char *filename,std::vector<int> *res_id,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px,py,pz;
	voropp_loop l1(this);
	int q,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px;y=p[s][sz*q+1]+py;z=p[s][sz*q+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				if(compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)){
if(c.surface_area() < 75){
c.draw_rgl(os,res_id->at(id[s][q]),x,y,z);}
}
			}
		}
	} while((s=l1.inc(px,py,pz))!=-1);
	os.close();
}

/** An overloaded version of draw_cells_pov() that computes the Voronoi
 * cells for the entire simulation region and saves the output in POV-Ray
 * format.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
void container_base<r_option>::draw_cells_pov(const char *filename) {
	draw_cells_pov(filename,ax,bx,ay,by,az,bz);
}

template<class r_option>
void container_base<r_option>::draw_cells_rgl(const char *filename,std::vector<int> *res_id) {
	draw_cells_rgl(filename,res_id,ax,bx,ay,by,az,bz);
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
template<class r_option>
void container_base<r_option>::compute_all_cells() {
	voronoicell c;
	int i,j,k,ijk=0,q;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(q=0;q<co[ijk];q++) compute_cell(c,i,j,k,ijk,q);
	}
}

/** Computes the Voronoi volumes for all the particles, and stores the results
 * according to the particle ID numbers in a floating point array that the user
 * has supplied. No bounds checking on the array is performed, so it is up to
 * the user to ensure that the array is large enough to store the computed
 * numbers.
 * \param[in] bb a pointer to an array to store the volumes. The volume of the
 *               particle with ID number n will be stored at bb[n]. */
template<class r_option>
void container_base<r_option>::store_cell_volumes(fpoint *bb) {
	voronoicell c;
	int i,j,k,ijk=0,q;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(q=0;q<co[ijk];q++) bb[id[ijk][q]]=compute_cell(c,i,j,k,ijk,q)?c.volume():0;
	}
}

/** Computes the local packing fraction at a point, by summing the volumes
 * of all particles within a test sphere, and dividing by the sum of their
 * Voronoi volumes that were previously computed using the store_cell_volumes()
 * function.
 * \param[in] bb a pointer to an array holding the Voronoi volumes of the
 *               particles.
 * \param[in] (cx,cy,cz) the center of the test sphere.
 * \param[in] r the radius of the test sphere. */
template<class r_option>
fpoint container_base<r_option>::packing_fraction(fpoint *bb,fpoint cx,fpoint cy,fpoint cz,fpoint r) {
	voropp_loop l1(this);
	fpoint px,py,pz,x,y,z,rsq=r*r,pvol=0,vvol=0;
	int q,s;
	s=l1.init(cx,cy,cz,r,px,py,pz);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px-cx;
			y=p[s][sz*q+1]+py-cy;
			z=p[s][sz*q+2]+pz-cz;
			if(x*x+y*y+z*z<rsq) {
				pvol+=radius.volume(s,q);
				vvol+=bb[id[s][q]];
			}
		}
	} while((s=l1.inc(px,py,pz))!=-1);
	return vvol>tolerance?pvol/vvol*4.1887902047863909846168578443726:0;
}

/** Computes the local packing fraction at a point, by summing the volumes of
 * all particles within test box, and dividing by the sum of their Voronoi
 * volumes that were previous computed using the store_cell_volumes() function.
 * \param[in] bb an array holding the Voronoi volumes of the particles.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box. */
template<class r_option>
fpoint container_base<r_option>::packing_fraction(fpoint *bb,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	voropp_loop l1(this);
	fpoint x,y,z,px,py,pz,pvol=0,vvol=0;
	int q,s;
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px;
			y=p[s][sz*q+1]+py;
			z=p[s][sz*q+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				pvol+=radius.volume(s,q);
				vvol+=bb[id[s][q]];
			}
		}
	} while((s=l1.inc(px,py,pz))!=-1);
	return vvol>tolerance?pvol/vvol*4.1887902047863909846168578443726:0;
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
template<class r_option>
fpoint container_base<r_option>::sum_cell_volumes() {
	voronoicell c;
	int i,j,k,ijk=0,q;fpoint vol=0;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(q=0;q<co[ijk];q++) if (compute_cell(c,i,j,k,ijk,q)) vol+=c.volume();
	}
	return vol;
}

/** Computes the Voronoi cells for all particles in the container, and for each
 * cell, outputs a line containing custom information about the cell structure.
 * The output format is specified using an input string with control sequences
 * similar to the standard C printf() routine. Full information about the
 * control sequences is available at http://math.lbl.gov/voro++/doc/custom.html
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] os an open output stream to write to. */
template<class r_option>
void container_base<r_option>::print_all_custom(const char *format,ostream &os) {
	int fp=0;

	// Check to see if the sequence "%n" appears in the format sequence
	while(format[fp]!=0) {
		if(format[fp]=='%') {
			fp++;
			if(format[fp]=='n') {

				// If a "%n" is detected, then we're going to
				// need neighbor information during the custom
				// output, so use the voronoicell_neighbor
				// class
				voronoicell_neighbor c;
				print_all_custom_internal(c,format,os);
				return;
			} else if(format[fp]==0) break;
		}
		fp++;
	}

	// No "%n" was detected, so we can just use the regular voronoicell
	// class without computing neighbor information.
	voronoicell c;
	print_all_custom_internal(c,format,os);
}

/** An overloaded version of print_all_custom() that prints to standard output.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics. */
template<class r_option>
void container_base<r_option>::print_all_custom(const char *format) {
	print_all_custom(format,cout);
}

/** An overloaded version of print_all_custom(), which outputs the result to a
 * particular file.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
void container_base<r_option>::print_all_custom(const char *format,const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all_custom(format,os);
	os.close();
}
/* Get residue level information */
template<class r_option>
		void container_base<r_option>::get_res_volume(
	vector<int>* res_type,
	vector<int>* res_id,
	vector<string>* at_type,
//Base Statistics
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
//Back Bone Statistics
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

	
	vector<int>  *parent,
	vector<int>  *child,
	vector<int>  *order

){

	fpoint x,y,z;
	voronoicell_neighbor c;	
//Set volume vector to zero
//res_type 0 = backbone, 1=G, 2=A, 3=U, 4=C 
for(int r=0; r<(res_id->back()+1);r++){
volume->push_back(0);
unpacked->push_back(0);
surface->push_back(0);
bvolume->push_back(0);
bunpacked->push_back(0);
bsurface->push_back(0);

orderG->push_back(0);
orderA->push_back(0);
orderU->push_back(0);
orderC->push_back(0);

borderG->push_back(0);
borderA->push_back(0);
borderU->push_back(0);
borderC->push_back(0);

surfG->push_back(0);
surfA->push_back(0);
surfU->push_back(0);
surfC->push_back(0);

bsurfG->push_back(0);
bsurfA->push_back(0);
bsurfU->push_back(0);
bsurfC->push_back(0);
//Initialize all atom types
vatom0->push_back(0);
vatom1->push_back(0);
vatom2->push_back(0);
vatom3->push_back(0);
vatom4->push_back(0);
vatom5->push_back(0);
vatom6->push_back(0);
vatom7->push_back(0);
vatom8->push_back(0);
vatom9->push_back(0);
vatom10->push_back(0);
vatom11->push_back(0);
vatom12->push_back(0);
vatom13->push_back(0);
vatom14->push_back(0);
vatom15->push_back(0);
vatom16->push_back(0);
vatom17->push_back(0);
vatom18->push_back(0);
vatom19->push_back(0);
vatom20->push_back(0);
vatom21->push_back(0);
vatom22->push_back(0);
vatom23->push_back(0);
vatom24->push_back(0);
vatom25->push_back(0);
vatom26->push_back(0);
atom0->push_back(0);
atom1->push_back(0);
atom2->push_back(0);
atom3->push_back(0);
atom4->push_back(0);
atom5->push_back(0);
atom6->push_back(0);
atom7->push_back(0);
atom8->push_back(0);
atom9->push_back(0);
atom10->push_back(0);
atom11->push_back(0);
atom12->push_back(0);
atom13->push_back(0);
atom14->push_back(0);
atom15->push_back(0);
atom16->push_back(0);
atom17->push_back(0);
atom18->push_back(0);
atom19->push_back(0);
atom20->push_back(0);
atom21->push_back(0);
atom22->push_back(0);
atom23->push_back(0);
atom24->push_back(0);
atom25->push_back(0);
atom26->push_back(0);
}
int i,j,k,ijk=0,q,fp;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) for(q=0;q<co[ijk];q++) {
		x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
		if(!compute_cell(c,i,j,k,ijk,q,x,y,z)) continue;
int atom_id = id[ijk][q];
vector<int> tempres;
vector<int> temporder;
int parent_id = res_id->at(atom_id);
//Get res_id
//printf("s:%f\n",c.surface_area());
vector<int> atom_ids;
vector<fpoint> surfAt;
c.output_face_areasV(&surfAt,&atom_ids);
//Identify unique residues that are associated with the particular residue
//Assign surface area to a particular residue
for(int jp=0; jp<atom_ids.size(); jp++){
if(atom_ids[jp]>=0){
//rid is the residue id
int rid1 = res_id->at(atom_ids[jp]);
int rid2 = res_id->at(atom_id);

if(surfAt[jp]>75){
surfAt[jp]=0;
}else{

if(rid1!=rid2){
int rtype=res_type->at(atom_ids[jp]);
int rtype2=res_type->at(atom_id);
if(rtype>0 && rtype2>0){
//If child is not already in the graph add it
int ingraph=0;
for(int jtmp=0; jtmp<tempres.size(); jtmp++){
if(tempres[jtmp]==rid1){ingraph=1;temporder[jtmp]=temporder[jtmp]++;}
}
if(ingraph==0){
tempres.push_back(rid1);
temporder.push_back(1);
}
}
}
}





if(rid1!=rid2){
int rtype=res_type->at(atom_ids[jp]);
int rtype2=res_type->at(atom_id);
//Only consider the base atoms
if(rtype2 > 0){ 
if(rtype==1){surfG->at(res_id->at(atom_id))=surfG->at(res_id->at(atom_id))+surfAt[jp];orderG->at(res_id->at(atom_id))++;surface->at(res_id->at(atom_id))=surface->at(res_id->at(atom_id))+surfAt[jp];}else{
if(rtype==2){surfA->at(res_id->at(atom_id))=surfA->at(res_id->at(atom_id))+surfAt[jp];orderA->at(res_id->at(atom_id))++;surface->at(res_id->at(atom_id))=surface->at(res_id->at(atom_id))+surfAt[jp];}else{
if(rtype==3){surfU->at(res_id->at(atom_id))=surfU->at(res_id->at(atom_id))+surfAt[jp];orderU->at(res_id->at(atom_id))++;surface->at(res_id->at(atom_id))=surface->at(res_id->at(atom_id))+surfAt[jp];}else{
if(rtype==4){surfC->at(res_id->at(atom_id))=surfC->at(res_id->at(atom_id))+surfAt[jp];orderC->at(res_id->at(atom_id))++;surface->at(res_id->at(atom_id))=surface->at(res_id->at(atom_id))+surfAt[jp];}else{
if(rtype==-1){bsurfG->at(res_id->at(atom_id))=bsurfG->at(res_id->at(atom_id))+surfAt[jp];borderG->at(res_id->at(atom_id))++;bsurface->at(res_id->at(atom_id))=bsurface->at(res_id->at(atom_id))+surfAt[jp];}else{
if(rtype==-2){bsurfA->at(res_id->at(atom_id))=bsurfA->at(res_id->at(atom_id))+surfAt[jp];borderA->at(res_id->at(atom_id))++;bsurface->at(res_id->at(atom_id))=bsurface->at(res_id->at(atom_id))+surfAt[jp];}else{
if(rtype==-3){bsurfU->at(res_id->at(atom_id))=bsurfU->at(res_id->at(atom_id))+surfAt[jp];borderU->at(res_id->at(atom_id))++;bsurface->at(res_id->at(atom_id))=bsurface->at(res_id->at(atom_id))+surfAt[jp];}else{
if(rtype==-4){bsurfC->at(res_id->at(atom_id))=bsurfC->at(res_id->at(atom_id))+surfAt[jp];borderC->at(res_id->at(atom_id))++;bsurface->at(res_id->at(atom_id))=bsurface->at(res_id->at(atom_id))+surfAt[jp];}else{
}
}
}
}
}
}
}
}


}
}

}

}

//Add temporary array elements to order, parent, and child
//If parent,child combination exists just add order

for(int jtmp=0; jtmp < tempres.size(); jtmp++){
int ingraph=0;

for(int ptmp=0; ptmp<parent->size(); ptmp++){
if(parent->at(ptmp)==parent_id && child->at(ptmp)==tempres[jtmp]){
order->at(ptmp)=temporder[jtmp]+order->at(ptmp);
ingraph=1;
}
}


if(ingraph==0){
parent->push_back(parent_id);
child->push_back(tempres[jtmp]);
order->push_back(temporder[jtmp]);
}


}

//Clear tmp vectors

tempres.clear();
temporder.clear();

if(c.surface_area()<75){
//Monster Conditional that fills atom vectors 
string aname=at_type->at(atom_id);
int restype=res_type->at(atom_id);
if(restype==1 || restype==-1){
if(aname.compare(" P  ")==0){atom0->at(res_id->at(atom_id))=1;vatom0->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" OP1")==0){atom1->at(res_id->at(atom_id))=1;vatom1->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" OP2")==0){atom2->at(res_id->at(atom_id))=1;vatom2->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O5'")==0){atom3->at(res_id->at(atom_id))=1;vatom3->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5'")==0){atom4->at(res_id->at(atom_id))=1;vatom4->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4'")==0){atom5->at(res_id->at(atom_id))=1;vatom5->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O4'")==0){atom6->at(res_id->at(atom_id))=1;vatom6->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C3'")==0){atom7->at(res_id->at(atom_id))=1;vatom7->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O3'")==0){atom8->at(res_id->at(atom_id))=1;vatom8->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2'")==0){atom9->at(res_id->at(atom_id))=1;vatom9->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O2'")==0){atom10->at(res_id->at(atom_id))=1;vatom10->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C1'")==0){atom11->at(res_id->at(atom_id))=1;vatom11->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N9 ")==0){atom12->at(res_id->at(atom_id))=1;vatom12->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4 ")==0){atom13->at(res_id->at(atom_id))=1;vatom13->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N3 ")==0){atom14->at(res_id->at(atom_id))=1;vatom14->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2 ")==0){atom15->at(res_id->at(atom_id))=1;vatom15->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N2 ")==0){atom16->at(res_id->at(atom_id))=1;vatom16->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N1 ")==0){atom17->at(res_id->at(atom_id))=1;vatom17->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C6 ")==0){atom18->at(res_id->at(atom_id))=1;vatom18->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O6 ")==0){atom19->at(res_id->at(atom_id))=1;vatom19->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5 ")==0){atom20->at(res_id->at(atom_id))=1;vatom20->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N7 ")==0){atom21->at(res_id->at(atom_id))=1;vatom21->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C8 ")==0){atom22->at(res_id->at(atom_id))=1;vatom22->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H1 ")==0){atom23->at(res_id->at(atom_id))=1;vatom23->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H8 ")==0){atom24->at(res_id->at(atom_id))=1;vatom24->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H21")==0){atom25->at(res_id->at(atom_id))=1;vatom25->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H22")==0){atom26->at(res_id->at(atom_id))=1;vatom26->at(res_id->at(atom_id))=c.volume();}
}else{                                                        
if(restype==2 || restype==-2){                                
if(aname.compare(" P  ")==0){atom0->at(res_id->at(atom_id))=1;vatom0->at(res_id->at(atom_id))=c.volume();}                             
if(aname.compare(" OP1")==0){atom1->at(res_id->at(atom_id))=1;vatom1->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" OP2")==0){atom2->at(res_id->at(atom_id))=1;vatom2->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O5'")==0){atom3->at(res_id->at(atom_id))=1;vatom3->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5'")==0){atom4->at(res_id->at(atom_id))=1;vatom4->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4'")==0){atom5->at(res_id->at(atom_id))=1;vatom5->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O4'")==0){atom6->at(res_id->at(atom_id))=1;vatom6->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C3'")==0){atom7->at(res_id->at(atom_id))=1;vatom7->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O3'")==0){atom8->at(res_id->at(atom_id))=1;vatom8->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2'")==0){atom9->at(res_id->at(atom_id))=1;vatom9->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O2'")==0){atom10->at(res_id->at(atom_id))=1;vatom10->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C1'")==0){atom11->at(res_id->at(atom_id))=1;vatom11->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N9 ")==0){atom12->at(res_id->at(atom_id))=1;vatom12->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4 ")==0){atom13->at(res_id->at(atom_id))=1;vatom13->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N3 ")==0){atom14->at(res_id->at(atom_id))=1;vatom14->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2 ")==0){atom15->at(res_id->at(atom_id))=1;vatom15->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N1 ")==0){atom16->at(res_id->at(atom_id))=1;vatom16->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C6 ")==0){atom17->at(res_id->at(atom_id))=1;vatom17->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N6 ")==0){atom18->at(res_id->at(atom_id))=1;vatom18->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5 ")==0){atom19->at(res_id->at(atom_id))=1;vatom19->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N7 ")==0){atom20->at(res_id->at(atom_id))=1;vatom20->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C8 ")==0){atom21->at(res_id->at(atom_id))=1;vatom21->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H2 ")==0){atom22->at(res_id->at(atom_id))=1;vatom22->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H8 ")==0){atom23->at(res_id->at(atom_id))=1;vatom23->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H61")==0){atom24->at(res_id->at(atom_id))=1;vatom24->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H62")==0){atom25->at(res_id->at(atom_id))=1;vatom25->at(res_id->at(atom_id))=c.volume();}
}else{
if(restype==-3 || restype==3){
//Uracil Atom Names
if(aname.compare(" P  ")==0){atom0->at(res_id->at(atom_id))=1;vatom0->at(res_id->at(atom_id))=c.volume();}                             
if(aname.compare(" OP1")==0){atom1->at(res_id->at(atom_id))=1;vatom1->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" OP2")==0){atom2->at(res_id->at(atom_id))=1;vatom2->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O5'")==0){atom3->at(res_id->at(atom_id))=1;vatom3->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5'")==0){atom4->at(res_id->at(atom_id))=1;vatom4->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4'")==0){atom5->at(res_id->at(atom_id))=1;vatom5->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O4'")==0){atom6->at(res_id->at(atom_id))=1;vatom6->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C3'")==0){atom7->at(res_id->at(atom_id))=1;vatom7->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O3'")==0){atom8->at(res_id->at(atom_id))=1;vatom8->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2'")==0){atom9->at(res_id->at(atom_id))=1;vatom9->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O2'")==0){atom10->at(res_id->at(atom_id))=1;vatom10->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C1'")==0){atom11->at(res_id->at(atom_id))=1;vatom11->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N1 ")==0){atom12->at(res_id->at(atom_id))=1;vatom12->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2 ")==0){atom13->at(res_id->at(atom_id))=1;vatom13->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C6 ")==0){atom14->at(res_id->at(atom_id))=1;vatom14->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O2 ")==0){atom15->at(res_id->at(atom_id))=1;vatom15->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N3 ")==0){atom16->at(res_id->at(atom_id))=1;vatom16->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4 ")==0){atom17->at(res_id->at(atom_id))=1;vatom17->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O4 ")==0){atom18->at(res_id->at(atom_id))=1;vatom18->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5 ")==0){atom19->at(res_id->at(atom_id))=1;vatom19->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H3 ")==0){atom20->at(res_id->at(atom_id))=1;vatom20->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H5 ")==0){atom21->at(res_id->at(atom_id))=1;vatom21->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H6 ")==0){atom22->at(res_id->at(atom_id))=1;vatom22->at(res_id->at(atom_id))=c.volume();}
}else{
//Cytosine Atom Names
if(restype==4 || restype==-4){
if(aname.compare(" P  ")==0){atom0->at(res_id->at(atom_id))=1;vatom0->at(res_id->at(atom_id))=c.volume();}                             
if(aname.compare(" OP1")==0){atom1->at(res_id->at(atom_id))=1;vatom1->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" OP2")==0){atom2->at(res_id->at(atom_id))=1;vatom2->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O5'")==0){atom3->at(res_id->at(atom_id))=1;vatom3->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5'")==0){atom4->at(res_id->at(atom_id))=1;vatom4->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4'")==0){atom5->at(res_id->at(atom_id))=1;vatom5->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O4'")==0){atom6->at(res_id->at(atom_id))=1;vatom6->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C3'")==0){atom7->at(res_id->at(atom_id))=1;vatom7->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O3'")==0){atom8->at(res_id->at(atom_id))=1;vatom8->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2'")==0){atom9->at(res_id->at(atom_id))=1;vatom9->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O2'")==0){atom10->at(res_id->at(atom_id))=1;vatom10->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C1'")==0){atom11->at(res_id->at(atom_id))=1;vatom11->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N1 ")==0){atom12->at(res_id->at(atom_id))=1;vatom12->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C2 ")==0){atom13->at(res_id->at(atom_id))=1;vatom13->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C6 ")==0){atom14->at(res_id->at(atom_id))=1;vatom14->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" O2 ")==0){atom15->at(res_id->at(atom_id))=1;vatom15->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N3 ")==0){atom16->at(res_id->at(atom_id))=1;vatom16->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C4 ")==0){atom17->at(res_id->at(atom_id))=1;vatom17->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" N4 ")==0){atom18->at(res_id->at(atom_id))=1;vatom18->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" C5 ")==0){atom19->at(res_id->at(atom_id))=1;vatom19->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H5 ")==0){atom20->at(res_id->at(atom_id))=1;vatom20->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H6 ")==0){atom21->at(res_id->at(atom_id))=1;vatom21->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H41")==0){atom22->at(res_id->at(atom_id))=1;vatom22->at(res_id->at(atom_id))=c.volume();}
if(aname.compare(" H42")==0){atom23->at(res_id->at(atom_id))=1;vatom22->at(res_id->at(atom_id))=c.volume();}
}
}
}
}
if(res_type->at(atom_id)>0){
volume->at(res_id->at(atom_id))=volume->at(res_id->at(atom_id))+c.volume();
}else{
bvolume->at(res_id->at(atom_id))=bvolume->at(res_id->at(atom_id))+c.volume();
}
}else{
if(res_type->at(atom_id)>0){
unpacked->at(res_id->at(atom_id))=unpacked->at(res_id->at(atom_id))+1;
}else{
bunpacked->at(res_id->at(atom_id))=bunpacked->at(res_id->at(atom_id))+1;
}

}
}
}

/** The internal part of the print_all_custom() routine, that can be called
 * with either a voronoicell class (if no neighbor computations are needed) or
 * with a voronoicell_neighbor class (if neighbor computations are needed).
 * \param[in,out] c a Voronoi cell object to use for the computation.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] os an open output stream to write to. */
template<class r_option>
template<class n_option>
void container_base<r_option>::print_all_custom_internal(voronoicell_base<n_option> &c,const char *format,ostream &os) {
	fpoint x,y,z;
	int i,j,k,ijk=0,q,fp;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) for(q=0;q<co[ijk];q++) {
		x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
		if(!compute_cell(c,i,j,k,ijk,q,x,y,z)) continue;
		fp=0;
		while(format[fp]!=0) {
			if(format[fp]=='%') {
				fp++;
				switch(format[fp]) {

					// Particle-related output
					case 'i': os << id[ijk][q];break;
					case 'x': os << x;break;
					case 'y': os << y;break;
					case 'z': os << z;break;
					case 'q': os << x << " " << y << " " << z;break;
					case 'r': radius.print(os,ijk,q,false);break;

					// Vertex-related output
					case 'w': os << c.p;break;
					case 'p': c.output_vertices(os);break;
					case 'P': c.output_vertices(os,x,y,z);break;
					case 'o': c.output_vertex_orders(os);
					case 'm': os << 0.25*c.max_radius_squared();break;

					// Edge-related output
					case 'g': os << c.number_of_edges();break;
					case 'E': os << c.total_edge_distance();break;
					case 'e': c.output_face_perimeters(os);break;

					// Face-related output
					case 's': os << c.number_of_faces();break;
					case 'F': os << c.surface_area();break;
					case 'A': c.output_face_freq_table(os);break;
					case 'a': c.output_face_orders(os);break;
					case 'f': c.output_face_areas(os);break;
					case 't': c.output_face_vertices(os);break;
					case 'l': c.output_normals(os);break;
					case 'n': c.output_neighbors(os);break;

					// Volume-related output
					case 'v': os << c.volume();break;
					case 'c': {
							  fpoint cx,cy,cz;
							  c.centroid(cx,cy,cz);
							  os << cx << " " << cy << " " << cz;
						  } break;
					case 'C': {
							  fpoint cx,cy,cz;
							  c.centroid(cx,cy,cz);
							  os << x+cx << " " << y+cy << " " << z+cz;
						  } break;

					// End-of-string reached
					case 0: fp--;break;

					// The percent sign is not part of a
					// control sequence
					default: os << '%' << format[fp];
				}
			} else os << format[fp];
			fp++;
		}
		os << "\n";
	}
}

/** The internal part of the print_all() and print_all_neighbor() routines, that
 * computes all of the Voronoi cells, and then outputs simple information about
 * them. The routine can be called with either a voronoicell class (if no
 * neighbor computations are needed) or with a voronoicell_neighbor class (if
 * neighbor computations are needed).
 * \param[in,out] c a Voronoi cell object to use for the computation.
 * \param[in] os an open output stream to write to. */
template<class r_option>
template<class n_option>
inline void container_base<r_option>::print_all_internal(voronoicell_base<n_option> &c,ostream &os) {
	fpoint x,y,z;
	int i,j,k,ijk=0,q;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(q=0;q<co[ijk];q++) {
			x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
			os << id[ijk][q] << " " << x << " " << y << " " << z;
			radius.print(os,ijk,q);
			if(compute_cell(c,i,j,k,ijk,q,x,y,z)) {
				os << " " << c.volume();
				c.output_neighbors(os,true);
				os << "\n";
			} else os << " 0\n";
		}
	}
}

/** Prints a list of all particle labels, positions, and Voronoi volumes to the
 * standard output.
 * \param[in] os the output stream to print to. */
template<class r_option>
void container_base<r_option>::print_all(ostream &os) {
	voronoicell c;
	print_all_internal(c,os);
}

/** An overloaded version of print_all(), which just prints to standard output. */
template<class r_option>
void container_base<r_option>::print_all() {
	voronoicell c;
	print_all_internal(c,cout);
}

/** An overloaded version of print_all(), which outputs the result to a particular
 * file.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
inline void container_base<r_option>::print_all(const char* filename) {
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all_internal(c,os);
	os.close();
}

/** Prints a list of all particle labels, positions, Voronoi volumes, and a list
 * of neighboring particles to an output stream.
 * \param[in] os the output stream to print to.*/
template<class r_option>
void container_base<r_option>::print_all_neighbor(ostream &os) {
	voronoicell_neighbor c;
	print_all_internal(c,os);
}

/** An overloaded version of print_all_neighbor(), which just prints to
 * standard output. */
template<class r_option>
void container_base<r_option>::print_all_neighbor() {
	voronoicell_neighbor c;
	print_all_internal(c,cout);
}

/** An overloaded version of print_all_neighbor(), which outputs the results to
 * a particular file.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
inline void container_base<r_option>::print_all_neighbor(const char* filename) {
	voronoicell_neighbor c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all_internal(c,os);
	os.close();
}

/** Initialize the Voronoi cell to be the entire container. For non-periodic
 * coordinates, this is set by the position of the walls. For periodic
 * coordinates, the space is equally divided in either direction from the
 * particle's initial position. That makes sense since those boundaries would
 * be made by the neighboring periodic images of this particle. It also applies
 * plane cuts made by any walls that have been added to the container.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (x,y,z) the position of the particle.
 * \return False if the plane cuts applied by walls completely removed the
 *         cell, true otherwise. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::initialize_voronoicell(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
	fpoint x1,x2,y1,y2,z1,z2;
	if(xperiodic) x1=-(x2=0.5*(bx-ax));else {x1=ax-x;x2=bx-x;}
	if(yperiodic) y1=-(y2=0.5*(by-ay));else {y1=ay-y;y2=by-y;}
	if(zperiodic) z1=-(z2=0.5*(bz-az));else {z1=az-z;z2=bz-z;}
	c.init(x1,x2,y1,y2,z1,z2);
	for(int j=0;j<wall_number;j++) if(!(walls[j]->cut_cell(c,x,y,z))) return false;
	return true;
}

/** This function tests to see if a given vector lies within the container
 * bounds and any walls.
 * \param[in] (x,y,z) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 *         outside. */
template<class r_option>
bool container_base<r_option>::point_inside(fpoint x,fpoint y,fpoint z) {
	if(x<ax||x>bx||y<ay||y>by||z<az||z>bz) return false;
	return point_inside_walls(x,y,z);
}

/** This function tests to see if a give vector lies within the walls that have
 * been added to the container, but does not specifically check whether the
 * vector lies within the container bounds.
 * \param[in] (x,y,z) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 *         outside. */
template<class r_option>
bool container_base<r_option>::point_inside_walls(fpoint x,fpoint y,fpoint z) {
	for(int j=0;j<wall_number;j++) if(!walls[j]->point_inside(x,y,z)) return false;
	return true;
}

/** This routine is a simpler alternative to compute_cell(), that constructs
 * the cell by testing over successively larger spherical shells of particles.
 * For a container that is homogeneously filled with particles, this routine
 * runs as fast as compute_cell(). However, it rapidly becomes inefficient
 * for cases when the particles are not homogeneously distributed, or where
 * parts of the container might not be filled. In that case, the spheres may
 * grow very large before being cut off, leading to poor efficiency.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \param[in] (x,y,z) the coordinates of the particle.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
bool container_base<r_option>::compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z) {

	// This length scale determines how large the spherical shells should
	// be, and it should be set to approximately the particle diameter
	const fpoint length_scale=1;
	fpoint x1,y1,z1,qx,qy,qz,lr=0,lrs=0,ur,urs,rs;
	int q,t;
	voropp_loop l(this);
	if(!initialize_voronoicell(c,x,y,z)) return false;

	// Now the cell is cut by testing neighboring particles in concentric
	// shells. Once the test shell becomes twice as large as the Voronoi
	// cell we can stop testing.
	radius.init(ijk,s);
	while(radius.cutoff(lrs)<c.max_radius_squared()) {
		ur=lr+0.5*length_scale;urs=ur*ur;
		t=l.init(x,y,z,ur,qx,qy,qz);
		do {
			for(q=0;q<co[t];q++) {
				x1=p[t][sz*q]+qx-x;y1=p[t][sz*q+1]+qy-y;z1=p[t][sz*q+2]+qz-z;
				rs=x1*x1+y1*y1+z1*z1;
				if(lrs-tolerance<rs&&rs<urs&&(q!=s||ijk!=t)) {
					if(!c.nplane(x1,y1,z1,radius.scale(rs,t,q),id[t][q])) return false;
				}
			}
		} while((t=l.inc(qx,qy,qz))!=-1);
		lr=ur;lrs=urs;
	}
	return true;
}

/** A overloaded version of compute_cell_sphere(), that sets up the x, y, and z
 * variables.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s) {
	fpoint x=p[ijk][sz*s],y=p[ijk][sz*s+1],z=p[ijk][sz*s+2];
	return compute_cell_sphere(c,i,j,k,ijk,s,x,y,z);
}

/** A overloaded version of compute_cell, that sets up the x, y, and z variables.
 * It can be run by the user, and it is also called multiple times by the
 * functions print_all(), store_cell_volumes(), and the output routines.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s) {
	fpoint x=p[ijk][sz*s],y=p[ijk][sz*s+1],z=p[ijk][sz*s+2];
	return  compute_cell(c,i,j,k,ijk,s,x,y,z);
}

/** This routine computes a Voronoi cell for a single particle in the
 * container. It can be called by the user, but is also forms the core part of
 * several of the main functions, such as store_cell_volumes(), print_all(),
 * and the drawing routines. The algorithm constructs the cell by testing over
 * the neighbors of the particle, working outwards until it reaches those
 * particles which could not possibly intersect the cell. For maximum
 * efficiency, this algorithm is divided into three parts. In the first
 * section, the algorithm tests over the blocks which are in the immediate
 * vicinity of the particle, by making use of one of the precomputed worklists.
 * The code then continues to test blocks on the worklist, but also begins to
 * construct a list of neighboring blocks outside the worklist which may need
 * to be test. In the third section, the routine starts testing these
 * neighboring blocks, evaluating whether or not a particle in them could
 * possibly intersect the cell. For blocks that intersect the cell, it tests
 * the particles in that block, and then adds the block neighbors to the list
 * of potential places to consider.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \param[in] (x,y,z) the coordinates of the particle.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
bool container_base<r_option>::compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z) {
	const fpoint boxx=(bx-ax)/nx,boxy=(by-ay)/ny,boxz=(bz-az)/nz;
	fpoint x1,y1,z1,qx=0,qy=0,qz=0;
	fpoint xlo,ylo,zlo,xhi,yhi,zhi,rs;
	int ci,cj,ck,di,dj,dk,dijk,ei,ej,ek,eijk,si,sj,sk,sijk;
	fpoint gxs,gys,gzs,*radp;
	int f,g,l;unsigned int q,*e;
	const unsigned int b1=1<<21,b2=1<<22,b3=1<<24,b4=1<<25,b5=1<<27,b6=1<<28;

	radius.init(ijk,s);

	// Initialize the Voronoi cell to fill the entire container
	if(!initialize_voronoicell(c,x,y,z)) return false;
	fpoint crs,mrs;

	int next_count=3,list_index=0,list_size=8;
	int count_list[]={7,11,15,19,26,35,45,59};

	// Test all particles in the particle's local region first
	for(l=0;l<s;l++) {
		x1=p[ijk][sz*l]-x;
		y1=p[ijk][sz*l+1]-y;
		z1=p[ijk][sz*l+2]-z;
		rs=radius.scale(x1*x1+y1*y1+z1*z1,ijk,l);
		if(!c.nplane(x1,y1,z1,rs,id[ijk][l])) return false;
	}
	l++;
	while(l<co[ijk]) {
		x1=p[ijk][sz*l]-x;
		y1=p[ijk][sz*l+1]-y;
		z1=p[ijk][sz*l+2]-z;
		rs=radius.scale(x1*x1+y1*y1+z1*z1,ijk,l);
		if(!c.nplane(x1,y1,z1,rs,id[ijk][l])) return false;
		l++;
	}

	// Now compute the maximum distance squared from the cell center to a
	// vertex. This is used to cut off the calculation since we only need
	// to test out to twice this range.
	mrs=c.max_radius_squared();

	// Now compute the fractional position of the particle within its
	// region and store it in (fx,fy,fz). We use this to compute an index
	// (si,sj,sk) of which subregion the particle is within.
	unsigned int m1,m2;
	fpoint fx=x-ax-boxx*i,fy=y-ay-boxy*j,fz=z-az-boxz*k;
	si=int(fx*xsp*fgrid);sj=int(fy*ysp*fgrid);sk=int(fz*zsp*fgrid);

	// The indices (si,sj,sk) tell us which worklist to use, to test the
	// blocks in the optimal order. But we only store worklists for the
	// eighth of the region where si, sj, and sk are all less than half the
	// full grid. The rest of the cases are handled by symmetry. In this
	// section, we detect for these cases, by reflecting high values of si,
	// sj, and sk. For these cases, a mask is constructed in m1 and m2
	// which is used to flip the worklist information when it is loaded.
	if(si>=hgrid) {
		gxs=fx;
		m1=127+(3<<21);si=fgrid-1-si;m2=1+(1<<21);
	} else {m1=m2=0;gxs=boxx-fx;}
	if(sj>=hgrid) {
		gys=fy;
		m1|=(127<<7)+(3<<24);sj=fgrid-1-sj;m2|=(1<<7)+(1<<24);
	} else gys=boxy-fy;
	if(sk>=hgrid) {
		gzs=fz;
		m1|=(127<<14)+(3<<27);sk=fgrid-1-sk;m2|=(1<<14)+(1<<27);
	} else gzs=boxz-fz;
	gxs*=gxs;gys*=gys;gzs*=gzs;

	// It's possible that a problem with the int() function could lead to
	// spurious values with particles lying on the boundaries of the regions.
	// In this section we correct for that.
	if(si<0) si=0;if(sj<0) sj=0;if(sk<0) sk=0;

	// Now compute which worklist we are going to use, and set radp and e to
	// point at the right offsets
	sijk=si+hgrid*(sj+hgrid*sk);
	radp=mrad+sijk*seq_length;
	e=(const_cast<unsigned int*> (wl))+sijk*seq_length;

	// Read in how many items in the worklist can be tested without having to
	// worry about writing to the mask
	f=e[0];g=0;
	do {

		// At the intervals specified by count_list, we recompute the
		// maximum radius squared
		if(g==next_count) {
			mrs=c.max_radius_squared();
			if(list_index!=list_size) next_count=count_list[list_index++];
		}

		// If mrs is less than the minimum distance to any untested
		// block, then we are done
		if(mrs<radius.cutoff(radp[g])) return true;
		g++;

		// Load in a block off the worklist, permute it with the
		// symmetry mask, and decode its position. These are all
		// integer bit operations so they should run very fast.
		q=e[g];q^=m1;q+=m2;
		di=q&127;di-=64;
		dj=(q>>7)&127;dj-=64;
		dk=(q>>14)&127;dk-=64;

		// Check that the worklist position is in range
		if(xperiodic) {if(di<-nx) continue;else if(di>nx) continue;}
		else {if(di<-i) continue;else if(di>=nx-i) continue;}
		if(yperiodic) {if(dj<-ny) continue;else if(dj>ny) continue;}
		else {if(dj<-j) continue;else if(dj>=ny-j) continue;}
		if(zperiodic) {if(dk<-nz) continue;else if(dk>nz) continue;}
		else {if(dk<-k) continue;else if(dk>=nz-k) continue;}

		// Call the compute_min_max_radius() function. This returns
		// true if the minimum distance to the block is bigger than the
		// current mrs, in which case we skip this block and move on.
		// Otherwise, it computes the maximum distance to the block and
		// returns it in crs.
		if(compute_min_max_radius(di,dj,dk,fx,fy,fz,gxs,gys,gzs,crs,mrs)) continue;

		// Now compute which region we are going to loop over, adding a
		// displacement for the periodic cases
		di+=i;dj+=j;dk+=k;
		if(xperiodic) {if(di<0) {qx=ax-bx;di+=nx;} else if(di>=nx) {qx=bx-ax;di-=nx;} else qx=0;}
		if(yperiodic) {if(dj<0) {qy=ay-by;dj+=ny;} else if(dj>=ny) {qy=by-ay;dj-=ny;} else qy=0;}
		if(zperiodic) {if(dk<0) {qz=az-bz;dk+=nz;} else if(dk>=nz) {qz=bz-az;dk-=nz;} else qz=0;}
		dijk=di+nx*(dj+ny*dk);

		// If mrs is bigger than the maximum distance to the block,
		// then we have to test all particles in the block for
		// intersections. Otherwise, we do additional checks and skip
		// those particles which can't possibly intersect the block.
		if(mrs>radius.cutoff(crs)) {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]+qy-y;
				z1=p[dijk][sz*l+2]+qz-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
			}
		} else {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]+qy-y;
				z1=p[dijk][sz*l+2]+qz-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(rs<mrs) {
					if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
				}
			}
		}
	} while(g<f);

	// If we reach here, we were unable to compute the entire cell using
	// the first part of the worklist. This section of the algorithm
	// continues the worklist, but it now starts preparing the mask that we
	// need if we end up going block by block. We do the same as before,
	// but we put a mark down on the mask for every block that's tested.
	// The worklist also contains information about which neighbors of each
	// block are not also on the worklist, and we start storing those
	// points in a list in case we have to go block by block.
	ci=xperiodic?nx:i;
	cj=yperiodic?ny:j;
	ck=zperiodic?nz:k;

	// Update the mask counter, and if it wraps around then reset the
	// whole mask; that will only happen once every 2^32 tries
	mv++;
	if(mv==0) {
		for(l=0;l<hxyz;l++) mask[l]=0;
		mv=1;
	}

	// Reset the block by block counters
	s_start=s_end=0;

	while(g<seq_length-1) {

		// At the intervals specified by count_list, we recompute the
		// maximum radius squared
		if(g==next_count) {
			mrs=c.max_radius_squared();
			if(list_index!=list_size) next_count=count_list[list_index++];
		}

		// If mrs is less than the minimum distance to any untested
		// block, then we are done
		if(mrs<radius.cutoff(radp[g])) return true;
		g++;

		// Load in a block off the worklist, permute it with the
		// symmetry mask, and decode its position. These are all
		// integer bit operations so they should run very fast.
		q=e[g];q^=m1;q+=m2;
		di=q&127;di-=64;
		dj=(q>>7)&127;dj-=64;
		dk=(q>>14)&127;dk-=64;

		// Compute the position in the mask of the current block. If
		// this lies outside the mask, then skip it. Otherwise, mark
		// it.
		ei=ci+di;
		ej=cj+dj;
		ek=ck+dk;
		if(ei<0) continue;else if(ei>=hx) continue;
		if(ej<0) continue;else if(ej>=hy) continue;
		if(ek<0) continue;else if(ek>=hz) continue;
		eijk=ei+hx*(ej+hy*ek);
		mask[eijk]=mv;

		// Call the compute_min_max_radius() function. This returns
		// true if the minimum distance to the block is bigger than the
		// current mrs, in which case we skip this block and move on.
		// Otherwise, it computes the maximum distance to the block and
		// returns it in crs.
		if(compute_min_max_radius(di,dj,dk,fx,fy,fz,gxs,gys,gzs,crs,mrs)) continue;

		// Now compute which region we are going to loop over, adding a
		// displacement for the periodic cases
		di+=i;dj+=j;dk+=k;
		if(xperiodic) {if(di<0) {qx=ax-bx;di+=nx;} else if(di>=nx) {qx=bx-ax;di-=nx;} else qx=0;}
		if(yperiodic) {if(dj<0) {qy=ay-by;dj+=ny;} else if(dj>=ny) {qy=by-ay;dj-=ny;} else qy=0;}
		if(zperiodic) {if(dk<0) {qz=az-bz;dk+=nz;} else if(dk>=nz) {qz=bz-az;dk-=nz;} else qz=0;}
		dijk=di+nx*(dj+ny*dk);

		// If mrs is bigger than the maximum distance to the block,
		// then we have to test all particles in the block for
		// intersections. Otherwise, we do additional checks and skip
		// those particles which can't possibly intersect the block.
		if(mrs>radius.cutoff(crs)) {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]+qy-y;
				z1=p[dijk][sz*l+2]+qz-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
			}
		} else {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]+qy-y;
				z1=p[dijk][sz*l+2]+qz-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(rs<mrs) {
					if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
				}
			}
		}

		// If there might not be enough memory on the list for these
		// additions, then add more
		if(s_end+18>s_size) add_list_memory();

		// Test the parts of the worklist element which tell us what
		// neighbors of this block are not on the worklist. Store them
		// on the block list, and mark the mask.
		if((q&b2)==b2) {
			if(ei>0) if(mask[eijk-1]!=mv) {mask[eijk-1]=mv;sl[s_end++]=ei-1;sl[s_end++]=ej;sl[s_end++]=ek;}
			if((q&b1)==0) if(ei<hx-1) if(mask[eijk+1]!=mv) {mask[eijk+1]=mv;sl[s_end++]=ei+1;sl[s_end++]=ej;sl[s_end++]=ek;}
		} else if((q&b1)==b1) {if(ei<hx-1) if(mask[eijk+1]!=mv) {mask[eijk+1]=mv;sl[s_end++]=ei+1;sl[s_end++]=ej;sl[s_end++]=ek;}}
		if((q&b4)==b4) {if(ej>0) if(mask[eijk-hx]!=mv) {mask[eijk-hx]=mv;sl[s_end++]=ei;sl[s_end++]=ej-1;sl[s_end++]=ek;}
			if((q&b3)==0) if(ej<hy-1) if(mask[eijk+hx]!=mv) {mask[eijk+hx]=mv;sl[s_end++]=ei;sl[s_end++]=ej+1;sl[s_end++]=ek;}
		} else if((q&b3)==b3) {if(ej<hy-1) if(mask[eijk+hx]!=mv) {mask[eijk+hx]=mv;sl[s_end++]=ei;sl[s_end++]=ej+1;sl[s_end++]=ek;}}
		if((q&b6)==b6) {if(ek>0) if(mask[eijk-hxy]!=mv) {mask[eijk-hxy]=mv;sl[s_end++]=ei;sl[s_end++]=ej;sl[s_end++]=ek-1;}
			if((q&b5)==0) if(ek<hz-1) if(mask[eijk+hxy]!=mv) {mask[eijk+hxy]=mv;sl[s_end++]=ei;sl[s_end++]=ej;sl[s_end++]=ek+1;}
		} else if((q&b5)==b5) if(ek<hz-1) if(mask[eijk+hxy]!=mv) {mask[eijk+hxy]=mv;sl[s_end++]=ei;sl[s_end++]=ej;sl[s_end++]=ek+1;}
	}

	// Do a check to see if we've reach the radius cutoff
	if(mrs<radius.cutoff(radp[g])) return true;

	// Update the mask counter, and if it has wrapped around, then
	// reset the mask
	fx+=boxx*ci;fy+=boxy*cj;fz+=boxz*ck;

	// We were unable to completely compute the cell based on the blocks in
	// the worklist, so now we have to go block by block, reading in items
	// off the list
	while(s_start!=s_end) {

		// If we reached the end of the list memory loop back to the
		// start
		if(s_start==s_size) s_start=0;

		// Read in a block off the list, and compute the upper and lower
		// coordinates in each of the three dimensions
		di=sl[s_start++];dj=sl[s_start++];dk=sl[s_start++];
		xlo=di*boxx-fx;xhi=xlo+boxx;
		ylo=dj*boxy-fy;yhi=ylo+boxy;
		zlo=dk*boxz-fz;zhi=zlo+boxz;

		// Carry out plane tests to see if any particle in this block
		// could possibly intersect the cell
		if(di>ci) {
			if(dj>cj) {
				if(dk>ck) {
					if(corner_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if(edge_z_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if(corner_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xlo,yhi,zhi,xhi,ylo,zlo)) continue;
				} else {
					if(edge_z_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if(edge_y_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(edge_y_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if(face_x_test(c,xlo,ylo,zlo,yhi,zhi)) continue;
				}
			}
		} else if(di<ci) {
			if(dj>cj) {
				if(dk>ck) {
					if(corner_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xhi,ylo,zhi,xlo,yhi,zlo)) continue;
				} else {
					if(edge_z_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if(corner_test(c,xhi,yhi,zlo,xlo,ylo,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xhi,yhi,zhi,xlo,ylo,zlo)) continue;
				} else {
					if(edge_z_test(c,xhi,yhi,zlo,xlo,ylo,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if(edge_y_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(edge_y_test(c,xhi,ylo,zhi,xlo,yhi,zlo)) continue;
				} else {
					if(face_x_test(c,xhi,ylo,zlo,yhi,zhi)) continue;
				}
			}
		} else {
			if(dj>cj) {
				if(dk>ck) {
					if(edge_x_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(edge_x_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if(face_y_test(c,xlo,ylo,zlo,xhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if(edge_x_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				} else if(dk<ck) {
					if(edge_x_test(c,xlo,yhi,zhi,xhi,ylo,zlo)) continue;
				} else {
					if(face_y_test(c,xlo,yhi,zlo,xhi,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if(face_z_test(c,xlo,ylo,zlo,xhi,yhi)) continue;
				} else if(dk<ck) {
					if(face_z_test(c,xlo,ylo,zhi,xhi,yhi)) continue;
				} else voropp_fatal_error("Compute cell routine revisiting central block, which should never\nhappen.",VOROPP_INTERNAL_ERROR);
			}
		}

		// Now compute the region that we are going to test over, and
		// set a displacement vector for the periodic cases
		if(xperiodic) {ei=i+di-nx;if(ei<0) {qx=ax-bx;ei+=nx;} else if(ei>=nx) {qx=bx-ax;ei-=nx;} else qx=0;} else ei=di;
		if(yperiodic) {ej=j+dj-ny;if(ej<0) {qy=ay-by;ej+=ny;} else if(ej>=ny) {qy=by-ay;ej-=ny;} else qy=0;} else ej=dj;
		if(zperiodic) {ek=k+dk-nz;if(ek<0) {qz=az-bz;ek+=nz;} else if(ek>=nz) {qz=bz-az;ek-=nz;} else qz=0;} else ek=dk;
		eijk=ei+nx*(ej+ny*ek);

		// Loop over all the elements in the block to test for cuts. It
		// would be possible to exclude some of these cases by testing
		// against mrs, but this will probably not save time.
		for(l=0;l<co[eijk];l++) {
			x1=p[eijk][sz*l]+qx-x;
			y1=p[eijk][sz*l+1]+qy-y;
			z1=p[eijk][sz*l+2]+qz-z;
			rs=radius.scale(x1*x1+y1*y1+z1*z1,eijk,l);
			if(!c.nplane(x1,y1,z1,rs,id[eijk][l])) return false;
		}

		// If there's not much memory on the block list then add more
		if((s_start<=s_end?s_size-s_end+s_start:s_end-s_start)<18) add_list_memory();

		// Test the neighbors of the current block, and add them to the
		// block list if they haven't already been tested
		dijk=di+hx*(dj+hy*dk);
		if(di>0) if(mask[dijk-1]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-1]=mv;sl[s_end++]=di-1;sl[s_end++]=dj;sl[s_end++]=dk;}
		if(dj>0) if(mask[dijk-hx]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-hx]=mv;sl[s_end++]=di;sl[s_end++]=dj-1;sl[s_end++]=dk;}
		if(dk>0) if(mask[dijk-hxy]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-hxy]=mv;sl[s_end++]=di;sl[s_end++]=dj;sl[s_end++]=dk-1;}
		if(di<hx-1) if(mask[dijk+1]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+1]=mv;sl[s_end++]=di+1;sl[s_end++]=dj;sl[s_end++]=dk;}
		if(dj<hy-1) if(mask[dijk+hx]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+hx]=mv;sl[s_end++]=di;sl[s_end++]=dj+1;sl[s_end++]=dk;}
		if(dk<hz-1) if(mask[dijk+hxy]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+hxy]=mv;sl[s_end++]=di;sl[s_end++]=dj;sl[s_end++]=dk+1;}
	}

	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is at a corner.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (xl,yl,zl) the relative coordinates of the corner of the block
 *                       closest to the cell center.
 * \param[in] (xh,yh,zh) the relative coordinates of the corner of the block
 *                       furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::corner_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint zl,fpoint xh,fpoint yh,fpoint zh) {
	if(c.plane_intersects_guess(xh,yl,zl,radius.cutoff(xl*xh+yl*yl+zl*zl))) return false;
	if(c.plane_intersects(xh,yh,zl,radius.cutoff(xl*xh+yl*yh+zl*zl))) return false;
	if(c.plane_intersects(xl,yh,zl,radius.cutoff(xl*xl+yl*yh+zl*zl))) return false;
	if(c.plane_intersects(xl,yh,zh,radius.cutoff(xl*xl+yl*yh+zl*zh))) return false;
	if(c.plane_intersects(xl,yl,zh,radius.cutoff(xl*xl+yl*yl+zl*zh))) return false;
	if(c.plane_intersects(xh,yl,zh,radius.cutoff(xl*xh+yl*yl+zl*zh))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on an edge which points along the x
 * direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
 *                    block.
 * \param[in] (yl,zl) the relative y and z coordinates of the corner of the
 *                    block closest to the cell center.
 * \param[in] (yh,zh) the relative y and z coordinates of the corner of the
 *                    block furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::edge_x_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint zl,fpoint x1,fpoint yh,fpoint zh) {
	if(c.plane_intersects_guess(x0,yl,zh,radius.cutoff(yl*yl+zl*zh))) return false;
	if(c.plane_intersects(x1,yl,zh,radius.cutoff(yl*yl+zl*zh))) return false;
	if(c.plane_intersects(x1,yl,zl,radius.cutoff(yl*yl+zl*zl))) return false;
	if(c.plane_intersects(x0,yl,zl,radius.cutoff(yl*yl+zl*zl))) return false;
	if(c.plane_intersects(x0,yh,zl,radius.cutoff(yl*yh+zl*zl))) return false;
	if(c.plane_intersects(x1,yh,zl,radius.cutoff(yl*yh+zl*zl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on an edge which points along the y
 * direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
 *                    block.
 * \param[in] (xl,zl) the relative x and z coordinates of the corner of the
 *                    block closest to the cell center.
 * \param[in] (xh,zh) the relative x and z coordinates of the corner of the
 *                    block furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::edge_y_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint zl,fpoint xh,fpoint y1,fpoint zh) {
	if(c.plane_intersects_guess(xl,y0,zh,radius.cutoff(xl*xl+zl*zh))) return false;
	if(c.plane_intersects(xl,y1,zh,radius.cutoff(xl*xl+zl*zh))) return false;
	if(c.plane_intersects(xl,y1,zl,radius.cutoff(xl*xl+zl*zl))) return false;
	if(c.plane_intersects(xl,y0,zl,radius.cutoff(xl*xl+zl*zl))) return false;
	if(c.plane_intersects(xh,y0,zl,radius.cutoff(xl*xh+zl*zl))) return false;
	if(c.plane_intersects(xh,y1,zl,radius.cutoff(xl*xh+zl*zl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on an edge which points along the z
 * direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the block.
 * \param[in] (xl,yl) the relative x and y coordinates of the corner of the
 *                    block closest to the cell center.
 * \param[in] (xh,yh) the relative x and y coordinates of the corner of the
 *                    block furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::edge_z_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint z0,fpoint xh,fpoint yh,fpoint z1) {
	if(c.plane_intersects_guess(xl,yh,z0,radius.cutoff(xl*xl+yl*yh))) return false;
	if(c.plane_intersects(xl,yh,z1,radius.cutoff(xl*xl+yl*yh))) return false;
	if(c.plane_intersects(xl,yl,z1,radius.cutoff(xl*xl+yl*yl))) return false;
	if(c.plane_intersects(xl,yl,z0,radius.cutoff(xl*xl+yl*yl))) return false;
	if(c.plane_intersects(xh,yl,z0,radius.cutoff(xl*xh+yl*yl))) return false;
	if(c.plane_intersects(xh,yl,z1,radius.cutoff(xl*xh+yl*yl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on a face aligned with the x direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] xl the minimum distance from the cell center to the face.
 * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
 *                    block.
 * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
 *                    block.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::face_x_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint z0,fpoint y1,fpoint z1) {
	if(c.plane_intersects_guess(xl,y0,z0,radius.cutoff(xl*xl))) return false;
	if(c.plane_intersects(xl,y0,z1,radius.cutoff(xl*xl))) return false;
	if(c.plane_intersects(xl,y1,z1,radius.cutoff(xl*xl))) return false;
	if(c.plane_intersects(xl,y1,z0,radius.cutoff(xl*xl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on a face aligned with the y direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] yl the minimum distance from the cell center to the face.
 * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
 *                    block.
 * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
 *                    block.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::face_y_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint z0,fpoint x1,fpoint z1) {
	if(c.plane_intersects_guess(x0,yl,z0,radius.cutoff(yl*yl))) return false;
	if(c.plane_intersects(x0,yl,z1,radius.cutoff(yl*yl))) return false;
	if(c.plane_intersects(x1,yl,z1,radius.cutoff(yl*yl))) return false;
	if(c.plane_intersects(x1,yl,z0,radius.cutoff(yl*yl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on a face aligned with the z direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] zl the minimum distance from the cell center to the face.
 * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
 *                    block.
 * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
 *                    block.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_base<r_option>::face_z_test(voronoicell_base<n_option> &c,fpoint x0,fpoint y0,fpoint zl,fpoint x1,fpoint y1) {
	if(c.plane_intersects_guess(x0,y0,zl,radius.cutoff(zl*zl))) return false;
	if(c.plane_intersects(x0,y1,zl,radius.cutoff(zl*zl))) return false;
	if(c.plane_intersects(x1,y1,zl,radius.cutoff(zl*zl))) return false;
	if(c.plane_intersects(x1,y0,zl,radius.cutoff(zl*zl))) return false;
	return true;
}

/** Creates a voropp_loop object, by setting the necessary constants about the
 * container geometry from a pointer to the current container class.
 * \param[in] q a pointer to the current container class. */
template<class r_option>
voropp_loop::voropp_loop(container_base<r_option> *q) : sx(q->bx-q->ax), sy(q->by-q->ay), sz(q->bz-q->az),
	xsp(q->xsp),ysp(q->ysp),zsp(q->zsp),
	ax(q->ax),ay(q->ay),az(q->az),
	nx(q->nx),ny(q->ny),nz(q->nz),nxy(q->nxy),nxyz(q->nxyz),
	xperiodic(q->xperiodic),yperiodic(q->yperiodic),zperiodic(q->zperiodic) {}

/** Initializes a voropp_loop object, by finding all blocks which are within a
 * given sphere. It calculates the index of the first block that needs to be
 * tested and sets the periodic displacement vector accordingly.
 * \param[in] (vx,vy,vz) the position vector of the center of the sphere.
 * \param[in] r the radius of the sphere.
 * \param[out] (px,py,pz) the periodic displacement vector for the first block
 *                        to be tested.
 * \return The index of the first block to be tested. */
inline int voropp_loop::init(fpoint vx,fpoint vy,fpoint vz,fpoint r,fpoint &px,fpoint &py,fpoint &pz) {
	ai=step_int((vx-ax-r)*xsp);
	bi=step_int((vx-ax+r)*xsp);
	if(!xperiodic) {
		if(ai<0) {ai=0;if(bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if(ai>=nx) ai=nx-1;}
	}
	aj=step_int((vy-ay-r)*ysp);
	bj=step_int((vy-ay+r)*ysp);
	if(!yperiodic) {
		if(aj<0) {aj=0;if(bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if(aj>=ny) aj=ny-1;}
	}
	ak=step_int((vz-az-r)*zsp);
	bk=step_int((vz-az+r)*zsp);
	if(!zperiodic) {
		if(ak<0) {ak=0;if(bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if(ak>=nz) ak=nz-1;}
	}
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	akp=kp=step_mod(k,nz);apz=pz=step_div(k,nz)*sz;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(ny+ajp-step_mod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
}

/** Initializes a voropp_loop object, by finding all blocks which overlap a given
 * rectangular box. It calculates the index of the first block that needs to be
 * tested and sets the periodic displacement vector (px,py,pz) accordingly.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box.
 * \param[out] (px,py,pz) the periodic displacement vector for the first block
 *                        to be tested.
 * \return The index of the first block to be tested. */
inline int voropp_loop::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax,fpoint &px,fpoint &py,fpoint &pz) {
	ai=step_int((xmin-ax)*xsp);
	bi=step_int((xmax-ax)*xsp);
	if(!xperiodic) {
		if(ai<0) {ai=0;if(bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if(ai>=nx) ai=nx-1;}
	}
	aj=step_int((ymin-ay)*ysp);
	bj=step_int((ymax-ay)*ysp);
	if(!yperiodic) {
		if(aj<0) {aj=0;if(bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if(aj>=ny) aj=ny-1;}
	}
	ak=step_int((zmin-az)*zsp);
	bk=step_int((zmax-az)*zsp);
	if(!zperiodic) {
		if(ak<0) {ak=0;if(bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if(ak>=nz) ak=nz-1;}
	}
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	akp=kp=step_mod(k,nz);apz=pz=step_div(k,nz)*sz;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(ny+ajp-step_mod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
}

/** Returns the next block to be tested in a loop, and updates the periodicity
 * vector if necessary.
 * \param[in,out] (px,py,pz) the current block on entering the function, which
 *                           is updated to the next block on exiting the
 *                           function. */
inline int voropp_loop::inc(fpoint &px,fpoint &py,fpoint &pz) {
	if(i<bi) {
		i++;
		if(ip<nx-1) {ip++;s++;} else {ip=0;s+=1-nx;px+=sx;}
		return s;
	} else if(j<bj) {
		i=ai;ip=aip;px=apx;j++;
		if(jp<ny-1) {jp++;s+=inc1;} else {jp=0;s+=inc1-nxy;py+=sy;}
		return s;
	} else if(k<bk) {
		i=ai;ip=aip;j=aj;jp=ajp;px=apx;py=apy;k++;
		if(kp<nz-1) {kp++;s+=inc2;} else {kp=0;s+=inc2-nxyz;pz+=sz;}
		return s;
	} else return -1;
}

/** Custom int function, that gives consistent stepping for negative numbers.
 * With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
 * With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1). */
inline int voropp_loop::step_int(fpoint a) {
	return a<0?int(a)-1:int(a);
}

/** Custom modulo function, that gives consistent stepping for negative
 * numbers. */
inline int voropp_loop::step_mod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
}

/** Custom integer division function, that gives consistent stepping for
 * negative numbers. */
inline int voropp_loop::step_div(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
}

/** Adds a wall to the container.
 * \param[in] w a wall object to be added.*/
template<class r_option>
void container_base<r_option>::add_wall(wall& w) {
	if(wall_number==current_wall_size) {
		current_wall_size*=2;
		if(current_wall_size>max_wall_size)
			voropp_fatal_error("Wall memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
		wall **pwall;
		pwall=new wall*[current_wall_size];
		for(int i=0;i<wall_number;i++) pwall[i]=walls[i];
		delete [] walls;walls=pwall;
	}
	walls[wall_number++]=&w;
}

/** Sets the radius of the jth particle in region i to r, and updates the
 * maximum particle radius.
 * \param[in] i the region of the particle to consider.
 * \param[in] j the number of the particle within the region.
 * \param[in] r the radius to set. */
inline void radius_poly::store_radius(int i,int j,fpoint r) {
	cc->p[i][4*j+3]=r;
	if(r>max_radius) max_radius=r;
}

/** Clears the stored maximum radius. */
inline void radius_poly::clear_max() {
	max_radius=0;
}

/** Imports a list of particles from an input stream for the monodisperse case
 * where no radius information is expected.
 * \param[in] is an input stream to read from. */
inline void radius_mono::import(istream &is) {
	int n;fpoint x,y,z;
	is >> n >> x >> y >> z;
	while(!is.eof()) {
		cc->put(n,x,y,z);
		is >> n >> x >> y >> z;
	}
}

/** Imports a list of particles from an input stream for the polydisperse case,
 * where both positions and particle radii are both stored.
 * \param[in] is an input stream to read from. */
inline void radius_poly::import(istream &is) {
	int n;fpoint x,y,z,r;
	is >> n >> x >> y >> z >> r;
	while(!is.eof()) {
		cc->put(n,x,y,z,r);
		is >> n >> x >> y >> z >> r;
	}
}

/** Initializes the radius_poly class for a new Voronoi cell calculation, by
 * computing the radial cut-off value, based on the current particle's radius
 * and the maximum radius of any particle in the packing.
 * \param[in] ijk the region to consider.
 * \param[in] s the number of the particle within the region. */
inline void radius_poly::init(int ijk,int s) {
	crad=cc->p[ijk][4*s+3];
	mul=1+(crad*crad-max_radius*max_radius)/((max_radius+crad)*(max_radius+crad));
	crad*=crad;
}

/** This routine is called when deciding when to terminate the computation of a
 * Voronoi cell. For the Voronoi radical tessellation for a polydisperse case,
 * this routine multiplies the cutoff value by the scaling factor that was
 * precomputed in the init() routine.
 * \param[in] lrs a cutoff radius for the cell computation.
 * \return The value scaled by the factor mul. */
inline fpoint radius_poly::cutoff(fpoint lrs) {
	return mul*lrs;
}

/** This routine is called when deciding when to terminate the computation of a
 * Voronoi cell. For the monodisperse case, this routine just returns the same
 * value that is passed to it.
 * \param[in] lrs a cutoff radius for the cell computation.
 * \return The same value passed to it. */
inline fpoint radius_mono::cutoff(fpoint lrs) {
	return lrs;
}

/** Prints the radius of particle, by just supplying a generic value of "s".
 * \param[in] os the output stream to write to.
 * \param[in] l the region to consider.
 * \param[in] c the number of the particle within the region. */
inline void radius_mono::rad(ostream &os,int l,int c) {
	os << "s";
}

/** Prints the radius of a particle to an open output stream.
 * \param[in] os the output stream to write to.
 * \param[in] l the region to consider.
 * \param[in] c the number of the particle within the region. */
inline void radius_poly::rad(ostream &os,int l,int c) {
	os << cc->p[l][4*c+3];
}

/** Returns the scaled volume of a particle, which is always set to 0.125 for
 * the monodisperse case where particles are taken to have unit diameter.
 * \param[in] ijk the region to consider.
 * \param[in] s the number of the particle within the region.
 * \return The cube of the radius of the particle, which is 0.125 in this case.
 */
inline fpoint radius_mono::volume(int ijk,int s) {
	return 0.125;
}

/** Returns the scaled volume of a particle.
 * \param[in] ijk the region to consider.
 * \param[in] s the number of the particle within the region.
 * \return The cube of the radius of the particle. */
inline fpoint radius_poly::volume(int ijk,int s) {
	fpoint a=cc->p[ijk][4*s+3];
	return a*a*a;
}

/** Scales the position of a plane according to the relative sizes
 * of the particle radii.
 * \param[in] rs the distance between the Voronoi cell and the cutting plane.
 * \param[in] t the region to consider
 * \param[in] q the number of the particle within the region.
 * \return The scaled position. */
inline fpoint radius_poly::scale(fpoint rs,int t,int q) {
	return rs+crad-cc->p[t][4*q+3]*cc->p[t][4*q+3];
}

/** Applies a blank scaling to the position of a cutting plane.
 * \param[in] rs the distance between the Voronoi cell and the cutting plane.
 * \param[in] t the region to consider
 * \param[in] q the number of the particle within the region.
 * \return The scaled position, which for this case, is equal to rs. */
inline fpoint radius_mono::scale(fpoint rs,int t,int q) {
	return rs;
}

/** Prints the radius of a particle to an open file stream.
 * \param[in] os an open file stream.
 * \param[in] ijk the region to consider.
 * \param[in] q the number of the particle within the region.
 * \param[in] later A boolean value to determine whether or not to write a
 *                  space character before the first entry. */
inline void radius_poly::print(ostream &os,int ijk,int q,bool later) {
	if(later) os << " ";
	os << cc->p[ijk][4*q+3];
}

/** This function is called during container construction. The routine scans
 * all of the worklists in the wl[] array. For a given worklist of blocks
 * labeled \f$w_1\f$ to \f$w_n\f$, it computes a sequence \f$r_0\f$ to
 * \f$r_n\f$ so that $r_i$ is the minimum distance to all the blocks
 * \f$w_{j}\f$ where \f$j>i\f$ and all blocks outside the worklist. The values
 * of \f$r_n\f$ is calculated first, as the minimum distance to any block in
 * the shell surrounding the worklist. The \f$r_i\f$ are then computed in
 * reverse order by considering the distance to \f$w_{i+1}\f$. */
template<class r_option>
void container_base<r_option>::initialize_radii() {
	const unsigned int b1=1<<21,b2=1<<22,b3=1<<24,b4=1<<25,b5=1<<27,b6=1<<28;
	const fpoint xstep=(bx-ax)/nx/fgrid;
	const fpoint ystep=(by-ay)/ny/fgrid;
	const fpoint zstep=(bz-az)/nz/fgrid;
	int i,j,k,lx,ly,lz,l=0,q;
	unsigned int *e=const_cast<unsigned int*> (wl);
	fpoint xlo,ylo,zlo,xhi,yhi,zhi,minr;fpoint *radp=mrad;
	unsigned int f;
	for(zlo=0,zhi=zstep,lz=0;lz<hgrid;zlo=zhi,zhi+=zstep,lz++) {
		for(ylo=0,yhi=ystep,ly=0;ly<hgrid;ylo=yhi,yhi+=ystep,ly++) {
			for(xlo=0,xhi=xstep,lx=0;lx<hgrid;xlo=xhi,xhi+=xstep,l++,lx++) {
				minr=large_number;
				for(q=e[0]+1;q<seq_length;q++) {
					f=e[q];
					i=(f&127)-64;
					j=(f>>7&127)-64;
					k=(f>>14&127)-64;
					if((f&b2)==b2) {
						compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i-1,j,k);
						if((f&b1)==0) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i+1,j,k);
					} else if((f&b1)==b1) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i+1,j,k);
					if((f&b4)==b4) {
						compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j-1,k);
						if((f&b3)==0) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j+1,k);
					} else if((f&b3)==b3) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j+1,k);
					if((f&b6)==b6) {
						compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k-1);
						if((f&b5)==0) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k+1);
					} else if((f&b5)==b5) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k+1);
				}
				q--;
				while(q>0) {
					radp[q]=minr;
					f=e[q];
					i=(f&127)-64;
					j=(f>>7&127)-64;
					k=(f>>14&127)-64;
					compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k);
					q--;
				}
				radp[0]=minr;
				e+=seq_length;
				radp+=seq_length;
			}
		}
	}
}

/** Computes the minimum distance from a subregion to a given block. If this distance
 * is smaller than the value of minr, then it passes
 * \param[in,out] minr a pointer to the current minimum distance. If the distance
 *                     computed in this function is smaller, then this distance is
 *                     set to the new one.
 * \param[out] (xlo,ylo,zlo) the lower coordinates of the subregion being
 *                           considered.
 * \param[out] (xhi,yhi,zhi) the upper coordinates of the subregion being
 *                           considered.
 * \param[in] (ti,tj,tk) the coordinates of the block. */
template<class r_option>
inline void container_base<r_option>::compute_minimum(fpoint &minr,fpoint &xlo,fpoint &xhi,fpoint &ylo,fpoint &yhi,fpoint &zlo,fpoint &zhi,int ti,int tj,int tk) {
	const fpoint boxx=(bx-ax)/nx,boxy=(by-ay)/ny,boxz=(bz-az)/nz;
	fpoint radsq,temp;
	if(ti>0) {temp=boxx*ti-xhi;radsq=temp*temp;}
	else if(ti<0) {temp=xlo-boxx*(1+ti);radsq=temp*temp;}
	else radsq=0;

	if(tj>0) {temp=boxy*tj-yhi;radsq+=temp*temp;}
	else if(tj<0) {temp=ylo-boxy*(1+tj);radsq+=temp*temp;}

	if(tk>0) {temp=boxz*tk-zhi;radsq+=temp*temp;}
	else if(tk<0) {temp=zlo-boxz*(1+tk);radsq+=temp*temp;}

	if(radsq<minr) minr=radsq;
}

/** This routine checks to see whether a point is within a particular distance
 * of a nearby region. If the point is within the distance of the region, then
 * the routine returns true, and computes the maximum distance from the point
 * to the region. Otherwise, the routine returns false.
 * \param[in] (di,dj,dk) the position of the nearby region to be tested,
 *                       relative to the region that the point is in.
 * \param[in] (fx,fy,fz) the displacement of the point within its region.
 * \param[in] (gxs,gys,gzs) the minimum squared distances from the point to the
 *                          sides of its region.
 * \param[out] crs a reference in which to return the maximum distance to the
 *                 region (only computed if the routine returns positive).
 * \param[in] mrs the distance to be tested.
 * \return False if the region is further away than mrs, true if the region in
 *         within mrs.*/
template<class r_option>
inline bool container_base<r_option>::compute_min_max_radius(int di,int dj,int dk,fpoint fx,fpoint fy,fpoint fz,fpoint gxs,fpoint gys,fpoint gzs,fpoint &crs,fpoint mrs) {
	fpoint xlo,ylo,zlo;
	const fpoint boxx=(bx-ax)/nx,boxy=(by-ay)/ny,boxz=(bz-az)/nz;
	const fpoint bxsq=boxx*boxx+boxy*boxy+boxz*boxz;
	if(di>0) {
		xlo=di*boxx-fx;
		crs=xlo*xlo;
		if(dj>0) {
			ylo=dj*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo+boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo+boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(2*xlo+boxx)+boxy*(2*ylo+boxy)+gzs;
			}
		} else if(dj<0) {
			ylo=(dj+1)*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo-boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo-boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(2*xlo+boxx)+boxy*(-2*ylo+boxy)+gzs;
			}
		} else {
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=gys+boxx*(2*xlo+boxx);
		}
	} else if(di<0) {
		xlo=(di+1)*boxx-fx;
		crs=xlo*xlo;
		if(dj>0) {
			ylo=dj*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo+boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo+boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(-2*xlo+boxx)+boxy*(2*ylo+boxy)+gzs;
			}
		} else if(dj<0) {
			ylo=(dj+1)*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo-boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo-boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(-2*xlo+boxx)+boxy*(-2*ylo+boxy)+gzs;
			}
		} else {
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=gys+boxx*(-2*xlo+boxx);
		}
	} else {
		if(dj>0) {
			ylo=dj*boxy-fy;
			crs=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=boxy*(2*ylo+boxy);
		} else if(dj<0) {
			ylo=(dj+1)*boxy-fy;
			crs=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=boxy*(-2*ylo+boxy);
		} else {
			if(dk>0) {
				zlo=dk*boxz-fz;crs=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;crs=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				crs=0;
				voropp_fatal_error("Min/max radius function called for central block, which should never\nhappen.",VOROPP_INTERNAL_ERROR);
			}
			crs+=gys;
		}
		crs+=gxs;
	}
	return false;
}

#include "worklist.cc"
