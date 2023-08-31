//NOTE: LBM
//D2Q9i model
//The top, bottom, inlet and outlet are all ZOU HE BC with slip/no-slip
//velocity inlet and pressure outlet with uniform flow all over
//inline oscillation
//height 20 D
//length 30 D

//NOTE: IBM
//srping mass system with Guo's explicit forcing scheme
//a four point stencil is used for the kernel function

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<omp.h>

#define num_threads 8

////////////////////////////// global variable //////////////////////////////////////////

//~~~~~~~~~~Eulerian Varaibles~~~~~~~~~~~
int i, j, k, p, m1=1, m2=201, n1=1, n2=3501, alpha=8, nihal=0, count=0;
double omega=1.0/0.6364, U=0.1;
double V_in=0.0, rho_in=0.0, rho_out=1.0;	//Re40
double t=0.0, dx=0.0, dy=0.0;

//~~~~~~~~~Lagrangian Varaibles~~~~~~~~~~
int n, P=1;
double Cx_ref[2] = { 101.0 , 101.0 };
double Cy_ref[2] = { 2251.0 , 2251.0 };
double Cx[2], Cy[2];
double axis[2] = { 25.0 , 12.5 }; 	// half length of major axis
double axis_ratio = 2.0;	// major axis/minor axis
double surface_area[2] = { 121.10560270162605 , 121.10560270162605 };	//2*M_PI*r for circle; perimeter for ellipse
double vol[2] = { 981.74770424681037 , 981.747704246810387 };
double mom_inertia[2] = { 195.3125 , 195.3125 };	// here only r^2 is there. (for ellipse (a2+b2)/4)
int N[2] = { 140 , 140 };	//N>[(perimeter of geometry)]
double psi[2] = { M_PI/4.0 , M_PI/4.0 };

double F_x[2], F_y[2], Tor[2], Uc[2], Vc[2], V_ang[2], grav[2];
double F_tot_x[2], F_tot_y[2], Tor_tot[2];
double Uc_old[2], Vc_old[2], V_ang_old[2];
double int_fx[2], int_fy[2], int_tor[2];
double F_pw[2], F_pp_x[2], F_pp_y[2];
double C_d[2], C_l[2];
double del_s[2], rho_avg[2];

double rho_p_by_f[2] = { 1.1 , 1.1 };
double gravity = 0.000162073;

//~~~~~~~~~~~~~~extras~~~~~~~~~~~~~~~~~~~
double sum=0.0, sumu=0.0, sumv=0.0;
double temp1=0.0, temp2=0.0, temp3=0.0, temp4=0.0, temp5=0.0;
double a=0.0, b=0.0;
double phi=0.0, theta=0.0;            
double error=0.01, e=1.0e-10;


////////////////////////////// data structure ///////////////////////////////////////////
struct value1{
	int e_x, e_y;
	int k_opp;
}vector[10];

struct value_1{
	double w;
}weight[10];

struct value2{
	double rho, u, v, vor, u_old, v_old;
	double fn_x, fn_y, mom_x, mom_y;
}node[500][4101];

struct value5{
	double flag;
}domain[500][4101][3];

struct value3{
	double f, f_0, f_tilde;
	double force;
}lattice[500][4101][10];

struct value4{
	double x, y, angle, x_ref, y_ref, fp_x, fp_y, u, v, rho, mom_x, mom_y;
	double r, V_tang;
}particle[3][200];

///////////////////////// defining velocity vectors /////////////////////////////////////
void velocity_vectors()
{
	vector[0].e_x = 0;	vector[0].e_y = 0;
	vector[1].e_x = 1;	vector[1].e_y = 0;
	vector[2].e_x = 0;	vector[2].e_y = 1;
	vector[3].e_x = -1;	vector[3].e_y = 0;
	vector[4].e_x = 0;	vector[4].e_y = -1;
	vector[5].e_x = 1;	vector[5].e_y = 1;
	vector[6].e_x = -1;	vector[6].e_y = 1;
	vector[7].e_x = -1;	vector[7].e_y = -1;
	vector[8].e_x = 1; 	vector[8].e_y = -1;
	
	vector[0].k_opp = 0;
	vector[1].k_opp = 3;
	vector[2].k_opp = 4;
	vector[3].k_opp = 1;
	vector[4].k_opp = 2;
	vector[5].k_opp = 7;
	vector[6].k_opp = 8;
	vector[7].k_opp = 5;
	vector[8].k_opp = 6;
}

////////////////////////// assigning weighting factors //////////////////////////////////
void weighting_factors()
{
	weight[0].w = 4.0/9.0;
	for(k=1;k<=4;++k)
		weight[k].w = 1.0/9.0;
	for(k=5;k<=8;++k)
		weight[k].w = 1.0/36.0;
}

/////////////////////// functions ///////////////////////////////////////////////////////
void grid_gen();
void body();
void initialising();
void initialising_pdf();
void lattice_force();
void collision();
void streaming();
void bounce_back();
void rho_momentum();
void interpolate_to_lag_nodes();
void compute_lag_node_force();
void spread_force_to_eul_nodes();
void u_v_calculator();
void force_coefficients();
void update_lag_node_position();
void printing_everything();
void time_files();

/////////////////////// assigning the domain ////////////////////////////////////////////
void grid_gen()
{
	FILE *dom;
	dom = fopen("domain.dat", "w");
	fprintf(dom, "ZONE\tI=%d\tJ=%d\n", (m2-m1+1), (n2-n1+1));
	for(j=n1;j<=n2;++j)
	for(i=m1;i<=m2;++i)
	{
		fprintf(dom, "%d\t%d\t\n", i, j);
	}
	fclose(dom);
}

////////////////////////// constructing body/////////////////////////////////////////////
void body()
{
	for(p=0;p<P;++p)
	{	
		for(n=1;n<=N[p];++n)
		{
			//calculate the angular orientation of the Lagrangian particles
			particle[p][n].angle = 2.0*M_PI*((double) n)/N[p];
			
			if(particle[p][n].angle>=(2.0*M_PI))
			{	particle[p][n].angle = particle[p][n].angle - 2.0*M_PI;	}
			else if(particle[p][n].angle<0.0)
			{	particle[p][n].angle = particle[p][n].angle + 2.0*M_PI;	}
			
		//	particle[p][n].x_ref = Cx_ref[p] + axis[p]*cos(particle[p][n].angle - psi[p])*cos(psi[p]) - (axis[p]/axis_ratio)*sin(particle[p][n].angle - psi[p])*sin(psi[p]);
		//	particle[p][n].y_ref = Cy_ref[p] + axis[p]*cos(particle[p][n].angle - psi[p])*sin(psi[p]) + (axis[p]/axis_ratio)*sin(particle[p][n].angle - psi[p])*cos(psi[p]);
			
			//calculate the position of the lagrangian nodes
			particle[p][n].x_ref = Cx_ref[p] + axis[p]*cos(particle[p][n].angle);
			particle[p][n].y_ref = Cy_ref[p] + (axis[p]/axis_ratio)*sin(particle[p][n].angle);
		
		//	particle[p][n].angle = 2.0*M_PI*((double) n)/N[p] + (psi[p]);
			
			//calculate the radial distance from the geometric centre
			particle[p][n].r = sqrt( pow((particle[p][n].x_ref - Cx_ref[p]),2.0) + pow((particle[p][n].y_ref - Cy_ref[p]),2.0) );
			a = particle[p][n].x_ref - Cx_ref[p];
			b = particle[p][n].y_ref - Cy_ref[p];

			//calculate (again) the angular orientation of the Lagrangian nodes from the geometric centre
			if(a==0.0000000000000000 && b==0.0000000000000000)
			{	a = 0.0000000000000001;
				b = 0.0000000000000001;
				phi = fabs(atan(fabs(b/a)));	}
			else
			{	
				phi = fabs(atan(fabs(b/a)));
				if(a>0.0)
				{
					if(b<0.0)
					{	phi = 2.0*M_PI - phi;	}
					else if(b==0.0)
					{	phi = 0.0;	}
					else if(b>0.0)
					{	phi = phi;	}
				}
				else if(a<0.0)
				{
					if(b>0.0)
					{	phi = M_PI - phi;	}
					else if(b<0.0)
					{	phi = M_PI + phi;	}
					else if(b==0.0)
					{	phi = M_PI;	}
				}
				else if(a==0.0)
				{
					if(b>0.0)
					{	phi = M_PI/2.0;	}
					else if(b<0.0)
					{	phi = 3.0*M_PI/2.0;	}
					else if(b==0.0)
					{	phi = 0.0000000000000001;	}
				}
			}
				
			particle[p][n].angle = phi;
		}
		
		//change the orientation of the body
		for(n=1;n<=N[p];++n)
		{
			particle[p][n].angle = particle[p][n].angle + psi[p];
			
			if(particle[p][n].angle>=(2.0*M_PI))
			{	particle[p][n].angle = particle[p][n].angle - 2.0*M_PI;	}
			else if(particle[p][n].angle<0.0)
			{	particle[p][n].angle = particle[p][n].angle + 2.0*M_PI;	}
			
			particle[p][n].x_ref = Cx_ref[p] + particle[p][n].r*cos(particle[p][n].angle);
			particle[p][n].y_ref = Cy_ref[p] + particle[p][n].r*sin(particle[p][n].angle);
			
		}
	}
	
	//determining which nodes are interior at each time step
	for(p=0;p<P;++p)
	{
		temp1=0.0;
		for(i=m1;i<=m2;++i)
		for(j=n1;j<=n2;++j)
		{
			temp1 = (pow(((i-Cx_ref[p])*cos(psi[p])+(j-Cy_ref[p])*sin(psi[p])),2.0))/(pow(axis[p],2.0)) + (pow(((i-Cx_ref[p])*sin(psi[p])-(j-Cy_ref[p])*cos(psi[p])),2.0))/(pow((axis[p]/axis_ratio),2.0));
			if(temp1<=1.0)
				domain[i][j][p].flag = 1.0;
			else
				domain[i][j][p].flag = 0.0;
		}
	}
	
	p=0;
	FILE *bdy1;
	bdy1 = fopen("particle_initial1.dat", "w");
	for(n=1;n<=N[p];++n)
		fprintf(bdy1, "%lf\t%lf\n", particle[p][n].x_ref, particle[p][n].y_ref);
	fclose(bdy1);
	p=1;
	FILE *bdy2;
	bdy2 = fopen("particle_initial2.dat", "w");
	for(n=1;n<=N[p];++n)
		fprintf(bdy2, "%lf\t%lf\n", particle[p][n].x_ref, particle[p][n].y_ref);
	fclose(bdy2);
	
	FILE *dom;
	dom = fopen("domain.dat", "w");
	fprintf(dom, "ZONE\tI=%d\tJ=%d\n", (m2-m1+1), (n2-n1+1));
	for(j=n1;j<=n2;++j)
	for(i=m1;i<=m2;++i)
	{
		fprintf(dom, "%d\t%d\t%f\n", i, j, domain[i][j][0].flag);
	}
	fclose(dom);
}

/////////////////////////////// initialising ////////////////////////////////////////////
void initialising()
{
	for(i=m1-1;i<=m2+1;++i)
	for(j=n1-1;j<=n2+1;++j)
	{
		node[i][j].rho = 1.0;
		node[i][j].u = 0.0;
		node[i][j].v = 0.0;
		node[i][j].vor = 0.0;
		node[i][j].mom_x = 0.0;
		node[i][j].mom_y = 0.0;
		node[i][j].fn_x = 0.0;
		node[i][j].fn_y = 0.0;
		for(k=0;k<=alpha;++k)
		{
			lattice[i][j][k].f = 0.0;
			lattice[i][j][k].f_tilde = 0.0;
			lattice[i][j][k].f_0 = 0.0;
		}
	}
	
	for(p=0;p<P;++p)
	{
		Cx[p] = Cx_ref[p];
		Cy[p] = Cy_ref[p];
		
		
		
		
		
		
		Uc[p]=0.0;
		Vc[p]=0.0;//(1.0/rho_p_by_f-1.0)*gravity;
		F_pw[p]=0.0;
		F_pp_x[p]=0.0;
		F_pp_y[p]=0.0;
		for(n=1;n<=N[p];++n)
		{
			particle[p][n].x = particle[p][n].x_ref;
			particle[p][n].y = particle[p][n].y_ref;
			particle[p][n].u = Uc[p];
			particle[p][n].v = Vc[p];
			particle[p][n].rho = 1.0;
			particle[p][n].mom_x = 0.0;
			particle[p][n].mom_y = 0.0;
			particle[p][n].fp_x = 0.0;
			particle[p][n].fp_y = 0.0;
		}
	}
	
//	FILE *para;
//	para = fopen("parabola.dat", "w");
//	for(j=n1;j<=n2;++j)
//		fprintf(para, "%lf\t%lf\n", dy*(j-n1), node[m1][j].u/U);
//	fclose(para);


/*	//....scanning....
	int points=0, p=0;

	double a=0.0, b=0.0, c=0.0, d=0.0;
	FILE *in;
	in = fopen("input1.dat", "r");
	points = m2*n2;
	printf("%d\n", points);
	for(p=1;p<=points;++p)
	{
		fscanf(in, "%d	%d	%lf	%lf	%lf	%lf\n", &i, &j, &a, &b, &c, &d);
		node[i][j].u = a;
		node[i][j].v = b;
		node[i][j].rho = c;

		node[i][j].vor = d;
	}
	fclose(in);
	FILE *in2;
	in2 = fopen("input2.dat", "r");
	points=m2*n2*(alpha+1);
	printf("%d\n", points);
	for(p=1;p<=points;++p)
	{
		fscanf(in2, "%d	%d	%d	%lf	%lf\n", &i, &j, &k, &a, &b);
		lattice[i][j][k].f_tilde = a;

		lattice[i][j][k].f_0 = b;
	}
	fclose(in2);
*/	
}

////////////////////////////initialising pdfs ///////////////////////////////////////////
void initialising_pdf()
{
	for(i=m1;i<=m2;++i)
	for(j=n1;j<=n2;++j)
	{
		temp3 = pow((node[i][j].u),2.0) + pow((node[i][j].v),2.0);
		for(k=0;k<=alpha;++k)
		{
			temp1 = vector[k].e_x*node[i][j].u + vector[k].e_y*node[i][j].v;
			lattice[i][j][k].f_0 = weight[k].w*( node[i][j].rho + 3.0*temp1 + (9.0/2.0)*temp1*temp1 - (3.0/2.0)*temp3 );
			lattice[i][j][k].f = lattice[i][j][k].f_0;
			lattice[i][j][k].f_tilde = lattice[i][j][k].f_0;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN BODY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
/////////////////////////////////////////////////////////////////////////////////////////
void main()
{
	velocity_vectors();
	weighting_factors();
	grid_gen();
	body();
	initialising();
	initialising_pdf();	//initialising the pdfs (f=f_eq)
	
	FILE *ptr3;
	ptr3 = fopen("error.dat", "w");
	
	FILE *part1, *part2;
	part1 = fopen("body1.dat", "w");
	part2 = fopen("body2.dat", "w");

//	struct timeval tv1, tv2;
//	gettimeofday(&tv1, NULL);

	//time loop starts
	do{
		t=t+1.0;

		//saving old data (u and v)
		omp_set_num_threads(num_threads);
		#pragma omp parallel
		{int i;
		#pragma omp for
		for(i=m1;i<=m2;++i)
		{int j;
		for(j=n1;j<=n2;++j)
		{   node[i][j].u_old = node[i][j].u;
			node[i][j].v_old = node[i][j].v;	}}}

			//lattice force computation
			lattice_force();
		
			//collision step
			collision();

			//streaming step
			streaming();

			//bounce back
			bounce_back();
			
			//calculating desnity
			rho_momentum();
			
			//interpolation velocity to lagrangian nodes
			interpolate_to_lag_nodes();
		
			//calculate force on lagrangian particles
			compute_lag_node_force();
			
			//spread the computed force to Eulerian nodes
			spread_force_to_eul_nodes();
		
			//calculating macroscopic variables
			u_v_calculator();
			
			//drag coefficient
			force_coefficients();
			
			// update the lagrangian particle position
			update_lag_node_position();
			
		p=0;
		fprintf(part1, "%1.1lf\t%1.9lf\t%1.9lf\t%1.9lf\t%1.9lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%1.9lf\n", t, Cx[p], Cy[p], Uc[p], Vc[p], sqrt(Uc[p]*Uc[p]+Vc[p]*Vc[p]), V_ang[p], psi[p], F_x[p], F_y[p], Tor[p], F_pw[p]);
		p=1;
		fprintf(part2, "%1.1lf\t%1.9lf\t%1.9lf\t%1.9lf\t%1.9lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%1.9lf\n", t, Cx[p], Cy[p], Uc[p], Vc[p], sqrt(Uc[p]*Uc[p]+Vc[p]*Vc[p]), V_ang[p], F_x[p], F_y[p], Tor[p], F_pw[p]);
		
		FILE *ptr9;
		ptr9 = fopen("Cd_Cl.dat","a+");	
		fprintf(ptr9, "%1.1lf\t%1.9lf\t%1.9lf\t%1.9lf\t%1.9lf\n", t, C_d[0], C_l[0], C_d[1], C_l[1]);
		fclose(ptr9);
		
		if((((int) (t))%500)==0)
		{	
			
			//.....calculating error......
			sumu=0.0; sumv=0.0; sum=0.0; error=0.0;
			omp_set_num_threads(num_threads);
			#pragma omp parallel
			{
				int i;
				#pragma omp for reduction (+:sumu,sumv)
				for(i=m1;i<=m2;++i)
				{int j;
				for(j=n1;j<=n2;++j)
				{
						sumu = sumu + pow((node[i][j].u-node[i][j].u_old),2.0);
						sumv = sumv + pow((node[i][j].v-node[i][j].v_old),2.0);
				}}
				sum=(sumu+sumv)/2.0;
			}
			error = sqrt(sum/((m2-m1+1)*(n2-n1+1)));
			fprintf(ptr3, "%0.1lf\t%1.16lf\n", t, error);
			printf("%0.1lf\t%1.16lf\n", t, error);
		
			//.........printing everything............
			printing_everything();
	
		//	time_files();
		
			FILE *dom;
			dom = fopen("domain.dat", "w");
			fprintf(dom, "ZONE\tI=%d\tJ=%d\n", (m2-m1+1), (n2-n1+1));
			for(j=n1;j<=n2;++j)
			for(i=m1;i<=m2;++i)
			{
				fprintf(dom, "%d\t%d\t%f\n", i, j, domain[i][j][0].flag);
			}
			fclose(dom);
			
			p=0;
			FILE *bdy2;
			bdy2 = fopen("particle.dat", "w");
			for(n=1;n<=N[p];++n)
				fprintf(bdy2, "%lf\t%lf\n", particle[p][n].x, particle[p][n].y);
			fclose(bdy2);
		}
		
	}while(Cy[0]>(axis[0]+0.1));//(t<10000.0);//(fabs(error)>e);//
	//time loop ends

//	gettimeofday(&tv2, NULL);
//	float sec = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);

//	FILE *time;
//	time = fopen("time.dat", "w");
//	fprintf(time, "time = %1.6f seconds\n", sec);
//	fclose(time);
	
	//.........printing everything............
	printing_everything();
	
/*	FILE *bdy3;
	bdy3 = fopen("particle.dat", "w");
	for(n=1;n<=N;++n)
		fprintf(bdy3, "%1.9lf\t%1.9lf\t%1.9lf\t%1.9lf\t%1.9lf\n", particle[n].x, particle[n].y, particle[n].u, particle[n].v, particle[n].angle/M_PI);
	fclose(bdy3);
*/	
	fclose(ptr3);
	fclose(part1);
	fclose(part2);
}

////////////////////// lattice force computation ////////////////////////////////////////
void lattice_force()
{
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
		int i;
		#pragma omp for
		for(i=m1;i<=m2;++i)
		{int j, k;
		for(j=n1;j<=n2;++j)
		{for(k=0;k<=alpha;++k)
		{
			lattice[i][j][k].force = (1.0-0.5*omega)*weight[k].w*( 3.0*(( vector[k].e_x-node[i][j].u )*node[i][j].fn_x + ( vector[k].e_y-node[i][j].v )*node[i][j].fn_y ) + 9.0*( (vector[k].e_x*node[i][j].u+vector[k].e_y*node[i][j].v)*(vector[k].e_x*node[i][j].fn_x+vector[k].e_y*node[i][j].fn_y) ) );
		}}
		}
	}
}

/////////////////////////////// collision ///////////////////////////////////////////////
void collision()
{
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
		int i;
		double temp1=0.0,temp3=0.0;
		#pragma omp for
		for(i=m1;i<=m2;++i)
		{int j, k;
		for(j=n1;j<=n2;++j)
		{
			temp3 = pow((node[i][j].u),2.0) + pow((node[i][j].v),2.0);
			for(k=0;k<=alpha;++k)
			{
				temp1 = vector[k].e_x*node[i][j].u + vector[k].e_y*node[i][j].v;
				lattice[i][j][k].f_0 = weight[k].w*( node[i][j].rho + 3.0*temp1 + (9.0/2.0)*temp1*temp1 - (3.0/2.0)*temp3 );
				lattice[i][j][k].f = omega*lattice[i][j][k].f_0 + (1.0-omega)*lattice[i][j][k].f_tilde + lattice[i][j][k].force;
			}
		}}
	}
}

///////////////////////////// streaming step ////////////////////////////////////////////
void streaming()
{
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
		int i;
		#pragma omp for
		for(i=m1;i<=m2;++i)
		{int j, k;
		for(j=n1;j<=n2;++j)
		{for(k=0;k<=alpha;++k)
			lattice[i][j][k].f_tilde = lattice[i-vector[k].e_x][j-vector[k].e_y][k].f;}}
	}
}

/////////////////////// begin bounce back ///////////////////////////////////////////////
void bounce_back()
{
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
	int i;
	#pragma omp for
	for(i=m1+1;i<=m2-1;++i)
	{
		//bottom wall (2,5,6)
	//	node[i][n1].u = node[i][n1+1].u;
	//	node[i][n1].u = (4.0*node[i][n1+1].u-node[i][n1+2].u)/3.0;
	//	node[i][n1].v = node[i][n1+1].v;
		node[i][n1].u = 0.0;
		node[i][n1].v = 0.0;
		node[i][n1].rho = node[i][n1].v+(lattice[i][n1][0].f_tilde+lattice[i][n1][1].f_tilde+lattice[i][n1][3].f_tilde+2.0*(lattice[i][n1][4].f_tilde+lattice[i][n1][7].f_tilde+lattice[i][n1][8].f_tilde));
		lattice[i][n1][2].f_tilde = lattice[i][n1][4].f_tilde+(2.0/3.0)*node[i][n1].v;
		lattice[i][n1][5].f_tilde = lattice[i][n1][7].f_tilde-0.5*(lattice[i][n1][1].f_tilde-lattice[i][n1][3].f_tilde)+(1.0/6.0)*node[i][n1].v+(1.0/2.0)*node[i][n1].u;
		lattice[i][n1][6].f_tilde = lattice[i][n1][8].f_tilde+0.5*(lattice[i][n1][1].f_tilde-lattice[i][n1][3].f_tilde)+(1.0/6.0)*node[i][n1].v-(1.0/2.0)*node[i][n1].u;
		
	//	lattice[i][n1][2].f_tilde = lattice[i][n2-1][2].f;
	//	lattice[i][n1][5].f_tilde = lattice[i-1][n2-1][5].f;
	//	lattice[i][n1][6].f_tilde = lattice[i+1][n2-1][6].f;
		
		//top wall (4,7,8)
		node[i][n2].u = 0.0;
	//	node[i][n2].u = (4.0*node[i][n2-1].u-node[i][n2-2].u)/3.0;
		node[i][n2].v = 0.0;
		node[i][n2].rho = -node[i][n2].v+(lattice[i][n2][0].f_tilde+lattice[i][n2][1].f_tilde+lattice[i][n2][3].f_tilde+2.0*(lattice[i][n2][2].f_tilde+lattice[i][n2][6].f_tilde+lattice[i][n2][5].f_tilde));
		lattice[i][n2][4].f_tilde = lattice[i][n2][2].f_tilde-(2.0/3.0)*node[i][n2].v;
		lattice[i][n2][7].f_tilde = lattice[i][n2][5].f_tilde+0.5*(lattice[i][n2][1].f_tilde-lattice[i][n2][3].f_tilde)-(1.0/6.0)*node[i][n2].v-(1.0/2.0)*node[i][n2].u;
		lattice[i][n2][8].f_tilde = lattice[i][n2][6].f_tilde-0.5*(lattice[i][n2][1].f_tilde-lattice[i][n2][3].f_tilde)-(1.0/6.0)*node[i][n2].v+(1.0/2.0)*node[i][n2].u;
	
	//	lattice[i][n2][4].f_tilde = lattice[i][n1+1][4].f;
	//	lattice[i][n2][7].f_tilde = lattice[i+1][n1+1][7].f;
	//	lattice[i][n2][8].f_tilde = lattice[i-1][n1+1][8].f;
	}
	}
	
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
	int j;
	#pragma omp for
	for(j=n1+1;j<=n2-1;++j)
	{
		//inlet (1,5,8) (velocity inlet)
	//	node[m1][j].u = node[m1+1][j].u;
	//	node[m1][j].v = node[m1+1][j].v;
		node[m1][j].u = 0.0;
		node[m1][j].v = 0.0;
		node[m1][j].rho = node[m1][j].u+(lattice[m1][j][0].f_tilde+lattice[m1][j][2].f_tilde+lattice[m1][j][4].f_tilde+2.0*(lattice[m1][j][3].f_tilde+lattice[m1][j][6].f_tilde+lattice[m1][j][7].f_tilde));
		lattice[m1][j][1].f_tilde = lattice[m1][j][3].f_tilde+(2.0/3.0)*node[m1][j].u;
		lattice[m1][j][5].f_tilde = lattice[m1][j][7].f_tilde-0.5*(lattice[m1][j][2].f_tilde-lattice[m1][j][4].f_tilde)+(1.0/6.0)*node[m1][j].u+(1.0/2.0)*node[m1][j].v;
		lattice[m1][j][8].f_tilde = lattice[m1][j][6].f_tilde+0.5*(lattice[m1][j][2].f_tilde-lattice[m1][j][4].f_tilde)+(1.0/6.0)*node[m1][j].u-(1.0/2.0)*node[m1][j].v;
		
		//outlet (3,6,7) (pressure outlet)
	//	node[m2][j].u = node[m2-1][j].u;
	//	node[m2][j].v = node[m2-1][j].v;
		node[m2][j].u = 0.0;
		node[m2][j].v = 0.0;
		node[m2][j].rho = -node[m2][j].u+(lattice[m2][j][0].f_tilde+lattice[m2][j][2].f_tilde+lattice[m2][j][4].f_tilde+2.0*(lattice[m2][j][1].f_tilde+lattice[m2][j][5].f_tilde+lattice[m2][j][8].f_tilde));
		lattice[m2][j][3].f_tilde = lattice[m2][j][1].f_tilde-(2.0/3.0)*node[m2][j].u;
		lattice[m2][j][7].f_tilde = lattice[m2][j][5].f_tilde+0.5*(lattice[m2][j][2].f_tilde-lattice[m2][j][4].f_tilde)-(1.0/6.0)*node[m2][j].u-(1.0/2.0)*node[m2][j].v;
		lattice[m2][j][6].f_tilde = lattice[m2][j][8].f_tilde-0.5*(lattice[m2][j][2].f_tilde-lattice[m2][j][4].f_tilde)-(1.0/6.0)*node[m2][j].u+(1.0/2.0)*node[m2][j].v;
	}
	}
	
	//corners
		//bottom left
			i=m1; j=n1;
			node[i][j].rho = node[i][j+1].rho;
		//	node[i][j].u = node[i][j-1].u;
		//	node[i][j].v = node[i][j-1].v;
			node[i][j].u = 0.0;
			node[i][j].v = 0.0;
			lattice[i][j][1].f_tilde = lattice[i][j][3].f_tilde+(2.0/3.0)*vector[1].e_x*node[i][j].u;
			lattice[i][j][2].f_tilde = lattice[i][j][4].f_tilde+(2.0/3.0)*vector[2].e_y*node[i][j].v;
			lattice[i][j][5].f_tilde = lattice[i][j][7].f_tilde+(1.0/6.0)*(vector[5].e_x*node[i][j].u+vector[5].e_y*node[i][j].v);
			lattice[i][j][6].f_tilde = (1.0/12.0)*(vector[6].e_x*node[i][j].u+vector[6].e_y*node[i][j].v);
			lattice[i][j][8].f_tilde = (1.0/12.0)*(vector[8].e_x*node[i][j].u+vector[8].e_y*node[i][j].v);
			lattice[i][j][0].f_tilde = node[i][j].rho - (lattice[i][j][1].f_tilde+lattice[i][j][2].f_tilde+lattice[i][j][3].f_tilde+lattice[i][j][4].f_tilde+lattice[i][j][5].f_tilde+lattice[i][j][6].f_tilde+lattice[i][j][7].f_tilde+lattice[i][j][8].f_tilde);
		
		//top left
			i=m1; j=n2;
			node[i][j].rho = node[i][j-1].rho;
		//	node[i][j].u = node[i][j-1].u;
		//	node[i][j].v = node[i][j-1].v;
			node[i][j].u = 0.0;
			node[i][j].v = 0.0;
			lattice[i][j][1].f_tilde = lattice[i][j][3].f_tilde+(2.0/3.0)*vector[1].e_x*node[i][j].u;
			lattice[i][j][4].f_tilde = lattice[i][j][2].f_tilde+(2.0/3.0)*vector[4].e_y*node[i][j].v;
			lattice[i][j][8].f_tilde = lattice[i][j][6].f_tilde+(1.0/6.0)*(vector[8].e_x*node[i][j].u+vector[8].e_y*node[i][j].v);
			lattice[i][j][5].f_tilde = (1.0/12.0)*(vector[5].e_x*node[i][j].u+vector[5].e_y*node[i][j].v);
			lattice[i][j][7].f_tilde = (1.0/12.0)*(vector[7].e_x*node[i][j].u+vector[7].e_y*node[i][j].v);
			lattice[i][j][0].f_tilde = node[i][j].rho - (lattice[i][j][1].f_tilde+lattice[i][j][2].f_tilde+lattice[i][j][3].f_tilde+lattice[i][j][4].f_tilde+lattice[i][j][5].f_tilde+lattice[i][j][6].f_tilde+lattice[i][j][7].f_tilde+lattice[i][j][8].f_tilde);
			
		//bottom right
			i=m2; j=n1;
			node[i][j].rho = node[i][j+1].rho;
		//	node[i][j].u = node[i][j-1].u;
		//	node[i][j].v = node[i][j-1].v;
			node[i][j].u = 0.0;
			node[i][j].v = 0.0;	
			lattice[i][j][3].f_tilde = lattice[i][j][1].f_tilde+(2.0/3.0)*vector[3].e_x*node[i][j].u;
			lattice[i][j][2].f_tilde = lattice[i][j][4].f_tilde+(2.0/3.0)*vector[2].e_y*node[i][j].v;
			lattice[i][j][6].f_tilde = lattice[i][j][8].f_tilde+(1.0/6.0)*(vector[6].e_x*node[i][j].u+vector[6].e_y*node[i][j].v);
			lattice[i][j][5].f_tilde = (1.0/12.0)*(vector[5].e_x*node[i][j].u+vector[5].e_y*node[i][j].v);
			lattice[i][j][7].f_tilde = (1.0/12.0)*(vector[7].e_x*node[i][j].u+vector[7].e_y*node[i][j].v);
			lattice[i][j][0].f_tilde = node[i][j].rho - (lattice[i][j][1].f_tilde+lattice[i][j][2].f_tilde+lattice[i][j][3].f_tilde+lattice[i][j][4].f_tilde+lattice[i][j][5].f_tilde+lattice[i][j][6].f_tilde+lattice[i][j][7].f_tilde+lattice[i][j][8].f_tilde);
		
		//top right
			i=m2; j=n2;
			node[i][j].rho = node[i][j-1].rho;
		//	node[i][j].u = node[i][j-1].u;
		//	node[i][j].v = node[i][j-1].v;
			node[i][j].u = 0.0;
			node[i][j].v = 0.0;
			lattice[i][j][3].f_tilde = lattice[i][j][1].f_tilde+(2.0/3.0)*vector[3].e_x*node[i][j].u;
			lattice[i][j][4].f_tilde = lattice[i][j][2].f_tilde+(2.0/3.0)*vector[4].e_y*node[i][j].v;
			lattice[i][j][7].f_tilde = lattice[i][j][5].f_tilde+(1.0/6.0)*(vector[7].e_x*node[i][j].u+vector[7].e_y*node[i][j].v);
			lattice[i][j][6].f_tilde = (1.0/12.0)*(vector[6].e_x*node[i][j].u+vector[6].e_y*node[i][j].v);
			lattice[i][j][8].f_tilde = (1.0/12.0)*(vector[8].e_x*node[i][j].u+vector[8].e_y*node[i][j].v);
			lattice[i][j][0].f_tilde = node[i][j].rho - (lattice[i][j][1].f_tilde+lattice[i][j][2].f_tilde+lattice[i][j][3].f_tilde+lattice[i][j][4].f_tilde+lattice[i][j][5].f_tilde+lattice[i][j][6].f_tilde+lattice[i][j][7].f_tilde+lattice[i][j][8].f_tilde);
	
}

/////////////////////////// calculating macro variables /////////////////////////////////
void rho_momentum()
{
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
		int i;
		double sum=0.0, sumu = 0.0, sumv = 0.0;
		#pragma omp for
		for(i=m1;i<=m2;++i)
		{int j, k;
		for(j=n1;j<=n2;++j)
		{
			sum=0.0; sumu=0.0; sumv=0.0;
			for(k=0;k<=alpha;++k)
			{
				sum = sum + lattice[i][j][k].f_tilde;
				sumu = sumu + vector[k].e_x*lattice[i][j][k].f_tilde;
				sumv = sumv + vector[k].e_y*lattice[i][j][k].f_tilde;
			}
			node[i][j].rho = sum;
			node[i][j].mom_x = sumu;
			node[i][j].mom_y = sumv;
		}}
	}
}

/////////////////// interpolation velocity to lagrangian nodes //////////////////////////
void interpolate_to_lag_nodes()
{
	int x1,y1;
	double dist_x=0.0, dist_y=0.0, wt_x=0.0, wt_y=0.0;
	
//	FILE *bdy3;
//	bdy3 = fopen("particle_vel.dat", "w");
	
	for(p=0;p<P;++p)
	{
		for(n=1;n<=N[p];++n)
		{
		//	particle[n].u = 0.0;
		//	particle[n].v = 0.0;
			particle[p][n].rho = 0.0;
			particle[p][n].mom_x = 0.0;
			particle[p][n].mom_y = 0.0;
		
			x1 = (int) (particle[p][n].x);
			y1 = (int) (particle[p][n].y);
		
		//	printf("%d\t%d\n", x1, y1);
	
			for(i=x1-3;i<=x1+3;++i)
			for(j=y1-3;j<=y1+3;++j)
			{
				dist_x = particle[p][n].x-i;
				dist_y = particle[p][n].y-j;
			
				if(fabs(dist_x)>=0.0 && fabs(dist_x)<=1.0)
				{	wt_x = (1.0-fabs(dist_x));	}
				else if(fabs(dist_x)>=1.0)
				{	wt_x = 0.0;	}
			
				if(fabs(dist_y)>=0.0 && fabs(dist_y)<=1.0)
				{	wt_y = (1.0-fabs(dist_y));	}
				else if(fabs(dist_y)>=1.0)
				{	wt_y = 0.0;	}
			
			//	particle[n].u = particle[n].u + node[i][j].u * wt_x * wt_y;
			//	particle[n].v = particle[n].v + node[i][j].v * wt_x * wt_y;
				particle[p][n].mom_x = particle[p][n].mom_x + node[i][j].mom_x * wt_x * wt_y;
				particle[p][n].mom_y = particle[p][n].mom_y + node[i][j].mom_y * wt_x * wt_y;
				particle[p][n].rho = particle[p][n].rho + node[i][j].rho * wt_x * wt_y;
			}
		
	//		fprintf(bdy3, "%lf\t%lf\n", particle[n].u, particle[n].v);
		}
	}
	
//	fclose(bdy3);
}

///////////////// calculate force on lagrangian particles ///////////////////////////////
void compute_lag_node_force()
{
//	area = 2.0*M_PI*r/N;
	
	for(p=0;p<P;++p)
	{
		omp_set_num_threads(num_threads);
		#pragma omp parallel
		{
			int n;
			#pragma omp for
			for(n=1;n<=N[p];++n)
			{
				particle[p][n].fp_x = 2.0*(particle[p][n].u-particle[p][n].mom_x);//*area;
				particle[p][n].fp_y = 2.0*(particle[p][n].v-particle[p][n].mom_y);//*area;
	
			//	particle[n].fp_x = -2.0*particle[n].u*particle[n].rho*area;
			//	particle[n].fp_y = -2.0*particle[n].v*particle[n].rho*area;
			}
		}
	}
	
//	FILE *force;
//	force = fopen("particle_force.dat", "w");
//	for(n=1;n<=N;++n)
//		fprintf(force, "%d\t%lf\t%lf\n", n, particle[n].fp_x, particle[n].fp_y);
//	fclose(force);
	
}

///////////// spread the computed force to Eulerian nodes ///////////////////////////////
void spread_force_to_eul_nodes()
{
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{int i;
	#pragma omp for
	for(i=m1;i<=m2;++i)
	{int j;
	for(j=n1;j<=n2;++j)
	{
		node[i][j].fn_x = 0.0;
		node[i][j].fn_y = 0.0;
	}}}
	
	int x1,y1;
	double dist_x=0.0, dist_y=0.0, wt_x=0.0, wt_y=0.0;
	
//	FILE *weight;
//	weight = fopen("interpolation_weight.dat", "w");
	
	for(p=0;p<P;++p)
	{
		del_s[p] = surface_area[p]/N[p];	//area of each arc
		
		for(n=1;n<=N[p];++n)
		{
			x1 = (int) (particle[p][n].x);
			y1 = (int) (particle[p][n].y);
	
			for(i=x1-3;i<=x1+3;++i)
			for(j=y1-3;j<=y1+3;++j)
			{
				dist_x = particle[p][n].x-i;
				dist_y = particle[p][n].y-j;
			
				if(fabs(dist_x)>=0.0 && fabs(dist_x)<=1.0)
				{	wt_x = (1.0-fabs(dist_x));	}
				else if(fabs(dist_x)>=1.0)
				{	wt_x = 0.0;	}
			
				if(fabs(dist_y)>=0.0 && fabs(dist_y)<=1.0)
				{	wt_y = (1.0-fabs(dist_y));	}
				else if(fabs(dist_y)>=1.0)
				{	wt_y = 0.0;	}
			
	//			fprintf(weight, "%d\t%d\t%d\t%lf\t%lf\n", n , x1, y1, wt_x, wt_y);
			
				node[i][j].fn_x = node[i][j].fn_x + particle[p][n].fp_x * wt_x * wt_y * del_s[p];
				node[i][j].fn_y = node[i][j].fn_y + particle[p][n].fp_y * wt_x * wt_y * del_s[p];
			}
		}
	}
	
//	fclose(weight);
}

//////////////////////// updating velocity //////////////////////////////////////////////
void u_v_calculator()
{
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
		int i;
		#pragma omp for
		for(i=m1;i<=m2;++i)
		{int j;
		for(j=n1;j<=n2;++j)
		{
				node[i][j].u = (node[i][j].mom_x + 0.5*node[i][j].fn_x);
				node[i][j].v = (node[i][j].mom_y + 0.5*node[i][j].fn_y);
		}}
	}
}

///////////////////////////calculating coefficient of drag //////////////////////////////
void force_coefficients()
{
/*
	for(p=0;p<P;++p)
	{
		double a=0.0, b=0.0;
		sum=0.0; sumu=0.0; sumv=0.0;
		temp1=0.0; temp2=0.0; temp3=0.0;
		phi=0.0; theta=0.0;
		for(i=((int) (Cx[p]-axis[p]));i<=((int) (Cx[p]+axis[p]));++i)
		for(j=((int) (Cy[p]-axis[p]));j<=((int) (Cy[p]+axis[p]));++j)
		if(domain[i][j][p].flag==1.0)
		{
			sumu = sumu + node[i][j].rho*(node[i][j].u-node[i][j].u_old);
			sumv = sumv + node[i][j].rho*(node[i][j].v-node[i][j].v_old);
			
			a=i-Cx[p];
			b=j-Cy[p];
			theta = fabs(atan(fabs((b)/(a))));
				if(a>0.0)
				{
					if(b<0.0)
					{	theta = 2.0*M_PI - theta;	}
					else if(b==0.0)
					{	theta = 0.0;	}
					else if(b>0.0)
					{	theta = theta;	}
				}
				else if(a<0.0)
				{
					if(b>0.0)
					{	theta = M_PI - theta;	}
					else if(b<0.0)
					{	theta = M_PI + theta;	}
					else if(b==0.0)
					{	theta = M_PI;	}
				}
				else if(a==0.0)
				{
					if(b>0.0)
					{	theta = M_PI/2.0;	}
					else if(b<0.0)
					{	theta = 3.0*M_PI/2.0;	}
					else if(b==0.0)
					{	theta = 0.0000000000000001;	}
				}
				
			a=node[i][j].u;
			b=node[i][j].v;
			phi = fabs(atan(fabs((b)/(a))));
				if(a>0.0)
				{
					if(b<0.0)
					{	phi = 2.0*M_PI - phi;	}
					else if(b==0.0)
					{	phi = 0.0;	}
					else if(b>0.0)
					{	phi = phi;	}
				}
				else if(a<0.0)
				{
					if(b>0.0)
					{	phi = M_PI - phi;	}
					else if(b<0.0)
					{	phi = M_PI + phi;	}
					else if(b==0.0)
					{	phi = M_PI;	}
				}
				else if(a==0.0)
				{
					if(b>0.0)
					{	phi = M_PI/2.0;	}
					else if(b<0.0)
					{	phi = 3.0*M_PI/2.0;	}
					else if(b==0.0)
					{	phi = 0.0000000000000001;	}
				}
			
			a=i-Cx[p];
			b=j-Cy[p];
			temp1 = sqrt( a*a+b*b );
			temp2 = M_PI/2.0 - phi + theta;
			temp3 = sqrt(node[i][j].u*node[i][j].u+node[i][j].v*node[i][j].v) - sqrt(node[i][j].u_old*node[i][j].u_old+node[i][j].v_old*node[i][j].v_old);
			sum = sum + node[i][j].rho*temp3*cos(temp2)*temp1;
		}
		int_fx[p] = sumu;
		int_fy[p] = sumv;
		int_tor[p] = sum;
	}
*/
	//added mass effect
	for(p=0;p<P;++p)
	{
		int_fx[p] = vol[p]*(Uc[p]-Uc_old[p]);
		int_fy[p] = vol[p]*(Vc[p]-Vc_old[p]);
		int_tor[p] = vol[p]*mom_inertia[p]*(V_ang[p]-V_ang_old[p]);
	}
	
	
	sum=0.0; sumu=0.0; sumv=0.0;
	for(p=0;p<P;++p)
	{	
		del_s[p] = surface_area[p]/N[p];
		
//		omp_set_num_threads(num_threads);
//		#pragma omp parallel
		{
//			int n;
			sum=0.0; sumu=0.0; sumv=0.0;
//			#pragma omp for reduction (+:sum,sumu,sumv)
			for(n=1;n<=N[p];++n)
			{
				sum = sum + particle[p][n].rho;
				sumu = sumu + particle[p][n].fp_x * del_s[p];
				sumv = sumv + particle[p][n].fp_y * del_s[p];
			}
			rho_avg[p] = sum/N[p];
			C_d[p] = -sumu;//(rho_avg*U*U*r);
			C_l[p] = -sumv;//(rho_avg*U*U*r);
		}
	}
}

///////////////////// update the lagrangian particle position ///////////////////////////
void update_lag_node_position()
{
	sum=0.0; sumu=0.0; sumv=0.0;
	temp1=0.0; temp2=0.0; temp3=0.0; temp4=0.0; temp5=0.0;
	a=0.0, b=0.0;
	
//	FILE *phi1;
//	phi1 = fopen("phi.dat", "w");
	
	//calculating the forces and torque
	for(p=0;p<P;++p)
	{
		del_s[p] = surface_area[p]/N[p];
		sum=0.0; sumu=0.0; sumv=0.0;
		temp1=0.0; temp2=0.0; temp3=0.0; temp4=0.0; temp5=0.0; a=0.0; b=0.0;
	//	omp_set_num_threads(num_threads);
	//	#pragma omp parallel
		{
	//		int n;
	//		#pragma omp for reduction (+:sum,sumu,sumv)
			for(n=1;n<=N[p];++n)
			{
				sumu = sumu + particle[p][n].fp_x*del_s[p];
				sumv = sumv + particle[p][n].fp_y*del_s[p];
				sum = sum + ( particle[p][n].fp_y*( (particle[p][n].x-Cx[p]) ) - particle[p][n].fp_x*( (particle[p][n].y-Cy[p]) ) )*del_s[p];
			}
			F_x[p] = -sumu;	// viscous focrce:x
			F_y[p] = -sumv;	// viscous focrce:y
			Tor[p] = -sum;
		}
	//	fclose(phi1);
	}
	

/*	//....particle wall collision....
	temp1=0.0; temp2=0.0; temp3=0.0;
	for(p=0;p<P;++p)
	{
		temp1 = 1.0*r[p]+1.0;
		temp2 = fabs(Cx[p]-((double) (m1)));
		temp3 = fabs(Cx[p]-((double) (m2)));
		if(temp2>temp1 || temp2>temp3)
		{	F_pw[p] = 0.0;	}
		else if(temp2<=temp1)
		{	F_pw[p] = 2.4*r[p]*r[p]*( 2.0*(pow((r[p]/temp2),14))-(pow((r[p]/temp2),8)) ) * (Cx[p]-((double) (m1)))/(r[p]*r[p]);	}
		else if(temp2<=temp1)
		{	F_pw[p] = 2.4*r[p]*r[p]*( 2.0*(pow((r[p]/temp2),14))-(pow((r[p]/temp2),8)) ) * (Cx[p]-((double) (m2)))/(r[p]*r[p]);	}
	}

	//.....particle particle collison.....
	temp1=0.0; temp2=0.0; temp3=0.0;
	sumu=0.0; sumv=0.0;
	temp1 = r[0]+r[1]+1.0;
	temp2 = sqrt((Cx[0]-Cx[1])*(Cx[0]-Cx[1])+(Cy[0]-Cy[1])*(Cy[0]-Cy[1]));
	if(temp2>temp1)
	{	for(p=0;p<P;++p)
		{	F_pp_x[p] = 0.0;
			F_pp_y[p] = 0.0;	}
	}
	else if(temp2<=temp1)
	{
		F_pp_x[0] = 2.4*r[0]*r[0]*( 2.0*(pow((2.0*r[0]/temp2),14))-(pow((2.0*r[0]/temp2),8)) ) * (Cx[0]-Cx[1])/(4.0*r[0]*r[0]);
		F_pp_y[0] = 2.4*r[0]*r[0]*( 2.0*(pow((2.0*r[0]/temp2),14))-(pow((2.0*r[0]/temp2),8)) ) * (Cy[0]-Cy[1])/(4.0*r[0]*r[0]);
		
		F_pp_x[1] = 2.4*r[1]*r[1]*( 2.0*(pow((2.0*r[1]/temp2),14))-(pow((2.0*r[1]/temp2),8)) ) * (Cx[1]-Cx[0])/(4.0*r[1]*r[1]);
		F_pp_y[1] = 2.4*r[1]*r[1]*( 2.0*(pow((2.0*r[1]/temp2),14))-(pow((2.0*r[1]/temp2),8)) ) * (Cy[1]-Cy[0])/(4.0*r[1]*r[1]);
	}
*/	
	
	for(p=0;p<P;++p)
	{	
		Uc_old[p] = Uc[p];
		Vc_old[p] = Vc[p];
		V_ang_old[p] = V_ang[p];
		
		grav[p] = (1.0/rho_p_by_f[p]-1.0)*(rho_p_by_f[p]*vol[p])*gravity;
		F_tot_x[p] = F_x[p] + int_fx[p];// + F_pw[p] + F_pp_x[p];
		F_tot_y[p] = F_y[p] + int_fy[p] + grav[p];// + F_pp_y[p];
		Tor_tot[p] = Tor[p] + int_tor[p];
		
		Uc[p] = Uc[p] + F_tot_x[p]/(rho_p_by_f[p]*vol[p]);
		Vc[p] = Vc[p] + F_tot_y[p]/(rho_p_by_f[p]*vol[p]);
		V_ang[p] = V_ang[p] + Tor_tot[p]/(rho_p_by_f[p]*vol[p]*mom_inertia[p]);
	
		Cx[p] = Cx[p] + Uc[p];
		Cy[p] = Cy[p] + Vc[p];
		
		psi[p] = psi[p] + V_ang[p];	//psi = orientation of particle

	//	omp_set_num_threads(num_threads);
	//	#pragma omp parallel
		{
	//		int n;
	//		#pragma omp for
			for(n=1;n<=N[p];++n)
			{
				particle[p][n].angle = particle[p][n].angle + V_ang[p];
				
				if(particle[p][n].angle>=(2.0*M_PI))
				{	particle[p][n].angle = particle[p][n].angle - 2.0*M_PI;	}
				else if(particle[p][n].angle<0.0)
				{	particle[p][n].angle = particle[p][n].angle + 2.0*M_PI;	}
				
				particle[p][n].x = Cx[p] + particle[p][n].r*cos(particle[p][n].angle);
				particle[p][n].y = Cy[p] + particle[p][n].r*sin(particle[p][n].angle);
				
				particle[p][n].V_tang = V_ang[p]*particle[p][n].r;
				particle[p][n].u = Uc[p] + particle[p][n].V_tang*cos(M_PI/2.0 + particle[p][n].angle);
				particle[p][n].v = Vc[p] + particle[p][n].V_tang*sin(M_PI/2.0 + particle[p][n].angle);
			}
		}
	}
	
	/*
	// flag denotes solid(=1.0) or fluid(=0.0) Eulerian nodes
	for(p=0;p<P;++p)
	{
		temp1=0.0;
		for(i=((int) (Cx[p]-1.1*axis[p]));i<=((int) (Cx[p]+1.1*axis[p]));++i)
		for(j=((int) (Cy[p]-1.1*axis[p]));j<=((int) (Cy[p]+1.1*axis[p]));++j)
		{
			temp1 = (pow(((i-Cx[p])*cos(psi[p])+(j-Cy[p])*sin(psi[p])),2.0))/(pow(axis[p],2.0)) + (pow(((i-Cx[p])*sin(psi[p])-(j-Cy[p])*cos(psi[p])),2.0))/(pow((axis[p]/axis_ratio),2.0));
			if(temp1<=1.0)
				domain[i][j][p].flag = 1.0;
			else
				domain[i][j][p].flag = 0.0;
		}
	}
	*/
	
//	FILE *bdy2;
//	bdy2 = fopen("particle.dat", "w");
//	for(n=1;n<=N;++n)
//		fprintf(bdy2, "%lf\t%lf\n", particle[n].x, particle[n].y);
//	fclose(bdy2);
}

///////////////////////////////// printing everything ///////////////////////////////////
void printing_everything()
{
//	dx = ((double) (m2/(n2-n1))) /(m2-m1);
//	dy = ((double) (n2/(n2-n1))) /(n2-n1);
	dx = 1.0/(2.0*axis[0]);
	dy = 1.0/(2.0*axis[0]);
	FILE *all, *pdf;
	all = fopen("all.dat","w");
	pdf = fopen("pdf.dat","w");
	fprintf(all, "ZONE\tI=%d\tJ=%d\n", (m2-m1+1), (n2-n1+1));
	for(j=n1;j<=n2;++j)
	for(i=m1;i<=m2;++i)
	{
			node[i][j].vor = (node[i+1][j].v-node[i-1][j].v)/dx - (node[i][j+1].u-node[i][j-1].u)/dy;
		fprintf(all, "%d\t%d\t%lf\t%lf\t%lf\t%lf\n", i, j, node[i][j].u/U, node[i][j].v/U, node[i][j].rho, node[i][j].vor/U);
		for(k=0;k<=alpha;++k)
			fprintf(pdf, "%d\t%d\t%d\t%lf\t%lf\n", i, j, k, lattice[i][j][k].f_tilde, lattice[i][j][k].f_0);
	}
	fclose(all);
	fclose(pdf);
}

void time_files()
{
/*	count=((int) t);
	if(count<=500001)
	{
		FILE *files[2001];
		FILE *files1[2001];
		FILE *files2[2001];
		if((count%20)==0)
		{
			char filename[100];
			char filename1[100];
			char filename2[100];
			sprintf(filename, "all_%d.dat", count);
			sprintf(filename1, "particle1_%d.dat", count);
			sprintf(filename2, "particle2_%d.dat", count);
			files[nihal] = fopen(filename, "w");
			files1[nihal] = fopen(filename1, "w");
			files2[nihal] = fopen(filename2, "w");
			fprintf(files[nihal], "ZONE T=\"%lf\"\tSTRANDID=1\tSOLUTIONTIME=%lf\tI=%d\tJ=%d\n", t, t, (m2-m1+1), (n2-n1+1));
			fprintf(files1[nihal], "ZONE T=\"%lf\"\tSTRANDID=1\tSOLUTIONTIME=%lf\n", t, t, (m2-m1+1), (n2-n1+1));
			fprintf(files2[nihal], "ZONE T=\"%lf\"\tSTRANDID=1\tSOLUTIONTIME=%lf\n", t, t, (m2-m1+1), (n2-n1+1));
			
			for(j=n1;j<=n2;++j)
			for(i=m1;i<=m2;++i)
			//	fprintf(files[nihal], "%d\t%d\t%1.9lf\t%1.9lf\t%1.9lf\n", i, j, node[i][j].u/U, node[i][j].v/U, node[i][j].vor/U);
				fprintf(files[nihal], "%d\t%d\t%1.9lf\n", i, j, node[i][j].vor/U);
			
			for(n=1;n<=63;++n)
			{	fprintf(files1[nihal], "%1.9lf\t%1.9lf\n", particle[0][n].x, particle[0][n].y);
				fprintf(files2[nihal], "%1.9lf\t%1.9lf\n", particle[1][n].x, particle[1][n].y);			}
		
			fclose(files[nihal]);
			fclose(files1[nihal]);
			fclose(files2[nihal]);
			nihal++;
		}
	} */
}
