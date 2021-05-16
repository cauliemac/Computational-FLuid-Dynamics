/*! \file exact_adiabatic.c
 * 
 * \brief Contains an exact Riemann solver for adiabatic hydro.
 *
 * This solver is borrowed from the 3D code of T.P. Downes and is
 * based on Sam Falle's notes for hydro Riemann solvers.  Note that we
 * have changed gamma_calc() to simplify its interaction with the rest
 * of the code a little.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "defs.h"
//#include "typedefs.h"
#include "functions.h"	


#define JMAX 20 //max amount of iterations for Newton-Ralphson iterator
#define gamma 5/3 //value of gamma for monatomic gas

/*! Newton-Raphson support function for calculating f and df
 *
 * @param[in] x Argument at which f and df are to be evaluated.
 * @param[out] f pointer to value of f
 * @param[out] df pointer to value of the derivative of f
 * @param[in] t1 constant for iteration over rarefactions
 * @param[in] t2 constant for iteration over rarefactions
 * @param[in] t3 constant for iteration over rarefactions
 * @param[in] gl Pointer to array of functions of \gamma in left state
 * @param[in] gr Pointer to array of functions of \gamma in right state
 */
void func_newt(double x,double *f, double *df, double t1, double t2, double t3, double gl, double gr)
{

	*f=t1*pow(x,gl)+t2*pow(x,gr)-t3;

	*df=gl*t1*pow(x,-gl)+gr*t2*pow(x,-gr);

}

/*! Newton-Raphson iterator
 *
 * @param[in] x Initial guess for solution
 * @param[in] xmin Lower bound for solution
 * @param[in] xmax Upper bound for solution
 * @param[in] t1 Constant for iteration
 * @param[in] t2 Constant for iteration
 * @param[in] t3 Constant for iteration
 * @param[in] gl Pointer to array of functions of \gamma in left state
 * @param[in] gr Pointer to array of functions of \gamma in right state
 * @param[in] newt_check A flag to check whether we converged.
 */

double newt(double x, double xmin, double xmax, double t1, double t2, double t3, double gl, double gr, int *newt_check)
{
	int j;
	double df,dx,f,rtn,greater;

	rtn=x;
	*newt_check=0;
	for (j=1;j<=JMAX;j++)
	{
		func_newt(rtn,&f,&df, t1, t2, t3, gl, gr);
		dx=f/df;
		if(dx>0.0) greater=rtn;
		rtn -= dx;

/* Correction for leaping to negative pressures */

		if(rtn<=0.0)  /* If the guess is wrong do a bisection */ 
			rtn=(rtn+dx)/2.0;

		if ((xmin-rtn)*(rtn-xmax) < 0.0)
			printf("Jumped out of brackets in RTNEWT");
		if (fabs(dx/rtn) < 1.0e-06) return rtn;
		if (fabs(rtn/t3) < 1.0e-10)// If the function is crap
		{
			*newt_check=1;
			rtn=greater;
			printf("Switching to bisection\n");
			return rtn;
		}
	}
	*newt_check=1;
	rtn=greater;
	return rtn;
}

/*! Function evaluation for bisection method iterator
 *
 * @param[in] x Initial guess for solution
 * @param[in] t1 Constant for iteration
 * @param[in] t2 Constant for iteration
 * @param[in] t3 Constant for iteration
 * @param[in] gl Pointer to array of functions of \gamma in left state
 * @param[in] gr Pointer to array of functions of \gamma in right state
 */

double func_bis(double x, double t1, double t2, double t3, double gl, double gr)
{

	double f;

	f=t1*pow(x,gl)+t2*pow(x,gr)-t3;

	return(f);

}

/*! Bisection method iterator
 *
 * @param[in] x1 Initial guess for solution
 * @param[in] x2 Initial guess for solution
 * @param[in] xacc Required accuracy in x
 * @param[in] t1 Constant for iteration
 * @param[in] t2 Constant for iteration
 * @param[in] t3 Constant for iteration
 * @param[in] gl Pointer to array of functions of \gamma in left state
 * @param[in] gr Pointer to array of functions of \gamma in right state
 */


double rtbis(double x1,double x2,double xacc, double t1, double t2, double t3, double gl, double gr)
{

	int j;
	double dx,f,fmid,xmid,rtb;

	f=func_bis(x1, t1, t2, t3, gl, gr);
	fmid=func_bis(x2, t1, t2, t3, gl, gr);
	if (f*fmid >= 0.0)
	{
		rtb=x2*1.0e+10;
		return (rtb);
	}
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++)
	{
		fmid=func_bis(xmid=rtb+(dx *= 0.5), t1, t2, t3, gl, gr);
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx/rtb) < xacc || fmid == 0.0) return rtb;
	}
	printf("Too many bisections in rtbis");
	return 0.0;
}


/*! Calculates the resolved state given the left and right states.
 *
 *  This is based on the notes from Sam Falle, used in the 3D hydro code of T Downes.  It
 *  is a non-linear adiabatic Riemann solver extended to allow for differences in the
 *  ratio of specific heats on either side of the interface.
 *
 * @param[in] left_state State to the left of the interface
 * @param[in] right_state State to the right of the interface
 * @param[in] perp Direction normal to interface
 * @param[out] resolved_state Pointer to the resolved state.
 */

//void adiflux(double left_state, double right_state, double dx, double dt, int perp, double *resolved_state)
void adiflux(cell_state temp_cell_state, int Left, int Right, double dx, double dt, interface_cell_state* riemann_cell_state)
{
	extern double CFL; //GAMMA[N_CHARGED_FLUIDS];
  /*
   * Local variables:-
   */
	int iter, newt_check, equal_gamma=1;
	double cl,cr,wleft,wright,u,v,w,c,ci,p,pi,gl,gr, t1, t2, t3;
	double rhol, ul, pl, rhor, ur, pr;
	double density, c_v;

/* First get our primitive variables */

	//rhol= (*temp_cell_state).Density[Left];//rho left (density ρ)
	rhol= temp_cell_state.Density[Left];//rho left (density ρ)
	ul  = temp_cell_state.Velocity[Left]/rhol;//speed left
	pl =  temp_cell_state.Pressure[Left];// pressure left

	rhor= temp_cell_state.Density[Right];//rho right (density ρ)
	ur  = temp_cell_state.Velocity[Right]/rhor;//speed right
	pr =  temp_cell_state.Pressure[Right];// pressure right

	//printf("\nadiflux function started");
	//TODO

  /***************************************************************************
   *                                                                         *
   *                         Linear Riemann solver                           *
   *                                                                         *
   ***************************************************************************/

/* Calculate various algebraic expressions involving the specific heat
 * capacity (for the neutrals, 'cos that's what this Riemann solver is
 * for). For the moment we're assuming this is uniform throughout
 * the domain, though most coding is in place to avoid this
 * assumption. */

  //no need for gamma_calc cos monatomic
  gl = gamma; //gamma left
  gr = gamma; //gamma right

  cl = sqrt(gl*pl*rhol);
  cr = sqrt(gr*pr*rhor);
  wleft = -cl;
  wright = cr;

  //find value of the pressure
  p = (wright*pl - wleft*pr + wright*wleft*(ur - ul))/(wright - wleft);
  //printf("\n\npressure is equal to %f ",p);
  //TODO remove this

  if(p<= 0.0)
  {
    //negative pressure cannot exist, so set it to a small number instead
	printf("Negative pressure in linear Riemann solver\n");
    printf("p = %f, pl = %f, pr = %f\n",p,pl,pr); 
    p = 1.0e-06;
  }

  /***************************************************************************
   *                                                                         *
   *                     End of Linear Riemann solver                        *
   *                                                                         *
   ***************************************************************************/

  /*
   *  Check if we need non-linear Riemann solver
   *  Will use the non_Linear Solver if the discrepancy between left and right
   *  states is greater than 10%
   */

  if((fabs(p-pl)>0.1*pl || fabs(p-pr)>0.1*pr))  
  {
     /***************************************************************************
      *                                                                         *
      *                       Non-linear Riemann solver                         *
      *                                                                         *
      ***************************************************************************/
	 	printf("\nnon-Linear Riemann function started\n");
		//TODO
		if((p < pl) && (p < pr))
		{

      /*
       * Non-iterative solution for two rarefactions
       */
			if(equal_gamma)// Note that gl and gr are the same by def
			{
				ci = cr/(rhor*pow(pr,gr));
				p = (0.5*gr*(ul - ur) + cr/rhor + cl/rhol)/(ci + cl/(rhol*pow(pl,gr)));

				if(p<0.0) // Gas is expanding to form a vacuum (bad news!)
				{
					printf("Negative pressures in rarefaction solver: %lf\n",pow(-p,gr));
					p=5.0e-11*(pr+pl); 
				}
				
				else 
					p = pow(p,gr);
			}
		}
    else
	{
      /* 
       * Iterative solution for shock/rarefaction
       *
       *        Initialise iteration
       */
      pi = p;
		wleft = -cl*wave(p,pl,gl);	//TODO FIX THE WAVE() FUNCTION
		wright = cr*wave(p,pr,gr);	//TODO FIX THE WAVE() FUNCTION
      iter = 0;
      p = (wright*pl - wleft*pr + wright*wleft*(ur - ul))/(wright - wleft);
		if(p<0.0) 
			p=1.0e-10; /* if <0 treat it as zero */
		pi = 0.5*p; /* Ensure we get into iteraction */

      /*
       *        End of initialisation
       */

      while(fabs((p-pi)/p) > 0.1 && iter < 20)
	  {
   		pi = p;
			wleft = -cl*wave(p,pl,gl);	//TODO FIX THE WAVE() FUNCTION
			wright = cr*wave(p,pr,gr);	//TODO FIX THE WAVE() FUNCTION
      	p = (wright*pl - wleft*pr + wright*wleft*(ur - ul))/(wright - wleft);
			if(p<0.0) 
				p=1.0e-10; /* if <0 treat it as zero */
			++iter;
      }
      if(iter >= 20)
	  {
		printf("Convergence failure in Riemann solver: pr=%e pl=%e p=%e\n", pr, pl, p);
      }

      /*
       * End of iterative solver
       */
    }
  }


  /*
   * Start calculation of resolved density and velocity
   */
  /* csos 
     if(p<1.0e-0){
     printf("resolved pressure < 0.0\n");
     }
  */

  wleft = -cl*wave(p,pl,gl);	//TODO FIX THE WAVE() FUNCTION
  wright = cr*wave(p,pr,gr);	//TODO FIX THE WAVE() FUNCTION
  u = (wright*ur - wleft*ul - pr + pl)/(wright -wleft);
  
  /*
   * Check for position of contact
   */
  if(u > 0.0)
  {
    /*
     * Contact is on right of interface
     */
    density = wleft*rhol/(wleft - (u -ul)*rhol);

	//g5 = gl;
	c_v = 1.0/(gl-1.0);
    
    /*
     * Check velocity of left waves
     */
    if(p < pl) 
    {
      /*
       * Left wave is rarefaction
       */
      c = sqrt(gl*p/(density));//speed of sound
      
      if((u-c) > 0.0) 
      {
	/*
	 * Tail of rarefaction is on right of interface
	 */
	cl = cl/rhol;
	if((ul-cl) < 0.0)
	{
	  /* 
	   * Head of rarefaction is on left of interface
	   * Left rarefaction spans interface
	   */
	  u = gl*(ul +gl*cl);
	  p = pl*pow((u/cl),gl);
	  density = gl*p/(u*u);
	}
	else
	{
	  /*
	   * Left rarefaction does not span interface
	   */
	  density = rhol;
	  u = ul;
	  p = pl;
	}
      }
    }
    else
    {
      /*
       * Left wave is shock
       */
      if((wleft/rhol+ul) > 0.0)
      {
	/*
	 * Shock is on right of interface
	 */
	density = rhol;
	u = ul;
	p = pl;
      }
    }
  }
  else
  {
    /* 
     * Contact is on left of interface
     */
    density = wright*rhor/(wright - (u - ur)*rhor);

	 //g5 = gr[5];
	 c_v = 1.0/(gr-1.0);
    
    /*
     * Check velocity of right waves
     */
    if(p < pr) 
    {
      /*
       * Right wave is rarefaction
       */
      c = sqrt(gr*p/(density));//speed of sound
      if((u+c) < 0.0)
      {
	/*
	 * Tail of rarefaction is on left of interface
	 */
	cr = cr/rhor;
	if((ur+cr) > 0.0)
	{
	  /*
	   * Head of rarefaction is on right of interface 
	   * Right rarefaction spans interface
	   */
	  u = gr*(ur - gr*cr);
	  p = pr*pow((-u/cr),gr);
	  density = gr*p/(u*u);
	}
	else
	{
	  /*
	   * Right rarefaction does not span interface
	   */
	  density = rhor;
	  u = ur;
	  p = pr;
	}
      }
    }
    else
    {
      /* 
       * Right wave is shock 
       */
      if((wright/rhor+ur) < 0.0)
      {
	/*
	 * Shock is on left of interface
	 */
	density = rhor;
	u = ur;
	p = pr;
      }
    }
  }

/* For very strong rarefactions we can end up with a negative density
 * (if we don't get convergence in our nonlinear solver).  If this
 * happens, just fix it to be the pressure (which is hopefully positive), but 
 * print an error */
  if(density < 0.0)
  {
	fprintf(stderr,"Had to fix density: %e -> %e\n", density, p);
	density = p;
  }
  if(isnan(p))// Checks if the pressure is Not A Number (NAN)
  {
	fprintf(stderr,"Got a nan for p\n");
	p = 0.5*(pl + pr);
	density = 0.5*(rhol + rhor);
	u = 0.5*(ul + ur);
  }

/* One thing for Sam's suggestion for artificial dissipation - the
 * resolved sound speed */

	if(u > 0.0)
		c=sqrt(gl*p/(density));
	else
		c=sqrt(gr*p/(density));

  /*
   *  End of resolved pressure, velocities and density calculation
   *
   * Calculate fluxes in resolved state
   */

	
	//resolved_state = [density,u,p,c_v];

	
	riemann_cell_state->Density = density;
	riemann_cell_state->Pressure = p;
	riemann_cell_state->Velocity = u;

	/*
	printf("inside adiabatic.c density = %f\n",density);
	printf("inside adiabatic.c pressure = %f\n",p);
	printf("inside adiabatic.c velocitty = %f\n",u);
	system("pause");
	*/
}