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
#include "defs.h"
#include "typedefs.h"
#include "functions.h"	

//#if ADIABATIC

#define JMAX 20
#define gamma = 5/3

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
void func_newt(REAL x,REAL *f, REAL *df, REAL t1, REAL t2, REAL t3, REAL *gl, REAL *gr)
{

	*f=t1*pow(x,gl[4])+t2*pow(x,gr[4])-t3;

	*df=gl[4]*t1*pow(x,-gl[3])+gr[4]*t2*pow(x,-gr[3]);	

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

REAL newt(REAL x, REAL xmin, REAL xmax, REAL t1, REAL t2, REAL t3, REAL *gl, REAL *gr, int *newt_check)
{
	int j;
	REAL df,dx,f,rtn,greater;

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
			tc_error("Jumped out of brackets in RTNEWT", 1);
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

REAL func_bis(REAL x, REAL t1, REAL t2, REAL t3, REAL *gl, REAL *gr)
{

	REAL f;

	f=t1*pow(x,gl[4])+t2*pow(x,gr[4])-t3;

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


REAL rtbis(REAL x1,REAL x2,REAL xacc, REAL t1, REAL t2, REAL t3, REAL *gl, REAL *gr)
{

	int j;
	REAL dx,f,fmid,xmid,rtb;

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
	tc_error("Too many bisections in rtbis", 1);
	return 0.0;
}

/*! Calculates various functions of gamma, the ratio of specific heats
 *
 * @param[in] gamma The value of c_p/c_v
 * @param[in] g Pointer to an array to contain the various functions of gamma
 */
void gamma_calc(REAL gamma, REAL *g)// just return all as the same, since we have a monatomic gas
{
	g[0]=gamma;
	g[1]=gamma;
	g[2]=gamma;
	g[3]=gamma;
	g[4]=gamma;
	g[5]=gamma;
	g[6]=gamma;
	g[7]=gamma;
	g[8]=gamma;
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

void adiflux(zone left_state, zone right_state, REAL dx, REAL dt, int perp, zone *resolved_state)
{
  
	extern REAL CFL, GAMMA[N_CHARGED_FLUIDS];
  /*
   * Local variables:-
   */
	int iter, newt_check, equal_gamma=1;
	REAL cl,cr,wleft,wright,u,v,w,c,ci,p,pi,g5, gl[9], gr[9], t1, t2, t3;
	REAL rhol, ul, pl, psil, rhor, ur, pr, psir;
	REAL density, psi, c_v;

/* First get our primitive variables */

	rhol=left_state.c[0];//rho left (density ρ)
	ul  =left_state.c[perp]/rhol;//speed left
	pl = pressure(left_state);
	psil =left_state.c[8];//PSI LEFT (ψ)

	rhor=right_state.c[0];//rho right (density ρ)
	ur  =right_state.c[perp]/rhor;//speed right
	pr = pressure(right_state);
	psir =right_state.c[8];//PSI RIGHT (ψ)

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

  gamma_calc(GAMMA[0], gl);
  gamma_calc(GAMMA[0], gr);
  //TODO use a #define to set GAMMA = to 5/3 as we are using a monotomic gas
  //or just set it to 1.6666 for ease and less accuracy

  cl = sqrt(gl[0]*pl*rhol);
  cr = sqrt(gr[0]*pr*rhor);
  wleft = -cl;
  wright = cr;

  //find value of the pressure
  p = (wright*pl - wleft*pr + wright*wleft*(ur - ul))/(wright - wleft);
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
//#if NONLINEAR_SOLVER
  if((fabs(p-pl)>0.1*pl || fabs(p-pr)>0.1*pr))  
  {
     /***************************************************************************
      *                                                                         *
      *                       Non-linear Riemann solver                         *
      *                                                                         *
      ***************************************************************************/

		if((p < pl) && (p < pr))
		{

      /*
       * Non-iterative solution for two rarefactions
       */
			if(equal_gamma)// Note that gl and gr are the same by def
			{
				ci = cr/(rhor*pow(pr,gr[4]));
				p = (0.5*gr[1]*(ul - ur) + cr/rhor + cl/rhol)/(ci + cl/(rhol*pow(pl,gr[4])));

				if(p<0.0) // Gas is expanding to form a vacuum (bad news!)
				{
					printf("Negative pressures in rarefaction solver: %lf\n",pow(-p,gr[5]));
					p=5.0e-11*(pr+pl); 
				}
				
				else 
					p = pow(p,gr[5]);
			}
			else // Have to do rarefactions with different gammas
			{

				/* Calculate the constants for the iteration */

				t1=2.0*cl/(rhol*gl[1]*pow(pl,gl[4]));//left going wave
				t2=2.0*cr/(rhor*gr[1]*pow(pr,gr[4]));//right going wave
				t3=ul-ur+2.0*((cr/(rhor*gr[1]))+(cl/(rhol*gl[1])));//if left going wave is heading right relative to gas

				p=mymin(pl,pr); // Simple first estimate for p

				p=newt(p,0.0,1.0e+10*pr, t1, t2, t3, gl, gr, &newt_check);
				if(newt_check) /* rtnewt has failed for function or #iter */
					p=rtbis(1.0e-10*pr,p,1.0e-04*pr, t1, t2, t3, gl, gr); /* Try bisection */

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
		wleft = -cl*wave(p,pl,gl);
		wright = cr*wave(p,pr,gr);
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
			wleft = -cl*wave(p,pl,gl);
			wright = cr*wave(p,pr,gr);
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
//#endif /* NONLINEAR_SOLVER */

  /*
   * Start calculation of resolved density and velocity
   */
  /* csos 
     if(p<1.0e-0){
     printf("resolved pressure < 0.0\n");
     }
  */

  wleft = -cl*wave(p,pl,gl);
  wright = cr*wave(p,pr,gr);
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
    v = vl;
	w = wl;
	Bx = Bxl;
	By = Byl;
	Bz = Bzl;
	psi = psil;

	g5 = gl[5];
	c_v = 1.0/(gl[0]-1.0);
    
    /*
     * Check velocity of left waves
     */
    if(p < pl) 
    {
      /*
       * Left wave is rarefaction
       */
      c = sqrt(gl[0]*p/(density));//speed of sound
      
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
	  u = gl[8]*(ul +gl[7]*cl);
	  p = pl*pow((u/cl),gl[5]);
	  density = gl[0]*p/(u*u);
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
    v = vr;
    w = wr;
	 Bx = Bxr;
	 By = Byr;
	 Bz = Bzr;
	 psi = psir;

	 g5 = gr[5];
	 c_v = 1.0/(gr[0]-1.0);
    
    /*
     * Check velocity of right waves
     */
    if(p < pr) 
    {
      /*
       * Right wave is rarefaction
       */
      c = sqrt(gr[0]*p/(density));//speed of sound
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
	  u = gr[8]*(ur - gr[7]*cr);
	  p = pr*pow((-u/cr),gr[5]);
	  density = gr[0]*p/(u*u);
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
	v = 0.5*(vl + vr);
	w = 0.5*(wl + wr);
	Bx = 0.5*(Bxl + Bxr);
	By = 0.5*(Byl + Byr);
	Bz = 0.5*(Bzl + Bzr);
	psi = 0.5*(psil + psir);
  }

/* One thing for Sam's suggestion for artificial dissipation - the
 * resolved sound speed */

	if(u > 0.0)
		c=sqrt(gl[0]*p/(density));
	else
		c=sqrt(gr[0]*p/(density));

  /*
   *  End of resolved pressure, velocities and density calculation
   *
   * Calculate fluxes in resolved state
   */

	(*resolved_state).c[0] = density;
	(*resolved_state).c[perp] = u;

	(*resolved_state).c[4] = p;
	(*resolved_state).c[8] = psi;
	(*resolved_state).c_v = c_v;
}
//#endif /* ADIABATIC */
