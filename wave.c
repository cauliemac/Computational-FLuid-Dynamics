/*
 * This function calculates the wave speed for a wave connecting states with
 * pressures pi, p ahead an behind respectively
 *
 * Arguments:-
 * Pressure in upstream state, pi
 * Pressure in downstream state, p
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "typedefs.h"
#include "functions.h"

/*! Finds wave speed between two pressure states.
 *
 * @param[in] p Pressure of one of the states
 * @param[in] pi Pressure of the other state
 * @param[in] g Pointer to array of functions of gamma
 */

double wave(double p,double pi, double gamma)
{

	int i;
  /*
   * Local variables
   */
	double ratio,waveval;

	if(p==1.0e-10)
		waveval = gamma;

	else{
		ratio = p/pi;
  
		if(fabs(ratio-1.0) < 1.0e-03) /* Use linear expression */
			waveval = 1.0 + 0.5*gamma*(ratio - 1.0);

		else{ /* Use non-linear expression */

			if(ratio >= 1.0) /* wave is a shock */
				waveval = sqrt(1.0 + gamma*(ratio - 1.0));

			else /* wave is a rarefaction */
				waveval = gamma*(1.0 - ratio)/(1.0 - pow(ratio,gamma));
		}
	}
	return waveval;
}
