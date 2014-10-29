#include <R.h>
#include <Rmath.h>
#include <stdlib.h>

/*******************************************************************************/
void firpsvms(double* u, int length, double* z)
{ 
int i;
for(i=0;i<length;i++)
{
*z = -2*(1-u[i]);
z++;
}
}

/*******************************************************************************/
void secpsvms(double* u, int length, double* z)
{ 
int i;
for(i=0;i<length;i++)
{
*z = 2;
z++;
}
}

/*******************************************************************************/

void anglepsvmstart(double* x, double* y, int* kkminus, int* nnobs, int* nnp, double* epsi, double* w, double* warmbeta, double* warminner, double* betaout, double* innerout)
{ 

int i,z,numiter,q,j;

int kminus = *kkminus;

int nobs = *nnobs;

int np = *nnp;

double epsilon= *epsi;  

double inner[nobs];

double betacheck1[((np+1)*kminus)];

double beta[((np+1)*kminus)];

double diff;

double partial, secpartial;

double temp;

double tempfd[nobs],tempsd[nobs];

/*- prepare -------------------------------------------------------------------*/

for (i=0;i<nobs;i++) {inner[i]=warminner[i];}

for (z=0;z<((np+1)*kminus);z++)	
	{
	betacheck1[z] = warmbeta[z];
	beta[z] = warmbeta[z];
	}

/*- update --------------------------------------------------------------------*/

for (numiter=0;numiter<500;numiter++)

{

	for (q=0;q<kminus;q++)
	{
	/*- update beta0 --------*/
	/* basically, update beta[q*(np+1)] */
	
		for (j=0;j<100;j++)
		{
		
		partial=0;
		firpsvms(inner,nobs,tempfd);
			for (i=0;i<nobs;i++)
			{
			partial += w[i]*tempfd[i]*y[(q*nobs+i)]; 
			}
		partial/=nobs;
		if (fabs(partial)<epsilon) {break;}
			
		secpartial=0;
		secpsvms(inner,nobs,tempsd);
			for (i=0;i<nobs;i++)
			{
			secpartial += w[i]*tempsd[i]*y[(q*nobs+i)]*y[(q*nobs+i)]; 
			}
		secpartial/=nobs;
		
		temp=partial/secpartial;
		beta[q*(np+1)] -= temp;
			for (i=0;i<nobs;i++)
			{
			inner[i] -= temp*y[(q*nobs+i)];
			}

		} /* for (j=0;j<100;j++) */
	
	/*- update beta0 --------*/

	}
	
	
/*- check difference and comapre to epsilon ------*/
	
	diff=0;
	for (q=0;q<kminus;q++)
	{
	diff+=fabs(beta[q*(np+1)]-betacheck1[q*(np+1)]);
	}

	if (diff<epsilon) {break;}
		else	{
			for (z=0;z<((np+1)*kminus);z++)	
				{
				betacheck1[z] = beta[z];
				}
			}

} /* for (numiter=0;numiter<500;numiter++) */

/*- update --------------------------------------------------------------------*/

/* report the result */

for (i=0;i<nobs;i++) {innerout[i]=inner[i];}

for (z=0;z<((np+1)*kminus);z++)	
	{
	betaout[z] = beta[z];
	}

} /* void anglepsvm */
