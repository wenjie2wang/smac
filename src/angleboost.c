#include <R.h>
#include <Rmath.h>
#include <stdlib.h>

/*******************************************************************************/
void firboost(double* u, int length, double* z)
{ 
int i;
for(i=0;i<length;i++)
{
*z = -exp(-u[i]);
z++;
}
}


/*******************************************************************************/
void secboost(double* u, int length, double* z)
{ 
int i;
for(i=0;i<length;i++)
{
*z = exp(-u[i]);
z++;
}
}


/*******************************************************************************/

void angleboost(double* x, double* y, int* kkminus, int* nnobs, int* nnp, double* lam, double* epsi, double* w, double* warmbeta, double* warminner, double* betaout, double* innerout)
{ 

int i,z,numiter,q,p,j;

int kminus = *kkminus;

int nobs = *nnobs;

int np = *nnp;

double lambda = *lam;

double epsilon= *epsi; /*line 41*/ 

double inner[nobs];

double betacheck1[((np+1)*kminus)];

double beta[((np+1)*kminus)];

double zeroinner[nobs];

double diff;

double partial, secpartial, zeropartial;

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
		firboost(inner,nobs,tempfd);
			for (i=0;i<nobs;i++)
			{
			partial += w[i]*tempfd[i]*y[(q*nobs+i)]; 
			}
		partial/=nobs;
		if (fabs(partial)<epsilon) {break;}
			
		secpartial=0;
		secboost(inner,nobs,tempsd);
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
	
	for (q=0;q<kminus;q++)
	{
	/*= update beta =========*/

		for (p=1;p<(np+1);p++)
		{
		/* basically, update beta[(q*(np+1)+p)] */
			
		zeropartial=0;
			for (i=0;i<nobs;i++)
			{
			zeroinner[i] = inner[i] - beta[(q*(np+1)+p)]*y[(q*nobs+i)]*x[((p-1)*nobs+i)];
			}
		firboost(zeroinner,nobs,tempfd);
			for (i=0;i<nobs;i++)
			{
			zeropartial += w[i]*tempfd[i]*y[(q*nobs+i)]*x[((p-1)*nobs+i)]; 
			}
		zeropartial/=nobs;
				
			/* soft update part */
			if (zeropartial > lambda)
				{
					if (beta[(q*(np+1)+p)] > 0)
						{
							for (i=0;i<nobs;i++)
							{
							inner[i] = zeroinner[i];
							}
							beta[(q*(np+1)+p)]=0;
						}
					
					for (j=0;j<100;j++)
					{

					partial=0;
					firboost(inner,nobs,tempfd);
						for (i=0;i<nobs;i++)
						{
						partial += w[i]*tempfd[i]*y[(q*nobs+i)]*x[((p-1)*nobs+i)]; 
						}
					partial=partial/nobs-lambda;
					if (fabs(partial)<epsilon) {break;}

					secpartial=0;
					secboost(inner,nobs,tempsd);
						for (i=0;i<nobs;i++)
						{
				secpartial += w[i]*tempsd[i]*y[(q*nobs+i)]*y[(q*nobs+i)]*x[((p-1)*nobs+i)]*x[((p-1)*nobs+i)]; 
						}
					secpartial/=nobs;
						
						temp=partial/secpartial;
						beta[(q*(np+1)+p)] -= temp;
							for (i=0;i<nobs;i++)
							{
							inner[i] -= temp*y[(q*nobs+i)]*x[((p-1)*nobs+i)];
							}
					
					} /* for (j=0;j<100;j++) */ 
					
				} /* if (zeropartial > lambda) */

			if (zeropartial < (-lambda))
				{
				
					if (beta[(q*(np+1)+p)] < 0)
						{
							for (i=0;i<nobs;i++)
							{
							inner[i] = zeroinner[i];
							}
							beta[(q*(np+1)+p)]=0;
						}
					
					for (j=0;j<100;j++)
					{

					partial=0;
					firboost(inner,nobs,tempfd);
						for (i=0;i<nobs;i++)
						{
						partial += w[i]*tempfd[i]*y[(q*nobs+i)]*x[((p-1)*nobs+i)]; 
						}
					partial=partial/nobs+lambda;
					if (fabs(partial)<epsilon) {break;}

					secpartial=0;
					secboost(inner,nobs,tempsd);
						for (i=0;i<nobs;i++)
						{
				secpartial += w[i]*tempsd[i]*y[(q*nobs+i)]*y[(q*nobs+i)]*x[((p-1)*nobs+i)]*x[((p-1)*nobs+i)]; 
						}
					secpartial/=nobs;
						
						temp=partial/secpartial;
						beta[(q*(np+1)+p)] -= temp;
							for (i=0;i<nobs;i++)
							{
							inner[i] -= temp*y[(q*nobs+i)]*x[((p-1)*nobs+i)];
							}
					
					} /* for (j=0;j<100;j++) */ 
	
				} /* if (zeropartial < (-lambda)) */

			if (fabs(zeropartial) < lambda)
				{
					for (i=0;i<nobs;i++)
					{
					inner[i] = zeroinner[i];
					}
					beta[(q*(np+1)+p)]=0;
				} /* if (fabs(zeropartial) < lambda) */
		
			/* soft update part */
			
		} /* for (p=1;p<(np+1);p++) */

	/*= update beta =========*/

	


	} /* for (q=0;q<kminus;q++) */
	
/*- check difference and comapre to epsilon ------*/
	
	diff=0;
	for (z=0;z<((kminus*(1+np)));z++)
	{
	diff+=fabs(beta[z]-betacheck1[z]);
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

} /* void angleboost */
