#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
double time =0.0;
double timeto =1.0e-07;


//double  speed,timeou;
double  dx,dt;
double u[102], flux[102];
struct initidata
 {  
 int nfreq,ntmaxi,cells;
 double cfl,domlen,wavespeed,timeou;
 		//double *U;
 };
struct initidata s1={10,1000000 ,100,0.9,1.0,1.0,0.5 };


//void initialcondition(double *domlen, int itest, int *cells) ;
void initialcondition(struct initidata s1) ;
void cflcondition(struct initidata s1);
void boundary(struct initidata s1);
void update(struct initidata s1);
void solution(struct initidata s1);
int riemann(double ul,double ur,double ustar);
void godunovflux(struct initidata s1);


/* *******************************************printf(" inticond%10.5f\n",dx);********** */
void main()
/* ***************************************************** */
  {
  int i,n;
 
 // double  speed, timeou;
  //itest=1;
  //cells=100;
  //speed=1.0;
 // timeou=1.0;
  //time =0.0;
  //timeto =1.0e-07;
    
//            initialcondition(&s1.domlen, itest, &s1.cells);

  printf("nfreq=%d ntmaxi=%d cells=%d cfl=%g domlen=%g wavespeed=%g timeou=%g\n",
					s1.nfreq,s1.ntmaxi,s1.cells,s1.cfl,s1.domlen,s1.wavespeed,s1.timeou);

	initialcondition (s1);  


           // printf("\n----------------------------------------\n");
            //printf("\n  time step n        time \n");   
           // printf("\n----------------------------------------\n");
           
	
          
 for(i=1;i <= s1.ntmaxi; i++){
// for(i=1;i <= s1.cells; i++){  
 // call function for calculating boundary,impose CFL condition and Calculate Godunov flux
     
		boundary(s1);
		printf("i=%d\ttime=%g\n", i, time);
		
	  cflcondition(s1);
	  time = time + dt;                                         
	  godunovflux(s1);     
	  update(s1);
  
		if(fabs(time-s1.timeou) < timeto)
		{
			
        solution(s1);	    
			printf("\n----------------------------------------\n");
			printf("\n  number of time steps = %5d\n",i);   
			printf("\n----------------------------------------\n");
            
			break;
		}
	} 

	  
}
/* ***************************************************** */
void initialcondition(struct initidata s1)
/* ***************************************************** */  
{
	int i,idimension=s1.cells+2;
	double  xleft, xpos, xright;
 
	dx = (s1.domlen)/(double)s1.cells;
 
//	printf("dx=%g\n",dx);    

	for(i=0; i < idimension; i++)
	{
		u[i]    = 0.0;
		flux[i] = 0.0;
	}
   
	xpos = -1.0;
 
	for(i=1; i <= s1.cells; i++)
	{
		xpos = xpos + 2.0/(double)(s1.cells);
		u[i] = exp(-8.0*xpos*xpos);
		
//		printf("u[%d]=%g\n",i,u[i]);
	}      
}                 
/* ***************************************************** */
void boundary(struct initidata s1)
/* ***************************************************** */
{
	u[0]=u[s1.cells];
	u[s1.cells+1]=u[1];

//	printf("u[0]=%g u[%d]=%g\n",u[0],s1.cells+1,u[s1.cells+1]);                    
}
/* ***************************************************** */
void cflcondition(struct initidata s1)
/* ***************************************************** */   
{ int i;
	double smax;
//	printf("\n s1.wavespeed= %g s1.cfl=%g dx=%g\n",s1.wavespeed,s1.cfl,dx);
	//smax = fabs(s1.wavespeed);
	smax=-1.0e+06;
	for(i=0;i <= s1.cells+1;i++)
	{
	if(fabs(u[i]) > smax)
	   smax=fabs(u[i]);
	}
	
	dt= s1.cfl*dx/smax;

//	printf("\n smax=%g dt=%g\n",smax,dt); for(i=1;i <= iter; i++){

	if((time+dt) > s1.timeou)
	{  
  	dt = s1.timeou - time;
	}
  return dt;                  
  
      
//	printf("\n dt=%g  time  = %10.5f\n", dt, time );   
}
/* ***************************************************** */
void update(struct initidata s1)
/* ***************************************************** */      
{
 int i;
 double deltax = dt/dx;

 printf("deltax=%g\n", deltax);
      
	for(i=1; i <= s1.cells; i++)
	{      
		u[i] = u[i] + deltax*(flux[i-1]-flux[i]);
//		printf("u[%d]=%g\n",i,u[i]);     
	}
}
/* ***************************************************** */
void solution(struct initidata s1)
/* ***************************************************** */          
{
	int i;
	double xpos;
	FILE * Output;
  Output = fopen("burger.out", "a");  
	for(i=1; i <= s1.cells; i++){
      
		xpos = ((double)i-0.5)*dx;
		fprintf(Output,"%g  %g\n",xpos,u[i]);        
	}
	 fclose(Output); 
}
/* ***************************************************** */  
void  godunovflux(struct initidata s1)
/* ***************************************************** */     
{
	int i;
	double  ul,ur,ustar;
          
//	printf("line=%d wavespeed\n",__LINE__,s1.wavespeed);
      
	for(i=0; i <= s1.cells; i++)
	{
	
	ul=u[i];
	ur=u[i+1];
	riemann (ul, ur,ustar);
	
	/*fright = 0.5*(1.0 + 2*s1.cfl)*(s1.wavespeed)*u[i];
	fleft = 0.5*(1.0 - 2*s1.cfl)*(s1.wavespeed)*u[i+1];

		flux[i] = fright + fleft;
		*/
//		printf ("flux[%d]=%g\n",i,flux[i]);
    flux[i]=0.5*(ustar*ustar);
    //printf ("flux[%d]=%g\n",i,flux[i]);
	}                     
}
/* ***************************************************** */                     
                     
int riemann(double ul,double ur,double ustar)
/* ***************************************************** */ 
{
    double s;
    
    if(ul > ur)
    {
      s=0.5*(ul+ur);
    
      if(s >= 0.0)
      ustar=ul;    
      else
      ustar=ur;
      }
    else{
       if(ul >= 0.0)
       ustar=ul;
       
       if(ur < 0.0)
       ustar=ur;
       
       if(ul < 0.0 && ur >=  0.0)
       ustar=0.0;
       }
       
 }
/* ***************************************************** */     
