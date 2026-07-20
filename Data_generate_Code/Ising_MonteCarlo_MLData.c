  /*  program to simulate a 2 dim ising model    */
  /*     using the Metropolis algorithm          */
  /* Witten by Debasih Banerjee, Jyotirmoy Bhattacharya (May 2021) */
/*****************************************************************/

/*****************************************************************/

//HEADER FILES

#include<stdio.h>
#include<math.h>



//GLOBAL CONSTANTS
#define L 40 /* Lattice Size */
#define N L*L
#define J 1 /* Nearest Neighbour (NN) Coupling strength */
#define G 0.0 /* Next to Nearest Neighbour (NNN) Coupling strength */
#define equil 50000
#define meas  10000
#define sweep 1000
#define TemC 2.269

//GLOBAL VARIABLES 
int spin[N],fx[N],bx[N],fy[N],by[N],dfx[N],dbx[N],dfy[N],dby[N];
float Tempu;
double beta;
float mag;
float magavg;
char datafilename[50],magfilename[50];	



/******************************************************/
/***Random Number generating functions***/

/*random no generator:ran3+rani3*/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

long idum=-137456;
static int inext,inextp;
static long ma[56];

void ran3i(long *idum)
{
  int iff;
  long mj,mk;
  int i,ii,k;
  
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
}

float ran3()
{
  //static int inext,inextp;
  //static long ma[56];
  //static int iff;
  long mj,mk;
  int i,ii,k;
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}



/******************************************************/





// FUNTION TO SET THE INITIAL SPIN CONFIGURATION
void setspin()
{
  int i;
  for(i=0;i<N;i++)
    spin[i]=1;
}


// FUNCTION TO DEFINE BOUNDARY CONDITIONS
void boundarycond()
{
  int i,x,y;
  i=0;
  for(y=0;y<L;y++){
    for(x=0;x<L;x++){
           i= y*L + x;
           // forward neighbours
           fy[i] = ((y+1)%L)*L + x;
           fx[i] = y*L + ((x+1)%L);
	   // backward neigbours
           by[i] = ((y-1+L)%L)*L + x;
           bx[i] = y*L + ((x-1+L)%L);
           
             
           //Diagonal neighbours             
           dfy[i] = ((y+1)%L)*L + ((x+1)%L);
           dfx[i] = ((y+1)%L)*L + ((x-1+L)%L);	  
           dby[i] = ((y-1+L)%L)*L + ((x-1+L)%L);
           dbx[i] = ((y-1+L)%L)*L + ((x+1)%L);
                     
     }
  }
}


//FUNCTION IMPLEMENTING THE METROPOLIS ALGORITHM
void metropolis()
{
  int h,i,hnn;
  double ei,ef,de;
  double prob,rand_no;
  

  //Initialize magnetization 
  mag=0.0;
  
  
  for(i=0;i<N;i++)
    {      
      //Neighbouring spins 
      h=spin[fx[i]]+spin[bx[i]]+spin[fy[i]]+spin[by[i]];
      hnn=spin[dfx[i]]+spin[dbx[i]]+spin[dfy[i]]+spin[dby[i]];
     
      //previous energy
      ei= -J*spin[i]*h-G*spin[i]*hnn; 
      //new energy after flip 
      ef= -J*(-spin[i])*h-G*(-spin[i])*hnn;
      //energy difference 
      de= ef-ei;
      
      //energy=lookup[-e+4]; Defining probability
      prob= exp(-beta*de);   
      
      //Spin flip decission
      if (de<0 || de==0) {spin[i]=-spin[i];}
      else{  
	    rand_no=ran3(); /*JB: This works more efficiently than Mersenne twister which was also experimented with */
	    if(rand_no<prob) spin[i]=-spin[i];
	  }
	  
	  
      //Calculate magnetization of the resulting configuration
      mag=mag+spin[i];
    }
   mag = mag/((double)N);

//Magnetization square and 4th power can be specified at this stage 
   
}



//MAIN PROGRAM RECORDING DATA FOR A PARTICULAR TEMPERATURE 
int recordDATA()
{ 

int i,k,kk;
double magsum=0.0;
FILE *fptr1,*fptr2,*fptrFD,*fptrDL;
 
// SET INITIAL SPIN CONFIGURATION 
  setspin();
  
// SET BOUNDARY CONDITION 
  boundarycond();
  


//Run metropolis till equilibriation time
for(k=0;k<equil;k++)
metropolis();
 
 
/***********************************/
//Run for data collections:
/***********************************/



// Open files for recording data
fptr1=fopen(magfilename,"a");
fptrFD=fopen("G0_dataCombi/fulldataJBset1.csv","a");
fptrDL=fopen("G0_dataCombi/datalabelJBset1.csv","a");
fptr2=fopen(datafilename,"a");
//Data collection run  
for(k=1;k<=meas;k++)
{ 
for(kk=0;kk<sweep;kk++)
{
      metropolis();
}
    
//Record spins  
for(i=0;i<N;i++)
{
if (i!=(N-1))
{fprintf(fptr2,"%d,",spin[i]);
fprintf(fptrFD,"%d,",spin[i]);}
else
{fprintf(fptr2,"%d\n",spin[i]);
fprintf(fptrFD,"%d\n",spin[i]);}
}


if (Tempu<=TemC) fprintf(fptrDL,"1\n");
else fprintf(fptrDL,"0\n");

//Record magnetization
if (k!=meas) fprintf(fptr1,"%f, ",mag);   
else fprintf(fptr1,"%f};\n",mag);
magsum=magsum+(pow(mag,2));   
}

magavg=magsum/((double)meas);


fclose(fptrFD);
fclose(fptrDL);
fclose(fptr1);
fclose(fptr2);
  
   
return 0;
}





//RUN THE MAIN PROGRAM FOR DIFERENT TEMPERATURES
int main()
{

ran3i(&idum);
int k;
double Tc,Ta;
FILE *fptr3,*fptr4,*fptr5,*fptr6,*fptr7;

  

/**************************************************************************/
/*JB: Input the value	 of beta instead of taking input from the file */
double Ts,Te,Ti;
printf("Input the range of value of T: \n Input starting value: ");
scanf("%le",&Ts);
printf("\n Input ending value: ");
scanf("%le",&Te);
printf("\n Input interval: ");
scanf("%le",&Ti);

/**************************************************************************/


//Prepare file to record magnetization data
sprintf(magfilename,"magDATAset1");
fptr3=fopen(magfilename,"w");
fclose(fptr3);

//Prepare full data file and labels file 
fptr5=fopen("G0_dataCombi/fulldataJBset1.csv","w");
fptr6=fopen("G0_dataCombi/datalabelJBset1.csv","w");
/*******/
// print temperuature on top of the combined datafile
for(Ta=Ts;Ta<=Te;Ta=Ta+Ti)
{fprintf(fptr5,"%f,",Ta);
}
fprintf(fptr5,"\n");
/*******/
fclose(fptr5);
fclose(fptr6);

//Prepare file to record average magnetization vs T data
fptr4=fopen("MagvsT_set1","w");
fprintf(fptr4,"MagvsT={");
//We are not closing it here.



//Run the main program for different values of temperature
k=1;
for(Tc=Ts;Tc<=Te;Tc=Tc+Ti)
{Tempu=Tc;
printf("\n %le", Tempu);
/**************************************************************************/
beta = 1/Tempu;
/**************************************************************************/
/*JB: File name creation using  T */
sprintf(datafilename,"G0_data/DATA_L40_M10000_T%f",Tempu);
printf("\n %s",datafilename);

fptr7=fopen(datafilename,"w");
fprintf(fptr7,"%f\n",Tempu);
fclose(fptr7);
/**************************************************************************/
// magnetization data file: Input for different temperature
fptr3=fopen(magfilename,"a");
fprintf(fptr3,"(*Magnetization data for T=%f*)\n",Tempu);
fprintf(fptr3,"magvals%d={",k);
k=k+1;
fclose(fptr3);
/**************************************************************************/
recordDATA();

//Record average magnetization vs T data
if (Tc<(Te-Ti)) fprintf(fptr4,"{%f, %f}, ",Tempu, magavg);
else fprintf(fptr4,"{%f, %f}};\n",Tempu, magavg);
}

fclose(fptr4);	


return 0;
}

