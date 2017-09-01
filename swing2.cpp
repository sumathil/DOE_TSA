#include <sys/time.h>
#include "iostream"
#include "iomanip"
#include "cmath"
#include <stdio.h>
//#include<omp.h>
using namespace std;

#define pi 3.14159265358979323846
// When the fault is occured
void Differential(float *deltapresent,float *omegapresent,float deltaprevious, float omegaprevious,float omega0,float c_h)
{
	float temp,temp1,ddeltapresent,domegapresent,ddeltaprevious,domegaprevious;
	ddeltaprevious =omegaprevious-omega0;
	temp=deltaprevious+(c_h)*(ddeltaprevious);
	domegaprevious =((pi*60)/5)* (0.8-0.65*sin(temp));
	temp1 = omegaprevious+(c_h*(domegaprevious));
	ddeltapresent = temp1-omegaprevious;
	*deltapresent = deltaprevious + (c_h/2)*(ddeltaprevious+ddeltapresent);
	domegapresent =((pi*60)/5)* (0.8-(0.65*sin(*deltapresent)));
	//domegapresent = 32-(173.68*sin(*deltapresent-(10*pi)/180));
	*omegapresent = omegaprevious+(c_h/2)*(domegaprevious+domegapresent);
}

//Once the fault is cleared
void Differentialpostfault(float *deltapresent,float *omegapresent,float deltaprevious,float omegaprevious,float omega0,float c_h)
{
	float temp,temp1,ddeltapresent,domegapresent,ddeltaprevious,domegaprevious;
	ddeltaprevious =omegaprevious-omega0;
	temp=deltaprevious+(c_h)*(ddeltaprevious);
	domegaprevious =((pi*60)/5)* (0.8-1.4625*sin(temp));
	temp1 = omegaprevious+(c_h*(domegaprevious));
	ddeltapresent = temp1-omegaprevious;
	*deltapresent = deltaprevious + (c_h/2)*(ddeltaprevious+ddeltapresent);
	domegapresent =((pi*60)/5)* (0.8-1.4625*sin(*deltapresent));
	*omegapresent = omegaprevious+(c_h/2)*(domegaprevious+domegapresent);
}


int main()
{
	struct timeval start,end;
	float tint,tfin,omega0,*omega,*delta,*a,f_h,dint;
	float et[110],faverage,fsum;
	cout<<"The initial time value is : "<<endl;
	cin>>tint;
	cout<<"The final time value is: "<<endl;
	cin>>tfin;
	cout<<"The fine grid step size value is: "<<endl;
	cin>>f_h;
	cout<<"Enter the intial value of delta in degrees: "<<endl;
	cin>>dint;
	int num_steps = ((tfin-tint)/f_h)+1;
	omega = new float[num_steps];
	delta = new float[num_steps];
	a = new float[num_steps];
	num_steps = num_steps -1;
	omega0=2*pi*60;
	omega[0]=omega0;
	delta[0]=(dint*pi)/180;
	cout<<"The value in radians is: "<<delta[0]<<endl;
	a[0]=tint;
	for(int k=0;k<110;k++)
	{
		gettimeofday(&start,NULL);
		for (int i=0;i<num_steps;i++)
		{
			a[i+1]=a[i]+f_h;
			//cout<<a[i+1]<<endl;
			if(a[i+1]<=0.8)
			{
				//cout<<a[i+1]<<endl;	//a[i+1]=a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
				Differential(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,f_h);
				//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
				//cout<<(delta[i+1]*180)/pi<<endl;

			}
			//cout<<"Break "<<endl;
			if(a[i+1]>0.8)
			{
				//cout<<"Break " <<endl;	//a[i+1]=a[i]+c_h;
				Differentialpostfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,f_h);
				//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
				//cout<<(delta[i+1]*180)/pi<<endl;
			}

		}

		//      #pragma omp parallel for private(tint,tfin) shared(delta,omega,fine_omega,fine_delta,fine_tempd,fine_tempo) num_threads(n)  
		gettimeofday(&end,NULL);
		et[k]=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
	}
	faverage=0;
	for(int i=0;i<110;i++)
	{
		faverage+=et[i];
	}
	fsum=faverage/110;
	cout<<"The execution time is "<<fsum<<" ms"<<endl;
	delete[] omega;
	delete[] delta;
	delete[] a;
}
