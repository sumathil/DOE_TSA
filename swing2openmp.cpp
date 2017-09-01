#include <sys/time.h>
#include "iostream"
#include "iomanip"
#include "cmath"
#include <stdio.h>
#include <omp.h>

using namespace std;

#define pi 3.14159265358979643846

// When fault occurs 
void Differentialduringfault(float *deltapresent,float *omegapresent,float deltaprevious,float omegaprevious,float omega0,float c_h)
{
	float temp,temp1,ddeltapresent,domegapresent,ddeltaprevious,domegaprevious;
	ddeltaprevious =omegaprevious-omega0;
	temp=deltaprevious+(c_h)*(ddeltaprevious);
	domegaprevious =((pi*60)/5)* (0.8-0.65*sin(temp));
	temp1 = omegaprevious+(c_h*(domegaprevious));
	ddeltapresent = temp1-omegaprevious;
	*deltapresent = deltaprevious + (c_h/2)*(ddeltaprevious+ddeltapresent);
	domegapresent =((pi*60)/5)* (0.8-(0.65*sin(*deltapresent)));
	*omegapresent = omegaprevious+(c_h/2)*(domegaprevious+domegapresent);

}

//When the fault is cleared
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
	float tint,tfin,omega0,*omega,*delta,*a,c_h,f_h,dint,*fine_delta,*fine_omega,*del_fine,*omega_fine;
	double t1,t2,t3,t4,et[110],tet[10],faverage,fsum;
	float *pred_delt,*corec_delt,*pred_omega,*corec_omega,delta_temp,omega_temp,*fine_tempd,*fine_tempo,*del_diff,*omega_diff,tempd,tempo;
	tint = 0.5;
	tfin = 102.9;
	dint = 26.39;
	cout<<"The coarse grid time step value is: "<<endl;
	cin>>c_h;
	cout<<"The fine grid step size value is: "<<endl;
	cin>>f_h;
	int num_steps = ((tfin-tint)/c_h)+1;
	int fine_size = ((tfin-tint)/f_h)+1;
	omega = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	delta = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	a = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	fine_delta =(float*)_mm_malloc((num_steps*sizeof(float)),32);
	fine_omega = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	del_fine = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	omega_fine = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	fine_tempd= (float*)_mm_malloc((fine_size*sizeof(float)),32);
	fine_tempo=(float*)_mm_malloc((fine_size*sizeof(float)),32);
	pred_delt= (float*)_mm_malloc((num_steps*sizeof(float)),32);
	pred_omega= (float*)_mm_malloc((num_steps*sizeof(float)),32);
	corec_delt=(float*)_mm_malloc((num_steps*sizeof(float)),32);
	corec_omega = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	del_diff =(float*)_mm_malloc((num_steps*sizeof(float)),32);
	omega_diff = (float*)_mm_malloc((num_steps*sizeof(float)),32);
	num_steps= num_steps-1;
	fine_size = fine_size -1;
	omega0=2*pi*60;
	omega[0]=omega0;
	delta[0]=(dint*pi)/180;
	cout<<"The value in radians is: "<<delta[0]<<endl;
	a[0]=tint;
	for(int k=0;k<2;k++)
	{
		tet[k]=0;
		gettimeofday(&start,NULL);
			if(k==0)
			{
				for(int i=0;i<num_steps;i++)
				{
					if(a[i]<0.8)
					{
						a[i+1]=a[i]+c_h;
						Differentialduringfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					}
					if(a[i]>=0.8)
					{
						a[i+1]=a[i]+c_h;
						Differentialpostfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					}
				}

			}
			else
			{
				for(int i=k;i<num_steps;i++)
				{
					delta[i]=corec_delt[i];
	                                omega[i]=corec_omega[i];

					if(a[i]<0.8)
					{
						a[i+1]=a[i]+c_h;
						Differentialduringfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					}
					if(a[i]>=0.8)
					{
						a[i+1] = a[i]+c_h;
						Differentialpostfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					}

				}
			}
		gettimeofday(&end,NULL);
		tet[k] = ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec));
		for(int i=0;i<110;i++)
		{
			t1 = omp_get_wtime();
			int n = 256;//omp_get_max_threads();
#pragma vector aligned
			{
#pragma omp parallel for private(tint,tfin) shared(delta,omega,fine_omega,fine_delta,fine_tempd,fine_tempo,a) num_threads(n)
				for(int j=0;j<num_steps;j++)
				{
					float fine_step;
					tint=a[j];
					tfin=a[j+1];
					fine_delta[j]=delta[j];
					fine_omega[j]=omega[j];
					int umax=round((tfin-tint)/f_h);
					if(a[j]<0.8)
					{
						for(int u=0;u<umax;u++)
						{
							fine_step = tint+f_h;
							Differentialduringfault(&fine_tempd[umax*j+u],&fine_tempo[umax*j+u],fine_delta[j],fine_omega[j],omega0,f_h);
							fine_delta[j]=fine_tempd[umax*j+u];
							fine_omega[j]=fine_tempo[umax*j+u];
							tint = fine_step;
						}
					}
					if(a[j]>=0.8)
					{
						for(int u=0;u<umax;u++)
						{
							fine_step = tint+f_h;
							Differentialpostfault(&fine_tempd[umax*j+u],&fine_tempo[umax*j+u],fine_delta[j],fine_omega[j],omega0,f_h);
							fine_delta[j]=fine_tempd[umax*j+u];
							fine_omega[j]=fine_tempo[umax*j+u];
							tint = fine_step;
						}

					}
					del_fine[j+1] = fine_delta[j];
					omega_fine[j+1]=fine_omega[j];
					del_diff[j]=del_fine[j+1]-delta[j+1];
					omega_diff[j]=omega_fine[j+1]-omega[j+1];
				}
			}
			t2 = omp_get_wtime();
			et[i]=(t2-t1)*1e6;
		}
		faverage = 0.0f;
		fsum = 0.0f;
		for(int j=10;j<110;j++)
		{
			faverage+=et[j];
		}
		fsum=faverage/100;
		cout<<"Fine grid OpenMP execution time is "<<fsum/1000<<" ms"<<endl;
		pred_delt[0]=del_fine[1];
		pred_omega[0]=omega_fine[1];
		gettimeofday(&start,NULL);
		for (int i=1;i<num_steps;i++)
		{
			a[i+1]=a[i]+c_h;
			if(a[i+1]<0.8)
			{
				Differentialduringfault(&pred_delt[i],&pred_omega[i],pred_delt[i-1],pred_omega[i-1],omega0,c_h);
			}
			if(a[i+1]>=0.8)
			{
				Differentialpostfault(&pred_delt[i],&pred_omega[i],pred_delt[i-1],pred_omega[i-1],omega0,c_h);
			}

		}
		for (int i=0;i<num_steps;i++)
		{
			corec_delt[i+1] = del_diff[i]+pred_delt[i];
			corec_omega[i+1] = omega_diff[i]+pred_omega[i];
		}
		for(int i=0;i<num_steps;i++)
		{
			if(abs(delta[i+1]-corec_delt[i+1])<1e-6)
			{
				cout<<"converged for "<<i+1<<"\telement"<<endl;
				break;
			}
		}
		gettimeofday(&end,NULL);
		tet[k]+=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec));
		cout<<"Sequential Execution time is "<< tet[k]/1000<<endl;
	}
	_mm_free(omega);
	_mm_free(delta);
	_mm_free(a);
	_mm_free(fine_delta);
	_mm_free(fine_omega);
	_mm_free(del_fine);
	_mm_free(omega_fine);
	_mm_free(fine_tempd);
	_mm_free(fine_tempo);
	_mm_free(pred_delt);
	_mm_free(pred_omega);
	_mm_free(corec_delt);
	_mm_free(corec_omega);
	_mm_free(del_diff);
	_mm_free(omega_diff);

}




