#include <sys/time.h>
#include "iostream"
#include "iomanip"
#include "cmath"
#include <stdio.h>
#include<omp.h>
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
	float tint,tfin,omega0,*omega,*delta,*a,c_h,f_h,dint,*fine_delta,*fine_omega,*del_fine,*omega_fine;
	double t1,t2,t3,t4,et[110],tet[10],faverage,fsum[10];
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
	omega = new float[num_steps];
	delta = new float[num_steps];
	a = new float[num_steps];
	fine_delta =  new float[num_steps];
	fine_omega = new float[num_steps];
	del_fine = new float[num_steps];
	omega_fine = new float[num_steps];
	fine_tempd= new float[fine_size];
	fine_tempo=new float[fine_size];
	pred_delt= new float[num_steps];
	pred_omega= new float[num_steps];
	corec_delt=new float[num_steps];
	corec_omega = new float[num_steps];
	del_diff = new float[num_steps];
	omega_diff = new float[num_steps];
	num_steps= num_steps-1;
	fine_size = fine_size -1;
	omega0=2*pi*60;
	omega[0]=omega0;
	delta[0]=(dint*pi)/180;
	cout<<"The value in radians is: "<<delta[0]<<endl;
	a[0]=tint;
	for(int k=0;k<2;k++)
	{       tet[k]=0;
		gettimeofday(&start,NULL);
		if(k==0)
		{
			for (int i=0;i<num_steps;i++)
			{
				if(a[i]<0.8)
				{
					a[i+1]=a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
					Differential(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
				}
				if(a[i]>=0.8)
				{
					a[i+1]=a[i]+c_h;
					Differentialpostfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
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
					a[i+1]=a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
					Differential(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
				}
				if(a[i]>=0.8)
				{
					a[i+1]=a[i]+c_h;
					Differentialpostfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
				}
			}
		}
		gettimeofday(&end,NULL);
		tet[k]=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec));
		//cout<<"Sequential Execution time is "<< tet[k]<<endl;
		for(int i=0;i<110;i++)
		{
			t1=omp_get_wtime();
			int n=omp_get_max_threads();
#pragma omp parallel for private(tint,tfin) shared(delta,omega,fine_omega,fine_delta,fine_tempd,fine_tempo,tempd,tempo) num_threads(n)  
			for(int j=0;j<num_steps;j++)
			{
				tint=a[j];
				tfin=a[j+1];
				fine_delta[j]=delta[j];
				fine_omega[j]=omega[j];
				int umax=round((tfin-tint)/f_h);
				if(a[j]<0.8)
				{

					for (int u=0;u<umax;u++)
					{
						float fine_step = tint+f_h;
						Differential(&fine_tempd[umax*j+u],&fine_tempo[umax*j+u],fine_delta[j],fine_omega[j],omega0,f_h);
						fine_delta[j]=fine_tempd[umax*j+u];
						fine_omega[j]=fine_tempo[umax*j+u];
					//	fine_delta[j]=tempd;
					//	fine_omega[j]=tempo;
						tint=fine_step;
						//printf("%lf\n , %d for time fine_step %lf\n",(fine_tempd[umax*j+u]*180)/pi,omp_get_thread_num(),fine_step);
					}
					//del_fine[j+1]=fine_delta[j];
					//omega_fine[j+1]=fine_omega[j];
					//del_diff[j]=del_fine[j+1]-delta[j+1];
					//omega_diff[j]=omega_fine[j+1]-omega[j+1];
				}
				if(a[j]>=0.8)
				{
					for (int u=0;u<umax;u++)
					{
						float fine_step = tint+f_h;
						Differentialpostfault(&fine_tempd[umax*j+u],&fine_tempo[umax*j+u],fine_delta[j],fine_omega[j],omega0,f_h);
						fine_delta[j]=fine_tempd[umax*j+u];
                                                fine_omega[j]=fine_tempo[umax*j+u];
						//fine_delta[j]=tempd;
						//fine_omega[j]=tempo;
						tint=fine_step;
						//printf("%lf\n , %d for time fine_step %lf\n",(fine_tempd[umax*j+u]*180)/pi,omp_get_thread_num(),fine_step);
					}
					//del_fine[j+1]=fine_delta[j];
					//omega_fine[j+1]=fine_omega[j];
					//del_diff[j]=del_fine[j+1]-delta[j+1];
					//omega_diff[j]=omega_fine[j+1]-omega[j+1];

				}
				del_fine[j+1]=fine_delta[j];
				omega_fine[j+1]=fine_omega[j];
				del_diff[j]=del_fine[j+1]-delta[j+1];
				omega_diff[j]=omega_fine[j+1]-omega[j+1];
			}
			t2=omp_get_wtime();
			et[i]=(t2-t1)*1e6;
			//cout<<"The elapsed time is "<<et[i]<<" us"<<endl;
		}
		faverage=0;
		for(int i=10;i<110;i++)
		{
			faverage+=et[i];
		}
		fsum[k]=faverage/100;
		cout<<"Fine grid OpenMP execution time is "<<fsum[k]/1000<<" ms"<<endl;
		pred_delt[0]=del_fine[1];
		pred_omega[0]=omega_fine[1];
		gettimeofday(&start,NULL);
		for (int i=1;i<num_steps;i++)
		{
			a[i+1]=a[i]+c_h;
			if(a[i+1]<0.8)
			{
				Differential(&pred_delt[i],&pred_omega[i],pred_delt[i-1],pred_omega[i-1],omega0,c_h);
			}
			if(a[i+1]>=0.8)
			{
				Differentialpostfault(&pred_delt[i],&pred_omega[i],pred_delt[i-1],pred_omega[i-1],omega0,c_h);
			}

		}
		/*		gettimeofday(&end,NULL);
				tet[k]+=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec));
				cout<<"Sequential Execution time is "<< tet[k]/1000<<endl;
				for(int z=0;z<110;z++)
				{
				t3=omp_get_wtime();
#pragma omp parallel for*/
		for (int i=0;i<num_steps;i++)
		{
			corec_delt[i+1] = del_diff[i]+pred_delt[i];
			corec_omega[i+1] = omega_diff[i]+pred_omega[i];
		}
		for(int i=0;i<num_steps;i++)
		{
			if(abs(delta[i+1]-corec_delt[i+1])<1e-10)
			{
				cout<<"converged for "<<i+1<<"\telement"<<endl;
				break;
			}
			/*else
			 *                           {
			 *                                                     cout<<"Not converged"<<endl;
			 *                                                                               }*/
		}
		//		t4=omp_get_wtime();
		//		et[z]=(t4-t3)*1e6;
		//		}
		/*	faverage=0;
			for(int i=10;i<110;i++)
			{
			faverage+=et[i];
			}
			fsum[k]+=faverage/100;
			cout<<"Total OpenMP execution time is "<<fsum[k]/1000<<" ms"<<endl;*/
		gettimeofday(&end,NULL);
		tet[k]+=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec));
		cout<<"Sequential Execution time is "<< tet[k]/1000<<endl;
	}
	delete[] omega;
	delete[] delta;
	delete[] a;
	delete[] fine_delta;
	delete[] fine_omega;
	delete[] del_fine;
	delete[] omega_fine;
	delete[] fine_tempd;
	delete[] fine_tempo;
	delete[] pred_delt;
	delete[] pred_omega;
	delete[] corec_delt;
	delete[] corec_omega;
	delete[] del_diff;
	delete[] omega_diff;
}
