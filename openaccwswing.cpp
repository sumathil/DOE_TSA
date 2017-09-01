#include <sys/time.h>
#include "iostream"
#include "iomanip"
#include "cmath"
#include <stdio.h>
#include "openacc.h"
using namespace std;

#define pi 3.14159265358979323846
// When the fault is occured
#pragma acc routine seq
void Differential(double *deltapresent,double *omegapresent,double deltaprevious, double omegaprevious,double omega0,double c_h)
{
	double temp,temp1,ddeltapresent,domegapresent,ddeltaprevious,domegaprevious;
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
#pragma acc routine seq
void Differentialpostfault(double *deltapresent,double *omegapresent,double deltaprevious, double omegaprevious,double omega0,double c_h)
{
	double temp,temp1,ddeltapresent,domegapresent,ddeltaprevious,domegaprevious;
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
	double tint,tfin,omega0,*omega,*delta,*a,c_h,f_h,dint,*fine_delta,*fine_omega,*del_fine,*omega_fine,t1,t2,et[110],tet[110],faverage,fsum[10], fMemcpy[110],fmem,fmemsum[10],cet[110];
	double *pred_delt,*corec_delt,*pred_omega,*corec_omega,delta_temp,omega_temp,*fine_tempd,*fine_tempo,*del_diff,*omega_diff;
	tint=0.5;
	tfin=102.9;
	dint=26.39;
	c_h = 0.01;
	f_h = 0.00001;
	/*cout<<"The initial time value is : "<<endl;
	cin>>tint;
	cout<<"The final time value is: "<<endl;
	cin>>tfin;
	cout<<"The coarse grid time step value is: "<<endl;
	cin>>c_h;
	cout<<"The fine grid step size value is: "<<endl;
	cin>>f_h;
	cout<<"Enter the intial value of delta in degrees: "<<endl;
	cin>>dint;*/
	int num_steps = ((tfin-tint)/c_h)+1;
	int fine_size = ((tfin-tint)/f_h)+1;
	omega = new double[num_steps];
	delta = new double[num_steps];
	a = new double[num_steps];
	fine_delta =  new double[num_steps];
	fine_omega = new double[num_steps];
	del_fine = new double[num_steps];
	omega_fine = new double[num_steps];
	fine_tempd= new double[fine_size];
	fine_tempo=new double[fine_size];
	pred_delt= new double[num_steps];
	pred_omega= new double[num_steps];
	corec_delt=new double[num_steps];
	corec_omega = new double[num_steps];
	del_diff = new double[num_steps];
	omega_diff = new double[num_steps];
	num_steps= num_steps-1;
	fine_size = fine_size -1;
	omega0=2*pi*60;
	omega[0]=omega0;
	delta[0]=(dint*pi)/180;
	cout<<"The value in radians is: "<<delta[0]<<endl;
	a[0]=tint;
	//gettimeofday(&start,NULL);
	for(int k=0;k<2;k++)
	{
		tet[k]=0;
		if(k==0)
		{	for (int j=0;j<110;j++)
			{	
				gettimeofday(&start,NULL);
				for (int i=0;i<num_steps;i++)
				{
					a[i+1]=a[i]+c_h;
					//cout<<a[i+1]<<endl;
					if(a[i+1]<=0.8)
					{
						//cout<<a[i+1]<<endl;	//a[i+1]=a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
						Differential(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
						//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;

					}
					//cout<<"Break "<<endl;
					if(a[i+1]>0.8)
					{
						//cout<<"Break " <<endl;	//a[i+1]=a[i]+c_h;
						Differentialpostfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
						//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
					}

				}
				gettimeofday(&end,NULL);
				tet[j]=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
			}
		}
		else
		{
			for(int i=1;i<num_steps;i++)
			{
				delta[i]=corec_delt[i];
				omega[i]=corec_omega[i];
				a[i+1]=a[i]+c_h;
				//cout<<a[i+1]<<endl;	
				if(a[i+1]<=0.8)
				{
					//cout<<a[i+1]<<endl;	//	a[i+1]=a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
					Differential(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					//	cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
				}
				if(a[i+1]>0.8)
				{
					//a[i+1]=a[i]+c_h;
					Differentialpostfault(&delta[i+1],&omega[i+1],delta[i],omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
				}

			}
		}

		//gettimeofday(&end,NULL);
		//tet[k]=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		//	t1=omp_get_wtime();
		//	int n = omp_get_max_threads();
		//	#pragma omp parallel for private(tint,tfin) shared(delta,omega,fine_omega,fine_delta,fine_tempd,fine_tempo) num_threads(n) 
		for(int i=0;i<110;i++)
		{
			gettimeofday(&start,NULL); 
#pragma acc parallel num_workers(4) vector_length(256) copyin(delta[0:num_steps+1],omega[0:num_steps+1],a[0:num_steps+1]) copy(fine_delta[0:num_steps+1],fine_omega[0:num_steps+1],fine_tempd[0:fine_size],fine_tempo[0:fine_size]) copyout(del_fine[0:num_steps+1],omega_fine[0:num_steps+1],del_diff[0:num_steps+1],omega_diff[0:num_steps+1])
			//gettimeofday(&end,NULL);
			//fMemcpy[i]= ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;			//#pragma acc parallel
			//gettimeofday(&start,NULL);
			{
#pragma acc loop independent gang worker
				for(int j=0;j<num_steps;j++)
				{
					if(a[j]<0.8)
					{
						tint=a[j];
						tfin=a[j+1];
						int umax=(tfin-tint)/f_h;
						//	printf("%d\n ",umax);
						fine_delta[j]=delta[j];
						fine_omega[j]=omega[j];
						int k = umax*j;
						//cout << " less than 0.8 fine_size" <<fine_size<< " umax = " << umax <<endl;
						for (int u=0;u<umax;u++)
						{
							double fine_step = tint+f_h;
							Differential(&fine_tempd[k+u],&fine_tempo[k+u],fine_delta[j],fine_omega[j],omega0,f_h);
							fine_delta[j]=fine_tempd[k+u];
							fine_omega[j]=fine_tempo[k+u];
							tint=fine_step;
							//              printf("%lf\n , %d for time fine_step %lf\n",(fine_delta[j]*180)/pi,omp_get_thread_num(),fine_step);


						}
						del_fine[j+1]=fine_delta[j];
						omega_fine[j+1]=fine_omega[j];
						del_diff[j]=del_fine[j+1]-delta[j+1];
						omega_diff[j]=omega_fine[j+1]-omega[j+1];
						//	cout<<"The fine grid Delta funtion is " <<(fine_delta[j]*180)/pi<<" for the time step "<< tfin<<endl;
					}
					//when fault is cleared
					if(a[j]>=0.8)
					{
						tint=a[j];
						tfin=a[j+1];
						int umax=((tfin-tint)/f_h);
						//printf("%d\n",umax);
						fine_delta[j]=delta[j];
						fine_omega[j]=omega[j];
						//int k = umax*j;
						//cout << "greater tha 0.8 fine_size" <<fine_size<< " umax = " << umax <<"tint " <<tint<< "tfin " << tfin <<"   " << (tfin-tint)/f_h <<endl;
						for (int u=0;u<umax;u++)
						{
							double fine_step = tint+f_h;
							Differentialpostfault(&fine_tempd[k+u],&fine_tempo[k+u],fine_delta[j],fine_omega[j],omega0,f_h);
							fine_delta[j]=fine_tempd[k+u];
							fine_omega[j]=fine_tempo[k+u];
							tint=fine_step;
							//                printf("%lf\n , %d for time fine_step %lf\n",(fine_tempd[umax*j+u]*180)/pi,omp_get_thread_num(),fine_step);

						}
						del_fine[j+1]=fine_delta[j];
						omega_fine[j+1]=fine_omega[j];
						del_diff[j]=del_fine[j+1]-delta[j+1];
						omega_diff[j]=omega_fine[j+1]-omega[j+1];
						//cout<<"The fine grid Delta funtion is " <<(fine_delta[j]*180)/pi<<" for the time step "<< tfin<<endl;
					}
					//      cout<<"fine_delta[j] "<<(temp3*180)/pi<<"for time step "<<tfin<<"for array value ["<<j<<"]"<<j<<endl;
					//cout<<"The fine grid Delta funtion is " <<temp3<<endl;
					//cout<<"The omega function of coarse grid is " <<temp4<<endl;
				}
				//gettimeofday(&end,NULL);
			}
			gettimeofday(&end,NULL);
			et[i]=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		}
		faverage=0;
		fmem=0;
		for(int i=10;i<110;i++)
		{
			faverage+=et[i];
			fmem+=tet[i];
		}
		fsum[k]=faverage/100;
		fmemsum[k]=fmem/100;
		cout<<"Coarse grid time is "<<fmemsum[k]<<endl;
		cout<<"Fine grid time is "<<fsum[k]<<" ms"<<endl;
		gettimeofday(&start,NULL);
		for(int j=0;j<110;j++)
		{
			cet[j]=0;
			pred_delt[0]=del_fine[1];
			pred_omega[0]=omega_fine[1];
			for (int i=0;i<num_steps;i++)
			{
				a[i+1]=a[i]+c_h;
				if(a[i+1]<=0.8)
				{
					//a[i+1]=a[i]+c_h;
					Differential(&pred_delt[i+1],&pred_omega[i+1],pred_delt[i],pred_omega[i],omega0,c_h);
				}
				if(a[i+1]>0.8)
				{
					//a[i+1]=a[i]+c_h;
					Differentialpostfault(&pred_delt[i+1],&pred_omega[i+1],pred_delt[i],pred_omega[i],omega0,c_h);
				}

			}
			for (int i=0;i<num_steps;i++)
			{
				corec_delt[i+1] = del_diff[i]+pred_delt[i];
				corec_omega[i+1] = omega_diff[i]+pred_omega[i];
				//cout<< "The corrected grid values are "<< (corec_delt[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<endl;
			}
			gettimeofday(&end,NULL);
			cet[j]=((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		}

		double sum=0;
		for (int j=10;j<110;j++)
		{
			sum+=cet[j];
			cout<<sum<<endl;
		}
		double avg = 0;
		avg =  sum/100;
		cout<<"The sequential execution time is "<<avg<<" ms"<<endl;
	}
	// cout<<"The elapsed time for "<<k<<" iteration is "<<tet[k]<<" ms"<<endl;

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
