#include "basic_include.h"

void DifferentialEquations(float *delta_present,float *omega_present,float *Edp_present, float *Eqp_present, float delta_previous,float omega_previous,float Edp_previous, float Eqp_previous, float dt, float omega0)
{

	float I_q, I_d, d_delta, d_omega, d_Edp, d_Eqp, Te;

	I_d=(Eb*cos(delta_previous) - Eqp_previous)/(x_dp+x_e);
	I_q=(Edp_previous+Eb*sin(delta_previous))/(x_qp+x_e);

	d_Eqp = 1/T_dp*(-Eqp_previous+ (x_d-x_dp)*I_d+Efd);
	*Eqp_present = Eqp_previous + dt*d_Eqp;

	d_Edp = 1/T_qp*(-Edp_previous-(x_q-x_qp)*I_q);	
	*Edp_present = Edp_previous +dt*d_Edp;

	Te = Eqp_previous*I_q+Edp_previous*I_d+(x_dp-x_qp)*I_d*I_q;

	d_omega= (1/(2*H))*(0.9- Te);
	*omega_present= omega_previous + dt*d_omega;

	d_delta = (omega_previous - omega0);
	*delta_present = delta_previous + dt*d_delta;
}

int main()
{
	struct timeval start, stop;
	//float *I_d, *I_q, *Edp, *Eqp, *omega, *delta, *Te;
	float tint, tfinal, faverage, fsum, et[110], cg_t, fg_t;
	cout<<"Enter the coarse grid step size: "<< endl;
	cin>>cg_t;
	cout<<"Enter the fine grid step size: "<<endl;
	cin>>fg_t;
	cout<< "Enter the initial time step: "<< endl;
	cin>> tint;
	cout<<"Enter the final simulation time: "<<endl;
	cin>> tfinal;
	int coarse_steps = (tfinal-tint)/cg_t;
	int fine_steps = (tfinal-tint)/fg_t;

	// Initialization of variable for coarse funtion
	float *Edp = new float[coarse_steps];
	float *Eqp= new float[coarse_steps];
	float *omega = new float[coarse_steps];
	float *delta = new float[coarse_steps];
	//	float *Te = new float[coarse_steps];


	// Initialization of variables for fine funtion

	float *fine_Edp = new float[coarse_steps];
	float *fine_Eqp = new float[coarse_steps];
	float *fine_omega = new float[coarse_steps];
	float *fine_delta = new float[coarse_steps];
	float *temp_Edp	= new float[fine_steps];
	float *temp_Eqp = new float[fine_steps];
	float *temp_delta = new float[fine_steps];
	float *temp_omega = new float[fine_steps]; 
	float *Edp_fine = new float[coarse_steps];
	float *Eqp_fine = new float[coarse_steps];
	float *omega_fine = new float[coarse_steps];
	float *delta_fine = new float[coarse_steps];

	// Predictor-Corrector initialization
	float *pred_Edp = new float[coarse_steps];
	float *pred_Eqp = new float[coarse_steps];
	float *pred_delta = new float[coarse_steps];
	float *pred_omega = new float[coarse_steps];
	float *correc_Edp = new float[coarse_steps];
	float *correc_Eqp = new float[coarse_steps];
	float *correc_delta = new float[coarse_steps];
	float *correc_omega = new float[coarse_steps];
	float *Edp_diff = new float[coarse_steps];
	float *Eqp_diff = new float[coarse_steps];
	float *delta_diff = new float[coarse_steps];
	float *omega_diff = new float[coarse_steps];

	// Initialization of initial conditions
	float omega0=0;
	omega[0]=omega0;
	delta[0]=(44.3265*pi)/180;
	Edp[0]= -0.4523;
	Eqp[0]= 1.1050;
	int kmax=8;
	float *a = new float[coarse_steps];
	a[0]=tint;
	//gettimeofday(&start,NULL);
	for (int k=0;k<kmax;k++)
	{
		if(k==0)
		{
			for (int i=0;i<coarse_steps;i++)
			{
				a[i+1] = a[i]+cg_t;
				DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],delta[i],omega[i],Edp[i], Eqp[i], cg_t, omega0);
				//				cout<< "The coarse grid values are "<< (delta[i+1]*180)/pi<<endl;	

			}
		}

		else
		{

			for (int i=k;i<coarse_steps;i++)
			{
				Edp[i] = correc_Edp[i];
				Eqp[i] = correc_Eqp[i];
				delta[i] = correc_delta[i];
				omega[i] = correc_omega[i];
				DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],delta[i],omega[i],Edp[i], Eqp[i], cg_t, omega0);
				// cout<< "The coarse grid values for else case are "<< (delta[i+1]*180)/pi<<endl;

			}

		}
		// Fine grid computation
		for(int j=0;j<coarse_steps;j++)
		{
			tint=a[j];
			tfinal=a[j+1];
			int umax=round((tfinal-tint)/fg_t);
			fine_delta[j]=delta[j];
			fine_omega[j]=omega[j];
			fine_Edp[j]=Edp[j];
			fine_Eqp[j]=Eqp[j];
			//cout<<"Initial Value assignment " <<(fine_delta[j]*180)/pi<<" ,"<<fine_omega[j]<<" ,"<<fine_Edp[j]<<" ,"<<fine_Eqp[j]<<" ,"<<umax<<" ,"<<tint<<" ,"<<tfinal<<endl;	
			for (int u=0;u<umax;u++)
			{
				double fine_step = tint+fg_t;
				/*
				//      cout<<"The omega function of fine grid is " <<temp4<<endl;*/
				DifferentialEquations(&temp_delta[umax*j+u],&temp_omega[umax*j+u],&temp_Edp[umax*j+u],&temp_Eqp[umax*j+u],fine_delta[j],fine_omega[j],fine_Edp[j],fine_Eqp[j],fg_t,omega0);
				fine_delta[j]=temp_delta[umax*j+u];
				fine_omega[j]=temp_omega[umax*j+u];
				fine_Edp[j]=temp_Edp[umax*j+u];
				fine_Eqp[j]=temp_Eqp[umax*j+u];
				tint=fine_step;
				//     printf("%lf\n , %d for time fine_step %lf\n",(fine_delta[j]*180)/pi,omp_get_thread_num(),fine_step);
				//   cout<<"The fine grid Delta funtion is " <<(fine_delta[j]*180)/pi<<" for the time step "<< tint<<endl;
				//cout<<(temp_delta[umax*j+u]*180)/pi<<endl;
				//               cout<<"The omega function of fine grid is " <<fine_omega[j]<<endl;

			}
			delta_fine[j+1]=fine_delta[j];
			omega_fine[j+1]=fine_omega[j];
			Edp_fine[j+1]=fine_delta[j];
			Eqp_fine[j+1]=fine_omega[j];
			delta_diff[j]=delta_fine[j+1]-delta[j+1];
			omega_diff[j]=omega_fine[j+1]-omega[j+1];
			Edp_diff[j]=Edp_fine[j+1]-Edp[j+1];
			Eqp_diff[j]=Eqp_fine[j+1]-Eqp[j+1];

			//cout<<"The fine grid Delta funtion is " <<(fine_delta[j]*180)/pi<<" for the time step "<< tfinal<<"for iteration "<<k<<endl;
		}


		pred_delta[0]=delta_fine[1];
		pred_omega[0]=omega_fine[1];
		pred_Edp[0]=Edp_fine[1];
		pred_Eqp[0]=Eqp_fine[1];

		for (int i=0;i<coarse_steps;i++)
		{
			//	Differential(&pred_delt[i+1],&pred_omega[i+1],pred_delt[i],pred_omega[i],omega0,c_h);
			DifferentialEquations(&pred_delta[i+1],&pred_omega[i+1], &pred_Edp[i+1], &pred_Eqp[i+1],pred_delta[i],pred_omega[i],pred_Edp[i], pred_Eqp[i], cg_t, omega0);
		}
		for (int i=0;i<coarse_steps;i++)
		{
			correc_delta[i+1] = delta_diff[i]+pred_delta[i];
			correc_omega[i+1] = omega_diff[i]+pred_omega[i];
			correc_Edp[i+1] = Edp_diff[i]+pred_Edp[i];
			correc_Eqp[i+1] = Eqp_diff[i]+pred_Eqp[i];
			//cout<< "The corrected grid values are "<< (corec_delt[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<endl;
		}

		for(int i=2;i<coarse_steps;i++)
		{
			if(abs(delta[i+1]-correc_delta[i+1])<1e-6)
			{
				if(i<1000)
				{
					cout<<"converged for "<<i+1<<"\telement"<<" for iteration: "<<k<<endl;
					break;
				}
			}
		}
		//cout<<"iteration converged is "<<k<<endl;

	}
	/*	for (int i=0;i<coarse_steps;i++)
		{
		cout<< "Delta values are "<< (delta[i]*180)/pi<<endl;                
		}
	 */


	/*gettimeofday(&stop,NULL);
	  et[0]=((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
	  cout<<"The execution time is "<<et[0]<<" ms"<<endl;
	 */

	delete[] Eqp;
	delete[] Edp;
	delete[] omega;
	delete[] delta;
	delete[] Eqp_fine;
	delete[] Edp_fine;
	delete[] omega_fine;
	delete[] delta_fine;
	delete[] temp_Eqp;
	delete[] temp_Edp;
	delete[] temp_omega;
	delete[] temp_delta;
	delete[] fine_Eqp;
	delete[] fine_Edp;
	delete[] fine_omega;
	delete[] fine_delta;
	delete[] Eqp_diff;
	delete[] Edp_diff;
	delete[] omega_diff;
	delete[] delta_diff;
	delete[] pred_Eqp;
	delete[] pred_Edp;
	delete[] pred_omega;
	delete[] pred_delta;
	delete[] correc_Eqp;
	delete[] correc_Edp;
	delete[] correc_omega;
	delete[] correc_delta;
}







