#include "basic_include.h"

void DifferentialEquations(float *delta_present,float *omega_present,float *Edp_present, float *Eqp_present, float delta_previous,float omega_previous,float Edp_previous, float Eqp_previous, float dt, float omega0, float x_e)
{
	float I_q_predic, I_d_predic, d_delta_predic, d_omega_predic, d_Edp_predic, d_Eqp_predic, Te_predic;
	float I_q_correc, I_d_correc, d_delta_correc, d_omega_correc, d_Edp_correc, d_Eqp_correc;
	float Edp_predic, Eqp_predic, delta_predic, omega_predic,Te;

	/* Predictor Code*/
	I_d_predic=(Eb*cos(delta_previous) - Eqp_previous)/(x_dp+x_e);
	I_q_predic=(Edp_previous+Eb*sin(delta_previous))/(x_qp+x_e);

	d_Eqp_predic = 1/T_dp*(-Eqp_previous+ (x_d-x_dp)*I_d_predic+Efd);
	Eqp_predic = Eqp_previous + dt*d_Eqp_predic;

	d_Edp_predic = 1/T_qp*(-Edp_previous-(x_q-x_qp)*I_q_predic);
	Edp_predic = Edp_previous +dt*d_Edp_predic;

	Te_predic = Eqp_previous*I_q_predic+Edp_previous*I_d_predic+(x_dp-x_qp)*I_d_predic*I_q_predic;

	d_omega_predic= (1/(2*H))*(0.9- Te_predic);
	omega_predic= omega_previous + dt*d_omega_predic;

	d_delta_predic = 2*pi*60*(omega_previous - omega0);
	delta_predic = delta_previous + dt*d_delta_predic;

	/*Corrector Code*/

	I_d_correc=(Eb*cos(delta_predic) - Eqp_predic)/(x_dp+x_e);
	I_q_correc=(Edp_predic+Eb*sin(delta_predic))/(x_qp+x_e);


	d_Eqp_correc = 1/T_dp*(-Eqp_predic+ (x_d-x_dp)*I_d_correc+Efd);
	*Eqp_present = Eqp_previous + (dt/2)*(d_Eqp_predic+d_Eqp_correc);

	d_Edp_correc = 1/T_qp*(-Edp_predic-(x_q-x_qp)*I_q_correc);
	*Edp_present = Edp_previous +(dt/2)*(d_Edp_predic+d_Edp_correc);

	Te = Eqp_predic*I_q_correc+Edp_predic*I_d_correc+(x_dp-x_qp)*I_d_correc*I_q_correc;
	d_omega_correc= (1/(2*H))*(0.9- Te);
	*omega_present= omega_previous + (dt/2)*(d_omega_predic+d_omega_correc);

	d_delta_correc = 2*pi*60*(omega_predic - omega0);
	*delta_present = delta_previous + (dt/2)*(d_delta_correc+d_delta_predic);


}

int main()
{
	struct timeval start, stop;
	//float *I_d, *I_q, *Edp, *Eqp, *omega, *delta, *Te;
	float tint, tfinal, faverage, fsum, et, cg_t, fg_t;
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
	delta[0]=(56.1201*pi)/180;
	Edp[0]= -0.4973;
	Eqp[0]= 1.0337;
	float x_edf, x_epf;
	x_edf = 0.9;
	x_epf = 0.5;
	int kmax=1;
	float *a = new float[coarse_steps];
	a[0]=tint;
	gettimeofday(&start,NULL);
	for (int k=0;k<kmax;k++)
	{
		if(k==0)
		{
			for (int i=0;i<coarse_steps;i++)
			{	
				a[i+1]=a[i]+cg_t;
				if(a[i+1]<=0.5)
				{
					DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],delta[i],omega[i],Edp[i], Eqp[i], cg_t, omega0,x_edf);
					cout<<(delta[i+1]*180)/pi<<endl;
				}
				if(a[i+1]>0.5)
				{
					DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],delta[i],omega[i],Edp[i], Eqp[i], cg_t, omega0,x_epf);
					cout<<(delta[i+1]*180)/pi<<endl;
				}

			}
		}

		else
		{

			cout<<"Iteration is "<<k<<endl;	 
			for (int i=1;i<coarse_steps;i++)
			{
				Edp[i] = correc_Edp[i];
				Eqp[i] = correc_Eqp[i];
				delta[i] = correc_delta[i];
				omega[i] = correc_omega[i];
				if(a[i+1]<=0.5)
				{

					DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],delta[i],omega[i],Edp[i], Eqp[i], cg_t, omega0, x_edf);
					//			cout<<"Iteration is "<<k<<endl;	 
					//		cout<< (delta[i+1]*180)/pi<<endl;
				}

				if(a[i+1]>0.5)
				{
					DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],delta[i],omega[i],Edp[i], Eqp[i], cg_t, omega0,x_epf);
					//                      cout<<(delta[i+1]*180)/pi<<endl;

				}

			}

		}
		//	cout<<"Fine grid computation" << endl;
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
			if (a[j+1]<=0.5)
			{
				for (int u=0;u<umax;u++)
				{
					double fine_step = tint+fg_t;
					/*
					//      cout<<"The omega function of fine grid is " <<temp4<<endl;*/
					DifferentialEquations(&temp_delta[umax*j+u],&temp_omega[umax*j+u],&temp_Edp[umax*j+u],&temp_Eqp[umax*j+u],fine_delta[j],fine_omega[j],fine_Edp[j],fine_Eqp[j],fg_t,omega0, x_edf);
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
			}

			if (a[j+1]>0.5)
			{
				for (int u=0;u<umax;u++)
				{
					double fine_step = tint+fg_t;
					/*
					//      cout<<"The omega function of fine grid is " <<temp4<<endl;*/
					DifferentialEquations(&temp_delta[umax*j+u],&temp_omega[umax*j+u],&temp_Edp[umax*j+u],&temp_Eqp[umax*j+u],fine_delta[j],fine_omega[j],fine_Edp[j],fine_Eqp[j],fg_t,omega0, x_epf);
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
			}



			delta_fine[j+1]=fine_delta[j];
			omega_fine[j+1]=fine_omega[j];
			Edp_fine[j+1]=fine_delta[j];
			Eqp_fine[j+1]=fine_omega[j];
			delta_diff[j]=delta_fine[j+1]-delta[j+1];
			omega_diff[j]=omega_fine[j+1]-omega[j+1];
			Edp_diff[j]=Edp_fine[j+1]-Edp[j+1];
			Eqp_diff[j]=Eqp_fine[j+1]-Eqp[j+1];

			//		cout<<"The fine grid Delta funtion is " <<(delta_fine[j]*180)/pi<<" for the time step "<< tfinal<<"for iteration "<<k<<endl;
		}


		pred_delta[k]=delta[k];
		pred_omega[k]=omega[k];
		pred_Edp[k]=Edp[k];
		pred_Eqp[k]=Eqp[k];

		for (int i=k;i<coarse_steps;i++)
		{
			if(a[i+1]<=0.5)
			{
				DifferentialEquations(&pred_delta[i+1],&pred_omega[i+1], &pred_Edp[i+1], &pred_Eqp[i+1],pred_delta[i],pred_omega[i],pred_Edp[i], pred_Eqp[i], cg_t, omega0,x_edf);

			}
			if(a[i+1]>0.5)
			{
				DifferentialEquations(&pred_delta[i+1],&pred_omega[i+1], &pred_Edp[i+1], &pred_Eqp[i+1],pred_delta[i],pred_omega[i],pred_Edp[i], pred_Eqp[i], cg_t, omega0,x_epf);

			}

		}
		for (int i=0;i<coarse_steps;i++)
		{
			correc_delta[i+1] = delta_diff[i]+delta[i];
			correc_omega[i+1] = omega_diff[i]+omega[i];
			correc_Edp[i+1] = Edp_diff[i]+Edp[i];
			correc_Eqp[i+1] = Eqp_diff[i]+Eqp[i];
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
		/*	for (int i=0;i<coarse_steps;i++)

			{
			cout<<(delta[i+1]*180)/pi<<endl;
			}
		 */

	}
	gettimeofday(&stop,NULL);
	et=((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
	cout<<"The execution time is "<<et<<" ms"<<endl;
	for (int i=0;i<coarse_steps;i++)

	{
		cout<<(delta[i]*180)/pi<<endl;                
	}



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






