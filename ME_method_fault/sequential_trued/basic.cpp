#include "basic_include.h"

void DifferentialEquations(double *delta_present,double *omega_present,double *Edp_present,double *Eqp_present,double *Te,double delta_previous,double omega_previous,double Edp_previous,double Eqp_previous,double dt,double omega0, double x_e)
{
	
	double I_q_predic, I_d_predic, d_delta_predic, d_omega_predic, d_Edp_predic, d_Eqp_predic, Te_predic;
	double I_q_correc, I_d_correc, d_delta_correc, d_omega_correc, d_Edp_correc, d_Eqp_correc;
	double Edp_predic, Eqp_predic, delta_predic, omega_predic;

/* Predictor Code*/
	I_d_predic=(1*cos(delta_previous) - Eqp_previous)/(x_dp+x_e);
        I_q_predic=(Edp_previous+1*sin(delta_previous))/(x_qp+x_e);
	
//	cout<<"Id and Iq are: "<<I_d_predic<<" ,  "<<I_q_predic<<endl;
	
	d_Eqp_predic = (1/T_dp)*(-Eqp_previous+(x_d-x_dp)*I_d_predic+Efd);
        Eqp_predic = Eqp_previous + dt*d_Eqp_predic;
//        cout<< "Eqp is : "<< Eqp_predic<< endl;

        d_Edp_predic = (1/T_qp)*(-Edp_previous-(x_q-x_qp)*I_q_predic);
        Edp_predic = Edp_previous +dt*d_Edp_predic;
//        cout<< "Edp is : "<< Edp_predic<< endl;

        Te_predic = (Eqp_previous*I_q_predic)+(Edp_previous*I_d_predic)+((x_dp-x_qp)*I_d_predic*I_q_predic);
//	cout<< "Te is : "<< Te_predic<< endl;
		
        d_omega_predic= (1/(2*H))*(Tm- Te_predic);
        omega_predic= omega_previous + dt*d_omega_predic;
//        cout<< "The omega values are "<< omega_predic<<endl;

        d_delta_predic = (2*pi*60)*(omega_previous - omega0);
        delta_predic = delta_previous + dt*d_delta_predic;
//       cout<< "The coarse grid values are "<< delta_predic <<endl;

/*Corrector Code*/

	I_d_correc=(1*cos(delta_predic) - Eqp_predic)/(x_dp+x_e);
	I_q_correc=(Edp_predic+1*sin(delta_predic))/(x_qp+x_e);
	
//	cout<<"Id and Iq are: "<<I_d_correc<<" ,  "<<I_q_correc<<endl;
	
	d_Eqp_correc = (1/T_dp)*(-Eqp_predic+(x_d-x_dp)*I_d_correc+Efd);
	*Eqp_present = Eqp_previous + (dt/2)*(d_Eqp_predic+d_Eqp_correc);
//	cout<< "Eqp is : "<< *Eqp_present<< endl;
	
	d_Edp_correc = (1/T_qp)*(-Edp_predic-(x_q-x_qp)*I_q_correc);	
	*Edp_present = Edp_previous +(dt/2)*(d_Edp_predic+d_Edp_correc);
//	cout<< "Edp is : "<< *Edp_present<< endl;

	*Te = (Eqp_predic*I_q_correc)+(Edp_predic*I_d_correc)+((x_dp-x_qp)*I_d_correc*I_q_correc);
//	cout<< "Te is : "<< *Te<< endl;
//	double temp = *Te;
//	cout<<"Temp is: "<< temp<< endl;
	d_omega_correc= (1/(2*H))*(Tm - *Te);
	*omega_present= omega_previous + (dt/2)*(d_omega_predic+d_omega_correc);
//	cout<< "The omega values are "<< *omega_present<<endl;

	d_delta_correc = (2*pi*60)*(omega_predic - omega0);
	*delta_present = delta_previous + (dt/2)*(d_delta_correc+d_delta_predic);
//	cout<< "The coarse grid values are "<< *delta_present <<endl;
}

int main()
{
	struct timeval start, stop;
	//double *I_d, *I_q, *Edp, *Eqp, *omega, *delta, *Te;
	double tint, tfinal, faverage, fsum, et[110], dt;
	cout<<" Enter the step size: "<< endl;
	cin>>dt;
	cout<< "Enter the initial time step: "<< endl;
	cin>> tint;
	cout<<"Enter the final simulation time: "<<endl;
	cin>> tfinal;
	int num_steps = (tfinal-tint)/dt;

	// Initialization of variable
	double *I_d=new double[num_steps];
	double *I_q=new double[num_steps];
	double *Edp = new double[num_steps];
	double *Eqp= new double[num_steps];
	double *omega = new double[num_steps];
	double *delta = new double[num_steps];
	double *Te = new double[num_steps];
	double *a = new double[num_steps];
	double x_edf, x_epf;
	x_edf = 0.9;
	x_epf = 0.5;
	a[0] = tint;
	// Initialization of initial conditions
	double omega0=0;
	omega[0]=omega0;
	delta[0]=(56.1201*pi)/180;
	cout<<"Delta in radians : " << delta[0]<<endl;
	Edp[0]= -0.4973;
	Eqp[0]= 1.0337;
	gettimeofday(&start,NULL);
	for (int i=0;i<num_steps;i++)
	{	
		a[i+1] = a[i]+dt;
		if(a[i+1]<=0.5)
		{
		DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],&Te[i],delta[i],omega[i],Edp[i], Eqp[i], dt, omega0, x_edf);
//		cout<< (delta[i+1]*180)/pi<<endl;	
		}
		if(a[i+1]>0.5)
		{
		DifferentialEquations(&delta[i+1],&omega[i+1], &Edp[i+1], &Eqp[i+1],&Te[i],delta[i],omega[i],Edp[i], Eqp[i], dt, omega0, x_epf);
		}
	}
	gettimeofday(&stop,NULL);
	et[0]=((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
	cout<<"The execution time is "<<et[0]<<" ms"<<endl;
	for (int i=0;i<num_steps;i++)
	{
	 cout<< (delta[i]*180)/pi<<endl;
	}
	


delete[] I_d;
delete[] I_q;
delete[] Eqp;
delete[] Edp;
delete[] omega;
delete[] delta;
delete[] Te;
delete[] a;
}
