#include "basic_include.h"

void DifferentialEquations(float *delta_present,float *omega_present,float *Edp_present, float *Eqp_present, float delta_previous,float omega_previous,float Edp_previous, float Eqp_previous, float dt, float omega0, const float x_e)
{
	float I_q_predic, I_d_predic, d_delta_predic, d_omega_predic, d_Edp_predic, d_Eqp_predic, Te_predic;
	float I_q_correc, I_d_correc, d_delta_correc, d_omega_correc, d_Edp_correc, d_Eqp_correc;
	float Edp_predic, Eqp_predic, delta_predic, omega_predic, Te;

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


__device__ void DifferentialEquationsGPU(float *delta_present,float *omega_present,float *Edp_present, float *Eqp_present, float delta_previous,float omega_previous,float Edp_previous, float Eqp_previous, float dt, float omega0, const float x_e)
{
	float I_q_predic, I_d_predic, d_delta_predic, d_omega_predic, d_Edp_predic, d_Eqp_predic, Te_predic;
	float I_q_correc, I_d_correc, d_delta_correc, d_omega_correc, d_Edp_correc, d_Eqp_correc;
	float Edp_predic, Eqp_predic, delta_predic, omega_predic, Te;

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

	d_delta_correc = 2*60*pi*(omega_predic - omega0);
	*delta_present = delta_previous + (dt/2)*(d_delta_correc+d_delta_predic);
}

__global__ void gpuparareal(float *g_delta,float *g_omega,float *g_Edp, float *g_Eqp, float *g_a,const float omega0,const float f_h,float *g_temp_delta,float *g_temp_omega,float *g_temp_Edp, float *g_temp_Eqp, float *g_delta_fine,float *g_omega_fine,float *g_Edp_fine, float *g_Eqp_fine,float *g_delta_diff,float *g_omega_diff,float *g_Edp_diff, float *g_Eqp_diff,int num_steps,const float c_h, const float x_edf, const float x_epf)
{
	const int idx = threadIdx.x + (blockIdx.x*blockDim.x);
	if(idx>=num_steps)
	{
		return;
	}
	float tempd,tempo,tempEdp,tempEqp,tint,tfin,fine_step,fine_tempd,fine_tempo, fine_tempEdp, fine_tempEqp;
	//      printf("g_a[%d] = %lf \n",idx,g_a[idx]);
	//      __syncthreads();
	tint = g_a[idx];
	tfin = g_a[idx+1];
	tempd = g_delta[idx];
	tempo = g_omega[idx];
	tempEdp = g_Edp[idx];
	tempEqp = g_Eqp[idx];
	bool flag = (g_a[idx]<=0.5);
	if(flag)
	{

		// printf("idx=%d  g_a = %lf g_a fin =%lf tint = %lf  tfinal = %lf \n\n",idx,g_a[idx],g_a[idx+1],tint,tfin);
		int umax = round((tfin-tint)/f_h);
		for (int u=0;u<umax;u++)
		{
			fine_step = tint+f_h;
			DifferentialEquationsGPU(&fine_tempd,&fine_tempo,&fine_tempEdp,&fine_tempEqp,tempd,tempo,tempEdp,tempEqp,f_h,omega0, x_edf);
			tempd=fine_tempd;
			tempo=fine_tempo;
			tempEdp=fine_tempEdp;
			tempEqp=fine_tempEqp;
			tint=fine_step;
			//printf("%f for time %f\n",(tempd*180)/pi,fine_step);
		}
		// printf("idx = %d The value of fine is %lf for time %lf\n",idx,(tempd*180/pi),tfin);
	}
	if(!flag)
	{

		// printf("idx=%d  g_a = %lf g_a fin =%lf tint = %lf  tfinal = %lf \n\n",idx,g_a[idx],g_a[idx+1],tint,tfin);
		int umax = round((tfin-tint)/f_h);
		for (int u=0;u<umax;u++)
		{
			fine_step = tint+f_h;
			DifferentialEquationsGPU(&fine_tempd,&fine_tempo,&fine_tempEdp,&fine_tempEqp,tempd,tempo,tempEdp,tempEqp,f_h,omega0, x_epf);
			tempd=fine_tempd;
			tempo=fine_tempo;
			tempEdp=fine_tempEdp;
			tempEqp=fine_tempEqp;
			tint=fine_step;
			//printf("%f for time %f\n",(tempd*180)/pi,fine_step);
		}
		// printf("idx = %d The value of fine is %lf for time %lf\n",idx,(tempd*180/pi),tfin);
	}

	g_delta_fine[idx+1]=tempd;
	g_omega_fine[idx+1]=tempo;
	g_Edp_fine[idx+1]=tempEdp;
	g_Eqp_fine[idx+1]=tempEqp;
	g_delta_diff[idx] = tempd - g_delta[idx+1];
	g_omega_diff[idx] = tempo - g_omega[idx+1];
	g_Edp_diff[idx] = tempEdp - g_Edp[idx+1];
	g_Eqp_diff[idx] = tempEqp - g_Eqp[idx+1];
	//      printf("idx = %d The value of fine is %lf for time %lf\n",idx,(tempd*180/pi),tfin);
	//      printf("The value of fine is %lf for time %lf\n",(tempd*180/pi),tfin);
}







__global__ void gpucorrection(float *d_delta_diff,float *d_omega_diff,float *d_Edp_diff, float *d_Eqp_diff,float *d_pred_delta,float *d_pred_omega,float *d_pred_Edp, float *d_pred_Eqp, float *d_correc_delta,float *d_correc_omega, float *d_correc_Edp, float *d_correc_Eqp, int num_steps)
{
	const int idx = threadIdx.x + (blockIdx.x*blockDim.x);
	if(idx>=num_steps)
	{
		return;
	}
	d_correc_delta[idx+1]=d_pred_delta[idx]+d_delta_diff[idx];
	d_correc_omega[idx+1]=d_pred_omega[idx]+d_omega_diff[idx];
	d_correc_Edp[idx+1] = d_Edp_diff[idx]+d_pred_Edp[idx];
	d_correc_Eqp[idx+1] = d_Eqp_diff[idx]+d_pred_Eqp[idx];
}




int main()

{
	cudaEvent_t kernel_start;
	cudaEvent_t kernel_stop;
	struct timeval start,stop;
	float tint, tfinal, faverage, fsum[10], et[110], cg_t, fg_t;
	float fElapsedTime;
	float fMemoryCopyTime[10];
	float fSequential_time[10];
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
	size_t num_steps_bytes_coarse = coarse_steps*sizeof(float);
	size_t num_steps_bytes_fine = fine_steps*sizeof(float);

	// Initialization of variable for coarse funtion
	float *h_Edp = new float[coarse_steps];
	float *h_Eqp= new float[coarse_steps];
	float *h_omega = new float[coarse_steps];
	float *h_delta = new float[coarse_steps];
	float *d_Edp,*d_Eqp,*d_omega,*d_delta; 

	// Initialization of variables for fine funtion

	float *h_fine_Edp = new float[coarse_steps];
	float *h_fine_Eqp = new float[coarse_steps];
	float *h_fine_omega = new float[coarse_steps];
	float *h_fine_delta = new float[coarse_steps];
	float *h_temp_Edp = new float[fine_steps];
	float *h_temp_Eqp = new float[fine_steps];
	float *h_temp_delta = new float[fine_steps];
	float *h_temp_omega = new float[fine_steps];
	float *h_Edp_fine = new float[coarse_steps];
	float *h_Eqp_fine = new float[coarse_steps];
	float *h_omega_fine = new float[coarse_steps];
	float *h_delta_fine = new float[coarse_steps];

	float *d_fine_Edp,*d_fine_Eqp,*d_fine_omega,*d_fine_delta;
	float *d_temp_Edp,*d_temp_Eqp,*d_temp_omega,*d_temp_delta;
	float *d_Edp_fine,*d_Eqp_fine,*d_omega_fine,*d_delta_fine;



	// Predictor-Corrector initialization
	float *h_pred_Edp = new float[coarse_steps];
	float *h_pred_Eqp = new float[coarse_steps];
	float *h_pred_delta = new float[coarse_steps];
	float *h_pred_omega = new float[coarse_steps];
	float *h_correc_Edp = new float[coarse_steps];
	float *h_correc_Eqp = new float[coarse_steps];
	float *h_correc_delta = new float[coarse_steps];
	float *h_correc_omega = new float[coarse_steps];
	float *h_Edp_diff = new float[coarse_steps];
	float *h_Eqp_diff = new float[coarse_steps];
	float *h_delta_diff = new float[coarse_steps];
	float *h_omega_diff = new float[coarse_steps];


	float *d_pred_Edp,*d_pred_Eqp,*d_pred_omega,*d_pred_delta;
	float *d_correc_Edp,*d_correc_Eqp,*d_correc_omega,*d_correc_delta;
	float *d_Edp_diff,*d_Eqp_diff,*d_omega_diff,*d_delta_diff;
	float *d_a;
	float x_edf, x_epf;
	x_edf = 0.9;
	x_epf = 0.5;
	float omega0=0;
	h_omega[0]=omega0;
	h_delta[0]=(56.1210*pi)/180;
	h_Edp[0]= -0.4973;
	h_Eqp[0]= 1.0337;
	int kmax=2;
	float *h_a = new float[coarse_steps];
	h_a[0]=tint;
	
	for (int i=0;i<coarse_steps;i++)
	{
	h_a[i+1] = h_a[i]+cg_t;
	}
		
	for (int k=0;k<kmax;k++)
	{	
		// Coarse Solver
		gettimeofday(&start,NULL);
		if(k==0)
		{	
			gettimeofday(&start,NULL);
			for (int i=0;i<coarse_steps;i++)
			{
				//h_a[i+1] = h_a[i]+cg_t;
				if(h_a[i+1]<=0.5)
				{
				DifferentialEquations(&h_delta[i+1],&h_omega[i+1], &h_Edp[i+1], &h_Eqp[i+1],h_delta[i],h_omega[i],h_Edp[i], h_Eqp[i], cg_t, omega0, x_edf);
				//    cout<< "The coarse grid values are "<< (h_delta[i+1]*180)/pi<<endl;
				}
				if(h_a[i+1]>0.5)
                                {
                                DifferentialEquations(&h_delta[i+1],&h_omega[i+1], &h_Edp[i+1], &h_Eqp[i+1],h_delta[i],h_omega[i],h_Edp[i], h_Eqp[i], cg_t, omega0,x_epf);
                                //    cout<< "The coarse grid values are "<< (h_delta[i+1]*180)/pi<<endl;
                                }

			}
				gettimeofday(&stop,NULL);
			//	fSequential_time[k] = ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		}

		else
		{
			gettimeofday(&start,NULL);
			for (int i=k;i<coarse_steps;i++)
			{
				h_Edp[i] = h_correc_Edp[i];
				h_Eqp[i] = h_correc_Eqp[i];
				h_delta[i] = h_correc_delta[i];
				h_omega[i] = h_correc_omega[i];
				//h_a[i+1]=h_a[i]+cg_t;
				if(h_a[i+1]<=0.5)
                                {
				DifferentialEquations(&h_delta[i+1],&h_omega[i+1], &h_Edp[i+1], &h_Eqp[i+1],h_delta[i],h_omega[i],h_Edp[i], h_Eqp[i], cg_t, omega0, x_edf);
				// cout<< "The coarse grid values for else case are "<< (delta[i+1]*180)/pi<<endl;
				}
				if(h_a[i+1]>0.5)
                                {
                                DifferentialEquations(&h_delta[i+1],&h_omega[i+1], &h_Edp[i+1], &h_Eqp[i+1],h_delta[i],h_omega[i],h_Edp[i], h_Eqp[i], cg_t, omega0,x_epf);
                                //    cout<< "The coarse grid values are "<< (h_delta[i+1]*180)/pi<<endl;
                                }

			}
			 gettimeofday(&stop,NULL);

		}
		gettimeofday(&stop,NULL);
		fSequential_time[k] = ((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		// Fine Solver
		CHECK(cudaEventCreate(&kernel_start));
		CHECK(cudaEventCreate(&kernel_stop));
		//Allocating memory on GPU for device variables

		CHECK(cudaMalloc((float**)&d_delta,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_omega,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_Eqp,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_Edp,num_steps_bytes_coarse+4));

		CHECK(cudaMalloc((float**)&d_temp_delta,num_steps_bytes_fine+4));
		CHECK(cudaMalloc((float**)&d_temp_omega,num_steps_bytes_fine+4));
		CHECK(cudaMalloc((float**)&d_temp_Edp,num_steps_bytes_fine+4));
		CHECK(cudaMalloc((float**)&d_fine_Eqp,num_steps_bytes_fine+4));

		CHECK(cudaMalloc((float**)&d_fine_delta,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_fine_omega,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_fine_Eqp,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_fine_Edp,num_steps_bytes_coarse+4));		

		CHECK(cudaMalloc((float**)&d_delta_fine,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_omega_fine,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_Eqp_fine,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_Edp_fine,num_steps_bytes_coarse+4));

		CHECK(cudaMalloc((float**)&d_a,num_steps_bytes_coarse+4));

		CHECK(cudaMalloc((float**)&d_delta_diff,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_omega_diff,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_Edp_diff,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_Eqp_diff,num_steps_bytes_coarse+4));

		CHECK(cudaMalloc((float**)&d_pred_delta,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_pred_omega,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_pred_Edp,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_pred_Eqp,num_steps_bytes_coarse+4));

		CHECK(cudaMalloc((float**)&d_correc_delta,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_correc_omega,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_correc_Edp,num_steps_bytes_coarse+4));
		CHECK(cudaMalloc((float**)&d_correc_Eqp,num_steps_bytes_coarse+4));
		//copying the data to device from host
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(d_delta,h_delta,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_omega,h_omega,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_Edp,h_Edp,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_Eqp,h_Eqp,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_a,h_a,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		gettimeofday(&stop,NULL);
		fMemoryCopyTime[k] = ((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;


		int ilen = 256;
		dim3 block (ilen,1,1);
		dim3 grid ((coarse_steps+block.x-1)/block.x,1,1);
		/*		cout << "1D Grid Dimension" << endl;
				cout << "\tNumber of Blocks along X dimension: " << grid.x << endl;
				cout << "1D Block Dimension" << endl;
				cout << "\tNumber of threads along X dimension: " << block.x << endl;
		 */
		et[0]=0;
		for(int i=0;i<110;i++)
		{	
			et[i]=0;
			CHECK(cudaEventRecord(kernel_start));
			gpuparareal<<<grid,block>>>(d_delta,d_omega,d_Edp,d_Eqp,d_a,omega0,fg_t,d_temp_delta,d_temp_omega,d_temp_Edp,d_temp_Eqp,d_delta_fine,d_omega_fine,d_Edp_fine, d_Eqp_fine, d_delta_diff,d_omega_diff,d_Edp_diff,d_Eqp_diff, coarse_steps,cg_t, x_edf,x_epf);
			CHECK(cudaEventRecord(kernel_stop));
			CHECK(cudaEventSynchronize(kernel_stop));
			CHECK(cudaEventElapsedTime(&fElapsedTime,kernel_start,kernel_stop));
			et[i]=fElapsedTime;
			//		cout<<"Et is "<< et<<endl;
		}


		/*CHECK(cudaMemcpy(h_temp_delta,d_temp_delta,num_steps_bytes_fine+4,cudaMemcpyDeviceToHost));
		  CHECK(cudaMemcpy(h_temp_omega,d_temp_omega,num_steps_bytes_fine+4,cudaMemcpyDeviceToHost));
		  CHECK(cudaMemcpy(h_temp_Edp,d_temp_Edp,num_steps_bytes_fine+4,cudaMemcpyDeviceToHost));
		  CHECK(cudaMemcpy(h_temp_Eqp,d_temp_Eqp,num_steps_bytes_fine+4,cudaMemcpyDeviceToHost));
		 */
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(h_delta_diff,d_delta_diff,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_omega_diff,d_omega_diff,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_Eqp_diff,d_Eqp_diff,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_Edp_diff,d_Edp_diff,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));

		CHECK(cudaMemcpy(h_delta_fine,d_delta_fine,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_omega_fine,d_omega_fine,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_Eqp_fine,d_Eqp_fine,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_Edp_fine,d_Edp_fine,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		gettimeofday(&stop,NULL);

		fMemoryCopyTime[k]+= ((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		//                cout<< "Memory transfer time =  " << fMemoryCopyTime[k] <<" ms"<<endl;		
		h_pred_delta[k]=h_delta[k];
		h_pred_omega[k]=h_omega[k];
		h_pred_Edp[k]=h_Edp[k];
		h_pred_Eqp[k]=h_Eqp[k];
		gettimeofday(&start,NULL);
		for (int i=k;i<coarse_steps;i++)
		{
			if(h_a[i+1]<=0.5)
                        {
                         
			//      Differential(&pred_delt[i+1],&pred_omega[i+1],pred_delt[i],pred_omega[i],omega0,c_h);
			DifferentialEquations(&h_pred_delta[i+1],&h_pred_omega[i+1], &h_pred_Edp[i+1], &h_pred_Eqp[i+1],h_pred_delta[i],h_pred_omega[i],h_pred_Edp[i], h_pred_Eqp[i], cg_t, omega0, x_edf);
			}
			
			if(h_a[i+1]>0.5)
                        {

                        //      Differential(&pred_delt[i+1],&pred_omega[i+1],pred_delt[i],pred_omega[i],omega0,c_h);
                        DifferentialEquations(&h_pred_delta[i+1],&h_pred_omega[i+1], &h_pred_Edp[i+1], &h_pred_Eqp[i+1],h_pred_delta[i],h_pred_omega[i],h_pred_Edp[i], h_pred_Eqp[i], cg_t, omega0, x_epf);
                        }
		}
		gettimeofday(&stop,NULL);
		fSequential_time[k] += ((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;

		/*	for (int i=0;i<coarse_steps;i++)
			{
			correc_delta[i+1] = delta_diff[i]+pred_delta[i];
			correc_omega[i+1] = omega_diff[i]+pred_omega[i];
			correc_Edp[i+1] = Edp_diff[i]+pred_Edp[i];
			correc_Eqp[i+1] = Eqp_diff[i]+pred_Eqp[i];
		//cout<< "The corrected grid values are "<< (corec_delt[i+1]*180)/pi<<" for time"<<a[i+1]<<"for array element "<<i<<endl;
		}
		 */

		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(d_delta_diff,h_delta_diff,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_omega_diff,h_omega_diff,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_Edp_diff,h_Edp_diff,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_Eqp_diff,h_Eqp_diff,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_pred_delta,h_pred_delta,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_pred_omega,h_pred_omega,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_pred_Edp,h_pred_Edp,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_pred_Eqp,h_pred_Eqp,num_steps_bytes_coarse+4,cudaMemcpyHostToDevice));
		gettimeofday(&stop,NULL);
		float time = 0;
		time = ((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		//      int ilen = 256;
		//      dim3 block (ilen,1,1);
		//    dim3 grid ((num_steps+block.x-1)/block.x,1,1);
		CHECK(cudaEventRecord(kernel_start));
		gpucorrection<<<grid,block>>>(d_delta_diff,d_omega_diff,d_Edp_diff, d_Eqp_diff, d_pred_delta,d_pred_omega,d_Edp,d_Eqp, d_correc_delta,d_correc_omega,d_correc_Edp,d_correc_Eqp,coarse_steps);
		CHECK(cudaEventRecord(kernel_stop));
		CHECK(cudaEventSynchronize(kernel_stop));
		CHECK(cudaEventElapsedTime(&fElapsedTime,kernel_start,kernel_stop));
		//	cout<<"Elapsed time is for correction is " <<fElapsedTime<<" ms"<<endl;
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(h_correc_delta,d_correc_delta,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_correc_omega,d_correc_omega,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_correc_Edp,d_correc_Edp,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_correc_Eqp,d_correc_Eqp,num_steps_bytes_coarse+4,cudaMemcpyDeviceToHost));
		gettimeofday(&stop,NULL);
		time+= ((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		// cout<<"Correction memory copy time is: "<<time<<" ms"<<endl;

		for(int i=k;i<coarse_steps;i++)
		{
			if(abs(h_delta[i+1]-h_correc_delta[i+1])<1e-6)
			{
				if(i<1000)
				{
				        cout<<"converged for "<<i+1<<"\telement"<<" for iteration: "<<k<<endl;
					break;
				}
			}
		}

		//		fSequential_time[k] += ((stop.tv_sec*1e6+stop.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		faverage = 0;	
		for(int i=10;i<110;i++)
		{
			faverage+=et[i];
		}
		fsum[k]=faverage/100;
		cout<<"The gpu execution time  for fine grid is: "<<fsum[k]<<" ms"<<"\t"<<"sequential time for coarse + predictor is: "<<fSequential_time[k]<<" ms"<<endl;
		cout<<"Memory transfer time for fine grid:  " << fMemoryCopyTime[k] <<" ms"<<"\t"<<"Memory transfer for Corrector is: "<<time<<" ms"<<endl;
		cout<<"Elapsed time is for correction on GPU is: " <<fElapsedTime<<" ms"<<endl;
		//cout<<"Correction memory copy time is: "<<time<<" ms"<<endl;
	}
/*		for (int i=0;i<coarse_steps;i++)
		{
		cout<< "Delta values are "<< (h_delta[i]*180)/pi<<endl;
		}
*/

	CHECK(cudaEventDestroy(kernel_start));
	CHECK(cudaEventDestroy(kernel_stop));
	CHECK(cudaFree(d_omega));
	CHECK(cudaFree(d_delta));
	CHECK(cudaFree(d_Edp));
	CHECK(cudaFree(d_Eqp));
	CHECK(cudaFree(d_temp_delta));
	CHECK(cudaFree(d_temp_omega));
	//CHECK(cudaFree(d_temp_Edp));
	//CHECK(cudaFree(d_temp_Eqp));
	CHECK(cudaFree(d_delta_fine));
	CHECK(cudaFree(d_omega_fine));
	CHECK(cudaFree(d_Edp_fine));
	CHECK(cudaFree(d_Eqp_fine));
	CHECK(cudaFree(d_a));
	CHECK(cudaFree(d_delta_diff));
	CHECK(cudaFree(d_omega_diff));
	CHECK(cudaFree(d_Edp_diff));
	CHECK(cudaFree(d_Eqp_diff));
	CHECK(cudaFree(d_fine_delta));
	CHECK(cudaFree(d_fine_omega));
	CHECK(cudaFree(d_fine_Edp));
	CHECK(cudaFree(d_fine_Eqp));
	CHECK(cudaFree(d_pred_delta));
	CHECK(cudaFree(d_pred_omega));
	CHECK(cudaFree(d_pred_Edp));
	CHECK(cudaFree(d_pred_Eqp));
	CHECK(cudaFree(d_correc_delta));
	CHECK(cudaFree(d_correc_omega));
	CHECK(cudaFree(d_correc_Edp));
	CHECK(cudaFree(d_correc_Eqp));
	CHECK(cudaDeviceReset());

	delete[] h_Eqp;
	delete[] h_Edp;
	delete[] h_omega;
	delete[] h_delta;
	delete[] h_Eqp_fine;
	delete[] h_Edp_fine;
	delete[] h_omega_fine;
	delete[] h_delta_fine;
	delete[] h_temp_Eqp;
	delete[] h_temp_Edp;
	delete[] h_temp_omega;
	delete[] h_temp_delta;
	delete[] h_fine_Eqp;
	delete[] h_fine_Edp;
	delete[] h_fine_omega;
	delete[] h_fine_delta;
	delete[] h_Eqp_diff;
	delete[] h_Edp_diff;
	delete[] h_omega_diff;
	delete[] h_delta_diff;
	delete[] h_pred_Eqp;
	delete[] h_pred_Edp;
	delete[] h_pred_omega;
	delete[] h_pred_delta;
	delete[] h_correc_Eqp;
	delete[] h_correc_Edp;
	delete[] h_correc_omega;
	delete[] h_correc_delta;
	delete[] h_a;
}
