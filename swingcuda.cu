#include <cuda_runtime.h>
#include <sys/time.h>
#include "iostream"
#include "iomanip"
#include "cmath"
#include <stdio.h>
using namespace std;
#define pi 3.14159265358979323846
#define CHECK(call) \
{                                                                        \
	const cudaError_t error = call;                                       \
	if (error != cudaSuccess)                                             \
	{                                                                     \
		printf("Error: %s:%d, ", __FILE__, __LINE__);                      \
		printf("code:%d, reason: %s\n", error, cudaGetErrorString(error)); \
		exit(1);                                                           \
	}                                                                     \
}


// When the fault is occured for GPU function
__device__ void Differentiald(double *deltapresent,double *omegapresent,double deltaprevious,double omegaprevious,double omega0,double c_h)
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



//Once the fault is cleared for GPU function

__device__ void Differentialpostfaultd(double *deltapresent,double *omegapresent,double deltaprevious,double omegaprevious,double omega0,double c_h)
{
	double temp,temp1,ddeltapresent,domegapresent,ddeltaprevious,domegaprevious;
	ddeltaprevious =omegaprevious-omega0;
	temp=deltaprevious+(c_h)*(ddeltaprevious);
	domegaprevious =((pi*60)/5)* (0.8-1.4625*sin(temp));
	temp1 = omegaprevious+(c_h*(domegaprevious));
	ddeltapresent = temp1-omegaprevious;
	*deltapresent = deltaprevious + (c_h/2)*(ddeltaprevious+ddeltapresent);
	domegapresent =((pi*60)/5)* (0.8-(1.4625*sin(*deltapresent)));
	*omegapresent = omegaprevious+(c_h/2)*(domegaprevious+domegapresent);
}



__global__ void gpuparareal(double *g_delta,double *g_omega,double *g_a,const double omega0,const double f_h,double *g_fine_tempd,double *g_fine_tempo,double *g_del_fine,double *g_omega_fine,double *g_diff_delta,double *g_diff_omega,int num_steps,int c_h)
{
	const int idx = threadIdx.x + (blockIdx.x*blockDim.x);
	if(idx>=num_steps)
	{
		return;
	}
	double tempd,tempo,tint,tfin,fine_step,fine_tempd,fine_tempo;
	tint = g_a[idx];
	tfin = g_a[idx+1];
	tempd = g_delta[idx];
	tempo = g_omega[idx];
	//printf("g_a[%d] = %lf \n",idx,g_a[idx]);
	//	__syncthreads();
	bool flag = (g_a[idx]<0.8);
	if(flag)
	{
		int umax = round((tfin-tint)/f_h);
		for (int u=0;u<umax;u++)
		{
			fine_step = tint+f_h;
			Differentiald(&fine_tempd,&fine_tempo,tempd,tempo,omega0,f_h);
			tempd=fine_tempd;
			tempo=fine_tempo;
			tint=fine_step;
		}
		//   printf("idx = %d The value of fine is %lf for time %lf\n",idx,(tempd*180/pi),tfin);
	}

	if(!flag)
	{
		int umax = round((tfin-tint)/f_h);
		for (int u=0;u<umax;u++)
		{
			fine_step = tint+f_h;
			Differentialpostfaultd(&fine_tempd,&fine_tempo,tempd,tempo,omega0,f_h);
			tempd=fine_tempd;
			tempo=fine_tempo;
			tint=fine_step;
		}
	}
	g_del_fine[idx+1]=tempd;
	g_omega_fine[idx+1]=tempo;
	g_diff_delta[idx] = tempd - g_delta[idx+1];
	g_diff_omega[idx] = tempo - g_omega[idx+1];
	//	printf("idx = %d The value of fine is %lf for time %lf\n",idx,(tempd*180/pi),tfin);
	//printf("The value of fine is %lf for time %lf\n",(tempd*180/pi),tfin);

}

__global__ void gpucorrection(double *d_diff_delta,double *d_diff_omega,double *d_pred_delta,double *d_pred_omega,double *d_corec_delta,double *d_corec_omega, int num_steps)
{
	const int idx = threadIdx.x + (blockIdx.x*blockDim.x);
	if(idx>=num_steps)
	{
		return;
	}
	d_corec_delta[idx+1]=d_pred_delta[idx]+d_diff_delta[idx];
	d_corec_omega[idx+1]=d_pred_omega[idx]+d_diff_omega[idx];
}



// When the fault is occured
void Differential(double *deltapresent,double *omegapresent,double deltaprevious,double omegaprevious,double omega0,double c_h)
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

void Differentialpostfault(double *deltapresent,double *omegapresent,double deltaprevious,double omegaprevious,double omega0,double c_h)
{
	double temp,temp1,ddeltapresent,domegapresent,ddeltaprevious,domegaprevious;
	ddeltaprevious =omegaprevious-omega0;
	temp=deltaprevious+(c_h)*(ddeltaprevious);
	domegaprevious =((pi*60)/5)* (0.8-1.4625*sin(temp));
	temp1 = omegaprevious+(c_h*(domegaprevious));
	ddeltapresent = temp1-omegaprevious;
	*deltapresent = deltaprevious + (c_h/2)*(ddeltaprevious+ddeltapresent);
	domegapresent =((pi*60)/5)* (0.8-(1.4625*sin(*deltapresent)));
	*omegapresent = omegaprevious+(c_h/2)*(domegaprevious+domegapresent);
}



int main()
{
	cudaEvent_t kernel_start;
	cudaEvent_t kernel_stop;
	struct timeval start,end;
	double tint,tfin,omega0;
	float fElapsedTime,faverage,fsum[10];
	double fMemoryCopyTime[10];
	double fSequential_time[10],tet[10];
	//host variables
	double *h_omega,*h_delta,*h_a,c_h,f_h,dint,*h_del_fine,*h_omega_fine,*h_diff_delta,*h_diff_omega;
	//device variables
	double *d_omega,*d_delta,*d_a,*d_del_fine,*d_omega_fine,*d_fine_tempd,*d_fine_tempo,*d_diff_delta,*d_diff_omega,*d_pred_delta,*d_corec_delta,*d_pred_omega,*d_corec_omega;
	double *h_pred_delta,*h_corec_delta,*h_pred_omega,*h_corec_omega,*h_fine_tempd,*h_fine_tempo,et[110];
	cout<<"The initial time value is : "<<endl;
	cin>>tint;
	cout<<"The final time value is: "<<endl;
	cin>>tfin;
	cout<<"The coarse grid time step value is: "<<endl;
	cin>>c_h;
	cout<<"The fine grid step size value is: "<<endl;
	cin>>f_h;
	cout<<"Enter the intial value of delta in degrees: "<<endl;
	cin>>dint;
	int num_steps = ((tfin-tint)/c_h)+1;
	cout<<"the number of steps for coarse : "<<num_steps<<endl;
	size_t num_steps_bytes_coarse = num_steps*sizeof(double);
	int fine_size = ((tfin-tint)/f_h)+1;
	cout<<"The number of steps for fine : "<<fine_size<<endl;
	size_t num_steps_bytes_fine = fine_size*sizeof(double);
	h_omega = new double[num_steps];
	h_delta = new double[num_steps];
	h_a = new double[num_steps];
	h_del_fine = new double[num_steps];
	h_omega_fine = new double[num_steps];
	h_fine_tempd= new double[fine_size];
	h_fine_tempo=new double[fine_size];
	h_pred_delta= new double[num_steps];
	h_pred_omega= new double[num_steps];
	h_corec_delta=new double[num_steps];
	h_corec_omega = new double[num_steps];
	h_diff_delta = new double[num_steps];
	h_diff_omega = new double [num_steps];
	omega0=2*pi*60;
	h_omega[0]=omega0;
	h_delta[0]=(dint*pi)/180;
	cout<<"The value in radians is: "<<h_delta[0]<<endl;
	h_a[0] =0;
	h_a[0]=tint;
	num_steps =num_steps - 1;
	fine_size =fine_size - 1;
	cout<<num_steps<<endl;	
	for(int k=0;k<2;k++)
	{
		fMemoryCopyTime[k]=0;	
		fSequential_time[k]=0;
		//	gettimeofday(&start,NULL);
		if(k==0)
		{ 	
			gettimeofday(&start,NULL);
			for (int i=0;i<num_steps;i++)
			{
				h_a[i+1]=h_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
				if(h_a[i+1]<=0.8)
				{
					//cout << "a= " <<h_a[i+1]<<__LINE__<<endl;
					//	h_a[i+1]=h_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
					Differential(&h_delta[i+1],&h_omega[i+1],h_delta[i],h_omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (h_delta[i+1]*180)/pi<<" for time"<<h_a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;
					//cout<<"break 2"<<endl;

				}
				if(h_a[i+1]>0.8)
				{
					//cout << "a= " <<h_a[i]<<__LINE__<<endl;
					//h_a[i+1]=h_a[i]+c_h;
					Differentialpostfault(&h_delta[i+1],&h_omega[i+1],h_delta[i],h_omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (h_delta[i+1]*180)/pi<<" for time"<<h_a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl; 
				}

			}
			gettimeofday(&end,NULL);
			fSequential_time[k] = ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
			cout<<" The Sequential Execution time is : "<<fSequential_time[k]<<" ms"<<endl;
		}
		else
		{
			gettimeofday(&start,NULL);
			for(int i=1;i<num_steps;i++)
			{
				h_delta[i]=h_corec_delta[i];
				h_omega[i]=h_corec_omega[i];
				h_a[i+1]=h_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
				if(h_a[i+1]<=0.8)
				{
					//	h_a[i+1]=h_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
					Differential(&h_delta[i+1],&h_omega[i+1],h_delta[i],h_omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (h_delta[i+1]*180)/pi<<" for time"<<h_a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl;

				}
				if(h_a[i+1]>0.8)
				{
					//h_a[i+1]=h_a[i]+c_h;
					Differentialpostfault(&h_delta[i+1],&h_omega[i+1],h_delta[i],h_omega[i],omega0,c_h);
					//cout<< "The coarse grid values are "<< (h_delta[i+1]*180)/pi<<" for time"<<h_a[i+1]<<"for array element "<<i<<"for k value "<<k<<endl; 
				}

			}
			gettimeofday(&end,NULL);
			fSequential_time[k] = ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
			cout<<" The Sequential Execution time is : "<<fSequential_time[k]<<" ms"<<endl;
		}

		//cout<<" The Sequential Execution time is : "<<fSequential_time<<" ms"<<endl;
		//      cudaEvent_t kernel_start;
		//      cudaEvent_t kernel_stop;

		cudaError_t cudaSetDevice(int device);
		cudaSetDevice(0);
		CHECK(cudaEventCreate(&kernel_start));
		CHECK(cudaEventCreate(&kernel_stop));
		//Allocating memory on GPU for device variables

		CHECK(cudaMalloc((double**)&d_delta,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_omega,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_fine_tempd,num_steps_bytes_fine+8));
		CHECK(cudaMalloc((double**)&d_fine_tempo,num_steps_bytes_fine+8));
		CHECK(cudaMalloc((double**)&d_del_fine,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_omega_fine,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_a,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_diff_delta,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_diff_omega,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_pred_delta,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_pred_omega,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_corec_delta,num_steps_bytes_coarse+8));
		CHECK(cudaMalloc((double**)&d_corec_omega,num_steps_bytes_coarse+8));
		//copying the data to device from host
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(d_delta,h_delta,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_omega,h_omega,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_a,h_a,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		gettimeofday(&end,NULL);
		fMemoryCopyTime[k] = ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		//Kernel call
		int ilen = 256;
		dim3 block (ilen,1,1);
		dim3 grid ((num_steps+block.x-1)/block.x,1,1);
		cout << "1D Grid Dimension" << endl;
		cout << "\tNumber of Blocks along X dimension: " << grid.x << endl;
		cout << "1D Block Dimension" << endl;
		cout << "\tNumber of threads along X dimension: " << block.x << endl;
		//kernel function
		et[0]=0;
		for(int i=0;i<110;i++)
		{
			CHECK(cudaEventRecord(kernel_start));
			gpuparareal<<<grid,block>>>(d_delta,d_omega,d_a,omega0,f_h,d_fine_tempd,d_fine_tempo,d_del_fine,d_omega_fine,d_diff_delta,d_diff_omega,num_steps,c_h);
			CHECK(cudaEventRecord(kernel_stop));
			CHECK(cudaEventSynchronize(kernel_stop));
			CHECK(cudaEventElapsedTime(&fElapsedTime,kernel_start,kernel_stop));
			et[i]=fElapsedTime;	
		}
		//	cout << "Kernel with Compiler Implementation = " << fElapsedTime << " msecs" << endl;
		//	gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(h_fine_tempd,d_fine_tempd,num_steps_bytes_fine+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_fine_tempo,d_fine_tempo,num_steps_bytes_fine+8,cudaMemcpyDeviceToHost));
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(h_diff_delta,d_diff_delta,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_diff_omega,d_diff_omega,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_del_fine,d_del_fine,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_omega_fine,d_omega_fine,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		gettimeofday(&end,NULL);
		fMemoryCopyTime[k]+= ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		cout<< "Memory transfer time =  " << fMemoryCopyTime[k] <<" ms"<<endl;
		h_pred_delta[0]=h_del_fine[1];
		h_pred_omega[0]=h_omega_fine[1];
		cout<<"Fine values are: "<<"\tdelta"<< (h_del_fine[1]*180/pi)<<"\t omega"<<h_omega_fine[1]<<endl;
		gettimeofday(&start,NULL);
		for (int i=1;i<num_steps;i++)
		{
			h_a[i+1]=h_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
			if(h_a[i+1]<=0.8)
			{
				//h_a[i+1]=h_a[i]+c_h;
				Differential(&h_pred_delta[i],&h_pred_omega[i],h_pred_delta[i-1],h_pred_omega[i-1],omega0,c_h);
				//cout<<"The predicted value is "<<(h_pred_delta[i]*180)/pi<<endl;
			}
			if(h_a[i+1]>0.8)
			{
				//h_a[i+1]=h_a[i]+c_h;
				Differentialpostfault(&h_pred_delta[i],&h_pred_omega[i],h_pred_delta[i-1],h_pred_omega[i-1],omega0,c_h);
				//cout<<"The predicted value is "<<(h_pred_delta[i]*180)/pi<<endl;
			}

		}
		gettimeofday(&end,NULL);
		fSequential_time[k] += ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(d_diff_delta,h_diff_delta,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_diff_omega,h_diff_omega,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_pred_delta,h_pred_delta,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_pred_omega,h_pred_omega,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		gettimeofday(&end,NULL);
		double time = 0;
		time = ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		//	int ilen = 256;
		//      dim3 block (ilen,1,1);
		//    dim3 grid ((num_steps+block.x-1)/block.x,1,1);
		CHECK(cudaEventRecord(kernel_start));
		gpucorrection<<<grid,block>>>(d_diff_delta,d_diff_omega,d_pred_delta,d_pred_omega,d_corec_delta,d_corec_omega,num_steps);
		CHECK(cudaEventRecord(kernel_stop));
		CHECK(cudaEventSynchronize(kernel_stop));
		CHECK(cudaEventElapsedTime(&fElapsedTime,kernel_start,kernel_stop));
		cout<<"Elapsed time is for correction is " <<fElapsedTime<<endl;
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(h_corec_delta,d_corec_delta,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_corec_omega,d_corec_omega,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		gettimeofday(&end,NULL);
		time+= ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		cout<<"Correction memory copy time is: "<<time<<" ms"<<endl;
		/*for (int i=0;i<num_steps;i++)
		  {
		  corec_delt[i+1] = h_diff_delta[i]+pred_delt[i];
		  corec_omega[i+1] = h_diff_omega[i]+pred_omega[i];
		//cout<< "The corrected grid values are "<< (corec_delt[i+1]*180)/pi<<" for time"<<h_a[i+1]<<"for array element "<<i<<endl;
		}*/

		fSequential_time[k] += ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		faverage=0;
		for(int i=10;i<110;i++)
		{
			faverage+=et[i];
		}
		fsum[k]=faverage/100;
		cout<<"The gpu execution time is "<<fsum[k]<<"\t"<<"sequential time is "<<fSequential_time[k]<<" ms"<<endl;
		tet[k]=fsum[k]+fSequential_time[k]+fMemoryCopyTime[k];
		cout<<"the elapsed time is "<<tet[k]<<" ms"<<endl;

		CHECK(cudaEventDestroy(kernel_start));
		CHECK(cudaEventDestroy(kernel_stop));
		CHECK(cudaFree(d_omega));
		CHECK(cudaFree(d_delta));
		CHECK(cudaFree(d_fine_tempd));
		CHECK(cudaFree(d_fine_tempo));
		CHECK(cudaFree(d_del_fine));
		CHECK(cudaFree(d_omega_fine));
		CHECK(cudaFree(d_a));
		CHECK(cudaFree(d_diff_delta));
		CHECK(cudaFree(d_diff_omega));
		CHECK(cudaFree(d_pred_delta));
		CHECK(cudaFree(d_pred_omega));
		CHECK(cudaFree(d_corec_delta));
		CHECK(cudaFree(d_corec_omega));
		CHECK(cudaDeviceReset());
	}
	/*	delete[] h_omega;
		delete[] h_delta;
		delete[] h_a;
		delete[] h_del_fine;
		delete[] h_omega_fine;
		delete[] h_fine_tempd;
		delete[] h_fine_tempo;
		delete[] h_pred_delta;
		delete[] h_pred_omega;
		delete[] h_corec_delta;
		delete[] h_corec_omega;
		delete[] h_diff_delta;
		delete[] h_diff_omega;*/



}







