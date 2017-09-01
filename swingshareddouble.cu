#include <cuda_runtime.h>
#include <sys/time.h>
#include "iostream"
#include "iomanip"
#include "cmath"
#include <stdio.h>
using namespace std;
#define pi 3.14159265358979323846
#define tile 256
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



__global__ void gpuparareal(double *g_delta,double *g_omega,double *g_a,const double omega0,const double f_h,double *g_fine_tempd,double *g_fine_tempo,double *g_del_fine,double *g_omega_fine,double *g_diff_delta,double *g_diff_omega,int num_steps, const double c_h,double *g_correct_delta)
{
	const int idx = threadIdx.x + (blockIdx.x*blockDim.x);
	const int tid = threadIdx.x;
	if(idx>=num_steps)
	{
		return;
	}
	double tempd,tempo,tint,tfin,fine_step;
	bool flag = (g_a[idx]<0.8);
	__shared__ double s_fine_tempo[tile];
	__shared__ double s_fine_tempd[tile];
	__shared__ double s_predict_delta[tile+2];
	__shared__ double s_predict_omega[tile+2];
	__shared__ double s_a[tile];
	__shared__ double s_correct_delta[tile+2];
	__shared__ double s_correct_omega[tile+2];
	s_a[tid]=g_a[idx];
	tint = g_a[idx];
	int condition = num_steps/tile;
	tfin = g_a[idx+1];
	tempd = g_delta[idx];
	tempo = g_omega[idx];
	for (int tw=0;tw<1;tw++)
	{
//	int tw =1;
		if(flag)
		{
			int umax = round((tfin-tint)/f_h);
			for (int u=0;u<umax;u++)
			{
				fine_step = tint+f_h;
				Differentiald(&s_fine_tempd[tid],&s_fine_tempo[tid],tempd,tempo,omega0,f_h);
				tempd=s_fine_tempd[tid];
				tempo=s_fine_tempo[tid];
				tint=fine_step;
			}
		}

		if(!flag)
		{
			int umax = round((tfin-tint)/f_h);
			for (int u=0;u<umax;u++)
			{
				fine_step = tint+f_h;
				Differentialpostfaultd(&s_fine_tempd[tid],&s_fine_tempo[tid],tempd,tempo,omega0,f_h);
				tempd=s_fine_tempd[tid];
				tempo=s_fine_tempo[tid];
				tint=fine_step;
			}
		}
		g_del_fine[idx+1]=s_fine_tempd[tid];
		g_omega_fine[idx+1]=s_fine_tempo[tid];
		g_diff_delta[idx]=s_fine_tempd[tid] - g_delta[idx+1];
		g_diff_omega[idx]=s_fine_tempo[tid] - g_omega[idx+1];
		//printf("idx = %d The value of fine is %f for time %f for local threadID %d for loop iteration %d condition %d\n",idx,(tempd*180/pi),tfin,tid,tw,condition);
		if(tid == 0 )
		{
			s_predict_delta[0] = s_fine_tempd[0];
			s_predict_omega[0] = s_fine_tempo[0];
		//	printf("The initial delta and omega are %f and %f\n", s_predict_delta[0],s_predict_omega[0]);
			s_a[1]=s_a[0]+c_h;
		//	printf("the value of time step is %f\n",s_a[1]);
			for (int i=1;i<tile;i++)
			{
				s_a[i+1]=s_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
				if(s_a[i+1]<=0.8)
				{
					Differentiald(&s_predict_delta[i],&s_predict_omega[i],s_predict_delta[i-1],s_predict_omega[i-1],omega0,c_h);
		//			printf("The predicted value is %f for s_a[i+1] is %f\n",(s_predict_delta[i]*180)/pi,s_a[i+1]);
				}
				if(s_a[i+1]>0.8)
				{
					Differentialpostfaultd(&s_predict_delta[i],&s_predict_omega[i],s_predict_delta[i-1],s_predict_omega[i-1],omega0,c_h);
		//			printf("The predicted value is %f for s_a[i+1] is %f\n",(s_predict_delta[i]*180)/pi,s_a[i+1]);
				}

			}
		}
		__syncthreads();
		s_correct_delta[tid] = s_predict_delta[tid]+g_diff_delta[tid];
		g_correct_delta[idx+1] = s_correct_delta[tid];
		s_correct_omega[tid] = s_predict_omega[tid]+g_diff_omega[tid];		
		//printf("The corrected value is %f for local thread ID %d\n",(s_correct_delta[tid]*180)/pi,tid);
		
	}
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
	float fMemoryCopyTime[10];
	float fSequential_time[10],tet[10];
	//host variables
	double *h_omega,*h_delta,*h_a,c_h,f_h,dint,*h_del_fine,*h_omega_fine,*h_diff_delta,*h_diff_omega,*h_fine_tempd,*h_fine_tempo,et[110],*h_correct_delta;
	//device variable
	double *d_omega,*d_delta,*d_a,*d_del_fine,*d_omega_fine,*d_fine_tempd,*d_fine_tempo,*d_diff_delta,*d_diff_omega,*d_correct_delta;
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
	h_diff_delta = new double[num_steps];
	h_diff_omega = new double [num_steps];
	h_correct_delta = new double [num_steps];
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
		gettimeofday(&start,NULL);
		for (int i=0;i<num_steps;i++)
		{
			h_a[i+1]=h_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
			if(h_a[i+1]<=0.8)
			{
				//cout << "a= " <<h_a[i+1]<<__LINE__<<endl;
				//h_a[i+1]=h_a[i]+c_h; //a[i] contains all the time step required for coarse grid calculation
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
		CHECK(cudaMalloc((double**)&d_correct_delta,num_steps_bytes_coarse+8));
		//copying the data to device from host
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(d_delta,h_delta,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_omega,h_omega,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_a,h_a,num_steps_bytes_coarse+8,cudaMemcpyHostToDevice));
		gettimeofday(&end,NULL);
		fMemoryCopyTime[k] = ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		//Kernel call
		//Kernel call
		int ilen = 128;
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
			gpuparareal<<<grid,block>>>(d_delta,d_omega,d_a,omega0,f_h,d_fine_tempd,d_fine_tempo,d_del_fine,d_omega_fine,d_diff_delta,d_diff_omega,num_steps,c_h,d_correct_delta);
			CHECK(cudaEventRecord(kernel_stop));
			CHECK(cudaEventSynchronize(kernel_stop));
			CHECK(cudaEventElapsedTime(&fElapsedTime,kernel_start,kernel_stop));
			et[i]=fElapsedTime;
		}
		CHECK(cudaMemcpy(h_fine_tempd,d_fine_tempd,num_steps_bytes_fine+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_fine_tempo,d_fine_tempo,num_steps_bytes_fine+8,cudaMemcpyDeviceToHost));
//		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(h_diff_delta,d_diff_delta,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_diff_omega,d_diff_omega,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_del_fine,d_del_fine,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(h_omega_fine,d_omega_fine,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		gettimeofday(&start,NULL);
		CHECK(cudaMemcpy(h_correct_delta,d_correct_delta,num_steps_bytes_coarse+8,cudaMemcpyDeviceToHost));
		//CHECK(cudaMemcpy(h_correct_omega,d_correct_omega,num_steps_bytes_coarse=4,cudaMemcpyDeviceToHost));
		gettimeofday(&end,NULL);
	/*	for(int i=1;i<=num_steps;i++)
		{
	//	cout<<(h_correct_delta[i]*180)/pi<<endl;
		}*/
		fMemoryCopyTime[k]+= ((end.tv_sec*1e6+end.tv_usec)-(start.tv_sec*1e6+start.tv_usec))/1000;
		cout<< "Memory transfer time =  " << fMemoryCopyTime[k] <<" ms"<<endl;
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
		CHECK(cudaDeviceReset());
	}
/*	delete[] h_omega;
	delete[] h_delta;
	delete[] h_a;
	delete[] h_del_fine;
	delete[] h_omega_fine;
	delete[] h_fine_tempd;
	delete[] h_fine_tempo;
	delete[] h_diff_delta;
	delete[] h_diff_omega;
*/
}























