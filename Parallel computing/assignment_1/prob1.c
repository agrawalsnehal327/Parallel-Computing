#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"
#include<math.h>

int main(int argc, char* argv[])
{
    MPI_Init(&argc,&argv);

    int mode; // to decide the execution of program
    int stencil;// stencil used(either 5 or 9)
    int seed;// the value needed to initialise the matrix
    int myrank,size_total;
    MPI_Status status;
    MPI_Request request;
    double maxTime;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&size_total);
    int n=atoi(argv[2]);
    int N,steps;
    N = sqrt((double)n);// Can be changed and argument can be given number like 512 instead of 512^2
    steps=atoi(argv[3]);
    seed=atoi(argv[4]);
    stencil=atoi(argv[5]);

    int t=0;
    double sTime,eTime,total_Time,time;
    int sends=0;

    double* N_matrix;
    N_matrix=(double*)malloc(N*N*sizeof(double));
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {

            srand(seed*(myrank+10));
	    N_matrix[i*N+j] = abs(rand()+(i*rand()+j*myrank))/100;
        }
    }
    int Px=atoi(argv[1]);//Number of process in each row of process grid
    int Py=size_total/Px;
    int Pi,Pj;
    Pi = myrank/Px; Pj=myrank%Px;
    int lflag=0,rflag=0,uflag=0,dflag=0;
    if(Pi>0)
    {
       	uflag=1;
    }
    if(Pi<Py-1)
    {
       	dflag=1;
    }
    if(Pj>0)
    {
       	lflag=1;
    }
    if(Pj<Px-1)
    {
	rflag=1;
    }
	double left[N],right[N],up[N],down[N];
    double send_left[N],send_right[N],send_up[N],send_down[N];
    double recv_left[N],recv_right[N],recv_up[N],recv_down[N];
	double s9_left1[N],s9_right1[N],s9_up1[N],s9_down1[N];
	double s9_left2[N],s9_right2[N],s9_up2[N],s9_down2[N];
	double s9_send_left1[N],s9_send_right1[N],s9_send_up1[N],s9_send_down1[N];
	double s9_send_left2[N],s9_send_right2[N],s9_send_up2[N],s9_send_down2[N];
	double s9_recv_left1[N],s9_recv_right1[N],s9_recv_up1[N],s9_recv_down1[N];
	double s9_recv_left2[N],s9_recv_right2[N],s9_recv_up2[N],s9_recv_down2[N];

    if(stencil==5){	
    	sTime = MPI_Wtime();
    	while(t<steps)
    	{
        	

        	
            	int position=0;
            	if(uflag) //if there is process above
            	{
                	for(int i=0;i<N;i++) //sending row 0
                	{
                    		MPI_Pack(N_matrix+i,1,MPI_DOUBLE,send_up,N*sizeof(double),&position,MPI_COMM_WORLD);
                	}
                	MPI_Isend(send_up,position,MPI_PACKED,myrank-Px,myrank,MPI_COMM_WORLD,&request);
                	sends++;
            	}
            	if(dflag) //if there is a process below  
            	{
                	position=0;
                	for(int i=0;i<N;i++) //sending row (N-1)
                	{
                    		MPI_Pack(N_matrix+(N-1)*N+i,1,MPI_DOUBLE,send_down,N*sizeof(double),&position,MPI_COMM_WORLD);
                	}
                	MPI_Isend(send_down,position,MPI_PACKED,myrank+Px,myrank,MPI_COMM_WORLD,&request);
                	sends++;
                
            	}
            	if(lflag) //if there is a process on left
            	{
                	position=0;
                	for(int i=0;i<N;i++) //sending column 0
                	{
                    		MPI_Pack(N_matrix+i*N,1,MPI_DOUBLE,send_left,N*sizeof(double),&position,MPI_COMM_WORLD);
                	}
                	MPI_Isend(send_left,position,MPI_PACKED,myrank-1,myrank,MPI_COMM_WORLD,&request);
                	sends++;
            	}
            	if(rflag) //if there is a process on right
            	{
                	position=0;
                	for(int i=0;i<N;i++) //sending column (n-1)
                	{
                    		MPI_Pack(N_matrix+(i*N+N-1),1,MPI_DOUBLE,send_right,N*sizeof(double),&position,MPI_COMM_WORLD);
                	}
                	sends++;
                	MPI_Isend(send_right,position,MPI_PACKED,myrank+1,myrank,MPI_COMM_WORLD,&request);
            	}
                //receiving data from neighbouring processes
            	if(uflag)
            	{
                	position=0;
                	MPI_Recv(recv_up,N*sizeof(double),MPI_PACKED,myrank-Px,myrank-Px,MPI_COMM_WORLD,&status);
                	for(int i=0;i<N;i++)
                	{
                    		MPI_Unpack(recv_up,N,&position,up+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                	}
            	}
            	if(dflag)
            	{
                	position=0;
                	MPI_Recv(recv_down,N*sizeof(double),MPI_PACKED,myrank+Px,myrank+Px,MPI_COMM_WORLD,&status);
                	for(int i=0;i<N;i++)
                	{
                    		MPI_Unpack(recv_down,N,&position,down+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                	}
            	}
            	if(lflag)
            	{
                	position=0;
                	MPI_Recv(recv_left,N*sizeof(double),MPI_PACKED,myrank-1,myrank-1,MPI_COMM_WORLD,&status);
                	for(int i=0;i<N;i++)
                	{
                    		MPI_Unpack(recv_left,N,&position,left+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                	}
            	}
            	if(rflag)
            	{
                	position=0;
                	MPI_Recv(recv_right,N*sizeof(double),MPI_PACKED,myrank+1,myrank+1,MPI_COMM_WORLD,&status);
                	for(int i=0;i<N;i++)
                	{
                    		MPI_Unpack(recv_right,N,&position,right+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                	}
            	}
        	
			if(Pi==Py-1 && Pj==0){ //For the process in bottom leftmost corner
				N_matrix[0*N+0]=(up[0]+N_matrix[0*N+1]+N_matrix[1*N+0])/3;
				N_matrix[0*N+N-1] = (up[N-1]+right[0]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
				N_matrix[(N-1)*N+0]=(N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/3;
				N_matrix[(N-1)*N+N-1]=(right[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/4;


				for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+up[i]+N_matrix[0*N+i])/stencil;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+N_matrix[(N-1)*N+i])/4;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+0])/4;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+right[i]+N_matrix[i*N+N-1])/stencil;
        		}

			}

			else if(Pi==Py-1 && Pj==Px-1){//For process bottom right most corner(Last process,P-1)
				N_matrix[0*N+0] = (up[0]+left[0]+N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/stencil;
        		N_matrix[0*N+N-1] = (up[N-1]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/4;
        		N_matrix[(N-1)*N+0] = (left[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/4;
        		N_matrix[(N-1)*N+N-1] = (N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/3;        

        		for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+up[i]+N_matrix[0*N+i])/stencil;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+N_matrix[(N-1)*N+i])/4;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+left[i]+N_matrix[i*N+0])/stencil;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1])/4;
        		}

			}

			else if(Pi==0 && Pj==0){//for process in top leftmost corner(Process 0)
				N_matrix[0*N+0] = (N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/3;
        		N_matrix[0*N+N-1] = (right[0]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/4;
        		N_matrix[(N-1)*N+0] = (down[0]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/4;
        		N_matrix[(N-1)*N+N-1] = (down[N-1]+right[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/stencil;        

        		for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+N_matrix[0*N+i])/4;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+down[i]+N_matrix[(N-1)*N+i])/stencil;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+0])/4;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+right[i]+N_matrix[i*N+N-1])/stencil;
        		}
			}

			else if(Pi==0 && Pj==Px-1){//For process top right most corner
				N_matrix[0*N+0] = (left[0]+N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/4;
        		N_matrix[0*N+N-1] = (N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/3;
        		N_matrix[(N-1)*N+0] = (down[0]+left[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/stencil;
        		N_matrix[(N-1)*N+N-1] = (down[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/4;        

        		for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+N_matrix[0*N+i])/4;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+down[i]+N_matrix[(N-1)*N+i])/stencil;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+left[i]+N_matrix[i*N+0])/stencil;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1])/4;
        		}
		
			}

			else if(Pi==0){//For other processes in top row
        
				N_matrix[0*N+0] = (left[0]+N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/4;
        		N_matrix[0*N+N-1] = (right[0]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/4;
        		N_matrix[(N-1)*N+0] = (down[0]+left[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/stencil;
        		N_matrix[(N-1)*N+N-1] = (down[N-1]+right[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/stencil;
        		for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+N_matrix[0*N+i])/4;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+down[i]+N_matrix[(N-1)*N+i])/stencil;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+left[i]+N_matrix[i*N+0])/stencil;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+right[i]+N_matrix[i*N+N-1])/stencil;
        		}
			}
			else if(Pi==Py-1){//for other processes in bottom row
				N_matrix[0*N+0] = (up[0]+left[0]+N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/stencil;
        		N_matrix[0*N+N-1] = (up[N-1]+right[0]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
        		N_matrix[(N-1)*N+0] = (left[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/4;
        		N_matrix[(N-1)*N+N-1] = (right[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/4;        

        		for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+up[i]+N_matrix[0*N+i])/stencil;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+N_matrix[(N-1)*N+i])/4;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+left[i]+N_matrix[i*N+0])/stencil;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+right[i]+N_matrix[i*N+N-1])/stencil;
        		}
			}

			else if(Pj==0){ //For other processes in leftmost column
				N_matrix[0*N+0] = (up[0]+N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/4;
        		N_matrix[0*N+N-1] = (up[N-1]+right[0]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
        		N_matrix[(N-1)*N+0] = (down[0]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/4;
        		N_matrix[(N-1)*N+N-1] = (down[N-1]+right[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/stencil;        

        		for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+up[i]+N_matrix[0*N+i])/stencil;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+down[i]+N_matrix[(N-1)*N+i])/stencil;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+0])/4;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+right[i]+N_matrix[i*N+N-1])/stencil;
        		}
			}

			else if(Pj==Px-1){//for other processes in right most column
				N_matrix[0*N+0] = (up[0]+left[0]+N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/stencil;
        		N_matrix[0*N+N-1] = (up[N-1]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/4;
        		N_matrix[(N-1)*N+0] = (down[0]+left[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/stencil;
        		N_matrix[(N-1)*N+N-1] = (down[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/4;        

        		for(int i=1;i<N-1;i++)
        		{
            		N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+up[i]+N_matrix[0*N+i])/stencil;
            		N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+down[i]+N_matrix[(N-1)*N+i])/stencil;
            		N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+left[i]+N_matrix[i*N+0])/stencil;
            		N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1])/4;
        		}

			}
			else{//For the remaining process
				N_matrix[0*N+0] = (up[0]+left[0]+N_matrix[0*N+1]+N_matrix[1*N+0]+N_matrix[0*N+0])/stencil;
				N_matrix[0*N+N-1] = (up[N-1]+right[0]+N_matrix[0*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
				N_matrix[(N-1)*N+0] = (down[0]+left[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-1)*N+0])/stencil;
				N_matrix[(N-1)*N+N-1] = (down[N-1]+right[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-1])/stencil;        

				for(int i=1;i<N-1;i++)
				{
						N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i-1]+N_matrix[1*N+i]+up[i]+N_matrix[0*N+i])/stencil;
						N_matrix[(N-1)*N+i]= (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i]+down[i]+N_matrix[(N-1)*N+i])/stencil;
						N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[i*N+1]+left[i]+N_matrix[i*N+0])/stencil;
						N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i-1)*N+N-1]+N_matrix[i*N+N-2]+right[i]+N_matrix[i*N+N-1])/stencil;
				}
			}
        	//at interior points(Common for all processes)
        	for(int i=1;i<N-1;i++)
        	{
            		for(int j=1;j<N-1;j++)
            		{
                		N_matrix[i*N+j]=(N_matrix[(i-1)*N+j]+N_matrix[(i+1)*N+j]+N_matrix[i*N+j-1]+N_matrix[i*N+j+1]+N_matrix[i*N+j])/stencil;
            		}
        	}
        	//incrementing the steps
        	t++;
    	}
    	eTime = MPI_Wtime();
        time = eTime-sTime;
    }
    
    
    if(stencil==9)
    {
			
    	sTime=MPI_Wtime();
    	while(t<steps)
    	{
    		
    		
    		
    			int position1=0;
    			int position2=0;
            		if(uflag) //if there is process above
            		{
                		for(int i=0;i<N;i++) //sending row 0
                		{
                    			MPI_Pack(N_matrix+i,1,MPI_DOUBLE,s9_send_up1,N*sizeof(double),&position1,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+N+i,1,MPI_DOUBLE,s9_send_up2,N*sizeof(double),&position2,MPI_COMM_WORLD);
                		}
                		MPI_Isend(s9_send_up1,position1,MPI_PACKED,myrank-Px,myrank,MPI_COMM_WORLD,&request);
                		sends++;
                		MPI_Isend(s9_send_up2,position2,MPI_PACKED,myrank-Px,myrank,MPI_COMM_WORLD,&request);
                		sends++;
            		}
            		if(dflag) //if there is a process below  
            		{
                		position1=0;
                		position2=0;
                		for(int i=0;i<N;i++) //sending row (N-1)
                		{
                    			MPI_Pack(N_matrix+(N-2)*N+i,1,MPI_DOUBLE,s9_send_down2,N*sizeof(double),&position2,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+(N-1)*N+i,1,MPI_DOUBLE,s9_send_down1,N*sizeof(double),&position1,MPI_COMM_WORLD);
                		}
                		MPI_Isend(s9_send_down1,position2,MPI_PACKED,myrank+Px,myrank,MPI_COMM_WORLD,&request);
                		sends++;
						MPI_Isend(s9_send_down2,position1,MPI_PACKED,myrank+Px,myrank,MPI_COMM_WORLD,&request);
                		sends++;     
            		}
            		if(lflag) //if there is a process on left
            		{
                		position1=0;
                		position2=0;
                		for(int i=0;i<N;i++) //sending column 0
                		{
                    			MPI_Pack(N_matrix+i*N,1,MPI_DOUBLE,s9_send_left1,N*sizeof(double),&position1,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+(i*N)+1,1,MPI_DOUBLE,s9_send_left2,N*sizeof(double),&position2,MPI_COMM_WORLD);

                		}
                		MPI_Isend(s9_send_left1,position1,MPI_PACKED,myrank-1,myrank,MPI_COMM_WORLD,&request);
                		sends++;
                		MPI_Isend(s9_send_left2,position2,MPI_PACKED,myrank-1,myrank,MPI_COMM_WORLD,&request);
                		sends++;
            		}
            		if(rflag) //if there is a process on right
            		{
                		position1=0;
                		position2=0;
                		for(int i=0;i<N;i++) //sending column (n-1)
                		{
                    			MPI_Pack(N_matrix+(i*N+N-2),1,MPI_DOUBLE,s9_send_right1,N*sizeof(double),&position1,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+(i*N+N-1),1,MPI_DOUBLE,s9_send_right2,N*sizeof(double),&position2,MPI_COMM_WORLD);
                		}
                		sends++;
                		MPI_Isend(s9_send_right1,position1,MPI_PACKED,myrank+1,myrank,MPI_COMM_WORLD,&request);
                		sends++;
                		MPI_Isend(s9_send_right2,position2,MPI_PACKED,myrank+1,myrank,MPI_COMM_WORLD,&request);
            		}
                        //receiving data from neighbouring processes
            		if(uflag)
            		{
                		position1=0;
                		position2=0;
                		MPI_Recv(s9_recv_up2,N*sizeof(double),MPI_PACKED,myrank-Px,myrank-Px,MPI_COMM_WORLD,&status);
                		MPI_Recv(s9_recv_up1,N*sizeof(double),MPI_PACKED,myrank-Px,myrank-Px,MPI_COMM_WORLD,&status);
                		for(int i=0;i<N;i++)
                		{
                		        MPI_Unpack(s9_recv_up2,N,&position2,s9_up2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			MPI_Unpack(s9_recv_up1,N,&position1,s9_up1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);

                		}
            		}
            		if(dflag)
            		{
                		position1=0;
                		position2=0;
                		MPI_Recv(s9_recv_down1,N*sizeof(double),MPI_PACKED,myrank+Px,myrank+Px,MPI_COMM_WORLD,&status);
                		MPI_Recv(s9_recv_down2,N*sizeof(double),MPI_PACKED,myrank+Px,myrank+Px,MPI_COMM_WORLD,&status);
                		for(int i=0;i<N;i++)
                		{
                			MPI_Unpack(s9_recv_down1,N,&position1,s9_down1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			MPI_Unpack(s9_recv_down2,N,&position2,s9_down2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			
                		}
            		}
            		if(lflag)
            		{
                		position1=0;
                		position2=0;
                		MPI_Recv(s9_recv_left2,N*sizeof(double),MPI_PACKED,myrank-1,myrank-1,MPI_COMM_WORLD,&status);
                		MPI_Recv(s9_recv_left1,N*sizeof(double),MPI_PACKED,myrank-1,myrank-1,MPI_COMM_WORLD,&status);
                		for(int i=0;i<N;i++)
                		{
                    			MPI_Unpack(s9_recv_left1,N,&position1,s9_left1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			MPI_Unpack(s9_recv_left2,N,&position2,s9_left2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                		}
            		}
            		if(rflag)
            		{
                		position1=0;
                		position2=0;
                		MPI_Recv(s9_recv_right1,N*sizeof(double),MPI_PACKED,myrank+1,myrank+1,MPI_COMM_WORLD,&status);
                		MPI_Recv(s9_recv_right2,N*sizeof(double),MPI_PACKED,myrank+1,myrank+1,MPI_COMM_WORLD,&status);
                		for(int i=0;i<N;i++)
                		{
                    			MPI_Unpack(s9_recv_right1,N,&position1,s9_right1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			MPI_Unpack(s9_recv_right2,N,&position2,s9_right2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                		}
            		}

			if(Pi==0 && Pj==0){//For top leftmost process
				N_matrix[0*N+0] = (N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/5;
				N_matrix[0*N+N-1] = (s9_right1[0]+s9_right2[0]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/7;
				N_matrix[(N-1)*N+0] =(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/7;
				N_matrix[(N-1)*N+N-1]=(s9_down1[N-1]+s9_down2[N-1]+s9_right1[N-1]+s9_right2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/stencil;
				
				
				
				N_matrix[0*N+1]=(N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+N_matrix[0*N+1])/6;
				N_matrix[0*N+N-2]=(N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+N_matrix[0*N+N-2])/7;
				N_matrix[(N-1)*N+1]=(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+N_matrix[(N-1)*N+1])/8;
				N_matrix[(N-1)*N+N-2]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+s9_right1[N-1]+N_matrix[(N-1)*N+N-2])/stencil;
				
				
				
				N_matrix[1*N+0]=(N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+s9_up1[0]+N_matrix[1*N+0])/6;
				N_matrix[1*N+N-1]=(s9_right1[0]+s9_right1[1]+N_matrix[0*N+N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/8;
				N_matrix[(N-2)*N+0]=(N_matrix[(N-1)*N+0]+s9_down1[0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/7;
				N_matrix[(N-2)*N+N-1]=(s9_right1[N-2]+s9_right2[N-2]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+s9_down1[N-1]+N_matrix[(N-2)*N+N-1])/stencil;
			
			
			
				N_matrix[1*N+1]=(N_matrix[0*N+1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1])/7;
				N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+s9_right1[1]+N_matrix[0*N+N-2]+N_matrix[1*N+N-2])/8;
				N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_down1[1]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/8;
				N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+s9_down1[N-2]+N_matrix[(N-2)*N+(N-1)]+s9_right1[N-2]+N_matrix[(N-2)*N+(N-2)])/stencil;
			
			
				for(int i=2;i<(N-2);i++)
				{
					N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+N_matrix[0*N+i])/7;
					N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+s9_down1[i]+s9_down2[i]+N_matrix[(N-1)*N+i])/stencil;
					N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+N_matrix[i*N+0])/7;
					N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+s9_right1[i]+s9_right2[i]+N_matrix[i*N+N-1])/stencil;
				
			
			
					N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/8;
					N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+s9_down1[i]+N_matrix[(N-2)*N+i])/stencil;
					N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+N_matrix[i*N+1])/8;
					N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1]+s9_right1[i])/stencil;

				}

			}

			else if(Pi==0 && Pj==Px-1){//For top rightmost process
				N_matrix[0*N+0] = (s9_left1[0]+s9_left2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/7;
				N_matrix[0*N+N-1] = (N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/5;
				N_matrix[(N-1)*N+0] =(s9_down1[0]+s9_down2[0]+s9_left1[N-1]+s9_left2[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/stencil;
				N_matrix[(N-1)*N+N-1]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/7;
					
					
					
				N_matrix[0*N+1]=(N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+s9_left1[0]+N_matrix[0*N+1])/7;
				N_matrix[0*N+N-2]=(N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+N_matrix[0*N+N-2])/6;
				N_matrix[(N-1)*N+1]=(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+s9_left1[N-1]+N_matrix[(N-1)*N+1])/stencil;
				N_matrix[(N-1)*N+N-2]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+N_matrix[(N-1)*N+N-2])/8;
					
					
					
				N_matrix[1*N+0]=(s9_left1[1]+s9_left2[1]+N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+N_matrix[1*N+0])/8;
				N_matrix[1*N+N-1]=(N_matrix[0*N+N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/6;
				N_matrix[(N-2)*N+0]=(s9_left1[N-2]+s9_left2[N-2]+N_matrix[(N-1)*N+0]+s9_down1[0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/stencil;
				N_matrix[(N-2)*N+N-1]=(N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+s9_down1[N-1]+N_matrix[(N-2)*N+N-1])/7;
				
				
				
				N_matrix[1*N+1]=(N_matrix[0*N+1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1]+s9_left1[1])/8;
				N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-2]+N_matrix[1*N+N-2])/7;
				N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_down1[1]+s9_left1[N-2]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/stencil;
				N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+s9_down1[N-2]+N_matrix[(N-2)*N+(N-1)]+N_matrix[(N-2)*N+(N-2)])/8;
				
				
				for(int i=2;i<(N-2);i++)
				{
					N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+N_matrix[0*N+i])/7;
					N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+s9_down1[i]+s9_down2[i]+N_matrix[(N-1)*N+i])/stencil;
					N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+s9_left1[i]+s9_left2[i]+N_matrix[i*N+0])/stencil;
					N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-1])/7;
					
				
				
					N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/8;
					N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+s9_down1[i]+N_matrix[(N-2)*N+i])/stencil;
					N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+s9_left1[i]+N_matrix[i*N+1])/stencil;
					N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1])/8;

					}
			}
			
			else if(Pi==Py-1 && Pj==0){//For bottom left most process
				N_matrix[0*N+0] = (s9_up1[0]+s9_up2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/7;
				N_matrix[0*N+N-1] = (s9_up1[N-1]+s9_up2[N-1]+s9_right1[0]+s9_right2[0]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
				N_matrix[(N-1)*N+0] =(N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/5;
				N_matrix[(N-1)*N+N-1]=(s9_right1[N-1]+s9_right2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/7;
				
				
				
				N_matrix[0*N+1]=(s9_up1[0]+s9_up2[0]+N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+N_matrix[0*N+1])/8;
				N_matrix[0*N+N-2]=(s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+s9_right1[0]+N_matrix[0*N+N-2])/stencil;
				N_matrix[(N-1)*N+1]=(N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+N_matrix[(N-1)*N+1])/6;
				N_matrix[(N-1)*N+N-2]=(N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+s9_right1[N-1]+N_matrix[(N-1)*N+N-2])/7;
				
				
				
				N_matrix[1*N+0]=(N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+N_matrix[1*N+0])/7;
				N_matrix[1*N+N-1]=(s9_right1[0]+s9_right1[1]+N_matrix[0*N+N-1]+s9_up1[N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/stencil;
				N_matrix[(N-2)*N+0]=(N_matrix[(N-1)*N+0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/6;
				N_matrix[(N-2)*N+N-1]=(s9_right1[N-2]+s9_right2[N-2]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+N_matrix[(N-2)*N+N-1])/8;
			
			
			
				N_matrix[1*N+1]=(N_matrix[0*N+1]+s9_up1[1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1])/8;
				N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+s9_right1[1]+N_matrix[0*N+N-2]+s9_up1[N-2]+N_matrix[1*N+N-2])/stencil;
				N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/7;
				N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+N_matrix[(N-2)*N+(N-1)]+s9_right1[N-2]+N_matrix[(N-2)*N+(N-2)])/8;
			
			
				for(int i=2;i<(N-2);i++)
				{
					N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+s9_up1[i]+s9_up2[i]+N_matrix[0*N+i])/stencil;
					N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i])/7;
					N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+N_matrix[i*N+0])/7;
					N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+s9_right1[i]+s9_right2[i]+N_matrix[i*N+N-1])/stencil;
				
			
			
					N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+s9_up1[i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/stencil;
					N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+N_matrix[(N-2)*N+i])/8;
					N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+N_matrix[i*N+1])/8;
					N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1]+s9_right1[i])/stencil;

				}
			}

			else if(Pi==Py-1 && Pj==Px-1){//For bottom rightmost process
				N_matrix[0*N+0] = (s9_up1[0]+s9_up2[0]+s9_left1[0]+s9_left2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/stencil;
				N_matrix[0*N+N-1] = (s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/7;
				N_matrix[(N-1)*N+0] =(s9_left1[N-1]+s9_left2[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/7;
				N_matrix[(N-1)*N+N-1]=(N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/5;
				
				
				
				N_matrix[0*N+1]=(s9_up1[0]+s9_up2[0]+N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+s9_left1[0]+N_matrix[0*N+1])/stencil;
				N_matrix[0*N+N-2]=(s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+N_matrix[0*N+N-2])/8;
				N_matrix[(N-1)*N+1]=(N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+s9_left1[N-1]+N_matrix[(N-1)*N+1])/7;
				N_matrix[(N-1)*N+N-2]=(N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+N_matrix[(N-1)*N+N-2])/6;
				
				
				
				N_matrix[1*N+0]=(s9_left1[1]+s9_left2[1]+N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+s9_up1[0]+N_matrix[1*N+0])/stencil;
				N_matrix[1*N+N-1]=(N_matrix[0*N+N-1]+s9_up1[N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/7;
				N_matrix[(N-2)*N+0]=(s9_left1[N-2]+s9_left2[N-2]+N_matrix[(N-1)*N+0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/8;
				N_matrix[(N-2)*N+N-1]=(N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+N_matrix[(N-2)*N+N-1])/6;
			
			
			
				N_matrix[1*N+1]=(N_matrix[0*N+1]+s9_up1[1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1]+s9_left1[1])/stencil;
				N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-2]+s9_up1[N-2]+N_matrix[1*N+N-2])/8;
				N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_left1[N-2]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/8;
				N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+N_matrix[(N-2)*N+(N-1)]+N_matrix[(N-2)*N+(N-2)])/7;
			
			
				for(int i=2;i<(N-2);i++)
				{
					N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+s9_up1[i]+s9_up2[i]+N_matrix[0*N+i])/stencil;
					N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i])/7;
					N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+s9_left1[i]+s9_left2[i]+N_matrix[i*N+0])/stencil;
					N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-1])/7;
				
			
			
					N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+s9_up1[i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/stencil;
					N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+N_matrix[(N-2)*N+i])/8;
					N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+s9_left1[i]+N_matrix[i*N+1])/stencil;
					N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1])/8;
				}

		}

		else if(Pi==0){//for other processes in top row of grid
			N_matrix[0*N+0] = (s9_left1[0]+s9_left2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/7;
        	N_matrix[0*N+N-1] = (s9_right1[0]+s9_right2[0]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/7;
        	N_matrix[(N-1)*N+0] =(s9_down1[0]+s9_down2[0]+s9_left1[N-1]+s9_left2[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/stencil;
        	N_matrix[(N-1)*N+N-1]=(s9_down1[N-1]+s9_down2[N-1]+s9_right1[N-1]+s9_right2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/stencil;
        	
        	
        	
        	N_matrix[0*N+1]=(N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+s9_left1[0]+N_matrix[0*N+1])/7;
        	N_matrix[0*N+N-2]=(N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+s9_right1[0]+N_matrix[0*N+N-2])/7;
        	N_matrix[(N-1)*N+1]=(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+s9_left1[N-1]+N_matrix[(N-1)*N+1])/stencil;
        	N_matrix[(N-1)*N+N-2]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+s9_right1[N-1]+N_matrix[(N-1)*N+N-2])/stencil;
        	
        	
        	
        	N_matrix[1*N+0]=(s9_left1[1]+s9_left2[1]+N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+N_matrix[1*N+0])/8;
        	N_matrix[1*N+N-1]=(s9_right1[0]+s9_right1[1]+N_matrix[0*N+N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/8;
        	N_matrix[(N-2)*N+0]=(s9_left1[N-2]+s9_left2[N-2]+N_matrix[(N-1)*N+0]+s9_down1[0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/stencil;
	  		N_matrix[(N-2)*N+N-1]=(s9_right1[N-2]+s9_right2[N-2]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+s9_down1[N-1]+N_matrix[(N-2)*N+N-1])/stencil;
	  	
	  	
	  	
	  		N_matrix[1*N+1]=(N_matrix[0*N+1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1]+s9_left1[1])/8;
	  		N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+s9_right1[1]+N_matrix[0*N+N-2]+N_matrix[1*N+N-2])/8;
	  		N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_down1[1]+s9_left1[N-2]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/stencil;
	  		N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+s9_down1[N-2]+N_matrix[(N-2)*N+(N-1)]+s9_right1[N-2]+N_matrix[(N-2)*N+(N-2)])/stencil;
	  	
	  	
	  		for(int i=2;i<(N-2);i++)
        	{
            	N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+N_matrix[0*N+i])/7;
            	N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+s9_down1[i]+s9_down2[i]+N_matrix[(N-1)*N+i])/stencil;
            	N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+s9_left1[i]+s9_left2[i]+N_matrix[i*N+0])/stencil;
            	N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+s9_right1[i]+s9_right2[i]+N_matrix[i*N+N-1])/stencil;
        	
	  	
	  	
	  			N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/8;
	  			N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+s9_down1[i]+N_matrix[(N-2)*N+i])/stencil;
	  			N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+s9_left1[i]+N_matrix[i*N+1])/stencil;
	  			N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1]+s9_right1[i])/stencil;

			}
		}

		else if(Pi==Py-1){//For other processes in bottom row 

			N_matrix[0*N+0] = (s9_up1[0]+s9_up2[0]+s9_left1[0]+s9_left2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/stencil;
        	N_matrix[0*N+N-1] = (s9_up1[N-1]+s9_up2[N-1]+s9_right1[0]+s9_right2[0]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
        	N_matrix[(N-1)*N+0] =(s9_left1[N-1]+s9_left2[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/7;
        	N_matrix[(N-1)*N+N-1]=(s9_right1[N-1]+s9_right2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/7;
        	
        	
        	
        	N_matrix[0*N+1]=(s9_up1[0]+s9_up2[0]+N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+s9_left1[0]+N_matrix[0*N+1])/stencil;
        	N_matrix[0*N+N-2]=(s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+s9_right1[0]+N_matrix[0*N+N-2])/stencil;
        	N_matrix[(N-1)*N+1]=(N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+s9_left1[N-1]+N_matrix[(N-1)*N+1])/7;
        	N_matrix[(N-1)*N+N-2]=(N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+s9_right1[N-1]+N_matrix[(N-1)*N+N-2])/7;
        	
        	
        	
        	N_matrix[1*N+0]=(s9_left1[1]+s9_left2[1]+N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+s9_up1[0]+N_matrix[1*N+0])/stencil;
        	N_matrix[1*N+N-1]=(s9_right1[0]+s9_right1[1]+N_matrix[0*N+N-1]+s9_up1[N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/stencil;
        	N_matrix[(N-2)*N+0]=(s9_left1[N-2]+s9_left2[N-2]+N_matrix[(N-1)*N+0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/8;
	  		N_matrix[(N-2)*N+N-1]=(s9_right1[N-2]+s9_right2[N-2]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+N_matrix[(N-2)*N+N-1])/8;
	  	
	  	
	  	
	  		N_matrix[1*N+1]=(N_matrix[0*N+1]+s9_up1[1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1]+s9_left1[1])/stencil;
	  		N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+s9_right1[1]+N_matrix[0*N+N-2]+s9_up1[N-2]+N_matrix[1*N+N-2])/stencil;
	  		N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_left1[N-2]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/8;
	  		N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+N_matrix[(N-2)*N+(N-1)]+s9_right1[N-2]+N_matrix[(N-2)*N+(N-2)])/8;
	  	
	  	
	  		for(int i=2;i<(N-2);i++)
        	{
            	N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+s9_up1[i]+s9_up2[i]+N_matrix[0*N+i])/stencil;
            	N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i])/7;
            	N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+s9_left1[i]+s9_left2[i]+N_matrix[i*N+0])/stencil;
            	N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+s9_right1[i]+s9_right2[i]+N_matrix[i*N+N-1])/stencil;
        	
	  	
	  	
	  			N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+s9_up1[i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/stencil;
	  			N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+N_matrix[(N-2)*N+i])/8;
	  			N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+s9_left1[i]+N_matrix[i*N+1])/stencil;
	  			N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1]+s9_right1[i])/stencil;

			}
		}
		else if(Pj==0){//For all process in columns 1
			N_matrix[0*N+0] = (s9_up1[0]+s9_up2[0]+s9_left1[0]+s9_left2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/stencil;
        	N_matrix[0*N+N-1] = (s9_up1[N-1]+s9_up2[N-1]+s9_right1[0]+s9_right2[0]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
        	N_matrix[(N-1)*N+0] =(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/7;
        	N_matrix[(N-1)*N+N-1]=(s9_down1[N-1]+s9_down2[N-1]+s9_right1[N-1]+s9_right2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/stencil;
        	
        	
        	
        	N_matrix[0*N+1]=(s9_up1[0]+s9_up2[0]+N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+s9_left1[0]+N_matrix[0*N+1])/stencil;
        	N_matrix[0*N+N-2]=(s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+s9_right1[0]+N_matrix[0*N+N-2])/stencil;
        	N_matrix[(N-1)*N+1]=(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+N_matrix[(N-1)*N+1])/8;
        	N_matrix[(N-1)*N+N-2]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+s9_right1[N-1]+N_matrix[(N-1)*N+N-2])/stencil;
        	
        	
        	
        	N_matrix[1*N+0]=(s9_left1[1]+s9_left2[1]+N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+s9_up1[0]+N_matrix[1*N+0])/stencil;
        	N_matrix[1*N+N-1]=(s9_right1[0]+s9_right1[1]+N_matrix[0*N+N-1]+s9_up1[N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/stencil;
        	N_matrix[(N-2)*N+0]=(N_matrix[(N-1)*N+0]+s9_down1[0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/7;
	  		N_matrix[(N-2)*N+N-1]=(s9_right1[N-2]+s9_right2[N-2]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+s9_down1[N-1]+N_matrix[(N-2)*N+N-1])/stencil;
	  	
	  	
	  	
	  		N_matrix[1*N+1]=(N_matrix[0*N+1]+s9_up1[1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1]+s9_left1[1])/stencil;
	  		N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+s9_right1[1]+N_matrix[0*N+N-2]+s9_up1[N-2]+N_matrix[1*N+N-2])/stencil;
	  		N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_down1[1]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/8;
	  		N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+s9_down1[N-2]+N_matrix[(N-2)*N+(N-1)]+s9_right1[N-2]+N_matrix[(N-2)*N+(N-2)])/stencil;
	  	
	  	
	  		for(int i=2;i<(N-2);i++)
        	{
            	N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+s9_up1[i]+s9_up2[i]+N_matrix[0*N+i])/stencil;
            	N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+s9_down1[i]+s9_down2[i]+N_matrix[(N-1)*N+i])/stencil;
            	N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+N_matrix[i*N+0])/7;
            	N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+s9_right1[i]+s9_right2[i]+N_matrix[i*N+N-1])/stencil;
        	
	  	
	  	
	  			N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+s9_up1[i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/stencil;
	  			N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+s9_down1[i]+N_matrix[(N-2)*N+i])/stencil;
	  			N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+N_matrix[i*N+1])/8;
	  			N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1]+s9_right1[i])/stencil;

			}
		}
		else if(Pj==Px-1){//for all processes in last column
			N_matrix[0*N+0] = (s9_up1[0]+s9_up2[0]+s9_left1[0]+s9_left2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/stencil;
        	N_matrix[0*N+N-1] = (s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/7;
        	N_matrix[(N-1)*N+0] =(s9_down1[0]+s9_down2[0]+s9_left1[N-1]+s9_left2[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/stencil;
        	N_matrix[(N-1)*N+N-1]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/7;
        	
        	
        	
        	N_matrix[0*N+1]=(s9_up1[0]+s9_up2[0]+N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+s9_left1[0]+N_matrix[0*N+1])/stencil;
        	N_matrix[0*N+N-2]=(s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+N_matrix[0*N+N-2])/8;
        	N_matrix[(N-1)*N+1]=(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+s9_left1[N-1]+N_matrix[(N-1)*N+1])/stencil;
        	N_matrix[(N-1)*N+N-2]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+N_matrix[(N-1)*N+N-2])/8;
        	
        	
        	
        	N_matrix[1*N+0]=(s9_left1[1]+s9_left2[1]+N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+s9_up1[0]+N_matrix[1*N+0])/stencil;
        	N_matrix[1*N+N-1]=(N_matrix[0*N+N-1]+s9_up1[N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/7;
        	N_matrix[(N-2)*N+0]=(s9_left1[N-2]+s9_left2[N-2]+N_matrix[(N-1)*N+0]+s9_down1[0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/stencil;
	  		N_matrix[(N-2)*N+N-1]=(N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+s9_down1[N-1]+N_matrix[(N-2)*N+N-1])/7;
	  	
	  	
	  	
	  		N_matrix[1*N+1]=(N_matrix[0*N+1]+s9_up1[1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1]+s9_left1[1])/stencil;
	  		N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+N_matrix[0*N+N-2]+s9_up1[N-2]+N_matrix[1*N+N-2])/8;
	  		N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_down1[1]+s9_left1[N-2]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/stencil;
	  		N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+s9_down1[N-2]+N_matrix[(N-2)*N+(N-1)]+N_matrix[(N-2)*N+(N-2)])/8;
	  	
	  	
	  		for(int i=2;i<(N-2);i++)
        	{
            	N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+s9_up1[i]+s9_up2[i]+N_matrix[0*N+i])/stencil;
            	N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+s9_down1[i]+s9_down2[i]+N_matrix[(N-1)*N+i])/stencil;
            	N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+s9_left1[i]+s9_left2[i]+N_matrix[i*N+0])/stencil;
            	N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-1])/7;
        	
	  	
	  	
	  			N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+s9_up1[i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/stencil;
	  			N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+s9_down1[i]+N_matrix[(N-2)*N+i])/stencil;
	  			N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+s9_left1[i]+N_matrix[i*N+1])/stencil;
	  			N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1])/8;

			}
		}
    	else{	//for all other process
    		N_matrix[0*N+0] = (s9_up1[0]+s9_up2[0]+s9_left1[0]+s9_left2[0]+N_matrix[0*N+1]+N_matrix[0*N+2]+N_matrix[1*N+0]+N_matrix[2*N+0]+N_matrix[0*N+0])/stencil;
        	N_matrix[0*N+N-1] = (s9_up1[N-1]+s9_up2[N-1]+s9_right1[0]+s9_right2[0]+N_matrix[0*N+N-2]+N_matrix[0*N+N-3]+N_matrix[2*N+N-1]+N_matrix[1*N+N-1]+N_matrix[0*N+N-1])/stencil;
        	N_matrix[(N-1)*N+0] =(s9_down1[0]+s9_down2[0]+s9_left1[N-1]+s9_left2[N-1]+N_matrix[(N-1)*N+1]+N_matrix[(N-1)*N+2]+N_matrix[(N-2)*N+0]+N_matrix[(N-3)*N+0]+N_matrix[(N-1)*N+0])/stencil;
        	N_matrix[(N-1)*N+N-1]=(s9_down1[N-1]+s9_down2[N-1]+s9_right1[N-1]+s9_right2[N-1]+N_matrix[(N-2)*N+N-1]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-1)*N+N-2]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-1])/stencil;
        	
        	
        	
        	N_matrix[0*N+1]=(s9_up1[0]+s9_up2[0]+N_matrix[0*N+2]+N_matrix[0*N+3]+N_matrix[1*N+1]+N_matrix[2*N+1]+N_matrix[0*N+0]+s9_left1[0]+N_matrix[0*N+1])/stencil;
        	N_matrix[0*N+N-2]=(s9_up1[N-1]+s9_up2[N-1]+N_matrix[0*N+N-3]+N_matrix[0*N+N-4]+N_matrix[1*N+N-2]+N_matrix[2*N+N-2]+N_matrix[0*N+N-1]+s9_right1[0]+N_matrix[0*N+N-2])/stencil;
        	N_matrix[(N-1)*N+1]=(s9_down1[0]+s9_down2[0]+N_matrix[(N-1)*N+2]+N_matrix[(N-1)*N+3]+N_matrix[(N-2)*N+1]+N_matrix[(N-3)*N+1]+N_matrix[(N-1)*N+0]+s9_left1[N-1]+N_matrix[(N-1)*N+1])/stencil;
        	N_matrix[(N-1)*N+N-2]=(s9_down1[N-1]+s9_down2[N-1]+N_matrix[(N-1)*N+N-3]+N_matrix[(N-1)*N+N-4]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-3)*N+N-2]+N_matrix[(N-1)*N+N-1]+s9_right1[N-1]+N_matrix[(N-1)*N+N-2])/stencil;
        	
        	
        	
        	N_matrix[1*N+0]=(s9_left1[1]+s9_left2[1]+N_matrix[2*N+0]+N_matrix[3*N+0]+N_matrix[1*N+1]+N_matrix[1*N+2]+N_matrix[0*N+0]+s9_up1[0]+N_matrix[1*N+0])/stencil;
        	N_matrix[1*N+N-1]=(s9_right1[0]+s9_right1[1]+N_matrix[0*N+N-1]+s9_up1[N-1]+N_matrix[1*N+N-2]+N_matrix[1*N+N-3]+N_matrix[2*N+N-1]+N_matrix[3*N+N-1]+N_matrix[1*N+N-1])/stencil;
        	N_matrix[(N-2)*N+0]=(s9_left1[N-2]+s9_left2[N-2]+N_matrix[(N-1)*N+0]+s9_down1[0]+N_matrix[(N-2)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-3)*N+0]+N_matrix[(N-4)*N+0]+N_matrix[(N-2)*N+0])/stencil;
	  		N_matrix[(N-2)*N+N-1]=(s9_right1[N-2]+s9_right2[N-2]+N_matrix[(N-3)*N+N-1]+N_matrix[(N-4)*N+N-1]+N_matrix[(N-2)*N+N-2]+N_matrix[(N-2)*N+N-3]+N_matrix[(N-1)*N+N-1]+s9_down1[N-1]+N_matrix[(N-2)*N+N-1])/stencil;
	  	
	  	
	  	
	  		N_matrix[1*N+1]=(N_matrix[0*N+1]+s9_up1[1]+N_matrix[1*N+2]+N_matrix[1*N+3]+N_matrix[2*N+1]+N_matrix[3*N+1]+N_matrix[1*N+0]+N_matrix[1*N+1]+s9_left1[1])/stencil;
	  		N_matrix[1*N+N-2]=(N_matrix[1*N+N-4]+N_matrix[1*N+N-3]+N_matrix[2*N+N-2]+N_matrix[3*N+N-2]+N_matrix[1*N+N-1]+s9_right1[1]+N_matrix[0*N+N-2]+s9_up1[N-2]+N_matrix[1*N+N-2])/stencil;
	  		N_matrix[(N-2)*N+1]=(N_matrix[(N-3)*N+1]+N_matrix[(N-4)*N+1]+N_matrix[(N-2)*N+2]+N_matrix[(N-2)*N+3]+N_matrix[(N-1)*N+1]+s9_down1[1]+s9_left1[N-2]+N_matrix[(N-2)*N+0]+N_matrix[(N-2)*N+1])/stencil;
	  		N_matrix[(N-2)*N+(N-2)]=(N_matrix[(N-3)*N+(N-2)]+N_matrix[(N-4)*N+(N-2)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-2)*N+(N-3)]+N_matrix[(N-1)*N+(N-2)]+s9_down1[N-2]+N_matrix[(N-2)*N+(N-1)]+s9_right1[N-2]+N_matrix[(N-2)*N+(N-2)])/stencil;
	  	
	  	
	  		for(int i=2;i<(N-2);i++)
        	{
            	N_matrix[0*N+i]= (N_matrix[0*N+i+1]+N_matrix[0*N+i+2]+N_matrix[0*N+i-1]+N_matrix[0*N+i-2]+N_matrix[1*N+i]+N_matrix[2*n+i]+s9_up1[i]+s9_up2[i]+N_matrix[0*N+i])/stencil;
            	N_matrix[(N-1)*N+i] = (N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+N_matrix[(N-1)*N+i-1]+N_matrix[(N-1)*N+i-2]+N_matrix[(N-1)*N+i+1]+N_matrix[(N-1)*N+i+2]+s9_down1[i]+s9_down2[i]+N_matrix[(N-1)*N+i])/stencil;
            	N_matrix[i*N+0]= (N_matrix[(i+1)*N+0]+N_matrix[(i+2)*N+0]+N_matrix[(i-1)*N+0]+N_matrix[(i-2)*N+0]+N_matrix[i*N+1]+N_matrix[i*N+2]+s9_left1[i]+s9_left2[i]+N_matrix[i*N+0])/stencil;
            	N_matrix[i*N+N-1]=(N_matrix[(i+1)*N+N-1]+N_matrix[(i+2)*N+N-2]+N_matrix[(i-1)*N+N-1]+N_matrix[(i-2)*N+N-1]+N_matrix[i*N+N-2]+N_matrix[i*N+N-3]+s9_right1[i]+s9_right2[i]+N_matrix[i*N+N-1])/stencil;
        	
	  	
	  	
	  			N_matrix[1*N+i]=(N_matrix[1*N+i+1]+N_matrix[1*N+i+2]+N_matrix[0*N+i]+s9_up1[i]+N_matrix[1*N+i-1]+N_matrix[1*N+i-2]+N_matrix[2*N+i]+N_matrix[3*N+i]+N_matrix[1*N+i])/stencil;
	  			N_matrix[(N-2)*N+i]=(N_matrix[(N-2)*N+(i+1)]+N_matrix[(N-2)*N+(i+2)]+N_matrix[(N-2)*N+(i-1)]+N_matrix[(N-2)*N+(i-2)]+N_matrix[(N-3)*N+i]+N_matrix[(N-4)*N+i]+N_matrix[(N-1)*N+i]+s9_down1[i]+N_matrix[(N-2)*N+i])/stencil;
	  			N_matrix[i*N+1]=(N_matrix[i*N+2]+N_matrix[i*N+3]+N_matrix[(i+1)*N+1]+N_matrix[(i+2)*N+1]+N_matrix[(i-1)*N+1]+N_matrix[(i-2)*N+1]+N_matrix[i*N+0]+s9_left1[i]+N_matrix[i*N+1])/stencil;
	  			N_matrix[i*N+N-2]=(N_matrix[(i+1)*N+(N-2)]+N_matrix[(i+2)*N+N-2]+N_matrix[i*N+N-3]+N_matrix[i*N+N-4]+N_matrix[(i-2)*N+N-2]+N_matrix[(i-1)*N+N-2]+N_matrix[i*N+N-2]+N_matrix[i*N+N-1]+s9_right1[i])/stencil;

			}
		}
        	//at interior points(for all process)
        	for(int i=2;i<N-2;i++)
        	{
            	for(int j=2;j<N-2;j++)
            	{
                	N_matrix[i*N+j]=(N_matrix[(i-1)*N+j]+N_matrix[(i-2)*N+j]+N_matrix[(i+2)*N+j]+N_matrix[(i+1)*N+j]+N_matrix[i*N+j-1]+N_matrix[i*N+j-2]+N_matrix[i*N+j+1]+N_matrix[i*N+j+2]+N_matrix[i*N+j])/stencil;
            	}
        	}
        	//incrementing the steps
        	t++;
    			
    	}
    	eTime = MPI_Wtime();
        time = eTime-sTime;
    }
    	
    MPI_Reduce (&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(myrank==0){
		printf("%lf\n",maxTime);    	
    }
    MPI_Finalize();
    return 0;
}
