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

	double s9_left1[N],s9_right1[N],s9_up1[N],s9_down1[N];
	double s9_left2[N],s9_right2[N],s9_up2[N],s9_down2[N];
	double s9_send_left1[N],s9_send_right1[N],s9_send_up1[N],s9_send_down1[N];
	double s9_send_left2[N],s9_send_right2[N],s9_send_up2[N],s9_send_down2[N];
	double s9_recv_left1[N],s9_recv_right1[N],s9_recv_up1[N],s9_recv_down1[N];
	double s9_recv_left2[N],s9_recv_right2[N],s9_recv_up2[N],s9_recv_down2[N];
	double send_leader_up1[Px][N],send_leader_up2[Px][N],send_leader_down1[Px][N],send_leader_down2[Px][N];
	double recv_leader_up1[Px][N],recv_leader_up2[Px][N],recv_leader_down1[Px][N],recv_leader_down2[Px][N];

    
    if(stencil==9)
    {
		int leader_rank = myrank - myrank%Px	// to know the leader rank in every node
    	sTime=MPI_Wtime();
    	while(t<steps)
    	{
    		
    		
    		
    			int position1=0;
    			int position2=0;
    			int position3=0;
    			int position4=0;

            		if(uflag) //if there is process above
            		{
                		for(int i=0;i<N;i++) //Packing row 0 & 1 
                		{
                    			MPI_Pack(N_matrix+i,1,MPI_DOUBLE,s9_send_up1,N*sizeof(double),&position1,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+N+i,1,MPI_DOUBLE,s9_send_up2,N*sizeof(double),&position2,MPI_COMM_WORLD);
                		}
						if(myrank%Px != 0) // all sending except leader
						{
                		MPI_Isend(s9_send_up1,position1,MPI_PACKED,leader_rank,4*myrank+0,MPI_COMM_WORLD,&request);
                		sends++;
                		MPI_Isend(s9_send_up2,position2,MPI_PACKED,leader_rank,4*myrank+1,MPI_COMM_WORLD,&request);
                		sends++;
						}
            		}
            		if(dflag) //if there is a process below  
            		{
                		position1=0;
                		position2=0;
						position3=0;
                		position4=0;
                		for(int i=0;i<N;i++) //Packing row (N-1) & (N-2)
                		{
                    			MPI_Pack(N_matrix+(N-2)*N+i,1,MPI_DOUBLE,s9_send_down2,N*sizeof(double),&position2,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+(N-1)*N+i,1,MPI_DOUBLE,s9_send_down1,N*sizeof(double),&position1,MPI_COMM_WORLD);
                		}

						if(myrank%Px != 0) // all sending except leader
						{
                		MPI_Isend(s9_send_down1,position2,MPI_PACKED,leader_rank,4*myrank+3,MPI_COMM_WORLD,&request);
                		sends++;
						MPI_Isend(s9_send_down2,position1,MPI_PACKED,leader_rank,4*myrank+2,MPI_COMM_WORLD,&request);
                		sends++;     
						}
            		}

					if(uflag && myrank == leader_rank) // Leader receiving top rows data from other processes
            		{
						for(int i=1; i<Px;i++)
						{
                		MPI_Recv(send_leader_up1[i],N*sizeof(double),MPI_PACKED,myrank+i,4*(myrank+i)+0,MPI_COMM_WORLD,&status);
                		MPI_Recv(send_leader_up2[i],N*sizeof(double),MPI_PACKED,myrank+i,4*(myrank+i)+1,MPI_COMM_WORLD,&status);
						}
            		}

					if(dflag && myrank == leader_rank) // Leader receiving bottom rows data from other processes
            		{
						for(int i=1; i<Px;i++)
						{
                		MPI_Recv(send_leader_down1[i],N*sizeof(double),MPI_PACKED,myrank+i,4*(myrank+i)+3,MPI_COMM_WORLD,&status);
                		MPI_Recv(send_leader_down2[i],N*sizeof(double),MPI_PACKED,myrank+i,4*(myrank+i)+2,MPI_COMM_WORLD,&status);
						}
            		}
					
					// All leaders packing their data for transmission
					if(myrank%Px == 0)
					//Packing row 0, 1, N-2 and N-1
                		{
						position1=0;
                		position2=0;
						position3=0;
                		position4=0;
							for(int i=0;i<N;i++) //Packing row 0 & 1 
                				{
                    			MPI_Pack(N_matrix+i,1,MPI_DOUBLE,send_leader_up1[0],N*sizeof(double),&position1,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+N+i,1,MPI_DOUBLE,send_leader_up2[0],N*sizeof(double),&position2,MPI_COMM_WORLD);
                		
                    			MPI_Pack(N_matrix+(N-2)*N+i,1,MPI_DOUBLE,send_leader_down2[0],N*sizeof(double),&position4,MPI_COMM_WORLD);
                    			MPI_Pack(N_matrix+(N-1)*N+i,1,MPI_DOUBLE,send_leader_down1[0],N*sizeof(double),&position3,MPI_COMM_WORLD);
                				}
							
                		}
					
					// Sending packed data at leader to Up process leader
					if(uflag && myrank == leader_rank)
					{
						for(int i=0; i<Px;i++)
						{		// Using tag 100, 101, 103, 104 for process 0 and so on
                			MPI_Isend(send_leader_up1[i],N,MPI_PACKED,myrank-Px,1000+4*Px*(myrank/Px)+4*i+0,MPI_COMM_WORLD,&request);
                			MPI_Isend(send_leader_up2[i],N,MPI_PACKED,myrank-Px,1000+4*Px*(myrank/Px)+4*i+1,MPI_COMM_WORLD,&request);
							
						}
					}

					// Sending packed data at leader to Down process
					if(dflag && myrank == leader_rank)
					{
						for(int i=0; i<Px;i++)
						{
                			MPI_Isend(send_leader_down1[i],N,MPI_PACKED,myrank+Px,1000+4*Px*(myrank/Px)+4*i+3,MPI_COMM_WORLD,&request);
                			MPI_Isend(send_leader_down2[i],N,MPI_PACKED,myrank+Px,1000+4*Px*(myrank/Px)+4*i+2,MPI_COMM_WORLD,&request);
							
						}
					}

					// Receiving packed data at leader from Up process
					if(uflag && myrank == leader_rank)
					{
						for(int i=0; i<Px;i++)
						{
                			MPI_Recv(recv_leader_up1[i],N*sizeof(double),MPI_PACKED,myrank-Px,1000+4*Px*(myrank/Px-1)+4*i+3,MPI_COMM_WORLD,&status);
                			MPI_Recv(recv_leader_up2[i],N*sizeof(double),MPI_PACKED,myrank-Px,1000+4*Px*(myrank/Px-1)+4*i+2,MPI_COMM_WORLD,&status);
						}
					}

					// Receiving packed data at leader from Down process
					if(dflag && myrank == leader_rank)
					{
						for(int i=0; i<Px;i++)
						{
                			MPI_Recv(recv_leader_down1[i],N*sizeof(double),MPI_PACKED,myrank+Px,1000+4*Px*(myrank/Px+1)+4*i+0,MPI_COMM_WORLD,&status);
                			MPI_Recv(recv_leader_down2[i],N*sizeof(double),MPI_PACKED,myrank+Px,1000+4*Px*(myrank/Px+1)+4*i+1,MPI_COMM_WORLD,&status);
						}
					}

					// Sending received data at leader to other processes on same node 0 to 1, 2, 3 etc

					if(myrank == leader_rank)
					{
						if (uflag && dflag)
						{
						for(int i=1; i<Px;i++)
						{
							MPI_Isend(recv_leader_up1[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+0,MPI_COMM_WORLD,&request);
							MPI_Isend(recv_leader_up2[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+1,MPI_COMM_WORLD,&request);
							MPI_Isend(recv_leader_down1[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+3,MPI_COMM_WORLD,&request);
							MPI_Isend(recv_leader_down2[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+2,MPI_COMM_WORLD,&request);
						}	

						}

						if (uflag == 0)		
						{
							for(int i=1; i<Px;i++)
							{
							MPI_Isend(recv_leader_down1[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+3,MPI_COMM_WORLD,&request);
							MPI_Isend(recv_leader_down2[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+2,MPI_COMM_WORLD,&request);
							}	
						}
						
						if (dflag == 0)		
						{
							for(int i=1; i<Px;i++)
							{
							MPI_Isend(recv_leader_up1[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+0,MPI_COMM_WORLD,&request);
							MPI_Isend(recv_leader_up2[i],N,MPI_PACKED,myrank+i,2000+4*Px*(myrank/Px)+4*i+1,MPI_COMM_WORLD,&request);
							}	
						}
						
					}

					// Receiving data at non leader nodes

					if(myrank != leader_rank)
					{
						if(uflag && dflag)
						{
						MPI_Recv(s9_recv_up1,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+0,MPI_COMM_WORLD,&status);
						MPI_Recv(s9_recv_up2,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+1,MPI_COMM_WORLD,&status);
						MPI_Recv(s9_recv_down1,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+3,MPI_COMM_WORLD,&status);
						MPI_Recv(s9_recv_down2,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+2,MPI_COMM_WORLD,&status);
						}	

						if (uflag == 0)
						{
							MPI_Recv(s9_recv_down1,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+3,MPI_COMM_WORLD,&status);
							MPI_Recv(s9_recv_down2,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+2,MPI_COMM_WORLD,&status);
						}
						
						if (dflag == 0)
						{
							MPI_Recv(s9_recv_up1,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+0,MPI_COMM_WORLD,&status);
							MPI_Recv(s9_recv_up2,N*sizeof(double),MPI_PACKED,leader_rank,2000+4*Px*(myrank/Px)+4*i+1,MPI_COMM_WORLD,&status);
						
						}
					}
					
					// Unpacking top row data at all process

					if(uflag)
            		{
						if(myrank != leader_rank)
						{
						position1=0;
                		position2=0;
						for(int i=0;i<N;i++)
                		{
                		        MPI_Unpack(s9_recv_up2,N,&position2,s9_up2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			MPI_Unpack(s9_recv_up1,N,&position1,s9_up1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);

                		}
						}
						if(myrank == leader_rank) // because we have not sent leaders 0th array so we unpack from leader's  received array directly
						{
						position1=0;
                		position2=0;
						for(int i=0;i<N;i++)
                		{
                		        MPI_Unpack(recv_leader_up2[0],N,&position2,s9_up2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			MPI_Unpack(recv_leader_up1[0],N,&position1,s9_up1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);

                		}
						}

					}
					// Unpacking bottom row data at all process
					if(dflag)
            		{
						if(myrank != leader_rank)
						{
                		position1=0;
                		position2=0;
						for(int i=0;i<N;i++)
                		{
                			MPI_Unpack(s9_recv_down1,N,&position1,s9_down1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    		MPI_Unpack(s9_recv_down2,N,&position2,s9_down2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    			
                		}
						}

						if(myrank == leader_rank) // because we have not sent leaders 0th array so we unpack from leader's received array directly
						{
						position1=0;
                		position2=0;
						for(int i=0;i<N;i++)
                		{
                		    MPI_Unpack(recv_leader_down1[0],N,&position1,s9_down1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    		MPI_Unpack(recv_leader_down2[0],N,&position2,s9_down2+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    		
                		}
						}


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
