#include "stdafx.h"
#include "Live.h";
#include "iostream"
#include "vector"
#include "limits"
#include "float.h"
#include "fstream"
#include "string"
#include "iomanip"
#include "cmath"
#include "conio.h"
#include <Windows.h>
using namespace std;

/*int T=85;
int N=5;
int M=32;
vector<long double>pi;
vector<vector<long double>>alpha;
vector<vector<long double>>beta;
vector<vector<long double>>gamma;
vector<vector<long double>>shi;
vector<vector<long double>>delta;
vector<vector<vector<long double>>>zeta;
vector<vector<long double>>A;
vector<vector<long double>>B;
long double p_star=__DBL_MIN__;
vector<long int>q_star;
vector<long int>observation(T+1);
long double prob_O_by_model_lambda = 0;
long double prev_p_star;

void resize()
{
	pi.resize(N+1);
	alpha.resize(T+1,vector<long double>(N+1));
	beta.resize(T+1,vector<long double>(N+1));
	gamma.resize(T+1,vector<long double>(N+1));
	shi.resize(T+1,vector<long double>(N+1));
	delta.resize(T+1,vector<long double>(N+1));
	A.resize(N+1,vector<long double>(N+1));
	B.resize(N+1,vector<long double>(M+1));
    q_star.resize(T+1);
	zeta.resize(T+1,vector<vector<long double>>(N+1,vector<long double>(N+1)));
}*/


int T = 85;
const int N = 5;
const int M = 32;
long double pi[N+1];
long double alpha[161][N+1];
long double beta[161][N+1];
long double gamma[161][N+1];
long double shi[161][N+1];
long double delta[161][N+1];
long double zeta[161][N+1][N+1];
long double A[N+1][N+1];
long double p_star;
long double prev_p_star;
long int q_star[161];
long double B[N+1][M+1];
long int observation[161];
long double prob_O_by_model_lambda = 0;
long double codebook[32][12];

long double ri[13];
long double ai[13];
long double ci[13];
long double raised_sine_ci[13];
int p = 12;
long double weight[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};


// Calculate CI by Raised Sine Window
void calculate_raised_sine_window()
{
	for(int i=1;i<=12;i++)
	{
	  long double pi = 3.141592654;
	  long double ans = 1 + 6 * sin((long double)pi*i/12);
      raised_sine_ci[i] = ci[i]*ans;
	}
}


// Calculate Ci value
void calculate_ci(){
	ci[0] = log(pow(ri[0],2));
	for(int i=1;i<=12;i++){
	    ci[i] = ai[i];
		long double sum = 0;
		for(int j=1;j<=i-1;j++){
		   sum = sum + ((double)j/i)*ci[j]*ai[i-j];
		}
		ci[i] = ci[i]+sum;
	}
}

//Calculate Ai value
void calculate_ai(){
  long double a[13][13],E[13],K[13];
  E[0]=ri[0];
  for(int i=1;i<=p;i++){
     K[i]=ri[i];
	 long double sum = 0;
	 for(int j=1;j<=i-1;j++){
	    sum = sum + a[j][i-1]*ri[i-j];
	 }
	 K[i] = (K[i]-sum)/E[i-1];
	 a[i][i]=K[i];
	 for(int k=1;k<=i-1;k++){
	    a[k][i] = a[k][i-1] - (K[i]*a[i-k][i-1]);
	 }
	 E[i] = (1-pow(K[i],2))*E[i-1];
   }
  for(int i=1;i<=p;i++){
    ai[i] = a[i][12];
  }

}

//Calculate Ri value
void calculate_ri(long double s[]){
	for(int j=1;j<=320;j++){
	  long double w = 0.54 - 0.46 * cos((2*3.1415926535*j)/319);
	  s[j] = s[j] * w; 
	}

	for(int i=0;i<=12;i++){
		long double sum = 0;
		 for(int j=1;j<=320-i;j++){
			sum = sum + (s[j]*s[j+i]);
		}
		ri[i]=sum;
	}
}

int minimum(int a,int b){
	return a<b ? a : b;
}
int maximum(int a,int b){
    return a<b ? b : a;
}

// Dc shift and Normalization
void calculate_dcshift(long double s[])
{
	long double count=0,sum=0;
	int mini = INT_MAX;
	int maxi = INT_MIN;
	 for(int i=1;i<=320;i++){
	   int data = s[i];
	   maxi = maximum(maxi,data);
	   mini = minimum(mini,data);
	   count++;
	   sum = sum + data;
	}

	sum = sum/count;


	 for(int i=1;i<=320;i++){
	   int data = s[i];
	   long double data_w = data + abs(sum);
	   s[i] = data_w;
	}



	long double avg = (abs(mini) + maxi) / (long double)2;
	
      for(int i=1;i<=320;i++){
	   long double data = s[i];
	   long double ans= (data/avg)*5000;
	   s[i]=data;
	}
}




void Forward_Procedure()
{
	prob_O_by_model_lambda = 0;
    for(int i=1;i<=N;i++){
        alpha[1][i]=pi[i]*B[i][observation[1]];
    }

    for(int t=1;t<=T-1;t++)
    {
        for(int j=1;j<=N;j++)
        {
            long double summation=0;
            for(int i=1;i<=N;i++)
            {
                summation = summation + alpha[t][i]*A[i][j];
            }
            alpha[t+1][j] = summation * B[j][observation[t+1]];
        }
    }
    
    for(int i=1;i<=N;i++)
    {
       prob_O_by_model_lambda = prob_O_by_model_lambda + alpha[T][i];
    }

}

void Backward_Procedure()
{

    for(int i=1;i<=N;i++)
    {
        beta[T][i]=1.0;
    }

     for(int t=T-1;t>=1;t--)
    {
        for(int i=1;i<=N;i++)
        {
            long double summation = 0;
            for(int j=1;j<=N;j++)
            {
                summation = summation + (A[i][j] * B[j][observation[t+1]] * beta[t+1][j]);   
            }
            beta[t][i]=summation;
        }
     }
}

void Forward_Backward_Procedure_Gamma()
{
  /*  for(int t=1;t<=T;t++)
    {
        long double summation = 0;
        for (int i = 1; i <= N; i++)
        {
             summation +=  alpha[t][i] * beta[t][i];
        }
        for (int j = 1; j<= N; j++)
        {
            gamma[t][j] = alpha[t][j] * beta[t][j] / summation;
        }
    } */
	for(int t=1;t<=T-1;t++){
		for(int i=1;i<=N;i++){
			long double sum=0;
			for(int j=1;j<=N;j++){
			sum+=zeta[t][i][j];
			}
			gamma[t][i]=sum;
		}
	}
}

void Veterbi_Algorithm()
{
    for(int i=1;i<=N;i++)
    {
       delta[1][i] = pi[i] * B[i][observation[1]];
       shi[1][i]=0;
    }


    for(int t=2;t<=T;t++)
    {
        for(int j=1;j<=N;j++)
        {
            long double maximum = 0;
            long int idx = 0;
            for(int i=1;i<=N;i++)
            {
               long double value = delta[t-1][i] * A[i][j];
               if (maximum < value) {
                     maximum = value;
                     idx=i;
                 } 
            }
            delta[t][j]=(B[j][observation[t]])*maximum;
            shi[t][j]=idx;
        }
    }

    q_star[T]=INT_MIN;
	long double max=0;
    for(int i=1;i<=N;i++)
    {
        if(max<delta[T][i]){
            max=delta[T][i];
            q_star[T]=i;
        }
    }
	p_star=max;


   for(int t=T-1;t>=1;t--){
       q_star[t] = shi[t+1][q_star[t+1]];
   }
}

void Zeta_Calculation()
{
    for(int t=1;t<=T-1;t++)
    {
        long double deno_summation=0.0;
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                deno_summation+=(alpha[t][i]*A[i][j]*B[j][observation[t+1]]*beta[t+1][j]);
				
            }
        }
		
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                zeta[t][i][j]=(alpha[t][i]*A[i][j]*B[j][observation[t+1]]*beta[t+1][j])/(deno_summation);
            }
        }
    }
}

void Reestimation_Problem()
{
	long double const threshold = 1e-30; 
	for(int i=1;i<=N;i++)
	{
		pi[i] = gamma[1][i];
	}
    for(int i=1;i<=N;i++)
    {
		//int mi=0;
       // long double max_value=DBL_MIN;
      //  long double adjust_sum=0;
        for(int j=1;j<=N;j++)
        {
            long double zeta_summation=0,gamma_summation=0;
            for(int t=1;t<=T-1;t++)
            {
                zeta_summation += zeta[t][i][j];
            }
            for(int t=1;t<=T-1;t++)
            {
                gamma_summation += gamma[t][i];
            }
			
            A[i][j] = (long double)(zeta_summation)/(gamma_summation);
			/*if(A[i][j]>max_value){
              max_value=A[i][j];
               mi=j;
            }    
            adjust_sum+=A[i][j];*/
        }
		//A[i][mi]+=(1-adjust_sum);

    }

    for(int j=1;j<=N;j++)
    {
		//int mi=0;
       // long double max_value=DBL_MIN;
        //long double adjust_sum=0;
		int count=0;
		long double max=0;
		int ind_j=0, ind_k=0;
        for(int k=1;k<=M;k++)
        {
            long double gamma_summation=0,deserving_symbol=0;
            for(int t=1;t<=T;t++)
            {
                 if(observation[t]==k)
                 {
                    deserving_symbol += gamma[t][j];
                 }
            }
            for(int t=1;t<=T;t++)
            {
                gamma_summation += gamma[t][j];
            }


            B[j][k] = (deserving_symbol)/(gamma_summation);
			
			if(B[j][k]>max){
				max=B[j][k];
				ind_j = j;
				ind_k = k;
			}
			
			if(B[j][k] == 0){
				B[j][k]=threshold;
				count++;
			}

        }
		 B[ind_j][ind_k] = max - count*threshold;
    }
}

void getObservation()
{
	std::ifstream seq("observation.txt");
	string data;
	int i=1;
	while(!seq.eof())
	{
		/*getline(seq,data);
		observation[i]=stoi(data);*/
		seq >> observation[i];
		i++;
	}
	//T=observation.size()-1;
	seq.close();
}

void getBmatrix()
{
	std::ifstream seq("B_matrix.txt");
	string data;
	int i=1;
	int j=1;
	while(!seq.eof())
	{
		/*getline(seq,data);
		long double a = stold(data);
		B[i][j]=a;*/
		seq >> B[i][j];
		//cout<<B[i][j]<<" ";
		if(j==32)
		{
			//cout<<endl;
		   j=1;
		   i++;
		}
		else j++;
	}
	seq.close();
}

void getAmatrix()
{
	std::ifstream seq("A_matrix.txt");
	string data;
	int i=1;
	int j=1;
	while(!seq.eof())
	{
		/*getline(seq,data);
		long double a = stold(data);
		A[i][j]=double(stold(data));*/
		seq >> long double(A[i][j]);
		if(j==5)
		{
		   j=1;
		   i++;
		}
		else j++;
	}
	seq.close();
}

void getpi()
{
	std::ifstream seq("pi.txt");
	string data;
	int i=1;
	while(!seq.eof())
	{
		seq >>pi[i];
		i++;
	}
	seq.close();
}


void getCodeBook()
{
   ifstream fp("codebook_my.txt");
   int i=0,j=0;
   while(!fp.eof())
   {
	   fp>>codebook[i][j];
	   j++;
	   if(j==12)
	   {
		   i++;
		   j=0;
	   }
   }
   fp.close();
}

int Thokura_Distance()
{
	long double ans = LDBL_MAX;
	int idx=-1;
	for(int i=0;i<32;i++)
	{
      long double dist=0;
	  for(int j=0;j<12;j++)
	  {
		  dist += (long double)(weight[j] * powl((raised_sine_ci[j+1] - codebook[i][j]),2));
	  }
	  if(dist<ans)
	  {
		  ans=dist;
		  idx=i;
	  }
	}
	//cout<<idx<<endl;
	return idx+1;
}

void generate_Observation_Sequence(string file)
{
   ifstream fp(file);
   long double temp[321];
   int i=1;
   int count=1;
   vector<int>ob;
   while(!fp.eof())
   {
	   fp>>temp[i++];
	   if(i==321)
	   {
		   count++;
		   if(count==161){
		     cout<<"Exceed"<<endl;
			 //exit(0);
			 break;
		   }
		  calculate_dcshift(temp);
	      calculate_ri(temp);
		  calculate_ai();
		  calculate_ci();
		  calculate_raised_sine_window();
		  ob.push_back(Thokura_Distance());
		  int k=1;
		 // for(int j=106;j<321;j++)
		//  {
			//  temp[k++]=temp[j];
		  //}
		  
		  i=k;
	   }
     }
       int size = ob.size();
	   for(int j=1;j<=ob.size();j++)
	   {
		   observation[j]=ob[j-1];
	   }
	   T=ob.size();
      fp.close();
}

long double B_temp[N+1][M+1];
long double A_temp[N+1][N+1];
long double prob[10];

void  clean()
{
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			A_temp[i][j]=0;
		}
	}
		for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			B_temp[i][j]=0;
		}
	}
}

void adding_values()
{
	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=N;val++)
		{
			A_temp[var][val] =  long double(A_temp[var][val] + A[var][val]);
		}
	}

	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=M;val++)
		{
			B_temp[var][val] += B[var][val];
		}
	}
}

void avg_values()
{
	for(int var=1;var<=N;var++)
    {
		for(int val=1;val<=N;val++)
		{
				A_temp[var][val] /= 20.0;
				A[var][val]=A_temp[var][val];
		}
	}

	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=M;val++)
		{
				B_temp[var][val] /= 20.0;
				B[var][val]=B_temp[var][val];
		}
	}
}

void store_in_file(string file)
{
	ofstream store(file);
	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=N;val++)
		{
				store << A[var][val] << " ";
		}
		store<<endl;
	}
	store<<endl<<endl<<endl;
	for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=M;val++)
		{
			store << B[var][val] << " ";
		}
		store<<endl;
	}
	store << prob_O_by_model_lambda;
	store.close();
}

void print()
{
  for(int var=1;var<=N;var++)
	{
		for(int val=1;val<=N;val++)
		{
			cout<<A[var][val]<<" ";
		}
		cout<<endl;
	}
  cout<<endl;
}

int correct=0,wrong=0;

void recognizer(int val)
{
	int idx=-1;
	long double maxi=LDBL_MIN;
	for(int i=0;i<=9;i++)
	{
		if(maxi<prob[i])
		{
			idx=i;
			maxi=prob[i];
		}
	}
	cout<<idx<<endl<<endl;
	if(val==idx)
	{
		correct++;
	}
	else{
		wrong++;
	}
}

void find_probability(int val)
{
	for(long long int k=0;k<=9;k++)
	{
		string file = "Data/234101029_E_"+to_string(k)+".txt";
		ifstream fp(file);
		int i=1,j=1;
		while(!fp.eof())
		{
			fp >> A[i][j];
			j++;
			if(j==6)
			{
				i++;
				j=1;
			}
			if(i==6){
			  break;
			}
		}
		string get;
		i=1,j=1;
		while(!fp.eof())
		{
			fp >> B[i][j];
			j++;
			if(j==33)
			{
				i++;
				j=1;
			}
			if(i==6)
			{
				break;
			}
		}
		Forward_Procedure();
	    prob[k] =  prob_O_by_model_lambda;
		fp.close();

	}
	recognizer(val);
}

void testing()
{
	for(long long int i=0;i<=9;i++)
	{
	   for(long long int k=21;k<=30;k++)
	   {
	     string file = "./Sample/234101029_E_"+to_string(i)+"_"+to_string(k)+".txt";
		 cout<<i<<"  "<<"  "<<k<<" -> ";
	     generate_Observation_Sequence(file);
		 find_probability(i); 
	   }
	   cout<<endl;
	}
	cout<<correct<<" "<<wrong<<endl;
}



void training()
{
	getCodeBook();
	for(long long int i=0;i<=9;i++)
	{
		for(long int j=0;j<1;j++)
		{
			clean();
			for(long long int k=1;k<=20;k++)
			{
				if(j==0){
				  getBmatrix();
	              getAmatrix();
                  getpi();
				}
				p_star=0,prev_p_star=-1;
				string file = "./Sample/234101029_E_"+to_string(i)+"_"+to_string(k)+".txt";
				generate_Observation_Sequence(file);
			    int r=1;
			    while(p_star > prev_p_star && r++ && r<1000){
					prev_p_star = p_star;
					Forward_Procedure();
					Backward_Procedure();
					Veterbi_Algorithm();
					Zeta_Calculation();
					Forward_Backward_Procedure_Gamma();
					Reestimation_Problem();
					//print();
			   }
			   adding_values();
			}
			avg_values();
		}
		string file = "./Data/234101029_E_"+to_string(i)+".txt";
		store_in_file(file);

	}
}

void live()
{
 short int  * t = StartRecord();
 ofstream f("./Live Testing/Sample.txt");
 for(int i=0;i<16025*3;i++)
 {
	 f << t[i] <<  endl;
 }
  f.close();
  generate_Observation_Sequence("./Live Testing/Sample.txt");
  find_probability(-2);
}


int main()
{
	cout<<"1.Training  \n2.File Testing \n3.Live Testing\n";
	int a;
	cin>>a;
	switch(a)
	{
	   case 1 : training(); break;
	   case 2 : getpi(); getCodeBook(); testing(); break;
	   case 3 : for(int i=0;i<=9;i++) { getpi(); getCodeBook(); cout<<"Say: "<<i<<endl; Sleep(1); live();} break;
	}
	getch();
	return 0;
}

/*
    
	// training();
	//find_probability();
	//getpi();
	//getCodeBook();
	//testing();
	cout<<"Difference -> ";
	cout<<prev_p_star-p_star<<endl;
    
	cout<<"Iteration "<<r-1<<"  ->   ";
	for(int i=1;i<=T;i++)
	{
		cout<<observation[i]<<" ";
	}
	cout<<endl<<endl<<endl;

    for(int t=1;t<=T;t++){
        cout<<q_star[t]<<" ";
    }

	cout<<endl<<"P_star ->";
	cout<<p_star;
    cout<<endl<<endl;*/
