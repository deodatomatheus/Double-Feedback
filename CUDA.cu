#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cstring>

using namespace std;

 struct matriz {
    
       int N1, M1;
       int total;
       double* ptr;
    
       matriz(int _N1, int _M1): N1(_N1), M1(_M1), total(_N1*_M1) {
        ptr = new double[_N1*_M1];
        }
    
        ~matriz() {
        delete[] ptr;
        }
    
        double* operator[](int i) {
        return &(ptr[i*M1]);
        }
        };

void hora(char data[80])
  { 
  time_t rawtime;
  struct tm * timeinfo;
  char data_now [80];
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  strftime (data_now,80,"%d-%m--%H-%M-%S",timeinfo);
  puts (data_now);
  sprintf(data,"%s",data_now);
  }

#define c 299792458
#define pi 3.1415926535897932384626433832795
  
#define tn 1e-12                            //Tempo de normalizacao do programa - unidade de t' (unidade de tempo do prog.) = tn segundos
#define dt 0.2 								//Passo temporal em unidades de t' 

#define Go1 3.2e3
#define No1 1.5e8
#define eps1 5.0e-7
#define gammap1 282.0e9  // gammap = (wox-woy)/2  //3.20347e9
#define gammas1 0.5e9
#define alpha1 3.0
#define kappa1 11e9
#define kappa2 11e9

#define wo 0.0 //2.5132741228718345907701147066236e15 //2.217e15

__device__ fun1 (double x1, double y1, double z1, double RE1delay,double RE1delay2) // cos1,cos2,sin1,sin2
{
       return (0.5*tn*((Go1*(z1-No1)/(1.0 + eps1*(x1*x1 + y1*y1))) - gammap1)*(x1-(alpha1*y1)) + 
       tn*kappa1*RE1delay+ tn*kappa2*RE1delay2);
}

__device__ fun2 (double x1, double y1, double z1, double IE1delay, double IE1delay2) // RE1delay,RE1delay2
{
       return (0.5*tn*((Go1*(z1-No1)/(1.0 + eps1*(x1*x1 + y1*y1))) - gammap1)*((alpha1*x1)+y1) + 
       tn*kappa1*IE1delay + tn*kappa2*IE1delay2);
}

__device__ fun3 (double x1, double y1, double z1,double J1)
{
       return (tn*(J1-(z1*gammas1)-((Go1*(z1-No1)/(1.0 + eps1*(x1*x1 + y1*y1))))*(x1*x1 + y1*y1)));
}

__global__ RK (,)
{
	double k11,k12,k13,k21,k22,k23,k31,k32,k33,k41,k42,k43;
	//ATENÇÃO: LAÇO "FOR" DO TRANSIENTE TRANSIENTE//
for (l=0;l<itransiente+1;l=l+1)
{
l1=l1+1;  
		
//MÉTODO RUNGE KUTTA 4ª ORDEM


        k11=dt*fun1(x1[0],x1[1],x1[2],RE1delay,RE1delay2); //
		k12=dt*fun2(x1[0],x1[1],x1[2],IE1delay,IE1delay2); //
		k13=dt*fun3(x1[0],x1[1],x1[2],J1); //

		k21=dt*fun1(x1[0]+k11/2.0,x1[1]+k12/2.0,x1[2]+k13/2.0,RE1delay,RE1delay2); //
		k22=dt*fun2(x1[0]+k11/2.0,x1[1]+k12/2.0,x1[2]+k13/2.0,IE1delay,IE1delay2); //
		k23=dt*fun3(x1[0]+k11/2.0,x1[1]+k12/2.0,x1[2]+k13/2.0,J1);

		k31=dt*fun1(x1[0]+k21/2.0,x1[1]+k22/2.0,x1[2]+k23/2.0,RE1delay,RE1delay2); //
		k32=dt*fun2(x1[0]+k21/2.0,x1[1]+k22/2.0,x1[2]+k23/2.0,IE1delay,IE1delay2); //
		k33=dt*fun3(x1[0]+k21/2.0,x1[1]+k22/2.0,x1[2]+k23/2.0,J1);

		k41=dt*fun1(x1[0]+k31, x1[1]+k32, x1[2]+k33,RE1delay,RE1delay2); //
		k42=dt*fun2(x1[0]+k31, x1[1]+k32, x1[2]+k33,IE1delay,IE1delay2); //
		k43=dt*fun3(x1[0]+k31, x1[1]+k32, x1[2]+k33,J1);

		/*Eaa=Ea; 	    Ea=x[0]; 	    Iaa=Ia; 	    Ia=I;*/
		x1[0]+=(k11+2*k21+2*k31+k41)/6.0;
		x1[1]+=(k12+2*k22+2*k32+k42)/6.0;
		x1[2]+=(k13+2*k23+2*k33+k43)/6.0;
   
	RE1[h11][0]=x1[0];
	IE1[h11][0]=x1[1];
	RE2[h12][0]=x1[0];
	IE2[h12][0]=x1[1];	
  
    RE1delay=RE1[h21][0];
    IE1delay=IE1[h21][0];
    RE1delay2=RE2[h22][0];
    IE1delay2=IE2[h22][0];
    
    h11=h11+1;
    h21=h21+1;
    h12=h12+1;
    h22=h22+1;

        if(h11==idelay1+1)
        {
                      h11=0;
        }
        
         if(h21==idelay1+1)
        {
                      h21=0;
        }
	
	 	if(h12==idelay2+1)
        {
                      h12=0;
        }
        
         if(h22==idelay2+1)
        {
                      h22=0;
        }	

    I1= x1[0]*x1[0]+x1[1]*x1[1];        
    t=t+dt;
 
    Ij=(I1+Ij_1*tfilter*exp(-1.0/tfilter))/tfilter;                                                          
    Ij_1=Ij;
    }

for (l=itransiente;l<itotal+1;l=l+1)
{
    
l1=l1+1;  
		
//MÉTODO RUNGE KUTTA 4ª ORDEM

     //******LASER 1******// k1[0]=k11,k1[1]=k12,k1[2]=k13, k-Tipo do k-Equação a que pertence
     
        //+eta*(x1[2]-x2[2]);

        k11=dt*fun1(x1[0],x1[1],x1[2],RE1delay,RE1delay2); //
		k12=dt*fun2(x1[0],x1[1],x1[2],IE1delay,IE1delay2); //
		k13=dt*fun3(x1[0],x1[1],x1[2],J1); //

		k21=dt*fun1(x1[0]+k11/2.0,x1[1]+k12/2.0,x1[2]+k13/2.0,RE1delay,RE1delay2); //
		k22=dt*fun2(x1[0]+k11/2.0,x1[1]+k12/2.0,x1[2]+k13/2.0,IE1delay,IE1delay2); //
		k23=dt*fun3(x1[0]+k11/2.0,x1[1]+k12/2.0,x1[2]+k13/2.0,J1);

		k31=dt*fun1(x1[0]+k21/2.0,x1[1]+k22/2.0,x1[2]+k23/2.0,RE1delay,RE1delay2); //
		k32=dt*fun2(x1[0]+k21/2.0,x1[1]+k22/2.0,x1[2]+k23/2.0,IE1delay,IE1delay2); //
		k33=dt*fun3(x1[0]+k21/2.0,x1[1]+k22/2.0,x1[2]+k23/2.0,J1);

		k41=dt*fun1(x1[0]+k31, x1[1]+k32, x1[2]+k33,RE1delay,RE1delay2); //
		k42=dt*fun2(x1[0]+k31, x1[1]+k32, x1[2]+k33,IE1delay,IE1delay2); //
		k43=dt*fun3(x1[0]+k31, x1[1]+k32, x1[2]+k33,J1);

		/*Eaa=Ea; 	    Ea=x[0]; 	    Iaa=Ia; 	    Ia=I;*/
		x1[0]+=(k11+2*k21+2*k31+k41)/6.0;
		x1[1]+=(k12+2*k22+2*k32+k42)/6.0;
		x1[2]+=(k13+2*k23+2*k33+k43)/6.0;
		
    
    
    
	RE1[h11][0]=x1[0];
	IE1[h11][0]=x1[1];
	RE2[h12][0]=x1[0];
	IE2[h12][0]=x1[1];	
  
    RE1delay=RE1[h21][0];
    IE1delay=IE1[h21][0];
    RE1delay2=RE2[h22][0];
    IE1delay2=IE2[h22][0];
    
    h11=h11+1;
    h21=h21+1;
    h12=h12+1;
    h22=h22+1;
    
	
        if(h11==idelay1+1)
        {
                      h11=0;
        }
        
         if(h21==idelay1+1)
        {
                      h21=0;
        }
	
	 	if(h12==idelay2+1)
        {
                      h12=0;
        }
        
         if(h22==idelay2+1)
        {
                      h22=0;
        }	
		


    I1= x1[0]*x1[0]+x1[1]*x1[1];
 //   N1 = x1[2];
 //   J11= J1;
	
	
	
    t=t+dt;
    
    Ij=(I1+Ij_1*tfilter*exp(-1.0/tfilter))/tfilter;
    
    Isum=Isum+I1;
    	                                                                            
     if(l1>2000000000&&((l1%jpp)==0))                        
     {                                                                                                                                             
     l1=jpp;                                                                                                                                       
     }                                                       
    
    //if (((long)(l1%jpp)==0)&&(l1>(int)((1000*1E-9/tn)/dt))) // && l1>1000ns(transiente) 
    //if (((long)(l1%jpp)==0)) // && l1>1000ns(transiente) 
    if (((long)(l1%jpp)==0)) 
			{
			    fprintf(p1,"%lf\t%lf\n",t*tn/1e-9,Ij); // Tempo plotado em nanosegundos

			}
    Ij_1=Ij;
    //fprintf(p1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x1[0]*x1[0]+x1[1]*x1[1],x2[0]*x2[0]+x2[1]*x2[1],x1[2],x2[2],t*tn/1e-9,((x1[0]*x1[0]+x1[1]*x1[1])-(x2[0]*x2[0]+x2[1]*x2[1])));

}
}

int main(int argc, char *argv[])
{

	char data [80];
	hora(data);
	
	// DEFINICOES DE COMANDOS NO WINDOWS
	char fazerdiretorio[] = "mkdir";
    char mover[] = "move";
    char copiar[]= "copy";
	
	//PONTEIROS PARA PEGAR NOME DO PROGRAMA
	char *ptr1, *ptr2, *ptr3, *ptr4;
	ptr2 = strtok (argv[0], "\\"); // pointer to first "token" in argv[0]
    while(ptr2 != NULL) { // while current is not null
        ptr1 = ptr2; // copy current to previous
        ptr2 = strtok(NULL, "\\"); // find next
    } 
	   
    ptr2=ptr1;
    ptr2 = strstr (ptr1, ".exe"); // pointer to first "token" in argv[0]
	strncpy(ptr2,".cpp",4); 
	
	//IMPRIMIR NOME DO PROGRAMA
	char formatNomePrograma[]="%s";
	char nome_programa[300];
	sprintf(nome_programa,formatNomePrograma,ptr1);	
	cout<< nome_programa << endl;

	char formatnome_programa2[]="%s-%s";
	char nome_programa2[300];
	sprintf(nome_programa2,formatnome_programa2,data,nome_programa);
		
	char fileNomePrograma[300];
	sprintf(fileNomePrograma,"%s %s %s",copiar,nome_programa,nome_programa2);
	system(fileNomePrograma);
    
    // CRIACAO DA PASTA 0
    strncpy(ptr2,"\0",1);
    char formatNomePasta[]="%s";
	char nome_pasta[200];
	sprintf(nome_pasta,formatNomePasta,ptr1);	
    
    char pasta0[200];   
    sprintf(pasta0,"%s_%s",data,nome_pasta);
    char file00[200];
    sprintf(file00,"%s %s",fazerdiretorio,pasta0);
    system(file00);
    
    //COPIAR UMA COPIA DO PROGRAMA PARA A PASTA 0 //

	char fileNomePrograma2[300];
	sprintf(fileNomePrograma2,"%s %s %s",mover,nome_programa2,pasta0);
	system(fileNomePrograma2);


// DECLARAÇÕES DAS VARIÁVEIS
double t,ttotalreal,ttotal,tmem,transiente,J1,Jo,RE1ret,IE1ret,RE2ret,IE2ret;
int i,j,k,m1,m2,h11,h21,h12,h22,lmax,l1,jpp,cont,tfilter;
int idelay1,idelay2,itransiente;
double itotal,l,Lfeed1,Lfeed2,Lrelativo,Tmed,Tvar,Tdesv,TmedProv,NLFF,ps_p_pontos;
double taudelay1,taudelay2;
//int itotal2;
double taudelay1in, taudelay1fim, deltataudelay1;
// DEFINIÇÕES DOS PARÂMETROS


ttotalreal = 2000000;                     //Tempo de iteração em nanosegundos
ps_p_pontos= 20000;                        //pegar pontos a cada ps_p_pontos pico segundos
tmem = 10000;                               // em unidades de tn (de acordo com a unidade temporal da série bruta calculada. Geralmente em ps)
transiente= 10000; 
tfilter = (int)(tmem/dt);
lmax= (int) tfilter+1;
itransiente=(int)(((transiente*1E-9/tn)/dt)+1);                               //  em unidades de tn   //  ESCOLHER 1ns      Fonte: Programa de John
jpp = (int) ((ps_p_pontos/dt)+1);                                        // para dt=0.5 ps, jpp=20 (pegar ponto a cada 10ps)

taudelay2 =50.0e-9; 									// Tempo de ida e volta na cavidade - Em unidades de segundos // 

taudelay1in = 0.4946*taudelay2;
taudelay1fim= 0.49461*taudelay2;
deltataudelay1=0.00002*taudelay2;

taudelay1 = taudelay1in;

//condições iniciais

double RE1i,RE2i,IE1i,IE2i,N1i,N2i;

RE1i= 1;
IE1i= 1.01;
N1i=  1e8;

//ATENÇÃO: Os valores abaixo NÃO SÃO ARBITRÁRIOS e devem ser calculados com os parâmetros acima

ttotal = (ttotalreal*1E-9)/tn;                  //Tempo de iteração em unidades de t'   
itotal = ceil((ttotal/dt)+1);             	        //Nº total de iterações  
//itotal2 = (int) (itotal);

idelay2 = (int) (((taudelay2/tn)/dt)+1);            //Nº de iterações para feedback laser1  

int l2=0;

FILE *p2;
char filename2[200]; 
char format2[]="%s-LS-filtrada-t1_%.4lf-Lrelativo.txt";
sprintf(filename2,format2,data,taudelay1*1e9);
p2=fopen(filename2,"w"); 

char fileMovp2[200];
sprintf(fileMovp2,"%s %s %s",mover,filename2,pasta0); 


FILE *p2_1;
char filename2_1[200]; 
char format2_1[]="%s-Ibar-filtrada-t1_%.4lf-Lrelativo.txt";
sprintf(filename2_1,format2_1,data,taudelay1*1e9);
p2_1=fopen(filename2_1,"w"); 


char fileMovp2_1[200];
sprintf(fileMovp2_1,"%s %s %s",mover,filename2_1,pasta0); 

//0.17655 - 0.00535
//for(taudelay1=29.527e-9;taudelay1<29.5275e-9;taudelay1=taudelay1+0.001e-9)
//for(taudelay1=26.3541e-9;taudelay1<26.354101e-9;taudelay1=taudelay1+0.00535e-9)
//for(taudelay1=5e-9;taudelay1<=taudelay2;taudelay1=taudelay1+1e-9)
//for(taudelay1=taudelay2*(0.5+pi/300);taudelay1<=taudelay2+1e-9;taudelay1=taudelay1+1e-9)
//{

double Ibar=0, Isum=0;
	
for(taudelay1=taudelay1in;taudelay1<taudelay1fim;taudelay1=taudelay1+deltataudelay1)
{
                                                             
l2=l2+1;

idelay1 = ceil(((taudelay1/tn)/dt)+1);            //Nº de iterações para feedback laser1

Lfeed1=c*taudelay1/2;
Lfeed2=c*taudelay2/2;

Lrelativo=Lfeed1/Lfeed2;  

int dif_idelay;
dif_idelay = idelay1-idelay2;

//cout << "idelay1= " << idelay1 <<endl;
//cout << "idelay2= " << idelay2 <<endl;
//cout << "Lrelativo= " <<Lrelativo  <<endl;



matriz RE1(idelay1+1,1); 
matriz IE1(idelay1+1,1); 
matriz RE2(idelay2+1,1); 
matriz IE2(idelay2+1,1);

double RE1delay;
double IE1delay;
double RE1delay2;
double IE1delay2;


double I1;

//double N1[lmax+1];
//double N2[lmax+1];
//double J11[lmax+1];
//double J22[lmax+1];

// DEFINIÇÕES DOS PARÂMETROS 1


double Nth=No1+(gammap1/(Go1));
double Jth=Nth*gammas1;

//Jo=8e8*1e9;
J1=1.02*Jth;
//cout << "Jo= " << Jo <<endl;
//cout << "Nth= " << Nth <<endl;
//cout << "Jth= " << Jth <<endl;
						
// INICIANDO O PROGRAMA

//data quando o programa começou a rodar ==> para incluir no nome do arquivo
time_t rawtime;
  struct tm * timeinfo;
  char data [80];

  time (&rawtime);
  timeinfo = localtime (&rawtime);

  strftime (data,80,"%d-%m--%H-%M-%S",timeinfo);
  puts (data);
 
  
FILE *p1,*p3;

char format1[]="%s-LS-filtrada-t1_%.4lf.txt";
char filename1[200];
sprintf(filename1,format1,data,taudelay1*1e9);
p1=fopen(filename1,"w"); 

//char s1[] = "Filtro_de_serie-Laser_solitario-parte2";
//char s2[100];
//sprintf(s2,filename1);

char c1[] = "Tlff-hist-dtlff-NumLffBursts-19-03-2018";
//char c2[200];
//sprintf(c2,"filtroLS-%s",filename1);

char d2[300],d3[300],d4[300],d5[300];
sprintf(d2,"TLFF-%s",filename1);
sprintf(d3,"Est-%s",filename1);
sprintf(d4,"dTLFF-%s",filename1);
sprintf(d5,"#LFFs-Burst-%s",filename1);



char fazerdiretorio[] = "mkdir";
char pasta[100];
sprintf(pasta,"tau1-%lf-%s",taudelay1*1e9,data); 

char mover[] = "move";
char del[] = "del";
char destino[100];
sprintf(destino,"%s",pasta);

char file0[100];
sprintf(file0,"%s %s",fazerdiretorio,pasta);

//char file1[500];
//sprintf(file1,"%s %s",s1,s2);

char file2[500];
sprintf(file2,"%s %s",c1,filename1);

char file6[200];
sprintf(file6,"%s %s %s",mover,filename1,destino);

char filedelp1[200];
sprintf(filedelp1,"%s %s /q",del,filename1);

char file7[200];
sprintf(file7,"%s %s %s",mover,d2,destino);

char file8[200];
sprintf(file8,"%s %s %s",mover,d3,destino);

char file9[200];
sprintf(file9,"%s %s %s",mover,d4,destino);

char file9_2[200];
sprintf(file9_2,"%s %s %s",mover,d5,destino);

char fileMovDestino[200];
sprintf(fileMovDestino,"%s %s %s",mover,destino,pasta0);

//cout << "oi5";

//double k11[3], k12[3], k13[3], k14[3];
double x1[3];


//condições iniciais

x1[0]=RE1i; x1[1]=IE1i; x1[2]=N1i;

for(m1=0;m1<idelay1+1;m1++)
{
	RE1[m1][0]=RE1i;
	IE1[m1][0]=IE1i;
}
m1=0;

for(m1=0;m1<idelay2+1;m1++)
{
   	RE2[m1][0]=RE1i;
	IE2[m1][0]=IE1i;
}

t=0.0;

RE1delay=RE1i;
IE1delay=IE1i;
RE1delay2=RE1i;
IE1delay2=IE1i;

    I1= x1[0]*x1[0]+x1[1]*x1[1];
 //   N1[lmax]= x1[2];
 //   J11[lmax]=0;
 //   N2[lmax]= x2[2];
 //   J22[lmax]=0;


h11=idelay1-1;
h21=0;

h12=idelay2-1;
h22=0;

l1=0;

double Ij,Ij_1;
Ij=I1;
Ij_1=I1;


Ibar=Isum/(itotal-itransiente);
cout<<"L1/L2: " <<Lrelativo <<"\t" <<"Ibar: " <<Ibar << endl;
fprintf(p2_1,"%lf\t%lf\n",Lrelativo,Ibar);
Ibar=0;
Isum=0;
//cout << "Calculo feito\n";

fclose(p1);

system(file0); // criar pasta

//cout << endl <<endl <<"ABRINDO PROGRAMA: " <<s1 <<endl;
//system(file1);

cout << endl <<endl <<"ABRINDO PROGRAMA: " <<c1 <<endl;
system(file2);
cout << endl;

//system(file4);
system(file6);

//system(filedelp1);

system(file7);

FILE *leitorEST;
leitorEST=fopen(d3,"r"); 

    for (cont=0;cont<1;cont++)
	{
	char ignore1[1024];
    fgets(ignore1, sizeof(ignore1), leitorEST);
    }
              
   	for (cont=0;cont<1;cont++)
	{
        fscanf(leitorEST,"%lf",&NLFF);      //  coluna 1 do arquivo 1       
		fscanf(leitorEST,"%lf",&Tmed);     //  coluna 2 do arquivo 1 
		fscanf(leitorEST,"%lf",&TmedProv);     //  coluna 2 do arquivo 1 
        fscanf(leitorEST,"%lf",&Tvar);      //  coluna 1 do arquivo 1       
		fscanf(leitorEST,"%lf",&Tdesv);
          
    }

fclose(leitorEST);

system(file8);
system(file9);
system(file9_2);
system(fileMovDestino);

fprintf(p2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Lrelativo,Tmed,TmedProv,Tvar,Tdesv,NLFF);

}

fclose(p2);
system(fileMovp2);

fclose(p2_1);
system(fileMovp2_1);


hora(data);
//fclose(p2);
system("PAUSE");
return 0 ;


}
     
