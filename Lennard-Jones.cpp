#include <iostream>
#include <fstream> 
#include <cmath>
#include <iomanip> 
#include <random>
#include <ctime>
#include <cstdlib>
#include <algorithm>

#define N 16
#define L 4
#define h 0.002
#define t 60
#define pi 3.141592

using namespace std;
void aceleracion(double ax[], double ay[], double x[], double y[]);
bool verificar_distancia(double x[], double y[], int index, double min_distance);


int main()
{
    int i,j,contador,contador2;
    double tx,x[N],y[N],vx[N],vy[N],ax[N],ay[N],wx[N],wy[N],ek,ep,et, dx, dy, r, T,aux,P,xi,yi,rf;
    srand(time(NULL));

    //Abrimos los ficheros
    ofstream posicion("pos.dat");
    ofstream enk("ek.dat");
    ofstream ent("et.dat");
    ofstream enp("ep.dat");
    ofstream hist("hist.dat");
    ofstream histx("histx.dat");
    ofstream histy("histy.dat");
    ofstream temp("temperatura.dat");
    ofstream desp("desp.dat");

    for ( i = 0; i < N; i++)
    {   /*do 
        {
            x[i]=L*1.*(rand()%(RAND_MAX))/RAND_MAX;
            y[i]=L*1.*(rand()%(RAND_MAX))/RAND_MAX;
        } while (!verificar_distancia(x, y, i, 0.18*L)); // Verificar la distancia con las partículas anteriores
        vx[i]=pow(1,rand()%(2))*4.*(rand()%(RAND_MAX))/RAND_MAX;
        vy[i]=pow(-1,rand()%(2))*sqrt(-pow(vx[i],2)+16);*/

        vx[i]=0.;
        vy[i]=0.;
        
    }
    x[1]=0.5;
    y[1]=0.5;
    x[2]=1.5;
    y[2]=0.5;
    x[3]=2.5;
    y[3]=0.5;
    x[4]=3.5;
    y[4]=0.5;
    x[5]=0.5;
    y[5]=1.5;
    x[6]=1.5;
    y[6]=1.5;
    x[7]=2.5;
    y[7]=1.5;
    x[8]=3.5;
    y[8]=1.5;
    x[9]=0.5;
    y[9]=2.5;
    x[10]=1.5;
    y[10]=2.5;
    x[11]=2.5;
    y[11]=2.5;
    x[12]=3.5;
    y[12]=2.5;  
    x[13]=0.5;
    y[13]=3.5;
    x[14]=1.5;
    y[14]=3.5;
    x[15]=2.5;
    y[15]=3.5;
    x[0]=3.5;
    y[0]=3.5;
    
    xi=x[6];
    yi=y[6];
    contador=0;
    contador2=0;
    T=0;
    P=0;


    aceleracion(ax,ay,x,y);

    for (tx=0;tx<t;tx=tx+h) //iniciamos las iter
    {
        //if(tx>=4)
    //{
        aux=0;
        //contador2=contador2+1;
       for ( i = 0; i < N; i++)
        {
            aux=aux+(pow(vx[i],2)+pow(vy[i],2))/2.;
            //hist << sqrt(pow(vx[i],2)+pow(vy[i],2)) << endl;
            //histx << vx[i] << endl;
            //histy << vy[i] << endl;
        }
        aux=aux/(N*1.);
        T=T+aux;
        
    //}
    if (fmod(contador, 10)==0 || tx==0)
    {
        //Energia
        ek=0.;
        ep=0.;
        et=0.;
        for ( i = 0; i < N; i++)
        {
            ek=ek+(pow(vx[i],2)+pow(vy[i],2))/2.;
        }

    for (i=0;i<N;i++)
    {
        for (j=0;j<i;j++)
        {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            r=sqrt(pow(dx,2)+pow(dy,2));
            if (abs(dx)<=(L/2.) && abs(dy)<=(L/2.))
            {
                ep=ep+4.*(1./pow(r,12)-1./pow(r,6));
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)>0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx),2)+pow((dy-L),2)),12)-1./pow(sqrt(pow((dx),2)+pow((dy-L),2)),6));
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)<0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx),2)+pow((dy+L),2)),12)-1./pow(sqrt(pow((dx),2)+pow((dy+L),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)>0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx-L),2)+pow((dy),2)),12)-1./pow(sqrt(pow((dx-L),2)+pow((dy),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)<0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx+L),2)+pow((dy),2)),12)-1./pow(sqrt(pow((dx+L),2)+pow((dy),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)>0.)
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx-L),2)+pow((dy-L),2)),12)-1./pow(sqrt(pow((dx-L),2)+pow((dy-L),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)<0.)
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx-L),2)+pow((dy+L),2)),12)-1./pow(sqrt(pow((dx-L),2)+pow((dy+L),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)<0. && (dy)>0.)
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx+L),2)+pow((dy-L),2)),12)-1./pow(sqrt(pow((dx+L),2)+pow((dy-L),2)),6));
            }
            else 
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx+L),2)+pow((dy+L),2)),12)-1./pow(sqrt(pow((dx+L),2)+pow((dy+L),2)),6));
            }
        }   
        for (j=i+1;j<N;j++)
        {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            r=sqrt(pow(dx,2)+pow(dy,2));
            if (abs(dx)<=(L/2.) && abs(dy)<=(L/2.))
            {
                ep=ep+4.*(1./pow(r,12)-1./pow(r,6));
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)>0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx),2)+pow((dy-L),2)),12)-1./pow(sqrt(pow((dx),2)+pow((dy-L),2)),6));
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)<0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx),2)+pow((dy+L),2)),12)-1./pow(sqrt(pow((dx),2)+pow((dy+L),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)>0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx-L),2)+pow((dy),2)),12)-1./pow(sqrt(pow((dx-L),2)+pow((dy),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)<0.)
            {
                ep=ep+4.*(1./pow(sqrt(pow((dx+L),2)+pow((dy),2)),12)-1./pow(sqrt(pow((dx+L),2)+pow((dy),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)>0.)
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx-L),2)+pow((dy-L),2)),12)-1./pow(sqrt(pow((dx-L),2)+pow((dy-L),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)<0.)
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx-L),2)+pow((dy+L),2)),12)-1./pow(sqrt(pow((dx-L),2)+pow((dy+L),2)),6));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)<0. && (dy)>0.)
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx+L),2)+pow((dy-L),2)),12)-1./pow(sqrt(pow((dx+L),2)+pow((dy-L),2)),6));
            }
            else 
            {
               ep=ep+4.*(1./pow(sqrt(pow((dx+L),2)+pow((dy+L),2)),12)-1./pow(sqrt(pow((dx+L),2)+pow((dy+L),2)),6));
            }
        } 
  
    }
    ep=ep/2.;
    et=ep+ek;
    enk << tx << "   "<< ek << endl; 
    enp << tx << "   "<< ep << endl;    
    ent << tx << "   "<< et << endl; 

            for (i=0;i<N;i++)
        {
            posicion << x[i] << ",   " << y[i] << endl;
        }
        /*if (contador>5000)
        {
           
        
        
for (i=1;i<N;i++)
    {
            dx=x[i]-x[0];
            dy=y[i]-y[0];
            r=sqrt(pow(dx,2)+pow(dy,2));
            if (abs(dx)<=(L/2.) && abs(dy)<=(L/2.))
            {
                hist << r << endl;
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)>0.)
            {
                hist << sqrt(pow((dx),2)+pow((dy-L),2)) << endl;
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)<0.)
            {
                hist << sqrt(pow((dx),2)+pow((dy+L),2)) << endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)>0.)
            {
                hist << sqrt(pow((dx-L),2)+pow((dy),2)) << endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)<0.)
            {
                hist << sqrt(pow((dx+L),2)+pow((dy),2)) << endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)>0.)
            {
               hist << sqrt(pow((dx-L),2)+pow((dy-L),2))<< endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)<0.)
            {
              hist << sqrt(pow((dx-L),2)+pow((dy+L),2))<< endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)<0. && (dy)>0.)
            {
               hist << sqrt(pow((dx+L),2)+pow((dy-L),2)) << endl;
            }
            else 
            {
               hist << sqrt(pow((dx+L),2)+pow((dy+L),2))<< endl;
            }
      
    }
    }*/
            posicion << endl;
            
    }

        //Calculamos la posición en t+h
        for(i=0;i<N;i++)
        {
            x[i]=x[i]+h*vx[i]+pow(h,2)*ax[i]/2;
            wx[i]=vx[i]+h*ax[i]/2;
            y[i]=y[i]+h*vy[i]+pow(h,2)*ay[i]/2;
            wy[i]=vy[i]+h*ay[i]/2;
        }   

                for ( i = 0; i < N; i++)
        {
            if (x[i]>L) 
            {
                x[i]=x[i]-L*1.;
                P=P+2.0*abs(vx[i]);
            }

            if (x[i]<0) 
            {
                x[i]=x[i]+L*1.;
                P=P+2.0*abs(vx[i]);
            }

            if (y[i]>L) 
            {
                y[i]=y[i]-L*1.;
                P=P+2.0*abs(vy[i]);
            }

            if (y[i]<0) 
            {
                y[i]=y[i]+L*1.;
                P=P+2.0*abs(vy[i]);
            }
        }

        //Calculamos la new aceleracion en t+h 

        aceleracion(ax,ay,x,y);

        //Calculamos las nuevas v en t+h
        for(i=0;i<N;i++)
        {
        vx[i]=wx[i]+h*(ax[i])/2;
        vy[i]=wy[i]+h*(ay[i])/2;
        }
       /* if (contador==5000)
        {
            for(i=0;i<N;i++)
        {
        vx[i]=vx[i]*3.0;
        vy[i]=vy[i]*3.0;
        }
        }*/

        if (contador==9999||contador==14999||contador==17499||contador==22499)
        {
        for(i=0;i<N;i++)
        {
        vx[i]=vx[i]*1.5;
        vy[i]=vy[i]*1.5;
        }
        temp << T/((contador-contador2)*1.) << endl;
        contador2=contador;
        T=0;
        }

        /*if (contador==59999||contador==89999||contador==119999||contador==29999||contador==149999)
        {
        for(i=0;i<N;i++)
        {
        vx[i]=vx[i]*1.1;
        vy[i]=vy[i]*1.1;
        }
        temp << T/((contador-contador2)*1.) << endl;
        contador2=contador;
        T=0;
        }*/
        
           /* dx=x[6]-xi;
            dy=y[6]-yi;
            r=sqrt(pow(dx,2)+pow(dy,2));
            if (abs(dx)<=(L/2.) && abs(dy)<=(L/2.))
            {
                desp << tx <<"    "<< pow(r,2) << endl;
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)>0.)
            {
                desp << tx <<"    "<< (pow((dx),2)+pow((dy-L),2)) << endl;
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)<0.)
            {
                desp << tx <<"    "<< (pow((dx),2)+pow((dy+L),2)) << endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)>0.)
            {
                desp << tx <<"    "<< (pow((dx-L),2)+pow((dy),2)) << endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)<0.)
            {
                desp << tx <<"    "<< (pow((dx+L),2)+pow((dy),2)) << endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)>0.)
            {
               desp << tx <<"    "<< (pow((dx-L),2)+pow((dy-L),2))<< endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)<0.)
            {
              desp << tx <<"    "<< (pow((dx-L),2)+pow((dy+L),2))<< endl;
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)<0. && (dy)>0.)
            {
               desp << tx <<"    "<< (pow((dx+L),2)+pow((dy-L),2)) << endl;
            }
            else 
            {
               desp << tx <<"    "<< (pow((dx+L),2)+pow((dy+L),2))<< endl;
            }*/
        contador=contador+1;
    }
    //temp << T/(contador2*1.); //<< "   " << P/(4.*L*50.);

    return 0;
}

void aceleracion (double ax[], double ay[], double x[], double y[])
{
    int i,j;
    double r, dx, dy;

    for (i=0;i<N;i++)
    {
        ax[i]=0;
        ay[i]=0;
    }
    
    for (i=0;i<N;i++)
    {
        for (j=0;j<i;j++)
        {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            r=sqrt(pow(dx,2)+pow(dy,2));
            if (abs(dx)<=(L/2.) && abs(dy)<=(L/2.))
            {
                ax[i]=ax[i]+48.*dx/pow(r,14)-24.*dx/pow(r,8);
                ay[i]=ay[i]+48.*dy/pow(r,14)-24.*dy/pow(r,8);
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)<0.)
            {
                ax[i]=ax[i]+24.*(2.*(dx/pow(sqrt(pow(dx,2)+pow((dy+L),2)),14))-(dx/pow(sqrt(pow(dx,2)+pow((dy+L),2)),8)));
                ay[i]=ay[i]+24.*(2.*((dy+L)/pow(sqrt(pow(dx,2)+pow((dy+L),2)),14))-((dy+L)/pow(sqrt(pow(dx,2)+pow((dy+L),2)),8)));
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),14))-((dx)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy-L)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),14))-((dy-L)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),14))-((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),14))-((dy)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)<0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),14))-((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),14))-((dy)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),14))-((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),14))-((dy-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)<0. && (dy)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),14))-((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy-L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),14))-((dy-L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)<0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),14))-((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy+L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),14))-((dy+L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),8))));
            }
            else
            {
                ax[i]=(ax[i]+24.*(2.*((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),14))-((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),14))-((dy+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),8))));
            }
        }   
        for (j=i+1;j<N;j++)
        {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            r=sqrt(pow(dx,2)+pow(dy,2));
            if (abs(dx)<=(L/2.) && abs(dy)<=(L/2.))
            {
                ax[i]=ax[i]+48.*dx/pow(r,14)-24.*dx/pow(r,8);
                ay[i]=ay[i]+48.*dy/pow(r,14)-24.*dy/pow(r,8);
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)<0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx)/pow(abs(sqrt(pow((dx),2)+pow((dy+L),2))),14))-((dx)/pow(abs(sqrt(pow((dx),2)+pow((dy+L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy+L)/pow(abs(sqrt(pow((dx),2)+pow((dy+L),2))),14))-((dy+L)/pow(abs(sqrt(pow((dx),2)+pow((dy+L),2))),8))));
            }
            else if (abs(dx)<=(L/2.) && abs(dy)>(L/2.) && (dy)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),14))-((dx)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy-L)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),14))-((dy-L)/pow(abs(sqrt(pow((dx),2)+pow((dy-L),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),14))-((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),14))-((dy)/pow(abs(sqrt(pow((dx-L),2)+pow((dy),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)<=(L/2.) && (dx)<0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),14))-((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),14))-((dy)/pow(abs(sqrt(pow((dx+L),2)+pow((dy),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),14))-((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),14))-((dy-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy-L),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)<0. && (dy)>0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),14))-((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy-L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),14))-((dy-L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy-L),2))),8))));
            }
            else if (abs(dx)>(L/2.) && abs(dy)>(L/2.) && (dx)>0. && (dy)<0.)
            {
                ax[i]=(ax[i]+24.*(2.*((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),14))-((dx-L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy+L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),14))-((dy+L)/pow(abs(sqrt(pow((dx-L),2)+pow((dy+L),2))),8))));
            }
            else
            {
                ax[i]=(ax[i]+24.*(2.*((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),14))-((dx+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),8))));
                ay[i]=(ay[i]+24.*(2.*((dy+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),14))-((dy+L)/pow(abs(sqrt(pow((dx+L),2)+pow((dy+L),2))),8))));
            }
        }   
    }
}

bool verificar_distancia(double x[], double y[],int j, double min_distance) {
    for (int q = 0; q < j; q++) {
        if (sqrt(pow((x[q]-x[j]),2)+pow((y[q]-y[j]),2)) < min_distance || sqrt(pow((abs(x[q]-x[j])-L),2)+pow((y[q]-y[j]),2)) < min_distance || sqrt(pow((x[q]-x[j]),2)+pow((abs(y[q]-y[j])-L),2)) < min_distance || sqrt(pow((abs(x[q]-x[j])-L),2)+pow((abs(y[q]-y[j])-L),2)) < min_distance ) {
            return false; // La distancia es menor que el umbral mínimo
        }
    }
    return true; // La distancia es mayor que el umbral mínimo para todas las partículas
}