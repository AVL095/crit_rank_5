#include <string.h>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <random> 

using namespace std;

void moodstatistic(int k,int *m,long int num,vector <double> &h);
void kruskalstatistic(int k,int *m,long int num,vector <double> &h);
void lemanstatistic(int k,int *m,long int num,vector <double> &h);
void rank_distribution(vector <double> h,vector <double> &pw,vector <double> &wrange);

///////////////////////////////////////////////////////////////////

void rank_distribution(vector <double> h,vector <double> &pw,vector <double> &wrange) {

   long int nn,nnew,num,j;
   double s;

   num=h.size();
   sort(h.begin(),h.end());
  
    nn=0;nnew=0;
    for (j=0;j<num;j++) {
        if (roundf(h[j]*100000)/100000==roundf(h[j+1]*100000)/100000) {
            nnew++;
        }
        else {
            pw.push_back(nnew+1);
            wrange.push_back(h[j]);
            nn++;
            nnew=0;
        }
    }
    s=0.;
    for(j=0;j<nn;j++) {
        s+= pw[j]/double(num); 
        pw[j]=s;
    }


}


//////////////////////////////////////////////////////////////

void moodstatistic(int k,int *m,long int num,vector <double> &h) {
    int i,ii,n;
    double z;
    vector<int>r;
    
    n=m[0]+m[1];
    for(i=0;i<n;i++) r.push_back(i+1);
     random_device rd;
     mt19937 g(rd());

    for(ii=0;ii<num;ii++) { 
      shuffle(r.begin(),r.end(), g); 
      z= 0;
      for(i=0;i<m[0];i++) z=z+pow((r[i]-(n+1.)/2.),2);
      h.push_back(z);
   }
   r.clear();
}

////////////////////////////////////////////////////////////////

void lemanstatistic(int k,int *m,long int num,vector <double> &h) {
    int i,ii,n;
    double z,r1,r2;
    vector <int> x;
    vector <int> y;
    vector <int> r;

    n=m[0]+m[1];
    for(i=0;i<n;i++) r.push_back(i+1);
     random_device rd;
     mt19937 g(rd());  


    for(ii=0;ii<num;ii++) { 
      shuffle(r.begin(),r.end(), g); 
      for (i=0;i<m[0];i++) x.push_back(r[i]);
      for (i=0;i<m[1];i++) y.push_back(r[i+m[0]]);
      sort(x.begin(),x.end());
      sort(y.begin(),y.end());
      r1=0.;r2=0.;
      for (i=0;i<m[0];i++) r1+=pow(x[i]-(i+1),2);
      for (i=0;i<m[1];i++) r2+=pow(y[i]-(i+1),2);
      z=(r1*1.0/m[1]+r2*1.0/m[0]+0.1666666666666)/(1.0*m[0]*m[1])-0.6666666666;
      h.push_back(z);
      x.clear();y.clear();
    }

    r.clear();
}


///////////////////////////////////////////////////////////////////////

void kruskalstatistic(int k,int *m,long int num,vector <double> &h) {
    int i,j,ii,n,km;
    vector <int> r;
    double s,rr,z;

    n=0;
    for(i=0;i<k;i++) n+=m[i];

    for(i=0;i<n;i++) r.push_back(i+1);
     random_device rd;
     mt19937 g(rd());

    for(ii=0;ii<num;ii++) {
       shuffle(r.begin(),r.end(), g);
       s=0.;km=0;
      for (j=0;j<k;j++) {
        rr=0.;
        for (i=0;i<m[j];i++) rr+=r[i+km];
          s+=rr*rr/m[j];
          km+=m[j];
        }
      z=12.*s/(n*(n+1.))-3.*(n+1.);
      h.push_back(z);
   }

    r.clear();

 }


/////////////////////////////////////////////////////

int main() {

  int i,k,*m;
  long int j,num,nn;
  string st,ff;
  vector <double>h;
  vector <double>pw;
  vector <double>wrange;
  double s;

  ifstream inp1("main.inp");
  inp1>>ff;
  inp1.close();
 
  ifstream inp("Inp/" + ff + ".inp");
  ofstream out("Out/" + ff + ".out");

  inp>>st;
  inp>>k;
  m=new int[k];
  inp>>st;
  for(i=0;i<k;i++) inp>>m[i];
  inp>>st;
  inp>>num;
  inp.close();

  out<<"Criterion:"<<ff<<endl;
  out<<"Variants="<<num<<endl;
  for(i=0;i<k;i++) out<<m[i]<<";";
  out <<endl;

  if(ff=="Kruskal_exact") kruskalstatistic(k,m,num,h);
  if(ff=="Leman_exact") lemanstatistic(k,m,num,h);
  if(ff=="Mood_exact") moodstatistic(k,m,num,h);
     
  rank_distribution(h,pw,wrange);

  nn=pw.size();
  out << "Size=" << nn << endl;
  for (j=0;j<nn;j++) out<<(j+1)<<";"<<setprecision(12)<<fixed<<wrange[j]<<";"<<pw[j]<<endl;
  out <<endl;
  out.close();

  h.clear();
  wrange.clear();
  pw.clear();

  return 0;
}