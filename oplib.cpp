//date      :  20160310,2019/12/15,20191218
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.3)

#include "chash.cpp"
#include "lu.c"


#define DEG 4096
#define K 170
#define T 170  //K/2
#define E 12

unsigned short c[DEG]={0};
unsigned short mat[K][M]={0};
//unsigned short g[K+1]={1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
unsigned short g[K+1]={1,1,0,1,1,0,0,1,1,0,1,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,1};
unsigned short syn[K]={0};
unsigned short BH[K*E][M]={0};
//={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}; //={1,5,0,1,7,3,15}; //={1,2,9,4,0,6,4}; // //
unsigned short P[M][M]={0},inv_P[M][M]={0},S[E*K][E*K]={0},inv_S[E*K][E*K]={0};
unsigned int nn=0;


typedef struct {
  unsigned short n;
  unsigned short a;
} oterm;

typedef struct{
  oterm t[DEG];
} OP;

typedef struct {
  unsigned int x[DEG];
} vec;

OP ss={0};




void random_permutation(unsigned short* a){
	int i,j,x;
	for(i = 0; i < N; i++){
		a[i] = i;
	}
	for(i = 0; i < N - 2; i++){
		j = (rand() % (N-1-i)) + i + 1;

		x = a[j];
		a[j] = a[i];
		a[i] = x;
	}
	if(a[N-1] == N-1){
		a[N-1] = a[N-2];
		a[N-2] = N - 1;
	}


}


unsigned short b2B(unsigned short b[E]){
  int i,j,k;
  unsigned short a=0;

  for(i=0;i<E;i++)
    a^=(b[i]<<i);

  return a;
}


unsigned short oinv(unsigned short a){
  int i;

  for(i=0;i<M;i++){
    if(gf[mlt(fg[a],i)]==1)
      return (unsigned short)i;
  }

}


unsigned short equ(unsigned short a,unsigned short b){
  int i;

  for(i=0;i<M;i++){
    if(gf[mlt(fg[a],fg[i])]==b){
      return i;
    }
  }

}


int deg(vec a){
  int i,n=0;

  for(i=0;i<DEG;i++){
//        printf("%d[%d],",a.x[i],i);
    if(a.x[i]>0)
      n=i;
  }
  
  if(n>0){
    return n;
  }else{
    return 0;
  }
//exit(1);
}


vec o2v(OP f){
  vec a={0};
  int i,count=0;


  for(i=0;i<DEG;i++){
    if(f.t[i].a>0)
      a.x[f.t[i].n]=f.t[i].a;
  }
  //  exit(1);

  return a;
}



OP v2o(vec a){
  int i,count=0;
  OP f={0};

  for(i=0;i<DEG;i++){
    if(a.x[i]>0){
      f.t[i].n=i;
      f.t[i].a=a.x[i];
    }
  }

  return f;
}


OP init_op(OP f){
  int i;

  for(i=0;i<deg(o2v(f))+1;i++){
    f.t[i].a=0;
    f.t[i].n=0;
}
    
  return f;
}


vec i2v(unsigned int n){
  vec v;
  int i;
  
  for(i=0;i<32;i++){
  v.x[i]=n%2;
  n=(n>>1);
  }
	
  return v;
}


unsigned int v2i(vec v){
  unsigned int d=0,i;

  for(i=0;i<32;i++){
    d=(d<<1);
    d^=v.x[i];
  }

  return d;
}


vec init_vec(vec a){
int i;

for(i=0;i<deg(a)+1;i++)
  a.x[i]=0;

return a;
}


vec Setvec(int n){
int i,a,b;
vec v={0};

/*
 c[i]=1;
 for(i=0;i<K;i++)
   c[i]=xor128()%2;
 printf("@\n");
 for(i=0;i<n;i++)
   printf("%d,",c[i]);
 printf("\n");
 //  exit(1);
 */
 for(i=0;i<K+1;i++)
   c[i]=g[i];
 for(i=0;i<n;i++){
  v.x[n-1-i]=c[i];
  printf("%d,",c[n-1-i]);
 }
 printf("\n");
 //   exit(1);
 
  return v;
}


void printpol(vec a){
  int i,n;

  n=deg(a);
  if(n<0){
    printf("baka\n");
    exit(1);
    }
  

  for(i=n;i>-1;i--){
    if(a.x[i]>0){
      printf("%u",a.x[i]);
      if(i>0)
	printf("x^%d",i);
    if(i>0)
      printf("+");
    }
  }
  printf("\n");

  return;
}


OP oadd(OP f,OP g){
  vec a={0},b={0},c={0};
  int i,k;
  OP h={0};

  a=o2v(f);
  b=o2v(g);
  if(deg(a)>deg(b)){
    k=deg(a)+1;
  }else{
    k=deg(b)+1;
}
  printf("k=%d\n",k);
  //  exit(1);
  for(i=0;i<k;i++)
    c.x[i]=a.x[i]^b.x[i];
  h=v2o(c);

  return h;
}


OP oterml(OP f,oterm t){
  int i;
  OP h={0};

  printf("oterml\n");
  printf("deg=%d\n",deg(o2v(f)));
  if(t.a>0){
  for(i=0;i<deg(o2v(f))+1;i++){
    h.t[i].n=f.t[i].n+t.n;
    h.t[i].a=gf[mlt(fg[f.t[i].a],fg[t.a])];
  printf("%dx%d,",h.t[i].a,h.t[i].n);
  }
  }else if(t.a==0){t.n=0;}
  printf("\n");
  //  exit(1);

  return h;
}


OP omul(OP f,OP g){
  int i,count=0;
  oterm t;
  OP h={0},e={0},r={0};


  for(i=0;i<deg(o2v(g))+1;i++){
    t=g.t[i];
    e=oterml(f,t);
    h=oadd(h,e);
  }

  return h;
}


oterm LT(OP f){
  int i;
  oterm t;

  for(i=0;i<deg(o2v(f))+1;i++){
    if(f.t[i].a>0){
      t.n=f.t[i].n;
      t.a=f.t[i].a;
    }
  }

  return t;
}


oterm LTdiv(OP f,oterm t){
  oterm tt,s;


  tt=LT(f);
  if(tt.n<t.n){
    s.n=0;
    s.a=0;
  } else if(tt.n==t.n){
    s.n=0;
    if(t.a>0 && tt.a>0)
      s.a=equ(t.a,tt.a);
    if(t.a==0 || tt.a==0){
      s.a=0;
      s.n=0;
      exit(1);
    }
     
    return s;
  }else if(tt.n>t.n && tt.a>0 && t.a>0){
    s.n=tt.n-t.n;
    s.a=equ(t.a,tt.a);
    return s;
  } else if(tt.n>t.n && tt.a==0){
    s.n=tt.n-t.n;
    exit(1);
    return s;
  }else if(tt.n>t.n && t.n>0 && t.a==0){
    s.n=tt.n-t.n;
    exit(1);
    
    //return s;
  }
    if(tt.n==2 && t.n==1){
      nn++;
      printf("tt=%d,%d %d\n",tt.n,t.n,nn);
       if(nn==3)
      exit(1);
    }
    if(tt.n==1 && t.n==1){
      printf("That's true\n");

      return s;
     
    }

 return s;
}


OP omod(OP f,OP g){
  int i=0,j,n;
  OP h,e;
  oterm a,b,c;

 
  printpol(o2v(f));
  printpol(o2v(g));
  if(deg(o2v(f))==0 || deg(o2v(g))==0){
    printf("baka^\n");
    exit(1);
  }
  printf("in omod\n");
  //exit(1);  

  a=LT(f);
  b=LT(g);
  n=a.n-b.n;
  if(a.a==0 || b.a==0)
    return f;
  while(deg(o2v(f))>=deg(o2v(g))){
    printf("n=======%u\n",deg(o2v(f)));
    printf("m=======%u\n",deg(o2v(g)));
    printf("in omod before LTdiv\n");
    printpol(o2v(f));
    printf("LT(g)=%dx^%d\n",b.a,b.n);
    printf("LT(f)=%dx^%d\n",a.a,a.n);
    c=LTdiv(f,b);
    printf("c=%dx^%d\n",c.a,c.n);
    if(c.a==0)
      exit(1);
    /*
    if(deg(o2v(f))==1 && deg(o2v(g))==1)
        exit(1);
    if(i>4096){
      printf("baka\n");
      exit(1);
    }
    */
    ss.t[i++]=c;
    h=oterml(g,c);
    printf("in omod after oterml\n");
    printpol(o2v(h));
    printpol(o2v(f));
    printpol(o2v(ss));
//exit(1);  
    f=oadd(h,f);
    //f=e;
    printpol(o2v(f));
    //    exit(1);
    a=LT(f);
    b=LT(g);
  }
  //exit(1);
  if(deg(o2v(f))==0)
  return f;
}


OP opow(OP f,int n){
  int i;
  OP g={0};

  memcpy(g.t,f.t,sizeof(f.t));
  for(i=1;i<n;i++)
    g=omul(g,f);

  return g;
}


OP opowmod(OP f,OP mod,int n){
  OP g;
  int i;
  
  g=omod(opow(f,n),mod);
  
  return g;
}


unsigned short trace(OP f,unsigned short x){
  int i;
  unsigned short u=0;

  for(i=0;i<deg(o2v(f))+1;i++)
    u^=gf[mlt(fg[f.t[i].a],mltn(f.t[i].n,fg[x]))];

  return u;
}


// invert of polynomial
OP inv(OP a,OP I){
  OP d,x={0},s={0},q={0},r={0},t={0};

  memcpy(d.t,I.t,deg(o2v(I))+1);
  //  d = I;
  //  x = {0};
  s.t[0].a = 1;
  while (deg(o2v(r)) > 0){
    memset(ss.t,0,DEG);
    r=omod(d , a);
    printpol(o2v(ss));
    //    exit(1);
    memcpy(q.t,ss.t,deg(o2v(ss))+1);
    // r = d % a;
    memcpy(d.t , a.t,deg(o2v(a))+1);
    memcpy(a.t , r.t,deg(o2v(r))+1);
    //a = r;
    t = oadd(x,omul(q , s));
    memcpy(x.t,s.t,deg(o2v(s))+1);
    //    x = s;
    memcpy(s.t,t.t,deg(o2v(t))+1);
    //    s = t;
  }
  //  gcd = d;  // $\gcd(a, n)$
  memset(ss.t,0,DEG);
  omod(I,d);
  
    memset(ss.t,0,DEG);
 return omod(oadd(x , I) , ss);
}


OP ToHorner(OP f){
  OP h;
  int i;
  vec v;

  v=o2v(f);


  return h;
}


OP vx(OP f,OP g){
  OP h={0};
  OP v[K]={0},vv={0};
  oterm a,b;
  int i;

  v[0].t[0].a=0;
  v[0].t[1].n=0;
  v[1].t[0].a=1;
  v[1].t[1].n=0;

  a=LT(f);
  b=LT(g);
  printf("in vx\n");
  //  exit(1);
  for(i=2;i<T+2;i++){
    printf("in vx befoe omod\n");
    printpol(o2v(f));
    printpol(o2v(g));
    //        exit(1);
    memset(ss.t,0,DEG);
  h=omod(f,g);
  printf("ss=");
  printpol(o2v(ss));
  printf("h=");
  printpol(o2v(h));
  //  exit(1);
  v[i]=oadd(v[i-2],omul(ss,v[i-1]));
  printf("-------");
  printpol(o2v(v[i]));
  //  exit(1);
  memset(f.t,0,sizeof(f.t));
  memcpy(f.t,g.t,sizeof(g.t));
  memset(g.t,0,sizeof(g.t));
  memcpy(g.t,h.t,sizeof(h.t));
  //  exit(1);
  if(deg(o2v(v[i]))==T){
        printf("-------");
	//  printpol(o2v(v[i]));
	//  exit(1);
    break;
  }
  }
  memcpy(vv.t,v[i].t,sizeof(v[i].t));
  //      exit(1);
  printf("%d\n",deg(o2v(vv)));
  //    exit(1);
  
  return vv;
}


vec genrandompol(int n){
  arrayul u;
  vec x={0};
  int i,j=0,k;

  //  x=init_vec(x);
  x.x[K+1]=1;
  i=n/N+1;
  while(j<i){
  seed();
  u=crand((unsigned char*)password);
  for(k=0;k<n;k++){
    x.x[j*N+k]=u.d[k];
    printf("%d,",u.d[k]);
  }
  j++;
  }
  //  printf("\n%d\n",deg(x));
  //  exit(1);  

  return x;
}


OP ogcd(OP f,OP g){
  OP h;
  oterm a,b;
  int i=0;

  a=LT(f);
  b=LT(g);
  printpol(o2v(f));
  printpol(o2v(g));
  //  exit(1);
  
  while(i<T){
    memset(ss.t,0,DEG);
    h=omod(f,g);
  printpol(o2v(h));
  //  exit(1);
  f=g;
  g=h;
  //  memcpy(f.t,g.t,sizeof(g.t));
  //memcpy(g.t,h.t,sizeof(h.t));
  a=LT(f);
  printf("%dx%d\n",a.a,a.n);
  b=LT(g);
  printf("%dx%d\n",b.a,b.n);
  i++;
  if(deg(o2v(f))==0 || deg(o2v(g))==0)
    break;
   }
  //  exit(1);

  return h;
}


OP bibun(vec a){
 OP w[DEG]={0},l={0},t={0};
 int i,j,k,n;
 vec tmp={0};
 
 
 n=deg(a);
 for(i=0;i<T;i++){
   w[i].t[0].a=a.x[i];
   w[i].t[0].n=0;
   w[i].t[1].a=1;
   w[i].t[1].n=1;
   printpol(o2v(w[i]));
 }
 //  exit(1);

 for(i=0;i<T;i++){
   tmp.x[0]=1;
   t=v2o(tmp);
   for(j=0;j<T;j++){
     if(i!=j)
       t=omul(t,w[j]);
   }
   printpol(o2v(t));
   //   exit(1);
   l=oadd(l,t);
    printpol(o2v(l));
    //    exit(1);
 }
 printf("l=");
 printpol(o2v(l));
 // printf("%d\n",oinv(2));
 //   exit(1);
 
 return l;
}


vec chen(OP f){
  vec e={0};
  int i,count=0,n,x=0;
  unsigned short y[256]={0},z;

//  e=o2v(f);
n=deg(o2v(f));
  for(x=0;x<M;x++){
    z=0;
    for(i=0;i<n+1;i++){
      z^=gf[mlt(mltn(f.t[i].n,fg[(unsigned short)x]),fg[f.t[i].a])];
    }
    if(z==0)
      e.x[count++]=(unsigned short)x;
    //    printf("%d\n",x);
  }

return e;
}


OP decode(OP f,OP s){
int i,j,k;
 OP r={0},h={0},w,e={0};
oterm t1,t2;
 vec x={0};

 printf("in decode\n");
 r=vx(f,s);
 printpol(o2v(r));
 //   exit(1);

 x=chen(r);
 for(i=0;i<T;i++)
   printf("x=%d ",x.x[i]);
 printf("\n");
 //    exit(1);
 w=bibun(x);
 printpol(o2v(w));
printf("@@@@@@@@@\n");


 h=ogcd(f,s);
 // exit(1);
t1=LT(r);
t2.a=t1.a;
t2.n=0;
w=oterml(w,t2);
printpol(o2v(w));
printpol(o2v(h));
 printf("%d\n",deg(x)+1);
 //  exit(1);
for(i=0;i<deg(x)+1;i++){
  e.t[i].a=gf[mlt(fg[trace(h,x.x[i])],oinv(trace(w,x.x[i])))];
e.t[i].n=x.x[i];
}

 for(i=0;i<T;i++)
   printf("%d ",trace(h,x.x[i]));
 printf("\n");
 for(i=0;i<T;i++)
   printf("%d ",oinv(trace(w,x.x[i])));
 printf("\n");

 
return e;
}


OP setpol(unsigned short f[],int n){
OP g;
vec a;

  memset(c,0,DEG);
  memcpy(c,f,n);
  a=Setvec(n);
  g=v2o(a);


return g;
}


void det(unsigned short g[K+1]){
  OP f,h,w;
  unsigned short cc[K+1]={0},d[2]={0};
  unsigned short HH[K][M]={0};
  int i,j,a,b;
  oterm t={0};
  vec e;
  
    for(i=0;i<K+1;i++)
      printf("%d ",cc[i]);
    printf("\n");

    //#pragma omp parallel for    
  for(i=0;i<M;i++){
    memcpy(cc,g,K+1);
    w=setpol(g,K+1);

    a=trace(w,i);

    
    printf("a=");
    printf("%d ",gf[oinv(a)]);
   
  printf("\n");

  
  for(j=0;j<K+1;j++)
    printf("%d ",cc[j]);
  printf("\n");
    
  cc[K]^=a;
  f=setpol(cc,K+1);
    printpol(o2v(f));
  //  exit(1);
     
  for(j=0;j<K+1;j++)
    printf("%d ",cc[j]);
  printf("\n");
    
  //    exit(1);
  //d[0]=i;
  //d[1]=1;
  h.t[0].a=i;
  h.t[0].n=0;
  h.t[1].a=1;
  h.t[1].n=1;

    memset(ss.t,0,DEG);
  omod(f,h);

  //  printpol(o2v(ss));
  //  exit(1);

  b=oinv(a);
  t.a=gf[b];
  t.n=0;
  //  printf("t=%d\n",t.a);
  //  exit(1);
  w=oterml(ss,t);
  //  printpol(o2v(w));
  e=o2v(w);
  for(j=0;j<K;j++){
    printf("%d,",e.x[K-1-j]);
    HH[j][i]=e.x[K-1-j];
  }
  printf("\n");
    
  //    if(i==1)
  //  exit(1);
  }
  for(i=0;i<K;i++){
    for(j=0;j<M;j++){
      mat[i][j]=HH[i][j];
           printf("%d,",mat[i][j]);
    }
    printf("\n");
  }
  
  //   exit(1);
}


void bdet(unsigned short g[K+1]){
  int i,j,k,l;


  det(g);
  for(i=0;i<M;i++){
    for(j=0;j<K;j++){
      l=mat[j][i];
      for(k=0;k<E;k++){
	BH[j*E+k][i]=l%2;
	l=(l>>1);
      }
    }
  }
  for(i=0;i<M;i++){
    for(j=0;j<E*K;j++)
      printf("%d,",BH[j][i]);
    printf("\n");
  }  
}


int main(int argc,char **argv){
  int i,j,k,l,c;
  unsigned long a,x,count=1;
  //  unsigned short cc[K]={0};
  unsigned short m[K],mm[T]={0};
  time_t timer;
  FILE *fp,*fq;

  unsigned short g2[7]={1,0,9,0,0,6,4};
  //  unsigned short s[K]={0}; //{4,12,7,8,11,13};

  unsigned short ee[10]={1,2,3,4,5,6,7,8,9,10};
  unsigned short zz[T]={0}; //{86,97,114,105,97,98,108,101,32,80,111,108,121,110,111,109};
  //  unsigned short zz[T]={10,97,114,105,97,98,108,101,32,80,111,108,121,110,111,109};
  int y;
  OP f,h,r,w;
  vec v;
  unsigned short d=0;


  srand(clock());
  //  unsigned short syn[K]={4,12,7,8,11,13};
  //  unsigned short g[K+1]={1,0,0,0,1,0,1};

  //  makegf(M);
  //  makefg(M);

  //  zz[0]=1;
  //zz[1]=2;
  //zz[2]=4;
  /*
  g[K]=1;
  for(i=0;i<K;i++)
    g[i]=xor128()%2;
  for(i=0;i<K+1;i++)
    printf("%d,",g[i]);
  */
  printf("\n");
  //  exit(1);
  w=setpol(g,K+1);
  printpol(o2v(w));

  bdet(g);
  //  makeS();
      exit(1);

  /*
  fp=fopen(argv[1],"rb");
  fq=fopen(argv[2],"wb");
  
  
  while((c=fread(zz,1,T,fp))>0){
  
  for(i=0;i<M;i++){
    d=trace(w,(unsigned short)i);
    printf("%d,",d);
    if(d==0){
      printf("%i bad trace 0\n",i);
      exit(1);
    }
  }
  printf("\n");
  //exit(1);
  */

  for(i=0;i<T;i++)
    zz[i]=i+1;
  
   det(g);
  for(i=0;i<K;i++){
    for(j=0;j<M;j++)
      printf("%d,",mat[i][j]);
    printf("\n");
  }
  //exit(1);

  printf("zz=");
  for(i=0;i<K;i++){
    syn[i]=0;

    for(j=0;j<T;j++){
      printf("%u,",zz[j]);
      syn[i]^=gf[mlt(fg[zz[j]],fg[mat[i][j]])];
    }
    //    printf("%d,",syn[i]);
  }
  printf("\n");
  //    exit(1);  
  
  f=setpol(syn,K);
  printpol(o2v(f));
  //  exit(1);
  r=decode(w,f);
  
  for(i=0;i<T;i++){
    mm[i]=r.t[i].a;
    printf("e=%d %d\n",r.t[i].a,r.t[i].n);
  }
  /*
  fwrite(mm,1,T,fq);
  }
  */
//  fclose(fp);
//  fclose(fq);
  
  // h2g(mat);
  //  exit(1);

  return 0;
}
