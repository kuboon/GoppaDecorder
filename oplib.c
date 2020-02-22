 
//date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.8)
// gcdの停止条件を修正した。vxの停止条件を修正した。
//date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.8)
// 0ベクトルが出ないように生成多項式のトレースチェックを入れた。
   
//date      :  20160310
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.1)

#include "chash.cpp"
#include "lu.c"


#define K 16*2
#define DEG 3*K
#define T K/2
#define E 8
#define D 256

unsigned char tmp[E*K][D]={0};
unsigned char pub[E*K][D]={0};
unsigned char BH[E*K][D]={0};
unsigned short c[2*K+1]={0};
unsigned short mat[K][D]={0};
unsigned short m2[K][D]={0};

static unsigned short g[K+1]={0};
  /*
  {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       0,1};
  */
  //  unsigned short g[K+1]={1,1,0,1,1,0,0,1,1,0,1};
//unsigned short g[K+1]={1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};

//,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
//unsigned short g[K+1]={1,0,0,0,1,0,1};
unsigned short syn[K]={0};
unsigned char A[D][D]={0};
//={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}; //={1,5,0,1,7,3,15}; //={1,2,9,4,0,6,4}; // //
unsigned short P[D]={0};
unsigned short inv_P[D]={0};
unsigned short uu;


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

typedef struct {
  OP u;
  OP v;
  OP d;
} EX;

typedef union {
  unsigned long long int u[4];
  unsigned short s[16];
} SU;

OP ss={0};
OP hh={0};
OP qq={0};
OP p0={0};


static void ginit(void){
  int i,j,count=0;
  unsigned short gg[K+1]={0};
  
  printf("in ginit\n");
  //   g[128]=g[126]=g[124]=g[120]=g[118]=g[117]=g[116]=g[115]=g[114]=g[113]=g[112]=g[110]=g[108]=g[106]=g[105]=g[104]=g[102]=g[101]=g[99]=g[98]=g[97]=g[94]=g[93]=g[89]=g[88]=g[85]=g[82]=g[81]=g[78]=g[77]=g[76]=g[74]=g[72]=g[69]=g[68]=g[66]=g[62]=g[61]=g[60]=g[57]=g[56]=g[54]=g[52]=g[51]=g[49]=g[48]=g[47]=g[46]=g[41]=g[40]=g[36]=g[35]=g[32]=g[30]=g[28]=g[27]=g[26]=g[24]=g[19]=g[17]=g[13]=g[10]=g[8]=g[6]=g[4]=g[3]=g[2]=g[1]=g[0]=1;
       
  g[K]=1;
  g[0]=1;
  while(count<K/2-1){
    printf("@\n");
    j=rand()%(K-1);
    if(j<K && j>0 && g[j]==0){
      g[j]=1;
      count++;
	}
  }
  
  
  for(i=0;i<K+1;i++)
    gg[i]=g[K-i];
  for(i=0;i<K+1;i++)
    g[i]=gg[i];
  
}




void random_permutation(unsigned short* a){
	int i,j,x;
	for(i = 0; i < D; i++){
		a[i] = i;
	}
	for(i = 0; i < D - 2; i++){
		j = (rand() % (D-1-i)) + i + 1;

		x = a[j];
		a[j] = a[i];
		a[i] = x;
	}
	if(a[D-1] == D-1){
		a[D-1] = a[D-2];
		a[D-2] = D - 1;
	}


}


void P2Mat(unsigned short P[D]){
int i,j;
	
  for(i=0;i<D;i++)
      A[i][P[i]]=1;
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
    if(gf[mlt(fg[a],fg[i])]==b)
      break;
  }
  return i;
}


int distance(OP f){
int i,j,k;

for(i=0;i<DEG;i++){
if(f.t[i].a>0)
  j=i;
}

return j;
}

int terms(OP f){
int i,count=0;

for(i=0;i<DEG;i++)
  if(f.t[i].a>0)
    count++;

return count;
}


int odeg(OP f){
int i,j=0,k;

for(i=0;i<terms(f)+1;i++){
  if(j<f.t[i].n && f.t[i].a>0)
    j=f.t[i].n;
}

return j;
}



int deg(vec a){
  int i,n=0;

  for(i=0;i<DEG;i++){
    if(a.x[i]>0)
      n=i;
  }
  
  if(n>0)
    return n;


}


vec o2v(OP f){
  vec a={0};
  int i,count=0;


  for(i=0;i<DEG;i++){
    if(f.t[i].a>0)
      a.x[f.t[i].n]=f.t[i].a;
  }

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


OP init_op(OP f){
  int i;

  for(i=0;i<deg(o2v(f))+1;i++){
    f.t[i].a=0;
    f.t[i].n=0;
}
    
  return f;
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


 for(i=0;i<n;i++){
  v.x[n-1-i]=c[i];
 }

 
  return v;
}



void oprintpol(OP f){
int i,j,k,n;

n=distance(f);
printf("n=%d\n",n);
printf("terms=%d\n",terms(f));
printf("deg=%d\n",odeg(f));
 
//exit(1);

for(i=n;i>-1;i--){
  if(f.t[i].a>0)
    printf("%ux^%u+",f.t[i].a,f.t[i].n);
}
//printf("\n");
//exit(1);
return;
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
  //  printf("\n");

  return;
}


OP oadd(OP f,OP g){
  vec a={0},b={0},c={0};
  int i,k;
  OP h={0};

  a=o2v(f);
  //exit(1);
  b=o2v(g);
  //  oprintpol((g));
  //  exit(1);
  if(deg(a)>=deg(b)){
    k=deg(a)+1;
  }else{

    k=deg(b)+1;

}
  for(i=0;i<k;i++)
    c.x[i]=a.x[i]^b.x[i];
  h=v2o(c);

  return h;
}


OP oterml(OP f,oterm t){
  int i,k;
  OP h={0};
  vec test;
  unsigned short n;

  k=distance(f);
  for(i=0;i<k+1;i++){

    h.t[i].n=f.t[i].n+t.n;
    h.t[i].a=gf[mlt(fg[f.t[i].a],fg[t.a])];

  }

  return h;
}


OP omul(OP f,OP g){
  int i,count=0,k;
  oterm t={0};
  OP h={0},e={0},r={0};

  if(deg(o2v(f))>deg(o2v(g))){
    k=deg(o2v(f));
  }else{
    k=deg(o2v(g));
  }

  for(i=0;i<k+1;i++){
    t=g.t[i];
    e=oterml(f,t);
    h=oadd(h,e);
  }

  printpol(o2v(h));
  printf(" debug======\n");
  //  exit(1);
  
  return h;
}


oterm LT(OP f){
  int i,k;
  oterm t={0};

  k=distance(f);
  for(i=0;i<k+1;i++){
    //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if(f.t[i].a>0){
      t.n=f.t[i].n;
      t.a=f.t[i].a;
    }
  }

  return t;
}


oterm LT2(OP f){
  int i,k;
  oterm t={0};

      t.n=f.t[0].n;
      t.a=f.t[0].a;

  return t;
}


oterm LTdiv(OP f,oterm t){
  oterm tt={0},s={0};

  tt=LT(f);
  if(tt.n<t.n){
    s.n=0;
    s.a=0;
  } else if(tt.n==t.n){
    s.n=0;
    s.a=equ(t.a,tt.a);
  }else if(tt.n>t.n){
    s.n=tt.n-t.n;
    s.a=equ(t.a,tt.a);
    //printf("%u\n",s.a);
  }else if(t.n==0 && t.a>0){
    s.a=gf[mlt(fg[tt.a],oinv(t.a))];
    s.n=tt.n;
  }

  return s;    
}


OP omod(OP f,OP g){
  int i=0,j,n,k;
  OP h={0},e={0};
  oterm a,b={0},c={0};
  
  
  n=LT(g).n;
  if(deg(o2v(f))==0){ 
    printf("baka^\n");
    //return f;
    exit(1);
  }
  if(deg(o2v(g))==0){
    printf("baka--\n");
    //return g;
         exit(1);
  }
  if(deg(o2v(f))<deg(o2v(g))){
    //    exit(1);
    return f;
  }

  //  printf("in omod\n");
  //exit(1);
  k=deg(o2v(g));
  b=LT(g);
  
  

    printf("b=========%dx^%d\n",b.a,b.n);
    oprintpol((g));
    printf("\nin omod g=============%d\n",deg(o2v(g)));
    while((deg(o2v(f))-deg(o2v(g)))>=0){
      printf("in!\n");
      //    exit(1);
      
      c=LTdiv(f,b);
      //   printf("c========%dx^%d\n",c.a,c.n);
      //    exit(1);
      //ss.t[i]=c;
      //i++;
      oprintpol((g));
      printf("\ng=================\n");
      
      h=oterml(g,c);
      oprintpol((h));
      printf("\n");
      printf("modh===================\n");

      oprintpol((f));
      printf("\nmodf===================\n");
      //     exit(1);
      
      f=oadd(f,h);
      if(odeg(f)>0)
      oprintpol((f));
      printf("\nff1=====================\n");
      if(deg(o2v(f))==0 || deg(o2v(f))>100){
	printf("baka500\n");
	//break;
      }


    }
    
    return f;
}


OP odiv(OP f,OP g){
  int i=0,j,n,k;
  OP h={0},e={0},tt={0},o={0};
  oterm a,b={0},c={0};

  o=f; 

  if(deg(o2v(f))==0 && g.t[0].a==0){ 
    printf("baka^\n");
    //return f;
        exit(1);
  }
  if(deg(o2v(g))==0 && g.t[0].a==0){
    printf("baka--\n");
    exit(1);
  }
    printf("in odiv\n");
    //oprintpol((g));
    
  //exit(1);
    k=deg(o2v(g));
    b=LT(g);
    oprintpol((g));
    printf("in odiv1 g===========%d %d\n",b.a,b.n);
  
    if(deg(o2v(f))<deg(o2v(g))){
  return f;
  //  a=LT(f);
    }
  printf("odiv in b=========%dx^%d\n",b.a,b.n);
oprintpol((g));
 printf("\nf===================\n");
 // e=omul(g,g);
 //oprintpol((e));
 //printf("\ng^2================\n");

 printf("\nin odiv2 g=============%d\n",deg(o2v(g)));
 i=0;
 while((deg(o2v(f))-deg(o2v(g)))>=0){
    //  printf("in!\n");
    //    exit(1);
    
    c=LTdiv(f,b);
    //    if(c.a>0){
    // printf("in odev c========%dx^%d\n",c.a,c.n);
    //    exit(1);
    tt.t[c.n]=c;
    
    //oprintpol((g));
    //printf("\ng=================\n");

    h=oterml(g,c);
    //oprintpol((h));
    //printf("\n");
    //printf("h===================\n");
    //oprintpol((f));
    //printf("\nf===================\n");
    //     exit(1);
    
      f=oadd(f,h);
    
      //oprintpol((f));
      //printf("\nff=====================\n");
      if(deg(o2v(f))>500){
	printf("baka500\n");
	exit(1);
      }
      
      if(c.n==0)
	break;
  }
  //exit(1);
  


    return tt;
}




  
  
OP opow(OP f,int n){
  int i;
  OP g={0};

  g=f;
  //memcpy(g.t,f.t,sizeof(f.t));
  
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

  for(i=0;i<deg(o2v(f))+1;i++){
    u^=gf[mlt(fg[f.t[i].a],mltn(f.t[i].n,fg[x]))];
  }
  
  return u;
}


OP coeff(OP f,unsigned short d){
  int i,j,k;
  vec a,b;
  
  for(i=0;i<deg(o2v(f))+1;i++)
    f.t[i].a=gf[mlt(fg[f.t[i].a],oinv(d))];


  return f;

}


// invert of polynomial
OP inv(OP a,OP n){
  OP d={0},x={0},s={0},q={0},r={0},t={0},u={0},v={0},w={0},tt={0},gcd={0};
  oterm b={0};
  vec vv={0},xx={0};

  
  if(deg(o2v(a))>deg(o2v(n))){
    printf("baka_i\n");
    exit(1);
  }
  //  v=n;
  printf("n=============\n");
  oprintpol((n));
  printf("\n");
  printf("a=============\n");
  oprintpol((a));
  printf("\n");
  //  exit(1);
  
  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (LT(a).n>0){

    r=omod(d ,a);
    q=odiv(d, a);
    
    d = a;
    a = r;
    t=oadd(x,omul(q,s));
    oprintpol((a));
    printf("\nin roop a==================%d\n",deg(o2v(a)));
    printf("\n");
    
    x = s;
    s = t;
  }
  // exit(1);
  if(a.t[0].a>0){
    d = a;
    a = r;
    oprintpol((a));
    printf("\nin roop a==================%d\n",deg(o2v(a)));
    printf("\n");
    
    x = s;
    s = t;

   oprintpol((d));
   printf("\nout1================\n");
     gcd = d;  // $\gcd(a, n)$
     oprintpol((gcd));
   printf("\n");
   oprintpol((n));
   printf("\n");
   printf("out2===============\n");
   
   printf("before odiv\n");
   w=n;

   b=LT(w);
   oprintpol((w));
   printf("\nw=======%d %d\n",b.a,b.n);
   //w=tt;
   v=oadd(x,n);
   oprintpol((v));
   printf("\n");
    printf("d==============\n");
  }
    oprintpol((w));
   printf("\n");
   printf("ss==============\n");
   //       exit(1);
   u=omod(v, w);
   u=odiv(u, d);
   //u=coeff(u,d.t[0].a);
   oprintpol((u));
   printf("\nu==================\n");
   //   exit(1); 
   return u;

}


OP ToHorner(OP f){
  OP h;
  int i;
  vec v;

  v=o2v(f);


  return h;
}


OP vx(OP f,OP g){
  OP h={0},ww={0};
  OP v[K]={0},vv={0};
  oterm a,b;
  int i,j;

  v[0].t[0].a=0;
  v[0].t[1].n=0;
  v[1].t[0].a=1;
  v[1].t[1].n=0;

  //printf("in vx\n");
  //  exit(1);
  i=2;
       for(i=2;i<T+2;i++){
	 // memset(ss.t,0,DEG*sizeof(ss));
    //    if(deg(o2v(f))>=deg(o2v(g)) && deg(o2v(g))>0){
    if(deg(o2v(vv))<T){

	h=omod(f,g);
      oprintpol((h));
      printf(" modh vx==============\n");
      ww=odiv(f,g);
    }
    
  printf("ww======= ");
  oprintpol((ww));
  printf("\n");
  v[i]=oadd(v[i-2],omul(ww,v[i-1]));
  printf("-------");
  memset(f.t,0,sizeof(f.t));
  f=g;
  //memcpy(f.t,g.t,sizeof(g.t));
  memset(g.t,0,sizeof(g.t));
  g=h;
  
  if(deg(o2v(v[i]))<=T){
    vv=v[i];
    printf("vv==");
    oprintpol((vv));
    printf("\n");
    ss=h;
    // oprintpol((h));
    printf(" ll========\n");
    //    return vv;
    
  } else {
    
    printf("-------");
    break;
  }
       }

  
  return vv;
}



vec genrandompol(int n){
  vec x={0};
  int i,j=0,k;

  x=init_vec(x);
  for(i=0;i<K+1;i++)
    x.x[i]=xor128()%2;
  
  return x;
}


OP gcd(OP a,OP b){
  OP r={0},s={0};

  s=a;
  
  oprintpol((a));
  printf(" a==============\n");
  oprintpol((b));
  printf(" b==============\n");
  // exit(1);
  /* ユークリッドの互除法 */
  r = omod(a , b);
  if(deg(o2v(r))>0){
  oprintpol((r));
  printf(" rrr==============\n");
  }
  if(odeg(r)==0)
    return r;
  // exit(1);
  while(deg(o2v(r))>=0){
    oprintpol((r));
    printf(" r2==============\n");

    a = b;
    b = r;
     if(deg(o2v(b))>0)
      r = omod(a ,b);
      //  } 
    oprintpol((a));
    printf(" aa==============\n");
    if(deg(o2v(r))==0 && r.t[0].a==1)
      return r;
    if(deg(o2v(r))==0 && r.t[0].a==0)
      return b;
  }
 }




unsigned int ipow(unsigned int q,unsigned int u){
unsigned int i,m=1;

  for(i=0;i<u;i++)
    m*=q;

  printf("in ipow====%d\n",m);

  return m;
}





OP benor(int ww,int nn){
  OP w={0},ff={0},f={0},tt={0};
  int flg=0,i,j,k;
  vec v={0};

  while(1){
    for(i=0;i<DEG;i++)
      v.x[i]=0;
    
    v.x[0]=1;
    v.x[nn]=1;
    //    v.x[6]=1;
    k=0;
    
    while(k<ww){
      j=xor128()%10;
      if(j!=0 && j!=10 && v.x[j]==0){
	v.x[j]=1;
	k++;
      }
    }
    
  printpol(v);
  printf(" v=============\n");
  tt=v2o(v);
  //  exit(1);
  j=0;
  //    tt=irr(1,10);
  for(i=1;i<nn/2+1;i++){
    f.t[0].a=1;
    f.t[0].n=1;
    f.t[1].a=1;
    f.t[1].n=ipow(2,i);
    
    ff=gcd(tt,f);
    if(deg(o2v(ff))>0){
      oprintpol((ff));
      printf(" lcm==============%d\n",i);
      break;
    }
    if(deg(o2v(ff))==0)
      j++;
    oprintpol((ff));
    printf(" ff=============\n");
    printf("j=%d\n",j);
    printpol(v);
    printf(" vv===================\n");
  }
  if(j==nn/2){
    oprintpol((tt));
    printf(" irr?============\n");
    return tt;
  }
  j=0;
  //exit(1);
    }

  
}



unsigned char chk(OP f){
  int i,j,flg=0;
  vec x={0};

  x=o2v(f);
  for(i=0;i<DEG;i++){
    if(x.x[i]>0){
      flg=1;
      return 1;
    }
  }
  if(flg==0)
    return 0;
  
}


OP ogcd(OP f,OP g){
  OP h={0},ww={0};
  oterm a,b;
  int i=0;


  //oprintpol((f));
  //oprintpol((g));
  //  exit(1);
  
  while(i<T && deg(o2v(f))>=deg(o2v(g))){
    memset(ss.t,0,DEG);
    h=omod(f,g);
    ww=odiv(f,g);
  f=g;
  g=h;
  i++;
   if(deg(o2v(g))<=T-1)
     break;
  }
  // exit(1);

    
  return h;
}


EX xgcd(OP f,OP g){
  OP h={0},ww={0},v[K*2]={0},u[K*2]={0};
  oterm a,b;
  int i=0,j,k;
  EX e={0};

  
  u[0].t[0].a=1;
  u[0].t[0].n=0;
  u[1].t[0].a=0;
  u[1].t[0].n=0;
  u[2].t[0].a=1;
  u[2].t[0].n=0;

  v[0].t[0].a=0;
  v[0].t[0].n=0;
  v[1].t[0].a=1;
  v[1].t[0].n=0;


  //oprintpol((f));
  //oprintpol((g));
  //  exit(1);

  i=1;
  k=0;
  //while(deg(o2v(f))>0){
  for(i=1;i<T;i++){
    if(deg(o2v(g))==0)
      break;
    h=omod(f,g);
    ww=odiv(f,g);

    v[i+1]=oadd(v[i-1],omul(ww,v[i]));
    u[i+1]=oadd(u[i-1],omul(ww,u[i]));
    printf("i+1=%d\n",i+1);
    f=g;
    g=h;
    //i++;
    //    if(deg(o2v(v[i]))==4)
    //break;
    }

  v[i]=odiv(v[i],h);
  u[i]=odiv(u[i],h);
  h.t[0].a=1;
  h.t[0].n=0;
  printf("i=%d\n",i);
  printpol(o2v(v[i]));
  printf(" v=============\n");
  printpol(o2v(u[i]));
  printf(" u=============\n");
  printpol(o2v(h));
  printf(" h=============\n");
  //exit(1);

  e.d=h;
  e.v=v[i];
  e.u=u[i];
    
   return e;
}


OP bibun(vec a){
  OP w[T]={0};
  OP l={0},t={0};
 int i,j,k,n;
 vec tmp={0};
 
 
 n=deg(a);
 printf("n=%d\n",n);
 if(n==0){
   printf("baka8\n");
   //  exit(1);
 }
 
 for(i=0;i<T;i++){
   w[i].t[0].a=a.x[i];
   w[i].t[0].n=0;
   w[i].t[1].a=1;
   w[i].t[1].n=1;
   //oprintpol((w[i]));
 }
 //  exit(1);

 for(i=0;i<T;i++){
   tmp.x[0]=1;
   t=v2o(tmp);
   for(j=0;j<T;j++){
     if(i!=j)
       t=omul(t,w[j]);
   }
   //oprintpol((t));
   if(deg(o2v(t))==0){
     printf("baka9\n");
     // exit(1);
   }
   l=oadd(l,t);
 }
 
 
  return l;
}


vec chen(OP f){
  vec e={0};
  int i,count=0,n,x=0;
  unsigned short y[256]={0},z;
  oterm d={0};
  OP g={0};

  //  exit(1);
//  e=o2v(f);
n=deg(o2v(f));
//exit(1);
  for(x=0;x<D;x++){
    z=0;
    for(i=0;i<n+1;i++){
       if(f.t[i].a>0)
      z^=gf[mlt(mltn(f.t[i].n,fg[x]),fg[f.t[i].a])];
    }
    if(z==0){
      e.x[count]=x;
      count++;
    }
    //    printf("%d\n",x);
  }
  //printpol(e);
  // exit(1);
 

return e;
}

OP decode(OP f,OP s){
int i,j,k;
 OP r={0},w={0},e={0},l={0};
 oterm t1,t2,d1,a,b;
 vec x={0};
 unsigned short d=0;
 OP h={0};
 EX hh={0};

 
 printf("in decode\n");
 oprintpol((s));
 printf("\nsyn===========\n");
 r=vx(f,s);
 //hh=xgcd(f,s);
 
 if(deg(o2v(r))==0){
   printf("baka12\n");
    exit(1);
 }
 k=0;
 // exit(1);
 x=chen(r);
 // exit(1);

 /*
 printpol(x);
 printf("\n");
 // exit(1);
 t1=LT(r);
 printf("a=%d,n=%d\n",t1.a,t1.n);
 // exit(1);
 t2.a=t1.a;
 t2.n=0;
 oprintpol((r));
 printf("\n");
 w=oterml(v2o(x),t2);
 oprintpol((w));
 printf("\n");
 //  exit(1);
 printf("%d %d\n",deg(o2v(r)),deg(o2v(w)));
 // exit(1);
 r=ogcd(r,w);
 printf("r=");
 oprintpol((r));
 printf("\n");
 exit(1);
 */

 
 for(i=0;i<T;i++){
  printf("x[%d]=1\n",x.x[i]);
  if(x.x[i]==0)
    k++;
  if(k>1){
    printf("baka0\n");
    exit(1);
  }
 }
 //   exit(1);
 
 //  printf("\n");
 
  printf("あっ、でる！\n");  
  //  exit(1);
 
  if(deg(o2v(r))<T){
    oprintpol((r));
   printf("baka5 deg(r)<T\n");
   exit(1);
 }
 
 
 // printpol(x);
 //   exit(1);

 w=bibun(x);
 //  w=oterml(w,d1);
 //oprintpol((w));
printf("@@@@@@@@@\n");
//exit(1);

 h=ogcd(f,s);
 //oprintpol((h));
 // exit(1);
t1=LT(r);
// printf("t1=%d %d\n",t1.a,t1.n);
// exit(1);
t2.a=t1.a;
t2.n=0;
/*
 if(deg(o2v(w))==0){
   printf("baka3");
   exit(1);
 }
*/
 if(deg(o2v(w))==0){
   oprintpol((w));
 }
l=oterml(w,t2);
/*
 if(deg(o2v(l))==0){
   printf("baka5\n");
   exit(1);
 }
*/
// printf("\n");
 
 for(i=0;i<T;i++){
   //printf("%d ",trace(l,x.x[i]));
   if(trace(l,x.x[i])==0 || trace(h,x.x[i])==0){
     oprintpol((l));
     oprintpol((h));
     if(trace(l,x.x[i])==0 && trace(h,x.x[i])>0){
     printf("l=%d is 0 deg(x)=%d\n",i,deg(x));
       }else if(trace(l,x.x[i])>0 && trace(h,x.x[i])==0){
	 printf("h=0\n");
     }else if(trace(l,x.x[i])==0 && trace(h,x.x[i])==0){
       printf("l=0 & h=0\n");
     }
     if(x.x[i]==1)
       printf("x=1\n");
     exit(1);
   }
 }
 // printf("\n\n");
 for(i=0;i<T;i++){
   //printf("%d ",trace(h,x.x[i]));
   if(trace(h,x.x[i])==0){
     oprintpol((l));
     oprintpol((h));
     if(deg(x)<T-1){
	 printf("h=%d %d is 0 too few error deg(x)=%d\n",i,x.x[i],deg(x));
	 exit(1);
     }
   }
 }

 printf("%d\n",deg(x)+1);
 
 //    exit(1);
for(i=0;i<deg(x)+1;i++){
  //  if(x.x[i]>0){
  e.t[i].a=gf[mlt(fg[trace(h,x.x[i])],oinv(trace(l,x.x[i])))];
e.t[i].n=x.x[i];

// }
 }

   for(i=0;i<T;i++)
    if(gf[trace(h,x.x[i])]==0)
      printf("h=0");
   //printf("\n");
  for(i=0;i<T;i++)
    if(gf[oinv(trace(l,x.x[i]))]==0)
      printf("l=0\n");
  //  printf("\n");

  // exit(1);

 /*
   if(x.x[i]==0){
     printf("baka!");
     exit(1);
   }
 */

  // printf("\n");

 
return e;
}




OP setpol(unsigned short *f,int n){
OP g;
vec a;
 int i;

 
 memset(c,0,sizeof(c));
   // for(i=0;i<K+1;i++)
   // c[i]=0;
 memcpy(c,f,sizeof(c));
 //  for(i=0;i<n;i++)
 //  c[i]=f[i];
  a=Setvec(n);
  g=v2o(a);


return g;
}




void det(unsigned short g[]){
  OP f,h={0},w,u;
  unsigned short cc[K+1]={0},d[2]={0};
  unsigned short **HH;
  int i,j,a,b;
  oterm t={0};
  vec e;

  HH=malloc(D*sizeof(unsigned short *));
  for(i=0;i<2*K+1;i++) {
	HH[i] = malloc(sizeof(unsigned short) * D);
}
  //    memcpy(cc,g,sizeof(g));
  for(i=0;i<K+1;i++){
      cc[i]=g[i];
    printf("%d,",g[i]);
  }
  printf("\n");
  //  exit(1);
    //    cc[i]=g[i];
  k=cc[K];
  w=setpol(g,K+1);

  //#pragma omp parallel for       
  for(i=0;i<D;i++){
  
  a=trace(w,i);
  cc[K]=k;
  
  cc[K]^=a;
  f=setpol(cc,K+1);
  
  h.t[0].a=i;
  h.t[0].n=0;
  h.t[1].a=1;
  h.t[1].n=1;

  OP ww={0};

  memset(ss.t,0,sizeof(ss.t));
  ww=odiv(f,h);

  b=oinv(a);
  t.a=gf[b];
  t.n=0;

  u=oterml(ww,t);
  e=o2v(u);
  
  #pragma omp parallel for 
  for(j=0;j<K;j++)
    HH[j][i]=e.x[K-1-j];  

  }
  
  
  //#pragma omp parallel for    
  for(i=0;i<K;i++){
    //#pragma omp parallel for
    for(j=0;j<D;j++){
    m2[i][j]=mat[i][j]=HH[i][j];
    printf("%d,",mat[i][j]);
    }
    printf("\n");
  }  
  // exit(1);
}


void det2(unsigned short g[]){
  OP f,h={0},w,u;
  unsigned short cc[K+1]={0},d[2]={0};
  unsigned short **HH;
  int i,j,a,b,aa=0;
  oterm t={0};
  vec e[2*K]={0};

  
  HH=malloc(D*sizeof(unsigned short *));
  for(i=0;i<2*K+1;i++) {
	HH[i] = malloc(sizeof(unsigned short) * D);
}
  //    memcpy(cc,g,sizeof(g));
  for(i=0;i<K/2+1;i++){
      cc[i]=g[i];
    printf("%d,",g[i]);
  }
  printf("\n");
  //exit(1);
    //    cc[i]=g[i];
  k=cc[K];
  w=setpol(g,K/2+1);
  oprintpol((w));
  printf("\nw=================\n");
  //  exit(1);
  //#pragma omp parallel for       
  for(i=0;i<D;i++){
  
  a=trace(w,gf[i]);
  //  cc[K]=k;
  
  //cc[K]^=a;
  //f=setpol(cc,K+1);
  
  h.t[0].a=gf[i];
  h.t[0].n=0;
  h.t[1].a=1;
  h.t[1].n=1;

  OP ww={0};

  //  memset(ss.t,0,sizeof(ss.t));
  // ww=odiv(f,h);

  b=oinv(a);
  // b=b-1;
  printf("g^-%d=%d\n",i,mlt(b,b));
  //for(j=0;j<K;j++)
  //  for(j=0;j<K;j++){

  /*
  e[0].x[i]=b;
  e[1].x[i]=mlt(b,fg[gf[i]^f.t[1].a]);
  */

  e[0].x[i]=mlt(b,b);
  aa=gf[i]; //gf[mltn(j,i)]^gf[mlt(j-1,mlt(7,7))];
  e[1].x[i]=mlt(mlt(b,b),fg[aa]);
  e[2].x[i]=mlt(mlt(b,b),fg[gf[mltn(2,fg[aa])]^gf[15]]);
   e[3].x[i]=mlt(mlt(b,b),fg[gf[mltn(3,fg[aa])]^gf[mlt(15,fg[aa])]]);
  /*
    aa=0;
  for(j=1;j<2*K+1;j++){
    aa^=gf[mltn(j,i)]^gf[mlt(j-1,mlt(7,7))];
    e[j].x[i]=mlt(mlt(b,b),fg[aa]);
    }
  t.n=0;
  */
  //  exit(1);
  //u=oterml(ww,t);
  //e=o2v(u);
  
  #pragma omp parallel for 
  for(j=0;j<2*K;j++)
    HH[j][i]=e[j].x[i]; //e.x[K-1-j];  

  }
  
  
  //#pragma omp parallel for    
  for(i=0;i<K;i++){
    //#pragma omp parallel for
    for(j=0;j<D;j++){
    m2[i][j]=mat[i][j]=HH[i][j];
    printf("%d,",mat[i][j]);
    }
    printf("\n");
  }  
  //  exit(1);
}



void bdet(){
  int i,j,k,l;
  unsigned char dd[E*K]={0};
  FILE *ff;
  

  ff=fopen("Hb.key","wb");
  

  for(i=0;i<D;i++){
    for(j=0;j<K;j++){
      l=mat[j][i];
      //#pragma omp parallel for 
      for(k=0;k<E;k++){
	BH[j*E+k][i]=l%2;
	l=(l>>1);
      }
    }
  }
  for(i=0;i<D;i++){
    //#pragma omp parallel for 
    for(j=0;j<E*K;j++){
      //  printf("%d,",BH[j][i]);
      dd[j]=BH[j][i];
    }
    fwrite(dd,1,E*K,ff);
    //printf("\n");
  }
  fclose(ff);
}
 

void pubkeygen(){
  int i,j,k,l;
  FILE *fp;
  unsigned char dd[E*K]={0};

  
  fp=fopen("pub.key","wb");
  for(i=0;i<E*K;i++){
    for(j=0;j<D;j++){
      //#pragma omp parallel for 
      for(k=0;k<E*K;k++){
	tmp[i][j]^=cl[i][k]&BH[k][j];
      }
    }
  }
  P2Mat(P);

  for(i=0;i<E*K;i++){
    //  for(j=0;j<D;j++){
    //#pragma omp parallel for 
      for(k=0;k<D;k++)
	pub[i][k]=tmp[i][P[k]];//&A[k][j];
      //    }
  }
  for(i=0;i<D;i++){
    //#pragma omp parallel for 
    for(j=0;j<E*K;j++){
     dd[j]=pub[j][i];
    }
    fwrite(dd,1,E*K,fp);
  }
  fclose(fp);  

}


void Pgen(){
  unsigned int i,j;
  FILE *fp;
  

  fp=fopen("P.key","wb");
  random_permutation(P);
  for(i=0;i<D;i++)
    inv_P[P[i]]=i;
  fwrite(P,2,D,fp);
  fclose(fp);
  
  for(i=0;i<D;i++)
    printf("%d,",inv_P[i]);
  printf("\n");
  

  fp=fopen("inv_P.key","wb");
  fwrite(inv_P,2,D,fp);
  fclose(fp);

}


void key2(unsigned short g[]){
  FILE *fp;
  unsigned short dd[K]={0};

  printf("鍵を生成中です。４分程かかります。\n");
  fp=fopen("H.key","wb");
  det(g);
  //  exit(1);
      for(i=0;i<D;i++){
	for(j=0;j<K;j++)
	  dd[j]=mat[j][i];
	fwrite(dd,2,K,fp);
	
      }
      fclose(fp);

}


void keygen(unsigned short *g){
  
  int i;
  FILE *fp;
  

  key2(g);
  printf("end of ky2\n");
  makeS();
  printf("end of S\n");
  bdet();
  printf("end of bdet\n");
  Pgen();
  printf("end of Pgen\n");
  pubkeygen();  


  
}


void encrypt(){


}


void decrypt(){



}



int isqrt(unsigned short u){
  int i,j,k;

  for(i=0;i<D;i++){
    if(gf[mlt(i,i)]==u)
      return i;
  }

}



OP osqrt(OP f,OP w){
  int i,j,k,jj;
  OP even={0},odd={0},h={0},r={0},ww={0};
  oterm o={0};
  vec v={0};

  j=0;
  jj=0;
  k=distance(f);
  for(i=0;i<k+1;i++){
    if(f.t[i].n%2==0){
      even.t[j].n=f.t[i].n/2;
      even.t[j++].a=gf[isqrt(f.t[i].a)];
      printf("a=%d %d\n",f.t[i].a,i);
    }
    if(f.t[i].n%2==1){
      odd.t[jj].n=(f.t[i].n-1)/2;
      odd.t[jj++].a=gf[isqrt(f.t[i].a)];
      printf(" odd %d\n",i);
    }
  }
  printpol(o2v(even));
  printf(" T0============\n");
  printpol(o2v(odd));
  printf(" T1============\n");
  // exit(1);

  
  k=distance(w);
  printf("%d\n",k);
  //exit(1);
  j=0;
  jj=0;
  for(i=0;i<k+1;i++){
    if(w.t[i].n%2==0){
      h.t[j].a=gf[isqrt(w.t[i].a)];
      h.t[j++].n=w.t[i].n/2;
      printf("wi==%d %d\n",(w.t[i].n/2),i);
    }
    if(w.t[i].n%2==1){
      r.t[jj].a=gf[isqrt(w.t[i].a)];
      r.t[jj++].n=(w.t[i].n-1)/2;

    }
  }
  printpol(o2v(r));
  printf(" sqrt(g0)=======\n");
  
  //  exit(1);
  r=inv(r,w);
  //exit(1);
  ww=omul(h,r);
  //h=om(ww,ww);
    ww=omod(ww,w);
    printpol(o2v(ww));
    printf(" w==============\n");
    //    exit(1);
  printpol(o2v(r));
  printf(" g1^-1=========\n");
  printpol(o2v(h));
  printf(" g0=========\n");
  //exit(1);
  printpol(o2v(ww));
  printf(" ww==========\n");
  //  exit(1);
  h=ww;
  if(odeg(omod(omul(h,ww),w))==1){
    ww=h;
    h=omod(oadd(even,omul(ww,odd)),w);
    return h;
  }else{
    printf("vaka\n");
    exit(1);
  }
  
  printpol(o2v(ww));
  printf(" w==========\n");

}


int main(int argc,char **argv){
  int i,j,k,l,c;
  unsigned long a,x,count=1;
  //  unsigned short cc[K]={0};
  unsigned short m[K],mm[T]={0},dd[K*D]={0};
  time_t timer;
  FILE *fp,*fq;

  unsigned short g2[7]={1,0,9,0,0,6,4};
  //  unsigned short s[K]={0}; //{4,12,7,8,11,13};
  unsigned short jj[T]={0};
  unsigned short ee[10]={1,2,3,4,5,6,7,8,9,10};
  unsigned short zz[D]={0};
  //  unsigned short zz[T]={10,97,114,105,97,98,108,101,32,80,111,108,121,110,111,109};
  int y,flg,o1=0;
  OP f={0},h={0},r={0},w={0},aa[8]={0},tt={0},ff={0};
  EX hh={0};
  vec v;
  unsigned short d=0;
  time_t t;
  unsigned short gg[K+1]={0};
  oterm rr={0};
  OP r1={0},r2={0},t1={0},t2={0},a1={0},b1={0},a2={0},b2={0};
  
  //  unsigned short syn[K]={4,12,7,8,11,13};
  //unsigned short g[K+1]={1,0,0,0,1,0,1};
  //unsigned short g[K/2+1]={1,0,1,1};
    //  unsigned short syn[K]={4,12,7,8,11,13};
  //unsigned short g[K+1]={1,0,0,0,1,0,1};
  //unsigned short g[K+1]={1,3,9,4,5,2,3,5,10};
  //unsigned short g[K+1]={2,14,5,4,4,0,15,12,10,1,14,7,1};
  //unsigned short g[K+1]={1,1,0,2,5,3,1,0,1};
  //unsigned short g[K+1]={1,1,2};
  //unsigned short g[K+1]={1,0,9,0,1};
  //  unsigned short g[K+1]={1,0,1};
  //unsigned short g[K+1]={1,1,1};
  //unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};
  unsigned short yy[5]={15,0,8,0,11};
  /*
  r1=setpol(yy,5);
  w=setpol(g,7);
  r2=osqrt(r1,w);
  printpol(o2v(r2));
  printf(" osqrt==========\n");
  */
  //exit(1);
  //  makegf(M);
  //  makefg(M);
  srand(clock()+time(&t));
  printf("@");
  ginit();
  /*  
  ff.t[0].a=1;
  ff.t[0].n=0;
  ff.t[1].a=1;
  ff.t[1].n=4;
  ff.t[2].a=1;
  ff.t[2].n=10;
  f.t[0].a=1;
  f.t[0].n=1;
  f.t[1].a=1;
  f.t[1].n=32;
  
    h=gcd(f,ff);
  oprintpol(h);
  printf("gcd=========\n");
  while(odeg(ff)<10){
  ff=benor(3,10);
  }
  oprintpol((ff));
  
  printf(" irr?=============\n");
    exit(1);
  */

  
    w=setpol(g,K+1);
    oprintpol(w);

    // exit(1);

  /*  
  //-------------２乗するとき外す
   w=setpol(g,K/2+1);
  oprintpol((w));
  printf("\n");
  //   exit(1);
  w=omul(w,w);
  oprintpol((w));
  printf("\n");
  //exit(1);
    
  v=o2v(w);
  for(i=0;i<K+1;i++){
    printf("%d,",v.x[K-i]);
    gg[K-i]=v.x[i];
  }
  printf("\n");
  //w=setpol(gg,K+1);
  //oprintpol((w));
  //printf("\n");
  //   exit(1);
  //--------------
  */
  
#pragma omp parallel for  
  for(i=0;i<D;i++){
    a=trace(w,i);
    if(a==0){
      printf("trace 0 @ %d\n",i);
         exit(1);
    }
  }
  printf("@");
  //keygen(g);
  key2(g);
  //exit(1);
  /*
  for(i=0;i<K;i++){
    for(j=0;j<M;j++){
      printf("%d,",mat[i][j]);
    }
    printf("\n");
  }
  //exit(1);
  */
  

  /*
  do{
  flg=0;
  v=genrandompol(K+1);
  for(i=0;i<K+1;i++)
    g[i]=(unsigned short)v.x[i];
    for(i=0;i<K+1;i++){
      w=setpol(g,K+1);
      a=trace(w,i);
    if(a==0)
      flg=1;
      }
  } while(deg(v)<K || flg==1);
  */

  //  w=setpol(g,K+1);
  printf("すげ、オレもうイキそ・・・\n");

  uu=0;
  #pragma omp parallel for
  for(i=0;i<D;i++){
  a=trace(w,i);
  if(a==0){
    printf("trace 0\n");
       exit(1);
  }
  }
  //oprintpol((w));
  //    exit(1);

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

  /*
  j=0;
  while(j<T){
    flg=0;
    l=xor128()%6688;
    printf("%d,",l);
       for(k=0;k<T;k++){
	 if(l==jj[k])
	   flg=1;
       }
       if(flg==0){
	 jj[j]=l;
	 zz[jj[j]]=1;
	 j++;
       }
  }
  */
  
  j=0;
  while(j<T*2){
    l=xor128()%D;
    printf("l=%d\n",l);
    if(0==zz[l]){
      zz[l]=1;
      j++;
    }
  }
  
   
  
  
  //  for(i=0;i<T;i++)
  //zz[i]=1;
  
  // det(g);

  //  exit(1);
  //x^4 + a^3*x^3 + (a + 1)*x^2 + (a + 1)*x + a^3 + a
  //a^3*x^3 + (a^3 + a)*x^2 + 1


  fq=fopen("H.key","rb");
  
  fread(dd,2,K*D,fq);
  #pragma omp parallel for
  for(i=0;i<D;i++){
  for(j=0;j<K;j++)
    mat[j][i]=dd[K*i+j];
    }
  
  
  for(j=0;j<D;j++){
    flg=0;
    for(i=0;i<K;i++){
      //printf("%d,",mat[i][j]);
      if(mat[i][j]>0)
	flg=1;
      //      printf("\n");
    }
    if(flg==0)
      printf("0 is %d\n",j);
  }
  
  //  exit(1);   

  
  printf("zz=");
  for(i=0;i<D;i++)
    printf("%d,",zz[i]);
  printf("\n");
  //    exit(1);
  //  
  for(i=0;i<K;i++){
    syn[i]=0;
    //#pragma omp parallel for
    for(j=0;j<D;j++){
      //   printf("%u,",zz[jj[j]]);
      syn[i]^=gf[mlt(fg[zz[j]],fg[mat[i][j]])];
    }
       printf("syn%d,",syn[i]);
  }
  printf("\n");
  //    exit(1);  
  for(i=0;i<K;i++)
    printf("mat[%d][1]=%d\n",i,mat[i][1]);
  printf("\n");
  //    exit(1);
  
  OP g1={0},ll={0},s={0};

  tt.t[0].n=1;
  tt.t[0].a=1;

  f=setpol(syn,K);
  printpol(o2v(f));
  printf(" syn=============\n");
  //   exit(1);
  
  ff=inv(f,w);
  printpol(o2v(ff));
  printf("locater==========\n");
  //exit(1);
  r2=oadd(ff,tt);
  printpol(o2v(r2));
  printf(" h+x==============\n");
  //  exit(1);
  g1=osqrt(r2,w);
    printpol(o2v(g1));
  printf(" g1!=========\n");
  //exit(1);
  hh=xgcd(w,g1);
  ff=omod(omul(hh.v,g1),w);
  printpol(o2v(ff));
  printf(" beta!=========\n");
  printpol(o2v(hh.v));
  printf(" alpha!=========\n");
  //exit(1);
  ll=oadd(omul(ff,ff),omul(tt,omul(hh.v,hh.v)));
  v=chen(ll);
  for(i=0;i<T*2;i++)
    printf("%d\n",v.x[i]);
   exit(1);
  /*  
  // r=vx(w,f);
  //s=omod(omul(r,g1),w);
  
  hh=xgcd(w,g1);
  printpol(o2v(r));
  printf(" vx===========\n");
  printpol(o2v(hh.v));
  printf(" hh.v==============\n");
  printpol(o2v(hh.u));
  printf(" hh.v*g1==============\n");
  
  ll=r;
  v=chen(ll);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  exit(1);
  
  s=inv(f,w);
  s=oadd(s,tt);
  printpol(o2v(s));
  printf(" s^-1+x===========\n");
  // exit(1);
  
  //r=vx(w,f);
   
  //g1=osqrt(s,w);
  //r=vx(w,g1);
  //hh=xgcd(g1,w);
  */
  /*
  s=vx(g1,hh.u);
  r=vx(g1,s);
  f=omod(hh.v,w);
  printpol(o2v(s));
  printf(" a?===============\n");
  printpol(o2v(r));
  printf(" b?============\n");
  printpol(o2v(f));
  printf(" v?============\n");
  printpol(o2v(hh.u));
  printf(" u?============\n");
  //exit(1);
  //ff=vx(w,f);
   //  printf(" beta^2============\n");
  ll=oadd(omul(tt,omul(s,s)),omul(r,r));
  v=chen(ll);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  exit(1);
  
  g1=osqrt(s,w);
  printpol(o2v(g1));
  printf(" s^-1+x==========\n");
  //exit(1);
  r=vx(w,g1);
  r1=omod(omul(g1,r),w);
  oprintpol(r);
  printf(" b?=========\n");
  printpol(o2v(r1));
  printf(" a?==========\n");
  //r=omod(omul(r1,g1),w);
  //exit(1);
  
  hh=xgcd(g1,w);
  //s=omod(omul(r1,g1),w);
  ll=omod(oadd(omul(hh.u,hh.u),omul(tt,omul(hh.v,hh.v))),w);
  
  printpol(o2v(ll));
  printf(" locater?============\n");
  v=chen(ll);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  exit(1);
  */
  /*
  printpol(o2v(hh.u));
  printf(" u=============\n");
  printpol(o2v(hh.v));
  printf(" v=============\n");
  printpol(o2v(hh.d));
  printf(" d=============\n");
  ll=oadd(omul(hh.v,hh.v),omul(omul(hh.u,hh.u),tt));
  */
  /*
  v=chen(r);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  //exit(1);
  
  if(deg(o2v(f))==0)
    exit(1);
  
  ff=inv(f,w);
  ff=oadd(tt,ff);
  oprintpol(ff);
  printf(" inv(s)+x\n");
  //exit(1);
  r=vx(w,f);
  oprintpol(r);
  printf(" locater========\n");
  h=omod(omul(f,r),w);
  r=omod(omul(ff,h),w);
  printpol(o2v(r));
  printf(" b^2============\n");
  //exit(1);
  //ff=oadd(r,omul(h,tt));
    //ff=omod(omul(ff,h),w);
  //  ff=(omul(ff,h),w);
  //printpol(o2v(ff));
  //printf(" a^2=============\n");
  //r=oadd(ff,omul(h,tt));
  //printpol(o2v(r));
  //printf(" locate2============\n");
  //  exit(1);
  //exit(1);
  h=vx(w,g1);
  f=omod(omul(h,g1),w);
  h=omod(omul(f,g1),w);
  ll=omod(oadd(omul(h,h),omul(omul(f,f),tt)),w);
  
  printpol(o2v(ll));
  printf(" locater?==========\n");
  //exit(1);
  v=chen(ll);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  exit(1);
  */
  /*
  printpol(o2v(ff));
  printf(" b?==========\n");
  //  exit(1);
  printf(" again a=======\n");
  ff=omod(omul(f,g1),w);
  printpol(o2v(ff));
  printf(" again b==========\n");
  exit(1);
  ff=omod(omul(ff,g1),w);
  ff=coeff(ff,9);
  printpol(o2v(ff));
  printf(" regain a==========\n");
  printpol(o2v(f));
  printf(" again a=======\n");
  exit(1);
  ll=oadd(ff,h);
  //  exit(1);
    printpol(o2v(ll));
  printf(" locater?=============\n");
  //  exit(1);
  v=chen(ll);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  exit(1);
  */
  /*
  r1.t[0].a=9;
  r1.t[0].n=0;
  r1.t[1].a=14;
  r1.t[1].n=1;
  r1.t[2].a=8;
  r1.t[2].n=2;
  r1.t[3].a=15;
  r1.t[3].n=3;
  r1.t[4].a=4;
  r1.t[4].n=4;
  r1.t[5].a=9;
  r1.t[5].n=5;
  r2.t[0].a=14;
  r2.t[0].n=0;
  r2.t[1].a=2;
  r2.t[1].n=1;
  r2.t[2].a=9;
  r2.t[2].n=2;
  r2.t[3].a=1;
  r2.t[3].n=1;
  r2.t[4].a=8;
  r2.t[4].n=4;
  r2.t[5].a=14;
  r2.t[5].n=5;
  b1=vx(w,r1);
  oprintpol(b1);
  printf(" b1=======\n");
  //  oprintpol(t1);
  //printf(" test=====\n");
  a1=omod(w,omul(r1,b1));
  oprintpol(a1);
  printf(" a1==========\n");
  t1=omod(oadd(omul(a1,a1),omul(omul(b1,b1),tt)),w);
  oprintpol(t1);
  printf(" rho1=======\n");
  b2=vx(w,r2);
  oprintpol(b2);
  printf(" b2=========\n");
  a2=omod(omul(r2,b2),w);
  oprintpol(a2);
  printf("a2============\n");
  t2=omod(w,oadd(omul(a2,a2),omul(omul(b2,b2),tt)));
  oprintpol(t2);
  printf(" rho2=======\n");
  v=chen(t1);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  //  exit(1);
  */
  //r=vx(w,f);
  /*
  r=vx(w,f);
    oprintpol(r);
  printf(" real=========\n");
  exit(1);
  h=omod(omul(r,f),w);
  printpol(o2v(h));
  printf(" const?=====\n");
  h=coeff(omul(r,f),h.t[0].a);
  printpol(o2v(h));
  
  ff=vx(w,r);
  printpol(o2v(ff));
  printf(" b^2=======\n");
  exit(1);
  v=chen(r);
  for(i=0;i<D;i++)
    printf("%d\n",v.x[i]);
  exit(1);
  */

  r=decode(w,f);
  
  
  for(i=0;i<T;i++){
    mm[i]=r.t[i].a;
    if(i==0){
     printf("e=%d %d %s\n",r.t[i].a,r.t[i].n,"う");
    }else{     printf("e=%d %d %s\n",r.t[i].a,r.t[i].n,"お");
}
    if(r.t[i].a==0){
      printf("------------------\n");
      
       printf("err=%d i=%d\n",o1,i);
        exit(1);
    }
  }
  printf("         っ！！\n");
    for(i=0;i<D;i++){
    if(zz[i]>0)
      o1++;
  }
  printf("err=%d\n",o1);


  return 0;
}


  //exit(1);


  //  return h;

     /*
  //実験中
  ff=inv(f,w);
  tt.t[1].n=1;
  tt.t[1].a=1;
  h=oadd(ff,tt);
  oprintpol((h));
  printf(" d^2==================\n");
  //    exit(1);
  //ff=setpol(xx,K);
  //oprintpol((ff));
  printf(" sqrt===============\n");
  r=vx(w,ff);
  h=omul(omul(r,r),h);
  f=omod(h,w);
  ff=omul(ss,ss);
  //h=oadd(r,ss);
  //  f=omod(r,w);
  oprintpol((f));
  printf(" B^2*(h+x) mod g==============\n");
  oprintpol((ff));
  printf(" alpha^2==============\n");
  exit(1);
  oprintpol((r));
  printf(" beta==============\n");
  oprintpol((ss));
  printf(" alpha==============\n");
  r=omul(omul(r,r),tt);
  oprintpol((r));
  printf(" x*b^2================\n");
  oprintpol((ff));
  printf(" aptha^2================\n");
  //  ff=oadd(ff,r);
  // oprintpol((ff));
  //printf(" rho================\n");
  #pragma omp parallel for
  for(i=0;i<D;i++)
    printf("x=%d %d\n",trace(r,i),i);
  exit(1);
  //ff=omod(ff,w);
  v=chen(ff);
  for(i=0;i<16;i++)
    printf("x[%d]=1\n",v.x[i]);
  exit(1);
  */  
    /*
  exit(1);
  oprintpol((ss));
  printf(" b^2==============\n");
  exit(1);
  
  h=inv(f,w);
  tt.t[1].n=1;
  tt.t[1].a=1;
  h=oadd(tt,f);
  oprintpol((h));
  printf(" d^2==================\n");
  //  exit(1);
  ff=setpol(xx,K);
  r=vx(w,ff);
  //  exit(1);
  oprintpol((r));
  r=omul(r,r);
  printf(" b=================\n");
  //  exit(1);
  ff=omul(r,tt);
  ss=omul(ss,ss);
  h=oadd(ff,ss);
  oprintpol((h));
  printf(" rho1=================\n");
  f=omod(h,w);
  //  f=omod(f,w);
  oprintpol((f));
  printf(" rho2================\n");
  v=chen(f);
  for(i=0;i<16;i++)
    printf("x[%d]=1\n",v.x[i]);
     exit(1);
    r=inv(f,w);
    ff.t[1].n=1;
    ff.t[1].a=1;
    ff.t[0].a=0;
    ff.t[0].n=0;
    r=oadd(r,ff);
    //    tt=omul(r,f);
  oprintpol((r));
  printf(" d^2=============\n");
  oprintpol((w));
  printf(" g================\n");
  exit(1);
  h.t[1].n=1;
  h.t[1].a=7;
  ff=vx(w,r);
  oprintpol((ff));
  printf(" beta^2==============\n");
  exit(1);
  oprintpol((omod(omul(r,f),w)));
  printf("\nmodg===============\n");
  oprintpol((tt));
  printf(" tt===================\n");
  */
  //  rr=LT(r);
  /*
  printf("c=%d\n",rr.a);
  ff.t[0].a=rr.a;
  ff.t[0].n=0;
  r=odiv(r,ff);
  */
  /*
  r=inv(f,w);
  oprintpol((r));
  printf(" h==========\n");
  ff=omod(omul(r,f),w);
  oprintpol((ff));
  printf(" omod=============\n");
  //  exit(1);
  //oprintpol((r));
  //printf("s^-1=============\n");
  oprintpol((r));
  printf(" d^2============\n");
  //  exit(1);
  h=vx(w,r);
  oprintpol((h));
  printf(" beta^2=================\n");
  //  exit(1);
  f=omod(omul(h,r),w);
  oprintpol((f));
  printf(" alpha=============\n");
  h=oadd(omul(h,ff),f);
  oprintpol((h));
  printf(" rho============\n");
  //exit(1);
  v=chen(h);
  //  exit(1);
  // v=chen(h);
  k=0;
 for(i=0;i<6;i++){
  printf("x[%d]=1\n",v.x[i]);
  if(v.x[i]==0)
    k++;
  if(k>1){
    printf("baka0\n");
    //  exit(1);
  }
 }
 //  hh=bibun(v);
 // oprintpol((hh));
 //printf("w=================\n");
  exit(1);
  v=chen(r);
  for(i=0;i<D;i++){
    if(v.x[i]>0)
    printf("%d %d\n",v.x[i],i);
  }
  // exit(1);
  */
  /*
  ff=inv(f,w);
  oprintpol((ff));
  printf("\nff3================\n");
  // exit(1);
 oprintpol((omod(omul(f,ff),w)));
  printf("\n");
  //   exit(1);
  */
