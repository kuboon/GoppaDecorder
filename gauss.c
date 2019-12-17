#include <stdio.h>
#include <stdlib.h>


unsigned char a[4][4]={{0,1,1,1},{1,0,1,1},{0,0,1,1},{1,0,0,1}}; //{{1,2,0,-1},{-1,1,2,0},{2,0,1,1},{1,-2,-1,1}}; //入力用の配列
unsigned char inv_a[4][4]; //ここに逆行列が入る
unsigned char buf; //一時的なデータを蓄える
int i,j,k; //カウンタ
int n=4;  //配列の次数

void main(){
  int i,j,k,l,ii;
  unsigned char b[4][4]={0};
  unsigned char c[4][4]={0};
  /*
  srand(clock());
  for(i=0;i<n;i++){
    for(j=0;j<n;j++)
      a[i][j]=rand()%2;
  }
  */
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      printf("%d,",a[i][j]);
      c[i][j]=a[i][j];
    }
    printf("\n");
  }

  
//単位行列を作る
for(i=0;i<n;i++){
 for(j=0;j<n;j++){
 inv_a[i][j]=(i==j)?1.0:0.0;
 }
}
 printf("\n");
 for(i=0;i<n;i++){
   for(j=0;j<n;j++)
     printf("%d,",a[i][j]);
   printf("\n");
 }
//掃き出し法
for(i=0;i<n;i++){
  if(a[i][i]==0){
  j=0;
  while(a[j][i]==0){
    buf=a[j++][i];
  }
  printf("j=%d\n",j);
  
  //  exit(1);
 for(k=0;k<n;k++){
 a[i][k]^=a[j][k];
 inv_a[i][k]^=inv_a[j][k];
 }
 for(j=0;j<n;j++){
   for(k=0;k<n;k++)
     printf("%d",a[j][k]);
   printf("\n");
 }
 for(j=0;j<n;j++){
   for(k=0;k<n;k++)
     printf("%d,",inv_a[j][k]);
   printf("\n");
 }
  }
  //  exit(1);
  
 if(a[i][i]==1){
 for(l=i+1;l<n;l++){
   if(a[l][i]==1){
     for(k=0;k<n;k++){
     a[l][k]^=a[i][k];
     inv_a[l][k]^=inv_a[i][k];
     }
   }
 }

 // printf("@%d\n",i);
 }
 printf("@i=%d\n",i);
 }
 for(i=1;i<n;i++){
   for(k=0;k<i;k++){
     if(a[k][i]==1){
       for(j=0;j<n;j++){
       // if(a[k][i]==1){
	 a[k][j]^=a[i][j];
	 inv_a[k][j]^=inv_a[i][j];
	 //}
     }
   }
   printf("i=%d\n",i);
   for(l=0;l<n;l++){
     for(ii=0;ii<n;ii++)
       printf("%d",a[l][ii]);
     printf("\n");
   }
   for(l=0;l<n;l++){
     for(ii=0;ii<n;ii++)
       printf("%d, ",inv_a[l][ii]);
     printf("\n");
   }
   
   }
 }
 for(i=0;i<n;i++){
   for(j=0;j<n;j++)
     printf("%d,",inv_a[i][j]);
   printf("\n");
 }
 //   exit(1);
 
//逆行列を出力
for(i=0;i<n;i++){
 for(j=0;j<n;j++){
  printf(" %d,",inv_a[i][j]);
 }
 printf("\n");
 }

// exit(1);
//検算
 for(i=0;i<n;i++){
   for(j=0;j<n;j++){
     for(k=0;k<n;k++){
       b[i][j]^=(c[i][k]&inv_a[k][j]);
     }
   }
 }
 for(i=0;i<n;i++){
   for(j=0;j<n;j++)
     printf("%d,",b[i][j]);
   printf("\n");
 }
 
}
/* 出力
2.000000 2.000000 -1.000000 3.000000
-4.000000 -5.000000 3.000000 -7.000000
3.000000 4.000000 -2.000000 5.000000
-7.000000 -8.000000 5.000000 -11.000000
*/
