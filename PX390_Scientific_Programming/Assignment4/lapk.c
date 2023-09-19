#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mkl_lapacke.h>


int parr(double* arr,int len){
  printf("START\n");
  for(int i =0; i<len;i++){
    printf("%lf\n",arr[i]);
  }
  printf("END\n");
  return 0;
}

int main(){
    double ar[] = {1.0,2.0};
    double mat[]= {9.2,1.0,12.0,4.0};
    int ipiv[2];
    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,2,1,mat,2,ipiv,ar,1);
    parr(ar,2);
    printf("%d\n",info);
    //printf("Hello");
    return 0;
}