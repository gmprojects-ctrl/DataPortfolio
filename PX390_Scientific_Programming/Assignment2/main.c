#include <stdio.h>
#include <math.h>
#include <stdlib.h>
void printmem(int *arr,int arr_length){
  for(int i =0; i<arr_length;i++){
    printf("%d ",arr[i]);
  }
}

int  *move_to_front(int *arr, int arr_length){
  int* temp = (int*) malloc(sizeof(int)*arr_length);
  temp[0] = arr[arr_length-1];
  for(int i=1; i<arr_length;i++){
    temp[i] = arr[i-1];
  }
  return temp;
}

int  *arrdiff(int *arr, int arr_length){
  int* temp = (int*) malloc(sizeof(int)*arr_length-1);
  for(int i=0; i <arr_length-1;i++){
    temp[i] = arr[i+1] -arr[i];
  }
  return temp;
}

int  *reverse_order(int *arr, int arr_length){
    int* temp = (int*) malloc(sizeof(int)*arr_length);
  for(int i=0;i<arr_length;i++){
    temp[i] = arr[arr_length-(i+1)];
  }
  return temp;
}

void mean_std(double *arr, long arrlen, double *mean, double *std){
  double cumsum =0;
  double var=0;
  for(long i=0; i<arrlen;i++){
    cumsum+=arr[i];
  }
  double t_mean = cumsum/ ((double) arrlen);
  for(long i=0; i<arrlen;i++){
    var+=(arr[i]-t_mean)*(arr[i]-t_mean);
  }
  double t_std = sqrt(var/((double)arrlen));
  *mean = t_mean;
  *std  = t_std;
}



void en_pot(double *posx, double *posy, double *posz,  long ncharges, double *res){
  long curr_ch = 0;
  double sum =0;
  for(int i=0;i<ncharges;i++){
      for(int j=0;j<ncharges;j++){
        if(curr_ch==j){
          continue;
        }else{
          double x = posx[j] - posx[i];
          double y = posy[j] - posy[i];
          double z = posz[j] - posz[i];
          /* printf("%lf %lf %lf \n",x,y,z); */
          double rep = x*x + y*y + z*z;
          sum+= 1/sqrt(rep);
        }
    }
    curr_ch+=1;
    }
  sum = sum/2;
  *res = sum;
}
  
  




int main(void) {
  /* double res=0;
  double posx[]={1.,4.,3.,6.};
  double posy[]={3.,5.,1.,12.};
  double posz[]={7.,2.,2.,1.};
  en_pot(posx, posy,posz,4, &res);
  printf("%lF",res);*/
  //int arr[] ={3,4,1,12,-6};
  //int arr2[] = {3,4,-6};
  //int *mtf = move_to_front(arr, 5);
  //int *ard = arrdiff(arr2, 3);
  //int *rev = reverse_order(arr, 5);
  //printmem(rev,5); 
  double mean,std;
  double data[] = {1.,234.,23.,4.,5.};
  mean_std(data,5,&mean,&std);
  printf("%lf %lf",mean,std);
}