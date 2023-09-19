#include <stdio.h>
#include <math.h>

void happy_meter(int size){
  char* smile = " :-)";
  printf("I am happy about this assignment");
  for(int i=0; i<size;i++){
    printf("%s",smile);
    
  }

  
}

int numfactors(long k){
  int count =0;
  for(int i=1;i<=k;i++){
    if(k%i==0){
      count++;
    }
    
  }
  return count;
}


int order_of_magnitude(double a){
  int c = floor(log10(a));
  return c;
}
int num_conseq_digits(long k){
  int exponent = floor(log10(k));
  int lc = k % 10;
  int count=0;
  int maxi = 0;
  /* printf("Exponent %d\n",exponent); */
  for(int i=2;i<=exponent+1;i++){
    int cc = ((k % (int) pow(10,i))-lc)/ (int)pow(10,i-1);
    /* printf("Last Character %d \t Current Character %d\n",lc,cc); */
    count++;
    if(cc-lc != 0){
      count=0;
    }
    /* printf("Count %d Maxi %d\n",count,maxi); */
    if(count > maxi) {
      maxi = count;
    }
    lc = cc;
    }
  return maxi+1;

  
}







int main(void) {
  printf("%d",num_conseq_digits(1443334));
  return 0;
}