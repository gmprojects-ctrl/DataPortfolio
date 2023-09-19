#include<stdio.h>
#include<math.h>

int main(){
    long a =0;
    long b =0;
    scanf("%ld %ld",&a,&b);
    printf("%ld",(long) fabs(a-b));
    return 0;
}