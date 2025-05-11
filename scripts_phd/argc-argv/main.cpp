#include <iostream>

void fun1(char***lol){}

int main(int argc, char **argv)
{

char *argv2[2];
char **argv3;
argv2[0]="lol";
argv2[1]="lol2";


char c1[10]={'lol'};

fun1(&argv3);

for (int i=0; i<argc; i++) std::cout<<argv[i]<<std::endl;
for (int i=0; i<argc; i++) std::cout<<argv2[i]<<std::endl;


return 0; 
}
