#include <stdio.h>

int main()
{
    double s = 0;
    double tmp;
    FILE * fp = fopen("time.info", "r");
    while(fscanf(fp, "%lf", &tmp)!=EOF)
    {
        s += tmp;
    }
    printf("%g\n", s/10.0);
    fclose(fp);
    return 0;
}
