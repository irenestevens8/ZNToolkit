/* program to initialize binary hs chromosome files to zeros */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int,char *[]);
int totestMode(void);

FILE *fp;
FILE *fTmReadp;
FILE *fTmWritep;

static char chrSet[40] = "";
static char chrBinFn[100] = "";
static char tmWFn[100] = "";
static char znt[10] = ".znt";
static char test[10] = ".test";

main (int argc,char *argv[])
{

	int nVal, testMode;
	long int n, chrLen;
	char chrSet[20];
	
	nVal = 0;
	
	testMode = 0;
	
	if(argc != 3)
		{
		printf("argc %d; init_hsGenome chrNo chrLen\n");
		exit(0);
		}

	strcpy(chrSet,argv[1]);
	chrLen = atol(argv[2]);

	strcpy(chrBinFn,chrSet);
	strcat(chrBinFn,znt);
	strcpy(tmWFn,chrSet);
	strcat(tmWFn,test);

	fp = fopen(chrBinFn,"wb");
	if (fp == NULL)
		{
		printf("\nOpen chrBinFn wb command failed.\n");
		exit(0);
		}
	
	for(n = 0; n < chrLen; n++)
		{
		putc(nVal,fp);
		}

	fclose(fp);
			
	if(testMode)
		{
		totestMode();
		}

	return 0;
}

totestMode(void)
{
	int linPtr;
	int bbase;
	
	fTmReadp = fopen(chrBinFn,"rb");
		if (fTmReadp == NULL)
		{
		printf("chrBinFn rb test command failed.\n");
		exit(0);
		}
		
	fTmWritep = fopen(tmWFn,"w");
		if (fTmWritep == NULL)
		{
		printf("fTmFn w test command failed.\n");
		exit(0);
		}

	linPtr = 0;
	
	fprintf(fTmWritep,">%s\n",chrBinFn);
	do
		{
		bbase = getc(fTmReadp);
		if (bbase == 0 || bbase == 1 || bbase == 2 || bbase == 3 || bbase == 4)
			{
			fprintf(fTmWritep,"%d",bbase);
			linPtr++;
			if (linPtr == 59)
				{
				fprintf(fTmWritep,"\n");
				linPtr = 0;
				}
			}
		}while(!feof(fTmReadp));

	fprintf(fTmWritep,"\n");
	
	fclose(fTmReadp);
	fclose(fTmWritep);

	return 0;
}
