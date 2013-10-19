/* derivative of mark_blast designed to match with .znt files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int,char *[]);

FILE *fp;
FILE *fDatap;

static char chrDtFn[100];
static char chrZntFn[100];
static char chrSet[20] =;
static char Outp[20] = "btpOut.";
static char znt[5] = ".znt2";


main (int argc,char *argv[])
{

	int m;
	int n;
	int bbChk;
	int readLen;
	int errPrt;
	unsigned short int bbase;
	unsigned short int bbMax;
	long int cp, chrLen, posMax;
	char chrVal[20];
	char orient;
	char seqRead[120];
	
	cp = 0;
	chrLen = 0;
	/* readLen = 50; */
	/* readLen = 40;  for tufts_052110 set */
	
	
		if (argc < 2) {
		printf("Usage: mark_bt16 [options]\n\n");
		printf("where:\n");
		printf("-chrNo ZNTFILE \t (MANDATORY)\n"); 
		printf("-outfile STR \t the name of the output file (MANDATORY)\n\n");
		printf("Options:\n");
		printf("-readlen INT \t specifies the average size of the DNA reads\n");
		printf("-format STR \t format of the read files: sam, bowtiesam, fastq. Default is sam\n");
		printf("-genome STR \t genome version: hsG37.2, hg19, Ensembl. Default is hsG37.2\n");

		printf("\n");
		exit(0);
	}
	
	
	//initialize variables with the arguments values
	
	if (exist_parameter(argc, argv, "-chrDtFn"))
		-chrDtFn  = get_parameter(argc, argv, "-chrDtFn");

	if (exist_parameter(argc, argv, "-chrZntFn"))
		-chrZntFn  = get_parameter(argc, argv, "-chrZntFn");
		
	if (exist_parameter(argc, argv, "-chrSet"))
		-chrSet  = get_parameter(argc, argv, "-chrSet");

	
	if (exist_parameter(argc, argv, "-readlen"))
		-readlen  = get_parameter(argc, argv, "-readlen"));



	fp = fopen(chrDtFn,"rb+");
	if (fp == NULL)
		{
		printf("\nOpen chrDtFn rb+ command failed - %s\n",chrDtFn);
		exit(0);
		}
	fseek(fp,0,2);
	chrLen = ftell(fp);
	printf("chrLen  %ld\n",chrLen);
	
	fDatap = fopen(chrZntFn,"r");
	if (fDatap == NULL)
		{
		printf("\nOpen chrZntFn r command failed - %s\n",chrZntFn);
		exit(0);
		}

	/* mark genome */

	m = 0;
	errPrt = 0;
	bbMax = 0;
	posMax = 0;
	bbase = 0;
	do
		{
		/* fscanf(fDatap," %s %ld %c",chrVal,&cp,&orient); */
		fscanf(fDatap," %s %ld %c %s",chrVal,&cp,&orient,seqRead);  /* seq field added 5/25/2010 */
		if(!strcmp(chrSet,chrVal))
			{
			if(m < 5)
				{
				printf("chrSet %s, chrVal %s, cp %ld, orient %c\n",chrSet,chrVal,cp,orient);
				m++;
				}
			if(cp >= chrLen)
				{
				printf("error - chrLen exceeded\n");
				exit(0);
				}

			if(orient == '+')
				{
				for(n = 0; n < readLen; n++)
					{
					fseek(fp,2 * (cp - 1) + (2 * n),0);  /* note increments of 2 */
					bbChk = fread(&bbase,sizeof(short int),1,fp);
					bbase = bbase + 1;

					if(bbase <= 65000)
						{
						fseek(fp,-2,1);
						bbChk = fwrite(&bbase,sizeof(short int),1,fp);
						}
					else
						{
						if(errPrt < 20)
							{
							printf("error - mark lim 65536 exceeded; bbase %d,  cp %ld\n",bbase,cp);
							errPrt++;
							}
						}
					if(bbase > bbMax)
						{
						bbMax = bbChk;
						posMax = cp;
						}
					}
				}
			else if(orient == '-')
				{
				for(n = 0; n < readLen; n++)
					{
					fseek(fp,2 * (cp + readLen - 1) + (2 * n),0);  /* note increments of 2 */
					bbChk = fread(&bbase,sizeof(short int),1,fp);
					bbase = bbase + 1;

					if(bbase <= 65000)
						{
						fseek(fp,-2,1);
						bbChk = fwrite(&bbase,sizeof(short int),1,fp);
						}
					else
						{
						if(errPrt < 20)
							{
							printf("error - mark lim 65536 exceeded; bbase %d,  cp %ld\n",bbase,cp + readLen);
							errPrt++;
							}
						}
					if(bbase > bbMax)
						{
						bbMax = bbChk;
						posMax = cp;
						}
					}
				}
			}
		}while(!feof(fDatap));
		
	printf("bbMax %d,  posMax %ld\n\n",bbMax,posMax);

	fclose(fp);
	fclose(fDatap);

	return 0;
}

