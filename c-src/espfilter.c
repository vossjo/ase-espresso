/****************************************************************************
 Copyright (C) 2013 SUNCAT
 This file is distributed under the terms of the
 GNU General Public License. See the file `COPYING'
 in the root directory of the present distribution,
 or http://www.gnu.org/copyleft/gpl.txt .


 espfilter increases the buffer size of a pipe if needed
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define bufstep (1024*1024)

void dump(char *b, size_t r)
{
    size_t i;
    
    for(i=0; i+1023<r; i+=1024)
	fwrite(b+i, 1, 1024, stdout);
    i = r%1024;
    if(i) fwrite(b+r-i, 1, i, stdout);
    fflush(stdout);
}

int main(int argc, char **argv)
{
    size_t bufsize=bufstep,currsize=bufstep,bytesread=0,i,nase;
    char *buf,*p;
    FILE *log;
    
    if(argc!=3)
    {
	fputs("Expecting number of atoms and log filename as parameters.\n", stderr);
	return 1;
    }
    
    nase = strtol(argv[1], NULL, 10)*2;
    
    log = fopen(argv[2], "a");
    if(!log)
    {
	fprintf(stderr, "Cannot open log file %s for writing.\n", argv[2]);
	return 1;
    }
    
    buf = (char *) malloc(bufstep);
    if(!buf)
    {
	fputs("Buffer memory allocation error!\n", stderr);
	return 1;
    }
    p = buf;

    while(fgets(p, 65536, stdin))
    {
	fputs(p, log);
	fflush(log);
	if(strstr(p,"   total energy")
	    || strstr(p,"     stopping")
	    || strstr(p,"     convergence NOT")) {
	    fputs(p, stderr);
	    fflush(stderr);
	}
	else if(strstr(p,"     Writing output data file"))
	{
	    dump(buf, bytesread+strlen(p));
	    p = buf;
	    bytesread = 0;
	    continue;
	}
	else if(strstr(p,"!ASE"))
	{
	    dump(buf, bytesread+strlen(p));
	    p = buf;
	    bytesread = 0;
	    for(i=0;i<nase;i++)
	    {
		fgets(buf, 65536, stdin);
		fputs(buf, stdout);
	    }
	    fflush(stdout);
	    continue;
	}
	i = strlen(p);
	p += i;
	bytesread += i;
	if(bytesread+65536>=currsize)
	{
	    currsize += bufstep;
	    buf = (char *) realloc(buf, currsize);
	    if(!buf)
	    {
		fputs("Buffer memory re-allocation error!\n", stderr);
		return 1;
	    }
	    p = buf+bytesread;
	}
    }
    
    if(bytesread>0) dump(buf, bytesread);

    fclose(log);

    return 0;
}
