/****************************************************************************
 Copyright (C) 2013 SUNCAT
 This file is distributed under the terms of the
 GNU General Public License. See the file `COPYING'
 in the root directory of the present distribution,
 or http://www.gnu.org/copyleft/gpl.txt .


 cubecutperiodic removes periodic image planes (x=1,y=1,z=1) from
 cube files produced by espresso's pp.x for analysis tools
 that expect input without these planes
 (e.g. the bader tool from the Henkelman group)
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#define nbuf 65536

void readerror(const char *s) {
    fprintf(stderr, "error reading from file %s\n", s);
    exit(4);
}

int main(int argc, char **argv) {
    char buf[nbuf],buf2[512],buf3[512],buf4[512];
    FILE *inp,*out;
    int i,j,k,natoms,nx,ny,nz,nx1,ny1,nz1,counter;

    if(argc!=3) {
	fprintf(stderr, "usage: %s in.cube out.cube\n", argv[0]);
	return 1;
    }
    
    inp = fopen(argv[1], "r");
    if(!inp) {
	fprintf(stderr, "cannot open cube input file %s for reading.\n", argv[1]);
	return 2;
    }

    out = fopen(argv[2], "w");
    if(!out) {
	fprintf(stderr, "cannot open cube output file %s for writing.\n", argv[2]);
	return 3;
    }
    
    for(i=0;i<3;i++) {
	if(!fgets(buf, nbuf, inp)) readerror(argv[1]);
	fputs(buf, out);
    }

    sscanf(buf, "%d", &natoms);
    
    if(!fgets(buf, nbuf, inp)) readerror(argv[1]);
    sscanf(buf, "%d%512s%512s%512s", &nx, buf2, buf3, buf4);
    fprintf(out, " %d %s %s %s\n", nx-1, buf2, buf3, buf4);
    
    if(!fgets(buf, nbuf, inp)) readerror(argv[1]);
    sscanf(buf, "%d%512s%512s%512s", &ny, buf2, buf3, buf4);
    fprintf(out, " %d %s %s %s\n", ny-1, buf2, buf3, buf4);
    
    if(!fgets(buf, nbuf, inp)) readerror(argv[1]);
    sscanf(buf, "%d%512s%512s%512s", &nz, buf2, buf3, buf4);
    fprintf(out, " %d %s %s %s\n", nz-1, buf2, buf3, buf4);

    for(i=0;i<natoms;i++) {
	if(!fgets(buf, nbuf, inp)) readerror(argv[1]);
	fputs(buf, out);
    }
    
    counter = 0;
    nx1 = nx-1;
    ny1 = ny-1;
    nz1 = nz-1;
    for(i=0;i<nx;i++) {
	for(j=0;j<ny;j++) {
	    for(k=0;k<nz;k++) {
		fscanf(inp, "%512s", buf2);
		if(i!=nx1 && j!=ny1 && k!= nz1) {
		    fprintf(out, " %s", buf2);
		}
		counter++;
		if(counter==6) {
		    fputs("\n", out);
		    counter = 0;
		}
	    }
	    if(counter!=0) {
		fputs("\n", out);
		counter = 0;
	    }
	}
    }

    fclose(out);
    fclose(inp);
    
    return 0;
}
