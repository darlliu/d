#!/bin/awk -f
BEGIN{
    FS="\t"
    tempfname="output.csv"
    program="./fgc.exe -b -o fgc_genesym_affyhg.bin"
    }
{

    print $2 | program
    close (program)
    getline result < tempfname
    system("rm " tempfname)
    print $1 "\t" result;
    }
