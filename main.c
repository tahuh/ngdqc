/*
 * main.c
 * NGS qc main file
 * compile with gcc -O3 -Wall -o ngsqc main.c -lz
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

const char *usage = "ngsqc [list of gastq files] > qc.txt";
int main(int argc, char **argv) 
{
    char *fname = NULL;
    gzFile fastq_fp = NULL;
    kseq_t *seq_ent = NULL;
    // char *seq = NULL;
    char *qual = NULL;
    char qbase = 0;
    unsigned long long base_cnt = 0;
    unsigned long long base_q20_cnt = 0;
    unsigned long long base_q30_cnt = 0;
    unsigned long long reads = 0;
    float q20 = 0.0;
    float q30 = 0.0;
    int l = 0;
    if(argc < 2) {
        fprintf(stderr, "%s\n", usage);
        return -1;
    }

    size_t i = 0;
    size_t j = 0;
    fprintf(stderr, "[ngsqc] Process trigger for %d files\n", argc - 1);
    fprintf(stdout,"filename\ttotal_bases\ttotal_reads\tq20_bases\tq30_bases\tq20_percentage\tq30_percentage\n");
    for(i = 1; i < argc; i++){
        // initialize all
        fname = NULL; fastq_fp = NULL; seq_ent = NULL; qual = NULL; // seq = NULL;
        qbase = 0; base_cnt = 0; base_q20_cnt = 0; base_q30_cnt = 0; reads = 0;
        q20 = 0.0; q30 = 0.0;
        l = 0;
        j = 0;
        fprintf(stderr, "[ngsqc] Processing file = %s\n", argv[i]);
        fname = argv[i];
        fastq_fp = gzopen(fname, "r");
        if(fastq_fp == NULL){
             fprintf(stderr, "[ngsqc] Failed to open file = %s\n", fname);
             return -1;
        }
        seq_ent = kseq_init(fastq_fp);
        if(seq_ent == NULL){
             fprintf(stderr, "[ngsqc] failed to init kseq. \n");
             gzclose(fastq_fp);
             return -1;
        }
        while((l = kseq_read(seq_ent)) >= 0) {
            reads++;
            //seq = seq_ent->seq.s;
            qual = seq_ent->qual.s;
            base_cnt += seq_ent->seq.l;
            for(j = 0; j < seq_ent->seq.l; j++){
                qbase = qual[j] - 33; // illumina phred quality
                if(qbase >= 20) {
                    base_q20_cnt++;
                }
                if(qbase >= 30){
                    base_q30_cnt++;
                }
            }
        }
        kseq_destroy(seq_ent); gzclose(fastq_fp);
        q20 = base_q20_cnt / (float)(base_cnt);
        q20 = 100 * q20;
        q30 = base_q30_cnt / (float)(base_cnt);
        q30 = 100 * q30;
        fprintf(stdout, "%s\t%llu\t%llu\t%llu\t%llu\t%f\t%f\n", fname, base_cnt, reads, base_q20_cnt,base_q30_cnt,q20,q30);
    }
    fprintf(stderr, "[ngsqc] Done!\n");
    return 0;
}
