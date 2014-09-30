#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <fitsio.h>

void cp_4r4(unsigned char *dest, unsigned char *orig, int nitem){
    
    int i;
    for(i = 0; i < nitem; i++, orig += 4, dest += 4){
        
        dest[0] = orig[3];
        dest[1] = orig[2];
        dest[2] = orig[1];
        dest[3] = orig[0];
    }
}

int main(int argc, char **argv){
    
    if(argc < 2){
        printf("Usage: ./readuv <filename>\n");
        return 0;
    }

    fitsfile *fptr;
    char card[FLEN_CARD];
    int status = 0, nkeys, i;
    
    fits_open_file(&fptr, argv[1], READONLY, &status);
    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    printf("Number of keys: %d\n", nkeys);
    
    for(i = 1; i <= nkeys; ++i){
        fits_read_record(fptr, i, card, &status);
        if(strstr(card, "HISTORY") != NULL) break;
        printf("%s\n", card);
    }

    int hdutype;
    fits_get_hdu_type(fptr, &hdutype, &status);
    printf("HDU type: %d\n", hdutype);
    assert(hdutype == IMAGE_HDU);
    
    int bitpix;
    fits_get_img_type(fptr, &bitpix, &status);
    printf("bitpix: %d\n", bitpix);
    assert(bitpix == FLOAT_IMG);

    int nnaxis, naxis[128];
    fits_read_key(fptr, TINT, "NAXIS", &nnaxis, NULL, &status);
    printf("NAXIS: %d\n", nnaxis);
    assert(nnaxis < 128);

    char tmp[128];
    for(i = 1; i <= nnaxis; ++i){
        sprintf(tmp, "NAXIS%d", i);
        fits_read_key(fptr, TINT, tmp, &naxis[i], NULL, &status);
        printf("NAXIS%d: %d\n", i, naxis[i]);
    }

    int gcount, pcount;
    fits_read_key(fptr, TINT, "GCOUNT", &gcount, NULL, &status);
    fits_read_key(fptr, TINT, "PCOUNT", &pcount, NULL, &status);
    printf("GCOUNT: %d\n", gcount);
    printf("PCOUNT: %d\n", pcount);
    assert(pcount == 9);

    float nulval = 0;
    int anynul;
    int ncols = pcount + naxis[2] * naxis[5];
    int nelem = naxis[2] * naxis[5];
    float *buf_l = (float *)malloc(sizeof(float) * ncols * gcount);

    FILE *fp = fopen("out.dat", "w");
    assert(fp != NULL);

    float *buf = (float *)malloc(sizeof(float) * ncols * gcount);
    fits_read_tblbytes(fptr, 1, 1, ncols * gcount * 4, (unsigned char *)buf, &status);
    fits_report_error(stdout, status);
    cp_4r4((unsigned char *)buf_l, (unsigned char *)buf, ncols * gcount);
    free(buf);

    int ifcount = naxis[5]; 
    int offset;

    float *pflt = buf_l;
    for(i = 0; i < gcount; ++i){
//        fits_read_grppar_flt(fptr, 0, 1 + i, pcount, buf, &status);
//        fits_read_img(fptr, TFLOAT, 1 + i * ncols, pcount + nelem, &nulval, buf, &anynul, &status);
//        fprintf(fp, "%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n", 
//                pflt[0], pflt[1], pflt[2], 
//                pflt[pcount + 0], pflt[pcount + 1], pflt[pcount + 2]);
        fprintf(fp, "%14.5e%14.5e%14.5e", pflt[0], pflt[1], pflt[2]);
        for(offset = pcount; offset < pcount + 3 * ifcount; offset += 3){
            fprintf(fp, "%14.5e%14.5e%14.5e", 
                pflt[offset + 0], pflt[offset + 1], pflt[offset + 2]);
        }
        fprintf(fp, "\n");
        pflt += ncols;
    }
    fclose(fp);

    fits_close_file(fptr, &status);

    free(buf_l);

    return 0;
}
