#include "linearwaves.h"
#include <ctype.h>

#define PRINT_DOUBLE(NAME,VAL) printf("\t%s = %lg\n",NAME,VAL)
#define PRINT_INT(NAME,VAL) printf("\t%s = %d\n",NAME,VAL)
#define PRINT_STR(NAME,VAL) printf("\t%s = %s\n",NAME,VAL)
#define FPRINT_DOUBLE(F,NAME,VAL) fprintf(f,"%s = %lg\n",NAME,VAL)
#define FPRINT_INT(F,NAME,VAL) fprintf(f,"%s = %d\n",NAME,VAL)
#define FPRINT_STR(F,NAME,VAL) fprintf(f,"%s = %s\n",NAME,VAL)

void set_var(char *name,int int_val, double double_val, int bool_val, char *str_val, Params *params) {
    // Ints
    if (strcmp(name,"n") == 0) {
        params->n = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"nphi") == 0) {
        params->nphi = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"nrhs") == 0) {
        params->nrhs = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"mstart") == 0) {
        params->mstart = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"mend") == 0) {
        params->mend = int_val;
        PRINT_INT(name,int_val);
    }
    // Doubles
    else if (strcmp(name,"h") == 0) {	
        params->h = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"mu") == 0) {	
        params->mu = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"delta") == 0) {	
        params->delta = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"nuindx") == 0) {	
        params->nuindx = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"eta") == 0) {	
        params->eta = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"alpha") == 0) {	
        params->alpha = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"omf") == 0) {	
        params->omf = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"f") == 0) {	
        params->f = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"sig0") == 0) {	
        params->sig0 = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"ieps") == 0) {	
        params->ieps = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"rmin") == 0) {	
        params->rmin = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"rmax") == 0) {	
        params->rmax = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"eps") == 0) {	
        params->eps = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    // Bools
    else if (strcmp(name,"iso") == 0) {	
        params->iso = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"fromfile") == 0) {	
        params->fromfile = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"readomega") == 0) {	
        params->readomega = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"pcorrect") == 0) {	
        params->pcorrect = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"zero_inner_bc") == 0) {	
        params->zero_inner_bc = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"zero_outer_bc") == 0) {	
        params->zero_outer_bc = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"simple_visc") == 0) {	
        params->simple_visc = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"indirect") == 0) {	
        params->simple_visc = bool_val;
        PRINT_STR(name,str_val);

    }
    // Strings
    else if (strcmp(name,"outputname") == 0) {	
        sprintf(params->outputname,"%s",str_val);
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"diskfile") == 0) {	
        sprintf(params->diskfile,"%s",str_val);
        PRINT_STR(name,str_val);

    }

    return;
}

void parse_argument(int argc, char *argv[], Params *params) {
    int j;
    unsigned int i;
    char name[100],strval[100];
    double dval;
    int ival;
    int bool_val;
    char testbool;


    for(j=0;j<argc;j++) {
        sscanf(argv[j],"%32[^=]=%s",name,strval);
        dval = atof(strval);
        ival = atoi(strval);
        testbool = toupper(strval[0]);
        if (testbool == 'Y') bool_val = TRUE;
        else bool_val = FALSE;
        for (i = 0; i<strlen(name); i++) name[i] = (char)tolower(name[i]);
        set_var(name,ival,dval,bool_val,strval,params);
    }



    return;
}

void read_param_file(char *fname, int argc, char *argv[], Params *params) {
    FILE *f;

    char tok[20] = "\t :=>";

    char line[100],name[100],strval[100];
    char *data;
    double temp;
    int status;
    int int_val;
    int bool_val;
    char testbool;
    unsigned int i;

    f= fopen(fname,"r");

    while (fgets(line,100,f)) {
       // printf("%s\n",line);
        status = sscanf(line,"%s",name);

      //  printf("%s\n",name);
        if (name[0] != '#' && status == 1) {
        
             data = line + (int)strlen(name);
             sscanf(data + strspn(data,tok),"%lf",&temp);
             sscanf(data + strspn(data,tok),"%s",strval);
             //printf("%lf\t%s\n",temp,strval);
            int_val = (int)temp;
            testbool = toupper(strval[0]);
            if (testbool == 'Y') bool_val = TRUE;
            else bool_val = FALSE;
            
            for (i = 0; i<strlen(name); i++) name[i] = (char)tolower(name[i]);
            
            set_var(name,int_val,temp,bool_val,strval,params);

        }
    }

    printf("Redefined on the command line:\n");
    if (argc > 0) {
        parse_argument(argc,argv,params);
    }


    params->readomega = params->readomega && params->pcorrect;
    params->dlr = log(params->rmax/params->rmin) / (double)params->n;
    params->eps2 = params->eps*params->eps;
    params->indirect = params->indirect;
    params->a = 1;

    return;
}


