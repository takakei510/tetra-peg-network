#include "config.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <limits.h>
static char *trim(char *s) { while (isspace((unsigned char)*s)) ++s; char *end=s+strlen(s); while(end>s && isspace((unsigned char)end[-1])) --end; *end=0; return s; }
static int fail(char *error,size_t n,int line,const char *message) { if(n) snprintf(error,n,"line %d: %s",line,message); return -1; }
int config_load(const char *path, Config *out, char *error, size_t error_size) {
    if(!path||!out) return fail(error,error_size,0,"null argument");
    FILE *fp=fopen(path,"r"); if(!fp) { if(error_size) snprintf(error,error_size,"cannot open %s",path); return -1; }
    Config c={0}; unsigned seen=0; char linebuf[512]; int line=0, result=0;
    while(fgets(linebuf,sizeof linebuf,fp)) {
        ++line;
        if(!strchr(linebuf,'\n') && !feof(fp)) { result=fail(error,error_size,line,"line too long"); break; }
        char *comment=strchr(linebuf,'#'); if(comment) *comment=0;
        char *s=trim(linebuf); if(!*s) continue;
        char *eq=strchr(s,'='); if(!eq) {result=fail(error,error_size,line,"expected key=value");break;}
        *eq=0; char *key=trim(s), *value=trim(eq+1), *end;
        unsigned bit; if(!strcmp(key,"dim")) bit=1; else if(!strcmp(key,"L")) bit=2; else if(!strcmp(key,"seed")) bit=4; else {result=fail(error,error_size,line,"unknown key");break;}
        if(seen&bit) {result=fail(error,error_size,line,"duplicate key");break;} seen|=bit;
        if(!*value||*value=='-') {result=fail(error,error_size,line,"invalid value");break;}
        errno=0; unsigned long long v=strtoull(value,&end,10);
        if(errno||*end||end==value||((bit!=4)&&v>INT_MAX)) {result=fail(error,error_size,line,"invalid number");break;}
        if(bit==1) c.dim=(int)v; else if(bit==2) c.L=(int)v; else c.seed=(uint64_t)v;
    }
    if(ferror(fp)) result=fail(error,error_size,line,"read error");
    fclose(fp);
    if(!result && seen!=7) result=fail(error,error_size,line,"dim, L and seed are required");
    if(!result && (c.dim!=2&&c.dim!=3)) result=fail(error,error_size,line,"dim must be 2 or 3");
    if(!result && c.L<2) result=fail(error,error_size,line,"L must be >= 2");
    if(!result) *out=c;
    return result;
}
