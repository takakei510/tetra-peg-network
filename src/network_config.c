#include "network_config.h"
#include "molecule.h"
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int fail(char *error, size_t size, int line, const char *message) {
    if(error && size) snprintf(error,size,"line %d: %s",line,message);
    return -1;
}
static char *trim(char *s) {
    while(isspace((unsigned char)*s)) ++s;
    char *end=s+strlen(s);
    while(end>s && isspace((unsigned char)end[-1])) --end;
    *end=0; return s;
}
static int unsigned_value(const char *s, uintmax_t limit, uintmax_t *out) {
    /* 整数は10進数字のみ。負数・小数・符号・後続文字を黙って変換しない。 */
    if(!*s) return -1;
    for(const char *p=s;*p;++p) if(!isdigit((unsigned char)*p)) return -1;
    errno=0; char *end; uintmax_t value=strtoumax(s,&end,10);
    if(errno || *end || value>limit) return -1;
    *out=value; return 0;
}

int network_config_validate(const NetworkConfig *c, char *error, size_t size) {
    if(!c) return fail(error,size,0,"null config");
    if(c->base.dim!=3) return fail(error,size,0,"AB-SAW requires dim=3");
    if(c->base.L<2) return fail(error,size,0,"L must be >= 2");
    if(!c->target_molecules || c->target_molecules>(size_t)INT_MAX ||
       c->target_molecules>SIZE_MAX/(ARM_COUNT*24))
        return fail(error,size,0,"target_molecules exceeds supported positive range");
    if(!c->arm_length) return fail(error,size,0,"arm_length must be >= 1");
    if(!isfinite(c->a_probability) || c->a_probability<0 || c->a_probability>1)
        return fail(error,size,0,"a_probability must be finite and in [0,1]");
    /* lattice_initと同じ格子点数上限。ここではowner配列を確保しない。 */
    size_t sites=1;
    for(int d=0;d<3;++d) {
        if(sites>SIZE_MAX/(size_t)c->base.L) return fail(error,size,0,"L cubed overflows");
        sites*=(size_t)c->base.L;
    }
    if(sites>(size_t)INT_MAX || sites>SIZE_MAX/sizeof(int) || sites>SIZE_MAX/sizeof(size_t))
        return fail(error,size,0,"lattice/index size exceeds implementation limit");
    /* 4腕と中心の必要点数を、オーバーフロー前に検査する。 */
    if(c->arm_length>(SIZE_MAX-1)/ARM_COUNT)
        return fail(error,size,0,"arm_length exceeds size limit");
    size_t path_points=ARM_COUNT*c->arm_length, per_molecule=1+path_points;
    if(per_molecule>SIZE_MAX/sizeof(size_t) ||
       c->target_molecules>SIZE_MAX/sizeof(Molecule) ||
       c->target_molecules>SIZE_MAX/path_points/sizeof(size_t))
        return fail(error,size,0,"trajectory allocation size overflows");
    /* 必要点数が全格子点を超える設定は生成前に拒否。成功保証ではない。 */
    if(c->target_molecules>sites/per_molecule)
        return fail(error,size,0,"required occupied sites exceed L cubed");
    return 0;
}

int network_config_load(const char *path, NetworkConfig *out, char *error, size_t size) {
    if(!path || !out) return fail(error,size,0,"null argument");
    FILE *fp=fopen(path,"r");
    if(!fp) return fail(error,size,0,"cannot open config");
    const char *keys[]={"dim","L","seed","target_molecules","arm_length","a_probability"};
    NetworkConfig config={0}; unsigned seen=0; char buffer[512]; int line=0,result=0;
    while(fgets(buffer,sizeof buffer,fp)) {
        ++line;
        if(!strchr(buffer,'\n') && !feof(fp)) {result=fail(error,size,line,"line too long");break;}
        /* 既存cfgと同じkey=value、#コメント、前後空白の書式。 */
        char *comment=strchr(buffer,'#');if(comment) *comment=0;
        char *s=trim(buffer);if(!*s) continue;
        char *eq=strchr(s,'=');if(!eq) {result=fail(error,size,line,"expected key=value");break;}
        *eq=0;char *key=trim(s),*value=trim(eq+1);size_t k;
        for(k=0;k<6;++k) if(!strcmp(key,keys[k])) break;
        if(k==6) {result=fail(error,size,line,"unknown key");break;}
        unsigned bit=1u<<k;
        if(seen&bit) {result=fail(error,size,line,"duplicate key");break;}seen|=bit;
        if(k==5) {
            errno=0;char *end;double probability=strtod(value,&end);
            if(errno || end==value || *end || !isfinite(probability)) {
                result=fail(error,size,line,"invalid probability");break;
            }
            config.a_probability=probability;
        } else {
            uintmax_t v,limit=k<2?INT_MAX:k==2?UINT64_MAX:SIZE_MAX;
            if(unsigned_value(value,limit,&v)) {result=fail(error,size,line,"invalid integer");break;}
            if(k==0) config.base.dim=(int)v;
            else if(k==1) config.base.L=(int)v;
            else if(k==2) config.base.seed=(uint64_t)v;
            else if(k==3) config.target_molecules=(size_t)v;
            else config.arm_length=(size_t)v;
        }
    }
    if(ferror(fp)) result=fail(error,size,line,"read error");
    if(fclose(fp) && !result) result=fail(error,size,line,"close error");
    /* 既定値を補わず、入力し忘れを必ず知らせる。 */
    if(!result && seen!=63) result=fail(error,size,line,"all six network config keys are required");
    if(!result) result=network_config_validate(&config,error,size);
    if(!result) *out=config;
    return result;
}
