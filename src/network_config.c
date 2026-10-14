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

/* 十進の入力値を整数演算で扱う。0.105等の境界で浮動小数誤差により1分子減らさない。 */
static int phi_budget(const char *text, size_t sites, size_t *budget, double *display) {
    const char *p=text; if(*p=='+') ++p;
    char digits[512]; size_t count=0; int fraction=0, dot=0;
    while(isdigit((unsigned char)*p) || *p=='.') {
        if(*p=='.') {if(dot) return -1;dot=1;}
        else {if(count>=sizeof digits) return -1;digits[count++]=*p;if(dot) ++fraction;}
        ++p;
    }
    if(!count) return -1;
    long exponent=0;
    if(*p=='e' || *p=='E') {
        ++p;char *end;errno=0;exponent=strtol(p,&end,10);
        if(errno || end==p || exponent < -10000 || exponent > 10000) return -1;
        p=end;
    }
    if(*p) return -1; /* hex浮動小数、NaN、inf、後続文字は受け付けない。 */
    size_t first=0;while(first<count && digits[first]=='0') ++first;
    if(first==count) return -1; /* 占有率0はこの生成目標では認めない。 */
    memmove(digits,digits+first,count-first);count-=first;
    long shift=exponent-fraction;
    while(count>1 && digits[count-1]=='0') {--count;++shift;}
    /* 1を超える入力を、doubleで1へ丸まる前に拒否する。 */
    if((long)count+shift>0 && !(count==1 && digits[0]=='1' && shift==0)) return -1;
    unsigned char product[544];size_t length=0;uintmax_t carry=0;
    /* 有効桁の整数×格子点数を十進の配列へ保存。桁あふれする乗算を避ける。 */
    for(size_t k=count;k>0;--k) {
        uintmax_t value=(uintmax_t)(digits[k-1]-'0')*sites+carry;
        product[length++]=(unsigned char)(value%10);carry=value/10;
    }
    while(carry) {product[length++]=(unsigned char)(carry%10);carry/=10;}
    long integer_digits=(long)length+shift;size_t result=0;
    if(integer_digits>0) {
        for(long k=(long)length-1;k>=-shift;--k) {
            unsigned digit=product[k];
            if(result>(SIZE_MAX-digit)/10) return -1;
            result=result*10+digit;
        }
    }
    if(result>sites) return -1;
    /* 小数部を捨てた占有点上限を返す。分子数への切り捨ては呼出側で行う。 */
    char *end;errno=0;double phi=strtod(text,&end);
    if(errno || *end || !isfinite(phi) || phi<=0 || phi>1) return -1;
    *budget=result;*display=phi;return 0;
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
    if(c->target_mode!=TARGET_MOLECULE_COUNT && c->target_mode!=TARGET_OCCUPIED_FRACTION)
        return fail(error,size,0,"unknown target mode");
    if(c->target_mode==TARGET_OCCUPIED_FRACTION &&
       (!isfinite(c->requested_phi) || c->requested_phi<=0 || c->requested_phi>1 || c->phi_site_budget>sites))
        return fail(error,size,0,"invalid resolved target_phi");
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
    if(c->target_mode==TARGET_OCCUPIED_FRACTION && c->target_molecules!=c->phi_site_budget/per_molecule)
        return fail(error,size,0,"resolved molecule count does not match target_phi");
    return 0;
}

int network_config_load(const char *path, NetworkConfig *out, char *error, size_t size) {
    if(!path || !out) return fail(error,size,0,"null argument");
    FILE *fp=fopen(path,"r");
    if(!fp) return fail(error,size,0,"cannot open config");
    const char *keys[]={"dim","L","seed","target_molecules","arm_length","a_probability","target_phi"};
    NetworkConfig config={0}; unsigned seen=0; char buffer[512]; int line=0,result=0;char phi_text[512]={0};
    while(fgets(buffer,sizeof buffer,fp)) {
        ++line;
        if(!strchr(buffer,'\n') && !feof(fp)) {result=fail(error,size,line,"line too long");break;}
        /* 既存cfgと同じkey=value、#コメント、前後空白の書式。 */
        char *comment=strchr(buffer,'#');if(comment) *comment=0;
        char *s=trim(buffer);if(!*s) continue;
        char *eq=strchr(s,'=');if(!eq) {result=fail(error,size,line,"expected key=value");break;}
        *eq=0;char *key=trim(s),*value=trim(eq+1);size_t k;
        for(k=0;k<7;++k) if(!strcmp(key,keys[k])) break;
        if(k==7) {result=fail(error,size,line,"unknown key");break;}
        unsigned bit=1u<<k;
        if(seen&bit) {result=fail(error,size,line,"duplicate key");break;}seen|=bit;
        if(k==6) {
            if(strlen(value)>=sizeof phi_text) {result=fail(error,size,line,"target_phi too long");break;}
            strcpy(phi_text,value);
        } else if(k==5) {
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
    /* 共通5キーと目標の片方だけが必要。両方指定で優先順位を作らない。 */
    if(!result && (seen&55)!=55) result=fail(error,size,line,"dim,L,seed,arm_length,a_probability are required");
    if(!result && (!(seen&(8|64)) || (seen&(8|64))==(8|64)))
        result=fail(error,size,line,"specify exactly one of target_molecules or target_phi");
    if(!result && (seen&64)) {
        /* まず1分子の設定で共通制約を検査し、安全なL^3と1+4lを求める。 */
        config.target_molecules=1;
        result=network_config_validate(&config,error,size);
        if(!result) {
            size_t sites=(size_t)config.base.L*config.base.L*config.base.L;
            result=phi_budget(phi_text,sites,&config.phi_site_budget,&config.requested_phi);
            if(result) result=fail(error,size,line,"target_phi must be a decimal in (0,1]");
            else {
                config.target_mode=TARGET_OCCUPIED_FRACTION;
                config.target_molecules=config.phi_site_budget/(1+ARM_COUNT*config.arm_length);
                if(!config.target_molecules) result=fail(error,size,line,"target_phi resolves to zero molecules");
            }
        }
    }
    if(!result) result=network_config_validate(&config,error,size);
    if(!result) *out=config;
    return result;
}
