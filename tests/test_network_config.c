#include "network_config.h"
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static const char *path="build/test_network_config.cfg";
static void write_config(const char *text) {
    FILE *fp=fopen(path,"w");assert(fp);assert(fputs(text,fp)>=0);assert(!fclose(fp));
}
static void invalid(const char *text) {
    write_config(text); NetworkConfig c,before;char error[256]={0};
    memset(&c,0xA5,sizeof c);memcpy(&before,&c,sizeof c);
    assert(network_config_load(path,&c,error,sizeof error));
    assert(error[0] && !memcmp(&c,&before,sizeof c)); /* 失敗時に呼出元の設定を壊さない。 */
}
int main(void) {
    NetworkConfig c;char error[256];
    assert(!network_config_load("configs/network/preview_N20_L16_seed3.cfg",&c,error,sizeof error));
    assert(c.base.dim==3 && c.base.L==16 && c.base.seed==3);
    assert(c.target_molecules==20 && c.arm_length==5 && c.a_probability==0.5);
    write_config(" # comment\n a_probability = 2.5e-1 # r\narm_length=3\ntarget_molecules=12\nseed=18446744073709551615\nL=16\ndim=3\n");
    assert(!network_config_load(path,&c,error,sizeof error));
    assert(c.base.seed==UINT64_MAX && c.target_molecules==12 && c.arm_length==3 && c.a_probability==0.25);
    /* r=0と1、seed=0は有効。必要点数が合っても生成成功を保証はしない。 */
    c.a_probability=0;assert(!network_config_validate(&c,error,sizeof error));
    c.a_probability=1;assert(!network_config_validate(&c,error,sizeof error));
    write_config("dim=3\nL=2\nseed=0\ntarget_molecules=1\narm_length=1\na_probability=0\n");
    assert(!network_config_load(path,&c,error,sizeof error) && c.base.seed==0);
    const char *rows[]={"dim=3\n","L=16\n","seed=3\n","target_molecules=20\n","arm_length=5\n","a_probability=0.5\n"};
    for(int missing=0;missing<6;++missing) {
        char text[256]={0};for(int k=0;k<6;++k)if(k!=missing)strcat(text,rows[k]);invalid(text);
    }
    const char *bad_rows[]={"dim=2\n","L=1\n","seed=-1\n","target_molecules=0\n",
        "target_molecules=-3\n","target_molecules=1.5\n","arm_length=0\n","arm_length=18446744073709551615\n",
        "a_probability=-0.1\n","a_probability=1.1\n","a_probability=nan\n","a_probability=inf\n",
        "a_probability=1e999\n","a_probability=0.5junk\n","seed=18446744073709551616\n",
        "L=2147483648\n","L=1291\n","target_molecules=196\n"};
    for(size_t k=0;k<sizeof bad_rows/sizeof *bad_rows;++k) {
        char text[512]={0};char key[64];const char *eq=strchr(bad_rows[k],'=');
        size_t length=(size_t)(eq-bad_rows[k]);memcpy(key,bad_rows[k],length);key[length]=0;
        for(int j=0;j<6;++j)if(strncmp(rows[j],key,length) || rows[j][length]!='=')strcat(text,rows[j]);
        strcat(text,bad_rows[k]);invalid(text);
    }
    const char *base="dim=3\nL=16\nseed=3\ntarget_molecules=20\narm_length=5\na_probability=0.5\n";
    char text[512];snprintf(text,sizeof text,"%sdim=3\n",base);invalid(text);
    snprintf(text,sizeof text,"%sunknown=1\n",base);invalid(text);
    snprintf(text,sizeof text,"%sbad-line\n",base);invalid(text);
    char long_line[1024];memset(long_line,'x',sizeof long_line-1);long_line[sizeof long_line-1]=0;invalid(long_line);
    /* 古い3キーの読込はそのまま。新形式を古いデモへ渡すと明確にエラー。 */
    Config old;assert(!config_load("configs/default.cfg",&old,error,sizeof error));
    assert(config_load("configs/network/preview_N20_L16_seed3.cfg",&old,error,sizeof error));
    assert(network_config_load("configs/default.cfg",&c,error,sizeof error));
    assert(network_config_load("build/nonexistent_network_config.cfg",&c,error,sizeof error));
    assert(!remove(path));
    puts("network config tests passed (values, required keys, ranges, overflow, legacy compatibility)");
}
