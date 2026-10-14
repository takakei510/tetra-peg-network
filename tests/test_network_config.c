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
    assert(c.target_mode==TARGET_MOLECULE_COUNT && c.requested_phi==0 && c.phi_site_budget==0);
    assert(!network_config_load("configs/network/preview_phi010_L16_seed3.cfg",&c,error,sizeof error));
    assert(c.target_mode==TARGET_OCCUPIED_FRACTION && c.requested_phi==0.1);
    assert(c.target_molecules==19 && c.phi_site_budget==409);
    /* 十進入力の切り捨て境界。doubleの誤差で5分子が4分子にならない。 */
    const char *phi_values[]={"0.105","0.104999999999999999","0.105000000000000001","1.05e-1","1",".105","+0.1050"};
    const size_t phi_counts[]={5,4,5,5,47,5,5};
    for(size_t k=0;k<sizeof phi_counts/sizeof *phi_counts;++k) {
        char text[256];snprintf(text,sizeof text,"dim=3\nL=10\nseed=3\narm_length=5\na_probability=0.5\ntarget_phi=%s\n",phi_values[k]);
        write_config(text);assert(!network_config_load(path,&c,error,sizeof error));
        assert(c.target_molecules==phi_counts[k]);
    }
    const char *bad_phi[]={"0","-0.1","1.1","1.000000000000000001","nan","inf","0x1p-1","1e999","1e-999","0.0001","0.1junk",".","1e","0.1.2"};
    for(size_t k=0;k<sizeof bad_phi/sizeof *bad_phi;++k) {
        char text[256];snprintf(text,sizeof text,"dim=3\nL=10\nseed=3\narm_length=5\na_probability=0.5\ntarget_phi=%s\n",bad_phi[k]);invalid(text);
    }
    invalid("dim=3\nL=16\nseed=3\narm_length=5\na_probability=0.5\ntarget_phi=0.1\ntarget_molecules=19\n");
    invalid("dim=3\nL=16\nseed=3\narm_length=5\na_probability=0.5\ntarget_phi=0.1\ntarget_phi=0.1\n");
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
