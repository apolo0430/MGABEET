#include "structures.h"
#include <memory.h>
#include <ctime>

pairing_t pairing;
element_t r1, r2;   // random value for Hash

/* System setup and Key Generation */
SET InitialSetup(securityParam p, attUniverseNum u)
{
    //std::cout << "[Initial Setup]" << std::endl;
    //std::cout << "start" << std::endl;
    SET set;

    pbc_param_t param;
    pbc_param_init_a_gen(param, 2 * p, 6 * p);
    pairing_init_pbc_param(pairing, param);

    /* temp value */
    element_t h, g1, g2;
    element_t alp, bet, a;

    element_init_G1(h, pairing);
    element_init_G1(g1, pairing);
    element_init_G1(g2, pairing);

    element_init_Zr(alp, pairing);
    element_init_Zr(bet, pairing);
    element_init_Zr(a, pairing);

    ////

    element_init_G1(set.sp.g, pairing);
    for (int i = 0; i < u; i++)
        element_init_G1(set.sp.F[i], pairing);
    element_init_G1(set.sp.u, pairing);
    element_init_G1(set.sp.v, pairing);
    element_init_G1(set.sp.w, pairing);
    element_init_G1(set.sp.ga, pairing);
    element_init_G1(set.msk.ha, pairing);
    element_init_G1(set.msk.hb, pairing);

    element_init_GT(set.sp.eg1h, pairing);
    element_init_GT(set.sp.eg2h, pairing);

    element_init_Zr(set.msk.k1, pairing);
    element_init_Zr(set.msk.k2, pairing);
    element_init_Zr(set.msk.k3, pairing);
    
    element_random(set.sp.g);
    element_random(h);
    for (int i = 0; i < u; i++)
        element_random(set.sp.F[i]);
    element_random(set.msk.k1);
    element_random(set.msk.k2);
    element_random(set.msk.k3);
    element_random(alp);
    element_random(bet);
    element_random(a);
    element_pow_zn(set.sp.u, set.sp.g, set.msk.k1);
    element_pow_zn(set.sp.v, set.sp.g, set.msk.k2);
    element_pow_zn(set.sp.w, set.sp.g, set.msk.k3);
    element_pow_zn(g1, set.sp.g, alp);
    element_pow_zn(g2, set.sp.g, bet);
    element_pow_zn(set.sp.ga, set.sp.g, a);
    element_pairing(set.sp.eg1h, g1, h);
    element_pairing(set.sp.eg2h, g2, h);

    element_pow_zn(set.msk.ha, h, alp);
    element_pow_zn(set.msk.hb, h, bet);

    ////

    element_clear(h);
    element_clear(g1);
    element_clear(g2);
    element_clear(alp);
    element_clear(bet);
    element_clear(a);

    //std::cout << "finish" << std::endl;

    // sample random values for hash function
    element_init_Zr(r1, pairing);
    element_init_Zr(r2, pairing);
    element_random(r1);
    element_random(r2);

    return set;
}

GSK SYS_FUNC::GroupKeyGen()
{
    //std::cout << "[Group Key Generation]" << std::endl;
    //std::cout << "start" << std::endl;

    GSK gsk;

    element_init_Zr(gsk.t1, pairing);
    element_init_Zr(gsk.t2, pairing);

    element_random(gsk.t1);
    element_random(gsk.t2);

    //std::cout << "finish" << std::endl;

    return gsk;
}

GW GA_FUNC::Join(SP sp, MSK msk, GSK gsk, memberAttSet s[MEM_ATT_LEN])
{
    //std::cout << "[Join]" << std::endl;
    //std::cout << "start" << std::endl;

    GW gw;

    /* temp value */
    element_t gat1;

    element_init_G1(gat1, pairing);

    ////
    
    element_init_G1(gw.hb_gat1, pairing);
    element_init_G1(gw.gt1, pairing);
    for (int i = 0; i < MEM_ATT_LEN; i++)
        element_init_G1(gw.fxt1[i], pairing);
    
    element_init_GT(gw.hS, pairing);
    element_init_GT(gw.hSt1, pairing);

    element_init_Zr(gw.t2, pairing);

    H1(gw.hS, s);
    element_pow_zn(gw.hSt1, gw.hS, gsk.t1);
    element_pow_zn(gat1, sp.ga, gsk.t1);
    element_mul(gw.hb_gat1, msk.hb, gat1);
    for (int i = 0; i < MEM_ATT_LEN; i++)
        element_pow_zn(gw.fxt1[i], sp.F[s[i]], gsk.t1);
    element_set(gw.t2, gsk.t2);

    ////

    element_clear(gat1);

    //std::cout << "finish" << std::endl;

    return gw;
}

MK GA_FUNC::Extract(SP sp, MSK msk, memberAttSet s[MEM_ATT_LEN])
{
    //std::cout << "[Extract]" << std::endl;
    //std::cout << "start" << std::endl;

    MK mk;

    /* temp value */
    element_t d, t, gad, gat;

    element_init_G1(gad, pairing);
    element_init_G1(gat, pairing);
    element_init_Zr(d, pairing);
    element_init_Zr(t, pairing);

    ////

    element_init_G1(mk.Kd, pairing);
    element_init_G1(mk.Kt, pairing);
    element_init_G1(mk.gd, pairing);
    element_init_G1(mk.gt, pairing);
    for (int i = 0; i < MEM_ATT_LEN; i++) {
        element_init_G1(mk.fxd[i], pairing);
        element_init_G1(mk.fxt[i], pairing);
    }

    extStart = clock();
    
    element_random(d);
    element_random(t);

    element_pow_zn(gad, sp.ga, d);
    element_mul(mk.Kd, msk.ha, gad);
    element_pow_zn(gat, sp.ga, t);
    element_mul(mk.Kt, msk.hb, gat);
    for (int i = 0; i < MEM_ATT_LEN; i++) {
        element_pow_zn(mk.fxd[i], sp.F[s[i]], d);
        element_pow_zn(mk.fxt[i], sp.F[s[i]], t);
    }
    element_pow_zn(mk.gd, sp.g, d);
    element_pow_zn(mk.gt, sp.g, t);

    extEnd = clock();

    ////

    extTime += extEnd - extStart;

    element_clear(gad);
    element_clear(gat);
    element_clear(d);
    element_clear(t);

    //std::cout << "finish" << std::endl;

    return mk;
}

/* Original Encryption and Decryption */
OCT MEM_FUNC::originalEnc(SP sp, GW gw, AP ap, MSG msg)
{
    //std::cout << "[Encryption]" << std::endl;
    //std::cout << "start" << std::endl;

    OCT ct;

    /* temp value */
    element_t eg1hs, H2m, eg2hs;
    element_t uT_vc6_w;
    element_t w[COLUMN];    // (s, y_1,..., y_n)
    element_t l[ROW];       // l_i = M_i x w
    element_t r[ROW];   // (r_1,..., r_l), r_0 = c6
    element_t mul, neg, one, T;

    element_init_GT(eg1hs, pairing);
    element_init_GT(H2m, pairing);
    element_init_GT(eg2hs, pairing);

    element_init_G1(uT_vc6_w, pairing);

    for (int i = 0; i < COLUMN; i++)
        element_init_Zr(w[i], pairing);
    for (int i = 0; i < ROW; i++)
        element_init_Zr(l[i], pairing);
    for (int i = 0; i < ROW; i++)
        element_init_Zr(r[i], pairing);
    element_init_Zr(mul, pairing);
    element_init_Zr(neg, pairing);
    element_init_Zr(one, pairing);
    element_init_Zr(T, pairing);
    
    ////

    element_init_G1(ct.c3, pairing);
    element_init_G1(ct.c8, pairing);
    for (int i = 0; i < ROW; i++) {
        element_init_G1(ct.c7[i][0], pairing);
        element_init_G1(ct.c7[i][1], pairing);
    }

    element_init_GT(ct.c1, pairing);
    element_init_GT(ct.c2, pairing);
    element_init_GT(ct.c4, pairing);
    element_init_GT(ct.c5, pairing);

    element_init_Zr(ct.c6, pairing);

    originEncStart = clock();

    for (int i = 0; i < COLUMN; i++)
        element_random(w[i]);
    for (int i = 0; i < ROW; i++) {
        element_set0(l[i]);

        for (int j = 0; j < COLUMN; j++) {
            element_mul(mul, ap.M[i][j], w[j]);
            element_add(l[i], l[i], mul);
        }
    }
    element_random(ct.c6);
    for (int i = 0; i < ROW; i++)
        element_random(r[i]);
    
    /* c1 */
    element_pow_zn(eg1hs, sp.eg1h, w[0]);
    element_mul(ct.c1, msg, eg1hs);
    /* c2 */
    H2(H2m, msg);
    element_pow_zn(eg2hs, sp.eg2h, w[0]);
    element_mul(ct.c2, H2m, eg2hs);
    /* c3 */
    element_pow_zn(ct.c3, sp.g, w[0]);
    /* c4 */
    element_pow_zn(ct.c4, gw.hS, ct.c6);
    /* c5 */
    element_pow2_zn(ct.c5, H2m, gw.t2, gw.hSt1, ct.c6);
    /* c7 */
    for (int i = 0; i < ROW; i++) {
        element_neg(neg, r[i]);

        element_pow2_zn(ct.c7[i][0], sp.ga, l[i], sp.F[ap.map[i]], neg);
        element_pow_zn(ct.c7[i][1], sp.g, r[i]);
    }
    /* c8 */
    element_set1(one);
    H3(T, ct.c1, ct.c2, ct.c4, ct.c7);
    element_pow3_zn(uT_vc6_w, sp.u, T, sp.v, ct.c6, sp.w, one);
    element_pow_zn(ct.c8, uT_vc6_w, w[0]);

    originEncEnd = clock();

    ////

    originEncTime += originEncEnd - originEncStart;

    element_clear(eg1hs);
    element_clear(H2m);
    element_clear(eg2hs);
    element_clear(uT_vc6_w);
    for (int i = 0; i < COLUMN; i++)
        element_clear(w[i]);
    for (int i = 0; i < ROW; i++)
        element_clear(l[i]);
    for (int i = 0; i < ROW; i++)
        element_clear(r[i]);
    element_clear(mul);
    element_clear(neg);
    element_clear(one);
    element_clear(T);

    //std::cout << "finish" << std::endl;

    return ct;
}

void MEM_FUNC::originalDec(SP sp, OCT ct, MK mk, AP ap)
{
    //std::cout << "[Decryption]" << std::endl;
    //std::cout << "start" << std::endl;

    /* temp value */
    element_t uT_vc6_w;
    element_t ec3uT_vc6_w, ec8g,
              ec3Kd, ec3Kt,
              eghas, eghbs,
              ecL, edK, ecL_edK,
              Hh, mh, H2mh,
              mul, prod;
    element_t T, one, zero;

    element_init_G1(uT_vc6_w, pairing);

    element_init_GT(ec3uT_vc6_w, pairing);
    element_init_GT(ec8g, pairing);
    element_init_GT(ec3Kd, pairing);
    element_init_GT(ec3Kt, pairing);
    element_init_GT(eghas, pairing);
    element_init_GT(eghbs, pairing);
    element_init_GT(ecL, pairing);
    element_init_GT(edK, pairing);
    element_init_GT(ecL_edK, pairing);
    element_init_GT(Hh, pairing);
    element_init_GT(mh, pairing);
    element_init_GT(H2mh, pairing);
    element_init_GT(mul, pairing);
    element_init_GT(prod, pairing);

    element_init_Zr(T, pairing);
    element_init_Zr(one, pairing);
    element_init_Zr(zero, pairing);
    
    ////

    originDecStart = clock();

    element_set1(one);
    H3(T, ct.c1, ct.c3, ct.c4, ct.c7);
    element_pow3_zn(uT_vc6_w, sp.u, T, sp.v, ct.c6, sp.w, one);
    element_pairing(ec3uT_vc6_w, ct.c3, uT_vc6_w);
    element_pairing(ec8g, ct.c8, sp.g);

    
    if (element_cmp(ec3uT_vc6_w, ec8g) == 1) {
        //std::cout << "CHECK FAIL" << std::endl;

        return;
    }

    element_pairing(ec3Kd, ct.c3, mk.Kd);
    element_set1(mul);
    for (int i = 0; i < ROW; i++) {
        element_pairing(ecL, ct.c7[i][0], mk.gd);
        if (ap.map[i] < MEM_ATT_LEN)
            element_pairing(edK, ct.c7[i][1], mk.fxd[ap.map[i]]);
        element_mul(ecL_edK, ecL, edK);
        if (i == 0)
            element_pow_zn(prod, ecL_edK, one);
        else
            element_pow_zn(prod, ecL_edK, zero);
        element_mul(mul, mul, prod);
    }
    element_div(eghas, ec3Kd, mul);

    element_pairing(ec3Kt, ct.c3, mk.Kt);
    element_set1(mul);
    for (int i = 0; i < ROW; i++) {
        element_pairing(ecL, ct.c7[i][0], mk.gt);
        if (ap.map[i] < MEM_ATT_LEN)
            element_pairing(edK, ct.c7[i][1], mk.fxt[ap.map[i]]);
        element_mul(ecL_edK, ecL, edK);
        if (i == 0)
            element_pow_zn(prod, ecL_edK, one);
        else
            element_pow_zn(prod, ecL_edK, zero);
        element_mul(mul, mul, prod);
    }
    element_div(eghbs, ec3Kt, mul);

    element_div(Hh, ct.c2, eghbs);
    element_div(mh, ct.c1, eghas);

    H2(H2mh, mh);
    /*
    if (element_cmp(Hh, H2mh) == 0)
        std::cout << "DECRYPTION SUCCEED" << std::endl;
    else
        std::cout << "DECRYPTION FAILED" << std::endl;
    */

    originDecEnd = clock();

    ////

    originDecTime += originDecEnd - originDecStart;

    element_clear(ec3uT_vc6_w);
    element_clear(ec8g);
    element_clear(ec3Kd);
    element_clear(ec3Kt);
    element_clear(eghas);
    element_clear(eghbs);
    element_clear(ecL);
    element_clear(edK);
    element_clear(ecL_edK);
    element_clear(Hh);
    element_clear(mh);
    element_clear(H2mh);
    element_clear(T);
    element_clear(one);
    element_clear(mul);
    element_clear(prod);
    element_clear(zero);

    //std::cout << "finish" << std::endl;
}

/* original equality test */
GTD GA_FUNC::originalGroupTrap(SP sp, GSK gsk)
{
    GTD gtd;

    /* temp value */
    element_t z, zt1;

    element_init_Zr(z, pairing);
    element_init_Zr(zt1, pairing);

    ////

    element_random(z);
    element_pow_zn(gtd.gtd1, sp.ga, z);
    element_mul(zt1, z, gsk.t1);
    element_pow_zn(gtd.gtd2, sp.ga, zt1);

    ////

    element_clear(z);

    return gtd;
}

void ATP_FUNC::originalGroupTest(OCT ct1, OCT ct2, GTD gtd)
{
    /* temp value */
    element_t ec5gtd1, ec4gtd2,
              c5_inv, c4_inv;

    element_init_GT(ec5gtd1, pairing);
    element_init_GT(ec4gtd2, pairing);
    element_init_GT(c5_inv, pairing);
    element_init_GT(c4_inv, pairing);

    ////

    element_mul(c5_inv, ct1.c5, ct2.c5);
    element_pairing(ec5gtd1, c5_inv, gtd.gtd1);
    element_mul(c4_inv, ct1.c4, ct2.c4);
    element_pairing(ec4gtd2, c4_inv, gtd.gtd2);

    /*
    if (element_cmp(ec5gtd1, ec4gtd2) == 0)
        std::cout << "SAME" << std::endl;
    else
        std::cout << "DIFFERENT" << std::endl;
    */

    ////

    element_clear(ec5gtd1);
    element_clear(ec4gtd2);
    element_clear(c5_inv);
    element_clear(c4_inv);
}

/* modified Encryption and Decryption */
MCT MEM_FUNC::modifiedEnc(SP sp, GW gw, AP ap, MSG msg)
{
    //std::cout << "[Encryption]" << std::endl;
    //std::cout << "start" << std::endl;

    MCT ct;

    /* temp value */
    element_t eg1hs, H2m;
    element_t uT_vc5_w;
    element_t w[COLUMN];    // (s, y_1,..., y_n)
    element_t l[ROW];       // l_i = M_i x w
    element_t r[ROW];   // (r_1,..., r_l), r_0 = c6
    element_t mul, neg, one, T;

    element_init_GT(eg1hs, pairing);
    element_init_GT(H2m, pairing);

    element_init_G1(uT_vc5_w, pairing);

    for (int i = 0; i < COLUMN; i++)
        element_init_Zr(w[i], pairing);
    for (int i = 0; i < ROW; i++)
        element_init_Zr(l[i], pairing);
    for (int i = 0; i < ROW; i++)
        element_init_Zr(r[i], pairing);
    element_init_Zr(mul, pairing);
    element_init_Zr(neg, pairing);
    element_init_Zr(one, pairing);
    element_init_Zr(T, pairing);
    
    ////

    element_init_G1(ct.c2, pairing);
    element_init_G1(ct.c7, pairing);
    for (int i = 0; i < ROW; i++) {
        element_init_G1(ct.c6[i][0], pairing);
        element_init_G1(ct.c6[i][1], pairing);
    }

    element_init_GT(ct.c1, pairing);
    element_init_GT(ct.c3, pairing);
    element_init_GT(ct.c4, pairing);

    element_init_Zr(ct.c5, pairing);

    modifiedEncStart = clock();

    for (int i = 0; i < COLUMN; i++)
        element_random(w[i]);
    for (int i = 0; i < ROW; i++) {
        element_set0(l[i]);

        for (int j = 0; j < COLUMN; j++) {
            element_mul(mul, ap.M[i][j], w[j]);
            element_add(l[i], l[i], mul);
        }
    }
    element_random(ct.c5);
    for (int i = 0; i < ROW; i++)
        element_random(r[i]);
    
    /* c1 */
    element_pow_zn(eg1hs, sp.eg1h, w[0]);
    element_mul(ct.c1, msg, eg1hs);
    /* c2 */
    element_pow_zn(ct.c2, sp.g, w[0]);
    /* c3 */
    element_pow_zn(ct.c3, gw.hS, ct.c5);
    /* c4 */
    H2(H2m, msg);
    element_pow2_zn(ct.c4, H2m, gw.t2, gw.hSt1, ct.c5);
    /* c6 */
    for (int i = 0; i < ROW; i++) {
        element_neg(neg, r[i]);

        element_pow2_zn(ct.c6[i][0], sp.ga, l[i], sp.F[ap.map[i]], neg);
        element_pow_zn(ct.c6[i][1], sp.g, r[i]);
    }
    /* c7 */
    element_set1(one);
    H3(T, ct.c1, ct.c2, ct.c3, ct.c4, ct.c6);
    element_pow3_zn(uT_vc5_w, sp.u, T, sp.v, ct.c5, sp.w, one);
    element_pow_zn(ct.c7, uT_vc5_w, w[0]);

    modifiedEncEnd = clock();

    ////

    modifiedEncTime += modifiedEncEnd - modifiedEncStart;

    element_clear(eg1hs);
    element_clear(H2m);
    element_clear(uT_vc5_w);
    for (int i = 0; i < COLUMN; i++)
        element_clear(w[i]);
    for (int i = 0; i < ROW; i++)
        element_clear(l[i]);
    for (int i = 0; i < ROW; i++)
        element_clear(r[i]);
    element_clear(mul);
    element_clear(neg);
    element_clear(one);
    element_clear(T);

    //std::cout << "finish" << std::endl;

    return ct;
}

void MEM_FUNC::modifiedDec(SP sp, MCT ct, MK mk, AP ap)
{
    //std::cout << "[Decryption]" << std::endl;
    //std::cout << "start" << std::endl;

    /* temp value */
    element_t uT_vc5_w;
    element_t ec2uT_vc5_w, ec7g,
              ec2Kd, eghas,
              ecL, edK, ecL_edK,
              mh, mul, prod;
    element_t T, one, zero;

    element_init_G1(uT_vc5_w, pairing);

    element_init_GT(ec2uT_vc5_w, pairing);
    element_init_GT(ec7g, pairing);
    element_init_GT(ec2Kd, pairing);
    element_init_GT(eghas, pairing);
    element_init_GT(ecL, pairing);
    element_init_GT(edK, pairing);
    element_init_GT(ecL_edK, pairing);
    element_init_GT(mh, pairing);
    element_init_GT(mul, pairing);
    element_init_GT(prod, pairing);

    element_init_Zr(T, pairing);
    element_init_Zr(one, pairing);
    element_init_Zr(zero, pairing);
    
    ////

    modifiedDecStart = clock();

    element_set1(one);
    H3(T, ct.c1, ct.c2, ct.c3, ct.c4, ct.c6);
    element_pow3_zn(uT_vc5_w, sp.u, T, sp.v, ct.c5, sp.w, one);
    element_pairing(ec2uT_vc5_w, ct.c2, uT_vc5_w);
    element_pairing(ec7g, ct.c7, sp.g);

    
    if (element_cmp(ec2uT_vc5_w, ec7g) == 1) {
        //std::cout << "CHECK FAIL" << std::endl;

        return;
    }

    element_pairing(ec2Kd, ct.c2, mk.Kd);
    element_set1(mul);
    for (int i = 0; i < ROW; i++) {
        element_pairing(ecL, ct.c6[i][0], mk.gd);
        if (ap.map[i] < MEM_ATT_LEN)
            element_pairing(edK, ct.c6[i][1], mk.fxd[ap.map[i]]);
        element_mul(ecL_edK, ecL, edK);
        if (i == 0)
            element_pow_zn(prod, ecL_edK, one);
        else
            element_pow_zn(prod, ecL_edK, zero);
        element_mul(mul, mul, prod);
    }
    element_div(eghas, ec2Kd, mul);

    element_div(mh, ct.c1, eghas);
    //std::cout << "DECRYPTION SUCCEED" << std::endl;

    modifiedDecEnd = clock();

    ////

    modifiedDecTime += modifiedDecEnd - modifiedDecStart;

    element_clear(ec2uT_vc5_w);
    element_clear(ec2Kd);
    element_clear(eghas);
    element_clear(ecL);
    element_clear(edK);
    element_clear(ecL_edK);
    element_clear(mh);
    element_clear(T);
    element_clear(one);
    element_clear(mul);
    element_clear(prod);
    element_clear(zero);

    //std::cout << "finish" << std::endl;
}

/* modified equality test */
GTD GA_FUNC::modifiedGroupTrap(SP sp, GSK gsk)
{
    GTD gtd;

    /* temp value */
    element_t z, zt1;

    element_init_Zr(z, pairing);
    element_init_Zr(zt1, pairing);

    ////

    element_random(z);
    element_pow_zn(gtd.gtd1, sp.ga, z);
    element_mul(zt1, z, gsk.t1);
    element_pow_zn(gtd.gtd2, sp.ga, zt1);

    ////

    element_clear(z);

    return gtd;
}

void ATP_FUNC::modifiedGroupTest(MCT ct1, MCT ct2, GTD gtd)
{
    /* temp value */
    element_t ec4gtd1, ec3gtd2,
              c4_inv, c3_inv;

    element_init_GT(ec4gtd1, pairing);
    element_init_GT(ec3gtd2, pairing);
    element_init_GT(c4_inv, pairing);
    element_init_GT(c3_inv, pairing);

    ////

    element_mul(c4_inv, ct1.c4, ct2.c4);
    element_pairing(ec4gtd1, c4_inv, gtd.gtd1);
    element_mul(c3_inv, ct1.c3, ct2.c3);
    element_pairing(ec3gtd2, c3_inv, gtd.gtd2);

    /*
    if (element_cmp(ec4gtd1, ec3gtd2) == 0)
        std::cout << "SAME" << std::endl;
    else
        std::cout << "DIFFERENT" << std::endl;
    */

    ////

    element_clear(ec4gtd1);
    element_clear(ec3gtd2);
    element_clear(c4_inv);
    element_clear(c3_inv);
}

/* Hash Functions */
void H1(element_s *result, int data[ROW])
{
    element_t tmpArr[ROW];
    int arrLen = 0, r1Len = element_length_in_bytes(r1);
    unsigned char *arrByte, r1Byte[r1Len], *inputByte;

    for (int i = 0; i < ROW; i++) {
        element_init_Zr(tmpArr[i], pairing);
        element_set_si(tmpArr[i], data[i]);

        arrLen += element_length_in_bytes(tmpArr[i]);
    }
    element_to_bytes(r1Byte, r1);

    arrByte = (unsigned char *)malloc(arrLen);

    int offset = 0;
    for (int i = 0; i < ROW; i++) {
        int len = element_length_in_bytes(tmpArr[i]);

        memcpy(arrByte + offset, tmpArr[i], len);
        offset += len;
    }

    inputByte = (unsigned char*)malloc(r1Len + arrLen);

    memcpy(inputByte, r1Byte, r1Len);
    memcpy(inputByte + r1Len, arrByte, arrLen);

    element_from_hash(result, inputByte, r1Len + arrLen);
}

void H2(element_s *result, element_t data)
{
    int arrLen = element_length_in_bytes(data);
    int r2Len = element_length_in_bytes(r2);
    unsigned char arrByte[arrLen], r2Byte[r2Len], *inputByte;

    element_to_bytes(arrByte, data);
    element_to_bytes(r2Byte, r2);
    inputByte = (unsigned char*)malloc(r2Len + arrLen);

    memcpy(inputByte, r2Byte, r2Len);
    memcpy(inputByte + r2Len, arrByte, arrLen);

    element_from_hash(result, inputByte, r2Len + arrLen);
}

void H3(element_s *result, element_t arr1, element_t arr2, element_t arr3, element_t arr4[ROW][2])
{
    int arr1Len = element_length_in_bytes(arr1),
        arr2Len = element_length_in_bytes(arr2),
        arr3Len = element_length_in_bytes(arr3),
        arr4Len = 0;
    unsigned char arr1Byte[arr1Len],
                  arr2Byte[arr2Len],
                  arr3Byte[arr3Len],
                  *arr4Byte,
                  *inputByte;
    
    for (int i = 0; i < ROW; i++)
        arr4Len += element_length_in_bytes(arr4[i][0]) + element_length_in_bytes(arr4[i][1]);
    arr4Byte = (unsigned char*)malloc(arr4Len);

    element_to_bytes(arr1Byte, arr1);
    element_to_bytes(arr2Byte, arr2);
    element_to_bytes(arr3Byte, arr3);

    int offset = 0;
    for (int i = 0; i < ROW; i++) {
        int arr4cLen = element_length_in_bytes(arr4[i][0]);
        int arr4dLen = element_length_in_bytes(arr4[i][1]);
        
        unsigned char arr4cByte[arr4cLen], arr4dByte[arr4dLen];

        element_to_bytes(arr4cByte, arr4[i][0]);
        element_to_bytes(arr4dByte, arr4[i][1]);

        memcpy(arr4Byte + offset, arr4cByte, arr4cLen);
        offset += arr4cLen;
        memcpy(arr4Byte + offset, arr4dByte, arr4dLen);
        offset += arr4dLen;
    }
    inputByte = (unsigned char*)malloc(arr1Len + arr2Len + arr3Len + arr4Len);

    memcpy(inputByte, arr1Byte, arr1Len);
    memcpy(inputByte + arr1Len, arr2Byte, arr2Len);
    memcpy(inputByte + arr2Len, arr3Byte, arr3Len);
    memcpy(inputByte + arr3Len, arr4Byte, arr4Len);

    element_from_hash(result, inputByte, arr1Len + arr2Len + arr3Len + arr4Len);
}

void H3(element_s *result, element_t arr1, element_t arr2, element_t arr3, element_t arr4, element_t arr5[ROW][2])
{
    int arr1Len = element_length_in_bytes(arr1),
        arr2Len = element_length_in_bytes(arr2),
        arr3Len = element_length_in_bytes(arr3),
        arr4Len = element_length_in_bytes(arr3),
        arr5Len = 0;
    unsigned char arr1Byte[arr1Len],
                  arr2Byte[arr2Len],
                  arr3Byte[arr3Len],
                  arr4Byte[arr4Len],
                  *arr5Byte,
                  *inputByte;
    
    for (int i = 0; i < ROW; i++)
        arr5Len += element_length_in_bytes(arr5[i][0]) + element_length_in_bytes(arr5[i][1]);
    arr5Byte = (unsigned char*)malloc(arr5Len);

    element_to_bytes(arr1Byte, arr1);
    element_to_bytes(arr2Byte, arr2);
    element_to_bytes(arr3Byte, arr3);
    element_to_bytes(arr4Byte, arr4);

    int offset = 0;
    for (int i = 0; i < ROW; i++) {
        int arr5cLen = element_length_in_bytes(arr5[i][0]);
        int arr5dLen = element_length_in_bytes(arr5[i][1]);
        
        unsigned char arr5cByte[arr5cLen], arr5dByte[arr5dLen];

        element_to_bytes(arr5cByte, arr5[i][0]);
        element_to_bytes(arr5dByte, arr5[i][1]);

        memcpy(arr5Byte + offset, arr5cByte, arr5cLen);
        offset += arr5cLen;
        memcpy(arr5Byte + offset, arr5dByte, arr5dLen);
        offset += arr5dLen;
    }
    inputByte = (unsigned char*)malloc(arr1Len + arr2Len + arr3Len + arr4Len + arr5Len);

    memcpy(inputByte, arr1Byte, arr1Len);
    memcpy(inputByte + arr1Len, arr2Byte, arr2Len);
    memcpy(inputByte + arr2Len, arr3Byte, arr3Len);
    memcpy(inputByte + arr3Len, arr4Byte, arr4Len);
    memcpy(inputByte + arr4Len, arr4Byte, arr5Len);

    element_from_hash(result, inputByte, arr1Len + arr2Len + arr3Len + arr4Len + arr5Len);
}