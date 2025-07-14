#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "pbc/pbc.h"
#include <iostream>

#define SECURITY_PARAM_BIT 128
// size of attribute universe
#define SU 10
// number of iterations
#define ITR 2
// number of members
#define MEM_NUM 5
// size of member's attribute set
#define MEM_ATT_LEN 5
#define ROW 5
#define COLUMN 5

/* pairing */
extern pairing_t pairing;

typedef int securityParam;
typedef int attUniverseNum;
typedef int memberAttSet;
typedef element_t MSG;
typedef double test;

// extract time
extern test extStart, extEnd, extTime;
// encryption time
extern test originEncStart, originEncEnd, originEncTime;
extern test modifiedEncStart, modifiedEncEnd, modifiedEncTime;
// decryption time
extern test originDecStart, originDecEnd, originDecTime;
extern test modifiedDecStart, modifiedDecEnd, modifiedDecTime;

typedef struct systemParam
{
    /* G element */
    element_t g, ga,
              eg1h, eg2h,
              u, v, w,
              F[SU];
} SP;

typedef struct masterSecretKey
{
    /* G element */
    element_t ha, hb,
              k1, k2, k3;
} MSK;

typedef struct setup
{
    SP sp;
    MSK msk;
} SET;

typedef struct groupSecretKey
{
    /* Zp element */
    element_t t1, t2;
} GSK;

typedef struct groupWarrantee
{
    /* G element */
    element_t hb_gat1, gt1;
    element_t fxt1[MEM_ATT_LEN];

    /* GT element */
    element_t hS, hSt1;

    /* Zp element */
    element_t t2;
} GW;

typedef struct memberKey
{
    /* G element */
    element_t Kd, Kt,
              gd, gt;
    element_t fxd[MEM_ATT_LEN], fxt[MEM_ATT_LEN];
} MK;

typedef struct groupTrapDoor
{
    /* G element */
    element_t gtd1, gtd2;
} GTD;

typedef struct originalCiphertext
{
    /* G element */
    element_t c3, c8;
    element_t c7[ROW][2];   // (c_i, d_i)

    /* GT element */
    element_t c1, c2, c4, c5;

    /* Zp element */
    element_t c6;
} OCT;

typedef struct modifiedCiphertext
{
    /* G element */
    element_t c2, c7;
    element_t c6[ROW][2];

    /* GT element */
    element_t c1, c3, c4;

    /* Zp element */
    element_t c5;
} MCT;

typedef struct accessPolicy
{
    element_t M[ROW][COLUMN];
    int map[ROW];
} AP;

/* Parties */
typedef struct groupAdmin
{
    GSK gsk;
    MSK msk;
} GA;

typedef struct member
{
    memberAttSet S[MEM_ATT_LEN];
    GW gw;
    MK mk;
    AP ap;
    MSG msg;
} MEM;

typedef struct authorizedThirdParty
{
    GTD gtd;

    OCT oct[MEM_NUM];
    MCT mct[MEM_NUM];
} ATP;

/* Functions */
SET InitialSetup(securityParam p, attUniverseNum u);

typedef class systemFunction
{
    public:
    GSK GroupKeyGen();
} SYS_FUNC;

typedef class groupAdminFunction
{
    public:
    GW Join(SP, MSK, GSK, memberAttSet[MEM_ATT_LEN]);
    MK Extract(SP, MSK, memberAttSet[MEM_ATT_LEN]);
    GTD originalGroupTrap(SP, GSK);
    GTD modifiedGroupTrap(SP, GSK);
} GA_FUNC;

typedef class memberFunction
{
    public:
    OCT originalEnc(SP, GW, AP, MSG);
    void originalDec(SP, OCT, MK, AP);

    MCT modifiedEnc(SP, GW, AP, MSG);
    void modifiedDec(SP, MCT, MK, AP);
} MEM_FUNC;

typedef class authorizedThirdPartyFunction
{
    public:
    void originalGroupTest(OCT, OCT, GTD);
    void modifiedGroupTest(MCT, MCT, GTD);
} ATP_FUNC;

/* Hashes */
void H1(element_s *, int[ROW]);
void H2(element_s *, element_t);
void H3(element_s *, element_t, element_t, element_t, element_t[ROW][2]);
void H3(element_s *, element_t, element_t, element_t, element_t, element_t[ROW][2]);

/* clear */
//void clearElements();

#endif