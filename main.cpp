#include "structures.h"
#include <random>
#include <algorithm>

SET set;

GA ga;
ATP atp;
MEM mem[MEM_NUM];

SYS_FUNC sysFunc;
GA_FUNC gaFunc;
MEM_FUNC memFunc;
ATP_FUNC atpFunc;

// extract time
test extStart, extEnd, extTime = 0;
// encryption time
test originEncStart, originEncEnd, originEncTime = 0;
test modifiedEncStart, modifiedEncEnd, modifiedEncTime = 0;
// decryption time
test originDecStart, originDecEnd, originDecTime = 0;
test modifiedDecStart, modifiedDecEnd, modifiedDecTime = 0;

std::random_device rd;
std::mt19937 mt(rd());

void setMembers();

int main()
{
    for (int i = 0; i < ITR; i++) {
        set = InitialSetup(SECURITY_PARAM_BIT, SU);

        ga.gsk = sysFunc.GroupKeyGen();

        setMembers();

        /* original */
        /* Encryption and Decryption*/
        for (int i = 0; i < MEM_NUM; i++)
            atp.oct[i] = memFunc.originalEnc(set.sp, mem[i].gw, mem[i].ap, mem[i].msg);

        for (int i = 0; i < MEM_NUM; i++) {
            for (int j = 0; j < MEM_NUM; j++) {
                if (i != j)
                    memFunc.originalDec(set.sp, atp.oct[i], mem[j].mk, mem[i].ap);
            }
        }

        /* Equality Test*/
        atp.gtd = gaFunc.originalGroupTrap(set.sp, ga.gsk);
        for (int i = 0; i < MEM_NUM; i++)
            for (int j = 0; j < i; j++)
                atpFunc.originalGroupTest(atp.oct[i], atp.oct[j], atp.gtd);
        
        /* modified */
        /* Encryption and Decryption*/
        for (int i = 0; i < MEM_NUM; i++)
            atp.mct[i] = memFunc.modifiedEnc(set.sp, mem[i].gw, mem[i].ap, mem[i].msg);

        for (int i = 0; i < MEM_NUM; i++) {
            for (int j = 0; j < MEM_NUM; j++) {
                if (i != j)
                    memFunc.modifiedDec(set.sp, atp.mct[i], mem[j].mk, mem[i].ap);
            }
        }

        /* Equality Test*/
        atp.gtd = gaFunc.modifiedGroupTrap(set.sp, ga.gsk);
        for (int i = 0; i < MEM_NUM; i++)
            for (int j = 0; j < i; j++)
                atpFunc.modifiedGroupTest(atp.mct[i], atp.mct[j], atp.gtd);
    }
    std::cout << "[SETTING]" << std::endl;
    std::cout << "SU: " << SU << std::endl;
    std::cout << "l:  " << ROW << std::endl;
    std::cout << std::endl;

    std::cout << "[Extract Time]" << std::endl;
    std::cout << ": " << extTime / (ITR * MEM_NUM) << std::endl;
    std::cout << std::endl;
    
    std::cout << "[Encryption Time]" << std::endl;
    std::cout << "Original: " << originEncTime / (ITR * MEM_NUM) << std::endl;
    std::cout << "Modified: " << modifiedEncTime / (ITR * MEM_NUM) << std::endl;
    std::cout << std::endl;

    std::cout << "[Decryption Time]" << std::endl;
    std::cout << "Original: " << originDecTime / (ITR * MEM_NUM * (MEM_NUM - 1)) << std::endl;
    std::cout << "Modified: " << modifiedDecTime / (ITR * MEM_NUM * (MEM_NUM - 1)) << std::endl;
    std::cout << std::endl;

    return 0;
}

void setMembers()
{
        for (MEM &m : mem) {
        std::uniform_int_distribution<int> distA(0, SU - 1);
        for (int i = 0; i < ROW; i++) {
            m.ap.map[i] = distA(mt);

            for (int j = 0; j < i; j++)
                if (m.ap.map[i] == m.ap.map[j]) {
                    i--;

                    break;
                }

            if (i < MEM_ATT_LEN)
                m.S[i] = m.ap.map[i];
        }

        std::sort(m.S, m.S + MEM_ATT_LEN);

        for (int i = 0; i < ROW; i++)
            for (int j = 0; j < ROW; j++) {
                element_init_Zr(m.ap.M[i][j], pairing);
                if (i == j)
                    element_set1(m.ap.M[i][j]);
                else
                    element_set0(m.ap.M[i][j]);
            }
        
        element_init_GT(m.msg, pairing);
        element_random(m.msg);

        m.gw = gaFunc.Join(set.sp, set.msk, ga.gsk, m.S);
        m.mk = gaFunc.Extract(set.sp, set.msk, m.S);
    }
}