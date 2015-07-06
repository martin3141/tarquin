#include "basis_set.hpp"
#include "common.hpp"
#include <stdexcept>

namespace tarquin
{
	void getCrCH2(std::vector<std::vector<double> >& doubmat);
	void getCrCH2_RT(std::vector<std::vector<double> >& doubmat);
	void getAla(std::vector<std::vector<double> >& doubmat);
	void getAsp(std::vector<std::vector<double> >& doubmat);
	void getCho_RT(std::vector<std::vector<double> >& doubmat);
	void getCr(std::vector<std::vector<double> >& doubmat);
	void getCr_RT(std::vector<std::vector<double> >& doubmat);
	void getGABA(std::vector<std::vector<double> >& doubmat);
	void getMEGAGABA(std::vector<std::vector<double> >& doubmat);
	void getGPC(std::vector<std::vector<double> >& doubmat);
	void getGlc(std::vector<std::vector<double> >& doubmat);
	void getGln(std::vector<std::vector<double> >& doubmat);
	void getGlu(std::vector<std::vector<double> >& doubmat);
	void getGlu_RT(std::vector<std::vector<double> >& doubmat);
	void getGua(std::vector<std::vector<double> >& doubmat);
	void getIns_RT(std::vector<std::vector<double> >& doubmat);
	void getIns(std::vector<std::vector<double> >& doubmat);
	void getLac_RT(std::vector<std::vector<double> >& doubmat);
	void getLac(std::vector<std::vector<double> >& doubmat);
	void getLip09(std::vector<std::vector<double> >& doubmat, double B0);
	void getLip13a(std::vector<std::vector<double> >& doubmat, double B0);
	void getLip13b(std::vector<std::vector<double> >& doubmat, double B0);
	void getLip20(std::vector<std::vector<double> >& doubmat, double B0);
	void getMM09(std::vector<std::vector<double> >& doubmat, double B0);
	void getMM12(std::vector<std::vector<double> >& doubmat, double B0);
	void getMM14(std::vector<std::vector<double> >& doubmat, double B0);
	void getMM17(std::vector<std::vector<double> >& doubmat, double B0);
	void getMM20(std::vector<std::vector<double> >& doubmat, double B0);
	void getMM30(std::vector<std::vector<double> >& doubmat, double B0);
	void getMM38(std::vector<std::vector<double> >& doubmat, double B0);
	void getNAA(std::vector<std::vector<double> >& doubmat);
	void getNAA_RT(std::vector<std::vector<double> >& doubmat);
	void getNAAG(std::vector<std::vector<double> >& doubmat);
	void getPCh(std::vector<std::vector<double> >& doubmat);
	void getScyllo(std::vector<std::vector<double> >& doubmat);
	void getTau(std::vector<std::vector<double> >& doubmat);
	void getGly(std::vector<std::vector<double> >& doubmat);
	void getCit(std::vector<std::vector<double> >& doubmat);
	void getGlth(std::vector<std::vector<double> >& doubmat);
	void getPCr(std::vector<std::vector<double> >& doubmat);
	void getPEth(std::vector<std::vector<double> >& doubmat);
	
    void getTNAA(std::vector<std::vector<double> >& doubmat);
    void getTCho(std::vector<std::vector<double> >& doubmat);
    
    // MEGA PRESS metabs
	void getMM09_MEGA(std::vector<std::vector<double> >& doubmat, double B0);
	void getGABA_A(std::vector<std::vector<double> >& doubmat);
	void getGABA_B(std::vector<std::vector<double> >& doubmat);
	void getNAA_MEGA(std::vector<std::vector<double> >& doubmat);
	void getGlx_A(std::vector<std::vector<double> >& doubmat);
	void getGlx_B(std::vector<std::vector<double> >& doubmat);
	void getGlx_C(std::vector<std::vector<double> >& doubmat);
	void getGlx_D(std::vector<std::vector<double> >& doubmat);

    // 31P metabs
    void get31P_ATP(std::vector<std::vector<double> >& doubmat);
    void get31P_GPC(std::vector<std::vector<double> >& doubmat);
    void get31P_GPE(std::vector<std::vector<double> >& doubmat);
    void get31P_NADH(std::vector<std::vector<double> >& doubmat);
    void get31P_PCH(std::vector<std::vector<double> >& doubmat);
    void get31P_PCR(std::vector<std::vector<double> >& doubmat);
    void get31P_PE(std::vector<std::vector<double> >& doubmat);
    void get31P_PI(std::vector<std::vector<double> >& doubmat);
}

void tarquin::getMetaboliteDescription(basis_vector_e metab, std::string& desc)
{
	switch( metab )
	{
		case BV_CRCH2:
			desc = "-CrCH2.csv";
			break;

		case BV_CRCH2_RT:
			desc = "-CrCH2.csv";
			break;

		case BV_ALA:
			desc = "Ala.csv";
			break;

		case BV_ASP:
			desc = "Asp.csv";
			break;
		
        case BV_CHO_RT:
            desc = "Cho.csv";
            break;

		case BV_CR: 
            desc = "Cr.csv";
            break;
		
        case BV_CR_RT:
			desc = "Cr.csv";
			break;
        
        case BV_PCR:
			desc = "PCr.csv";
			break;
        
        case BV_PETH:
			desc = "PEth.csv";
			break;

		case BV_GABA:
			desc = "GABA.csv";
			break;

		case BV_MEGAGABA:
			desc = "GABA.csv";
			break;

		case BV_GPC:
			desc = "GPC.csv";
			break;

		case BV_GLC:
			desc = "Glc.csv";
			break;

		case BV_GLN:
			desc = "Gln.csv";
			break;

        case BV_GLTH:
			desc = "Glth.csv";
			break;

		case BV_GLU_RT:
			desc = "Glu.csv";
			break;

		case BV_GLU:
			desc = "Glu.csv";
			break;

		case BV_GUA:
			desc = "Gua.csv";
			break;

		case BV_INS:
			desc = "Ins.csv";
			break;
		
        case BV_INS_RT:
			desc = "Ins.csv";
			break;

		case BV_LAC:
			desc = "Lac.csv";
			break;

		case BV_LAC_RT:
			desc = "Lac.csv";
			break;

		case BV_LIP09:
			desc = "Lip09.csv";
			break;

		case BV_LIP13A:
			desc = "Lip13a.csv";
			break;

		case BV_LIP13B:
			desc = "Lip13b.csv";
			break;

		case BV_LIP20:
			desc = "Lip20.csv";
			break;

		case BV_MM09:
			desc = "MM09.csv";
			break;

		case BV_MM12:
			desc = "MM12.csv";
			break;

		case BV_MM14:
			desc = "MM14.csv";
			break;

		case BV_MM17:
			desc = "MM17.csv";
			break;

		case BV_MM20:
			desc = "MM20.csv";
			break;
        
        case BV_MM30:
			desc = "MM30.csv";
			break;
        
        case BV_MM38:
			desc = "MM38.csv";
			break;

		case BV_NAA:
			desc = "NAA.csv";
			break;

		case BV_NAA_RT:
			desc = "NAA.csv";
			break;

		case BV_NAAG:
			desc = "NAAG.csv";
			break;

		case BV_PCH:
			desc = "PCh.csv";
			break;

		case BV_SCYLLO:
			desc = "Scyllo.csv";
			break;

		case BV_TAU:
			desc = "Tau.csv";
			break;
        
        case BV_GLY:
			desc = "Gly.csv";
			break;

        case BV_CIT:
			desc = "Cit.csv";
			break;
        
        case BV_TCHO:
			desc = "TCho.csv";
			break;

        case BV_TNAA:
			desc = "TNAA.csv";
			break;

        case BV_GABA_A:
            desc = "GABA_A.csv";
			break;

        case BV_GABA_B:
            desc = "GABA_B.csv";
			break;

        case BV_NAA_MEGA:
            desc = "NAA.csv";
			break;

        case BV_GLX_A:
            desc = "Glx_A.csv";
			break;

        case BV_GLX_B:
            desc = "Glx_B.csv";
			break;

        case BV_GLX_C:
            desc = "Glx_C.csv";
			break;

        case BV_GLX_D:
            desc = "Glx_D.csv";
			break;
        
        case BV_MM09_MEGA:
			desc = "MM09.csv";
			break;
    
        case BV_31P_ATP:
            desc = "APT.csv";
            break;

        case BV_31P_GPC:
            desc = "GPC.csv";
            break;

        case BV_31P_GPE:
            desc = "GPE.csv";
            break;

        case BV_31P_NADH:
            desc = "NADH.csv";
            break;

        case BV_31P_PCH:
            desc = "PCH.csv";
            break;

        case BV_31P_PCR:
            desc = "PCr.csv";
            break;

        case BV_31P_PE:
            desc = "PE.csv";
            break;

        case BV_31P_PI:
            desc = "PI.csv";
            break;
		
        default:
			throw std::runtime_error("unknown basis vector enumeration value");
	}
}


void tarquin::getMetaboliteMatrix(
		basis_vector_e metab,
		std::vector<std::vector<double> >& doubmat,
		double tf
		)
{
	switch( metab )
	{
		case BV_CRCH2:
			getCrCH2(doubmat);
			break;
            
		case BV_CRCH2_RT:
			getCrCH2_RT(doubmat);
			break;

		case BV_ALA:
			getAla(doubmat);
			break;

		case BV_ASP:
			getAsp(doubmat);
			break;

		case BV_CHO_RT:
			getCho_RT(doubmat);
			break;

		case BV_CR:
			getCr(doubmat);
			break;

		case BV_CR_RT:
			getCr_RT(doubmat);
			break;

        case BV_PCR:
			getPCr(doubmat);
			break;

        case BV_PETH:
			getPEth(doubmat);
			break;

		case BV_GABA:
			getGABA(doubmat);
			break;
		
        case BV_MEGAGABA:
			getMEGAGABA(doubmat);
			break;

		case BV_GPC:
			getGPC(doubmat);
			break;

		case BV_GLC:
			getGlc(doubmat);
			break;

		case BV_GLN:
			getGln(doubmat);
			break;
        
        case BV_GLTH:
			getGlth(doubmat);
			break;

		case BV_GLU:
			getGlu(doubmat);
			break;

		case BV_GLU_RT:
			getGlu_RT(doubmat);
			break;

		case BV_GUA:
			getGua(doubmat);
			break;

		case BV_INS:
			getIns(doubmat);
			break;
		
        case BV_INS_RT:
			getIns_RT(doubmat);
			break;

		case BV_LAC:
			getLac(doubmat);
			break;

		case BV_LAC_RT:
			getLac_RT(doubmat);
			break;

		case BV_LIP09:
			getLip09(doubmat, tf);
			break;

		case BV_LIP13A:
			getLip13a(doubmat, tf);
			break;

		case BV_LIP13B:
			getLip13b(doubmat, tf);
			break;

		case BV_LIP20:
			getLip20(doubmat, tf);
			break;

		case BV_MM09:
			getMM09(doubmat, tf);
			break;

		case BV_MM12:
			getMM12(doubmat, tf);
			break;

		case BV_MM14:
			getMM14(doubmat, tf);
			break;

		case BV_MM17:
			getMM17(doubmat, tf);
			break;

		case BV_MM20:
			getMM20(doubmat, tf);
			break;
        
        case BV_MM30:
			getMM30(doubmat, tf);
			break;
        
        case BV_MM38:
			getMM38(doubmat, tf);
			break;

		case BV_NAA:
			getNAA(doubmat);
			break;

		case BV_NAA_RT:
			getNAA_RT(doubmat);
			break;

		case BV_NAAG:
			getNAAG(doubmat);
			break;

		case BV_PCH:
			getPCh(doubmat);
			break;

		case BV_SCYLLO:
			getScyllo(doubmat);
			break;

		case BV_TAU:
			getTau(doubmat);
			break;

		case BV_GLY:
			getGly(doubmat);
			break;
       
        case BV_TCHO:
			getTCho(doubmat);
			break;

        case BV_TNAA:
			getTNAA(doubmat);
			break;

		case BV_CIT:
			getCit(doubmat);
			break;

        case BV_GABA_A:
            getGABA_A(doubmat);
			break;

        case BV_GABA_B:
            getGABA_B(doubmat);
			break;

        case BV_NAA_MEGA:
            getNAA_MEGA(doubmat);
			break;

        case BV_GLX_A:
            getGlx_A(doubmat);
			break;

        case BV_GLX_B:
            getGlx_B(doubmat);
			break;

        case BV_GLX_C:
            getGlx_C(doubmat);
			break;
        
        case BV_GLX_D:
            getGlx_D(doubmat);
			break;
    
        case BV_MM09_MEGA:
			getMM09_MEGA(doubmat, tf);
			break;
        
        case BV_31P_ATP:
			get31P_ATP(doubmat);
			break;

        case BV_31P_GPC:
			get31P_GPC(doubmat);
			break;

        case BV_31P_GPE:
			get31P_GPE(doubmat);
			break;

        case BV_31P_NADH:
			get31P_NADH(doubmat);
			break;

        case BV_31P_PCH:
			get31P_PCH(doubmat);
			break;

        case BV_31P_PCR:
			get31P_PCR(doubmat);
			break;

        case BV_31P_PE:
			get31P_PE(doubmat);
			break;

        case BV_31P_PI:
			get31P_PI(doubmat);
			break;

		default:
			throw std::runtime_error("unknown basis vector enumeration value");
	}
}

void tarquin::getMM38(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	doubvec[5] = -pow(0.26 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	doubvec[5] = -pow(0.15 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	doubvec[5] = -pow(0.15 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);

	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1.7;
	doubvec[2] = 0.5;
	doubvec[3] = 3.85;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1.2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.82;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 3;
	doubvec[1] = 0.9;
	doubvec[2] = 0.5;
	doubvec[3] = 3.77;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}


void tarquin::getMM20(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(8, d_nan);
	std::vector<double> doubvec(8, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1307;
	doubvec[5] = -pow(0.15 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 2325;
	doubvec[5] = -pow(0.2 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1.33;
	doubvec[2] = 0.5;
	doubvec[3] = 2.08;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 0.33;
	doubvec[2] = 0.5;
	doubvec[3] = 2.25;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 0.33;
	doubvec[2] = 0.5;
	doubvec[3] = 1.95;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 0.4;
	doubvec[2] = 0.5;
	doubvec[3] = 3;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getMM17(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1307;
	doubvec[5] = -pow(0.15 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 1.67;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getMM14(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1679;
	doubvec[5] = -pow(0.17 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 1.43;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getMM12(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1307;
	doubvec[5] = -pow(0.15 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 1.21;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getMM09(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1139;
	doubvec[5] = -pow(0.14 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 0.91;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getMM30(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1139;
	doubvec[5] = -pow(0.14 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 3.0;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getLip20(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(7, d_nan);
	std::vector<double> doubvec(7, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1307;
	doubvec[5] = -pow(0.15 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 2325;
	doubvec[5] = -pow(0.2 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1.33;
	doubvec[2] = 0.5;
	doubvec[3] = 2.04;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 0.67;
	doubvec[2] = 0.5;
	doubvec[3] = 2.25;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 0.87;
	doubvec[2] = 0.5;
	doubvec[3] = 2.8;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getLip13b(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 460.4;
	doubvec[5] = -pow(0.089 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 1.28;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getLip13a(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1307;
	doubvec[5] = -pow(0.15 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 1.28;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getLip09(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1139;
	doubvec[5] = -pow(0.14 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 0.89;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getTau(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(8, d_nan);
	std::vector<double> doubvec(8, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.4206;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.4206;
	doubvec[4] = -12.438;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.2459;
	doubvec[4] = 6.742;
	doubvec[5] = 6.403;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.2459;
	doubvec[4] = 6.464;
	doubvec[5] = 6.792;
	doubvec[6] = -12.93;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getScyllo(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 6;
	doubvec[2] = 0.5;
	doubvec[3] = 3.34;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGly(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.548;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getCit(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 2.6735;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
    doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 2.5265;
	doubvec[4] = -16.1;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getNAAG(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 2.042;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getLac(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(8, d_nan);
	std::vector<double> doubvec(8, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.0974;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.3142;
	doubvec[4] = 6.933;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.3142;
	doubvec[4] = 6.933;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.3142;
	doubvec[4] = 6.933;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getLac_RT(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(8, d_nan);
	std::vector<double> doubvec(8, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.0974;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.3142;
	doubvec[4] = 6.7;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.3142;
	doubvec[4] = 6.7;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.3142;
	doubvec[4] = 6.7;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
}


void tarquin::getIns(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.5217;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.0538;
	doubvec[4] = 2.889;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.5217;
	doubvec[4] = 0;
	doubvec[5] = 3.006;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.6144;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 9.997;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.269;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 9.485;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.6144;
	doubvec[4] = 9.998;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 9.482;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getIns_RT(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.5217;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.0538;
	doubvec[4] = 2.889;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.5217;
	doubvec[4] = 0;
	doubvec[5] = 3.006;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.6144;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 9.997;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.269;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 9.485;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.6144;
	doubvec[4] = 9.998;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 9.482;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGua(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.78;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGlu(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(9, d_nan);
	std::vector<double> doubvec(9, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.7433;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.0375;
	doubvec[4] = 7.331;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.12;
	doubvec[4] = 4.651;
	doubvec[5] = -14.849;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.3378;
	doubvec[4] = 0;
	doubvec[5] = 6.413;
	doubvec[6] = 8.478;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.352;
	doubvec[4] = 0;
	doubvec[5] = 8.406;
	doubvec[6] = 6.875;
	doubvec[7] = -15.915;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGlu_RT(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(9, d_nan);
	std::vector<double> doubvec(9, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.7433;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.0375;
	doubvec[4] = 7.331;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.12;
	doubvec[4] = 4.651;
	doubvec[5] = -14.849;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.3378;
	doubvec[4] = 0;
	doubvec[5] = 6.413;
	doubvec[6] = 8.478;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.352;
	doubvec[4] = 0;
	doubvec[5] = 8.406;
	doubvec[6] = 6.875;
	doubvec[7] = -15.915;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGlth(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(9, d_nan);
	std::vector<double> doubvec(9, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 4;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 5;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
    doubvec[0] = 6;
	doubvec[1] = 6;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.769;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
    doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.5608;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.9264;
	doubvec[4] = 7.09;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.9747;
	doubvec[4] = 4.71;
	doubvec[5] = -14.06;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
    doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.769;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.159;
	doubvec[4] = 6.34;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.146;
	doubvec[4] = 6.36;
	doubvec[5] = -15.48;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.51;
	doubvec[4] = 0;
	doubvec[5] = 6.7;
	doubvec[6] = 7.6;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.56;
	doubvec[4] = 0;
	doubvec[5] = 7.6;
	doubvec[6] = 6.7;
	doubvec[7] = -15.92;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGln(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(9, d_nan);
	std::vector<double> doubvec(9, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.753;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.129;
	doubvec[4] = 5.847;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.109; //corrected in v 4.3.5 from 2.19
	doubvec[4] = 6.5;
	doubvec[5] = -14.504;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.432;
	doubvec[4] = 0;
	doubvec[5] = 9.165;
	doubvec[6] = 6.324;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.454;
	doubvec[4] = 0;
	doubvec[5] = 6.347;
	doubvec[6] = 9.209;
	doubvec[7] = -15.371;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
}


void tarquin::getGlc(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(11, d_nan);
	std::vector<double> doubvec(11, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 7;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 5.216;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.519;
	doubvec[4] = 3.8;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.698;
	doubvec[4] = 0;
	doubvec[5] = 9.6;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.395;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 9.4;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.822;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 9.9;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.826;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 1.5;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 7;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.749;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 6;
	doubvec[9] = -12.1;
	doubvec[10] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGABA(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.0128;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.0128;
	doubvec[4] = -12.021;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.889;
	doubvec[4] = 5.372;
	doubvec[5] = 10.578;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.889;
	doubvec[4] = 7.127;
	doubvec[5] = 6.982;
	doubvec[6] = -13.121;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.284;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 7.755;
	doubvec[7] = 6.173;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.284;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 7.432;
	doubvec[7] = 7.933;
	doubvec[8] = -10.744;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getMEGAGABA(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.0128;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.0128;
	doubvec[4] = -12.021;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 1.889;
	doubvec[4] = 5.372;
	doubvec[5] = 10.578;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 1.889;
	doubvec[4] = 7.127;
	doubvec[5] = 6.982;
	doubvec[6] = -13.121;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 2.284;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 7.755;
	doubvec[7] = 6.173;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 2.284;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 7.432;
	doubvec[7] = 7.933;
	doubvec[8] = -10.744;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getAsp(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(7, d_nan);
	std::vector<double> doubvec(7, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.8914;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.8011;
	doubvec[4] = 3.647;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.6533;
	doubvec[4] = 9.107;
	doubvec[5] = -17.426;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getCrCH2(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = -2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.913;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getCrCH2_RT(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = -2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.922;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}


void tarquin::getAla(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(8, d_nan);
	std::vector<double> doubvec(8, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.7746;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.4667;
	doubvec[4] = 7.234;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.4667;
	doubvec[4] = 7.234;
	doubvec[5] = -14.366;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 1.4667;
	doubvec[4] = 7.234;
	doubvec[5] = -14.366;
	doubvec[6] = -14.366;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
}




void tarquin::getGPC(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 4;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 7;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 8;
	doubvec[1] = 4;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 9;
	doubvec[1] = 4;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 9;
	doubvec[2] = 0.5;
	doubvec[3] = 3.212;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.605;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.672;
	doubvec[4] = -14.78;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.903;
	doubvec[4] = 5.77;
	doubvec[5] = 4.53;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.871;
	doubvec[4] = 0;
	doubvec[5] = 5.77;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.946;
	doubvec[4] = 0;
	doubvec[5] = 4.53;
	doubvec[6] = 0;
	doubvec[7] = -14.78;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubvec[6] = 0;
	doubvec[7] = 6.03;
	doubvec[8] = 6.03;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 6;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.312;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.312;
	doubvec[4] = -9.32;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 7;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.659;
	doubvec[4] = 3.1;
	doubvec[5] = 5.9;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 7;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.659;
	doubvec[4] = 5.9;
	doubvec[5] = 3.1;
	doubvec[6] = -9.32;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 8;
	doubvec[1] = 0;
	doubvec[2] = 1;
	doubvec[3] = 0;
	doubvec[4] = 2.67;
	doubvec[5] = 2.67;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 9;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 6.03;
	doubvec[5] = 6.03;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 0;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getPCh(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 7;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 9;
	doubvec[2] = 0.5;
	doubvec[3] = 3.208;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.2805;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.2805;
	doubvec[4] = -14.89;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.641;
	doubvec[4] = 2.284;
	doubvec[5] = 7.326;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 5;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.641;
	doubvec[4] = 7.231;
	doubvec[5] = 2.235;
	doubvec[6] = -14.19;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 6;
	doubvec[1] = 0;
	doubvec[2] = 1;
	doubvec[3] = 0;
	doubvec[4] = 2.68;
	doubvec[5] = 2.772;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 7;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 6.298;
	doubvec[5] = 6.249;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 0;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getCho_RT(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 9;
	doubvec[2] = 0.5;
	doubvec[3] = 3.185;
	doubvec[4] = 0;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 0.57;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 4.054;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 4.054;
	doubvec[4] = -14.1;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.501;
	doubvec[4] = 3.14;
	doubvec[5] = 7.011;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.501;
	doubvec[4] = 6.979;
	doubvec[5] = 3.168;
	doubvec[6] = -14.07;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 2.572;
	doubvec[5] = 2.681;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getCr(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 3.027;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.913;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getCr_RT(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 3.027;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.922;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}


void tarquin::getPCr(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 3.029;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0.5;
	doubvec[3] = 3.930;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getNAA(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(7, d_nan);
	std::vector<double> doubvec(7, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 2.008;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.3817;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.6727;
	doubvec[4] = 3.861;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.4863;
	doubvec[4] = 9.821;
	doubvec[5] = -15.592;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);

	// print the doubmat
	/*for (vector< vector<double> >::size_type u = 0; u < naa_doubmat.size(); u++) {
		for (vector<double>::size_type v = 0; v < naa_doubmat[u].size(); v++) {
			cout << naa_doubmat[u][v] << " ";
		}
		cout << endl;
	}*/
}

void tarquin::getNAA_RT(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(7, d_nan);
	std::vector<double> doubvec(7, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 3;
	doubvec[2] = 0;
	doubvec[4] = 0.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 2.008;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 4.3817;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 3;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.6827;
	doubvec[4] = 3.65;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 4;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.4863;
	doubvec[4] = 9.821;
	doubvec[5] = -15.592;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);

	// print the doubmat
	/*for (vector< vector<double> >::size_type u = 0; u < naa_doubmat.size(); u++) {
		for (vector<double>::size_type v = 0; v < naa_doubmat[u].size(); v++) {
			cout << naa_doubmat[u][v] << " ";
		}
		cout << endl;
	}*/
}


void tarquin::getTNAA(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 2.008;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getTCho(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 9;
	doubvec[2] = 0.5;
	doubvec[3] = 3.212;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getGABA_A(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 30;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 0.4;
	doubvec[2] = 0.5;
	doubvec[3] = 3.04;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::getGABA_B(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 20;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 0.4;
	doubvec[2] = 0.5;
	doubvec[3] = 2.95;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::getNAA_MEGA(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = -3;
	doubvec[2] = 0.5;
	doubvec[3] = 2.0;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::getGlx_A(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.299;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::getGlx_B(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.400;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::getGlx_C(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.707;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::getGlx_D(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.789;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getMM09_MEGA(std::vector<std::vector<double> >& doubmat, double B0)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 0;
	//doubvec[5] = 1139;
	doubvec[5] = -pow(0.07 * B0 * 1.0e-6 * M_PI / 2.0, 2.0) / log(0.5);
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 3;
	doubvec[2] = 0.5;
	doubvec[3] = 0.91;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_ATP(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(7, d_nan);
	std::vector<double> doubvec(7, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = -7.616;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
    doubvec[0] = 1;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = -16.26;
	doubvec[4] = 16.3;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
    doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = -2.6;
	doubvec[4] = 0;
	doubvec[5] = 16.1;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_GPC(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 2.93;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_GPE(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.5;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_NADH(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = -8.3;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_PCH(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 6.243;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_PCR(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_PE(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 6.7222;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}
void tarquin::get31P_PI(std::vector<std::vector<double> >& doubmat)
{
    double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(6, d_nan);
	std::vector<double> doubvec(6, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 4.8161;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
}

void tarquin::getPEth(std::vector<std::vector<double> >& doubmat)
{
	double d_nan = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> nan_doubvec(10, d_nan);
	std::vector<double> doubvec(10, d_nan);
	doubmat.push_back(nan_doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 2;
	doubvec[2] = 0;
	doubvec[4] = 2.5;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubmat.push_back(nan_doubvec);
	doubvec = nan_doubvec;
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.9765;
	doubvec[4] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 1;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.9765;
	doubvec[4] = -14.56;
	doubvec[5] = 0;
	doubmat.push_back(doubvec);
	doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.216;
	doubvec[4] = 3.182;
	doubvec[5] = 7.204;
	doubvec[6] = 0;
	doubmat.push_back(doubvec);
    doubvec[0] = 2;
	doubvec[1] = 1;
	doubvec[2] = 0.5;
	doubvec[3] = 3.216;
	doubvec[4] = 6.716;
	doubvec[5] = 2.980;
	doubvec[6] = -14.71;
	doubvec[7] = 0;
	doubmat.push_back(doubvec);
    doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 7.228;
	doubvec[5] = 7.088;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 0;
	doubmat.push_back(doubvec);
    doubvec[0] = 2;
	doubvec[1] = 0;
	doubvec[2] = 0.5;
	doubvec[3] = 0;
	doubvec[4] = 0.464;
	doubvec[5] = 0.588;
	doubvec[6] = 0;
	doubvec[7] = 0;
	doubvec[8] = 0;
	doubvec[9] = 0;
	doubmat.push_back(doubvec);


}
