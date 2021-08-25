#include "read_pp.h"

//  read bulk local pseudopotential "blps" in the Unified
// int Pseudopot_upf::read_pseudo_blps(ifstream &ifs)
// {
//     if(!SCAN_BEGIN(ifs,"START")) WARNING_QUIT("read_pp_blps","Find no PP_HEADER");
//     skip_comment(ifs);

//     float max_g;
//     ifs >> max_g;

//     char comment[200];
//     int nline = 0;
//     while(ifs.getline(comment,sizeof(comment)))
//     {
//         if(comment != "1000"){
//             ++nline;
//         }
//         else{
//             break;
//         }
//     }
//     this->mesh = nline * 3;
//     double dg = max_g / (mesh - 1);

//     // begin read vlocal
//     delete[] r;
//     delete[] vloc;
//     this->r = new double[mesh];
//     this->vloc = new double[mesh];
//     skip_comment(ifs);
//     ifs >> max_g;
//     for(int i = 0;i < mesh;++i)
//     {
//         r[i] = dg * i;
//         ifs >> this->vloc[i];
//     }
//     return 0;
// }

// void Pseudopot_upf::skip_comment(ifstream &ifs)
// {
//     char comment[200];
//     while(ifs.getline(comment,sizeof(comment)))
//     {
//         if(comment != "END COMMENT") break;
//     }
//     int temp1;
//     int temp2;
//     ifs >> temp1 >> temp2;
// }

int Pseudopot_upf::read_pseudo_blps(ifstream &ifs)
{
    // double bohr2a = 0.529177249;
    this->nlcc = false;
    this->tvanp = false;


    this->nbeta = 0;
    delete[] kkbeta;
    delete[] lll;
    this->kkbeta = new int[1];
    this->lll = new int[1];
    this->beta.create(1, 1);
    this->dion.create(1, 1);

    this->nwfc = 0;
    delete[] nn;
    delete[] jchi;
    delete[] jjj;
    this->nn = new int[1];
    this->jchi = new double[1];
    this->jjj = new double[1];
    ZEROS(nn, 1);
    ZEROS(jchi, 1);
    ZEROS(jjj, 1);

    ifs >> this->psd;
    // if(!SCAN_BEGIN(ifs,"BLPS")) WARNING_QUIT("read_pp_blps","Find no PP_HEADER");
    ifs.ignore(300, '\n');

    double zatom;
    double zion;
    ifs >> zatom >> zion;
    this->zp = static_cast<int>(zion);
    ifs.ignore(300, '\n');

    int pspcod, pspxc, lloc, r2well;
    ifs >> pspcod >> pspxc >> this->lmax >> lloc >> this->mesh >> r2well;

    if(pspxc == 2)
    {
        this->dft[0] = "LDA";
        this->dft[1] = "LDA";
    }
    else if (pspxc == 11)
    {
        this->dft[0] = "GGA";
        this->dft[1] = "GGA";
    }
    this->dft[2] = "NOGX";
    this->dft[3] = "NOGC";

    ifs.ignore(300, '\n');
    ifs.ignore(300, '\n');
    ifs.ignore(300, '\n');
    ifs.ignore(300, '\n');
    ifs.ignore(300, '\n');

    assert(mesh > 0);

    delete[] r;
    delete[] rab;
    delete[] vloc;
    this->r = new double[mesh]; // Bohr
    this->rab = new double[mesh];
    this->vloc = new double[mesh]; // Hatree
    ZEROS(r,mesh);
    ZEROS(rab,mesh);
    ZEROS(vloc,mesh);
    int num;
    for(int i = 0;i < mesh; ++i)
    {
        ifs >> num >> this->r[i] >> this->vloc[i];
        this->vloc[i] = this->vloc[i]*2;
    }
    rab[0] = r[1] - r[0];
    for(int i = 1; i < mesh - 1; ++i)
    {
        rab[i] = (r[i+1] - r[i-1])/2.0;
    }
    rab[mesh-1] = r[mesh-1] - r[mesh-2];

    delete[] rho_at;
    this->rho_at = new double[mesh];
    ZEROS(rho_at,mesh);
    // double charge = zion/(4.0/3.0*3.1415926535*r[mesh-1]*r[mesh-1]*r[mesh-1]);
    // double charge = 1;
    double charge = zion/r[mesh-1];
    for(int i = 0;i < mesh; ++i)
    {
        rho_at[i] = charge;
    }
    // cout<<"mesh="<<this->mesh<<endl;
    return 0;
}

//parameters
//read_pp.h <--> blps_real
//nv             -
//psd            head
//pp_type(NC or US) -
//tvanp          False
//nlcc           False
//dft            pspxc 2->lda, 11->gga
//zp             zion
//etotps         -
//ecutwfc        -
//ecutrho        -
//lmax           lmax
//mesh           mmax
//nwfc           -
//nbeta          -
//els            -
//lchi           -
//oc             -

//rab            rab[ir]=(r[ir+1]-r[ir-1])/2.0
//rho_atc(nonlocal) -
//vloc
//chi            -
//rho_at         -

//lll            -
//kkbeta         -
//beta           -
//dion           -

//nn             -
//jchi           -
//jjj            -
//nd             -