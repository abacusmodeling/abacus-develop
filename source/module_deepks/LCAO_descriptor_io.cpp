//wenfei 2021-11-17
#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "npy.hpp"
//============================
//DeePKS Part 3
//subroutines that deals with io as well as interface with libtorch, libnpy
//============================

//calculates descriptors from projected density matrices
void LCAO_Deepks::cal_descriptor(void)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_descriptor");

    //init pdm_tensor and d_tensor
    torch::Tensor tmp;

    //if pdm_tensor and d_tensor is not empty, clear it !!
    if (!this->d_tensor.empty())
    {
        this->d_tensor.erase(this->d_tensor.begin(), this->d_tensor.end());
    }
    if (!this->pdm_tensor.empty())
    {
        this->pdm_tensor.erase(this->pdm_tensor.begin(), this->pdm_tensor.end());
    }

    for (int inl = 0;inl < this->inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        tmp = torch::ones({ nm, nm }, torch::TensorOptions().dtype(torch::kFloat64));
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                tmp.index_put_({m1, m2}, this->pdm[inl][m1 * nm + m2]);
            }
        }
        //torch::Tensor tmp = torch::from_blob(this->pdm[inl], { nm, nm }, torch::requires_grad());
        tmp.requires_grad_(true);
        this->pdm_tensor.push_back(tmp);
        this->d_tensor.push_back(torch::ones({ nm }, torch::requires_grad(true)));
    }

    //cal d_tensor
    for (int inl = 0;inl < inlmax;++inl)
    {
        torch::Tensor vd;
        std::tuple<torch::Tensor, torch::Tensor> d_v(this->d_tensor[inl], vd);
        //d_v = torch::symeig(pdm_tensor[inl], /*eigenvalues=*/true, /*upper=*/true);
        d_v = torch::linalg::eigh(pdm_tensor[inl], /*uplo*/"U");
        d_tensor[inl] = std::get<0>(d_v);
    }
    return;
}

//calculates gradient of descriptors from gradient of projected density matrices
void LCAO_Deepks::cal_gvx(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks","cal_gvx");
    //preconditions
    this->cal_gvdm(nat);
    if(!gdmr_vector.empty())
    {
        gdmr_vector.erase(gdmr_vector.begin(),gdmr_vector.end());
    }

    //gdmr_vector : nat(derivative) * 3 * inl(projector) * nm * nm
    if(GlobalV::MY_RANK==0)
    {
        //make gdmx as tensor
        int nlmax = this->inlmax/nat;
        for (int nl=0;nl<nlmax;++nl)
        {
            std::vector<torch::Tensor> bmmv;
            for (int ibt=0;ibt<nat;++ibt)
            {
                std::vector<torch::Tensor> xmmv;
                for (int i=0;i<3;++i)
                {
                    std::vector<torch::Tensor> ammv;
                    for (int iat=0; iat<nat; ++iat)
                    {
                        int inl = iat*nlmax + nl;
                        int nm = 2*this->inl_l[inl]+1;
                        std::vector<double> mmv;
                        for (int m1=0;m1<nm;++m1)
                        {
                            for(int m2=0;m2<nm;++m2)
                            {
                                if(i==0) mmv.push_back(this->gdmx[ibt][inl][m1*nm+m2]);
                                if(i==1) mmv.push_back(this->gdmy[ibt][inl][m1*nm+m2]);
                                if(i==2) mmv.push_back(this->gdmz[ibt][inl][m1*nm+m2]);
                            }
                        }//nm^2
                        torch::Tensor mm = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64) ).reshape({nm, nm});    //nm*nm
                        ammv.push_back(mm);
                    }
                    torch::Tensor amm = torch::stack(ammv, 0);  //nat*nm*nm
                    xmmv.push_back(amm);
                }
                torch::Tensor bmm = torch::stack(xmmv, 0);  //3*nat*nm*nm
                bmmv.push_back(bmm); 
            }
            this->gdmr_vector.push_back(torch::stack(bmmv, 0)); //nbt*3*nat*nm*nm
        }
        assert(this->gdmr_vector.size()==nlmax);

        //einsum for each inl: 
        //gdmr_vector : b:nat(derivative) * x:3 * a:inl(projector) * m:nm * n:nm
        //gevdm_vector : a:inl * v:nm (descriptor) * m:nm (pdm, dim1) * n:nm (pdm, dim2)
        //gvx_vector : b:nat(derivative) * x:3 * a:inl(projector) * m:nm(descriptor)
        std::vector<torch::Tensor> gvx_vector;
        for (int nl = 0;nl<nlmax;++nl)
        {
            gvx_vector.push_back(at::einsum("bxamn, avmn->bxav", {this->gdmr_vector[nl], this->gevdm_vector[nl]}));
        }
        
        // cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
        // concatenate index a(inl) and m(nm)
        this->gvx_tensor = torch::cat(gvx_vector, -1);

        assert(this->gvx_tensor.size(0) == nat);
        assert(this->gvx_tensor.size(1) == 3);
        assert(this->gvx_tensor.size(2) == nat);
        assert(this->gvx_tensor.size(3) == this->des_per_atom);
    }

    return;
}

void LCAO_Deepks::cal_gvx_k(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks","cal_gvx");
    //preconditions
    this->cal_gvdm(nat);

    if(!gdmr_vector.empty())
    {
        gdmr_vector.erase(gdmr_vector.begin(),gdmr_vector.end());
    }
    
    if(GlobalV::MY_RANK==0)
    {
        //make gdmx as tensor
        int nlmax = this->inlmax/nat;
        for (int nl=0;nl<nlmax;++nl)
        {
            std::vector<torch::Tensor> bmmv;
            for (int ibt=0;ibt<nat;++ibt)
            {
                std::vector<torch::Tensor> xmmv;
                for (int i=0;i<3;++i)
                {
                    std::vector<torch::Tensor> ammv;
                    for (int iat=0; iat<nat; ++iat)
                    {
                        int inl = iat*nlmax + nl;
                        int nm = 2*this->inl_l[inl]+1;
                        std::vector<double> mmv;
                        for (int m1=0;m1<nm;++m1)
                        {
                            for(int m2=0;m2<nm;++m2)
                            {
                                if(i==0) mmv.push_back(this->gdmx[ibt][inl][m1*nm+m2]);
                                if(i==1) mmv.push_back(this->gdmy[ibt][inl][m1*nm+m2]);
                                if(i==2) mmv.push_back(this->gdmz[ibt][inl][m1*nm+m2]);
                            }
                        }//nm^2
                        torch::Tensor mm = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64) ).reshape({nm, nm});    //nm*nm
                        ammv.push_back(mm);
                    }
                    torch::Tensor amm = torch::stack(ammv, 0);  //nat*nm*nm
                    xmmv.push_back(amm);
                }
                torch::Tensor bmm = torch::stack(xmmv, 0);  //3*nat*nm*nm
                bmmv.push_back(bmm); 
            }
            this->gdmr_vector.push_back(torch::stack(bmmv, 0)); //nbt*3*nat*nm*nm
        }
        assert(this->gdmr_vector.size()==nlmax);

        //einsum for each inl: 
        std::vector<torch::Tensor> gvx_vector;
        for (int nl = 0;nl<nlmax;++nl)
        {
            gvx_vector.push_back(at::einsum("bxamn, avmn->bxav", {this->gdmr_vector[nl], this->gevdm_vector[nl]}));
        }//
        
        // cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
        this->gvx_tensor = torch::cat(gvx_vector, -1);

        assert(this->gvx_tensor.size(0) == nat);
        assert(this->gvx_tensor.size(1) == 3);
        assert(this->gvx_tensor.size(2) == nat);
        assert(this->gvx_tensor.size(3) == this->des_per_atom);
    }

    return;
}

void LCAO_Deepks::load_model(const string& model_file)
{
    ModuleBase::TITLE("LCAO_Deepks", "load_model");

    try
	{
        this->module = torch::jit::load(model_file);
    }
    catch (const c10::Error& e)

	{
        std::cerr << "error loading the model" << std::endl;
        return;
    }
	return;
}

//obtain from the machine learning model dE_delta/dDescriptor
void LCAO_Deepks::cal_gedm(const int nat)
{
    //using this->pdm_tensor
    ModuleBase::TITLE("LCAO_Deepks", "cal_gedm");

    //forward
    std::vector<torch::jit::IValue> inputs;
    //input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(this->d_tensor, /*dim=*/0).reshape({ nat, this->des_per_atom }));
    std::vector<torch::Tensor> ec;
    ec.push_back(module.forward(inputs).toTensor());    //Hartree
    this->E_delta = ec[0].item().toDouble() * 2;//Ry; *2 is for Hartree to Ry
    
    //cal gedm
    std::vector<torch::Tensor> gedm_shell;
    gedm_shell.push_back(torch::ones_like(ec[0]));
    this->gedm_tensor = torch::autograd::grad(ec, this->pdm_tensor, gedm_shell, /*retain_grad=*/true, /*create_graph=*/false, /*allow_unused=*/true);

    //gedm_tensor(Hartree) to gedm(Ry)
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                int index = m1 * nm + m2;
                //*2 is for Hartree to Ry
                this->gedm[inl][index] = this->gedm_tensor[inl].index({ m1,m2 }).item().toDouble() * 2;
            }
        }
    }
    return;
}

void LCAO_Deepks::cal_gedm_k(const int nat)
{
    //using this->pdm_tensor
    ModuleBase::TITLE("LCAO_Deepks", "cal_gedm");

    //forward
    std::vector<torch::jit::IValue> inputs;
    //input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(this->d_tensor, /*dim=*/0).reshape({ nat, this->des_per_atom }));
    std::vector<torch::Tensor> ec;
    ec.push_back(module.forward(inputs).toTensor());    //Hartree
    this->E_delta = ec[0].item().toDouble() * 2;//Ry; *2 is for Hartree to Ry
    
    //cal gedm
    std::vector<torch::Tensor> gedm_shell;
    gedm_shell.push_back(torch::ones_like(ec[0]));
    this->gedm_tensor = torch::autograd::grad(ec, this->pdm_tensor, gedm_shell, /*retain_grad=*/true, /*create_graph=*/false, /*allow_unused=*/true);

    //gedm_tensor(Hartree) to gedm(Ry)
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                int index = m1 * nm + m2;
                //*2 is for Hartree to Ry
                this->gedm[inl][index] = this->gedm_tensor[inl].index({ m1,m2 }).item().toDouble() * 2;
            }
        }
    }
    return;
}

//dDescriptor / dprojected density matrix
void LCAO_Deepks::cal_gvdm(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gvdm");
    if(!gevdm_vector.empty())
    {
        gevdm_vector.erase(gevdm_vector.begin(),gevdm_vector.end());
    }
    //cal gevdm(d(EigenValue(D))/dD)
    int nlmax = inlmax/nat;
    for (int nl=0;nl<nlmax;++nl)
    {
        std::vector<torch::Tensor> avmmv;
        for (int iat = 0;iat<nat;++iat)
        {
            int inl = iat*nlmax+nl;
            int nm = 2*this->inl_l[inl]+1;
            //repeat each block for nm times in an additional dimension
            torch::Tensor tmp_x = this->pdm_tensor[inl].reshape({nm, nm}).unsqueeze(0).repeat({nm, 1, 1});
            //torch::Tensor tmp_y = std::get<0>(torch::symeig(tmp_x, true));
            torch::Tensor tmp_y = std::get<0>(torch::linalg::eigh(tmp_x, "U"));
            torch::Tensor tmp_yshell = torch::eye(nm, torch::TensorOptions().dtype(torch::kFloat64));
            std::vector<torch::Tensor> tmp_rpt;     //repeated-pdm-tensor (x)
            std::vector<torch::Tensor> tmp_rdt; //repeated-d-tensor (y)
            std::vector<torch::Tensor> tmp_gst; //gvx-shell
            tmp_rpt.push_back(tmp_x);
            tmp_rdt.push_back(tmp_y);
            tmp_gst.push_back(tmp_yshell);
            std::vector<torch::Tensor> tmp_res;
            tmp_res = torch::autograd::grad(tmp_rdt, tmp_rpt, tmp_gst, false, false, /*allow_unused*/true); //nm(v)**nm*nm
            avmmv.push_back(tmp_res[0]);
        }
        torch::Tensor avmm = torch::stack(avmmv, 0); //nat*nv**nm*nm
        this->gevdm_vector.push_back(avmm);
    }
    assert(this->gevdm_vector.size() == nlmax);
    return;
}

void LCAO_Deepks::print_H_V_delta(void)
{
    ModuleBase::TITLE("LCAO_Deepks", "print_H_V_delta");

    ofstream ofs;
    stringstream ss;

    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "H_V_delta.dat";

    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }

    ofs << "E_delta(Ry) from deepks model: " << this->E_delta << std::endl;
    ofs << "E_delta(eV) from deepks model: " << this->E_delta * ModuleBase::Ry_to_eV << std::endl;
    ofs << "H_delta(Ry)(gamma only)) from deepks model: " << std::endl;

    for (int i = 0;i < GlobalV::NLOCAL;++i)
    {
        for (int j = 0;j < GlobalV::NLOCAL;++j)
        {
            ofs<< std::setw(12)<< this->H_V_delta[i * GlobalV::NLOCAL + j] << " ";
        }
        ofs << std::endl;
    }
    ofs << "H_delta(eV)(gamma only)) from deepks model: " << std::endl;

    for (int i = 0;i < GlobalV::NLOCAL;++i)
    {
        for (int j = 0;j < GlobalV::NLOCAL;++j)
        {
            ofs<< std::setw(12)<< this->H_V_delta[i * GlobalV::NLOCAL + j] *ModuleBase::Ry_to_eV<< " ";
        }
        ofs << std::endl;
    }

    GlobalV::ofs_running << " H_delta has been printed to " << ss.str() << std::endl;
    return;
}


void LCAO_Deepks::print_F_delta(const string& fname, const UnitCell_pseudo &ucell)
{
    ModuleBase::TITLE("LCAO_Deepks", "print_F_delta");

    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"<< fname ;

    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }

    ofs << "F_delta(Hatree/Bohr) from deepks model: " << std::endl;
    ofs << std::setw(12) << "type" << std::setw(12) << "atom" << std::setw(15) << "dF_x" << std::setw(15) << "dF_y" << std::setw(15) << "dF_z" << std::endl;

    for (int it = 0;it < ucell.ntype;++it)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;++ia)
        {
            int iat = ucell.itia2iat(it, ia);
            ofs << std::setw(12) << ucell.atoms[it].label << std::setw(12) << ia
                << std::setw(15) << this->F_delta(iat, 0) / 2 << std::setw(15) << this->F_delta(iat, 1) / 2
                << std::setw(15) << this->F_delta(iat, 2) / 2 << std::endl;
        }
    }

    ofs << "F_delta(eV/Angstrom) from deepks model: " << std::endl;
    ofs << std::setw(12) << "type" << std::setw(12) << "atom" << std::setw(15) << "dF_x" << std::setw(15) << "dF_y" << std::setw(15) << "dF_z" << std::endl;

    for (int it = 0;it < ucell.ntype;++it)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;++ia)
        {
            int iat = ucell.itia2iat(it, ia);
            ofs << std::setw(12) << ucell.atoms[it].label << std::setw(12)
                << ia << std::setw(15) << this->F_delta(iat, 0) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A
                << std::setw(15) << this->F_delta(iat, 1) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A
                << std::setw(15) << this->F_delta(iat, 2) * ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A << std::endl;
        }
    }

    GlobalV::ofs_running << " F_delta has been printed to " << ss.str() << std::endl;
    ofs.close();

    return;
}


void LCAO_Deepks::save_npy_d(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_d");
    //save descriptor in .npy format
    vector<double> npy_des;
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2*inl_l[inl] + 1;
        for(int im=0;im<nm;im++)
        {
            npy_des.push_back(this->d_tensor[inl].index({im}).item().toDouble());
        }
    }
    const long unsigned dshape[] = {(long unsigned) nat, (long unsigned) this->des_per_atom };
    if (GlobalV::MY_RANK == 0)
    {
        npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, dshape, npy_des);
    }
    return;
}


void LCAO_Deepks::save_npy_e(const double &e, const std::string &e_file)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_e");
    //save e_base
    const long unsigned eshape[] = { 1 };
    vector<double> npy_e;
    npy_e.push_back(e);
    npy::SaveArrayAsNumpy(e_file, false, 1, eshape, npy_e);
    return;
}

void LCAO_Deepks::save_npy_f(const ModuleBase::matrix &f, const std::string &f_file, const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_f");
    //save f_base
    //caution: unit: Rydberg/Bohr
    const long unsigned fshape[] = {(long unsigned) nat, 3 };
    vector<double> npy_f;
    for (int iat = 0;iat < nat;++iat)
    {
        for (int i = 0;i < 3;i++)
        {
            npy_f.push_back(f(iat, i));
        }
    }
    npy::SaveArrayAsNumpy(f_file, false, 2, fshape, npy_f);
    return;
}

void LCAO_Deepks::save_npy_gvx(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_npy_gvx");
    //save grad_vx.npy (when  force label is in use)
    //unit: /Bohr
    const long unsigned gshape[] = {(long unsigned) nat, 3, nat, this->des_per_atom};
    vector<double> npy_gvx;
    for (int ibt = 0;ibt < nat;++ibt)
    {
        for (int i = 0;i < 3;i++)
        {
            for (int iat = 0;iat < nat;++iat)
            {
                for(int p=0;p<this->des_per_atom;++p)
                {
                    npy_gvx.push_back(this->gvx_tensor.index({ ibt, i, iat, p }).item().toDouble());
                }
            }
        }
    }
    npy::SaveArrayAsNumpy("grad_vx.npy", false, 4, gshape, npy_gvx);
    return;
}

#endif