//This file contains interfaces with libtorch,
//including loading of model and calculating gradients
//as well as subroutines that prints the results for checking

//The file contains 8 subroutines:
//1. cal_descriptor : obtains descriptors which are eigenvalues of pdm
//      by calling torch::linalg::eigh
//2. check_descriptor : prints descriptor for checking
//3. cal_gvx : gvx is used for training with force label, which is gradient of descriptors, 
//      calculated by d(des)/dX = d(pdm)/dX * d(des)/d(pdm) = gdmx * gvdm
//      using einsum
//4. check_gvx : prints gvx into gvx.dat for checking
//5. cal_gvdm : d(des)/d(pdm)
//      calculated using torch::autograd::grad
//6. load_model : loads model for applying V_delta
//7. cal_gedm : calculates d(E_delta)/d(pdm)
//      this is the term V(D) that enters the expression H_V_delta = |alpha>V(D)<alpha|
//      caculated using torch::autograd::grad
//8. check_gedm : prints gedm for checking

#ifdef __DEEPKS

#include "LCAO_deepks.h"

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

void LCAO_Deepks::check_descriptor(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "check_descriptor");
    ofstream ofs("descriptor.dat");
    ofs<<std::setprecision(12);
    for(int ia=0;ia<nat;ia++)
    {
        for(int inl=0;inl<inlmax/nat;inl++)
        {
            int nm = 2*inl_l[inl]+1;
            for(int im=0;im<nm;im++)
            {
                const int ind=ia*inlmax/nat+inl;
                ofs << std::setprecision(10) << d_tensor[ind].index({im}).item().toDouble() << " ";
            }
            ofs << std::endl;
        }   
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

void LCAO_Deepks::check_gvx(const int nat)
{
    std::stringstream ss;
    ofstream ofs_x;
    ofstream ofs_y;
    ofstream ofs_z;

    ofs_x<<std::setprecision(12);
    ofs_y<<std::setprecision(12);
    ofs_z<<std::setprecision(12);

    for(int ia=0;ia<nat;ia++)
    {
        ss.str("");
        ss<<"gvx_"<<ia<<".dat";
        ofs_x.open(ss.str().c_str());
        ss.str("");
        ss<<"gvy_"<<ia<<".dat";
        ofs_y.open(ss.str().c_str());
        ss.str("");
        ss<<"gvz_"<<ia<<".dat";
        ofs_z.open(ss.str().c_str());
        for(int ib=0;ib<nat;ib++)
        {
            for(int inl=0;inl<inlmax/nat;inl++)
            {
                int nm = 2*inl_l[inl]+1;
                {
                    const int ind=ib*inlmax/nat+inl;
                    ofs_x << std::setprecision(10) << gvx_tensor.index({ia,0,ib,inl}).item().toDouble() << " ";
                    ofs_y << std::setprecision(10) << gvx_tensor.index({ia,1,ib,inl}).item().toDouble() << " ";
                    ofs_z << std::setprecision(10) << gvx_tensor.index({ia,2,ib,inl}).item().toDouble() << " ";
                }
            }
            ofs_x << std::endl;
            ofs_y << std::endl;
            ofs_z << std::endl;
        }
        ofs_x.close();
        ofs_y.close();
        ofs_z.close();        
    }
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
    inputs.push_back(torch::cat(this->d_tensor, 0).reshape({ nat, this->des_per_atom }));
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

void LCAO_Deepks::check_gedm()
{
    ofstream ofs("gedm.dat");
    for(int inl=0;inl<inlmax;inl++)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                int index = m1 * nm + m2;
                //*2 is for Hartree to Ry
                ofs << this->gedm[inl][index] << " ";
            }
        }   
        ofs << std::endl;     
    }
}

#endif