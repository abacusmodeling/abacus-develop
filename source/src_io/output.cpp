/* output.cpp */

#include "output.h"

void output::printrm(ofstream &ofs,const string &s, const matrix &m, const double &limit) const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    ofs << "\n " << s << " nr=" << b1 << " nc=" << b2 ;
    if (b1*b2 == 0) return;
    for (int i = 0;i < b1;i++)
    {
//        ofs<<"\n row = "<<i;
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) ofs << "\n ";
			
			if ( abs(m(i,j)) > limit) 
			{
				ofs << setprecision(6) << setw(12) << m(i,j);
            }
			else 
			{
				ofs << setw(12) << "0";
			}
        }
    }
	ofs << endl;
    return;
}

void output::printrm(const string &s, const matrix &m, const double &limit) const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    cout << "\n " << s << " nr=" << b1 << " nc=" << b2 ;
    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
		//cout << "\n row=" << i;
        for (int j = 0;j < b2;j++)
        {
//            if (j % 4 == 0) cout << "\n ";
			if (j % 8 == 0) cout << "\n ";
//            if ( abs(m(i,j)) > limit) cout << setprecision(15) << setw(20) << m(i,j);
            if ( abs(m(i,j)) > limit) cout << setprecision(6) << setw(12) << m(i,j);
//            else cout<<setw(20)<<"0";
            else cout<<setw(12)<<"0";
        }
//		cout << "\n";
    }

	cout << endl;
    return;
}

void output::printcm_norm(ofstream &ofs, const string &s, const ComplexMatrix &m, const double &limit)const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    ofs << "\n" << s << " nr=" <<  b1 << " nc=" << b2 ;

    if (b1*b2 == 0) return;
    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) ofs << "\n ";
            double norm = ( conj( m(i, j) ) * m(i, j) ).real();
            if ( abs(norm) > limit) ofs<< setw(12) << sqrt(norm);
            else ofs<<setw(12)<<"0";
        }
    }
    return;
}

void output::printcm_norm(const string &s, const ComplexMatrix &m, const double &limit)const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    cout << "\n" << s << " nr=" <<  b1 << " nc=" << b2 ;

    if (b1*b2 == 0) return;
    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) cout << "\n ";
            double norm = ( conj( m(i, j) ) * m(i, j) ).real();
            if ( abs(norm) > limit) cout<< setw(12) << sqrt(norm);
            else cout<<setw(12)<<"0";
        }
    }
    return;
}


void output::printcm(ofstream &ofs, const string &s, const ComplexMatrix &m) const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    ofs << "\n" << s << " " <<  b1 << " " << b2 ;

    if (b1*b2 == 0) return;

    ofs.setf(ios::scientific, ios::floatfield);
    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j%2 == 0) ofs << "\n ";
            ofs << setw(20) << m(i, j).real() << setw(20) << m(i, j).imag();
        }
    }
    ofs.unsetf(ios::scientific);
    return;
}//end print cm

void output::printcm(const string &s, const ComplexMatrix &m) const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    cout << "\n" << s << " nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 4 == 0) cout << "\n ";
            cout << setw(18) << m(i, j);
        }
    }
    return;
}

void output::printcm_real_limit_hermit(const string &s, const ComplexMatrix &m,const double &limit) const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    cout << "\n" << s << "  nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) cout << "\n ";
            const double x = m(i,j).real();
            if ( abs(x) < limit || i < j) cout<<setw(12)<<"0";
            else cout << setw(12) << x;
        }
    }
    return;
}

void output::printcm_real(const string &s, const ComplexMatrix &m,const double &limit) const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    cout << "\n " << s << "  nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) cout << "\n ";
            const double x = m(i,j).real();
            if ( abs(x) < limit ) cout<<setw(12)<<"0";
            else cout << setw(12) << x;
        }
    }
    return;
}



void output::printcm_imag(const string &s, const ComplexMatrix &m,const double &limit) const
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    cout << "\n " << s << "  nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) cout << "\n ";
            const double x = m(i,j).imag();
            if ( abs(x) < limit ) cout<<setw(12)<<"0";
            else cout << setw(12) << x;
        }
    }
    return;
}





void output::printr3_d(ofstream &ofs, const string &s,const realArray &u) const
{
    const int b1 = u.getBound1();
    const int b2 = u.getBound2();
    const int b3 = u.getBound3();
    ofs << "\n\n " << s << "  b1 = " << b1
    << " b2 = " << b2 << " b3 = " << b3;

    if (u.getSize() == 0) return;

    for (int i = 0;i < b1;i++)
        for (int j = 0;j < b2;j++)
            for (int k = 0;k < b3;k++)
            {
                if (k % 4 == 0)
                {
                    ofs << "\n";
                }
                ofs << setw(18) << u(i, j, k);
            }
    return;
}//end printr3_d

void output::printr4_d(ofstream &ofs, const string &s,const realArray &u) const
{
    const int b1 = u.getBound1();
    const int b2 = u.getBound2();
    const int b3 = u.getBound3();
    const int b4 = u.getBound4();
    ofs << "\n\n " << s << "  b1 = " << b1 << " b2 = " << b2
    << " b3 = " << b3 << " b4 = " << b4;

    if (u.getSize() == 0) return;

    for (int i = 0;i < b1;i++)
        for (int j = 0;j < b2;j++)
            for (int k = 0;k < b3;k++)
                for (int m = 0;m < b4;m++)
                {
                    if (m % 4 == 0)
                    {
                        ofs << "\n";
                    }
                    ofs << setw(15) << u(i, j, k, m);
                }
}//end print4_d

void output::printM3(ofstream &ofs,const string &description, const Matrix3 &m)const
{
    ofs << " " << description << endl;
	ofs << setiosflags(ios::showpos);
    ofs << " " << setw(20) << m.e11 << setw(20) << m.e12 << setw(20) << m.e13 
	<< "\n " << setw(20) << m.e21 << setw(20) << m.e22 << setw(20) << m.e23 
	<< "\n " << setw(20) << m.e31 << setw(20) << m.e32 << setw(20) << m.e33 << endl;
	ofs << resetiosflags(ios::showpos);
    return;
}

void output::printM3(const string &description, const Matrix3 &m)const
{
    cout << "\n " << description << endl;
    cout << setw(20) << m.e11 << setw(20) << m.e12 << setw(20) << m.e13 << "\n"
         << setw(20) << m.e21 << setw(20) << m.e22 << setw(20) << m.e23 << "\n"
         << setw(20) << m.e31 << setw(20) << m.e32 << setw(20) << m.e33 << endl;
    return;
}




