#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>


int main()
{
    double
        kp = 5.0e4, // M-1 s-1
        km = 2.0e-8, // s-1
        kn = 5.0e-5, // M-1 s-1
        k2 = 0.0,
        nc = 2.0,
        n2 = 2.0;
    double dt = 1.0;

    unsigned int maxlen = 600000;
    unsigned int maxtime = static_cast<int>(80.0*3600.0/dt);
    std::vector<double> f(maxlen, 0.0);
    std::vector<double> P(maxtime, 0.0);
    std::vector<double> M(maxtime, 0.0);
    std::vector<double> Q(maxtime, 0.0);
    M.at(0) = 1.0e-8;
    P.at(0) = M.at(0) / 5000.0;
	f.at(5000) = P.at(0);
    double m0 = 1.0e-6 - M.at(0);
    double m = m0;
    double ko = 0.0*kp*m0;
    for (unsigned int timestep = 1; timestep < maxtime; ++timestep)
    {
        double sec_nucl_surf = 0.0;
        // for (unsigned int i = (int) nc; i < maxlen; ++i)
        // {
        //     sec_nucl_surf += static_cast<double>(i)*f.at(i);
        // }

        double totalAbove = 0.0;
        auto f_old = f;
        for (unsigned int j = (maxlen-2); j >= (int) nc; j--)
        {
            double dj = static_cast<double>(j);
            double fj_less_1 = f_old.at(j-1);
            double fj_more_1 = f_old.at(j+1);
            double fj = f_old.at(j);
            totalAbove += fj_more_1;
            double dfj_dt =
                (2.0*m*kp*fj_less_1) +
                (-2.0*m*kp*fj) +
                (2*ko*fj_more_1) +
                (-2*ko*fj) +
                (-km*(dj - 1.0)*fj) +
                (2.0*km*totalAbove);

            if (j == (int) n2) { dfj_dt += (k2*std::pow(m, n2)*sec_nucl_surf); }
            if (j == (int) nc) { dfj_dt += (kn*std::pow(m, nc)); }

            f.at(j) = f_old.at(j) + dfj_dt*dt;
        }

        for (unsigned int i = (int) nc; i < maxlen; ++i)
        {
            P.at(timestep) += f.at(i);
            M.at(timestep) += static_cast<double>(i)*f.at(i);
            Q.at(timestep) += std::pow(static_cast<double>(i),2.0)*f.at(i);
        }

        m = m0 - M.at(timestep);

        if (timestep%100 == 0){ std::cout << static_cast<double>(timestep)*dt/3600.0 << std::endl; }
    }

    FILE* outfile = fopen("tmoments_sigmoidal_fig5_0pct.dat","w");
    for (unsigned int timestep = 0; timestep < maxtime; ++timestep)
    {
        // std::cout << P.at(timestep) << " " <<  M.at(timestep) << std::endl;
        fprintf(outfile, "%e, %e, %e, %e\n", static_cast<double>(timestep)*dt/3600.0, P.at(timestep), M.at(timestep), Q.at(timestep));
    }

    fclose(outfile);

    return 0;
}
