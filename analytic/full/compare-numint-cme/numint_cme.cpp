#include <algorithm>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>


int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		printf("Error: insufficient arguments.\nSyntax: numint_cme.x max_length time filename_output\n");
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}

	char *endptr;
	long int ml = strtol(argv[1], &endptr, 10);
	if (!*argv[1] || *endptr)
	{
		std::cerr << "Error: Invalid length provided: " << argv[1] << std::endl;
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}
	
	char *endptr2;
	double t = strtod(argv[2], &endptr2);
	if (!*argv[2] || *endptr2)
	{
		std::cerr << "Error: Invalid time provided: " << argv[2] << std::endl;
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}
	
	char *endptr3;
	double dt = strtod(argv[3], &endptr3);
	if (!*argv[3] || *endptr3)
	{
		std::cerr << "Error: Invalid time provided: " << argv[2] << std::endl;
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}
	
	FILE* outfile = fopen(argv[4],"w");
	if (outfile == nullptr)
	{
		printf("Error: unable to open output file for writing.\n");
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}

    double 	a 	= 0.0,				// alpha
			kp 	= 1.+a,				// elongation, k_plus
			kd 	= 1.,				// defragmentation
			km 	= 1.,				// fragmentation, k_minus
			kn 	= 1.,				// primary nucleation
			nc 	= 2.,				// critical nucleus size
			n2 	= 2.,				// critical nucleus size for 2ary nucleation
			nm = std::min(nc, n2); 	// smallest segment length to keep track of

    unsigned int 	maxlen 		= ml,								// corresponds to total number of monomer units L
    				maxtime 	= static_cast<int>(t/dt);			// time for which to integrate

	std::vector<double> f(maxlen, 0.0);		// current (continuously updated) concentrations of segments of each length

	std::vector<double> mon(maxtime, 0.0);	// concentration of monomers over time
    std::vector<double> P(maxtime, 0.0);	// 0th moment of f(j) over time
    std::vector<double> M(maxtime, 0.0);	// 1st moment
    std::vector<double> Q(maxtime, 0.0);	// 2nd moment

	M.at(0) 	= 1.;
    P.at(0) 	= M.at(0);
    
	double 	m0 	= M.at(0),			// initial concentration of unbound monomers
			m 	= m0;				// current concentration of unbound monomers
    

	for (unsigned int elem = 0; elem < f.size(); ++elem)
	{
		if (elem == 0)	std::cout << m;
		else			std::cout << f.at(elem);
		std::cout << " ";
	} std::cout << "|| total bound conc = " << std::accumulate(f.begin(), f.end(), 0.0) << std::endl;

	for (unsigned int timestep = 1; timestep < maxtime; ++timestep)
    {
        double 	total_above 	= 0.0;		// total concentration of segments with length greater than j
        auto 	f_old 			= f;		// copy concentrations for use in calculating rates, since f will be overwritten at each timestep
		
		// iterate over segment lengths
        for (
				unsigned int j = maxlen;
				j >= (unsigned int) nm;
				j--
			)
        {
            double 	j_d 			= static_cast<double>(j),
					fj_less_1		= f_old.at(j-2),			// population with one monomer less than current length
					fj				= f_old.at(j-1);				// population with current length
			double& fj_new			= f.at(j-1);
 
			double 	fj_more_1 		= 0.;
			if (j < maxlen)
			{
				fj_more_1			= f_old.at(j);			// population with one monomer more than current length
				total_above 		+= fj_more_1;				// update running total
			}
			
			// increase in population due to elongation of segments shorter by 1 monomer
			// don't need special case for j = nc because f defined to be zero for j < nc
			fj_new += dt*(2.*m*kp*fj_less_1);

			// decrease in population due to elongation of segments of current length (by one monomer)
			fj_new -= dt*(2.*m*kp*fj);

			/* DISSOCIATION ------- REMOVED FROM THE SCHEME OF COHEN, VENDRUSCOLO, et al. 2011 
			 *
			 * TECHNICALLY THEIR METHOD OVERCOUNTS REMOVAL OF A SINGLE MONOMER, SINCE IT IS INCLUDED IN BOTH
			 * THE FRAGMENTATION AND DISSOCIATION TERMS
			 *
			 * I BELIEVE THEY DO THIS BECAUSE THEY WANT TO BE ABLE TO TREAT DISSOCIATION SEPARATELY (I.E.
			 * WITH A DIFFERENT, MOST LIKELY LARGER, RATE CONSTANT) WHICH MAKES THE CONTRIBUTION FROM THE
			 * DISSOCIATION TERMS NEGLIGIBLE
				
				// increase in population due to dissociation of segments of length up to one monomer greater than current length
				// only if the current length is less than the maximum segment length
				double ko = 1.0;
				if (j < (maxlen-1)) fj_new += dt*(2.*ko*fj_more_1);

				// decrease in population due to dissociation of segments of current length
				// only if the current length is greater than the minimum segment length (because dissociation of minimum-length segments has already been treated in the fragmentation term)
				if (j > nc) fj_new -= dt*(2.*ko*fj);
			*/
			
			// decrease in population due to fragmentation of segments of current length
			fj_new -= dt*(km*j_d*fj);

			// increase in population to fragmentation of segments of greater length
			fj_new += dt*(2.*km*total_above);

			// "de-fragmentation" term in CME
			// iterate over length of short component
			// and find the length of the corresponding long component needed to form a new segment of length j
			for (	
					unsigned int l_short = nc;
					l_short <= (unsigned int) std::round(j/2);
					++l_short
				)
			{
				unsigned int 	l_long 		= j - l_short;
				double 			pop_change 	= dt*(kd*f_old.at(l_short)*f_old.at(l_long));

				fj_new 						+= pop_change;
				f.at(l_short-1)				-= pop_change;
				f.at(l_long-1)				-= pop_change;
			}

            if (j == (unsigned int) nc)
			{
				fj_new += dt*(kn*std::pow(m, nc));
			}
        }

		for (unsigned int elem = 0; elem < f.size(); ++elem)
		{
			if (elem == 0)	std::cout << m;
			else			std::cout << f.at(elem);
			std::cout << " ";
		} std::cout << "|| total bound conc = " << std::accumulate(f.begin(), f.end(), 0.0) << std::endl;

        for (unsigned int i = (int) nm; i < maxlen; ++i)
        {
			double i_d 		= static_cast<double>(i);
            P.at(timestep) 	+= f.at(i);
            M.at(timestep) 	+= f.at(i)*i_d;
            Q.at(timestep) 	+= f.at(i)*i_d*i_d;
        }

        m 					= m0 - M.at(timestep);
		mon.at(timestep) 	= m;

        //if (timestep%100 == 0){ std::cout << static_cast<double>(timestep)*dt << std::endl; }
    }

    for (unsigned int timestep = 0; timestep < maxtime; ++timestep)
    {
        fprintf(outfile,
				"%e, %e, %e, %e, %e, %e\n",
				static_cast<double>(timestep)*dt,
				mon.at(timestep),
				P.at(timestep),
				M.at(timestep),
				m0/(P.at(timestep) + mon.at(timestep)),
				Q.at(timestep));
    }

    fclose(outfile);


    return 0;
}
