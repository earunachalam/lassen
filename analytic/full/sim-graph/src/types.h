#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include <string>
#include <vector>

typedef std::vector<unsigned int> ustate;

class summableEdge
{
	protected:
		long double base_rateconst;
		long double rateconst;
		std::string color;

	public:
		summableEdge(long double base_rateconst = 0, std::string color = "#00f"): base_rateconst(base_rateconst), rateconst(base_rateconst), color(color) {};

	    summableEdge operator+=(const summableEdge& rhs)
	    {
	        this->rateconst = (this->rateconst + rhs.rateconst);
			this->base_rateconst = std::numeric_limits<long double>::quiet_NaN();
	        return *this;
	    }

	    void s_base_rateconst(const long double& base_rateconst)
	    {
	    	this->base_rateconst = base_rateconst;
	    }

	    long double g_base_rateconst()
	    {
	    	return this->base_rateconst;
	    }

	    void s_rateconst(const long double& rateconst)
	    {
			this->rateconst = rateconst;
	    }

	    long double g_rateconst()
	    {
	    	return this->rateconst;
	    }

	    void s_color(const std::string& color)
	    {
	    	this->color = color;
	    }

	    std::string g_color()
	    {
	    	return this->color;
	    }
};

#endif
