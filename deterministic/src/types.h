#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include <string>
#include <vector>

typedef std::vector<unsigned int> ustate;

class activeEdge
{
	protected:
		long double m_base_rateconst;
		long double m_rateconst;
		std::string m_color;

	public:
		activeEdge(long double base_rateconst = 0, std::string color = "#00f"): m_base_rateconst(base_rateconst), m_rateconst(base_rateconst), m_color(color) {};

	    activeEdge operator+=(const activeEdge& rhs)
	    {
	        this->m_rateconst = (this->m_rateconst + rhs.m_rateconst);
			this->m_base_rateconst = std::numeric_limits<long double>::quiet_NaN();
	        return *this;
	    }

	    void s_base_rateconst(const long double& base_rateconst)
	    {
	    	m_base_rateconst = base_rateconst;
	    }

	    long double g_base_rateconst()
	    {
	    	return m_base_rateconst;
	    }

	    void s_rateconst(const long double& rateconst)
	    {
			m_rateconst = rateconst;
	    }

	    long double g_rateconst()
	    {
	    	return m_rateconst;
	    }

	    void s_color(const std::string& color)
	    {
	    	m_color = color;
	    }

	    std::string g_color()
	    {
	    	return m_color;
	    }
};

#endif
