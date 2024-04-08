#include "ideas/LegacyToJSON.h"
#include <fstream>
#include <vector>
#include <stack>
#include <stdexcept>
#include <string>

namespace loki
{

namespace {

struct Lines
{
	Lines(std::istream& is);
	struct Line
	{
		Line(std::istream& is);
		bool empty() { return key.empty() && value.empty(); }
		std::string::size_type indent;
		std::string key;
		json_type value;
	private:
		static json_type parse_word(const std::string& str)
		{
			if (str.empty())
			{
				return json_type{};
			}
			std::string s(str);
			s = s.substr(0,str.find_last_not_of(" \n\r\t")+1);
			const std::string::size_type indent = s.find_first_not_of(" \n\r\t");
			s = s.substr(indent);
			try {
				return json_type::parse(s.begin(),s.end());
			}
			catch(std::exception& exc)
			{
				return s;
			}
		}
	};
	std::vector<Line> m_lines;
	using const_iterator = std::vector<Line>::const_iterator;
};


Lines::Line::Line(std::istream& is)
{
	std::string str;
	while (!std::getline(is,str).eof())
	{
		// remove "*whitespace '%' *characters" from the end of the line
		const std::string::size_type comment_pos = str.find('%');
		str = str.substr(0,comment_pos);
		str = str.substr(0,str.find_last_not_of(" \n\r\t")+1);
		if (!str.empty())
		{
			// skip line with only whitespace and/or a comment
			break;
		}
	}
	if (is.eof()) return;
	indent = str.find_first_not_of(" \n\r\t");
	str = str.substr(indent);
	std::stringstream ss(str);
	std::string word;
	ss >> word;
	if (word.back()==':')
	{
		key = word.substr(0,word.size()-1);
		if (key.empty())
		{
			std::runtime_error("Empty key found.");
		}
		if (ss.tellg()>=0)
		{
			std::string v(str.begin()+ss.tellg(),str.end());
			value = parse_word(v);
		}
	}
	else if (word=="-")
	{
		if (ss.tellg()>=0)
		{
			std::string v(str.begin()+ss.tellg(),str.end());
			value = parse_word(v);
		}
	}
	else
	{
		throw std::runtime_error("While reading '" + ss.str() + "': expected '<key>: [value]' or '- <value>'.");
	}
}

Lines::Lines(std::istream& is)
{
	for (Line line(is);!line.empty();line=Line(is))
	{
		m_lines.push_back(line);
	}
}

void readSection(json_type& parent, Lines::const_iterator& it, Lines::const_iterator end)
{
	Lines::const_iterator begin(it);
	while (it != end && it->indent>=begin->indent)
	{
		//std::cout << "INDENT = " << it->indent << ", KEY = " << it->key << ", VALUE = " << it->value << std::endl;
		if (it->value.empty())
		{
			json_type& sec = parent[it->key];
			auto old = it++;
			if (it==end || it->indent==old->indent)
			{
	//			sec = json_type(R"([])");
				continue;
			}
			readSection( sec, it, end);
		}
		else if (it->key.empty())
		{
			parent.push_back(it->value);
			++it;
		}
		else
		{
			parent[it->key] = it->value;
			++it;
		}
	}
}

} // namespace

json_type legacyToJSON(std::istream& is)
try 
{
	Lines lines(is);
	json_type root;
	Lines::const_iterator it = lines.m_lines.begin();
	readSection(root,it,lines.m_lines.end());
	return root;
}
catch (std::exception& exc)
{
	throw std::runtime_error(std::string("Error parsing stream: ") + exc.what());
}

json_type legacyToJSON(const std::filesystem::path& fname)
{
	std::ifstream ifs(fname);
	if (!ifs)
	{
		throw std::runtime_error("Error opening file '" + fname.generic_string() + "' for reading.");
	}
	return legacyToJSON(ifs);
}

} // namespace loki
