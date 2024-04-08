/** \file
 *
 *  Convert legacy LoKI-B input files to JSON.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Jan van Dijk
 *  \date   8 April 2024
 */

#include "ideas/LegacyToJSON.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace loki
{

namespace {

/** All lines in the input file that are non-empty after removing
 *  comments (starting with a %-sign) are stored in this class,
 *  which sets up a vector of Line objects. Valid lines in LoKI-B
 *  input lines consist of a number of whitespace characters,
 *  stored in Line member 'indent', followed by
 *
 *    1) <key>: [value]
 *    2) - value.
 *
 *  A key, followed by a value represents a regular key-value pair;
 *  the value can be an integer or floating point number, a literal
 *  true of false or a string. Note that the value can contain whitespace,
 *  everything after <key>: is considered to be part of the value,
 *  but leading and trailing whitespace (and comments) are trimmed.
 *  As an example, the line
 *  \verbatim
      pair:  A B % comment\endverbatim
 *  will result in value "A B".
 *  When [value] is omitted, a key-object pair will be inserted.
 *
 *  When the value is missing, the value is made a JSON object or array.
 *  The object/array contant is formed by the subsequent lines with a
 *  indentation level higher than the one that has the key. An array
 *  is created if these subsequent lines are of the form "- value",
 *  an object if these lines are of the form "<key>: [value]".
 *
 *  Finally, values are parsed mostly as in JSON: the literals true and
 *  false result in a boolean, literals like 42 as integers, 42.0 and 1e3
 *  as floating point numbers. (Only) if these conversions fail, a string
 *  value is produced, as for the value "A B".
 *
 *  The function readSection produces JSON output from a Lines object.
 *
 *  \sa readSection
 */
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
