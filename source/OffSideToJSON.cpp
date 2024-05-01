/** \file
 *
 *  Conversion of 'off-side rule' files to JSON.
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

#include "LoKI-B/OffSideToJSON.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <list>

namespace loki
{

namespace {

/** All lines in the input file that are non-empty after removing
 *  comments (starting with a %-sign) are stored in this class,
 *  which sets up a list of Line objects. Valid lines in LoKI-B
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
 *  The object/array content is formed by the subsequent lines with an
 *  indentation level higher than the one that has the key. An array
 *  is created if these subsequent lines are of the form "- value",
 *  an object if these lines are of the form "<key>: [value]".
 *
 *  Finally, values are parsed mostly as in JSON: the literals true and
 *  false result in a boolean, literals like 42 as integers, 42.0 and 1e3
 *  as floating point numbers. (Only) if these conversions fail, a string
 *  value is produced, as for the value "A B".
 *
 *  \sa readSection
 *
 *  \author Jan van Dijk
 *  \date   8 April 2024
 */
class Lines
{
public:
	Lines(std::istream& is);
	class Line
	{
	public:
		Line(std::istream& is, unsigned& line_number);
		bool empty() { return key().empty() && value().empty(); }
		void throw_line_exception(const std::string& msg) const;
		std::string::size_type indent() const { return m_indent; }
		unsigned line_number() const { return m_line_number; }
		const std::string& raw_string() const { return m_raw_string; }
		const std::string& key() const { return m_key; }
		const json_type& value() const { return m_value; }
	private:
		std::string::size_type m_indent;
		unsigned m_line_number;
		std::string m_raw_string;
		std::string m_key;
		json_type m_value;
		static json_type parse_word(const std::string& str);
	};
	using container_type = std::list<Line>;
	using const_iterator = container_type::const_iterator;
	const_iterator begin() const { return m_lines.begin(); }
	const_iterator end() const { return m_lines.end(); }
private:
	container_type m_lines;
};


Lines::Line::Line(std::istream& is, unsigned& line_number)
{
	m_indent = -1;
	m_line_number = -1;
	std::string str;
	while (!std::getline(is,str).eof())
	{
		m_raw_string = str;
		++line_number;
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
	if (is.eof())
	{
		return;
	}
	m_line_number = line_number;
	m_indent = str.find_first_not_of(" \t");
	if (str.find('\t')<m_indent)
	{
		throw_line_exception("No tab characters are allowed in leading whitespace.");
	}
	str = str.substr(m_indent);
	std::stringstream ss(str);
	std::string word;
	ss >> word;
	if (word.back()==':')
	{
		m_key = word.substr(0,word.size()-1);
		if (m_key.empty())
		{
			throw_line_exception("Empty key found.");
		}
		if (ss.tellg()>=0)
		{
			std::string v(str.begin()+ss.tellg(),str.end());
			m_value = parse_word(v);
		}
	}
	else if (word=="-")
	{
		if (ss.tellg()>=0)
		{
			std::string v(str.begin()+ss.tellg(),str.end());
			m_value = parse_word(v);
		}
	}
	else
	{
		throw_line_exception("expected '<key>: [value]' or '- <value>'.");
	}
}

void Lines::Line::throw_line_exception(const std::string& msg) const
{
	throw std::runtime_error("line " + std::to_string(m_line_number) + ": [" + m_raw_string + "]: " + msg);
}

json_type Lines::Line::parse_word(const std::string& str)
{
	if (str.empty())
	{
		return json_type{};
	}
	std::string s(str);
	s = s.substr(0,str.find_last_not_of(" \n\r\t")+1);
	const std::string::size_type non_ws_index = s.find_first_not_of(" \n\r\t");
	s = s.substr(non_ws_index);
	try {
		return json_type::parse(s.begin(),s.end());
	}
	catch(std::exception& exc)
	{
		return s;
	}
}

Lines::Lines(std::istream& is)
{
	unsigned line_number = 0;
	for (Line line(is,line_number);!line.empty();line=Line(is,line_number))
	{
		m_lines.push_back(line);
	}
}

/** Populate JSON object \a parent with the content that is contained in the
 *  line-range [\a it, \a end). Lines are processed (and \a it is advanced)
 *  as long as it!=end and the indentation level of the line is equal to that
 *  of the initial line that is read. A line is characterized by its key and value
 *  members, at least one of which is nonempty. The three possible cases are:
 *  must be considered are then:
 *
 *   - When the key and value are both non-empty, a simple member of the form
 *       'key': value
 *     is added to \a parent. If parent is null, it is changed into an object
 *     first.
 *   - When the key is empty, we assume that \a parent is an array. We simply
 *     push back the value to the array. Of parent is null, it is changed into
 *     an array first.
 *   - When the value is empty, a member of the form
 *       'key': null
 *     is added to \a parent. Initially, the object is empty. Subsequently,
 *     readSection is called recursively to populate the object.
 *
 *  \author Jan van Dijk
 *  \date   8 April 2024
 */
void readSection(json_type& parent, Lines::const_iterator& it, Lines::const_iterator end)
{
	Lines::const_iterator begin(it);
	while (it != end && it->indent()==begin->indent())
	{
		//std::cout << "INDENT = " << it->indent() << ", KEY = " << it->key() << ", VALUE = " << it->value() << std::endl;
		if (it->value().empty())
		{
			if (parent.contains(it->key()))
			{
				it->throw_line_exception("Duplicate key '" + it->key() + "' found.");
			}
			json_type& sec = parent[it->key()];
			auto old = it++;
			if (it==end || it->indent()==old->indent())
			{
				/// \todo should sec become an array or an object?
				// sec = json_type(R"([])");
				continue;
			}
			readSection( sec, it, end);
		}
		else if (it->key().empty())
		{
			parent.push_back(it->value());
			++it;
		}
		else
		{
			if (parent.contains(it->key()))
			{
				it->throw_line_exception("Duplicate key '" + it->key() + "' found.");
			}
			parent[it->key()] = it->value();
			++it;
		}
	}
}

} // namespace

json_type offSideToJSON(std::istream& is)
try
{
	Lines lines(is);
	json_type root;
	Lines::const_iterator it = lines.begin();
	readSection(root,it,lines.end());
	/* If there are still lines left, the indentation level
	 * if it is different from that of any of its ancestors.
	 */
	if (it!=lines.end())
	{
		it->throw_line_exception("Unexpected indentation level.");
	}
	return root;
}
catch (std::exception& exc)
{
	throw std::runtime_error(std::string("Error parsing stream: ") + exc.what());
}

json_type offSideToJSON(const std::filesystem::path& fname)
{
	std::ifstream ifs(fname);
	if (!ifs)
	{
		throw std::runtime_error("Error opening file '" + fname.generic_string() + "' for reading.");
	}
	return offSideToJSON(ifs);
}

} // namespace loki
