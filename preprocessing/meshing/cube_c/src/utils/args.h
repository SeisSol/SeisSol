/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * TODO place this in the submodule folder
 */

#ifndef UTILS_ARGS_H
#define UTILS_ARGS_H

#include <getopt.h>
#include <algorithm>
#include <map>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace utils
{

/**
 * Parses command line arguments
 * @todo comment functions
 */
class Args
{
private:
	struct optionInfo {
		std::string longOption;		// We need a copy here to get the const char* correct
		/** Name of the value in the help command */
		std::string value;
		std::string description;
		bool required;
	};

	struct additionalOptionInfo {
		std::string name;
		std::string description;
		bool required;
	};

	/**
	 * Convert a long option into an "argument" that is shown in the help message
	 */
	struct valueConvert {
		void operator()(char& c)
		{
			c = toupper(static_cast<unsigned char>(c));
			switch (c) {
			case '-':
				c = '_';
				break;
			}
		}
	};

	/** Program description (can be empty) */
	const std::string m_description;
	/** Automatically add help option */
	const bool m_addHelp;

	/** The command line arguments */
	std::vector<struct option> m_options;
	/**
	 * Additional information for the command line arguments
	 * required to generate help information
	 */
	std::vector<optionInfo> m_optionInfo;

	std::vector<additionalOptionInfo> m_additionalOptionInfo;

	/** Maps from short option to index in m_options */
	std::map<char, size_t> m_short2option;

	/** Contains the arguments after parse was called */
	std::map<std::string, std::string> m_arguments;

	/**
	 * Contains additional arguments after parse was called
	 * @todo Find a better name
	 */
	std::map<std::string, std::string> m_additionalArguments;

public:
	enum Argument
	{
		Required = required_argument,
		No = no_argument,
		Optional = optional_argument
	};

	enum Result
	{
		Success = 0,
		Error,
		/** Help message printed */
		Help
	};

public:
	Args(const std::string &description = "", bool addHelp = true)
		: m_description(description),
		  m_addHelp(addHelp)
	{
	}

	void addOption(const std::string &longOption,
			char shortOption = 0,
			const std::string &description = "",
			Argument argument = Required,
			bool required = true)
	{
		std::string value;

		if (shortOption)
			m_short2option[shortOption] = m_options.size();

		if (argument != No) {
			value = longOption;
			for_each(value.begin(), value.end(), valueConvert());
		}

		struct optionInfo i = {longOption, value, description, required};
		m_optionInfo.push_back(i);

		struct option o = {m_optionInfo.back().longOption.c_str(), argument, 0, shortOption};
		m_options.push_back(o);
	}

	void addAdditionalOption(const std::string &name,
			const std::string &description = "",
			bool required = true)
	{
		if (!m_additionalOptionInfo.empty()) {
			if (required && !m_additionalOptionInfo.back().required)
				// After one optional argument there can only be more optional arguments
				return;
		}

		struct additionalOptionInfo i = {name, description, required};
		m_additionalOptionInfo.push_back(i);
	}

	/**
	 * @return True of options are successfully parsed, false otherwise
	 */
	Result parse(int argc, char* const* argv, bool printHelp = true)
	{
		if (m_addHelp)
			addOption("help", 'h', "Show this help message", No, false);

		std::ostringstream shortOptions;
		for (std::vector<struct option>::const_iterator i = m_options.begin();
			i != m_options.end(); i++) {
			if (i->val != 0) {
				shortOptions << static_cast<char>(i->val);
				switch (i->has_arg)
				{
				case required_argument:
					shortOptions << ':';
					break;
				case optional_argument:
					shortOptions << "::";
					break;
				}
			}
		}

		// Add null option
		struct option o = {0, 0, 0, 0};
		m_options.push_back(o);

		// Update const char* in m_options
		for (size_t i = 0; i < m_optionInfo.size(); i++)
			m_options[i].name = m_optionInfo[i].longOption.c_str();

		while (true) {
			int optionIndex = 0;

			int c = getopt_long(argc, argv, shortOptions.str().c_str(),
				&m_options[0], &optionIndex);

			if (c < 0)
				break;

			switch (c) {
			case '?':
				if (printHelp)
					helpMessage(argv[0], std::cerr);
				return Error;
			case 0:
				// Nothing to do
				break;
			default:
				optionIndex = m_short2option.at(c);
			}

			if (optarg == 0L)
				m_arguments[m_options[optionIndex].name] = "";
			else
				m_arguments[m_options[optionIndex].name] = optarg;
		}

		if (m_addHelp && isSet("help")) {
			if (printHelp)
				helpMessage(argv[0]);
			return Help;
		}

		for (std::vector<optionInfo>::const_iterator i = m_optionInfo.begin();
			i != m_optionInfo.end(); i++) {
			if (i->required && !isSet(i->longOption)) {
				if (printHelp) {
					std::cerr << argv[0] << ": option --" << i->longOption << " is required" << std::endl;
					helpMessage(argv[0], std::cerr);
				}
				return Error;
			}
		}

		// Parse additional options and check if all required options are set
		int i;
		for (i = 0; i < argc-optind; i++)
			m_additionalArguments[m_additionalOptionInfo[i].name] = argv[i+optind];
		if (static_cast<size_t>(i) < m_additionalOptionInfo.size()) {
			if (m_additionalOptionInfo[i].required) {
				if (printHelp) {
					std::cerr << argv[0] << ": option <" << m_additionalOptionInfo[i].name << "> is required" << std::endl;
					helpMessage(argv[0], std::cerr);
				}
				return Error;
			}
		}

		return Success;
	}

	bool isSet(const std::string &option) const
	{
		return m_arguments.find(option) != m_arguments.end();
	}

	bool isSetAdditional(const std::string &option) const
	{
		return m_additionalArguments.find(option) != m_additionalArguments.end();
	}

	template<typename T>
	T getArgument(const std::string &option)
	{
		std::istringstream ss(m_arguments.at(option));

		T result;
		ss >> result;

		return result;
	}

	template<typename T>
	T getArgument(const std::string &option, T defaultArgument)
	{
		if (!isSet(option))
			return defaultArgument;

		return getArgument<T>(option);
	}


	template<typename T>
	T getAdditionalArgument(const std::string &option)
	{
		std::istringstream ss(m_additionalArguments.at(option));

		T result;
		ss >> result;

		return result;
	}

	template<typename T>
	T getAdditionalArgument(const std::string &option, T defaultArgument)
	{
		if (!isSetAdditional(option))
			return defaultArgument;

		return getAdditionalArgument<T>(option);
	}

	void helpMessage(const char* prog, std::ostream &out = std::cout)
	{
		// First line with all short options
		out << "Usage: " << prog;
		for (size_t i = 0; i < m_options.size()-1; i++) {
			out << ' ';

			if (!m_optionInfo[i].required)
				out << '[';

			if (m_options[i].val != 0)
				out << '-' << static_cast<char>(m_options[i].val);
			else
				out << "--" << m_options[i].name;

			argumentInfo(i, out);

			if (!m_optionInfo[i].required)
				out << ']';
		}
		for (size_t i = 0; i < m_additionalOptionInfo.size(); i++) {
			out << ' ';

			if (!m_additionalOptionInfo[i].required)
				out << '[';

			out << '<' << m_additionalOptionInfo[i].name << '>';

			if (!m_additionalOptionInfo[i].required)
				out << ']';

		}
		out << std::endl;

		// General program description
		if (!m_description.empty())
			out << std::endl << m_description << std::endl;

		// Arguments
		out << std::endl << "arguments:" << std::endl;
		for (size_t i = 0; i < m_additionalOptionInfo.size(); i++) {
			out << "  <" << m_additionalOptionInfo[i].name << '>';

			// Number of characters used for the option
			size_t length = 4 + m_additionalOptionInfo[i].name.size();

			if (length >= 30) {
				out << std::endl;
				out << std::setw(30) << ' ';
			} else
				out << std::setw(30-length) << ' ';

			out << m_additionalOptionInfo[i].description << std::endl;
		}

		// Optional arguments
		out << std::endl << "optional arguments:" << std::endl;
		for (size_t i = 0; i < m_options.size()-1; i++) {
			out << "  ";

			// Number of characters used for the option
			size_t length = 2;

			if (m_options[i].val != 0) {
				out << '-' << static_cast<char>(m_options[i].val);
				out << ", ";
				length += 4;
			}

			out << "--" << m_options[i].name;
			length += m_optionInfo[i].longOption.size() + 2;
			length += argumentInfo(i, out);

			if (length >= 30) {
				out << std::endl;
				out << std::setw(30) << ' ';
			} else
				out << std::setw(30-length) << ' ';

			out << m_optionInfo[i].description << std::endl;
		}
	}

private:
	/**
	 * Writes the argument information to out
	 *
	 * @param i The index of the option for which the argument should be generated
	 * @return The number if characters written
	 */
	size_t argumentInfo(size_t i, std::ostream &out)
	{
		switch (m_options[i].has_arg) {
		case required_argument:
			out << ' ' << m_optionInfo[i].value;
			return m_optionInfo[i].value.size() + 1;
		case optional_argument:
			out << " [" << m_optionInfo[i].value << ']';
			return m_optionInfo[i].value.size() + 3;
		}

		return 0;
	}
};

template<> inline
std::string utils::Args::getArgument(const std::string &option)
{
	return m_arguments.at(option);
}

template<> inline
std::string utils::Args::getAdditionalArgument(const std::string &option)
{
	return m_additionalArguments.at(option);
}

}

#endif // UTILS_ARGS_H
